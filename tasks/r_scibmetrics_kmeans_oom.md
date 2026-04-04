# scib-metrics Scalability Research: JAX KMeans OOM and Alternatives

**Date**: 2026-04-02
**Context**: scib-metrics v0.5.6 Benchmarker on 244K cells x 9 methods, 80GB H100 GPU.

---

## 1. Why Does JAX KMeans Need 8.3TB for 244K Cells?

### Root Cause: O(n_cells x n_clusters x n_features) 3D Tensor

The memory explosion is in `_kmeans_step` inside `scib_metrics/utils/_kmeans.py` (lines 152-170). The centroid update computes:

```python
new_centroids = (
    jnp.sum(
        jnp.where(
            # Shape: (n_cells, n_clusters, n_features)
            new_labels[:, jnp.newaxis, jnp.newaxis]
            == jnp.arange(self.n_clusters)[jnp.newaxis, :, jnp.newaxis],
            X[:, jnp.newaxis, :],
            0.0,
        ),
        axis=0,
    )
    / counts
)
```

This creates a **dense 3D intermediate tensor** of shape `(n_cells, n_clusters, n_features)` in float32 (4 bytes). Instead of using scatter/segment-sum operations (which would be O(n_cells x n_features)), it broadcasts the full data matrix against all clusters.

### Memory Math

For 244K cells with various label counts and embedding dimensions:

| n_clusters | n_features | Single Tensor | Notes |
|-----------|------------|---------------|-------|
| 50 | 50 | 2.3 GB | Tight fit on GPU with overhead |
| 100 | 50 | 4.5 GB | Still possible |
| 100 | 100 | 9.1 GB | Over budget with JIT overhead |
| 200 | 50 | 9.1 GB | Common for Cell Painting |
| 200 | 200 | 36.4 GB | Half an H100 for ONE tensor |
| 500 | 50 | 22.7 GB | High-cardinality labels |

**But the 8.3TB figure is not from a single tensor.** It arises because:

1. **XLA JIT compilation** of `jax.lax.while_loop` traces the full computation graph. The while_loop body contains both this 3D tensor AND distance computations.
2. **XLA buffer assignment** must account for all live buffers simultaneously during compilation planning.
3. **`_initialize_plus_plus` uses `jax.lax.scan`** over `n_clusters-1` steps, each computing `cdist(n_local_trials, n_cells)`. The scan stores all intermediate states, creating O(n_clusters x n_local_trials x n_cells) intermediates.
4. **Convergence bug (still present in v0.5.6)**: The while_loop condition uses `jnp.logical_or(cond1, cond2)` instead of `jnp.logical_and`, meaning KMeans always runs all `max_iter=300` iterations regardless of convergence. XLA must plan memory for the full unrolled loop.

The 8.3TB figure comes from XLA's buffer planner trying to allocate for the compiled HLO graph, which compounds the 3D tensor across loop iterations plus k-means++ initialization buffers.

---

## 2. Upstream Status

### Relevant GitHub Issues

| Issue | Status | Summary |
|-------|--------|---------|
| [#43](https://github.com/YosefLab/scib-metrics/issues/43) | Closed | KMeans++ memory -- fixed by switching map vs vmap (v0.0.7) |
| [#50](https://github.com/YosefLab/scib-metrics/issues/50) | Closed | Large memory issue for ARI/NMI with KMeans -- partial fix |
| [#85](https://github.com/YosefLab/scib-metrics/issues/85) | Closed | Large-scale benchmarking slow on 892K cells |
| [#114](https://github.com/YosefLab/scib-metrics/issues/114) | Closed | KMeans results differ from sklearn (convergence + init bugs, fixed v0.4.1) |
| [#186](https://github.com/YosefLab/scib-metrics/issues/186) | Closed | GPU memory not de-allocated on 2.7M cells (closed with gc.collect suggestion) |

### Latest Version (0.5.9, Feb 2026)

Versions 0.5.7-0.5.9 contain **no changes to KMeans**. Changes were: kBET numpy compat (0.5.9), Leiden seed fix (0.5.8), arpack solver selection (0.5.7). The fundamental O(n x k x d) 3D tensor in centroid update **remains unfixed on main**. The code on GitHub main is identical to our v0.5.6.

### Conclusion

The OOM is a known architectural limitation that has not been addressed upstream. The maintainers' approach is to suggest Leiden clustering or subsampling.

---

## 3. Alternatives for NMI/ARI Computation at Scale

### Option A: Use Leiden Instead of KMeans (Recommended -- Simplest)

`nmi_ari_cluster_labels_leiden` is already implemented in scib-metrics and works on the **kNN graph** (sparse, O(n x k_neighbors) memory) rather than the full embedding matrix.

- **Memory**: O(n x k_neighbors) for the graph. 244K x 90 neighbors x 8 bytes = ~176 MB.
- **Scalability**: Proven at 900K+ cells in the scib-metrics large-scale tutorial.
- **Semantics**: Different from KMeans (graph-based vs. centroid-based), but this is what the original scIB paper (Luecken et al. 2022) used (Louvain, which Leiden replaces).
- **Integration**: Already in Benchmarker. Set `nmi_ari_cluster_labels_kmeans=False, nmi_ari_cluster_labels_leiden=True` in BioConservation.

### Option B: sklearn MiniBatchKMeans (Drop-in KMeans Replacement)

Monkeypatch `_compute_clustering_kmeans` with sklearn's MiniBatchKMeans:

```python
from sklearn.cluster import MiniBatchKMeans
import scib_metrics.metrics._nmi_ari as _nmi_ari

def _compute_clustering_kmeans_patched(X, n_clusters):
    return MiniBatchKMeans(
        n_clusters=n_clusters, batch_size=10000, n_init=10, random_state=0
    ).fit_predict(X)

_nmi_ari._compute_clustering_kmeans = _compute_clustering_kmeans_patched
```

- **Memory**: O(batch_size x n_features + n_clusters x n_features). ~4 MB for batch=10K, d=50.
- **Speed**: Much faster than full KMeans for large n.
- **Downside**: Approximate, but fine for NMI/ARI benchmarking.

### Option C: FAISS KMeans (GPU-Optimized)

```python
import faiss

def _compute_clustering_kmeans_faiss(X, n_clusters):
    X = np.ascontiguousarray(X, dtype=np.float32)
    kmeans = faiss.Kmeans(X.shape[1], n_clusters, niter=20, gpu=True, verbose=False)
    kmeans.train(X)
    _, labels = kmeans.index.search(X, 1)
    return labels.ravel()
```

- **Memory**: FAISS pages data to GPU as needed. Handles millions of points.
- **Speed**: ~17x faster than sklearn on GPU.
- **Downside**: Extra dependency (but already used for kNN in this pipeline).

### Option D: Pre-compute Labels Externally

Compute KMeans labels separately, skip KMeans in Benchmarker, and compute NMI/ARI manually:

```python
from sklearn.metrics.cluster import normalized_mutual_info_score, adjusted_rand_score

labels_pred = my_scalable_kmeans(X, n_clusters)
nmi = normalized_mutual_info_score(true_labels, labels_pred, average_method="arithmetic")
ari = adjusted_rand_score(true_labels, labels_pred)
```

Set `nmi_ari_cluster_labels_kmeans=False` in Benchmarker and merge results post-hoc.

### Recommendation

**Option A (Leiden)** for simplicity and alignment with original scIB methodology. If KMeans NMI/ARI comparability is specifically needed, **Option B (MiniBatchKMeans monkeypatch)** is the best balance of correctness and effort.

---

## 4. Which Other Metrics Scale Poorly to 200K+?

| Metric | Uses JAX? | Memory Complexity | 244K Status | 1M+ Risk |
|--------|----------|-------------------|-------------|----------|
| **KMeans NMI/ARI** | Yes (full JAX) | O(n x k x d) | **OOM on GPU** | OOM on CPU too |
| **Silhouette label/batch** | Yes (chunked cdist) | O(chunk x n) per chunk | OK (250 MB/chunk with chunk_size=256) | Slow (~4K chunks) but feasible |
| **iLISI/cLISI** | Yes (vmap over cells) | O(n x k_neighbors) | OK at 244K (~176 MB) | **OOM at 2.7M** (issue #186) |
| **kBET** | Yes (vmap, per-label) | O(n_label x k_neighbors) | OK (processes per-label subsets) | OK (subsets keep it bounded) |
| **Graph connectivity** | No (igraph) | O(edges) | OK | OK |
| **PCR comparison** | No (sklearn PCA) | O(n x d) | OK | OK |
| **Leiden NMI/ARI** | No (igraph) | O(edges) | OK | OK |

**Key finding**: For 244K cells, only KMeans NMI/ARI is problematic. At >1M cells, iLISI/cLISI would also become problematic.

---

## 5. What Do Other Papers/Benchmarks Do at Scale?

### Original scIB Paper (Luecken et al. 2022, Nature Methods)

- Used **Louvain clustering** (not KMeans) for NMI/ARI.
- scib-metrics switched to KMeans as an "optimization" which backfired.
- scIB's `metrics()` function has a `subsample` parameter (default 0.5) that subsamples cells.

### scib-metrics Large-Scale Tutorial

- Uses ~900K cells (lung atlas), **subsampled to 500K**.
- Uses default BioConservation (KMeans enabled).
- Uses FAISS for kNN on GPU.
- Reports 20 minutes total on RTX 3090.
- Likely runs KMeans on CPU by default (tutorial conda env does not set up JAX GPU).

### Deep Learning Integration Benchmark (2025, Genome Biology)

- Large-scale benchmarks typically subsample to 50K-200K.
- Use the scIB framework metrics.

### Summary

**Most benchmarks at >200K cells either subsample or use CPU-based clustering.** GPU JAX KMeans in scib-metrics is not designed for >100K cells with the current 3D tensor approach.

---

## 6. Correct Long-Term Fix

The proper fix (for an upstream PR to scib-metrics) would replace the 3D broadcast centroid update with `jax.ops.segment_sum`:

```python
# Current: O(n_cells x n_clusters x n_features) memory
new_centroids = jnp.sum(jnp.where(...3D tensor...), axis=0) / counts

# Fixed: O(n_cells x n_features) memory
new_centroids = jax.ops.segment_sum(X, new_labels, num_segments=self.n_clusters) / counts
```

This reduces memory from O(n x k x d) to O(n x d) + O(k x d), making GPU KMeans feasible at any scale. The convergence bug (`logical_or` vs `logical_and` in the while_loop condition) should also be reported/fixed.
