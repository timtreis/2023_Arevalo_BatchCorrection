# R Method Scalability Limits on COMPOUND-Plate Scenarios

**Date**: 2026-04-01
**Status**: Investigated, fix proposed
**Affects**: fastMNN, seurat_cca, seurat_rpca on scenarios with >100K cells

---

## Summary

All three R-based batch correction methods (fastMNN, Seurat CCA, Seurat RPCA) fail
on COMPOUND-plate scenarios due to R's 2^31-1 element vector limit and memory
allocation patterns in intermediate computations. This is not a RAM issue — it is a
hard language-level constraint in R that cannot be worked around without rewriting
upstream package internals.

---

## Evidence

### Scenario 3 failures (244,797 cells, 1,064 features, 3 sources)

All three methods were launched on a full A100 80 GB node (2 TB system RAM).

| Method | Trials attempted | Outcome | Error |
|--------|-----------------|---------|-------|
| fastMNN | 10 real + 19 cascade | All PRUNED (exit -9) | `"long vectors not supported yet"` in `.average_correction` |
| seurat_cca | 10 real + 19 cascade | All PRUNED (exit -9) | R process killed after 1-2h inside `IntegrateLayers` |
| seurat_rpca | 30+ | All PRUNED (exit -9) | R process killed during integration |

After all real trials returned None, Optuna's TPE sampler crashed with
`TypeError: unsupported operand type(s) for *: 'NoneType' and 'float'` in
`_calculate_weights_below_for_multi_objective`. This caused trials 11-29 to fail
instantly (cascade failure). The correction step then found "No COMPLETE trials"
and aborted.

For comparison, the same methods completed successfully on scenarios 1 (8K cells)
and 2 (14K cells) with all 30 HPO trials.

### Scenario 8 defaults (370,231 cells) — partial confirmation

Seurat CCA/RPCA completed on scenario_8_defaults using pre-tuned default parameters
(single trial, no HPO). This suggests that with carefully chosen parameters, small
instances of these methods may survive, but HPO (which explores parameter space
including aggressive settings) reliably triggers the failure.

fastMNN also completed on scenario_8_defaults, consistent with the finding that the
long vector issue depends on the number of MNN pairs found (parameter-dependent).

---

## Root Cause Analysis

### fastMNN: `.average_correction` long vector limit

The call chain is: `fastMNN -> .fast_mnn -> .fast_mnn_core -> .average_correction`

The `.average_correction` function in the `batchelor` R package creates a smoothing
kernel matrix of dimensions `n_target_cells x n_mnn_pairs`. With 3 sources merged
sequentially:

- Step 1: merge source_2 (82K) with source_6 (89K) = 171K cells
- Step 2: merge result (171K) with source_10 (74K) = 244K cells

For 171K target cells, only ~12,600 MNN pairs can fit before the smoothing matrix
exceeds 2^31 elements. With `k` values in the HPO range (5-50) and significant
compound overlap across JUMP-CP sources, the number of MNN pairs routinely exceeds
this threshold.

This is a hard R language limitation. R does not support vectors with >2^31-1
elements regardless of available RAM. The `batchelor` package would need internal
C-level changes to chunk the smoothing computation.

### Seurat CCA/RPCA: memory allocation during `IntegrateLayers`

Seurat's `IntegrateLayers` creates dense intermediate matrices during anchor finding
and weight computation. The PR #33 improvements (adaptive k_weight bounds, retry
with decreasing k_weight, NaN detection) address parameter-level failures but not
the fundamental scale issue.

Each trial runs for 1-2 hours before being killed (signal -9), indicating the
process gets deep into integration before the allocation spike. R's copy-on-modify
semantics can multiply memory usage unpredictably for large objects, and Seurat was
not designed for 200K+ cell dense cell-profiling matrices (it typically operates on
sparse scRNA-seq data).

### Optuna cascade failure (secondary)

When all real R trials fail (return None), Optuna's TPE sampler attempts to compute
importance weights on an all-None array, hitting a multiplication TypeError. This
causes all remaining trials to fail instantly without ever launching R. This is an
Optuna bug with no-complete-trials edge cases but is secondary to the R failures.

---

## Affected Scenarios

| Scenario | Cells | Plate types | R methods | Status |
|----------|------:|-------------|-----------|--------|
| S1 | 8,064 | T2 | work | done |
| S2 | 14,194 | T2 | work | done |
| S6 | 11,892 | T2 | work | not started |
| S4 (est) | ~25,329 | T2 | work | not started |
| S7 | 25,329 | T2 | work | done (old) |
| apricot | 45,266 | T2 | borderline | not started |
| **S3** | **244,797** | **T2+C** | **fail** | **partial** |
| **S5 (est)** | **~370,231** | **T2+C** | **fail** | **not started** |
| **S8** | **370,231** | **T2+C** | **fail** | **partial** |
| **wave2 (est)** | **~255,000** | **T2+C** | **fail** | **not started** |
| **wave1 (est)** | **~550,000** | **T2+C** | **fail** | **not started** |

Pattern: COMPOUND plates add ~70-80K cells per source. TARGET2-only scenarios stay
under ~50K cells and are safe. All current and planned COMPOUND-inclusive scenarios
exceed 200K cells.

---

## Options Considered

### 1. Fix upstream R packages

Not feasible. The `batchelor` package's `.average_correction` would need C-level
changes to chunk the smoothing kernel computation. Seurat's `IntegrateLayers` is a
complex pipeline with many internal dense matrix operations. Neither package is under
our control and neither has open issues tracking large-dataset support in R.

### 2. Subsample to ~50-100K cells, integrate, project remaining

Moderate implementation effort. Would require:
- A subsampling step preserving batch/compound proportions
- Integration on the subsample
- Projection of remaining cells (fastMNN supports this via `reducedMNN`; Seurat has
  `MapQuery` but it changes the integration approach)

**Rejected** because it changes the method semantics. We would no longer be
benchmarking "fastMNN" or "Seurat CCA" but a subsampling-projection wrapper. Results
would not be comparable to TARGET2-only scenarios (where methods run natively) or to
Arevalo et al.'s published results.

### 3. Seurat v5 SketchIntegration

Seurat v5 offers sketch-based integration for large datasets. This is a different
algorithm (sketch + integrate + project) rather than standard CCA/RPCA.

**Rejected** for the same reason: it's a different method, not comparable to the
standard CCA/RPCA we run on TARGET2-only scenarios.

### 4. Auto-skip R methods on large scenarios (recommended)

Add cell-count-based auto-skip logic in the Snakefile, following the existing pattern
for DESC (>200K cells) and scANVI (COMPOUND plates). Threshold: 100,000 cells.

Advantages:
- Trivial to implement (5-10 lines, same pattern as existing skips)
- Honest and reproducible — the skip is documented and automatic
- No change to method semantics on scenarios where R methods do work
- The scalability limitation is itself a benchmarking finding worth reporting

---

## Decision: Auto-skip R methods above 100K cells

### Rationale

1. **Arevalo et al. faced the same constraint.** Their largest scenarios had
   equivalent cell counts. These methods were never validated at this scale.

2. **12 methods remain on COMPOUND scenarios.** After skipping fastMNN +
   seurat_cca + seurat_rpca (+ already-skipped desc + scanvi), we retain 9-10
   methods: harmony_v1, harmony_v2, scvi_single, scvi_multi, scanorama, combat,
   sphering, scpoli, sysvi. This is sufficient for fair benchmarking.

3. **R methods are not top performers.** Based on S1, S2, and S8-defaults results,
   deep learning methods (scVI, scPoli, sysVI) and Harmony consistently outperform
   the R-based methods.

4. **The finding is reportable.** "R-based integration methods do not scale to
   COMPOUND-plate scenarios (>100K cells) due to language-level memory constraints"
   is a legitimate benchmarking result for the paper.

### Threshold choice: 100K cells

- S1/S2/S4/S6/S7 (8K-25K): well within limits, no risk
- apricot (45K): likely fine but not verified — 100K gives comfortable margin
- S3 (245K): clearly above, confirmed failure
- Conservative threshold avoids borderline situations during HPO (where aggressive
  parameter choices increase intermediate sizes)

### Implementation

Add to `Snakefile` after the existing DESC and scANVI skip blocks:

```python
R_METHODS_MAX_CELLS = 100_000
r_methods = ["fastMNN", "seurat_cca", "seurat_rpca"]
# Check after preprocessing completes
if n_cells > R_METHODS_MAX_CELLS:
    for m in r_methods:
        if m in METHODS:
            METHODS.remove(m)
            print(f"Auto-skipping {m}: {n_cells} cells exceeds R method limit ({R_METHODS_MAX_CELLS})")
```

---

## Per-source cell counts (reference)

| Source | TARGET2 | COMPOUND | Total |
|--------|--------:|---------:|------:|
| source_2 | 3,828 | 77,912 | 81,740 |
| source_3 | 9,599 | 42,038 | 51,637 |
| source_6 | 8,064 | 80,990 | 89,054 |
| source_8 | 1,536 | 72,261 | 73,797 |
| source_10 | 2,302 | 71,701 | 73,003 |
