import os
import sys
import logging
import anndata as ad
import numpy as np
import pandas as pd
import time as _time

logger = logging.getLogger(__name__)

# Minimum batch coverage for a label to be used in semi-supervised training.
# Labels appearing in fewer batches are marked "Unknown" (unlabeled).
# See scANVI (Xu et al. 2021) and scPoli (De Donno et al. 2023) for details
# on semi-supervised handling of unlabeled cells.
SEMISUP_MIN_BATCHES = 5
SEMISUP_UNLABELED = "Unknown"


def coarsen_labels(adata, label_key, batch_key, min_batches=SEMISUP_MIN_BATCHES):
    """Mark rare labels as unlabeled for semi-supervised methods.

    Compounds appearing in fewer than `min_batches` distinct batches provide
    insufficient cross-batch signal for semi-supervised learning. They are
    relabeled as "Unknown" so that scANVI/scPoli treat them as unlabeled
    during training.  All cells are retained; only the training label changes.

    Parameters
    ----------
    adata
        AnnData object (modified in place).
    label_key
        Column in ``adata.obs`` with the original labels.
    batch_key
        Column(s) in ``adata.obs`` used as the batch indicator.
        If a list, the first element is used to count batch coverage.
    min_batches
        Minimum number of distinct batches a label must appear in.

    Returns
    -------
    n_kept, n_total
        Number of labels kept and total number of original labels.
    """
    if isinstance(batch_key, list):
        batch_col = batch_key[0]
    elif "," in batch_key:
        batch_col = batch_key.split(",")[0]
    else:
        batch_col = batch_key

    obs = adata.obs
    n_batches = obs[batch_col].nunique()
    min_batches = max(3, min(min_batches, n_batches))
    batches_per_label = obs.groupby(label_key, observed=True)[batch_col].nunique()
    keep_labels = set(batches_per_label[batches_per_label >= min_batches].index)

    n_total = batches_per_label.shape[0]
    n_kept = len(keep_labels)

    if n_kept < n_total:
        mask = ~obs[label_key].isin(keep_labels)
        n_relabeled = mask.sum()
        # Use a plain string column (not categorical) to avoid issues with
        # scvi-tools / scarches category handling.
        adata.obs[label_key] = obs[label_key].astype(str)
        adata.obs.loc[mask, label_key] = SEMISUP_UNLABELED
        logger.info(
            "Label coarsening: %d/%d labels kept (≥%d batches), "
            "%d/%d cells marked as '%s'",
            n_kept, n_total, min_batches,
            n_relabeled, adata.n_obs, SEMISUP_UNLABELED,
        )
    else:
        logger.info(
            "Label coarsening: all %d labels appear in ≥%d batches, no filtering needed",
            n_total, min_batches,
        )

    return n_kept, n_total


# Re-export from optuna_utils so existing imports keep working
from optuna_utils import save_optuna_results  # noqa: F401

# Maximum cells for lightweight HPO evaluation. Above this threshold, we
# stratified-subsample to keep all scib-metrics (silhouette, kNN, etc.)
# feasible. 30K gives ~2-3 min per evaluation with full metrics.
_LIGHTWEIGHT_MAX_CELLS = 30_000

# Track whether warmup has been done this process
_WARMUP_DONE = False


def _stratified_subsample(adata, keys, max_cells, seed=42):
    """Subsample adata preserving proportions of batch/label combinations."""
    rng = np.random.RandomState(seed)
    if adata.n_obs <= max_cells:
        return adata

    # Group by the combination of keys
    group_col = adata.obs[keys[0]].astype(str)
    for k in keys[1:]:
        group_col = group_col + "||" + adata.obs[k].astype(str)

    indices = []
    groups = group_col.unique()
    # Proportional allocation per group
    for g in groups:
        g_idx = np.where(group_col.values == g)[0]
        n_sample = max(1, int(len(g_idx) * max_cells / adata.n_obs))
        n_sample = min(n_sample, len(g_idx))
        indices.append(rng.choice(g_idx, size=n_sample, replace=False))

    indices = np.concatenate(indices)
    # Trim or pad to exactly max_cells
    if len(indices) > max_cells:
        indices = rng.choice(indices, size=max_cells, replace=False)

    return adata[np.sort(indices)].copy()


def _enable_jax_gpu():
    """Try to enable JAX GPU backend for faster scib-metrics computation.

    Falls back to CPU if GPU init fails (e.g. cuDNN version mismatch between
    the container and the host driver injected via --nv).
    """
    try:
        os.environ.pop("JAX_PLATFORMS", None)
        os.environ["JAX_PLATFORMS"] = "cuda,cpu"
        import jax
        jax.config.update("jax_platforms", "cuda,cpu")
        devices = jax.devices()
        has_gpu = any(d.platform == "gpu" for d in devices)
        if has_gpu:
            # Smoke-test: run a trivial op to verify cuDNN actually works.
            import jax.numpy as jnp
            _ = float(jnp.ones(1).sum())
            print(f"JAX GPU enabled: {devices}", file=sys.__stderr__, flush=True)
        return has_gpu
    except Exception as exc:
        # cuDNN mismatch or other GPU failure — fall back to CPU-only JAX.
        print(f"JAX GPU failed ({exc}), falling back to CPU",
              file=sys.__stderr__, flush=True)
        os.environ["JAX_PLATFORMS"] = "cpu"
        import jax
        jax.config.update("jax_platforms", "cpu")
        return False


def warmup_benchmark(batch_key, label_key):
    """Run a tiny dummy benchmark to warm up JAX JIT cache and scib-metrics.

    Call this once before the Optuna trial loop. It triggers all JIT compilations,
    import-time work, and GPU kernel caching so that actual trials run fast.
    """
    global _WARMUP_DONE
    if _WARMUP_DONE:
        return
    _WARMUP_DONE = True

    t0 = _time.time()
    print("WARMUP: running dummy benchmark to pre-compile JAX kernels...",
          file=sys.__stderr__, flush=True)

    _enable_jax_gpu()

    # Create tiny dummy data with same structure as real data
    n = 200
    dummy = ad.AnnData(
        X=np.random.randn(n, 50).astype(np.float32),
        obs=pd.DataFrame({
            batch_key: np.random.choice(["b1", "b2"], n),
            label_key: np.random.choice(["l1", "l2", "l3"], n),
        }),
    )

    # Run the full benchmark pipeline to trigger all JIT compilations
    scib_benchmark_embedding(dummy, batch_key, label_key, lightweight=True)

    t1 = _time.time()
    print(f"WARMUP: done in {t1-t0:.1f}s", file=sys.__stderr__, flush=True)


def scib_benchmark_embedding(
    adata: ad.AnnData,
    batch_key: str,
    label_key: str,
    lightweight: bool = False,
) -> float:
    from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

    if lightweight and adata.n_obs > _LIGHTWEIGHT_MAX_CELLS:
        print(
            f"SUBSAMPLE: {adata.n_obs} -> {_LIGHTWEIGHT_MAX_CELLS} cells",
            file=sys.__stderr__, flush=True,
        )
        adata = _stratified_subsample(
            adata, [batch_key, label_key], _LIGHTWEIGHT_MAX_CELLS,
        )
        print(
            f"SUBSAMPLE DONE: {adata.n_obs} cells",
            file=sys.__stderr__, flush=True,
        )

    # Enable GPU for JAX-based metrics
    if lightweight:
        _enable_jax_gpu()

    adata.obsm["trial"] = adata.X

    if lightweight:
        # Fast metrics for HPO. Disabled: kBET (repeated sampling), KMeans and
        # Leiden clustering (expensive on 30K cells), isolated_labels (uses
        # silhouette internally). Kept: silhouette_label, silhouette_batch/bras,
        # iLISI, cLISI, graph_connectivity, PCR — 6 metrics covering both
        # batch correction and bio conservation.
        bio = BioConservation(
            nmi_ari_cluster_labels_kmeans=False,
            nmi_ari_cluster_labels_leiden=False,
            isolated_labels=False,
        )
        batch = BatchCorrection(
            kbet_per_label=False,
        )
        # Increase silhouette chunk_size so the full distance matrix is computed
        # in one pass instead of 117 small chunks (default 256). On GPU or with
        # 30K subsampled cells, one 30K×30K matrix fits easily in memory (~7 GB).
        import scib_metrics.utils._silhouette as _sil
        _orig_silhouette = _sil.silhouette_samples
        def _fast_silhouette(*args, chunk_size=256, **kwargs):
            return _orig_silhouette(*args, chunk_size=max(chunk_size, adata.n_obs), **kwargs)
        _sil.silhouette_samples = _fast_silhouette
    else:
        bio = BioConservation()
        batch = BatchCorrection()

    # silence output
    sys.stdout = open(os.devnull, "w")

    bm = Benchmarker(
        adata=adata,
        batch_key=batch_key,
        label_key=label_key,
        embedding_obsm_keys=["trial"],
        bio_conservation_metrics=bio,
        batch_correction_metrics=batch,
    )
    bm.benchmark()
    df = bm.get_results(min_max_scale=False)

    # restore output
    sys.stdout.close()
    sys.stdout = sys.__stdout__

    return df.loc["trial"][["Batch correction", "Bio conservation"]].values
