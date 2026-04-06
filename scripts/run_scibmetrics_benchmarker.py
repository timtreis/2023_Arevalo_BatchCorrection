import sys
import logging
import numpy as np
import scanpy as sc
import pandas as pd

import scib_metrics.metrics._nmi_ari as _nmi_ari

from scib_metrics.benchmark import Benchmarker, BioConservation
from metrics.scib import (
    _ensure_inchikey,
    _merge_with_duplication,
    _load_opentargets_moa_info,
    _load_repurposinghub_moa_info,
    _load_repurposinghub_target_info,
)

import warnings

warnings.simplefilter(action="ignore", category=FutureWarning)

logger = logging.getLogger(__name__)


# --- GPU-accelerated KNN via FAISS ---

def _make_faiss_gpu_knn():
    """Create a FAISS GPU brute-force KNN function for scib-metrics.

    Uses IndexFlatL2 (exact, no approximation) to avoid the duplicate-index
    problem that HNSW can produce, which corrupts downstream metrics.
    """
    import faiss
    from scib_metrics.nearest_neighbors import NeighborsResults

    n_gpus = faiss.get_num_gpus()
    if n_gpus == 0:
        logger.warning("No GPU available for FAISS, falling back to pynndescent")
        return None

    def faiss_knn(X: np.ndarray, k: int) -> NeighborsResults:
        X = np.ascontiguousarray(X, dtype=np.float32)
        index = faiss.IndexFlatL2(X.shape[1])
        try:
            res = faiss.StandardGpuResources()
            index = faiss.index_cpu_to_gpu(res, 0, index)
            logger.info("FAISS: using GPU KNN")
        except (RuntimeError, AssertionError) as e:
            logger.warning(f"FAISS GPU KNN failed ({e}), using CPU")
        index.add(X)
        distances_sq, indices = index.search(X, k)
        return NeighborsResults(
            indices=indices,
            distances=np.sqrt(np.maximum(distances_sq, 0.0)),
        )

    return faiss_knn


try:
    _faiss_knn = _make_faiss_gpu_knn()
except ImportError:
    logger.warning("faiss not installed, falling back to pynndescent")
    _faiss_knn = None


# --- KMeans monkeypatch: FAISS GPU KMeans replaces JAX KMeans ---
# JAX KMeans creates an O(n_cells * n_clusters * n_features) 3D tensor in the
# centroid update step, causing 8+ TB allocations at 244K cells.
# FAISS GPU KMeans runs entirely on GPU and completes in seconds.

def _compute_clustering_kmeans_patched(X, n_clusters):
    import faiss
    X = np.ascontiguousarray(X, dtype=np.float32)
    try:
        kmeans = faiss.Kmeans(X.shape[1], n_clusters, niter=20, gpu=True, seed=0)
        kmeans.train(X)
    except (RuntimeError, AssertionError):
        logger.warning("FAISS GPU KMeans failed, falling back to CPU")
        kmeans = faiss.Kmeans(X.shape[1], n_clusters, niter=20, gpu=False, seed=0)
        kmeans.train(X)
    _, labels = kmeans.index.search(X, 1)
    return labels.ravel()

_nmi_ari._compute_clustering_kmeans = _compute_clustering_kmeans_patched


def run_scibmetrics_benchmarker(
    adata_path, output_path, batch_key, eval_keys, methods
):
    adata = sc.read_h5ad(adata_path)

    if isinstance(batch_key, list) and len(batch_key) == 1:
        batch_key = batch_key[0]

    # keys are whitespace separated
    eval_keys = eval_keys.split(" ")

    results = []

    # iterate over keys and benchmark embeddings
    for eval_key in eval_keys:

        eval_key_function_mapping = {
            "Metadata_OT_MOA": _load_opentargets_moa_info,
            "Metadata_DRH_MOA": _load_repurposinghub_moa_info,
            "Metadata_DRH_TARGET": _load_repurposinghub_target_info,
        }

        if eval_key in eval_key_function_mapping:
            meta = eval_key_function_mapping[eval_key]()
            adata_copy = _ensure_inchikey(adata.copy())
            adata_for_eval = _merge_with_duplication(adata_copy, meta)
            adata_for_eval = adata_for_eval[~adata_for_eval.obs[eval_key].isna()].copy()
        else:
            adata_for_eval = adata.copy()

        if eval_key not in adata_for_eval.obs.columns:
            raise ValueError(f"Eval key '{eval_key}' not in metadata")

        # With high-cardinality labels (e.g. 82K compounds), some metrics are
        # degenerate and must be skipped:
        # - cLISI: random embeddings score ~0.97, masking real signal. Also OOMs
        #   on GPU (jax.vmap over 244K cells × 82K label bincount).
        # - kmeans NMI/ARI: scib-metrics sets n_clusters = n_labels, so 82K
        #   clusters on 244K cells (~3 cells/cluster) makes KMeans unstable.
        #   Results are dominated by initialization noise, not bio conservation.
        n_labels = adata_for_eval.obs[eval_key].nunique()
        high_cardinality = n_labels > 500

        if high_cardinality:
            logger.info(
                f"High-cardinality labels ({n_labels}): disabling cLISI and "
                f"kmeans NMI/ARI (degenerate with n_clusters={n_labels})"
            )

        bm = Benchmarker(
            adata_for_eval,
            batch_key=batch_key,
            label_key=eval_key,
            embedding_obsm_keys=methods.split(" "),
            bio_conservation_metrics=BioConservation(
                clisi_knn=not high_cardinality,
                nmi_ari_cluster_labels_kmeans=not high_cardinality,
            ),
            n_jobs=-1,
        )
        if _faiss_knn is not None:
            bm.prepare(neighbor_computer=_faiss_knn)
        bm.benchmark()

        df = bm.get_results(min_max_scale=False)

        # ugly code to convert to tidy format so we can properly save to parquet
        # TODO(ttreis): maybe refactor this to be more readable
        metric_types = df.loc["Metric Type"].to_frame(name="Metric Type").reset_index()
        metric_types.columns = ["Metric", "Metric_Type"]
        df2 = df.drop("Metric Type").reset_index()
        df2 = df2.rename(columns={"Embedding": "Method"})
        df_long = df2.melt(id_vars=["Method"], var_name="Metric", value_name="Value")
        df_tidy = df_long.merge(metric_types, on="Metric", how="left")
        df_tidy = df_tidy.query("Metric_Type != 'Aggregate score'")
        df_tidy["eval_key"] = eval_key

        results.append(df_tidy)

    # make pretty
    df_tidy = pd.concat(results)
    df_tidy = df_tidy[["eval_key", "Method", "Metric", "Value", "Metric_Type"]]
    df_tidy.columns = df_tidy.columns.str.replace(r"\s+", "_", regex=True).str.lower()
    cols_to_modify = ["method", "metric", "metric_type"]
    df_tidy[cols_to_modify] = df_tidy[cols_to_modify].map(
        lambda col: col.replace(" ", "_").lower() if isinstance(col, str) else col
    )

    df_tidy.reset_index(drop=True).to_parquet(output_path, index=False)


if __name__ == "__main__":
    adata_path = sys.argv[1]
    output_path = sys.argv[2]
    batch_key = sys.argv[3]
    eval_keys = sys.argv[4]
    methods = sys.argv[5]
    run_scibmetrics_benchmarker(adata_path, output_path, batch_key, eval_keys, methods)
