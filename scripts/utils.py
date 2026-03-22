import os
import sys
import logging
import anndata as ad

logger = logging.getLogger(__name__)

def scib_benchmark_embedding(
    adata: ad.AnnData,
    batch_key: str,
    label_key: str,
    lightweight: bool = False,
) -> float:
    from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

    adata.obsm["trial"] = adata.X

    if lightweight:
        # Fast proxy metrics for HPO — skip expensive kNN/clustering steps
        bio = BioConservation(
            isolated_labels=False,
            nmi_ari_cluster_labels_kmeans=False,
            nmi_ari_cluster_labels_leiden=False,
            silhouette_label=True,
            clisi_knn=False,
        )
        # Handle API differences across scib-metrics versions:
        # newer versions use 'bras', older use 'silhouette_batch'
        import inspect
        batch_params = inspect.signature(BatchCorrection).parameters
        silhouette_key = "bras" if "bras" in batch_params else "silhouette_batch"
        batch = BatchCorrection(
            **{silhouette_key: True},
            ilisi_knn=False,
            kbet_per_label=False,
            graph_connectivity=False,
            pcr_comparison=True,
        )
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

