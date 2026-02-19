import os
import sys
import logging
import anndata as ad
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

logger = logging.getLogger(__name__)

def scib_benchmark_embedding(
    adata: ad.AnnData,
    batch_key: str,
    label_key: str,
) -> float:
    adata.obsm["trial"] = adata.X

    # silence output
    sys.stdout = open(os.devnull, "w")

    bm = Benchmarker(
        adata=adata,
        batch_key=batch_key,
        label_key=label_key,
        embedding_obsm_keys=["trial"],
        bio_conservation_metrics=BioConservation(),
        batch_correction_metrics=BatchCorrection(),
    )
    bm.benchmark()
    df = bm.get_results(min_max_scale=False)

    # restore output
    sys.stdout.close()
    sys.stdout = sys.__stdout__

    return df.loc["trial"][["Batch correction", "Bio conservation"]].values

