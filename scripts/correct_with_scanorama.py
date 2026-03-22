import argparse
import logging
import scanpy as sc
from preprocessing import io
import anndata as ad

logger = logging.getLogger(__name__)


def _load_hyperparameters(parameter_path):
    import pandas as pd
    params = pd.read_csv(parameter_path)
    params = params.sort_values("total", ascending=False).iloc[0].to_dict()
    if params["state"] != "COMPLETE":
        raise ValueError("Optimization did not complete successfully")
    knn = int(params["params_knn"])
    sigma = params["params_sigma"]
    alpha = params["params_alpha"]
    print(f"\nUsing tuned params: knn={knn}, sigma={sigma}, alpha={alpha}\n")
    return knn, sigma, alpha


def correct_with_scanorama(parquet_path, batch_key, output_path, parameter_path=None, smoketest=False):
    """Scanorama correction on raw data"""
    import scanorama

    adata = io.to_anndata(parquet_path)

    knn, sigma, alpha = (20, 15.0, 0.10)  # defaults
    if parameter_path is not None:
        knn, sigma, alpha = _load_hyperparameters(parameter_path)

    # Sort adata based on batch_key (needed for scanorama)
    adata = adata[adata.obs.sort_values(by=batch_key).index]

    def split_adata_by_col(adata, col):
        """Split the AnnData object by a column"""
        splits = []
        for value in adata.obs[col].cat.categories:
            mask = adata.obs[col] == value
            splits.append(adata[mask].copy())
        return splits

    # batch-correct
    adatas_by_source = split_adata_by_col(adata, batch_key)
    scanorama.integrate_scanpy(adatas_by_source, knn=knn, sigma=sigma, alpha=alpha)
    corrected_adata = ad.concat(adatas_by_source)

    # write to parquet as old pipeline
    meta = corrected_adata.obs.reset_index(drop=True).copy()
    vals = corrected_adata.obsm["X_scanorama"]
    features = [f"scanorama_{i}" for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)


def correct_with_scanorama_pca(parquet_path, batch_key, output_path, smoketest=False):
    """Scanorama correction using PCA components"""
    adata = io.to_anndata(parquet_path)

    # Sort adata based on batch_key (needed for scanorama)
    adata = adata[adata.obs.sort_values(by=batch_key).index]

    # TODO(ttreis): Do these numbers make sense?
    n_comps = 2 if smoketest else 50

    sc.pp.pca(adata, svd_solver="arpack", n_comps=n_comps)
    sc.external.pp.scanorama_integrate(
        adata, key=batch_key, basis="X_pca", adjusted_basis="X_scanorama"
    )

    # write to parquet as old pipeline
    meta = adata.obs.reset_index(drop=True).copy()
    vals = adata.obsm["X_scanorama"]
    features = [f"scanorama_pca_{i}" for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Perform Scanorama correction on data."
    )

    parser.add_argument(
        "--mode",
        required=True,
        choices=["scanorama", "scanorama_pca"],
        help="Correction mode to use: 'scanorama' or 'scanorama_pca'.",
    )
    parser.add_argument(
        "--input_data", required=True, help="Path to the input data in Parquet format."
    )
    parser.add_argument("--batch_key", required=True, help="Batch key.")
    parser.add_argument("--parameter_path", default=None, help="Path to Optuna parameter CSV.")
    parser.add_argument(
        "--output_path", required=True, help="Path to save the corrected data."
    )
    parser.add_argument(
        "--smoketest",
        action="store_true",
        help="Run a smoketest with reduced computation for testing purposes.",
    )

    args = parser.parse_args()

    if args.mode == "scanorama":
        correct_with_scanorama(
            parquet_path=args.input_data,
            batch_key=args.batch_key,
            output_path=args.output_path,
            parameter_path=args.parameter_path,
            smoketest=args.smoketest,
        )
    elif args.mode == "scanorama_pca":
        correct_with_scanorama_pca(
            parquet_path=args.input_data,
            batch_key=args.batch_key,
            output_path=args.output_path,
            smoketest=args.smoketest,
        )
    else:
        raise ValueError("Invalid mode. Choose either 'scanorama' or 'scanorama_pca'.")
