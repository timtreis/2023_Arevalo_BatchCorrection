import os
import sys
import logging
import argparse
import pandas as pd
import anndata as ad
import optuna
import numpy as np
import scanorama

from preprocessing import io
from utils import scib_benchmark_embedding, save_optuna_results

logger = logging.getLogger(__name__)

# Number of PCs for HPO trials. Reduces 1040 features to 50 dims, making
# scanorama correction much faster. The final correction uses full features.
_HPO_N_PCS = 50


def objective(
    trial,
    adata: ad.AnnData,
    batch_key: str,
    label_key: str,
    smoketest: bool = False,
):
    sys.stdout = open(os.devnull, "w")

    knn = trial.suggest_int("knn", 5, 100)
    sigma = trial.suggest_float("sigma", 5.0, 150.0)
    alpha = trial.suggest_float("alpha", 0.05, 0.50)

    adata = adata[adata.obs.sort_values(by=batch_key).index].copy()

    def split_adata_by_col(adata, col):
        splits = []
        for value in adata.obs[col].cat.categories:
            mask = adata.obs[col] == value
            splits.append(adata[mask].copy())
        return splits

    adatas_by_source = split_adata_by_col(adata, batch_key)
    scanorama.integrate_scanpy(
        adatas_by_source, knn=knn, sigma=sigma, alpha=alpha
    )
    corrected_adata = ad.concat(adatas_by_source)

    vals = corrected_adata.obsm["X_scanorama"]
    features = [f"scanorama_{i}" for i in range(vals.shape[1])]
    integrated_adata = ad.AnnData(
        X=pd.DataFrame(vals, columns=features, index=corrected_adata.obs_names),
        obs=corrected_adata.obs.copy(),
    )

    batch, bio = scib_benchmark_embedding(
        adata=integrated_adata,
        batch_key=batch_key,
        label_key=label_key,
        lightweight=True,
    )

    sys.stdout.close()
    sys.stdout = sys.__stdout__

    return batch, bio


def optimize_scanorama(
    input_path: str,
    batch_key: str,
    label_key: str,
    n_trials: int,
    output_path: str,
    smoketest: bool = False,
):
    if smoketest:
        n_trials = 2
    adata = io.to_anndata(input_path)

    # PCA reduction for HPO speed — scanorama on 50 PCs instead of 1040 features.
    from sklearn.decomposition import PCA
    n_pcs = min(_HPO_N_PCS, adata.X.shape[1])
    logger.info("Computing PCA (%d components) for HPO...", n_pcs)
    X = adata.X if isinstance(adata.X, np.ndarray) else adata.X.toarray()
    pca = PCA(n_components=n_pcs, random_state=42)
    X_pca = pca.fit_transform(X)
    logger.info("PCA done. Variance explained: %.1f%%", 100 * pca.explained_variance_ratio_.sum())

    adata_pca = ad.AnnData(X=X_pca, obs=adata.obs.copy())

    from utils import warmup_benchmark
    warmup_benchmark(batch_key, label_key)

    study = optuna.create_study(directions=["maximize", "maximize"], sampler=optuna.samplers.TPESampler(seed=42))
    study.optimize(
        lambda trial: objective(trial, adata_pca, batch_key, label_key, smoketest),
        n_trials=n_trials,
    )

    save_optuna_results(study, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Use Optuna to tune hyperparameters for Scanorama."
    )
    parser.add_argument("--input_data", required=True, help="Path to input data.")
    parser.add_argument("--batch_key", required=True, help="Batch key.")
    parser.add_argument("--label_key", required=True, help="Label key.")
    parser.add_argument("--n_trials", required=True, help="How many trials to run.")
    parser.add_argument("--output_path", required=True, help="Where to save results.")
    parser.add_argument(
        "--smoketest", action="store_true", help="Run a smoketest with limited trials"
    )
    args = parser.parse_args()

    optimize_scanorama(
        input_path=args.input_data,
        batch_key=args.batch_key,
        label_key=args.label_key,
        n_trials=int(args.n_trials),
        output_path=args.output_path,
        smoketest=args.smoketest,
    )
