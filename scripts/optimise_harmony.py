import os
import sys
import logging
import argparse
import pandas as pd
import anndata as ad
import optuna
import numpy as np

from preprocessing import io
from harmonypy import run_harmony
from utils import scib_benchmark_embedding, save_optuna_results, _stratified_subsample

logger = logging.getLogger(__name__)

# Number of PCs for HPO trials. Reduces 1040 features to 50 dims, making
# harmony correction ~20× faster. The final correction uses full features.
_HPO_N_PCS = 50


def objective(
    trial,
    pca_mat: np.ndarray,
    obs: pd.DataFrame,
    batch_key: str,
    label_key: str,
    smoketest: bool = False,
):
    # Silence output
    sys.stdout = open(os.devnull, "w")

    # Optimize Harmony-specific hyperparameters
    sigma  = trial.suggest_float("sigma", 0.05, 1.0)
    theta  = trial.suggest_float("theta", 0.1, 5.0)
    lamb   = trial.suggest_float("lamb", 0.1, 5.0)
    nclust = trial.suggest_int("nclust", 2, 500)
    tau    = trial.suggest_float("tau", 0.0, 1.0)

    # Use GPU if available (harmonypy >= 0.2 supports device parameter)
    import inspect
    harmony_kwargs = dict(
        data_mat=pca_mat,
        meta_data=obs,
        vars_use=batch_key,
        theta=theta,
        lamb=lamb,
        sigma=sigma,
        nclust=nclust,
        tau=tau,
        random_state=0,
    )
    if "device" in inspect.signature(run_harmony).parameters:
        import torch
        harmony_kwargs["device"] = "cuda" if torch.cuda.is_available() else "cpu"

    ho = run_harmony(**harmony_kwargs)

    # harmonypy v1 returns (d, N), v2 returns (N, d)
    vals = ho.result()
    if vals.shape[0] != len(obs):
        vals = vals.T
    features = [f"harmony_{i}" for i in range(vals.shape[1])]
    integrated_adata = ad.AnnData(
        X=pd.DataFrame(vals, columns=features, index=obs.index),
        obs=obs.copy()
    )

    # Evaluate the integration using scIB metrics
    batch, bio = scib_benchmark_embedding(
        adata=integrated_adata,
        batch_key=batch_key,
        label_key=label_key,
        lightweight=True,
    )

    # Restore output
    sys.stdout.close()
    sys.stdout = sys.__stdout__

    return batch, bio

def optimize_harmony(
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

    # Subsample for HPO — relative rankings are preserved; correction uses full data.
    HPO_MAX_CELLS = 100_000
    if adata.n_obs > HPO_MAX_CELLS:
        _key = batch_key if isinstance(batch_key, str) else batch_key[0]
        adata = _stratified_subsample(adata, [_key], HPO_MAX_CELLS)
        print(f"HPO subsample: {adata.n_obs} cells")

    # PCA reduction for HPO speed — harmony on 50 PCs instead of 1040 features.
    # The tuned hyperparameters transfer to the full-feature correction.
    from sklearn.decomposition import PCA
    n_pcs = min(_HPO_N_PCS, adata.X.shape[1])
    logger.info("Computing PCA (%d components) for HPO...", n_pcs)
    X = adata.X if isinstance(adata.X, np.ndarray) else adata.X.toarray()
    pca = PCA(n_components=n_pcs, random_state=42)
    pca_mat = pca.fit_transform(X)
    obs = adata.obs.copy()
    logger.info("PCA done. Variance explained: %.1f%%", 100 * pca.explained_variance_ratio_.sum())

    # Warm up JAX JIT cache before trials so compilation cost is paid once
    from utils import warmup_benchmark
    warmup_benchmark(batch_key, label_key)

    study = optuna.create_study(directions=["maximize", "maximize"], sampler=optuna.samplers.TPESampler(seed=42))
    study.optimize(lambda trial: objective(trial, pca_mat, obs, batch_key, label_key, smoketest), n_trials=n_trials)

    save_optuna_results(study, output_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Use Optuna to tune hyperparameters for Harmony integration.")
    parser.add_argument("--input_data", required=True, help="Path to input data.")
    parser.add_argument("--batch_key", required=True, help="Batch key.")
    parser.add_argument("--label_key", required=True, help="Label key.")
    parser.add_argument("--n_trials", required=True, help="How many trials to run.")
    parser.add_argument("--output_path", required=True, help="Where to save the optimal parameter set.")
    parser.add_argument("--smoketest", action="store_true", help="Run a smoketest with limited epochs")
    args = parser.parse_args()

    optimize_harmony(
        input_path=args.input_data,
        batch_key=args.batch_key,
        label_key=args.label_key,
        n_trials=int(args.n_trials),
        output_path=args.output_path,
        smoketest=args.smoketest,
    )
