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
from utils import scib_benchmark_embedding

logger = logging.getLogger(__name__)

def objective(
    trial,
    adata: ad.AnnData, 
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

    # Ensure the data matrix is a dense numpy array
    if not isinstance(adata.X, np.ndarray):
        data_mat = adata.X.toarray()
    else:
        data_mat = adata.X
    meta_data = adata.obs.copy()

    ho = run_harmony(
        data_mat=data_mat,
        meta_data=meta_data,
        vars_use=batch_key,
        theta=theta,
        lamb=lamb,
        sigma=sigma,
        nclust=nclust,
        tau=tau,
        random_state=0
    )

    # Retrieve the integrated data (Z_corr has shape [d, N]) and transpose to [N, d]
    vals = ho.result().T
    features = [f"harmony_{i}" for i in range(vals.shape[1])]
    integrated_adata = ad.AnnData(
        X=pd.DataFrame(vals, columns=features, index=adata.obs_names),
        obs=adata.obs.copy()
    )

    # Evaluate the integration using scIB metrics
    batch, bio = scib_benchmark_embedding(
        adata=integrated_adata,
        batch_key=batch_key, 
        label_key=label_key
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

    study = optuna.create_study(directions=["maximize", "maximize"])
    study.optimize(lambda trial: objective(trial, adata.copy(), batch_key, label_key, smoketest), n_trials=n_trials)

    df = study.trials_dataframe()
    df = df.rename(columns={"values_0": "batch", "values_1": "bio"})
    df["total"] = 0.6 * df["bio"] + 0.4 * df["batch"]
    df = df.sort_values("total", ascending=False)
    df.to_csv(output_path, index=False)

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
