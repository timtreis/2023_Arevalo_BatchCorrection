import os
import sys
import logging
import argparse
from preprocessing import io
import pandas as pd
import anndata as ad
import optuna
import scvi

from utils import scib_benchmark_embedding

logger = logging.getLogger(__name__)

def objective(
    trial,
    adata: ad.AnnData, 
    batch_key: str, 
    label_key: str, 
    multiple_covariates: bool = False,
    smoketest: bool = False,
):

    # Silence output during training and evaluation
    sys.stdout = open(os.devnull, "w")

    # Optimize scVI hyperparameters:
    n_hidden = trial.suggest_int("n_hidden", 64, 256, step=64)
    n_latent = trial.suggest_int("n_latent", 10, 100)
    n_layers = trial.suggest_int("n_layers", 1, 3)
    dropout_rate = trial.suggest_float("dropout_rate", 0.0, 0.5)
    # learning_rate = trial.suggest_float("learning_rate", 1e-4, 1e-2, log=True)
    n_epochs = 2 if smoketest else 50

    # Preprocess data: subtract the minimum to ensure non-negative values (if needed)
    min_value = adata.X.min()
    adata.X -= min_value

    if multiple_covariates:
        if isinstance(batch_key, list) and len(batch_key) == 1:
            batch_key = batch_key[0]

        batch_key = batch_key.split(",")

        if isinstance(batch_key, list):
            actual_batch_key = batch_key[0]
            assert isinstance(actual_batch_key, str)

            categorical_covariate_keys = batch_key[1:]
            if isinstance(categorical_covariate_keys, str):
                categorical_covariate_keys = [categorical_covariate_keys]
        else:
            actual_batch_key = batch_key
            categorical_covariate_keys = [None]
    else:
        actual_batch_key = batch_key

    if multiple_covariates:
        scvi.model.SCVI.setup_anndata(
            adata,
            batch_key=actual_batch_key,
            categorical_covariate_keys=categorical_covariate_keys,
            labels_key=label_key,
        )
    else:
        scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key, labels_key=label_key)

    vae = scvi.model.SCVI(
        adata,
        n_hidden=n_hidden,
        n_latent=n_latent,
        n_layers=n_layers,
        dropout_rate=dropout_rate,
    )

    vae.train(
        max_epochs=n_epochs,
        early_stopping=True,
        early_stopping_monitor="elbo_validation",
    )

    vals = vae.get_latent_representation()
    features = [f"scvi_{i}" for i in range(vals.shape[1])]
    integrated_adata = ad.AnnData(
        X=pd.DataFrame(vals, columns=features, index=adata.obs_names),
        obs=adata.obs.copy()
    )

    batch_score, bio_score = scib_benchmark_embedding(
        adata=integrated_adata,
        batch_key=actual_batch_key,
        label_key=label_key,
    )

    sys.stdout.close()
    sys.stdout = sys.__stdout__

    return batch_score, bio_score

def optimize_scvi(
    input_path: str,
    batch_key: str,
    label_key: str,
    n_trials: int,
    output_path: str,
    multiple_covariates: bool = False,
    smoketest: bool = False,
):
    if smoketest:
        n_trials = 2
    adata = io.to_anndata(input_path)
    print(multiple_covariates)
    study = optuna.create_study(directions=["maximize", "maximize"])
    study.optimize(lambda trial: objective(trial, adata.copy(), batch_key, label_key, multiple_covariates, smoketest), n_trials=n_trials)

    df = study.trials_dataframe()
    df = df.rename(columns={"values_0": "batch", "values_1": "bio"})
    df["total"] = 0.6 * df["bio"] + 0.4 * df["batch"]
    df = df.sort_values("total", ascending=False)
    df.to_csv(output_path, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Use Optuna to tune hyperparameters for scPoli.")

    parser.add_argument("--input_data", required=True, help="Path to input data.")
    parser.add_argument("--batch_key", required=True, help="Batch key.")
    parser.add_argument("--label_key", required=True, help="Label key.")
    parser.add_argument("--multi", action="store_true", help="Use one or multiple categorical covariates")
    parser.add_argument("--n_trials", required=True, help="How many trials to run.")
    parser.add_argument("--output_path", required=True, help="Where to save the optimal parameter set.")
    parser.add_argument("--smoketest", action="store_true", help="Run a smoketest with limited epochs")

    args = parser.parse_args()

    optimize_scvi(
        input_path=args.input_data,
        batch_key=args.batch_key,
        label_key=args.label_key,
        multiple_covariates=args.multi,
        n_trials=int(args.n_trials),
        output_path=args.output_path,
        smoketest=args.smoketest,
    )
