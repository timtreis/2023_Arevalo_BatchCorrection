import os
import sys
import logging
import argparse
from preprocessing import io
import pandas as pd
import anndata as ad
import optuna
import scvi
import torch.distributions as dist
dist.Distribution.set_default_validate_args(False)    # disable global validation

from utils import scib_benchmark_embedding

logger = logging.getLogger(__name__)

def objective(
    trial,
    adata: ad.AnnData,
    batch_key: str,
    label_key: str,
    params: dict,
    multiple_covariates: bool = False,
    smoketest: bool = False,
):

    # silence stdout
    sys.stdout = open(os.devnull, "w")

    # ----------  fixed scVI architecture ----------
    n_hidden = params["params_n_hidden"]
    n_latent = params["params_n_latent"]
    n_layers = params["params_n_layers"]
    dropout_rate = params["params_dropout_rate"]

    # ----------  scANVI hyper-parameters to tune ----------
    classification_ratio = trial.suggest_int("classification_ratio", 20, 80)
    n_epochs_kl_warmup = trial.suggest_int("n_epochs_kl_warmup", 5, 50)
    linear_classifier = trial.suggest_categorical("linear_classifier", [True, False])

    n_epochs = 2 if smoketest else 50

    # make sure counts stay ≥0
    adata.X -= adata.X.min()

    # -------------  anndata setup -------------
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
    # -------------  pre-train SCVI -------------
    vae = scvi.model.SCVI(
        adata,
        n_hidden=n_hidden,
        n_latent=n_latent,
        n_layers=n_layers,
        dropout_rate=dropout_rate,
    )
    vae.train(max_epochs=n_epochs, early_stopping=True, early_stopping_monitor="elbo_validation")

    # -------------  transfer to scANVI -------------
    scanvi = scvi.model.SCANVI.from_scvi_model(
        vae, 
        unlabeled_category="Unknown",
        linear_classifier=linear_classifier,
    )

    plan_kwargs = {
        "classification_ratio": classification_ratio,
        "n_epochs_kl_warmup":   n_epochs_kl_warmup,
    }

    scanvi.train(
        max_epochs=n_epochs,
        early_stopping=True,
        early_stopping_monitor="elbo_validation",
        plan_kwargs=plan_kwargs,
    )

    # -------------  evaluation -------------
    vals = scanvi.get_latent_representation()
    features = [f"scanvi_{i}" for i in range(vals.shape[1])]
    integrated = ad.AnnData(
        X=pd.DataFrame(vals, columns=features, index=adata.obs_names),
        obs=adata.obs.copy(),
    )

    batch_score, bio_score = scib_benchmark_embedding(
        adata=integrated,
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
    scvi_params_path: str,
    n_trials: int,
    output_path: str,
    multiple_covariates: bool = False,
    smoketest: bool = False,
):
    if smoketest:
        n_trials = 2

    # load scVI parameters
    params = pd.read_csv(scvi_params_path)
    params = params.sort_values("total", ascending=False).iloc[0].to_dict()

    if params["state"] != "COMPLETE":
        raise ValueError("Optimization did not complete successfully")

    print("\nUsing pre-optimized scVI hyperparameters:\n")
    for k, v in params.items():
        if "params_" in k:
            print(f"- {k}: {v}")

    adata = io.to_anndata(input_path)

    study = optuna.create_study(directions=["maximize", "maximize"])
    study.optimize(lambda trial: objective(trial, adata.copy(), batch_key, label_key, params, multiple_covariates, smoketest), n_trials=n_trials)

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
    parser.add_argument("--scvi_params_path", required=True, help="Optimal hyperparameters for the scVI model.")
    parser.add_argument("--n_trials", required=True, help="How many trials to run.")
    parser.add_argument("--output_path", required=True, help="Where to save the optimal parameter set.")
    parser.add_argument("--smoketest", action="store_true", help="Run a smoketest with limited epochs")

    args = parser.parse_args()

    optimize_scvi(
        input_path=args.input_data,
        batch_key=args.batch_key,
        label_key=args.label_key,
        multiple_covariates=args.multi,
        scvi_params_path=args.scvi_params_path,
        n_trials=int(args.n_trials),
        output_path=args.output_path,
        smoketest=args.smoketest,
    )
