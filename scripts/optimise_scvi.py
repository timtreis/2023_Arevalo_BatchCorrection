import os
import sys
import logging
import argparse
from preprocessing import io
import pandas as pd
import anndata as ad
import optuna
import scvi

from utils import scib_benchmark_embedding, save_optuna_results, _stratified_subsample

logger = logging.getLogger(__name__)

def objective(
    trial,
    adata: ad.AnnData,
    batch_key: str,
    label_key: str,
    multiple_covariates: bool = False,
    smoketest: bool = False,
    gene_likelihood: str = "zinb",
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

    # For ZINB/NB likelihoods, shift data to non-negative (required by count models).
    # For normal likelihood, no shift needed (Gaussian handles negative values).
    if gene_likelihood != "normal":
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

    # labels_key is intentionally omitted: scVI is unsupervised and does not
    # use labels during training.  Registering a high-cardinality label column
    # (e.g. 82 K compounds) wastes memory without any benefit.
    if multiple_covariates:
        scvi.model.SCVI.setup_anndata(
            adata,
            batch_key=actual_batch_key,
            categorical_covariate_keys=categorical_covariate_keys,
        )
    else:
        scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key)

    vae = scvi.model.SCVI(
        adata,
        n_hidden=n_hidden,
        n_latent=n_latent,
        n_layers=n_layers,
        dropout_rate=dropout_rate,
        gene_likelihood=gene_likelihood,
    )

    # Gradient clipping stabilizes training with gene_likelihood="normal",
    # which can produce NaN gradients on certain architectures (deep + wide).
    # ZINB is inherently more stable due to the log-link function.
    train_kwargs = dict(
        max_epochs=n_epochs,
        early_stopping=True,
        early_stopping_monitor="validation_loss",
    )
    if gene_likelihood == "normal":
        train_kwargs["gradient_clip_val"] = 1.0

    try:
        vae.train(**train_kwargs)
        vals = vae.get_latent_representation()
    except (ValueError, RuntimeError) as e:
        # NaN in encoder output or CUDA errors — skip this trial.
        sys.stdout.close()
        sys.stdout = sys.__stdout__
        print(f"Trial failed during training: {e}")
        return None, None

    prefix = "scvi_normal" if gene_likelihood == "normal" else "scvi"
    features = [f"{prefix}_{i}" for i in range(vals.shape[1])]
    integrated_adata = ad.AnnData(
        X=pd.DataFrame(vals, columns=features, index=adata.obs_names),
        obs=adata.obs.copy()
    )

    batch_score, bio_score = scib_benchmark_embedding(
        adata=integrated_adata,
        batch_key=actual_batch_key,
        label_key=label_key,
        lightweight=True,
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
    gene_likelihood: str = "zinb",
):
    if smoketest:
        n_trials = 2
    adata = io.to_anndata(input_path)
    print(multiple_covariates)

    # Subsample for HPO — relative rankings are preserved; correction uses full data.
    HPO_MAX_CELLS = 100_000
    if adata.n_obs > HPO_MAX_CELLS:
        _key = batch_key if isinstance(batch_key, str) else batch_key[0]
        adata = _stratified_subsample(adata, [_key], HPO_MAX_CELLS)
        print(f"HPO subsample: {adata.n_obs} cells")

    from utils import warmup_benchmark
    warmup_benchmark(batch_key, label_key)

    study = optuna.create_study(directions=["maximize", "maximize"], sampler=optuna.samplers.TPESampler(seed=42))
    study.optimize(lambda trial: objective(trial, adata.copy(), batch_key, label_key, multiple_covariates, smoketest, gene_likelihood), n_trials=n_trials)

    save_optuna_results(study, output_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Use Optuna to tune hyperparameters for scPoli.")

    parser.add_argument("--input_data", required=True, help="Path to input data.")
    parser.add_argument("--batch_key", required=True, help="Batch key.")
    parser.add_argument("--label_key", required=True, help="Label key.")
    parser.add_argument("--multi", action="store_true", help="Use one or multiple categorical covariates")
    parser.add_argument("--n_trials", required=True, help="How many trials to run.")
    parser.add_argument("--output_path", required=True, help="Where to save the optimal parameter set.")
    parser.add_argument("--smoketest", action="store_true", help="Run a smoketest with limited epochs")
    parser.add_argument("--gene_likelihood", default="zinb", choices=["zinb", "nb", "normal"],
                        help="Observation model likelihood (default: zinb). Use 'normal' for continuous features.")

    args = parser.parse_args()

    optimize_scvi(
        input_path=args.input_data,
        batch_key=args.batch_key,
        label_key=args.label_key,
        multiple_covariates=args.multi,
        n_trials=int(args.n_trials),
        output_path=args.output_path,
        smoketest=args.smoketest,
        gene_likelihood=args.gene_likelihood,
    )
