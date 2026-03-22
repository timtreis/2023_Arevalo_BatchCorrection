import os
import sys
import logging
import argparse
import tempfile
import pandas as pd
import anndata as ad
import optuna
from desc import train, scale_bygroup

from preprocessing import io
from utils import scib_benchmark_embedding

logger = logging.getLogger(__name__)


def objective(
    trial,
    adata: ad.AnnData,
    batch_key: str,
    label_key: str,
    smoketest: bool = False,
):
    sys.stdout = open(os.devnull, "w")

    n_neighbors = trial.suggest_int("n_neighbors", 10, 100)
    louvain_resolution = trial.suggest_float("louvain_resolution", 0.2, 2.0)
    hidden_dim = trial.suggest_int("hidden_dim", 32, 256, step=32)
    tol = trial.suggest_float("tol", 0.001, 0.1, log=True)

    dims = [adata.shape[1], hidden_dim, 32]

    scale_bygroup(adata, batch_key, max_value=None)

    adata = train(
        adata,
        dims=dims,
        tol=tol,
        n_neighbors=n_neighbors,
        batch_size=128 if smoketest else 1024,
        louvain_resolution=louvain_resolution,
        save_encoder_weights=False,
        save_dir=tempfile.TemporaryDirectory().name,
        do_tsne=False,
        use_GPU=True,
        GPU_id=0,
        num_Cores=1,
        use_ae_weights=False,
        do_umap=False,
    )

    vals = adata.obsm[f"X_Embeded_z{louvain_resolution}"]
    features = [f"desc_{i}" for i in range(vals.shape[1])]
    integrated_adata = ad.AnnData(
        X=pd.DataFrame(vals, columns=features, index=adata.obs_names),
        obs=adata.obs.copy(),
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


def optimize_desc(
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
    study.optimize(
        lambda trial: objective(trial, adata.copy(), batch_key, label_key, smoketest),
        n_trials=n_trials,
    )

    df = study.trials_dataframe()
    df = df.rename(columns={"values_0": "batch", "values_1": "bio"})
    df["total"] = 0.6 * df["bio"] + 0.4 * df["batch"]
    df = df.sort_values("total", ascending=False)
    df.to_csv(output_path, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Use Optuna to tune hyperparameters for DESC."
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

    optimize_desc(
        input_path=args.input_data,
        batch_key=args.batch_key,
        label_key=args.label_key,
        n_trials=int(args.n_trials),
        output_path=args.output_path,
        smoketest=args.smoketest,
    )
