import logging
import argparse
from pathlib import Path
from typing import List, Union, Optional, Literal
from scarches.models.scpoli import scPoli
from preprocessing import io
from utils import coarsen_labels
import scanpy as sc
import anndata as ad
import pandas as pd

logger = logging.getLogger(__name__)


def correct_with_scpoli(
    dframe_path: str,
    batch_key: Union[List[str], str],  # we're using Python 3.9
    label_key: str,
    parameter_path: str,
    output_path: str,
    preproc: Optional[Literal["pca"]] = None,
    smoketest: bool = False,
    model_dir: Optional[str] = None,
    **kwargs,
):
    """scPoli correction from https://www.nature.com/articles/s41592-023-02035-2"""

    # load hyperparameters
    params = pd.read_csv(parameter_path)
    params = params.sort_values("total", ascending=False).iloc[0].to_dict()

    if params["state"] != "COMPLETE":
        raise ValueError("Optimization did not complete successfully")

    alpha_epoch_anneal = params["params_alpha_epoch_anneal"]
    embedding_dims = params["params_embedding_dims"]
    eta = params["params_eta"]
    latent_dim = params["params_latent_dim"]
    hidden_layer_sizes = [
        int(params[f"params_layer_{i}_size"]) for i in range(params["params_num_layers"])
    ]
    pretrain_to_train_ratio = params["params_pretrain_to_train_ratio"]

    # scPoli disables early stopping during pretraining, so pretraining epochs
    # must be finite. Use the HPO-tuned ratio to split a fixed budget for
    # pretraining, then train until early stopping.
    total_pretrain_budget = 4 if smoketest else 400
    n_pretrain_epochs = 2 if smoketest else int(total_pretrain_budget * pretrain_to_train_ratio)
    # Cap training epochs at 600.  scPoli's EarlyStopping with reduce_lr=True
    # can extend effective patience to 500+ epochs past the minimum on large
    # datasets, wasting hours of compute.  Observed best training epochs:
    # S1=84, S2=209, S3=188, S4=101.  600 gives ~3x headroom over the worst
    # case while preventing the runaway (S5 reached 700+ with 999999 cap).
    n_train_epochs = 2 if smoketest else 600

    print("\nUsing the following hyperparameters:")
    print(f"- alpha_epoch_anneal: {alpha_epoch_anneal}")
    print(f"- embedding_dims: {embedding_dims}")
    print(f"- eta: {eta}")
    print(f"- latent_dim: {latent_dim}")
    print(f"- hidden_layer_sizes: {hidden_layer_sizes}")
    print(f"- pretrain_to_train_ratio: {pretrain_to_train_ratio}")
    print(f"- n_pretrain_epochs: {n_pretrain_epochs}")
    print(f"- n_train_epochs: {n_train_epochs} (early stopping)\n")

    if isinstance(batch_key, list) and len(batch_key) == 1 and "," in batch_key[0]:
        batch_key = batch_key[0].split(",")

    if isinstance(batch_key, str):
        batch_key = batch_key.split(",") if "," in batch_key else [batch_key]

    adata = io.to_anndata(dframe_path)
    meta = adata.obs.reset_index(drop=True).copy()

    # Mark rare compounds as unlabeled for semi-supervised training
    coarsen_labels(adata, label_key, batch_key)

    if preproc == "pca":
        logger.info("Applying PCA preprocessing with Scanpy")
        sc.pp.pca(adata, svd_solver="arpack", n_comps=50)
        # Create new adata object with PCA-transformed data
        adata_for_training = ad.AnnData(adata.obsm["X_pca"].copy())
        adata_for_training.obs = adata.obs.copy()
        adata = adata_for_training
        logger.info("PCA completed. Shape of PCA-transformed data: %s", adata.X.shape)


    model = scPoli(
        adata=adata,
        condition_keys=batch_key,
        cell_type_keys=label_key,
        hidden_layer_sizes=hidden_layer_sizes,
        latent_dim=latent_dim,
        embedding_dims=embedding_dims,
        recon_loss="mse",
    )

    model.train(
        n_epochs=n_train_epochs,
        pretraining_epochs=n_pretrain_epochs,
        use_early_stopping=True,
        reload_best=True,
        alpha_epoch_anneal=alpha_epoch_anneal,
        eta=eta,
    )

    model.model.eval()

    if model_dir is not None:
        Path(model_dir).mkdir(parents=True, exist_ok=True)
        model.save(model_dir, overwrite=True)

    vals = model.get_latent(adata, mean=True)
    features = [f"scpoli_{i}" for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="scPoli correction from https://www.nature.com/articles/s41592-023-02035-2"
    )

    parser.add_argument("--input_data", required=True, help="Path to the input data")
    parser.add_argument(
        "--batch_key",
        nargs="+",
        required=True,
        help="Batch key(s), provide multiple keys separated by commas",
    )
    parser.add_argument("--label_key", required=True, help="Label key")
    parser.add_argument("--parameter_path", required=True, help="Path to the parameter file")
    parser.add_argument(
        "--output_path", required=True, help="Path to save the output data"
    )
    parser.add_argument(
        "--preproc",
        type=str,
        choices=["pca"],
        default=None,
        help="Preprocessing method to apply. Choices: 'pca' or None",
    )
    parser.add_argument(
        "--smoketest", action="store_true", help="Run a smoketest with limited epochs"
    )
    parser.add_argument(
        "--model_dir", default=None,
        help="If set, save the trained scPoli model to this directory for downstream reference mapping.",
    )

    args = parser.parse_args()

    correct_with_scpoli(
        dframe_path=args.input_data,
        batch_key=args.batch_key,
        label_key=args.label_key,
        parameter_path=args.parameter_path,
        output_path=args.output_path,
        preproc=args.preproc,
        smoketest=args.smoketest,
        model_dir=args.model_dir,
    )
