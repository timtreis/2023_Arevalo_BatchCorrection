import argparse
import logging
from scvi.external import SysVI
from preprocessing import io
import pandas as pd

logger = logging.getLogger(__name__)


def correct_with_sysvi(
    dframe_path: str,
    batch_key: list[str] | str,
    label_key: str,
    parameter_path: str,
    output_path: str,
    smoketest: bool = False,
):
    """sysVI correction from https://docs.scvi-tools.org/en/1.4.1/user_guide/models/sysvi.html"""

    # load hyperparameters
    params = pd.read_csv(parameter_path)
    params = params.sort_values("total", ascending=False).iloc[0].to_dict()

    if params["state"] != "COMPLETE":
        raise ValueError("Optimization did not complete successfully")

    dropout_rate = params["params_dropout_rate"]
    kl_weight = params["params_kl_weight"]
    n_hidden = params["params_n_hidden"]
    n_latent = params["params_n_latent"]
    n_layers = params["params_n_layers"]
    n_prior_components = params["params_n_prior_components"]
    prior = params["params_prior"]
    use_z_distance_cycle_weight = params["params_use_z_distance_cycle_weight"]
    z_distance_cycle_weight = params["params_z_distance_cycle_weight"]

    n_epochs = 2 if smoketest else 999999

    print("\nUsing the following hyperparameters:")
    print(f"- dropout_rate: {dropout_rate}")
    print(f"- kl_weight: {kl_weight}")
    print(f"- n_hidden: {n_hidden}")
    print(f"- n_latent: {n_latent}")
    print(f"- n_layers: {n_layers}")
    print(f"- n_prior_components: {n_prior_components}")
    print(f"- prior: {prior}")
    print(f"- use_z_distance_cycle_weight: {use_z_distance_cycle_weight}")
    print(f"- z_distance_cycle_weight: {z_distance_cycle_weight}")
    print(f"- n_epochs: {n_epochs}\n")

    adata = io.to_anndata(dframe_path)
    meta = adata.obs.reset_index(drop=True).copy()

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

    SysVI.setup_anndata(
        adata,
        batch_key=actual_batch_key,
        categorical_covariate_keys=categorical_covariate_keys,
    )
    vae = SysVI(
        adata,
        n_hidden=n_hidden,
        n_latent=n_latent,
        n_layers=n_layers,
        dropout_rate=dropout_rate,
        prior=prior,
        n_prior_components=n_prior_components,
    )
    vae.view_anndata_setup(adata=adata)
    vae.train(
        max_epochs=n_epochs,
        early_stopping=True,
        early_stopping_monitor="validation_loss",
        plan_kwargs={
            "kl_weight": kl_weight,
            "z_distance_cycle_weight": z_distance_cycle_weight if use_z_distance_cycle_weight else 0,
        }
    )

    vals = vae.get_latent_representation(adata=adata)
    features = [f"sysvi_{i}" for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform sysVI correction on data.")

    parser.add_argument("--input_data", required=True, help="Path to input data")
    parser.add_argument(
        "--batch_key",
        required=True,
        nargs="+",
        help="Batch key(s); provide multiple keys separated by spaces or as a comma-separated string.",
    )
    parser.add_argument("--label_key", required=True, help="Label key")
    parser.add_argument("--parameter_path", required=True, help="Path to the parameter file")
    parser.add_argument("--output_path", required=True, help="Path to save corrected data")
    parser.add_argument(
        "--smoketest",
        action="store_true",
        help="Run a smoketest with limited epochs (sets max_epochs to 2)",
    )

    args = parser.parse_args()

    correct_with_sysvi(
        dframe_path=args.input_data,
        batch_key=args.batch_key,
        label_key=args.label_key,
        parameter_path=args.parameter_path,
        output_path=args.output_path,
        smoketest=args.smoketest,
    )
