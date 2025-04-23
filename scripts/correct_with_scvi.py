import logging
import argparse
import scvi
import pandas as pd
from preprocessing import io

logger = logging.getLogger(__name__)


def correct_with_scvi(
    dframe_path: str,
    batch_key: str | list,
    label_key: str,
    parameter_path: str,
    output_path: str,
    multiple_covariates: bool = False,
    smoketest: bool = False,
):
    """scVI correction"""

    # load hyperparameters
    params = pd.read_csv(parameter_path)
    params = params.sort_values("total", ascending=False).iloc[0].to_dict()

    if params["state"] != "COMPLETE":
        raise ValueError("Optimization did not complete successfully")

    dropout_rate = params["params_dropout_rate"]
    n_hidden = params["params_n_hidden"]
    n_latent = params["params_n_latent"]
    n_layers = params["params_n_layers"]

    n_epochs = 2 if smoketest else 999999

    print("\nUsing the following hyperparameters:")
    print(f"- dropout_rate: {dropout_rate}")
    print(f"- n_hidden: {n_hidden}")
    print(f"- n_latent: {n_latent}")
    print(f"- n_layers: {n_layers}")
    print(f"- n_epochs: {n_epochs}\n")

    adata = io.to_anndata(dframe_path)
    meta = adata.obs.reset_index(drop=True).copy()

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
    vae.view_anndata_setup(adata=adata)
    vae.train(
        max_epochs=n_epochs,
        early_stopping=True,
        early_stopping_monitor="elbo_validation",
    )

    vals = vae.get_latent_representation()
    features = [f"scvi_{i}" for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform scVI correction on data.")
    parser.add_argument("--input_data", required=True, help="Path to input data")
    parser.add_argument("--batch_key", required=True, help="Batch key")
    parser.add_argument("--label_key", required=True, help="Label key")
    parser.add_argument("--multi", action="store_true", help="Use one or multiple categorical covariates")
    parser.add_argument("--parameter_path", required=True, help="Path to the parameter file")
    parser.add_argument("--output_path", required=True, help="Path to save corrected data")
    parser.add_argument("--smoketest", action="store_true", help="Run a smoketest with limited epochs")

    args = parser.parse_args()

    correct_with_scvi(
        dframe_path=args.input_data,
        batch_key=args.batch_key,
        label_key=args.label_key,
        multiple_covariates=args.multi,
        parameter_path=args.parameter_path,
        output_path=args.output_path,
        smoketest=args.smoketest,
    )
