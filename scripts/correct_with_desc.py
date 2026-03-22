import logging
import argparse
from desc import train, scale_bygroup
from preprocessing import io
import tempfile

logger = logging.getLogger(__name__)


def _load_hyperparameters(parameter_path):
    import pandas as pd
    params = pd.read_csv(parameter_path)
    params = params.sort_values("total", ascending=False).iloc[0].to_dict()
    if params["state"] != "COMPLETE":
        raise ValueError("Optimization did not complete successfully")
    n_neighbors = int(params["params_n_neighbors"])
    louvain_resolution = params["params_louvain_resolution"]
    hidden_dim = int(params["params_hidden_dim"])
    tol = params["params_tol"]
    print(f"\nUsing tuned params: n_neighbors={n_neighbors}, resolution={louvain_resolution}, hidden_dim={hidden_dim}, tol={tol}\n")
    return n_neighbors, louvain_resolution, hidden_dim, tol


def correct_with_desc(
    dframe_path: str, batch_key: str, output_path: str, parameter_path=None, smoketest=False
):
    """DESC correction"""

    adata = io.to_anndata(dframe_path)
    meta = adata.obs.reset_index(drop=True).copy()

    res = 1.0
    hidden_dim = 128
    tol = 0.05
    n_neighbors = 64

    if parameter_path is not None:
        n_neighbors, res, hidden_dim, tol = _load_hyperparameters(parameter_path)

    scale_bygroup(adata, batch_key, max_value=None)

    # Adjust parameters based on smoketest flag
    batch_size = 1024
    num_Cores = 100  # it'll grab what it can, snakemake scales it down
    if smoketest:
        batch_size = 128
        num_Cores = 1
        n_neighbors = min(n_neighbors, 10)

    adata = train(
        adata,
        dims=[adata.shape[1], hidden_dim, 32],
        tol=tol,
        n_neighbors=n_neighbors,
        batch_size=batch_size,
        louvain_resolution=res,
        save_encoder_weights=False,
        save_dir=tempfile.TemporaryDirectory().name,
        do_tsne=False,
        use_GPU=True,
        GPU_id=0,
        num_Cores=num_Cores,
        use_ae_weights=False,
        do_umap=False,
    )

    vals = adata.obsm[f"X_Embeded_z{res}"]

    features = [f"desc_{i}" for i in range(vals.shape[1])]
    io.merge_parquet(meta, vals, features, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform DESC correction on data.")

    parser.add_argument(
        "--input_data",
        required=True,
        help="Path to the input data file (Parquet format).",
    )
    parser.add_argument(
        "--batch_key",
        required=True,
        help="Batch key(s), separated by commas if multiple.",
    )
    parser.add_argument("--parameter_path", default=None, help="Path to Optuna parameter CSV.")
    parser.add_argument(
        "--output_path", required=True, help="Path to save the corrected output data."
    )
    parser.add_argument(
        "--smoketest",
        action="store_true",
        help="Run a smoketest with reduced parameters.",
    )

    args = parser.parse_args()

    correct_with_desc(
        dframe_path=args.input_data,
        batch_key=args.batch_key,
        output_path=args.output_path,
        parameter_path=args.parameter_path,
        smoketest=args.smoketest,
    )
