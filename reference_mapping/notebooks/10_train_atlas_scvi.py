# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: RefMap scVI (GPU)
#     language: python
#     name: refmap-scvi
# ---

# %% [markdown]
# # 10 — scVI atlas
#
# Two paths:
#
# 1. **Pipeline path (preferred)**: re-run the patched
#    `methods_scvi_single` rule against `scenario_<refscenario>` so the model
#    artifact lands at `outputs/<scenario>/mad_int_featselect_scvi_single_model/`.
#    This notebook then symlinks it into `models/`.
#
# 2. **Notebook path (fallback)**: re-train inline using the HPO best params,
#    skipping the Snakemake plumbing.
#
# Path 1 is faster per re-run and matches the figures shown in the paper, so
# default to it. Path 2 is for when the pipeline isn't conveniently runnable.

# %%
from pathlib import Path
import sys

sys.path.insert(0, str(Path.cwd().parent))
from src.paths import MODEL_OUT, scenario_model_dir, scenario_input_parquet, scenario_optuna_csv

# %% [markdown]
# ## Parameters

# %%
QUERY_SOURCE = "source_5"
REFERENCE_SCENARIO = "scenario_5"  # the scenario whose pipeline output we use
ATLAS_NAME = f"scvi_atlas_excl_{QUERY_SOURCE}"

# %% [markdown]
# ## Path 1 — symlink the pipeline-trained model
#
# Run AFTER `snakemake methods_scvi_single` has produced the model dir.

# %%
src_model_dir = scenario_model_dir(REFERENCE_SCENARIO, "scvi_single")
dst = MODEL_OUT / ATLAS_NAME

if src_model_dir.exists():
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.symlink_to(src_model_dir)
    print(f"symlinked {dst} -> {src_model_dir}")
else:
    print(f"NO pipeline model at {src_model_dir} — falling through to inline training")

# %% [markdown]
# ## Path 2 — inline training (fallback)
#
# Only runs if `dst` does not exist after Path 1.

# %%
if not dst.exists():
    import scvi
    import pandas as pd
    from src.data_io import load_parquet_as_anndata, split_by_source

    SOURCE_COL = "Metadata_Source"
    BATCH_KEY = "Metadata_Source"

    adata = load_parquet_as_anndata(scenario_input_parquet(REFERENCE_SCENARIO))
    ref, _ = split_by_source(adata, SOURCE_COL, QUERY_SOURCE)

    params = pd.read_csv(scenario_optuna_csv(REFERENCE_SCENARIO, "scvi_single"))
    p = params.sort_values("total", ascending=False).iloc[0].to_dict()

    # Pipeline shifts to non-negative for ZINB; mirror that here.
    ref.X -= ref.X.min()
    scvi.model.SCVI.setup_anndata(ref, batch_key=BATCH_KEY)
    vae = scvi.model.SCVI(
        ref,
        n_hidden=int(p["params_n_hidden"]),
        n_latent=int(p["params_n_latent"]),
        n_layers=int(p["params_n_layers"]),
        dropout_rate=float(p["params_dropout_rate"]),
        gene_likelihood="zinb",
    )
    vae.train(max_epochs=400, early_stopping=True)
    vae.save(str(dst), overwrite=True, save_anndata=True)
    print(f"inline-trained scVI atlas to {dst}")
