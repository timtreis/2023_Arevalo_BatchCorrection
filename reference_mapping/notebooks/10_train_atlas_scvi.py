# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: RefMap scVI (GPU)
#     language: python
#     name: refmap-scvi
# ---

# %% [markdown]
# # 10 — scVI atlas
#
# Train a scVI reference atlas on scenario_5 (all 5 sources).
# Uses HPO best params from `optuna_scvi_single.csv`.
# Saves model + reference latent embedding so nb11 can map queries
# and nb40 can compute joint reference+query metrics.

# %%
from pathlib import Path
import sys

sys.path.insert(0, str(Path.cwd().parent))
from src.paths import MODEL_OUT, scenario_input_parquet, scenario_optuna_csv

# %% [markdown]
# ## Parameters

# %%
QUERY_SOURCE = ""          # empty = use all sources; "source_8" = exclude that source
REFERENCE_SCENARIO = "scenario_5"
BATCH_KEY = "Metadata_Source"
LABEL_KEY = "Metadata_JCP2022"

ATLAS_NAME = (
    f"scvi_atlas_{REFERENCE_SCENARIO}_full"
    if not QUERY_SOURCE
    else f"scvi_atlas_excl_{QUERY_SOURCE}"
)
atlas_dir = MODEL_OUT / ATLAS_NAME
atlas_dir.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# ## Load reference data

# %%
from src.data_io import load_parquet_as_anndata, split_by_source
import json

adata = load_parquet_as_anndata(scenario_input_parquet(REFERENCE_SCENARIO))
ref = adata if not QUERY_SOURCE else split_by_source(adata, BATCH_KEY, QUERY_SOURCE)[0]
print(f"reference: {ref.shape}")

# HPO uses gene_likelihood="zinb" (default), which requires non-negative data.
# Shift X so the minimum is 0 and save the offset for query notebooks.
X_min = float(ref.X.min())
ref.X = ref.X - X_min
(atlas_dir / "X_min.json").write_text(json.dumps({"X_min": X_min}))
print(f"shifted X by {X_min:.4f} (saved to atlas_dir/X_min.json)")

# %% [markdown]
# ## Load HPO best params

# %%
import pandas as pd

params = pd.read_csv(scenario_optuna_csv(REFERENCE_SCENARIO, "scvi_single"))
p = params.sort_values("total", ascending=False).iloc[0].to_dict()
print("Best HPO params:")
for k, v in p.items():
    if k.startswith("params_"):
        print(f"  {k}: {v}")

# %% [markdown]
# ## Train scVI atlas

# %%
from scvi.model import SCVI

SCVI.setup_anndata(ref, batch_key=BATCH_KEY)

model = SCVI(
    ref,
    n_hidden=int(p["params_n_hidden"]),
    n_latent=int(p["params_n_latent"]),
    n_layers=int(p["params_n_layers"]),
    dropout_rate=float(p["params_dropout_rate"]),
    gene_likelihood="zinb",
)
model.train(
    max_epochs=400,
    early_stopping=True,
    early_stopping_monitor="elbo_validation",
)
model.save(str(atlas_dir), overwrite=True)
print(f"saved scVI atlas to {atlas_dir}")

# %% [markdown]
# ## Save reference embedding for nb40 joint metrics

# %%
import anndata as ad
import numpy as np

ref_emb = model.get_latent_representation()
print(f"reference embedding shape: {ref_emb.shape}")

ref_out = ad.AnnData(obs=ref.obs[[BATCH_KEY, LABEL_KEY]])
ref_out.obsm["X_scvi_mapped"] = ref_emb
ref_out.write_h5ad(atlas_dir / "reference.h5ad")
print(f"saved reference embedding to {atlas_dir / 'reference.h5ad'}")
