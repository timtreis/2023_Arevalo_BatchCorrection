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
#     display_name: RefMap Symphony (GPU)
#     language: python
#     name: refmap-symphony
# ---

# %% [markdown]
# # 30 — Train Harmony reference + Symphony compression
#
# Symphony is a linear projection on top of a Harmony-corrected reference.
# We compute the Harmony correction here, then build a Symphony index that
# new query batches can be projected onto without re-running Harmony.
#
# Cheapest paradigm of the three — start here as a smoke test.

# %%
from pathlib import Path
import sys

sys.path.insert(0, str(Path.cwd().parent))
from src.paths import DATA_OUT, MODEL_OUT
import anndata as ad
import numpy as np
import scanpy as sc

# %% [markdown]
# ## Parameters

# %%
QUERY_SOURCE = ""          # empty = full S5 atlas; "source_8" = excl. that source
REFERENCE_SCENARIO = "scenario_5"
BATCH_KEY = "Metadata_Source"
N_PCS = 50

ATLAS_DIR = (
    MODEL_OUT / f"symphony_atlas_{REFERENCE_SCENARIO}_full"
    if not QUERY_SOURCE
    else MODEL_OUT / f"symphony_atlas_excl_{QUERY_SOURCE}"
)
ATLAS_DIR.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# ## Load reference

# %%
if QUERY_SOURCE:
    ref = ad.read_h5ad(DATA_OUT / "query_arms" / f"{QUERY_SOURCE}_reference.h5ad")
else:
    from src.data_io import load_parquet_as_anndata
    from src.paths import scenario_input_parquet
    ref = load_parquet_as_anndata(scenario_input_parquet(REFERENCE_SCENARIO))
print(ref.shape)

# %% [markdown]
# ## Scale → PCA → Harmony
#
# These three steps are in one cell intentionally: scale must run before PCA,
# and PCA must run before Harmony. Running them separately on a live kernel
# risks using stale in-memory state from a previous run.

# %%
import symphonypy as sp

# Scale stores mean+std in ref.var — required by symphonypy for query projection.
sc.pp.scale(ref, zero_center=True)
assert "mean" in ref.var.columns and "std" in ref.var.columns, \
    "sc.pp.scale did not write mean/std to var"

sc.pp.pca(ref, n_comps=N_PCS, svd_solver="arpack")

sp.pp.harmony_integrate(ref, key=BATCH_KEY, ref_basis_source="X_pca")
print("X_pca_harmony shape:", ref.obsm["X_pca_harmony"].shape)
print("var cols:", ref.var.columns.tolist())

# %% [markdown]
# ## Persist atlas

# %%
ref.write_h5ad(ATLAS_DIR / "reference.h5ad")
print(f"wrote atlas to {ATLAS_DIR}")
