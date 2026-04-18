# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
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
QUERY_SOURCE = "source_5"
BATCH_KEY = "Metadata_Source"
N_PCS = 50

REF_PATH = DATA_OUT / "query_arms" / f"{QUERY_SOURCE}_reference.h5ad"
ATLAS_DIR = MODEL_OUT / f"symphony_atlas_excl_{QUERY_SOURCE}"
ATLAS_DIR.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# ## Load reference

# %%
ref = ad.read_h5ad(REF_PATH)
print(ref.shape)

# %% [markdown]
# ## PCA + Harmony

# %%
sc.pp.pca(ref, n_comps=N_PCS, svd_solver="arpack")

# %%
import symphonypy as sp

# Harmony correction in symphonypy lives under sp.pp.harmony_integrate; it
# wraps harmonypy and stores the corrected embedding in `ref.obsm["X_pca_harmony"]`
# alongside the components needed to project queries.
sp.pp.harmony_integrate(ref, key=BATCH_KEY, ref_basis_source="X_pca")
print("X_pca_harmony shape:", ref.obsm["X_pca_harmony"].shape)

# %% [markdown]
# ## Persist atlas

# %%
ref.write_h5ad(ATLAS_DIR / "reference.h5ad")
print(f"wrote atlas to {ATLAS_DIR}")
