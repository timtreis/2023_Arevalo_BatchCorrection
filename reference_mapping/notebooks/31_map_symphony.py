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
# # 31 — Map a query arm onto the Symphony atlas (papermill-parameterized)

# %%
from pathlib import Path
import sys

sys.path.insert(0, str(Path.cwd().parent))
from src.paths import DATA_OUT, MODEL_OUT
import anndata as ad

# %% tags=["parameters"]
# --- PARAMETERS (overridden by papermill) ---
QUERY_SOURCE = "source_5"
T_ARM = "tplus"  # one of: tplus, tminus, tminus_matched
BATCH_KEY = "Metadata_Source"
# --- END PARAMETERS ---

# %% [markdown]
# ## Load atlas + query

# %%
atlas_dir = MODEL_OUT / f"symphony_atlas_excl_{QUERY_SOURCE}"
ref = ad.read_h5ad(atlas_dir / "reference.h5ad")
query = ad.read_h5ad(DATA_OUT / "query_arms" / f"{QUERY_SOURCE}_{T_ARM}.h5ad")
print(f"ref: {ref.shape}, query: {query.shape}")

# %% [markdown]
# ## Symphony projection
#
# `symphonypy` exposes `sp.tl.map_embedding` which projects new cells onto the
# Harmony-corrected PCA basis using the reference's centroids/sigmas.

# %%
import symphonypy as sp
import scanpy as sc
import time

# Project onto reference PCA
t0 = time.perf_counter()
sp.tl.map_embedding(adata_query=query, adata_ref=ref, key=BATCH_KEY)
elapsed = time.perf_counter() - t0
print(f"mapping took {elapsed:.1f}s")

# The mapped embedding lands in query.obsm — store under a stable name
query.obsm["X_symphony_mapped"] = query.obsm["X_pca_harmony"]

# %% [markdown]
# ## Persist mapped query + timing

# %%
out_dir = DATA_OUT / "mapped"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / f"{QUERY_SOURCE}_{T_ARM}_symphony.h5ad"
query.write_h5ad(out_path)

import json
(out_dir / f"{QUERY_SOURCE}_{T_ARM}_symphony_timings.json").write_text(
    json.dumps({"mapping_seconds": elapsed})
)
print(f"wrote {out_path}")
