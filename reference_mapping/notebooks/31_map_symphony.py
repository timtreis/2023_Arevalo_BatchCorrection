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
# # 31 — Map a query arm onto the Symphony atlas (papermill-parameterized)

# %%
from pathlib import Path
import sys

sys.path.insert(0, str(Path.cwd().parent))
from src.paths import DATA_OUT, MODEL_OUT
import anndata as ad

# %% tags=["parameters"]
# --- PARAMETERS (overridden by papermill) ---
QUERY_SOURCE = "source_8"
T_ARM = "tplus"  # one of: tplus, tminus, tminus_matched
REFERENCE_ATLAS = "excl_source_8"   # "excl_{source}" for Exp1; "scenario_5_full" for Exp2
BATCH_KEY = "Metadata_Source"
# --- END PARAMETERS ---

# %% [markdown]
# ## Load atlas + query

# %%
atlas_dir = MODEL_OUT / f"symphony_atlas_{REFERENCE_ATLAS}"
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

# symphonypy needs mean/std in ref.var to z-score the query before PCA projection.
# If the atlas was built without sc.pp.scale (old nb30 run), compute them now
# from the raw reference and add them in-memory — no need to rebuild the atlas.
if "mean" not in ref.var.columns or "std" not in ref.var.columns:
    print("mean/std missing from ref.var — computing from raw reference")
    raw_ref = ad.read_h5ad(DATA_OUT / "query_arms" / f"{QUERY_SOURCE}_reference.h5ad")
    sc.pp.scale(raw_ref, zero_center=True)
    ref.var["mean"] = raw_ref.var["mean"].values
    ref.var["std"] = raw_ref.var["std"].values
    del raw_ref
    print("done")

# Project onto reference PCA
t0 = time.perf_counter()
sp.tl.map_embedding(adata_query=query, adata_ref=ref, key=BATCH_KEY, use_genes_column=None)
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
