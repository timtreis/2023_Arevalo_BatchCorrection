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
# # 11 — Map a query arm onto the scVI atlas (papermill-parameterized)

# %%
from pathlib import Path
import sys
import time
import json

sys.path.insert(0, str(Path.cwd().parent))
from src.paths import DATA_OUT, MODEL_OUT
import anndata as ad
import scvi

# %% tags=["parameters"]
# --- PARAMETERS (overridden by papermill) ---
QUERY_SOURCE = "source_5"
T_ARM = "tplus"
MAX_EPOCHS = 20
# --- END PARAMETERS ---

# %% [markdown]
# ## Load atlas + query

# %%
atlas_dir = MODEL_OUT / f"scvi_atlas_excl_{QUERY_SOURCE}"
query = ad.read_h5ad(DATA_OUT / "query_arms" / f"{QUERY_SOURCE}_{T_ARM}.h5ad")
# scVI ZINB expects non-negative; mirror pipeline shift.
query.X -= query.X.min()
print(f"atlas: {atlas_dir}, query: {query.shape}")

# %% [markdown]
# ## scArches surgery via SCVI.load_query_data

# %%
scvi.model.SCVI.prepare_query_anndata(query, str(atlas_dir))
model_query = scvi.model.SCVI.load_query_data(query, str(atlas_dir))

t0 = time.perf_counter()
model_query.train(
    max_epochs=MAX_EPOCHS,
    plan_kwargs={"weight_decay": 0.0},
)
elapsed = time.perf_counter() - t0
print(f"mapping took {elapsed:.1f}s")

# %%
query.obsm["X_scvi_mapped"] = model_query.get_latent_representation()

# %% [markdown]
# ## Persist

# %%
out_dir = DATA_OUT / "mapped"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / f"{QUERY_SOURCE}_{T_ARM}_scvi.h5ad"
query.write_h5ad(out_path)
(out_dir / f"{QUERY_SOURCE}_{T_ARM}_scvi_timings.json").write_text(
    json.dumps({"mapping_seconds": elapsed, "max_epochs": MAX_EPOCHS})
)
print(f"wrote {out_path}")
