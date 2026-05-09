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
# > **⚠ This notebook does not work with new (unseen) query batches.**
# >
# > `SysVI.load_query_data` hard-codes `transfer_batch=False`, raising
# > `ValueError` when query contains any batch not in the reference registry.
# > sysVI has no scArches surgery path for new batches.
# >
# > Kept as **scaffolding only**. For real reference mapping use **scPoli**
# > (`21_map_scpoli.ipynb`) or **Symphony** (`31_map_symphony.ipynb`).
# >
# > See `tasks/lessons.md` → "SysVI cannot do reference mapping with new
# > (unseen) batch categories (2026-04-19)" for the full diagnosis.

# %% [markdown]
# # 11 — Map a query arm onto the sysVI atlas (papermill-parameterized)
#
# Uses scArches surgery: freeze reference encoder, fine-tune query-specific encoder
# for a small number of epochs. No X shift needed (sysVI uses Gaussian likelihood).

# %%
from pathlib import Path
import sys
import time
import json

sys.path.insert(0, str(Path.cwd().parent))
from src.paths import DATA_OUT, MODEL_OUT
import anndata as ad
from scvi.external import SysVI

# %% tags=["parameters"]
# --- PARAMETERS (overridden by papermill) ---
QUERY_SOURCE = "source_8"
T_ARM = "tplus"
MAX_EPOCHS = 20
# --- END PARAMETERS ---

# %% [markdown]
# ## Load atlas + query

# %%
atlas_dir = MODEL_OUT / f"sysvi_atlas_excl_{QUERY_SOURCE}"
query = ad.read_h5ad(DATA_OUT / "query_arms" / f"{QUERY_SOURCE}_{T_ARM}.h5ad")
print(f"atlas: {atlas_dir}, query: {query.shape}")

# %% [markdown]
# ## scArches surgery via SysVI.load_query_data

# %%
SysVI.prepare_query_anndata(query, str(atlas_dir))
model_query = SysVI.load_query_data(query, str(atlas_dir))

t0 = time.perf_counter()
model_query.train(
    max_epochs=MAX_EPOCHS,
    plan_kwargs={"weight_decay": 0.0},
)
elapsed = time.perf_counter() - t0
print(f"mapping took {elapsed:.1f}s")

# %%
query.obsm["X_sysvi_mapped"] = model_query.get_latent_representation()

# %% [markdown]
# ## Persist

# %%
out_dir = DATA_OUT / "mapped"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / f"{QUERY_SOURCE}_{T_ARM}_sysvi.h5ad"
query.write_h5ad(out_path)
(out_dir / f"{QUERY_SOURCE}_{T_ARM}_sysvi_timings.json").write_text(
    json.dumps({"mapping_seconds": elapsed, "max_epochs": MAX_EPOCHS})
)
print(f"wrote {out_path}")
