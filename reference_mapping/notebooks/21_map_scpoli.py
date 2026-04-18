# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: RefMap scPoli (GPU)
#     language: python
#     name: refmap-scpoli
# ---

# %% [markdown]
# # 21 — Map a query arm onto the scPoli atlas (papermill-parameterized)

# %%
from pathlib import Path
import sys
import time
import json

sys.path.insert(0, str(Path.cwd().parent))
from src.paths import DATA_OUT, MODEL_OUT
import anndata as ad
from scarches.models.scpoli import scPoli

# %% tags=["parameters"]
# --- PARAMETERS (overridden by papermill) ---
QUERY_SOURCE = "source_5"
T_ARM = "tplus"
N_EPOCHS = 50
PRETRAIN_EPOCHS = 25
# --- END PARAMETERS ---

# %% [markdown]
# ## Load atlas + query

# %%
atlas_dir = MODEL_OUT / f"scpoli_atlas_excl_{QUERY_SOURCE}"
query = ad.read_h5ad(DATA_OUT / "query_arms" / f"{QUERY_SOURCE}_{T_ARM}.h5ad")
print(f"atlas: {atlas_dir}, query: {query.shape}")

# %% [markdown]
# ## scArches load_query_data + fine-tune

# %%
model_query = scPoli.load_query_data(adata=query, reference_model=str(atlas_dir))

t0 = time.perf_counter()
model_query.train(
    n_epochs=N_EPOCHS,
    pretraining_epochs=PRETRAIN_EPOCHS,
    eta=10,
)
elapsed = time.perf_counter() - t0
print(f"mapping took {elapsed:.1f}s")

# %%
model_query.model.eval()
query.obsm["X_scpoli_mapped"] = model_query.get_latent(query, mean=True)

# %% [markdown]
# ## Persist

# %%
out_dir = DATA_OUT / "mapped"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / f"{QUERY_SOURCE}_{T_ARM}_scpoli.h5ad"
query.write_h5ad(out_path)
(out_dir / f"{QUERY_SOURCE}_{T_ARM}_scpoli_timings.json").write_text(
    json.dumps({"mapping_seconds": elapsed, "n_epochs": N_EPOCHS})
)
print(f"wrote {out_path}")
