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
# # 11 — Map a query arm onto the scVI atlas (papermill-parameterized)
#
# Uses scArches-style `load_query_data` to fine-tune new batch embeddings
# while freezing the reference encoder. No prototype mechanism (unlike scPoli).

# %%
from pathlib import Path
import sys
import time
import json
import numpy as np

sys.path.insert(0, str(Path.cwd().parent))
from src.paths import DATA_OUT, MODEL_OUT
import anndata as ad
from scvi.model import SCVI

# %% tags=["parameters"]
# --- PARAMETERS (overridden by papermill) ---
QUERY_SOURCE = "source_5"
T_ARM = "tplus"
REFERENCE_ATLAS = "scenario_5_full"   # "excl_{source}" or "scenario_5_full"
N_EPOCHS = 200
# --- END PARAMETERS ---

# %% [markdown]
# ## Load atlas + query

# %%
import pandas as pd
import scipy.sparse as sp
import torch

atlas_dir = MODEL_OUT / f"scvi_atlas_{REFERENCE_ATLAS}"
query = ad.read_h5ad(DATA_OUT / "query_arms" / f"{QUERY_SOURCE}_{T_ARM}.h5ad")
print(f"atlas: {atlas_dir}, query shape: {query.shape}")

# Zero-pad features present in the reference atlas but missing from the query,
# then apply the non-negative scalar shift used during atlas training.
# Padding is in raw (pre-shift) space so missing features get 0 (= mean for
# MAD-normalised data). We build the padded AnnData manually to preserve obs.
pt_data = torch.load(atlas_dir / "model.pt", map_location="cpu", weights_only=False)
ref_var_names = list(pt_data["var_names"])
X_min = json.loads((atlas_dir / "X_min.json").read_text())["X_min"]

X_raw = query.X.toarray() if sp.issparse(query.X) else np.array(query.X, dtype=np.float32)
if set(ref_var_names) - set(query.var_names):
    missing = [v for v in ref_var_names if v not in query.var_names]
    print(f"zero-padding {len(missing)} features absent from query (raw space)")
    orig_idx = {v: i for i, v in enumerate(query.var_names)}
    X_full = np.zeros((query.shape[0], len(ref_var_names)), dtype=np.float32)
    for new_i, v in enumerate(ref_var_names):
        if v in orig_idx:
            X_full[:, new_i] = X_raw[:, orig_idx[v]]
    # X_full[:, missing positions] stays 0 (raw-space mean for MAD data)
else:
    orig_idx = {v: i for i, v in enumerate(query.var_names)}
    X_full = np.stack([X_raw[:, orig_idx[v]] for v in ref_var_names], axis=1)

query = ad.AnnData(
    X=X_full - X_min,
    obs=query.obs.copy(),
    var=pd.DataFrame(index=ref_var_names),
    uns=query.uns.copy(),
)
print(f"query after padding + shift: {query.shape}, X_min={X_min:.4f}")

# %% [markdown]
# ## scArches load_query_data + fine-tune

# %%
model_query = SCVI.load_query_data(
    adata=query,
    reference_model=str(atlas_dir),
    inplace_subset_query_vars=True,
)

t0 = time.perf_counter()
model_query.train(
    max_epochs=N_EPOCHS,
    plan_kwargs={"weight_decay": 0.0},
    early_stopping=True,
    early_stopping_monitor="elbo_validation",
)
elapsed = time.perf_counter() - t0
print(f"mapping took {elapsed:.1f}s")

# %%
query.obsm["X_scvi_mapped"] = model_query.get_latent_representation()

# %% [markdown]
# ## Persist mapped query + query model

# %%
out_dir = DATA_OUT / "mapped"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / f"{QUERY_SOURCE}_{T_ARM}_scvi.h5ad"
query.write_h5ad(out_path)
(out_dir / f"{QUERY_SOURCE}_{T_ARM}_scvi_timings.json").write_text(
    json.dumps({"mapping_seconds": elapsed, "n_epochs": N_EPOCHS})
)
print(f"wrote {out_path}")

model_save_path = out_dir / f"{QUERY_SOURCE}_{T_ARM}_scvi_model"
model_query.save(str(model_save_path), overwrite=True)
print(f"saved query model to {model_save_path}")
