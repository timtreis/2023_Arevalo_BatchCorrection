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
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scarches.models.scpoli import scPoli

# %% tags=["parameters"]
# --- PARAMETERS (overridden by papermill) ---
QUERY_SOURCE = "source_8"
T_ARM = "tplus"
REFERENCE_ATLAS = "excl_source_8"   # "excl_{source}" for Exp1; "scenario_5_full" for Exp2
REFERENCE_SCENARIO = "scenario_5"   # used to look up shared compounds for labeled_indices
N_EPOCHS = 50
PRETRAIN_EPOCHS = 0
# --- END PARAMETERS ---

# %% [markdown]
# ## Load atlas + query

# %%
atlas_dir = MODEL_OUT / f"scpoli_atlas_{REFERENCE_ATLAS}"
query = ad.read_h5ad(DATA_OUT / "query_arms" / f"{QUERY_SOURCE}_{T_ARM}.h5ad")
print(f"atlas: {atlas_dir}, query shape: {query.shape}")

# Zero-pad features present in the reference atlas but missing from the query.
# Build padded AnnData manually to preserve obs metadata (ad.concat drops columns).
var_names_path = atlas_dir / "var_names.txt"
if var_names_path.exists():
    ref_var_names = list(np.loadtxt(var_names_path, dtype=str))
    X_raw = query.X.toarray() if sp.issparse(query.X) else np.array(query.X, dtype=np.float32)
    orig_idx = {v: i for i, v in enumerate(query.var_names)}
    missing = [v for v in ref_var_names if v not in orig_idx]
    if missing:
        print(f"zero-padding {len(missing)} features absent from query")
        X_full = np.zeros((query.shape[0], len(ref_var_names)), dtype=np.float32)
        for new_i, v in enumerate(ref_var_names):
            if v in orig_idx:
                X_full[:, new_i] = X_raw[:, orig_idx[v]]
    else:
        X_full = np.stack([X_raw[:, orig_idx[v]] for v in ref_var_names], axis=1)
    query = ad.AnnData(
        X=X_full,
        obs=query.obs.copy(),
        var=pd.DataFrame(index=ref_var_names),
        uns=query.uns.copy(),
    )
    print(f"query after padding: {query.shape}")
else:
    print("WARNING: var_names.txt not found — skipping zero-padding")

# %% [markdown]
# ## scArches load_query_data + fine-tune

# %%
import pandas as pd
import numpy as np
import json
from src.paths import scenario_input_parquet
from src.target2 import load_manifest

# Read atlas training config to determine which wells should be labeled.
# If T2-only prototypes were used, only T2-matching wells are labeled.
config_path = atlas_dir / "training_config.json"
if config_path.exists():
    atlas_config = json.loads(config_path.read_text())
    t2_only_prototypes = atlas_config.get("t2_only_prototypes", False)
else:
    t2_only_prototypes = False

if t2_only_prototypes:
    t2_ids = load_manifest(DATA_OUT / "target2_compounds.json")
    labeled_mask = query.obs["Metadata_JCP2022"].isin(t2_ids)
    # Add _cell_type col matching training convention (T2 id or "BACKGROUND")
    query.obs["_cell_type"] = query.obs["Metadata_JCP2022"].astype(str).where(
        labeled_mask, "BACKGROUND"
    )
    print(f"T2-only labeled_indices: {labeled_mask.sum()}/{len(query)} wells")
else:
    ref_compounds = set(
        pd.read_parquet(scenario_input_parquet(REFERENCE_SCENARIO), columns=["Metadata_JCP2022"])
        ["Metadata_JCP2022"].unique()
    )
    labeled_mask = query.obs["Metadata_JCP2022"].isin(ref_compounds)
    print(f"All-compound labeled_indices: {labeled_mask.sum()}/{len(query)} wells")

labeled_indices = list(np.where(labeled_mask)[0])

# scPoli requires dense input; densify if zero-padding produced a sparse matrix.
if sp.issparse(query.X):
    query.X = query.X.toarray()

model_query = scPoli.load_query_data(
    adata=query,
    reference_model=str(atlas_dir),
    labeled_indices=labeled_indices,
)

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

# Save fine-tuned query model so mapping can be reproduced without retraining
model_save_path = out_dir / f"{QUERY_SOURCE}_{T_ARM}_scpoli_model"
model_query.save(str(model_save_path), overwrite=True)
print(f"saved query model to {model_save_path}")
