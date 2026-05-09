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
# # 01 — Prepare query arms (T+, T-, T-_matched)
#
# Build three AnnData objects for the held-out query source by stripping
# Target2 in different ways. See `src/query_arms.py` for the design rationale.

# %%
from pathlib import Path
import sys

sys.path.insert(0, str(Path.cwd().parent))
from src.paths import DATA_OUT, scenario_input_parquet
from src.target2 import load_manifest
from src.query_arms import make_tplus, make_tminus, make_tminus_matched
from src.data_io import load_parquet_as_anndata

# %% [markdown]
# ## Parameters

# %% tags=["parameters"]
QUERY_SOURCE = "source_5"
QUERY_SCENARIO = "scenario_wave2"  # parquet containing the query source (used if QUERY_PARQUET is empty)
QUERY_PARQUET = ""  # if non-empty, load directly from this path instead of QUERY_SCENARIO
REFERENCE_SCENARIO = "scenario_5"  # full S5; written as reference h5ad for downstream
COMPOUND_COL = "Metadata_JCP2022"
SOURCE_COL = "Metadata_Source"
SEED = 42

# %%
# Resolve query parquet path
from pathlib import Path as _Path
_query_pq = _Path(QUERY_PARQUET) if QUERY_PARQUET else scenario_input_parquet(QUERY_SCENARIO)

# %% [markdown]
# ## Load preprocessed (pre-correction) data and split off the query source

# %%
query_full = load_parquet_as_anndata(_query_pq)
print(f"loaded {query_full.shape} from {_query_pq}")
query_full.obs[SOURCE_COL].value_counts()

# %%
query = query_full[query_full.obs[SOURCE_COL] == QUERY_SOURCE].copy()
print(f"query ({QUERY_SOURCE}): {query.shape}")

ref = load_parquet_as_anndata(scenario_input_parquet(REFERENCE_SCENARIO))
print(f"reference: {ref.shape}")

# %% [markdown]
# ## Build the three arms

# %%
target2 = load_manifest(DATA_OUT / "target2_compounds.json")
print(f"{len(target2)} Target2 compounds in manifest")

# %%
tplus = make_tplus(query)
tminus = make_tminus(query, target2, COMPOUND_COL)
swap_log = DATA_OUT / "query_arms" / f"{QUERY_SOURCE}_swap_ids.json"
tminus_matched = make_tminus_matched(
    query, target2, COMPOUND_COL, seed=SEED, swap_log_path=swap_log
)

print(f"T+         : {tplus.shape}")
print(f"T-         : {tminus.shape}")
print(f"T-_matched : {tminus_matched.shape}")

# %% [markdown]
# ## Persist h5ads

# %%
out_dir = DATA_OUT / "query_arms"
out_dir.mkdir(parents=True, exist_ok=True)

tplus.write_h5ad(out_dir / f"{QUERY_SOURCE}_tplus.h5ad")
tminus.write_h5ad(out_dir / f"{QUERY_SOURCE}_tminus.h5ad")
tminus_matched.write_h5ad(out_dir / f"{QUERY_SOURCE}_tminus_matched.h5ad")

# Also persist the reference for downstream atlas builds.
ref.write_h5ad(DATA_OUT / "query_arms" / f"{QUERY_SOURCE}_reference.h5ad")
print("wrote arm h5ads")
