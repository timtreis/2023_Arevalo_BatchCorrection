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
from src.data_io import load_parquet_as_anndata, split_by_source

# %% [markdown]
# ## Parameters

# %%
# Acute first-pass: source_5 is a CV8000 confocal source matched to wave1
# microscope family and has Target2 plates.
QUERY_SOURCE = "source_5"
# Reference scenario whose preprocessed parquet contains the query source
# alongside the rest of wave1.
REFERENCE_SCENARIO = "scenario_5"  # contains source_5 + wave1 sources; verify
COMPOUND_COL = "Metadata_JCP2022"
SOURCE_COL = "Metadata_Source"
SEED = 42

# %% [markdown]
# ## Load preprocessed (pre-correction) data and split off the query source

# %%
input_pq = scenario_input_parquet(REFERENCE_SCENARIO)
adata = load_parquet_as_anndata(input_pq)
print(f"loaded {adata.shape} from {input_pq}")
adata.obs[SOURCE_COL].value_counts()

# %%
ref, query = split_by_source(adata, SOURCE_COL, QUERY_SOURCE)
print(f"reference: {ref.shape}, query: {query.shape}")

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
