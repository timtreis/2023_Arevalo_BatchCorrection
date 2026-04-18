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
# # 00 — Derive JUMP-Target2 manifest
#
# Read pipeline metadata, intersect compound sets across sources for the
# `TARGET2` plate type, cache the resulting compound IDs to JSON.
#
# **First run**: column-name guesses may be wrong. The helper raises with
# the available `Metadata_*` columns listed — paste the right ones below.

# %%
from pathlib import Path
import sys

sys.path.insert(0, str(Path.cwd().parent))
from src.target2 import derive_target2_manifest, cache_manifest
from src.paths import PIPELINE_OUT, DATA_OUT, scenario_method_parquet

# %% [markdown]
# ## Pick a metadata source
#
# Any pipeline parquet works — they all carry the same `Metadata_*` columns.
# We use the largest available scenario for maximum source coverage.

# %%
# Edit if the chosen scenario doesn't have all sources.
META_SOURCE = scenario_method_parquet("scenario_5", "scvi_single")
assert META_SOURCE.exists(), f"missing {META_SOURCE} — pick another scenario"
META_SOURCE

# %% [markdown]
# ## Inspect available `Metadata_*` columns
#
# Run this if column-name guesses fail.

# %%
import pandas as pd
head = pd.read_parquet(META_SOURCE).head(0)
[c for c in head.columns if c.startswith("Metadata_")]

# %% [markdown]
# ## Derive the manifest
#
# Override `plate_col` / `source_col` / `compound_col` if defaults don't
# match what you found above.

# %%
target2 = derive_target2_manifest(META_SOURCE)
print(f"Found {len(target2)} Target2 compounds")
sorted(target2)[:5]

# %% [markdown]
# ## Cache to JSON

# %%
cache_path = DATA_OUT / "target2_compounds.json"
cache_manifest(target2, cache_path)
print(f"Wrote {cache_path}")
