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
# # 40 — Metrics for one (query, arm, paradigm)
#
# Computes compound mAP, MOA mAP, kBET, iLISI on the mapped query embedding,
# masking out:
#   - Target2 compounds (always)
#   - swap_in compounds (only when arm == tminus_matched, to ensure a
#     fair compound population for the comparison)
#
# Appends one row to `results/acute/<query_source>_metrics.csv`.

# %%
from pathlib import Path
import sys
import pandas as pd
import anndata as ad

sys.path.insert(0, str(Path.cwd().parent))
from src.paths import DATA_OUT, MODEL_OUT, RESULTS_OUT
from src.target2 import load_manifest
from src.query_arms import load_swap_ids
from src import metrics as m

# %% tags=["parameters"]
# --- PARAMETERS (overridden by papermill) ---
QUERY_SOURCE = "source_5"
T_ARM = "tplus"
PARADIGM = "symphony"  # symphony | scvi | scpoli
COMPOUND_COL = "Metadata_JCP2022"
MOA_COL = "Metadata_MoA"  # verify; may be Metadata_target / Metadata_pert_iname
BATCH_KEY = "Metadata_Source"
# --- END PARAMETERS ---

EMBEDDING_KEY = f"X_{PARADIGM}_mapped"

# %% [markdown]
# ## Load mapped query and apply evaluation mask

# %%
query = ad.read_h5ad(DATA_OUT / "mapped" / f"{QUERY_SOURCE}_{T_ARM}_{PARADIGM}.h5ad")
target2 = load_manifest(DATA_OUT / "target2_compounds.json")

eval_mask = ~query.obs[COMPOUND_COL].isin(target2)
if T_ARM == "tminus_matched":
    swap_in = load_swap_ids(DATA_OUT / "query_arms" / f"{QUERY_SOURCE}_swap_ids.json")
    eval_mask &= ~query.obs[COMPOUND_COL].isin(swap_in)

eval_query = query[eval_mask].copy()
print(f"evaluating on {eval_query.shape[0]}/{query.shape[0]} cells")

# %% [markdown]
# ## Per-cell metrics on mapped embedding

# %%
emb = eval_query.obsm[EMBEDDING_KEY]
compound_map = m.compound_map(emb, eval_query.obs[COMPOUND_COL])
moa_map_val = m.moa_map(emb, eval_query.obs[MOA_COL])
print(f"compound_map={compound_map:.4f}  moa_map={moa_map_val:.4f}")

# %% [markdown]
# ## Joint reference + query batch-mixing metrics

# %%
# Find the matching reference embedding next to the atlas.
atlas_dir = MODEL_OUT / f"{PARADIGM}_atlas_excl_{QUERY_SOURCE}"
ref_h5ad = atlas_dir / "reference.h5ad"
if ref_h5ad.exists():
    ref = ad.read_h5ad(ref_h5ad)
    # Use the same key on the ref. Symphony stores X_pca_harmony; scVI/scPoli
    # need an explicit ref-embedding key — skip joint metrics if missing.
    ref_key_candidates = [EMBEDDING_KEY, f"X_{PARADIGM}_ref", "X_pca_harmony"]
    ref_key = next((k for k in ref_key_candidates if k in ref.obsm), None)
    if ref_key is None:
        kbet = ilisi = float("nan")
        print(f"reference embedding missing; skip kBET/iLISI for {PARADIGM}")
    else:
        joint = ad.concat(
            [
                ad.AnnData(X=ref.obsm[ref_key], obs=ref.obs[[BATCH_KEY]]),
                ad.AnnData(X=eval_query.obsm[EMBEDDING_KEY], obs=eval_query.obs[[BATCH_KEY]]),
            ]
        )
        kbet = m.kbet_score(joint.X, joint.obs[BATCH_KEY])
        ilisi = m.ilisi_score(joint.X, joint.obs[BATCH_KEY])
        print(f"kBET={kbet:.4f}  iLISI={ilisi:.4f}")
else:
    kbet = ilisi = float("nan")
    print(f"no reference at {ref_h5ad}; skip kBET/iLISI")

# %% [markdown]
# ## Append to CSV

# %%
row = {
    "query_source": QUERY_SOURCE,
    "t_arm": T_ARM,
    "paradigm": PARADIGM,
    "n_cells_eval": int(eval_query.shape[0]),
    "compound_map": compound_map,
    "moa_map": moa_map_val,
    "kbet": kbet,
    "ilisi": ilisi,
}

csv_path = RESULTS_OUT / "acute" / f"{QUERY_SOURCE}_metrics.csv"
csv_path.parent.mkdir(parents=True, exist_ok=True)
header = not csv_path.exists()
pd.DataFrame([row]).to_csv(
    csv_path, mode="a", header=header, index=False
)
print(f"appended row to {csv_path}")
