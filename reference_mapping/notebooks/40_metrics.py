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
# # 40 — Metrics for one (query, arm, paradigm)
#
# Computes compound ASW, batch ASW, kBET, iLISI on the mapped query embedding,
# masking out:
#   - Target2 compounds (always)
#   - swap_in compounds (only when arm == tminus_matched, for a fair comparison)
#
# Appends one row to `results/acute/<query_source>_metrics.csv`.

# %%
from pathlib import Path
import sys
import pandas as pd
import numpy as np
import anndata as ad

sys.path.insert(0, str(Path.cwd().parent))
from src.paths import DATA_OUT, MODEL_OUT, RESULTS_OUT
from src.target2 import load_manifest
from src.query_arms import load_swap_ids
from src import metrics as m

# %% tags=["parameters"]
# --- PARAMETERS (overridden by papermill) ---
QUERY_SOURCE = "source_8"
T_ARM = "tplus"
PARADIGM = "symphony"  # symphony | scpoli | scvi
REFERENCE_ATLAS = "excl_source_8"   # "excl_{source}" for Exp1; "scenario_5_full" for Exp2
COMPOUND_COL = "Metadata_JCP2022"
BATCH_KEY = "Metadata_Source"
N_NEIGHBORS = 50
# --- END PARAMETERS ---

# %%
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
# ## Compound coherence (ASW over compound identity)

# %%
emb = eval_query.obsm[EMBEDDING_KEY]
compound_asw_val = m.compound_asw(emb, eval_query.obs[COMPOUND_COL])
print(f"compound_asw={compound_asw_val:.4f}")

# %% [markdown]
# ## Joint reference + query batch-mixing metrics

# %%
atlas_dir = MODEL_OUT / f"{PARADIGM}_atlas_{REFERENCE_ATLAS}"
ref_h5ad = atlas_dir / "reference.h5ad"

kbet = ilisi = batch_asw_val = precision_at_k = float("nan")

if ref_h5ad.exists():
    ref = ad.read_h5ad(ref_h5ad)
    ref_key_candidates = [EMBEDDING_KEY, f"X_{PARADIGM}_ref", "X_pca_harmony"]
    ref_key = next((k for k in ref_key_candidates if k in ref.obsm), None)

    if ref_key is None:
        print(f"reference embedding missing; skip batch metrics for {PARADIGM}")
    else:
        # Build joint embedding: reference + eval query
        joint = ad.concat(
            [
                ad.AnnData(
                    X=ref.obsm[ref_key],
                    obs=ref.obs[[BATCH_KEY, COMPOUND_COL]],
                ),
                ad.AnnData(
                    X=eval_query.obsm[EMBEDDING_KEY],
                    obs=eval_query.obs[[BATCH_KEY, COMPOUND_COL]],
                ),
            ]
        )
        joint_emb = joint.X if isinstance(joint.X, np.ndarray) else joint.X.toarray()

        # Batch ASW (no kNN needed)
        batch_asw_val = m.batch_asw(
            joint_emb, joint.obs[COMPOUND_COL], joint.obs[BATCH_KEY]
        )
        print(f"batch_asw={batch_asw_val:.4f}")

        # kBET and iLISI via pre-computed kNN
        from scib_metrics.nearest_neighbors import pynndescent
        nn = pynndescent(joint_emb, n_neighbors=N_NEIGHBORS)
        kbet = m.kbet_score(nn, joint.obs[BATCH_KEY])
        ilisi = m.ilisi_score(nn, joint.obs[BATCH_KEY])
        print(f"kBET={kbet:.4f}  iLISI={ilisi:.4f}")

        # Cross-source compound precision@k
        from sklearn.neighbors import NearestNeighbors
        K = 10
        ref_emb_arr = ref.obsm[ref_key]
        query_emb_arr = eval_query.obsm[EMBEDDING_KEY]
        nbrs = NearestNeighbors(n_neighbors=K, metric="cosine").fit(ref_emb_arr)
        _, nn_indices = nbrs.kneighbors(query_emb_arr)
        ref_cpds = ref.obs[COMPOUND_COL].astype(str).values
        query_cpds = eval_query.obs[COMPOUND_COL].astype(str).values
        match = ref_cpds[nn_indices] == query_cpds[:, None]
        precision_at_k = float(match.mean())
        print(f"precision@{K}={precision_at_k:.4f}")
else:
    print(f"no reference at {ref_h5ad}; skip batch metrics")

# %% [markdown]
# ## Append to CSV

# %%
row = {
    "query_source": QUERY_SOURCE,
    "t_arm": T_ARM,
    "paradigm": PARADIGM,
    "n_cells_eval": int(eval_query.shape[0]),
    "compound_asw": compound_asw_val,
    "batch_asw": batch_asw_val,
    "kbet": kbet,
    "ilisi": ilisi,
    "precision_at_10": precision_at_k,
}

csv_path = RESULTS_OUT / "acute" / f"{QUERY_SOURCE}_metrics.csv"
csv_path.parent.mkdir(parents=True, exist_ok=True)
header = not csv_path.exists()
pd.DataFrame([row]).to_csv(
    csv_path, mode="a", header=header, index=False
)
print(f"appended row to {csv_path}")
