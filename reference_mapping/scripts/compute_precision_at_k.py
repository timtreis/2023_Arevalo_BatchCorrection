"""Compute precision@k at k=10,20,50,100 for all mapped query arms × methods.

Two variants per (k, arm, method):
  - all:    all non-TARGET2 eval cells (same mask as nb40)
  - shared: further restricted to compounds present in BOTH query and reference

Outputs: results/precision_at_k.csv
"""
import sys
from pathlib import Path
import json

import numpy as np
import pandas as pd
import anndata as ad
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.paths import DATA_OUT, MODEL_OUT, RESULTS_OUT
from src.target2 import load_manifest
from src.query_arms import load_swap_ids

KS = [10, 20, 50, 100]
QUERY_SOURCE = "source_5"
PARADIGMS = ["symphony", "scvi", "scpoli"]
ARMS = ["tplus", "tminus", "tminus_matched"]
COMPOUND_COL = "Metadata_JCP2022"

EMBEDDING_KEYS = {
    "symphony": "X_pca_harmony",
    "scvi":     "X_scvi_mapped",
    "scpoli":   "X_scpoli_mapped",
}

REF_KEY_FALLBACKS = {
    "symphony": ["X_pca_harmony", "X_pca_reference"],
    "scvi":     ["X_scvi_mapped"],
    "scpoli":   ["X_scpoli_mapped"],
}

target2 = load_manifest(DATA_OUT / "target2_compounds.json")

rows = []

for paradigm in PARADIGMS:
    atlas_dir = MODEL_OUT / f"{paradigm}_atlas_scenario_5_full"
    ref_h5ad = atlas_dir / "reference.h5ad"
    if not ref_h5ad.exists():
        print(f"  SKIP {paradigm}: no reference.h5ad at {ref_h5ad}")
        continue

    ref = ad.read_h5ad(ref_h5ad)
    ref_key = next(
        (k for k in REF_KEY_FALLBACKS[paradigm] if k in ref.obsm), None
    )
    if ref_key is None:
        print(f"  SKIP {paradigm}: no embedding key in reference")
        continue

    ref_emb = ref.obsm[ref_key]
    ref_cpds = ref.obs[COMPOUND_COL].astype(str).values
    ref_cpd_set = set(ref_cpds)
    print(f"\n[{paradigm}] ref: {ref_emb.shape}, {len(ref_cpd_set)} unique compounds")

    # Build kNN index over reference (max k = 100)
    max_k = max(KS)
    nbrs = NearestNeighbors(n_neighbors=max_k, metric="cosine").fit(ref_emb)

    for arm in ARMS:
        h5ad = DATA_OUT / "mapped" / f"{QUERY_SOURCE}_{arm}_{paradigm}.h5ad"
        if not h5ad.exists():
            print(f"  SKIP {paradigm}/{arm}: no mapped h5ad")
            continue

        query = ad.read_h5ad(h5ad)
        emb_key = EMBEDDING_KEYS[paradigm]

        # Eval mask: remove TARGET2 (and swap_in for tminus_matched)
        eval_mask = ~query.obs[COMPOUND_COL].isin(target2)
        if arm == "tminus_matched":
            swap_in = load_swap_ids(DATA_OUT / "query_arms" / f"{QUERY_SOURCE}_swap_ids.json")
            eval_mask &= ~query.obs[COMPOUND_COL].isin(swap_in)

        eval_q = query[eval_mask]
        query_cpds_all = eval_q.obs[COMPOUND_COL].astype(str).values
        query_emb_all = eval_q.obsm[emb_key]

        # Shared mask: further restrict to compounds in both query and reference
        shared_mask = np.isin(query_cpds_all, list(ref_cpd_set))
        query_emb_shared = query_emb_all[shared_mask]
        query_cpds_shared = query_cpds_all[shared_mask]

        print(f"  {arm}: {len(query_cpds_all)} all-eval, "
              f"{shared_mask.sum()} shared-compound cells "
              f"({len(set(query_cpds_shared))} unique compounds)")

        # Single kNN query for all cells (max_k neighbors), slice per k
        _, nn_all = nbrs.kneighbors(query_emb_all)
        _, nn_shared = nbrs.kneighbors(query_emb_shared)

        row = {
            "query_source": QUERY_SOURCE,
            "t_arm": arm,
            "paradigm": paradigm,
            "n_eval_all": int(len(query_cpds_all)),
            "n_eval_shared": int(shared_mask.sum()),
            "n_shared_compounds": int(len(set(query_cpds_shared))),
        }

        for k in KS:
            match_all = (ref_cpds[nn_all[:, :k]] == query_cpds_all[:, None]).mean()
            match_shared = (ref_cpds[nn_shared[:, :k]] == query_cpds_shared[:, None]).mean()
            row[f"precision_at_{k}_all"] = float(match_all)
            row[f"precision_at_{k}_shared"] = float(match_shared)

        rows.append(row)
        print(f"    p@10 all={row['precision_at_10_all']:.2e}  "
              f"shared={row['precision_at_10_shared']:.2e}")

out = pd.DataFrame(rows)
out_path = RESULTS_OUT / "precision_at_k.csv"
out.to_csv(out_path, index=False)
print(f"\nWrote {out_path}")
print(out.to_string())
