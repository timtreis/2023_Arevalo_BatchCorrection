"""Generate integration UMAPs for all completed (query_source, paradigm, arm) combinations.

For each combination:
  1. Load reference embedding + mapped query embedding
  2. Concatenate and SHUFFLE rows so no single source dominates the z-order
  3. Compute UMAP on the joint embedding
  4. Plot 2 panels: colored by Metadata_Source (batch) and by ref/query split

Outputs: results/umaps/{query}_{paradigm}_{arm}.png

Run with: pixi run -e symphony python scripts/generate_umaps.py
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.paths import DATA_OUT, MODEL_OUT, RESULTS_OUT

UMAP_OUT = RESULTS_OUT / "umaps"
UMAP_OUT.mkdir(parents=True, exist_ok=True)

BATCH_KEY = "Metadata_Source"
COMPOUND_COL = "Metadata_JCP2022"

# Embedding key per paradigm in the reference h5ad and query h5ad
REF_EMB_KEYS = {
    "symphony": "X_pca_harmony",
    "scvi":     "X_scvi_mapped",
    "scpoli":   "X_scpoli_mapped",
}
QUERY_EMB_KEY = {
    "symphony": "X_symphony_mapped",
    "scvi":     "X_scvi_mapped",
    "scpoli":   "X_scpoli_mapped",
}
ATLAS_NAMES = {
    "source_5": "scenario_5_full",
}

SEED = 42


def load_joint(query_source: str, paradigm: str, arm: str) -> ad.AnnData | None:
    atlas_name = ATLAS_NAMES[query_source]
    ref_path = MODEL_OUT / f"{paradigm}_atlas_{atlas_name}" / "reference.h5ad"
    query_path = DATA_OUT / "mapped" / f"{query_source}_{arm}_{paradigm}.h5ad"

    if not ref_path.exists():
        print(f"  skip: no reference at {ref_path}")
        return None
    if not query_path.exists():
        print(f"  skip: no mapped query at {query_path}")
        return None

    ref = ad.read_h5ad(ref_path)
    query = ad.read_h5ad(query_path)

    ref_key = REF_EMB_KEYS[paradigm]
    q_key = QUERY_EMB_KEY[paradigm]

    if ref_key not in ref.obsm:
        print(f"  skip: {ref_key} missing from reference obsm ({list(ref.obsm.keys())})")
        return None
    if q_key not in query.obsm:
        print(f"  skip: {q_key} missing from query obsm ({list(query.obsm.keys())})")
        return None

    ref_emb = ref.obsm[ref_key]
    q_emb = query.obsm[q_key]

    # Build joint obs table
    def _get_col(adata, col):
        return adata.obs[col].astype(str).values if col in adata.obs else np.full(len(adata), "?")

    n_ref = len(ref_emb)
    obs = pd.DataFrame({
        BATCH_KEY: np.concatenate([_get_col(ref, BATCH_KEY), _get_col(query, BATCH_KEY)]),
        COMPOUND_COL: np.concatenate([_get_col(ref, COMPOUND_COL), _get_col(query, COMPOUND_COL)]),
        "dataset": ["reference"] * n_ref + ["query"] * len(q_emb),
    })

    joint = ad.AnnData(
        X=np.vstack([ref_emb, q_emb]).astype(np.float32),
        obs=obs,
    )

    # SHUFFLE so z-order doesn't bias visualization
    rng = np.random.default_rng(SEED)
    idx = rng.permutation(len(joint))
    return joint[idx].copy()


def compute_umap(joint: ad.AnnData) -> ad.AnnData:
    joint.obsm["X_emb"] = joint.X if isinstance(joint.X, np.ndarray) else joint.X.toarray()
    sc.pp.neighbors(joint, use_rep="X_emb", n_neighbors=15, random_state=SEED)
    sc.tl.umap(joint, random_state=SEED)
    return joint


def plot_umap(joint: ad.AnnData, title: str, out_path: Path) -> None:
    n_sources = joint.obs[BATCH_KEY].nunique()
    palette_source = sc.pl.palettes.default_20[:n_sources]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel 1: colored by source
    sc.pl.umap(
        joint, color=BATCH_KEY, ax=axes[0], show=False,
        title="Batch (source)", frameon=False, palette=palette_source,
        size=3, alpha=0.5,
    )
    # Panel 2: colored by ref/query
    sc.pl.umap(
        joint, color="dataset", ax=axes[1], show=False,
        title="Reference vs Query", frameon=False,
        palette={"reference": "#aaaaaa", "query": "#e63946"},
        size=3, alpha=0.5,
    )

    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


def main():
    # All 3 arms × all 3 methods for source_5
    arms = ["tplus", "tminus", "tminus_matched"]
    queries = ["source_5"]
    paradigms = ["symphony", "scvi", "scpoli"]

    for query in queries:
        for paradigm in paradigms:
            for arm in arms:
                tag = f"{query}_{paradigm}_{arm}"
                out_path = UMAP_OUT / f"{tag}.png"
                if out_path.exists():
                    print(f"[{tag}] already exists — skip")
                    continue
                print(f"[{tag}] loading joint embedding ...")
                joint = load_joint(query, paradigm, arm)
                if joint is None:
                    continue
                print(f"[{tag}] computing UMAP on {len(joint)} cells ...")
                joint = compute_umap(joint)
                title = f"{query} → {paradigm} atlas | arm={arm}"
                plot_umap(joint, title, out_path)
                print(f"[{tag}] done")


if __name__ == "__main__":
    main()
