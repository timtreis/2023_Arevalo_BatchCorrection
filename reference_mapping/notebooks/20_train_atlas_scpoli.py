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
# # 20 — scPoli atlas
#
# Mirror of `10_train_atlas_scvi`: prefer the pipeline-trained model from
# `methods_scpoli`; fall back to inline training using the HPO best params
# if the pipeline artifact is missing.
#
# **Open question (plan §17 Q4)**: confirm scPoli was trained with Target2
# as a `condition_key` for prototype anchoring — if not, retraining with
# the right condition_key may matter for the T+ vs T-_matched contrast.

# %%
from pathlib import Path
import sys

sys.path.insert(0, str(Path.cwd().parent))
from src.paths import DATA_OUT, MODEL_OUT, scenario_model_dir, scenario_input_parquet, scenario_optuna_csv

# %% [markdown]
# ## Parameters

# %% tags=["parameters"]
QUERY_SOURCE = ""          # empty = use all sources; "source_8" = exclude that source
REFERENCE_SCENARIO = "scenario_5"
# Limits prototypes to ~292 T2 compounds instead of 30K+. Set False only for debugging.
T2_ONLY_PROTOTYPES = True
ATLAS_NAME = (
    f"scpoli_atlas_{REFERENCE_SCENARIO}_full"
    if not QUERY_SOURCE
    else f"scpoli_atlas_excl_{QUERY_SOURCE}"
)

# %% [markdown]
# ## Path 1 — symlink the pipeline-trained model

# %%
src_model_dir = scenario_model_dir(REFERENCE_SCENARIO, "scpoli")
dst = MODEL_OUT / ATLAS_NAME

if src_model_dir.exists():
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.symlink_to(src_model_dir)
    print(f"symlinked {dst} -> {src_model_dir}")
else:
    print(f"NO pipeline model at {src_model_dir} — falling through to inline training")

# %% [markdown]
# ## Path 2 — inline training (fallback)

# %%
if not dst.exists():
    from scarches.models.scpoli import scPoli
    import pandas as pd
    import json
    from src.data_io import load_parquet_as_anndata, split_by_source
    from src.target2 import load_manifest

    SOURCE_COL = "Metadata_Source"
    BATCH_KEY = "Metadata_Source"
    LABEL_KEY = "Metadata_JCP2022"

    adata = load_parquet_as_anndata(scenario_input_parquet(REFERENCE_SCENARIO))
    ref = adata if not QUERY_SOURCE else split_by_source(adata, SOURCE_COL, QUERY_SOURCE)[0]

    params = pd.read_csv(scenario_optuna_csv(REFERENCE_SCENARIO, "scpoli"))
    p = params.sort_values("total", ascending=False).iloc[0].to_dict()
    hidden = [int(p[f"params_layer_{i}_size"]) for i in range(int(p["params_num_layers"]))]

    # Limit prototypes to T2 compounds only to avoid 30K+ prototype explosion.
    # Non-T2 wells are grouped into a single "BACKGROUND" prototype.
    if T2_ONLY_PROTOTYPES:
        t2_ids = load_manifest(DATA_OUT / "target2_compounds.json")
        ref.obs["_cell_type"] = ref.obs[LABEL_KEY].astype(str).where(
            ref.obs[LABEL_KEY].isin(t2_ids), "BACKGROUND"
        )
        effective_cell_type_keys = "_cell_type"
        print(f"T2-only prototypes: {len(t2_ids)} T2 + 1 BACKGROUND group")
    else:
        effective_cell_type_keys = LABEL_KEY
        print(f"All-compound prototypes: {ref.obs[LABEL_KEY].nunique()} unique — may CPU-hang")

    model = scPoli(
        adata=ref,
        condition_keys=[BATCH_KEY],
        cell_type_keys=effective_cell_type_keys,
        hidden_layer_sizes=hidden,
        latent_dim=int(p["params_latent_dim"]),
        embedding_dims=int(p["params_embedding_dims"]),
        recon_loss="mse",
    )
    model.train(
        n_epochs=600,
        pretraining_epochs=int(400 * float(p["params_pretrain_to_train_ratio"])),
        use_early_stopping=True,
        reload_best=True,
        alpha_epoch_anneal=int(p["params_alpha_epoch_anneal"]),
        eta=float(p["params_eta"]),
    )
    dst.mkdir(parents=True, exist_ok=True)
    model.save(str(dst), overwrite=True)

    # Save training config so nb21 can read it for labeled_indices logic
    (dst / "training_config.json").write_text(json.dumps({
        "t2_only_prototypes": T2_ONLY_PROTOTYPES,
        "effective_cell_type_col": effective_cell_type_keys,
    }))
    print(f"inline-trained scPoli atlas to {dst}")

    # Save reference embedding so nb40 can build joint ref+query metrics
    import anndata as ad
    import numpy as np

    # Save var_names so nb21 can zero-pad the query for missing reference features
    np.savetxt(dst / "var_names.txt", ref.var_names, fmt="%s")
    print(f"saved {len(ref.var_names)} var_names to {dst / 'var_names.txt'}")
    model.model.eval()
    ref_emb = model.get_latent(ref, mean=True)
    ref_out = ad.AnnData(obs=ref.obs[[BATCH_KEY, LABEL_KEY]])
    ref_out.obsm["X_scpoli_mapped"] = ref_emb
    ref_out.write_h5ad(dst / "reference.h5ad")
    print(f"saved reference embedding ({ref_emb.shape}) to {dst / 'reference.h5ad'}")
