# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
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
from src.paths import MODEL_OUT, scenario_model_dir, scenario_input_parquet, scenario_optuna_csv

# %% [markdown]
# ## Parameters

# %%
QUERY_SOURCE = "source_5"
REFERENCE_SCENARIO = "scenario_5"
ATLAS_NAME = f"scpoli_atlas_excl_{QUERY_SOURCE}"

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
    from src.data_io import load_parquet_as_anndata, split_by_source

    SOURCE_COL = "Metadata_Source"
    BATCH_KEY = "Metadata_Source"
    LABEL_KEY = "Metadata_JCP2022"  # compound id; verify

    adata = load_parquet_as_anndata(scenario_input_parquet(REFERENCE_SCENARIO))
    ref, _ = split_by_source(adata, SOURCE_COL, QUERY_SOURCE)

    params = pd.read_csv(scenario_optuna_csv(REFERENCE_SCENARIO, "scpoli"))
    p = params.sort_values("total", ascending=False).iloc[0].to_dict()
    hidden = [int(p[f"params_layer_{i}_size"]) for i in range(int(p["params_num_layers"]))]

    model = scPoli(
        adata=ref,
        condition_keys=[BATCH_KEY],
        cell_type_keys=LABEL_KEY,
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
    print(f"inline-trained scPoli atlas to {dst}")
