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
#     display_name: RefMap scVI (GPU)
#     language: python
#     name: refmap-scvi
# ---

# %% [markdown]
# > **⚠ sysVI cannot do reference mapping with new (unseen) batches.**
# >
# > `SysVI.load_query_data` hard-codes `transfer_batch=False`, which raises
# > `ValueError` when any query batch is missing from the reference registry.
# > This is architectural: sysVI learns system-specific transformations at
# > train time and has no surgery path for new batches.
# >
# > This notebook (and `11_map_sysvi.ipynb`) is kept as **scaffolding only**.
# > For reference mapping with new batches use **scPoli** (`20_train_atlas_scpoli`)
# > or **Symphony** (`30_train_atlas_harmony_symphony`).
# >
# > sysVI is appropriate only for in-sample batch correction (all batches seen
# > at train time). See `tasks/lessons.md` → "SysVI cannot do reference mapping
# > with new (unseen) batch categories (2026-04-19)" for the full diagnosis.

# %% [markdown]
# # 10 — sysVI atlas
#
# Two paths:
#
# 1. **Pipeline path (preferred)**: re-run the patched `methods_sysvi` rule against
#    `scenario_<refscenario>` so the model artifact lands at
#    `outputs/<scenario>/mad_int_featselect_sysvi_model/`.
#    This notebook then symlinks it into `models/`.
#
# 2. **Notebook path (fallback)**: re-train inline using the HPO best params.
#
# sysVI uses a Gaussian likelihood (MSE reconstruction) — no non-negative shift needed.

# %%
from pathlib import Path
import sys

sys.path.insert(0, str(Path.cwd().parent))
from src.paths import MODEL_OUT, scenario_model_dir, scenario_input_parquet, scenario_optuna_csv

# %% [markdown]
# ## Parameters

# %%
QUERY_SOURCE = "source_8"
REFERENCE_SCENARIO = "scenario_5"
BATCH_KEY = "Metadata_Source"
SOURCE_COL = "Metadata_Source"
ATLAS_NAME = f"sysvi_atlas_excl_{QUERY_SOURCE}"

# %% [markdown]
# ## Path 1 — symlink the pipeline-trained model

# %%
src_model_dir = scenario_model_dir(REFERENCE_SCENARIO, "sysvi")
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
#
# Only runs if `dst` does not exist after Path 1.

# %%
if not dst.exists():
    import pandas as pd
    from scvi.external import SysVI
    from src.data_io import load_parquet_as_anndata, split_by_source

    adata = load_parquet_as_anndata(scenario_input_parquet(REFERENCE_SCENARIO))
    ref, _ = split_by_source(adata, SOURCE_COL, QUERY_SOURCE)
    print(f"reference: {ref.shape}")

    params = pd.read_csv(scenario_optuna_csv(REFERENCE_SCENARIO, "sysvi"))
    p = params.sort_values("total", ascending=False).iloc[0].to_dict()

    print("\nBest HPO params:")
    for k, v in p.items():
        if k.startswith("params_"):
            print(f"  {k}: {v}")

    SysVI.setup_anndata(ref, batch_key=BATCH_KEY)
    vae = SysVI(
        ref,
        n_hidden=int(p["params_n_hidden"]),
        n_latent=int(p["params_n_latent"]),
        n_layers=int(p["params_n_layers"]),
        dropout_rate=float(p["params_dropout_rate"]),
        prior=p["params_prior"],
        n_prior_components=int(p["params_n_prior_components"]),
    )
    vae.train(
        max_epochs=400,
        early_stopping=True,
        early_stopping_monitor="validation_loss",
        plan_kwargs={
            "kl_weight": float(p["params_kl_weight"]),
            "z_distance_cycle_weight": (
                float(p["params_z_distance_cycle_weight"])
                if p["params_use_z_distance_cycle_weight"]
                else 0.0
            ),
        },
    )
    dst.parent.mkdir(parents=True, exist_ok=True)
    vae.save(str(dst), overwrite=True, save_anndata=True)
    print(f"inline-trained sysVI atlas to {dst}")
