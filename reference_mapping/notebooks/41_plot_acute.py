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
# # 41 — Plot acute results
#
# Reads `results/acute/<query>_metrics.csv` and produces the two-panel
# grouped bar chart used as the acute figure.

# %%
from pathlib import Path
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path.cwd().parent))
from src.paths import RESULTS_OUT

# %%
QUERY_SOURCE = "source_5"
csv_path = RESULTS_OUT / "acute" / f"{QUERY_SOURCE}_metrics.csv"
df = pd.read_csv(csv_path)
df

# %% [markdown]
# ## Two-panel grouped bar
#
# Order paradigms by ascending complexity, arms by tplus → tminus_matched → tminus
# so the eye reads "with anchors" → "matched control" → "no anchors".

# %%
paradigm_order = ["symphony", "scvi", "scpoli"]
arm_order = ["tplus", "tminus_matched", "tminus"]

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
for ax, metric, title in zip(
    axes,
    ["compound_map", "moa_map"],
    ["Compound mAP", "MOA mAP"],
):
    sns.barplot(
        data=df,
        x="paradigm",
        y=metric,
        hue="t_arm",
        order=paradigm_order,
        hue_order=arm_order,
        ax=ax,
    )
    ax.set_title(title)
    ax.set_xlabel("")
    ax.legend(title="Query arm", loc="best", frameon=False)

fig.suptitle(f"Reference mapping — query={QUERY_SOURCE}")
fig.tight_layout()

out_path = RESULTS_OUT / "acute" / f"{QUERY_SOURCE}_plot.pdf"
fig.savefig(out_path)
print(f"wrote {out_path}")
