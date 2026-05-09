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

# %% tags=["parameters"]
QUERY_SOURCE = "source_8"

# %%
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

metrics = [
    ("batch_asw",       "Batch ASW"),
    ("kbet",            "kBET"),
    ("ilisi",           "iLISI"),
    ("precision_at_10", "Precision@10"),
    ("compound_asw",    "Compound ASW"),
]

# Filter to paradigms/arms that are present to avoid seaborn warnings on missing groups
df = df[df["paradigm"].isin(paradigm_order) & df["t_arm"].isin(arm_order)]

fig, axes = plt.subplots(1, len(metrics), figsize=(4 * len(metrics), 4))
for ax, (metric, title) in zip(axes, metrics):
    present = df[metric].notna().any()
    if not present:
        ax.set_visible(False)
        continue
    sns.barplot(
        data=df,
        x="paradigm",
        y=metric,
        hue="t_arm",
        order=paradigm_order,
        hue_order=arm_order,
        ax=ax,
        errwidth=1.2,
    )
    ax.set_title(title)
    ax.set_xlabel("")
    if ax is not axes[0]:
        ax.get_legend().remove()
    else:
        ax.legend(title="Query arm", loc="best", frameon=False, fontsize=8)

fig.suptitle(f"Reference mapping — query={QUERY_SOURCE}", fontsize=11)
fig.tight_layout()

out_path = RESULTS_OUT / "acute" / f"{QUERY_SOURCE}_plot.pdf"
fig.savefig(out_path, bbox_inches="tight")
print(f"wrote {out_path}")
