"""Compile all per-row metric CSVs into summary tables and bar charts.

Run after all nb40 runs complete:
  pixi run -e symphony python scripts/generate_report.py
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.paths import RESULTS_OUT

ACUTE = RESULTS_OUT / "acute"
METRIC_COLS = ["batch_asw", "kbet", "ilisi", "precision_at_10", "compound_asw"]
METRIC_LABELS = {
    "batch_asw":       "Batch ASW",
    "kbet":            "kBET",
    "ilisi":           "iLISI",
    "precision_at_10": "Precision@10",
    "compound_asw":    "Compound ASW",
}
PARADIGM_ORDER = ["symphony", "scvi", "scpoli"]
ARM_ORDER = ["tplus", "tminus_matched", "tminus"]
ARM_LABELS = {"tplus": "T+", "tminus_matched": "T−_matched", "tminus": "T−"}


def load_metrics(query_source: str) -> pd.DataFrame | None:
    csv = ACUTE / f"{query_source}_metrics.csv"
    if not csv.exists():
        print(f"  missing: {csv}")
        return None
    df = pd.read_csv(csv)
    # De-duplicate: keep last row per (paradigm, t_arm) combo in case of reruns
    df = df.drop_duplicates(subset=["paradigm", "t_arm"], keep="last").reset_index(drop=True)
    return df


def make_pivot(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        for col in METRIC_COLS:
            if col in row and pd.notna(row[col]):
                rows.append({
                    "paradigm": row["paradigm"],
                    "arm": row["t_arm"],
                    "metric": col,
                    "value": row[col],
                })
    long = pd.DataFrame(rows)
    if long.empty:
        return long
    pivot = long.pivot_table(index=["paradigm", "arm"], columns="metric", values="value")
    pivot = pivot.reindex(
        pd.MultiIndex.from_product([PARADIGM_ORDER, ARM_ORDER], names=["paradigm", "arm"]),
        fill_value=np.nan,
    ).dropna(how="all")
    pivot = pivot[[c for c in METRIC_COLS if c in pivot.columns]]
    return pivot


def plot_metrics(df: pd.DataFrame, query_source: str, out_path: Path) -> None:
    present_metrics = [m for m in METRIC_COLS if m in df.columns and df[m].notna().any()]
    n = len(present_metrics)
    if n == 0:
        print("  no metrics to plot")
        return

    fig, axes = plt.subplots(1, n, figsize=(4 * n, 4), sharey=False)
    if n == 1:
        axes = [axes]

    df_filt = df[
        df["paradigm"].isin(PARADIGM_ORDER) & df["t_arm"].isin(ARM_ORDER)
    ].copy()
    df_filt["arm_label"] = df_filt["t_arm"].map(ARM_LABELS)
    arm_label_order = [ARM_LABELS[a] for a in ARM_ORDER]

    for ax, metric in zip(axes, present_metrics):
        sns.barplot(
            data=df_filt,
            x="paradigm",
            y=metric,
            hue="arm_label",
            order=PARADIGM_ORDER,
            hue_order=arm_label_order,
            ax=ax,
            errwidth=1.2,
        )
        ax.set_title(METRIC_LABELS[metric])
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=20)
        if ax is not axes[0]:
            ax.get_legend().remove()
        else:
            ax.legend(title="Query arm", loc="best", frameon=False, fontsize=8)

    exp_label = "Exp1 (LOSO source_8)" if query_source == "source_8" else "Exp2 (wave2 source_5)"
    fig.suptitle(f"Reference mapping metrics — {exp_label}", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  wrote {out_path}")


def main():
    for query_source in ["source_5"]:
        print(f"\n=== {query_source} ===")
        df = load_metrics(query_source)
        if df is None or df.empty:
            continue

        # Pivot summary table
        pivot = make_pivot(df)
        csv_out = RESULTS_OUT / f"{query_source}_summary.csv"
        pivot.to_csv(csv_out)
        print(f"  pivot table ({pivot.shape}) → {csv_out}")
        print(pivot.to_string())

        # Bar chart
        plot_metrics(df, query_source, RESULTS_OUT / f"{query_source}_barplot.pdf")



if __name__ == "__main__":
    main()
