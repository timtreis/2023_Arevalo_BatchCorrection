"""Lightweight Optuna utilities. No pandas/anndata/scib-metrics dependencies.

This module can be used in containers that only have Python + optuna (e.g., r.sif).
For the full benchmark utilities (scib_benchmark_embedding), use utils.py instead.
"""

import csv


def save_optuna_results(study, output_path: str) -> None:
    """Save Optuna study results in a standardized CSV format.

    Output columns: number, batch, bio, params_*, state, total
    This matches the format used by R optimization scripts.
    """
    trials = study.trials
    if not trials:
        return

    # Collect all param names across trials
    all_param_names = sorted({k for t in trials for k in t.params})

    header = ["number", "batch", "bio"]
    header += [f"params_{k}" for k in all_param_names]
    header += ["state", "total"]

    rows = []
    for t in trials:
        batch = t.values[0] if t.values else ""
        bio = t.values[1] if t.values and len(t.values) > 1 else ""
        total = 0.6 * bio + 0.4 * batch if isinstance(batch, (int, float)) and isinstance(bio, (int, float)) else ""
        row = [t.number, batch, bio]
        row += [t.params.get(k, "") for k in all_param_names]
        row += [t.state.name, total]
        rows.append(row)

    # Sort by total descending
    rows.sort(key=lambda r: r[-1] if isinstance(r[-1], (int, float)) else float("-inf"), reverse=True)

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)
