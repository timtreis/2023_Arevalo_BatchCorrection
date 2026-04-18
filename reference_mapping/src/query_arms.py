"""Construct the three query arms for a held-out source.

T+   : all compounds in the query retained (Target2 included).
T-   : Target2 compounds removed; query is smaller.
T-_matched : Target2 compounds removed AND a size-matched random sample of
             non-Target2 compounds is downsampled to keep the same compound
             count as T-. This isolates "presence of the matched-anchor
             compound *class*" from "raw query size", defusing the
             reviewer concern that any T+ vs T- difference is a sample-size
             effect (Option A in the implementation plan).

Convention: the compound id column is whatever was used to derive the
Target2 manifest — pass it explicitly via `compound_col`.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import anndata as ad


def make_tplus(adata: ad.AnnData) -> ad.AnnData:
    """All compounds retained; baseline arm."""
    return adata.copy()


def make_tminus(
    adata: ad.AnnData, target2_ids: set[str], compound_col: str
) -> ad.AnnData:
    """Drop all Target2 compounds. Query shrinks."""
    mask = ~adata.obs[compound_col].isin(target2_ids)
    return adata[mask].copy()


def make_tminus_matched(
    adata: ad.AnnData,
    target2_ids: set[str],
    compound_col: str,
    seed: int = 42,
    swap_log_path: str | Path | None = None,
) -> ad.AnnData:
    """Size-matched negative control for the Target2-presence question.

    Strategy (Option A from the implementation plan):
      1. Drop Target2 compounds entirely.
      2. From the remaining non-Target2 compounds, randomly DESIGNATE
         |target2_ids| compounds as the "swap" set and drop them too.
      3. Result is a query that has the same compound count and roughly
         the same cell count as T- — i.e. the difference between T+ and
         T-_matched isolates the *content* of the Target2 set (a known
         shared pharmacological panel) from raw size shrinkage.

    The swap_log_path argument records the chosen swap-out IDs so the
    metrics notebook can mask them when computing fair comparisons.
    """
    non_target2 = adata[~adata.obs[compound_col].isin(target2_ids)]
    unique_non_target2 = np.array(non_target2.obs[compound_col].unique())

    n_swap = len(target2_ids)
    if len(unique_non_target2) < n_swap:
        raise ValueError(
            f"Not enough non-Target2 compounds ({len(unique_non_target2)}) "
            f"to size-match against Target2 ({n_swap})."
        )

    rng = np.random.default_rng(seed)
    swap_out = set(rng.choice(unique_non_target2, size=n_swap, replace=False))

    if swap_log_path is not None:
        Path(swap_log_path).parent.mkdir(parents=True, exist_ok=True)
        Path(swap_log_path).write_text(json.dumps(sorted(swap_out)))

    keep_mask = ~adata.obs[compound_col].isin(target2_ids | swap_out)
    return adata[keep_mask].copy()


def load_swap_ids(swap_log_path: str | Path) -> set[str]:
    return set(json.loads(Path(swap_log_path).read_text()))
