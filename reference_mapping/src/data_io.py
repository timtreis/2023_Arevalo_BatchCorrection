"""Loaders that bridge the pipeline's parquet outputs and AnnData notebooks."""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import anndata as ad
import numpy as np


def load_parquet_as_anndata(
    parquet_path: str | Path,
    feature_prefix: str | None = None,
) -> ad.AnnData:
    """Load a pipeline parquet and split metadata vs feature columns.

    The pipeline parquets store both metadata (`Metadata_*` columns) and
    feature values (everything else, optionally namespaced via prefix like
    `scvi_*`, `scpoli_*`, `harmony_*`).

    feature_prefix
        If set, only columns starting with this prefix are treated as
        feature values. Otherwise we keep every non-Metadata column.
    """
    df = pd.read_parquet(parquet_path)
    meta_cols = [c for c in df.columns if c.startswith("Metadata_")]

    if feature_prefix is not None:
        feat_cols = [c for c in df.columns if c.startswith(feature_prefix)]
    else:
        feat_cols = [c for c in df.columns if c not in meta_cols]

    if not feat_cols:
        raise ValueError(
            f"No feature columns found in {parquet_path} "
            f"(prefix={feature_prefix!r}). Columns: {list(df.columns)[:10]}..."
        )

    X = df[feat_cols].to_numpy(dtype=np.float32)
    obs = df[meta_cols].copy().reset_index(drop=True)
    var = pd.DataFrame(index=feat_cols)

    return ad.AnnData(X=X, obs=obs, var=var)


def split_by_source(
    adata: ad.AnnData,
    source_col: str,
    query_source: str,
) -> tuple[ad.AnnData, ad.AnnData]:
    """Return (reference, query) where query == one source, reference == rest."""
    is_query = adata.obs[source_col] == query_source
    return adata[~is_query].copy(), adata[is_query].copy()
