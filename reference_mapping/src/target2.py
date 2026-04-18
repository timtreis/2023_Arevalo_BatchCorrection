"""Derive the JUMP-Target2 compound manifest from pipeline metadata.

JUMP-Target2 is the consortium-defined set of ~302 compounds present in
every source by design. Used as the matched-anchor set for cross-source
reference mapping experiments.

Derivation strategy: read the metadata embedded in the corrected parquet
(or any pipeline parquet — they share the same `Metadata_*` columns), filter
to plate type TARGET2, and intersect compound sets across sources so we
end up with compounds that are truly shared everywhere.

Column names below are best-guess and need verification against the actual
parquet on first run. Override by passing explicit kwargs.
"""
from __future__ import annotations
import json
from pathlib import Path
import pandas as pd


DEFAULT_PLATE_COL = "Metadata_PlateType"
DEFAULT_SOURCE_COL = "Metadata_Source"
DEFAULT_COMPOUND_COL = "Metadata_JCP2022"  # JUMP CP 2022 compound ID
DEFAULT_TARGET2_VALUE = "TARGET2"


def derive_target2_manifest(
    metadata_parquet: str | Path,
    plate_col: str = DEFAULT_PLATE_COL,
    source_col: str = DEFAULT_SOURCE_COL,
    compound_col: str = DEFAULT_COMPOUND_COL,
    target2_value: str = DEFAULT_TARGET2_VALUE,
    columns: list[str] | None = None,
) -> set[str]:
    """Return the set of compound IDs shared by every source's TARGET2 plates.

    Parameters
    ----------
    metadata_parquet
        Path to any pipeline parquet — only the Metadata_* columns are read.
    plate_col, source_col, compound_col
        Override if the parquet uses different column names. The function
        validates that all three columns exist and raises with a helpful
        listing of available `Metadata_*` columns if not.
    target2_value
        The value of `plate_col` that identifies TARGET2 plates.
    columns
        Optional explicit list of columns to read. If None, reads only
        the three required columns.
    """
    metadata_parquet = Path(metadata_parquet)
    cols = columns or [plate_col, source_col, compound_col]

    try:
        meta = pd.read_parquet(metadata_parquet, columns=cols)
    except (KeyError, ValueError) as exc:
        head = pd.read_parquet(metadata_parquet).head(0)
        avail = [c for c in head.columns if c.startswith("Metadata_")]
        raise KeyError(
            f"Could not read columns {cols} from {metadata_parquet}. "
            f"Available Metadata_* columns: {avail}"
        ) from exc

    target2_plates = meta[meta[plate_col] == target2_value]
    if target2_plates.empty:
        plate_values = meta[plate_col].unique().tolist()
        raise ValueError(
            f"No rows with {plate_col} == {target2_value!r}. "
            f"Distinct plate types present: {plate_values}"
        )

    by_source = target2_plates.groupby(source_col)[compound_col].apply(set)
    target2_compounds = set.intersection(*by_source.tolist())

    if len(target2_compounds) < 250:
        raise ValueError(
            f"Expected ~302 Target2 compounds, got {len(target2_compounds)}. "
            "Verify column choice and source coverage."
        )

    return target2_compounds


def cache_manifest(compounds: set[str], cache_path: str | Path) -> None:
    """Persist a derived manifest as JSON for fast re-use."""
    cache_path = Path(cache_path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(sorted(compounds)))


def load_manifest(cache_path: str | Path) -> set[str]:
    return set(json.loads(Path(cache_path).read_text()))
