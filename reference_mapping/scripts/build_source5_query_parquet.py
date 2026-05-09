"""Subset scenario_source5/mad_int.parquet to scenario_5's 998 feature set.

The scenario_source5 pipeline preprocesses source_5 plates in isolation
(own negcon stats, own variant-feature selection), then MAD→INT normalises.
The result has ~4400 raw features. We subset to the exact 998 features used
by the scenario_5 atlas so the query shares the reference feature space with
real values (no zero-padding).

Output: reference_mapping/data/source5_query_mad_int_s5feats.parquet
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
SCENARIO5_FEATSELECT = ROOT / "outputs" / "scenario_5" / "mad_int_featselect.parquet"
SOURCE5_MAD_INT = ROOT / "outputs" / "scenario_source5" / "mad_int.parquet"
OUT = ROOT / "reference_mapping" / "data" / "source5_query_mad_int_s5feats.parquet"

def main():
    print("Loading scenario_5 feature list …")
    s5 = pd.read_parquet(SCENARIO5_FEATSELECT, columns=["Metadata_Source"])
    ref_cols = pd.read_parquet(SCENARIO5_FEATSELECT).columns.tolist()
    meta_cols = [c for c in ref_cols if c.startswith("Metadata_")]
    feat_cols = [c for c in ref_cols if not c.startswith("Metadata_")]
    print(f"  scenario_5: {len(feat_cols)} feature columns, {len(meta_cols)} metadata columns")

    print("Loading scenario_source5/mad_int.parquet …")
    q = pd.read_parquet(SOURCE5_MAD_INT)
    q_feat_cols = [c for c in q.columns if not c.startswith("Metadata_")]
    q_meta_cols = [c for c in q.columns if c.startswith("Metadata_")]
    print(f"  source5 mad_int: {len(q_feat_cols)} feature columns, {len(q)} rows")

    missing = [f for f in feat_cols if f not in q.columns]
    present = [f for f in feat_cols if f in q.columns]
    print(f"  Overlap: {len(present)}/{len(feat_cols)} reference features present in source5")
    if missing:
        print(f"  Missing {len(missing)} features — will zero-pad:")
        for f in missing[:10]:
            print(f"    {f}")
        if len(missing) > 10:
            print(f"    ... and {len(missing)-10} more")

    # Build output with all reference features, zero for any missing
    out_df = q[q_meta_cols].copy()
    for f in feat_cols:
        if f in q.columns:
            out_df[f] = q[f].values
        else:
            out_df[f] = 0.0

    # Reorder to match reference column order
    out_df = out_df[meta_cols + feat_cols]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_parquet(OUT, index=False)
    print(f"\nWrote {OUT}  shape={out_df.shape}")
    print(f"Sources: {out_df['Metadata_Source'].unique().tolist()}")
    print(f"Plate types: {out_df['Metadata_PlateType'].unique().tolist() if 'Metadata_PlateType' in out_df.columns else 'N/A'}")


if __name__ == "__main__":
    main()
