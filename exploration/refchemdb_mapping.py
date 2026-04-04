"""
RefChemDB ↔ JUMP-CP Compound Mapping (POC: TARGET2)

Maps JUMP-CP compounds (identified by InChIKey) into the RefChemDB universe
(identified by DTXSID) using PubChem as a bridge.

RefChemDB provides expert-curated chemical-target annotations with evidence
support scores from 16+ independent sources. Unlike the Drug Repurposing Hub's
broad MOA categories, RefChemDB provides granular target-level annotations
(gene symbol + mode: agonist/antagonist/inhibitor/etc.) with quantified evidence.

Pipeline:
  1. Load TARGET2 compound InChIKeys from JUMP-CP metadata
  2. Batch-query PubChem: InChIKey → CID → DTXSID
  3. Download RefChemDB (EPA, Judson et al. 2019)
  4. Join on DTXSID
  5. Analyse coverage, evidence quality, and annotation distribution
"""

import json
import time
import urllib.request
import urllib.error
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO = Path(__file__).resolve().parent.parent
INPUTS = REPO / "inputs" / "metadata"
CACHE_DIR = REPO / "exploration" / "cache"
CACHE_DIR.mkdir(exist_ok=True)

REFCHEMDB_URL = (
    "https://gaftp.epa.gov/Comptox/NCCT_Publication_Data/FriedmanPaul_K/"
    "CompTox-ToxCast_EDSPsteroidogenesis_prioritization/misc/"
    "refChemDB_step_3_2018-05-30.txt"
)
REFCHEMDB_CACHE = CACHE_DIR / "refchemdb_step3.tsv"
PUBCHEM_CACHE = CACHE_DIR / "pubchem_inchikey_to_dtxsid.csv"


# ===========================================================================
# Step 1: Load TARGET2 compounds
# ===========================================================================
def load_target2_inchikeys() -> pd.DataFrame:
    """Return unique TARGET2 compounds with JCP2022 and InChIKey."""
    compound = pd.read_csv(INPUTS / "compound.csv.gz")
    plate = pd.read_csv(INPUTS / "plate.csv.gz")
    well = pd.read_csv(INPUTS / "well.csv.gz")

    t2_plates = plate.loc[
        plate["Metadata_PlateType"] == "TARGET2", "Metadata_Plate"
    ]
    t2_wells = well[well["Metadata_Plate"].isin(t2_plates)]
    t2_jcp = t2_wells["Metadata_JCP2022"].unique()

    t2 = compound[compound["Metadata_JCP2022"].isin(t2_jcp)].drop_duplicates(
        "Metadata_InChIKey"
    )
    # Drop the DMSO placeholder
    t2 = t2[t2["Metadata_InChIKey"].notna()].reset_index(drop=True)
    return t2[["Metadata_JCP2022", "Metadata_InChIKey"]]


# ===========================================================================
# Step 2: PubChem bridge — InChIKey → CID → DTXSID
# ===========================================================================
def _pubchem_inchikey_to_cid_batch(
    inchikeys: list[str], batch_size: int = 100
) -> dict[str, int | None]:
    """Batch-query PubChem: InChIKey → CID via POST."""
    mapping = {}
    for i in range(0, len(inchikeys), batch_size):
        batch = inchikeys[i : i + batch_size]
        body = "inchikey=" + ",".join(batch)
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/cids/JSON"
        req = urllib.request.Request(
            url, data=body.encode(), method="POST",
            headers={"Content-Type": "application/x-www-form-urlencoded"},
        )
        try:
            resp = json.loads(urllib.request.urlopen(req, timeout=30).read())
            for cid in resp["IdentifierList"]["CID"]:
                # PubChem returns 0 for not-found in batch mode
                if cid != 0:
                    mapping[batch[len(mapping) - sum(1 for v in mapping.values() if v is None)]] = cid
        except Exception:
            # Fall back to individual queries for this batch
            for ik in batch:
                mapping[ik] = _pubchem_inchikey_to_cid_single(ik)
                time.sleep(0.25)
        time.sleep(0.5)  # rate limit
    return mapping


def _pubchem_inchikey_to_cid_single(inchikey: str) -> int | None:
    """Single InChIKey → CID lookup."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON"
    try:
        resp = json.loads(urllib.request.urlopen(url, timeout=10).read())
        return resp["IdentifierList"]["CID"][0]
    except Exception:
        return None


def _pubchem_cids_to_dtxsid(cids: list[int]) -> dict[int, str | None]:
    """Map CIDs to DTXSIDs via PubChem synonym lookup."""
    mapping = {}
    for cid in cids:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
        try:
            resp = json.loads(urllib.request.urlopen(url, timeout=10).read())
            syns = resp["InformationList"]["Information"][0]["Synonym"]
            dtxsids = [s for s in syns if s.startswith("DTXSID")]
            mapping[cid] = dtxsids[0] if dtxsids else None
        except Exception:
            mapping[cid] = None
        time.sleep(0.25)  # rate limit
    return mapping


def build_inchikey_to_dtxsid(inchikeys: list[str]) -> pd.DataFrame:
    """
    Map InChIKeys → DTXSIDs via PubChem (with caching).

    Returns DataFrame with columns: Metadata_InChIKey, pubchem_cid, dtxsid
    """
    if PUBCHEM_CACHE.exists():
        cached = pd.read_csv(PUBCHEM_CACHE)
        cached_iks = set(cached["Metadata_InChIKey"])
        remaining = [ik for ik in inchikeys if ik not in cached_iks]
        if not remaining:
            print(f"Using cached PubChem mapping ({len(cached)} compounds)")
            return cached[cached["Metadata_InChIKey"].isin(inchikeys)]
        print(f"Cache has {len(cached)} compounds, querying {len(remaining)} new ones")
    else:
        cached = pd.DataFrame(
            columns=["Metadata_InChIKey", "pubchem_cid", "dtxsid"]
        )
        remaining = list(inchikeys)

    if remaining:
        # Step 2a: InChIKey → CID (individual queries — batch POST unreliable for InChIKeys)
        print(f"Querying PubChem for {len(remaining)} InChIKeys → CIDs...")
        ik_to_cid = {}
        for i, ik in enumerate(remaining):
            ik_to_cid[ik] = _pubchem_inchikey_to_cid_single(ik)
            if (i + 1) % 50 == 0:
                print(f"  ... {i + 1}/{len(remaining)} CIDs resolved")
            time.sleep(0.25)

        found_cids = {ik: cid for ik, cid in ik_to_cid.items() if cid is not None}
        print(
            f"  CID resolution: {len(found_cids)}/{len(remaining)} "
            f"({len(remaining) - len(found_cids)} not found)"
        )

        # Step 2b: CID → DTXSID
        if found_cids:
            print(f"Querying PubChem for {len(found_cids)} CIDs → DTXSIDs...")
            cid_to_dtxsid = _pubchem_cids_to_dtxsid(list(found_cids.values()))
            print(
                f"  DTXSID resolution: "
                f"{sum(1 for v in cid_to_dtxsid.values() if v is not None)}"
                f"/{len(found_cids)}"
            )
        else:
            cid_to_dtxsid = {}

        # Build new rows
        new_rows = []
        for ik in remaining:
            cid = ik_to_cid.get(ik)
            dtxsid = cid_to_dtxsid.get(cid) if cid else None
            new_rows.append(
                {"Metadata_InChIKey": ik, "pubchem_cid": cid, "dtxsid": dtxsid}
            )
        new_df = pd.DataFrame(new_rows)

        # Merge with cache and save
        cached = pd.concat([cached, new_df], ignore_index=True)
        cached.to_csv(PUBCHEM_CACHE, index=False)
        print(f"Saved PubChem cache ({len(cached)} total compounds)")

    return cached[cached["Metadata_InChIKey"].isin(inchikeys)].reset_index(
        drop=True
    )


# ===========================================================================
# Step 3: Download and load RefChemDB
# ===========================================================================
def load_refchemdb() -> pd.DataFrame:
    """Download (or load cached) RefChemDB step-3 data."""
    if not REFCHEMDB_CACHE.exists():
        print(f"Downloading RefChemDB from EPA...")
        urllib.request.urlretrieve(REFCHEMDB_URL, REFCHEMDB_CACHE)
        print(f"Saved to {REFCHEMDB_CACHE}")

    df = pd.read_csv(REFCHEMDB_CACHE, sep="\t", engine="python")
    # The 'count' column is the evidence support score: number of independent
    # sources confirming the chemical-target interaction. Higher = more reliable.
    # count=0 entries are artifacts (pmid/count columns may be swapped in some rows).
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0).astype(int)
    print(
        f"RefChemDB: {len(df):,} records, "
        f"{df['dsstox_substance_id'].nunique():,} compounds, "
        f"{df['target'].nunique():,} targets"
    )
    return df


# ===========================================================================
# Step 4: Join and analyse
# ===========================================================================
def map_and_analyse(
    t2: pd.DataFrame,
    bridge: pd.DataFrame,
    refdb: pd.DataFrame,
    min_evidence: int = 2,
):
    """Join TARGET2 compounds to RefChemDB and print analysis.

    Parameters
    ----------
    min_evidence
        Minimum evidence support score (number of independent sources).
        RefChemDB records with count < min_evidence are excluded.
        Recommended: 2 (confirmed by ≥2 independent sources).
    """
    # Merge TARGET2 → bridge → RefChemDB
    t2_bridge = t2.merge(bridge, on="Metadata_InChIKey", how="left")

    print("\n" + "=" * 70)
    print("BRIDGE COVERAGE (InChIKey → PubChem → DTXSID)")
    print("=" * 70)
    n_total = len(t2_bridge)
    n_cid = t2_bridge["pubchem_cid"].notna().sum()
    n_dtxsid = t2_bridge["dtxsid"].notna().sum()
    print(f"  TARGET2 compounds:        {n_total}")
    print(f"  Found in PubChem (CID):   {n_cid} ({100*n_cid/n_total:.1f}%)")
    print(f"  Have DTXSID:              {n_dtxsid} ({100*n_dtxsid/n_total:.1f}%)")

    # Evidence filtering
    print(f"\n" + "=" * 70)
    print(f"EVIDENCE FILTERING (min_evidence={min_evidence})")
    print("=" * 70)
    print(f"\n  RefChemDB evidence support distribution:")
    for thresh in [1, 2, 3, 5, 10]:
        sub = refdb[refdb["count"] >= thresh]
        print(
            f"    count ≥ {thresh:>2d}: {len(sub):>7,} records, "
            f"{sub['dsstox_substance_id'].nunique():>5,} cpds, "
            f"{sub['target'].nunique():>5,} targets"
        )

    refdb_filtered = refdb[refdb["count"] >= min_evidence].copy()
    print(
        f"\n  Using count ≥ {min_evidence}: "
        f"{len(refdb_filtered):,} records, "
        f"{refdb_filtered['dsstox_substance_id'].nunique():,} compounds"
    )

    # Join — both unfiltered and filtered
    t2_ref_all = t2_bridge.merge(
        refdb, left_on="dtxsid", right_on="dsstox_substance_id", how="inner"
    )
    t2_ref = t2_bridge.merge(
        refdb_filtered, left_on="dtxsid", right_on="dsstox_substance_id", how="inner"
    )

    n_all = t2_bridge["dtxsid"].isin(refdb["dsstox_substance_id"]).sum()
    n_filt = t2_bridge["dtxsid"].isin(refdb_filtered["dsstox_substance_id"]).sum()
    print(f"\n  TARGET2 compounds in RefChemDB (any evidence):    {n_all} ({100*n_all/n_total:.1f}%)")
    print(f"  TARGET2 compounds in RefChemDB (count ≥ {min_evidence}):     {n_filt} ({100*n_filt/n_total:.1f}%)")
    print(f"  Annotation records (any):     {len(t2_ref_all):,}")
    print(f"  Annotation records (filtered): {len(t2_ref):,}")

    print("\n" + "=" * 70)
    print(f"REFCHEMDB ANNOTATION COVERAGE (count ≥ {min_evidence})")
    print("=" * 70)

    if len(t2_ref) == 0:
        print("  No records after filtering!")
        return t2_ref

    # Unique targets per compound
    targets_per_cpd = (
        t2_ref.groupby("Metadata_InChIKey")["target"].nunique().describe()
    )
    print(f"\n  Targets per compound:")
    print(f"    mean: {targets_per_cpd['mean']:.1f}")
    print(f"    median: {targets_per_cpd['50%']:.0f}")
    print(f"    max: {targets_per_cpd['max']:.0f}")

    # Mode distribution
    print(f"\n  Mode distribution:")
    mode_counts = t2_ref["mode"].value_counts()
    for mode, cnt in mode_counts.head(10).items():
        print(f"    {mode}: {cnt:,}")

    # Top targets
    print(f"\n  Top 20 targets (by compound count):")
    target_cpd_counts = (
        t2_ref.groupby("target")["Metadata_InChIKey"]
        .nunique()
        .sort_values(ascending=False)
    )
    for target, cnt in target_cpd_counts.head(20).items():
        # Also show evidence range for this target
        evidence = t2_ref[t2_ref["target"] == target]["count"]
        print(f"    {target:<12s}: {cnt:>3d} cpds  (evidence: {evidence.min()}-{evidence.max()})")

    print("\n" + "=" * 70)
    print("COMPARISON WITH DRUG REPURPOSING HUB")
    print("=" * 70)

    drh = pd.read_csv(INPUTS / "repurposinghub_inchi_to_moa_moaexploded.csv")
    t2_iks = set(t2["Metadata_InChIKey"])

    drh_t2 = drh[drh["Metadata_InChIKey"].isin(t2_iks)]
    drh_cpds = drh_t2["Metadata_InChIKey"].nunique()
    drh_moas = drh_t2["Metadata_DRH_MOA"].dropna().nunique()

    ref_cpds = t2_ref["Metadata_InChIKey"].nunique()
    ref_targets = t2_ref["target"].nunique()

    print(f"\n  {'':30s} {'DRH':>10s} {'RefChemDB':>12s}")
    print(f"  {'─' * 55}")
    print(f"  {'Compounds covered':30s} {drh_cpds:>10d} {ref_cpds:>12d}")
    print(f"  {'Coverage of 302 T2 cpds':30s} {100*drh_cpds/302:>9.1f}% {100*ref_cpds/302:>11.1f}%")
    print(f"  {'Annotation categories':30s} {drh_moas:>10d} {'MOAs':<5s} {ref_targets:>5d} {'targets'}")
    print(f"  {'Total records':30s} {len(drh_t2):>10d} {len(t2_ref):>12d}")
    print(f"  {'Evidence level':30s} {'curated':>10s} {'count≥'+str(min_evidence):>12s}")

    # Overlap analysis
    both = t2_iks & set(drh_t2["Metadata_InChIKey"]) & set(t2_ref["Metadata_InChIKey"])
    drh_only = set(drh_t2["Metadata_InChIKey"]) - set(t2_ref["Metadata_InChIKey"])
    ref_only = set(t2_ref["Metadata_InChIKey"]) - set(drh_t2["Metadata_InChIKey"])

    print(f"\n  Overlap:")
    print(f"    Both DRH + RefChemDB:    {len(both)}")
    print(f"    DRH only:                {len(drh_only)}")
    print(f"    RefChemDB only:          {len(ref_only)}")

    # Target-level grouping potential (for mAP evaluation)
    print("\n" + "=" * 70)
    print(f"TARGET-LEVEL GROUPING POTENTIAL (for mAP, count ≥ {min_evidence})")
    print("=" * 70)

    print("\n  RefChemDB target grouping:")
    for min_cpds in [2, 3, 5, 10]:
        usable_targets = target_cpd_counts[target_cpd_counts >= min_cpds]
        usable_cpds = (
            t2_ref[t2_ref["target"].isin(usable_targets.index)]
            ["Metadata_InChIKey"].nunique()
        )
        print(
            f"    ≥{min_cpds} cpds/target: "
            f"{len(usable_targets):>4d} targets, "
            f"{usable_cpds:>3d} compounds evaluable"
        )

    # Compare with DRH MOA grouping
    print("\n  DRH MOA grouping (for comparison):")
    drh_moa_counts = (
        drh_t2.dropna(subset=["Metadata_DRH_MOA"])
        .groupby("Metadata_DRH_MOA")["Metadata_InChIKey"]
        .nunique()
        .sort_values(ascending=False)
    )
    for min_cpds in [2, 3, 5, 10]:
        usable_moas = drh_moa_counts[drh_moa_counts >= min_cpds]
        usable_cpds = (
            drh_t2[drh_t2["Metadata_DRH_MOA"].isin(usable_moas.index)]
            ["Metadata_InChIKey"].nunique()
        )
        print(
            f"    ≥{min_cpds} cpds/MOA:    "
            f"{len(usable_moas):>4d} MOAs,    "
            f"{usable_cpds:>3d} compounds evaluable"
        )

    # Save the joined data for downstream use
    output_path = CACHE_DIR / "target2_refchemdb_annotations.csv"
    t2_ref.to_csv(output_path, index=False)
    print(f"\n  Saved joined data to {output_path}")

    # Also save a summary at different evidence thresholds
    print("\n" + "=" * 70)
    print("SENSITIVITY ANALYSIS: Coverage vs Evidence Threshold")
    print("=" * 70)
    for ev_thresh in [1, 2, 3, 5]:
        sub = t2_bridge.merge(
            refdb[refdb["count"] >= ev_thresh],
            left_on="dtxsid",
            right_on="dsstox_substance_id",
            how="inner",
        )
        n_cpds = sub["Metadata_InChIKey"].nunique()
        n_targets = sub["target"].nunique()
        # Evaluable: targets with ≥3 compounds
        tpc = sub.groupby("target")["Metadata_InChIKey"].nunique()
        eval_targets = (tpc >= 3).sum()
        eval_cpds = sub[sub["target"].isin(tpc[tpc >= 3].index)]["Metadata_InChIKey"].nunique()
        print(
            f"  count ≥ {ev_thresh}: {n_cpds:>3d} cpds, {n_targets:>4d} targets "
            f"→ {eval_targets:>3d} evaluable targets ({eval_cpds:>3d} cpds with ≥3/target)"
        )

    return t2_ref


# ===========================================================================
# Main
# ===========================================================================
def main():
    print("=" * 70)
    print("RefChemDB ↔ JUMP-CP Mapping POC (TARGET2 compounds)")
    print("=" * 70)

    # Step 1
    print("\n--- Step 1: Load TARGET2 compounds ---")
    t2 = load_target2_inchikeys()
    print(f"Loaded {len(t2)} TARGET2 compounds")

    # Step 2
    print("\n--- Step 2: Bridge InChIKey → DTXSID via PubChem ---")
    bridge = build_inchikey_to_dtxsid(t2["Metadata_InChIKey"].tolist())

    # Step 3
    print("\n--- Step 3: Load RefChemDB ---")
    refdb = load_refchemdb()

    # Step 4
    print("\n--- Step 4: Join and analyse ---")
    t2_ref = map_and_analyse(t2, bridge, refdb)

    return t2_ref


if __name__ == "__main__":
    main()
