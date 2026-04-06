# Annotation Database Comparison for JUMP-CP Cell Painting

Which compound annotation database provides the most reliable and informative biological labels for evaluating batch correction?

## Metrics

We evaluate each database on **5 objective criteria**:

### M1: Coverage
- **Raw coverage**: % of JUMP compounds with >= 1 annotation
- **Evaluable coverage**: % in groups with >= 3 compounds (required for mAP)
- **Stratified**: Separate for TARGET2 (302 compounds) vs COMPOUND (~116K)

### M2: Confidence
- **Confidence-stratified coverage**: Coverage at different quality thresholds
- **Evidence depth**: N independent sources per annotation
- Each database's confidence metric is normalized to 0-1 for comparison

### M3: Annotation Quality
- **Polypharmacology**: % of compounds with multiple targets/MOAs
- **Primary vs all**: Is there a clear primary annotation?
- **Group size distribution**: Are groups well-sized for mAP evaluation?

### M4: Morphological Predictivity (the ground truth test)
- **Label-mAP**: Do compounds with the same label cluster in Cell Painting embeddings?
- **mAP by label level**: Gene target vs MOA — which level predicts morphology best?
- **mAP by confidence tier**: Does filtering to high-confidence improve prediction?

### M5: Practical
- **Identifier bridging loss**: % of compounds lost during ID mapping
- **License**: Open vs restricted
- **Freshness**: Last update date

## Notebook Structure

| Notebook | Purpose |
|----------|---------|
| `00_compound_registry.ipynb` | Load JUMP compounds, create master list |
| `01_map_drh.ipynb` | Map to Drug Repurposing Hub |
| `02_map_chembl_curated.ipynb` | Map to ChEMBL DRUG_MECHANISM |
| `03_map_chembl_bioactivity.ipynb` | Map to ChEMBL bioactivity (with thresholds) |
| `04_map_opentargets.ipynb` | Map to Open Targets |
| `05_map_refchemdb.ipynb` | Map to RefChemDB |
| `10_coverage_comparison.ipynb` | Side-by-side coverage at all thresholds |
| `11_agreement_analysis.ipynb` | Inter-database concordance |
| `12_morphological_predictivity.ipynb` | mAP using each database's labels |
| `13_confidence_vs_map.ipynb` | Does higher confidence = better prediction? |
| `14_summary_figures.ipynb` | Publication-ready figures |

## Standardized Output Schema

Every mapping notebook produces `cache/{db}_annotations.parquet` with identical columns:

| Column | Type | Description |
|--------|------|-------------|
| `Metadata_InChIKey` | str | Chemical identifier |
| `Metadata_JCP2022` | str | JUMP compound ID |
| `source` | str | Database name |
| `label_target` | str | Gene symbol (e.g., EGFR) |
| `label_moa` | str | MOA (e.g., EGFR inhibitor) |
| `confidence` | float | Normalized 0-1 |
| `evidence_count` | int | N independent sources |
| `is_primary` | bool | Primary/curated annotation? |

Plus `cache/{db}_summary.parquet` with the M1-M5 metrics, and `cache/{db}_morphological_predictivity.parquet` with per-method mAP results.

## Confidence Normalization

| Database | Raw metric | Normalization |
|----------|-----------|---------------|
| DRH | clinical_phase | Launched=1.0, P3=0.75, P2=0.6, P1=0.5, Preclinical=0.25 |
| ChEMBL curated | max_phase (0-4) | phase / 4 |
| ChEMBL bioactivity | pChEMBL | min(pChEMBL / 9, 1.0) |
| Open Targets | association_score | Use directly (0-1) |
| RefChemDB | support count | min(count / 10, 1.0) |

## How to Reproduce

```bash
# Run from repo root
cd exploration/annotation_comparison

# 1. Create compound registry
jupyter execute 00_compound_registry.ipynb

# 2. Map each database (independent, can run in parallel)
jupyter execute 01_map_drh.ipynb
jupyter execute 02_map_chembl_curated.ipynb
# ... etc

# 3. Cross-database comparison (requires all mapping notebooks to have run)
jupyter execute 10_coverage_comparison.ipynb
jupyter execute 11_agreement_analysis.ipynb
jupyter execute 12_morphological_predictivity.ipynb
jupyter execute 14_summary_figures.ipynb
```
