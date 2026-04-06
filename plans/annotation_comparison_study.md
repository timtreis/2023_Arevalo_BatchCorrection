# Annotation Database Comparison Study: Plan

**Date**: 2026-04-05
**Goal**: Provide a data-driven, reproducible justification for which annotation database(s) to use when evaluating batch correction quality on JUMP-CP Cell Painting data. Publishable as supplementary analysis + online resource.

---

## 1. Research Question

> Given the JUMP-CP compound library, which annotation database provides the most reliable and informative biological labels for evaluating batch correction?

"Reliable" = high confidence, inter-database agreement, stable over time.
"Informative" = labels that actually predict morphological similarity in Cell Painting data.

The key insight from the literature: **coverage alone is insufficient**. A database with 62% coverage but noisy labels may be worse than one with 10% coverage but highly accurate labels. We need to measure both.

---

## 2. What Are We Evaluating?

### 2.1 Label Types (from specific to broad)

| Level | Example | Expected morphological coherence | Available in |
|-------|---------|--------------------------------|--------------|
| **Gene target** | EGFR, BRAF, HDAC1 | Low-moderate (~50% of genes produce detectable signatures; Rohban 2017) | ChEMBL, DRH, DrugBank, RefChemDB |
| **Target family** | RTKs, HDACs, tubulin | Moderate-high (natural grouping level for Cell Painting; Pahl & Ziegler 2023) | ChEMBL (target_type), UniProt keywords |
| **MOA** | kinase inhibitor, HDAC inhibitor, tubulin destabilizer | Highest (10 morphological subprofiles identified; Pahl & Ziegler 2023) | DRH, ChEMBL DRUG_MECHANISM |
| **Therapeutic area** | oncology, anti-inflammatory | Too broad for morphological grouping | DRH, DrugBank |

**For this study**: We evaluate at **gene target** and **MOA** levels. Target family is derived from gene target via UniProt protein family classification (computed, not from a single database).

### 2.2 The Ultimate Test

A database's labels are "good" for our purpose if **compounds sharing a label have similar morphological profiles after batch correction**. This is directly measurable via mAP. The database whose labels yield the highest mAP on well-corrected data provides the most morphologically informative grouping.

But there's a circularity: we need good labels to evaluate batch correction, and we need good batch correction to evaluate labels. **Solution**: Use a well-understood scenario (S1 or S2, single source or simple multi-source, where batch correction quality is established) as the test bed for comparing annotation databases.

---

## 3. Objective Evaluation Criteria

### 3.1 Coverage (Quantity)

| Metric | What it measures | How to compute |
|--------|-----------------|----------------|
| **Raw coverage** | % of JUMP compounds with >= 1 annotation | `n_mapped / n_total` |
| **Evaluable coverage** | % of compounds in groups with >= 3 members | Filter to groups with >= 3, count compounds |
| **Coverage by plate type** | Separate for TARGET2 (302) vs COMPOUND (~116K) | Stratify |
| **Coverage depth** | Distribution of annotations per compound | Histogram of n_labels per compound |

### 3.2 Confidence (Quality)

| Metric | What it measures | How to compute |
|--------|-----------------|----------------|
| **Confidence-stratified coverage** | Coverage at different quality thresholds | ChEMBL: confidence=9 vs >=7 vs any. RefChemDB: support>=5 vs >=2 vs >=1. DRH: all curated (no tiers). |
| **Evidence depth** | How much supporting evidence per annotation | ChEMBL: n_assays per compound-target pair. RefChemDB: support count. |
| **Primary vs all targets** | Fraction annotated with a single clear primary target vs polypharmacology | Count compounds with 1 vs 2+ targets per database |

### 3.3 Inter-Database Agreement (Consistency)

| Metric | What it measures | How to compute |
|--------|-----------------|----------------|
| **Pairwise concordance** | For compounds in both DB_A and DB_B, % agreeing on primary target | Jaccard similarity of target sets per compound, averaged |
| **N-way consensus** | Compounds where >= 3 databases agree | Intersection analysis |
| **Discordance rate** | Compounds where databases actively disagree (different primary targets) | Flag conflicts |

Key reference: Southan et al. 2013 found only 11.5% three-way consensus for approved drug targets across ChEMBL, DrugBank, TTD. We should expect similar.

### 3.4 Morphological Predictivity (The Ground Truth Test)

| Metric | What it measures | How to compute |
|--------|-----------------|----------------|
| **Label-stratified mAP** | Do compounds with the same label cluster morphologically? | copairs mAP using each database's labels as `pos_sameby` |
| **mAP by annotation level** | Which granularity best predicts morphology? | Compare mAP at gene target vs target family vs MOA |
| **mAP by confidence tier** | Do high-confidence annotations predict better? | Compare mAP using only high-confidence vs all annotations |
| **Label enrichment** | Fraction of label groups with mAP > random | Count groups with p < 0.05 |

This is the definitive criterion: **the database whose labels best predict morphological similarity wins**.

### 3.5 Practical Considerations

| Metric | What it measures |
|--------|-----------------|
| **Identifier bridging loss** | % of compounds lost during ID mapping (InChIKey -> DB-specific ID) |
| **License compatibility** | Open vs restricted |
| **Update frequency** | Last update date; is it maintained? |
| **Multi-label handling** | Fraction of compounds with multi-label annotations; impact on mAP |

---

## 4. Databases to Compare

| Database | Short name | Identifier | Confidence metric | Already in codebase |
|----------|-----------|------------|-------------------|-------------------|
| Drug Repurposing Hub | `drh` | InChIKey (direct) | None (all curated) | Yes |
| ChEMBL DRUG_MECHANISM | `chembl_curated` | InChIKey (direct) | max_phase, action_type | Partial |
| ChEMBL bioactivity | `chembl_bioactivity` | InChIKey (direct) | confidence score (0-9), pChEMBL | Partial |
| Open Targets | `opentargets` | ChEMBL ID -> InChIKey | Association score (0-1) | Yes (sparse) |
| RefChemDB | `refchemdb` | DTXSID (2-hop bridge) | Support count | Evaluated, rejected |
| DrugBank | `drugbank` | InChIKey | pharmacological_action flag | No |

---

## 5. Study Design

### 5.1 Central Compound Registry

Start from the JUMP-CP compound list as the single source of truth:

```
inputs/metadata/compound.csv.gz
  -> Metadata_JCP2022 (compound ID)
  -> Metadata_InChIKey (chemical identifier)
  -> ~116K compounds total
  -> 302 TARGET2 + ~116K COMPOUND
```

Map this INTO each database's space. Never the other way around. This ensures:
- Every database is evaluated on the same compound set
- Coverage is directly comparable
- No database-specific biases in compound selection

### 5.2 Notebook Structure

```
exploration/annotation_comparison/
    00_compound_registry.ipynb          # Load JUMP compounds, create master list
    01_map_drh.ipynb                    # Map to DRH, extract labels at all levels
    02_map_chembl_curated.ipynb         # Map to ChEMBL DRUG_MECHANISM
    03_map_chembl_bioactivity.ipynb     # Map to ChEMBL bioactivity (with thresholds)
    04_map_opentargets.ipynb            # Map to Open Targets
    05_map_refchemdb.ipynb              # Map to RefChemDB (document the failure)
    06_map_drugbank.ipynb               # Map to DrugBank (if license allows)
    10_coverage_comparison.ipynb        # Side-by-side coverage at all thresholds
    11_agreement_analysis.ipynb         # Inter-database concordance
    12_morphological_predictivity.ipynb # mAP using each database's labels (THE key notebook)
    13_confidence_vs_map.ipynb          # Does higher confidence = better morphological prediction?
    14_summary_figures.ipynb            # Publication-ready figures
    README.md                           # Motivation, methods, key findings
```

Each mapping notebook (01-06) follows the same template:

```python
# 1. Load JUMP compound registry
compounds = load_jump_compounds()  # InChIKey -> JCP2022

# 2. Load database
db = load_database(name)

# 3. Map: InChIKey join (or bridge for RefChemDB)
mapped = compounds.merge(db, on="Metadata_InChIKey", how="inner")

# 4. Report coverage
report_coverage(mapped, compounds, by=["plate_type"])

# 5. Extract labels at each granularity
labels_target = extract_target_labels(mapped)      # Gene-level
labels_family = extract_family_labels(mapped)       # Target family
labels_moa = extract_moa_labels(mapped)             # MOA-level

# 6. Apply confidence filters (database-specific)
for threshold in confidence_thresholds:
    labels_filtered = filter_by_confidence(mapped, threshold)
    report_coverage(labels_filtered, compounds)

# 7. Save standardized output
save_annotations(mapped, f"exploration/annotation_comparison/cache/{db_name}_annotations.parquet")
```

Standardized output schema (same for ALL databases):
```
Metadata_InChIKey     | str   | Chemical identifier
Metadata_JCP2022      | str   | JUMP compound ID
source                | str   | Database name
label_target          | str   | Gene symbol (e.g., EGFR)
label_family          | str   | Protein family (e.g., RTK) — derived via UniProt
label_moa             | str   | MOA (e.g., EGFR inhibitor) — if available
confidence            | float | Normalized 0-1 confidence (database-specific mapping)
evidence_count        | int   | Number of independent evidence sources
is_primary            | bool  | Whether this is the primary/curated annotation
```

### 5.3 Confidence Normalization

To compare across databases, normalize each database's confidence metric to a common 0-1 scale:

| Database | Raw metric | Normalization to 0-1 |
|----------|-----------|---------------------|
| DRH | None (all curated) | 1.0 for all |
| ChEMBL curated | max_phase (0-4) | phase / 4 |
| ChEMBL bioactivity | pChEMBL value | min(pChEMBL / 9, 1.0) — 9 = 1 nM |
| Open Targets | association_score (0-1) | Use directly |
| RefChemDB | support count (1-N) | min(count / 10, 1.0) |
| DrugBank | pharmacological_action (yes/no) | 1.0 if yes, 0.5 if unknown, 0.0 if no |

### 5.4 The Morphological Predictivity Test

This is the core experiment. For each database × confidence threshold × label level:

```python
# Load a well-corrected embedding (e.g., S1 or S2, best method)
adata = load_scenario("scenario_1", method="scpoli")  # or best method

# For each database's labels:
for db_name, annotations in all_annotations.items():
    # Merge labels into adata
    adata_labeled = merge_annotations(adata, annotations)

    for label_level in ["label_target", "label_family", "label_moa"]:
        for conf_threshold in [0.0, 0.25, 0.5, 0.75, 1.0]:
            subset = adata_labeled[adata_labeled["confidence"] >= conf_threshold]

            # Compute mAP
            map_result = copairs_map(
                subset,
                pos_sameby=[label_level],
                neg_diffby=[label_level],
                batch_col="Metadata_Source",
            )

            results.append({
                "database": db_name,
                "label_level": label_level,
                "confidence_threshold": conf_threshold,
                "n_compounds": subset.n_obs,
                "n_groups": subset[label_level].nunique(),
                "mean_map": map_result["mean_map"],
                "frac_retrievable": map_result["frac_p"],
            })
```

**Key analysis**: Plot mAP vs confidence threshold for each database. The database where mAP increases most steeply with confidence has the most informative confidence metric. The database with the highest mAP at any threshold has the best labels for our purpose.

---

## 6. Expected Findings (Hypotheses)

1. **DRH will have the highest mAP per compound on TARGET2** because its curated MOA labels are the natural grouping for well-characterized drugs
2. **ChEMBL bioactivity will have lowest mAP at low thresholds but competitive mAP at pChEMBL >= 7** because unfiltered bioactivity is noisy but high-affinity targets are real
3. **MOA-level labels will outperform gene-level labels** because multiple compounds targeting different genes in the same pathway produce similar morphology (Pahl & Ziegler 2023)
4. **Inter-database agreement will be low (~15-30%)** for target-level but higher (~50-70%) for MOA-level, because MOA is more stable than specific target assignment
5. **The confidence-mAP curve will be steepest for ChEMBL** because it has the most granular confidence metric — filtering from noisy (all bioactivity) to clean (confidence=9, pChEMBL>=7) should dramatically improve predictivity
6. **For COMPOUND plates, only ChEMBL bioactivity will have enough coverage** to compute meaningful mAP

---

## 7. Paper Narrative

### The argument:

1. "We need biological labels to evaluate whether batch correction preserves functional signal"
2. "Multiple databases exist, each with different coverage, curation, and label types"
3. "Rather than choosing arbitrarily, we systematically compared N databases on M criteria"
4. "We introduce a key criterion: **morphological predictivity** — do the labels actually predict Cell Painting similarity?"
5. "We find that [DRH MOA / ChEMBL target family / ...] provides the best trade-off of coverage, confidence, and morphological predictivity"
6. "We use a two-tier strategy: [Tier 1] for primary evaluation, [Tier 2] for extended COMPOUND evaluation"

### Novel contribution:

No published study has systematically compared annotation databases specifically for Cell Painting morphological evaluation. The "morphological predictivity" criterion — measuring whether labels predict embedding similarity — is a direct, task-specific evaluation that goes beyond generic database comparisons (Southan 2013, Isigkeit 2022).

### Figures:

1. **UpSet plot**: Compound overlap across databases (how many compounds are annotated in 1, 2, 3+ databases)
2. **Coverage curves**: Evaluable compounds vs confidence threshold, per database (line plot, one line per DB)
3. **mAP comparison**: Bar chart of mAP per database × label level (the main result)
4. **Confidence-mAP curves**: mAP vs confidence threshold per database (shows value of confidence filtering)
5. **Agreement heatmap**: Pairwise concordance matrix between databases
6. **Venn/UpSet for label agreement**: For shared compounds, overlap in assigned labels

---

## 8. Online Resource Structure

```
github.com/[repo]/annotation-comparison/
    README.md                  # Motivation + key findings + how to reproduce
    notebooks/
        00-06_mapping.ipynb    # One per database
        10-14_analysis.ipynb   # Comparison analyses
    data/
        compound_registry.parquet      # JUMP compound master list
        drh_annotations.parquet        # Standardized format
        chembl_curated_annotations.parquet
        chembl_bioactivity_annotations.parquet
        opentargets_annotations.parquet
        refchemdb_annotations.parquet
        agreement_matrix.parquet       # Inter-database concordance
        morphological_predictivity.parquet  # mAP results
    figures/
        coverage_comparison.pdf
        map_comparison.pdf
        confidence_curves.pdf
        agreement_heatmap.pdf
```

All annotation files use the same standardized schema (Section 5.2). Researchers can swap in their own compound list and re-run the mapping notebooks.

---

## 9. Implementation Timeline

| Phase | Work | Est. effort |
|-------|------|-------------|
| **Phase 0** | Download/prepare all databases | 1 day |
| **Phase 1** | Notebooks 00-06: Map JUMP compounds into each database | 2-3 days |
| **Phase 2** | Notebook 10: Coverage comparison | 0.5 day |
| **Phase 3** | Notebook 11: Agreement analysis | 1 day |
| **Phase 4** | Notebook 12: Morphological predictivity (requires S1/S2 embeddings) | 1-2 days |
| **Phase 5** | Notebook 13-14: Confidence analysis + figures | 1 day |
| **Total** | | ~6-8 days |

---

## 10. References

1. Corsello SM et al. The Drug Repurposing Hub. *Nat Med* 23, 405-408 (2017).
2. Zdrazil B et al. The ChEMBL Database in 2023. *NAR* 52, D1180-D1192 (2024).
3. Southan C et al. Comparing ChEMBL, DrugBank, HMDB and TTD. *Mol Inf* 32, 881-897 (2013).
4. Isigkeit L et al. Consensus Compound/Bioactivity Dataset. *J Med Chem* 65, 7602-7614 (2022).
5. Pahl A & Ziegler S. Morphological Subprofile Analysis. *Cell Chem Biol* 30, 839-851 (2023).
6. Rohban MH et al. Systematic Morphological Profiling. *eLife* 6, e24060 (2017).
7. Chandrasekaran SN et al. JUMP-CP. *Nat Methods* 21, 1114-1121 (2024).
8. Pabon NA et al. Reference Chemical Workflows. *Toxicol Sci* 152, 323-339 (2019).
9. Santos R et al. Comprehensive Map of Molecular Drug Targets. *Nat Rev Drug Discov* 16, 19-34 (2017).
10. Paluczak A et al. Integrative Annotation Workflows. *Database* 2025, baaf081 (2025).
