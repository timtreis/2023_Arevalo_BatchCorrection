# JUMP-CP Dataset Metadata Research

**Date**: 2026-03-29
**Purpose**: Document the experimental design of JUMP-CP to identify meaningful batch correction scenarios beyond the Arevalo et al. originals.

---

## 1. Data Generating Centers (Sources)

The JUMP-CP dataset was co-produced by 12 unique institutions across 14 source identifiers. Source_7 and source_13 are from the same institution. Source_4 is the Broad Institute (generates ORF/CRISPR data, minimal compound screening).

**Consortium partners**: Broad Institute (source_4), Amgen, AstraZeneca, Bayer, Biogen, Eisai, Janssen, Merck KGaA, Pfizer, Servier, Takeda, Ksilink.

The exact mapping of source number → institution is not publicly documented (except source_4 = Broad).

### Source roles and dataset codes (from Figure 1, Chandrasekaran et al.)

| Source | Wave | Primary dataset | Notes |
|--------|------|----------------|-------|
| source_1 | Wave 1 | cpg0016 (compound) | 1536-well plates |
| source_2 | Wave 1 | cpg0016 (compound) | |
| source_3 | Wave 1 | cpg0016 (compound) | |
| source_4 | Neither | cpg0000 (ORF), cpg0001/cpg0002 (compound), cpg0016 (CRISPR) | Broad Institute. Only 367 compounds. Main role is ORF/CRISPR. |
| source_5 | Wave 2 | cpg0016 (compound) | |
| source_6 | Wave 1 | cpg0016 (compound) | |
| source_7 | Bridge | cpg0016 (compound) | Overlaps both waves. 6,061 compounds. Same institution as source_13. |
| source_8 | Wave 1 | cpg0016 (compound) | |
| source_9 | Wave 2 | cpg0016 (compound) | 1536-well plates |
| source_10 | Wave 1 | cpg0016 (compound) | |
| source_11 | Wave 2 | cpg0016 (compound) | |
| source_13 | Neither | cpg0016 (CRISPR) | Same institution as source_7. No compounds. |
| source_15 | Wave 1 | cpg0016 (compound) | |

**Dataset codes:**
- **cpg0016** — main compound screening dataset (both waves) + CRISPR (source_13)
- **cpg0000** — ORF overexpression (source_4 + others)
- **cpg0001, cpg0002** — additional Broad compound datasets (source_4)
- **cpg0012** — Broad's own separate dataset

---

## 2. Microscope Systems

Five microscope systems were used, falling into two imaging modalities:

### Confocal (7 sources)

| Microscope | Manufacturer | Sources |
|------------|-------------|---------|
| CV8000 (CellVoyager) | Yokogawa | source_2, source_5, source_6, source_10 |
| CV7000 (CellVoyager) | Yokogawa | source_7, source_13 |
| ImageXpress Micro Confocal | Molecular Devices | source_8 |

### Widefield (6 sources)

| Microscope | Manufacturer | Sources | Notes |
|------------|-------------|---------|-------|
| Opera Phenix | PerkinElmer | source_1, source_3, source_4, source_9 | source_1 and source_9 use **1536-well** plates |
| Operetta | PerkinElmer | source_11 | |

### Microscope metadata in the pipeline

`preprocessing/io.py` contains a `MICRO_CONFIG` mapping and `add_microscopy_info()` function that adds `Metadata_Microscope` to all processed data. The microscope names are loaded from the JUMP datasets repository (`microscope_config.csv`).

### Microscope distribution in current scenarios

| Scenario | Sources | Microscopes used |
|----------|---------|-----------------|
| 1 | source_6 | CV8000 only |
| 2-3 | source_2, 6, 10 | CV8000 only (all confocal, same platform) |
| 4-5 | source_2, 3, 6, 8, 10 | CV8000 (3) + Opera Phenix (1) + ImageXpress (1) = 3 microscope types |
| 6 | source_2, 6 | CV8000 only |
| 7-8 | source_2, 3, 6, 8, 10 | Same as 4-5 |

**Key finding**: Scenarios 1-3 and 6 are all CV8000-only. Microscope heterogeneity only enters in scenarios 4-5 and 7-8. Sources not used in any scenario: 1, 4, 5, 7, 9, 11, 12, 13.

---

## 3. Well Plate Format

Most sources use **384-well plates** (16 rows × 24 columns), but two sources use **1536-well plates** (48 rows × 32 columns):

| Format | Layout | Sources |
|--------|--------|---------|
| 384-well | 16 rows × 24 columns (A01–P24) | source_2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13 |
| 1536-well | 48 rows × 32 columns (A01–AP32) | **source_1** (1472 wells used), **source_9** (full 1536) |

Both 1536-well sources use the Opera Phenix (widefield). This plate format difference is an additional technical confounder:
- Smaller wells → different cell seeding densities
- Different imaging field sizes / sites per well
- Potential edge effects differ between formats

**Neither source_1 nor source_9 is used in any current scenario (1-8).** This means plate format has not been tested as a batch factor yet — but it would be confounded with source if these sources are added.

---

## 4. Plate Types

| Plate type | Description | Usage |
|------------|-------------|-------|
| **TARGET2** | 306 positive control compounds, present in every batch across all sources. Primary anchor for cross-batch calibration. | Scenarios 1-8 |
| **COMPOUND** | Production plates with chemical compound perturbations (116,750 unique compounds total). | Scenarios 3, 5, 8 |
| **TARGET1** | Alternative positive control plate set | Available but rarely used |
| **ORF** | Open Reading Frame overexpression plates (12,602 genes) | Not in current scenarios |
| **CRISPR** | CRISPR knockout plates (7,975 genes) | Not in current scenarios |
| **COMPOUND_EMPTY** | Empty/control wells on compound plates | Quality control |

**Key finding**: Only TARGET2 and COMPOUND are used in current scenarios. ORF and CRISPR perturbation types represent entirely different modalities that could be explored.

---

## 4b. Compound Wave Structure

The COMPOUND plates are organized into **two production waves** with partial overlap between them. Two additional sources (source_4, source_7) are not part of either wave but share some compounds with both.

### Wave assignments

| Wave | Sources | Compounds per source | Replication |
|------|---------|---------------------|-------------|
| **Wave 1** | source_1, source_2, source_3, source_6, source_8, source_10, source_15 | ~11,000–12,700 each | Each partner exchanged with 4 others → compounds in up to 5 sites |
| **Wave 2** | source_5, source_9, source_11 | ~12,000 each | Each partner exchanged with the other 2 → compounds in up to 3 sites |
| **Neither** | source_4 (Broad Institute) | 367 | ORF/CRISPR focus, minimal compound screening |
| **Neither** | source_7 | 6,061 | Partial overlap with both waves (especially Wave 1 via source_3: 2,884 shared) |
| **Neither** | source_13 | 0 compounds | Same institution as source_7; generated CRISPR data |

### Pairwise compound overlap (from `compound_source.csv.gz`)

Most compounds appear in only 1 source (107,590 / 116,753 = 92%). The overlap structure:

```
Compounds by source count:
  1 source:  107,590
  2 sources:   8,209
  3 sources:     820
  4 sources:     119
  5 sources:      15
```

**Highest overlap pairs:**
- source_3 – source_7: 2,884 (source_7's primary overlap is with Wave 1)
- source_5 – source_8: 827 (cross-wave overlap)
- source_3 – source_9: 820 (cross-wave overlap)
- source_8 – source_11: 629 (cross-wave overlap)
- source_5 – source_6: 412 (cross-wave overlap)

**Within Wave 1 overlap** is relatively sparse (typically 50–400 shared compounds per pair). The design is NOT that all Wave 1 sources screen the same compounds — rather, each partner nominated their own ~12K and exchanged with 4 of the other 6 Wave 1 partners. So within-wave overlap is partial, not complete.

**Within Wave 2 overlap** is similarly sparse (55–629 per pair among sources 5, 9, 11).

**Cross-wave overlap** exists via bridge compounds scattered across pairs, typically 50–800 per cross-wave pair.

### Implications for batch correction and reference atlas design

1. **Compound replication is sparse**: Even within a wave, most compounds appear in only 1-2 sources. The "5-site replication" applies to compounds that were *exchanged*, not to each partner's nominated set.
2. **TARGET2 is the only universal anchor**: The 306 TARGET2 compounds are present in ALL sources regardless of wave. This is the primary evaluation tool for cross-source batch correction.
3. **Within-wave correction**: Still benefits from partial compound overlap (~100-400 shared compounds per pair within Wave 1) for mAP evaluation, beyond just TARGET2.
4. **Cross-wave correction**: Harder but not impossible — bridge compounds provide some cross-wave signal (50-800 per pair).
5. **Source_7 is a natural bridge**: Its heavy overlap with source_3 (2,884 compounds) makes it a potential cross-wave connector, despite not formally belonging to either wave.
6. **Reference atlas strategy**:
   - Building within Wave 1 (7 sources, 3 microscope types) gives the most diverse reference with reasonable internal overlap
   - Wave 2 (3 sources) can serve as a held-out query set
   - Source_7 can test bridging between waves
7. **Evaluation power varies by pair**: Some pairs share 2,884 compounds (great for mAP), others share only 2 (essentially TARGET2-only evaluation). Scenario design must account for this asymmetry.

---

## 5. Additional Metadata Columns

### Currently used in scenarios

- `Metadata_Source` — data generating center (batch_key in scenarios 2-8)
- `Metadata_Batch` — experimental batch within a source (batch_key in scenarios 1, 6-8)
- `Metadata_JCP2022` — compound/gene identifier (label_key, eval_key)
- `Metadata_DRH_MOA` — Mechanism of Action from Drug Repurposing Hub (eval_key in scenario 8)
- `Metadata_PlateType` — plate type category

### Available but not yet exploited

- `Metadata_Microscope` — microscope system name (added by preprocessing)
- `Metadata_Plate` — physical plate identifier
- `Metadata_Well` — well position (e.g., A01)
- `Metadata_Row` / `Metadata_Column` — row/column within plate
- `Metadata_PertType` — perturbation type (trt, poscon, negcon)
- `Metadata_InChIKey` / `Metadata_InChI` — chemical structure
- `Metadata_Symbol` — gene name (for ORF/CRISPR)
- `Metadata_NCBI_Gene_ID` — gene identifier
- `Anomaly` — quality flags (dye, segmentation, infection, none)
- `Density` — cell density (20%, 50%, 80%, 100%, 120%)
- `Polybrene` — transfection reagent (Present/Absent)

---

## 6. Analysis of Existing Scenarios

### Dimensions that vary across scenarios

| Dim | S1 | S2 | S3 | S4 | S5 | S6 | S7 | S8 | Apricot |
|-----|----|----|----|----|----|----|----|----|---------|
| # Sources | 1 | 3 | 3 | 5 | 5 | 2 | 5 | 5 | 11 |
| Sources | 6 | 2,6,10 | 2,6,10 | 2,3,6,8,10 | 2,3,6,8,10 | 2,6 | 2,3,6,8,10 | 2,3,6,8,10 | 1-8,10,11,13 |
| Plate types | T2 | T2 | T2+C | T2 | T2+C | T2 | T2 | T2+C | T2 |
| Batch key | Batch | Source | Source | Source | Source | Src+Batch | Src+Batch | Src+Batch | Src+Batch |
| Eval keys | JCP | JCP | JCP | JCP | JCP | JCP+MOA | JCP+MOA | JCP+MOA | JCP+MOA |
| Waves | W1 | W1 | W1 | W1 | W1 | W1 | W1 | W1 | W1+W2+bridge |
| Microscopes | CV8000 | CV8000 | CV8000 | 3 types | 3 types | CV8000 | 3 types | 3 types | 5 types |
| Plate format | 384 | 384 | 384 | 384 | 384 | 384 | 384 | 384 | 384+1536 |

### What the existing scenarios test

**Tested:**
- Intra-source batch correction (S1: batches within source_6)
- Same-microscope cross-source (S2-3: CV8000 only)
- Mixed-microscope cross-source (S4-5, S7-8: 3 microscope types)
- TARGET2-only vs COMPOUND+TARGET2 (S2→S3, S4→S5, S7→S8)
- Simple vs hierarchical batch key (S4→S7, S5→S8)
- MOA-based evaluation (S6-8)
- Small vs large scale (S1→S8)

**NOT tested by any existing scenario:**
- Widefield-only sources (all scenarios with >1 source use CV8000-dominated sets)
- Cross-wave correction (S1-8 are all Wave 1; only Apricot mixes waves)
- 1536-well plate sources (only Apricot includes S1/S9)
- Source_5 (CV8000, Wave 2 — never used despite being same microscope as S2,6,10)
- Source_7 as bridge between waves
- Source_11, source_15 in controlled settings (only in Apricot)
- Within-wave-2 correction
- Microscope type as explicit batch key
- ORF or CRISPR perturbation types

### Gaps in coverage

1. **Wave 2 is invisible**: No scenario tests within-Wave-2 or cross-wave correction in a controlled way
2. **Widefield is underrepresented**: source_3 is the only widefield source in S4-8, always mixed with confocal
3. **Same-microscope comparisons only exist for CV8000**: No Opera Phenix-only or CV7000-only scenarios
4. **Plate format is untested**: 1536-well sources (S1, S9) only appear in Apricot where they're confounded with everything else
5. **No microscope-type batch key**: Microscope type is always confounded with source — never isolated as the batch factor

---

## 7. Comprehensive List of Possible New Scenarios

### Notation
- **T2** = TARGET2, **C** = COMPOUND, **Src** = Metadata_Source, **Batch** = Metadata_Batch
- Microscope types: CV8 = CV8000, CV7 = CV7000, OP = Opera Phenix, IX = ImageXpress, OPR = Operetta

---

### A. Same-microscope scenarios (isolate lab effects from microscope effects)

| ID | Sources | Microscope | Wave | Plates | Batch key | Question |
|----|---------|-----------|------|--------|-----------|----------|
| A1 | 2, 5, 6, 10 | CV8000 (all) | W1+W2 | T2 | Src | Cross-wave, same microscope. Does wave matter when microscope is held constant? |
| A2 | 2, 6, 10 | CV8000 (all) | W1 only | T2+C | Src | Same as S3 but note: this IS scenario 3. Included for completeness. |
| A3 | 2, 5, 6, 10 | CV8000 (all) | W1+W2 | T2+C | Src | A1 + COMPOUND plates. Cross-wave with compound overlap for evaluation. |
| A4 | 1, 3, 4, 9 | Opera Phenix (all) | Mixed | T2 | Src | Widefield-only. Includes 384+1536 well plates. |
| A5 | 3, 4 | Opera Phenix (all) | W1+neither | T2 | Src | Widefield-only, 384-well only (exclude 1536). |
| A6 | 1, 3, 9 | Opera Phenix (all) | W1+W2 | T2 | Src | Widefield-only, cross-wave. Mixed plate format (S1,S9=1536, S3=384). |
| A7 | 7, 13 | CV7000 (all) | bridge+neither | T2 | Src | Same-institution, same microscope. Minimal scenario — do identical setups still have batch effects? Note: S13 has no compounds, only CRISPR. |
| A8 | 5, 9, 11 | CV8000+OP+OPR | W2 only | T2 | Src | Within-Wave-2 only. 3 different microscopes. |
| A9 | 5, 9, 11 | CV8000+OP+OPR | W2 only | T2+C | Src | A8 + COMPOUND plates. Exploits within-Wave-2 compound overlap. |

### B. Cross-microscope scenarios (isolate microscope effects)

| ID | Sources | Microscopes | Wave | Plates | Batch key | Question |
|----|---------|------------|------|--------|-----------|----------|
| B1 | 2, 3 | CV8000 + OP | W1 | T2 | Src | Minimal confocal vs widefield (2 sources). |
| B2 | 2, 6, 3, 15 | CV8000(2) + OP(2) | W1 | T2 | Src | Balanced confocal vs widefield, 2 per type. |
| B3 | 2, 6, 10, 3, 15 | CV8000(3) + OP(2) | W1 | T2 | Src | Slightly unbalanced. 5 sources, all Wave 1, 384-well. |
| B4 | 2, 6, 8 | CV8000(2) + IX(1) | W1 | T2 | Src | Cross-confocal-platform. Different manufacturers, same modality. |
| B5 | 2, 6, 10, 7 | CV8000(3) + CV7000(1) | W1+bridge | T2 | Src | Same manufacturer (Yokogawa), different generation. |
| B6 | 2, 6, 10, 8, 7 | CV8000(3)+IX(1)+CV7000(1) | W1+bridge | T2 | Src | All confocal, 3 different platforms. |
| B7 | 3, 1, 9 | OP(3) | W1+W2 | T2 | Src | Same microscope but mixed plate format (384+1536). Isolates plate format effect. |

### C. Wave-focused scenarios (isolate wave/compound-exchange effects)

| ID | Sources | Microscopes | Wave | Plates | Batch key | Question |
|----|---------|------------|------|--------|-----------|----------|
| C1 | All Wave 1: 1,2,3,6,8,10,15 | CV8000+OP+IX | W1 | T2 | Src | Full Wave 1. 7 sources, 3 microscope types, mixed plate format. |
| C2 | All Wave 1: 1,2,3,6,8,10,15 | CV8000+OP+IX | W1 | T2+C | Src | C1 + COMPOUND. Richest within-wave compound overlap for mAP. |
| C3 | W1(384 only): 2,3,6,8,10,15 | CV8000+OP+IX | W1 | T2 | Src | Full Wave 1 minus S1 (1536-well). Clean plate format. |
| C4 | W1(384 only): 2,3,6,8,10,15 | CV8000+OP+IX | W1 | T2+C | Src+Batch | C3 + compounds + hierarchical batch. The "clean Wave 1 stress test". |
| C5 | All Wave 2: 5, 9, 11 | CV8000+OP+OPR | W2 | T2 | Src | Same as A8. Within-Wave-2. |
| C6 | W1 subset + W2: 2,6,10 + 5,9,11 | CV8000+OP+OPR | W1+W2 | T2 | Src | Cross-wave. 6 sources. Compound overlap enables evaluation. |
| C7 | W1 subset + W2 + bridge: 2,6,10 + 5,9,11 + 7 | 4 types | all | T2 | Src | Cross-wave with bridge source. 7 sources. |
| C8 | W1 subset + bridge: 2,3,6,10 + 7 | CV8000+OP+CV7000 | W1+bridge | T2+C | Src | Test S7 as bridge. S3-S7 share 2,884 compounds. |

### D. Scale scenarios

| ID | Sources | Microscopes | Wave | Plates | Batch key | Question |
|----|---------|------------|------|--------|-----------|----------|
| D1 | All 384-well: 2,3,4,5,6,7,8,10,11,13,15 | 5 types | all | T2 | Src+Batch | Maximum diversity, controlled plate format. 11 sources. |
| D2 | All sources: 1,2,3,4,5,6,7,8,9,10,11,13,15 | 5 types | all | T2 | Src+Batch | True maximum. 13 sources. Mixed plate format. Essentially scenario_apricot + S9,S15. |
| D3 | All compound sources: 1,2,3,5,6,7,8,9,10,11,15 | 5 types | all | T2+C | Src+Batch | D2 minus non-compound sources (S4, S13). With compounds for richer eval. |

### E. Hierarchical batch key scenarios (re-runs of existing with different batch keys)

| ID | Base | Batch key change | Question |
|----|------|-----------------|----------|
| E1 | S2 sources (2,6,10) | Src+Batch (was just Src) | Does hierarchical batch key improve same-microscope correction? |
| E2 | S4 sources (2,3,6,8,10) | Src only (was Src+Batch in S7) | Is hierarchical batch key necessary, or does Src alone suffice? |
| E3 | Any multi-source | Metadata_Microscope | Can microscope type serve as the batch key? Tests whether microscope is the primary batch axis. |

### F. Reference mapping scenarios (train/test splits)

| ID | Train sources | Test source(s) | Microscope test | Question |
|----|--------------|----------------|-----------------|----------|
| F1 | 2, 6, 10 (CV8000, W1) | 8 (IX, W1) | Cross-platform confocal | Does same-wave, cross-microscope mapping work? |
| F2 | 2, 6, 10 (CV8000, W1) | 5 (CV8000, W2) | Same microscope, cross-wave | Does cross-wave mapping work when microscope is held constant? |
| F3 | 2, 6, 10 (CV8000, W1) | 3 (OP, W1) | Cross-modality, same wave | Confocal reference, widefield query. |
| F4 | 2, 3, 6, 8, 10 (W1, 384) | 5, 9, 11 (W2) | Cross-wave at scale | Full Wave 1 reference, Wave 2 query. |
| F5 | 2, 3, 6, 8, 10 (W1, 384) | 7 (bridge) | Bridge source mapping | Does S7's compound overlap help mapping? |
| F6 | 2, 3, 6, 10, 15 (W1, 384) | 1 (W1, 1536) | Within-wave, cross-plate-format | Same wave but different plate format. |

### G. Perturbation modality scenarios

| ID | Sources | Pert type | Plates | Question |
|----|---------|-----------|--------|----------|
| G1 | 4 + others with ORF | ORF | ORF plates | Do batch correction methods work for genetic perturbations? |
| G2 | 4, 13 + others | CRISPR | CRISPR plates | Same question for CRISPR knockouts. |
| G3 | 4 + compound sources | COMPOUND + ORF | Mixed | Cross-modality integration. Can methods align chemical and genetic perturbation spaces? |

**Note**: G1-G3 require pipeline changes to support ORF/CRISPR plate types and different label keys (gene symbols instead of JCP2022 compound IDs).

---

### Summary: scenario design dimensions

| Dimension | Values | Tested by existing | Not yet tested |
|-----------|--------|-------------------|----------------|
| Microscope homogeneity | Same / same-type / mixed | S1-3 (same), S4-8 (mixed) | Same-type widefield, CV7000-only |
| Wave | W1 / W2 / cross-wave | S1-8 (W1 only) | W2-only, cross-wave, bridge |
| Plate format | 384 / 1536 / mixed | S1-8 (384 only) | 1536-only, mixed |
| Plate types | T2 / T2+C / ORF / CRISPR | T2, T2+C | ORF, CRISPR |
| # Sources | 1-13 | 1,2,3,5 | 4+, widefield subset |
| Batch key | Src / Batch / Src+Batch | All three | Metadata_Microscope |
| Eval key | JCP / JCP+MOA | Both | Gene-level (for ORF/CRISPR) |
| Reference mapping | Train/test split | None | All of F1-F6 |

---

## 8. Scenario Assessment: Scientific Utility, Reference Mapping, and ML-Readiness

### Assessment criteria

Each scenario is scored on five axes (High / Medium / Low):

1. **Scientific utility**: Does the community need this? Does it answer a question nobody else has answered?
2. **Reference mapping suitability**: Would the corrected embedding be a good reference atlas for mapping new data into?
3. **Embedding quality (expected)**: Given the difficulty, how clean will the resulting embedding be? Easier corrections → higher quality.
4. **ML-readiness**: Is the embedding useful as features for downstream ML (compound activity prediction, MOA classification, hit calling)?
5. **Feasibility**: Can we run it with the current pipeline without major changes?

### Assessment key

- **Ref-map quality** depends on: (a) correction quality, (b) chemical diversity in the embedding, (c) microscope diversity (more = more generalizable but harder), (d) compound overlap for evaluation
- **ML-readiness** depends on: (a) embedding quality, (b) number of labeled compounds, (c) balanced representation across batches, (d) perturbation diversity
- **Practical note**: T2+C scenarios are always better for ML than T2-only (117K vs 306 compounds), but are much larger and harder to correct

---

### Compound replication analysis

A critical factor for both evaluation quality (mAP needs compounds in multiple sources) and ML-readiness (replicated compounds provide cross-validated features). Analysis of `compound_source.csv.gz`:

**Key finding: compound replication is extremely sparse.** No compound appears in all sources, all Wave 1 sources, or even all CV8000 sources. TARGET2 (306 compounds, present in ALL sources by design) is the only universal anchor.

#### Compound replication tally per scenario (COMPOUND plates only, excludes 306 TARGET2)

All scenarios additionally have 306 TARGET2 compounds shared across ALL sources.

| Scenario | Plates | #Src | in 1 | in 2 | in 3 | in 4 | in 5 | Total | % repl. |
|----------|--------|------|------|------|------|------|------|-------|---------|
| **Existing** | | | | | | | | | |
| S1 | T2 | 1 | 12,443 | - | - | - | - | 12,443 | 0% |
| S2 | T2 | 3 | 36,260 | 76 | - | - | - | 36,336 | 0.2% |
| S3 | T2+C | 3 | 36,260 | 76 | - | - | - | 36,336 | 0.2% |
| S4/S7 | T2 | 5 | 58,530 | 521 | 28 | - | - | 59,079 | 0.9% |
| S5/S8 | T2+C | 5 | 58,530 | 521 | 28 | - | - | 59,079 | 0.9% |
| S6 | T2 | 2 | 24,381 | 26 | - | - | - | 24,407 | 0.1% |
| Apricot | T2 | 11 | 86,560 | 6,751 | 550 | 84 | 7 | 93,952 | 7.9% |
| **A. Same-microscope** | | | | | | | | | |
| A1 | T2 | 4 | 46,726 | 843 | - | - | - | 47,569 | 1.8% |
| A3 | T2+C | 4 | 46,726 | 843 | - | - | - | 47,569 | 1.8% |
| A4 | T2 | 4 | 32,863 | 1,368 | 22 | - | - | 34,253 | 4.1% |
| A5 | T2 | 2 | 10,944 | 228 | - | - | - | 11,172 | 2.0% |
| A6 | T2 | 3 | 32,938 | 1,174 | 4 | - | - | 34,116 | 3.5% |
| A7 | T2 | 2 | 6,061 | - | - | - | - | 6,061 | 0% |
| A8/C5 | T2 | 3 | 35,781 | 241 | 1 | - | - | 36,023 | 0.7% |
| A9 | T2+C | 3 | 35,781 | 241 | 1 | - | - | 36,023 | 0.7% |
| **B. Cross-microscope** | | | | | | | | | |
| B1 | T2 | 2 | 22,491 | 266 | - | - | - | 22,757 | 1.2% |
| B2 | T2 | 4 | 47,508 | 314 | - | - | - | 47,822 | 0.7% |
| B3 | T2 | 5 | 58,823 | 604 | 28 | - | - | 59,455 | 1.1% |
| B4 | T2 | 3 | 36,576 | 34 | - | - | - | 36,610 | 0.1% |
| B5 | T2 | 4 | 41,498 | 444 | 29 | - | - | 41,971 | 1.1% |
| B6 | T2 | 5 | 53,659 | 469 | 29 | - | - | 54,157 | 0.9% |
| B7 | T2 | 3 | 32,938 | 1,174 | 4 | - | - | 34,116 | 3.5% |
| **C. Wave-focused** | | | | | | | | | |
| C1 | T2 | 7 | 82,639 | 796 | 30 | 1 | - | 83,466 | 1.0% |
| C2 | T2+C | 7 | 82,639 | 796 | 30 | 1 | - | 83,466 | 1.0% |
| C3 | T2 | 6 | 71,002 | 620 | 28 | - | - | 71,650 | 0.9% |
| C4 | T2+C | 6 | 71,002 | 620 | 28 | - | - | 71,650 | 0.9% |
| C6 | T2 | 6 | 67,834 | 2,305 | 78 | - | - | 70,217 | 3.4% |
| C7 | T2 | 7 | 72,769 | 2,727 | 164 | 6 | - | 75,666 | 3.8% |
| **C8** | **T2+C** | **5** | **46,835** | **2,798** | **321** | **28** | - | **49,982** | **6.3%** |
| **D. Scale** | | | | | | | | | |
| D1 | T2 | 11 | 88,075 | 6,366 | 535 | 78 | 6 | 95,060 | 7.3% |
| D2 | T2 | 13 | 107,590 | 8,209 | 820 | 119 | 15 | 116,753 | 7.9% |
| D3 | T2+C | 11 | 107,656 | 8,222 | 728 | 82 | 8 | 116,696 | 7.7% |
| **E. Batch key** | | | | | | | | | |
| E3 | T2 | 5 | 58,530 | 521 | 28 | - | - | 59,079 | 0.9% |

**Notes:**
- "in N" = number of COMPOUND-plate compounds appearing in exactly N sources within that scenario
- "% repl." = percentage of compounds in 2+ sources (evaluable via mAP beyond TARGET2)
- T2-only scenarios: evaluation uses 306 TARGET2 compounds regardless of COMPOUND replication
- T2+C scenarios: evaluation uses 306 TARGET2 + replicated COMPOUND compounds
- Bolded C8 has the best replication of any single-scenario design (6.3%), driven by S3↔S7 overlap

#### Actual well-level replicate counts per compound (from `well.csv.gz`)

The "in N sources" table above counts *source-level* presence. But within each source, compounds have multiple well-level replicates. This table shows actual replicate counts — the number of wells per compound across all sources in each scenario.

**TARGET2 compounds** (302 compounds, present in all sources by design):

| Scenario | # Sources | Median reps/cpd | Range |
|----------|----------|----------------|-------|
| S3 | 3 | 37 | 25–2,399 |
| S5/S8 | 5 | 66 | 54–4,255 |
| A3 | 4 | 61 | 49–3,935 |
| A9 (Wave 2) | 3 | 67 | 64–4,288 |
| C4 (W1 384) | 6 | 66 | 54–4,255 |
| C8 (W1+bridge) | 5 | 69 | 57–4,447 |
| D3 (all cpd) | 11 | 140 | 128–8,991 |

TARGET2 compounds are heavily replicated across all scenarios — this is by design. They always provide a well-powered evaluation anchor.

**COMPOUND-plate compounds** (the bulk of each T2+C scenario):

| Scenario | # Cpds | Median reps | 1x | 2x | 3x | 4-5x | 6-10x | >10x |
|----------|--------|------------|-----|-----|-----|------|-------|------|
| S3 | 82,152 | **2** | 14,185 | 40,688 | 24,194 | 2,686 | 377 | 22 |
| S5/S8 | 82,288 | **4** | 129 | 1,829 | 35,112 | 42,545 | 2,223 | 450 |
| A3 | 112,332 | **2** | 23,876 | 59,824 | 23,489 | 2,744 | 2,375 | 24 |
| A9 | 35,677 | **7** | 1,839 | 21 | 6 | 4,701 | 29,092 | 18 |
| C4 | 82,288 | **4** | 129 | 1,088 | 11,302 | 62,352 | 6,764 | 653 |
| C8 | 85,565 | **3** | 792 | 24,409 | 34,352 | 21,688 | 2,722 | 1,602 |
| D3 | 115,732 | **5** | 1 | 151 | 1,373 | 69,122 | 41,948 | 3,137 |

**Key observations:**

1. **S3 and A3 are dominated by 1-2x replicates** — most compounds have only 1-2 wells. This is because within Wave 1, each compound was typically screened at 1 source with ~2 wells per source. A3 has more total compounds (4 sources) but same median of 2.

2. **S5/S8 and C4 jump to median 4x** — adding source_3 (Opera Phenix) and source_8 (ImageXpress) brings compounds that each have ~4 wells (their own source's replicates). The within-source replication is the main driver, not cross-source overlap.

3. **A9 (Wave 2) has the highest median at 7x** — despite only 3 sources, Wave 2's exchange design gives each compound ~7 replicates. Wave 2 compounds are better replicated per compound than Wave 1.

4. **C8 has a long tail of well-replicated compounds** — 1,602 compounds with >10 replicates (the S3↔S7 overlap), plus good median of 3x across the bulk. Best combination of breadth and depth.

5. **D3 is the most uniformly replicated** — median 5x, only 1 compound with 1 replicate. At full scale, nearly every compound has usable replication. But correction quality at 11 sources is uncertain.

6. **For ML: median 4+ replicates is the practical threshold** for reliable per-compound features. This makes C4, A9, and D3 the most ML-ready scenarios. S3 and A3 (median 2) will have noisier per-compound features.

#### Reference mapping: compound overlap between train and test

| ID | Mapping | Train cpds | Test cpds | Shared | % test shared | Total eval (shared+T2) |
|----|---------|-----------|----------|--------|--------------|----------------------|
| F1 | CV8000 → IX | 36,336 | 12,211 | **8** | **0.1%** | 314 |
| F2 | CV8000 W1 → CV8000 W2 | 36,336 | 12,000 | 767 | 6.4% | 1,073 |
| F3 | CV8000 → Opera Phenix | 36,336 | 11,033 | 486 | 4.4% | 792 |
| **F4** | **Full W1 → Full W2** | **59,079** | **36,023** | **5,086** | **14.1%** | **5,392** |
| **F5** | **W1 → Bridge (S7)** | **59,079** | **6,061** | **2,975** | **49.1%** | **3,281** |
| F6 | W1 384 → W1 1536 | 59,455 | 11,999 | 181 | 1.5% | 487 |

**Notable:** F5 (mapping source_7) has 49.1% compound overlap with the W1 training set — S7 is genuinely a bridge source with strong compound overlap, making it the best-powered single-source mapping test.

**Key implications:**

1. **F1 is nearly impossible to evaluate with compounds** — only 8 shared compounds between CV8000 and ImageXpress. Evaluation relies almost entirely on TARGET2 (306 compounds). This limits confidence in cross-platform mapping quality.

2. **F4 is the best reference mapping scenario** — 5,086 shared compounds (14.1% of test set) plus 306 TARGET2 gives ~5,400 evaluable compounds. This is enough for meaningful mAP.

3. **C8 has the best compound replication of any single scenario** — 6.3% replicated, driven by the source_3↔source_7 overlap (2,884 compounds). This makes it the best choice for within-scenario compound-level evaluation.

4. **The replication sparsity means mAP evaluation is always anchored primarily on TARGET2.** For COMPOUND-level evaluation, only scenarios with source_7 (bridge) or cross-wave overlap provide meaningful numbers of replicated compounds.

5. **For ML-readiness, replication matters less than total compound count.** A compound in only 1 source still gets batch-corrected features. The question is whether those features are *reliable* — which is what TARGET2 evaluation answers. A model trained on C4's 71K compounds will have broader coverage than A3's 48K, even if most are unreplicated.

---

### Existing scenarios (S1-8, Apricot)

| ID | Utility | Ref-map | Embed quality | ML-ready | Feasibility | Notes |
|----|---------|---------|---------------|----------|-------------|-------|
| S1 | Med | Low | High | Low | Done | Single source, intra-batch only. Good baseline but too narrow for reference. |
| S2 | Med | Low | High | Low | Done | Same microscope, T2 only. Clean embedding but only 306 compounds. |
| S3 | Med | Med | High | Med | Running | Same microscope, T2+C. Good ML features within CV8000 labs, but no microscope diversity for generalizable reference. |
| S4 | Med | Low | Med | Low | Pending | 3 microscope types, T2 only. Tests mixed-microscope but too few compounds for ML. |
| S5 | High | Med | Med | Med | Pending | 3 microscope types, T2+C. Key benchmark scenario from Arevalo. Good compound diversity. |
| S6 | Med | Low | High | Low | Not run | Hierarchical batch, only 2 sources. |
| S7 | High | Low | Med | Low | Not run | 5 sources, hierarchical, T2 only. Key benchmark. |
| S8 | High | Med | Low-Med | Med | Partial | The hardest existing scenario. 370K cells, 3 microscope types, T2+C, MOA eval. Best existing benchmark but correction quality is uncertain. |
| Apricot | High | Med-High | Low | Med | Not run | 11 sources, 5 microscope types, T2. Most diverse but correction will be hard. Would be a strong reference if correction succeeds. |

---

### A. Same-microscope scenarios

| ID | Utility | Ref-map | Embed quality | ML-ready | Feasibility | Assessment |
|----|---------|---------|---------------|----------|-------------|------------|
| A1 | **High** | Med | High | Low | Easy | **Cross-wave with controlled microscope — novel.** Directly tests whether wave membership matters beyond microscope effects. Only scenario to isolate this. T2-only limits ML use. |
| A2 | — | — | — | — | — | Duplicate of S3. Skip. |
| A3 | **High** | **High** | High | **High** | Easy | **Best same-microscope reference atlas candidate.** 4 CV8000 sources, cross-wave, T2+C. Clean correction expected (same microscope). ~48K compounds with cross-wave overlap for eval. High ML-readiness: large, clean, well-evaluated. |
| A4 | Med | Low | Med | Low | Easy | Widefield-only but mixed plate format (384+1536) confounds interpretation. S4 has only 367 compounds. |
| A5 | Low | Low | High | Low | Easy | Only 2 sources, very small. Too narrow. |
| A6 | Med | Low | Med | Low | Easy | Widefield cross-wave, but mixed plate format and small. |
| A7 | Low | Low | Unknown | Low | Med | Same institution — interesting curiosity but not useful as reference. S13 has no compounds. |
| A8 | **High** | Med | Med | Low | Easy | **Wave 2 standalone — novel.** 3 sources, 3 different microscopes. Direct comparison with S4 (Wave 1, same setup). Critical for testing wave effects. |
| A9 | **High** | Med | Med | Med | Easy | A8 + compounds. Better ML-readiness but Wave 2 has fewer shared compounds → weaker evaluation. |

---

### B. Cross-microscope scenarios

| ID | Utility | Ref-map | Embed quality | ML-ready | Feasibility | Assessment |
|----|---------|---------|---------------|----------|-------------|------------|
| B1 | Med | Low | Med | Low | Easy | Minimal confocal vs widefield. Too small for reference use. |
| B2 | Med | Med | Med | Low | Easy | Balanced 2v2 confocal/widefield. Clean design but T2-only. |
| B3 | Med | Med | Med | Low | Easy | 5 sources all Wave 1 — similar to S4 but with S15 instead of S8. Not novel enough. |
| B4 | Med | Med | Med | Low | Easy | Cross-confocal-platform. Useful niche question but small. |
| B5 | Med | Med | Med | Low | Easy | Yokogawa generations. Niche. |
| B6 | **High** | Med | Med | Low | Easy | **All confocal, 3 platforms, 5 sources.** Clean modality-controlled design. Good for confocal-only reference. |
| B7 | **High** | Low | Med | Low | Easy | **Plate format isolation.** Same microscope (OP), same modality (widefield), only difference is 384 vs 1536 well. Novel and important. Not useful as reference but answers a key technical question. |

---

### C. Wave-focused scenarios

| ID | Utility | Ref-map | Embed quality | ML-ready | Feasibility | Assessment |
|----|---------|---------|---------------|----------|-------------|------------|
| C1 | **High** | Med-High | Med | Low | Easy | **Full Wave 1, T2.** 7 sources, 3 microscope types. The "complete Wave 1 reference" without compound complexity. |
| C2 | **High** | **High** | Med | **High** | Med | **Full Wave 1, T2+C. Best candidate for a general-purpose reference atlas.** 7 sources, 3 microscope types, maximum compound diversity within one wave. Richest compound overlap for evaluation. The flagship scenario for the paper. Very large — compute-intensive. |
| C3 | **High** | **High** | Med-High | Low | Easy | C1 minus 1536-well source. Cleaner than C1, 6 sources, still 3 microscope types. Good reference, just T2-only. |
| C4 | **High** | **High** | Med | **High** | Med | **C3 + compounds + hierarchical batch. The premium reference atlas.** Clean plate format, 6 Wave 1 sources, maximum compound diversity, hierarchical correction. If correction quality is good, this is the best embedding for ML. |
| C5 | **High** | Low | Med | Low | Easy | Same as A8. Important for wave comparison but too small for general reference. |
| C6 | **High** | Med-High | Low-Med | Med | Med | **Cross-wave, 6 sources, T2.** Key test of whether cross-wave correction works. Lower expected quality than within-wave, but if successful, more generalizable reference. |
| C7 | **High** | Med-High | Low-Med | Med | Med | C6 + bridge source S7. Tests whether bridge compounds improve cross-wave correction. Novel. |
| C8 | **High** | Med | Med | Med | Med | S7 as bridge with compounds. The S3-S7 overlap (2884 compounds) makes this a strong evaluation scenario. |

---

### D. Scale scenarios

| ID | Utility | Ref-map | Embed quality | ML-ready | Feasibility | Assessment |
|----|---------|---------|---------------|----------|-------------|------------|
| D1 | **High** | **High** | Low | Med | Hard | **Max diversity, clean plate format.** 11 sources, all 384-well, T2. If correction works, the most generalizable reference possible without plate format confound. Hard to correct. |
| D2 | Med | Med | Low | Med | Hard | All 13 sources including 1536-well. Plate format confound weakens it vs D1. Diminishing returns. |
| D3 | **High** | **High** | Low | **High** | Hard | **D1 + compounds.** The theoretical maximum ML-ready embedding. 11 sources, ~117K compounds. But correction at this scale may fail. |

---

### E. Batch key variants

| ID | Utility | Ref-map | Embed quality | ML-ready | Feasibility | Assessment |
|----|---------|---------|---------------|----------|-------------|------------|
| E1 | Med | Low | High | Low | Easy | Incremental over S2/S3. Useful for methods paper completeness. |
| E2 | Med | Low | Med | Low | Easy | Ablation study. Useful for methods paper completeness. |
| E3 | **High** | Med | Unknown | Med | Easy | **Microscope as batch key — novel.** Tests a fundamentally different batch correction axis. If microscope is the primary driver, this should outperform source-based correction. Important methodological finding. |

---

### F. Reference mapping scenarios

| ID | Utility | Ref-map | Embed quality | ML-ready | Feasibility | Assessment |
|----|---------|---------|---------------|----------|-------------|------------|
| F1 | **High** | Med | Med | Low | Med | **Cross-platform confocal mapping.** CV8000 ref → IX query. Same wave, same modality, different manufacturer. **BUT only 8 shared compounds (0.1%)** — evaluation relies entirely on TARGET2. Confidence in results will be limited. |
| F2 | **High** | **High** | High | Med | Med | **Cross-wave, same microscope.** CV8000 W1 → S5 W2. 767 shared compounds (6.4% of test) + 306 TARGET2 ≈ 1,073 evaluable compounds. Best controlled test of cross-wave generalization. |
| F3 | **High** | **High** | Med | Med | Med | **Cross-modality mapping.** Confocal ref → widefield query. 486 shared compounds (4.4%) + 306 TARGET2 ≈ 792 evaluable. The hardest but most important test. |
| F4 | **High** | **High** | Med | **High** | Hard | **Full-scale cross-wave mapping.** W1 ref, W2 query. **5,086 shared compounds (14.1%)** + 306 TARGET2 ≈ 5,400 evaluable. By far the best-powered ref-map test. Definitive. |
| F5 | **High** | Med | Med | Med | Med | Bridge source mapping. S7 has 2,884 compounds shared with S3. Novel angle on bridge role. |
| F6 | Med | Med | Med | Low | Med | Within-wave, cross-plate-format. Niche. |

---

### G. Perturbation modality scenarios

| ID | Utility | Ref-map | Embed quality | ML-ready | Feasibility | Assessment |
|----|---------|---------|---------------|----------|-------------|------------|
| G1 | **High** | Med | Unknown | Med | **Hard** | ORF perturbations. Novel for batch correction benchmarking. But requires pipeline changes (different label key, different eval logic). |
| G2 | **High** | Med | Unknown | Med | **Hard** | CRISPR perturbations. Same novelty and difficulty as G1. |
| G3 | **Very High** | **High** | Unknown | **High** | **Very Hard** | **Cross-modality integration (compound + ORF).** The holy grail — a unified embedding of chemical and genetic perturbations. Would be a landmark result. But requires significant pipeline work and novel evaluation approaches. |

---

### Top-tier recommendations by use case

#### For the batch correction methods paper (C1 priority: reproduce + extend Arevalo)
1. **S1-S5** (existing, in progress) — reproduce Arevalo
2. **A1** — cross-wave, same microscope (novel finding)
3. **A8/C5** — Wave 2 standalone (novel, mirrors S4 design)
4. **B7** — plate format isolation (novel technical finding)
5. **E3** — microscope as batch key (novel methodology)

#### For the best general-purpose reference atlas
1. **C4** — Wave 1, 384-well, 6 sources, T2+C, hierarchical batch. 71K compounds, 0.9% replicated + 306 TARGET2. Best quality/diversity tradeoff.
2. **A3** — CV8000 only, 4 sources (cross-wave), T2+C. 48K compounds, 1.8% replicated. Highest expected embedding quality. Best if downstream users also have CV8000 data.
3. **C8** — W1+bridge (S7), 5 sources, T2+C. 50K compounds, **6.3% replicated** (best single-scenario replication). Source_7's 2,884-compound overlap with S3 gives the strongest compound-level evaluation.
4. **C2** — Full Wave 1 including S1 (1536-well). 83K compounds but plate format confound.
5. **D1** — 11 sources, T2 only. Most generalizable if correction succeeds, but only 306 evaluable compounds.

#### For ML-ready embeddings (compound activity prediction, hit calling)
1. **C4** — 6 sources, 71K compounds (620 replicated in 2+ sources). Clean design, best breadth.
2. **A3** — 4 sources, 48K compounds (843 replicated). Easiest correction → most reliable per-compound features.
3. **C8** — 5 sources, 50K compounds (3,147 replicated). Best replication for cross-validated ML.
4. **D3** — 11 sources, 117K compounds (9,040 replicated). Maximum breadth but hardest correction — features may be noisy.

#### For reference mapping validation
1. **F4** → train C3/C4 (full W1), query A8 (full W2). **5,086 shared compounds** — by far the best-powered evaluation. The definitive cross-wave test.
2. **F2** → train A3 (CV8000 W1), query S5 (CV8000 W2). 767 shared compounds + 306 TARGET2. Same microscope, cross-wave — cleanest isolation of wave effect.
3. **F3** → train A3 (CV8000 W1), query S3 (OP W1). 486 shared compounds + 306 TARGET2. Cross-modality, same wave — tests microscope generalization.
4. **F1** → train CV8000, query S8 (IX). **Only 8 shared compounds** — essentially TARGET2-only evaluation (306). Results will be suggestive but not definitive. Consider dropping or reframing as "zero-shot mapping."

#### Replication-aware strategy note
The replication sparsity has a key implication: **mAP evaluation is fundamentally TARGET2-limited for most scenarios.** The 306 TARGET2 compounds, not COMPOUND plate overlap, drive evaluation confidence. This means:
- T2-only scenarios aren't as evaluation-limited as they appear (you're already using ~90%+ of the evaluable compounds)
- T2+C scenarios mainly add value for **ML breadth** (more compound features), not for evaluation power
- The exception is **C8** (source_7 bridge: 3,147 replicated) and **F4** (cross-wave: 5,086 shared) — these are the only scenarios where COMPOUND plate replication materially improves evaluation beyond TARGET2

---

## 9. Key Findings from Arevalo et al. (2024)

- Microscope type was identified as a **dominant batch factor** in the paper
- Scenario complexity increases from 1→8 primarily by adding more sources (= more microscope types)
- The paper found that method rankings were **not stable across scenarios** — methods that performed well on simple scenarios sometimes performed poorly on complex ones
- TARGET2 compounds serve as shared anchors across all sources, enabling cross-source evaluation
- The 306 TARGET2 compounds have known MOA annotations, enabling biological evaluation

---

## References

- Arevalo et al. (2024). "Evaluating batch correction methods for image-based cell profiling." Nature Communications 15, 6698.
- Chandrasekaran et al. (2023). "JUMP Cell Painting dataset." bioRxiv.
- JUMP-CP Consortium: https://jump-cellpainting.broadinstitute.org/
- JUMP datasets repo: https://github.com/jump-cellpainting/datasets
