# Scenario Registry

Structured documentation of all batch correction scenarios for the JUMP-CP integration benchmark. Each scenario card describes what it tests, why it matters, and what makes it hard.

**Background research**: See `tasks/jump_metadata_research.md` for full metadata details, compound wave structure, and microscope specifications.

---

## Card format

Each scenario is documented with:

- **Config**: sources, plate types, batch key, eval key
- **Microscopes**: which systems and modalities are present
- **Compounds**: how many, replication structure, evaluation power
- **Challenges**: what makes this scenario hard for batch correction
- **Use cases**: what the corrected embedding is good for
- **Status**: not started / running / done

---

## Existing scenarios (Arevalo et al. reproduction)

### S1 — Intra-source batch correction

| Field | Value |
|-------|-------|
| Sources | source_6 (1 source) |
| Plate types | TARGET2 |
| Batch key | Metadata_Batch |
| Eval key | Metadata_JCP2022 |
| Wave | Wave 1 |
| Microscopes | CV8000 (confocal) |
| Plate format | 384-well |
| Compounds | 302 TARGET2 (median 21 reps/cpd within source) |
| Evaluation | 302 TARGET2 compounds |

**Challenges**: Minimal — single microscope, single source, only batch-to-batch variation within one lab. The easiest scenario.

**Use cases**: Baseline. How well can methods remove within-lab plate/batch effects? Not useful as a reference atlas (too narrow).

**Status**: Done (2026-03-28). PR #29.

---

### S2 — Same-microscope cross-source (TARGET2)

| Field | Value |
|-------|-------|
| Sources | source_2, source_6, source_10 (3 sources) |
| Plate types | TARGET2 |
| Batch key | Metadata_Source |
| Eval key | Metadata_JCP2022 |
| Wave | Wave 1 only |
| Microscopes | CV8000 (confocal) — all identical |
| Plate format | 384-well |
| Compounds | 302 TARGET2 (median 37 reps/cpd) |
| Evaluation | 302 TARGET2 compounds |

**Challenges**: Low — all sources use the same microscope (CV8000), so batch effects are limited to lab-specific protocol differences (reagent lots, cell passage, image acquisition settings).

**Use cases**: Tests whether lab-specific effects are removable when microscope type is controlled. Not useful for ML (only 302 compounds).

**Status**: Done (2026-03-28). PR #30.

---

### S3 — Same-microscope cross-source (TARGET2 + COMPOUND)

| Field | Value |
|-------|-------|
| Sources | source_2, source_6, source_10 (3 sources) |
| Plate types | TARGET2 + COMPOUND |
| Batch key | Metadata_Source |
| Eval key | Metadata_JCP2022 |
| Wave | Wave 1 only |
| Microscopes | CV8000 (confocal) — all identical |
| Plate format | 384-well |
| Compounds | 82,152 COMPOUND (median 2 reps, 14K at 1x, 41K at 2x) + 302 TARGET2 |
| Cross-source overlap | 76 compounds in 2+ sources (0.2%) |
| Evaluation | 302 TARGET2 + 76 replicated COMPOUND = 378 evaluable |

**Challenges**: Low correction difficulty (same microscope). Larger data volume than S2. Most compounds have only 1-2 replicates.

**Use cases**: First ML-usable embedding (82K compounds), but per-compound features are noisy (median 2 reps). Good for within-CV8000 applications.

**Status**: Running (2026-03-29).

---

### S4 — Mixed-microscope cross-source (TARGET2)

| Field | Value |
|-------|-------|
| Sources | source_2, source_3, source_6, source_8, source_10 (5 sources) |
| Plate types | TARGET2 |
| Batch key | Metadata_Source |
| Eval key | Metadata_JCP2022 |
| Wave | Wave 1 only |
| Microscopes | CV8000 ×3, Opera Phenix ×1, ImageXpress ×1 (3 types, 2 modalities) |
| Plate format | 384-well |
| Compounds | 302 TARGET2 (median 66 reps/cpd) |
| Evaluation | 302 TARGET2 compounds |

**Challenges**: Medium — introduces microscope heterogeneity. Confocal-dominated (3 CV8000 vs 1 widefield + 1 other confocal). Batch effects now include microscope type, not just lab protocols.

**Use cases**: Tests whether methods handle microscope diversity. Arevalo found rankings shift at this transition. Not useful for ML (T2 only).

**Status**: Pending.

---

### S5 — Mixed-microscope cross-source (TARGET2 + COMPOUND)

| Field | Value |
|-------|-------|
| Sources | source_2, source_3, source_6, source_8, source_10 (5 sources) |
| Plate types | TARGET2 + COMPOUND |
| Batch key | Metadata_Source |
| Eval key | Metadata_JCP2022 |
| Wave | Wave 1 only |
| Microscopes | CV8000 ×3, Opera Phenix ×1, ImageXpress ×1 (3 types, 2 modalities) |
| Plate format | 384-well |
| Compounds | 82,288 COMPOUND (median 4 reps, 129 at 1x, 35K at 3x, 42K at 4-5x) + 302 TARGET2 |
| Cross-source overlap | 549 compounds in 2+ sources (0.9%) |
| Evaluation | 302 TARGET2 + 549 replicated COMPOUND = 851 evaluable |

**Challenges**: Medium-hard — microscope diversity + large data. This is the key Arevalo benchmark where rankings begin to diverge from simpler scenarios.

**Use cases**: Primary Arevalo benchmark with compounds. Reasonable ML embedding (82K compounds, median 4 reps). The within-source replication is decent.

**Status**: Pending.

---

### S6-S8 — Hierarchical batch key variants

S6 (2 sources, T2, hierarchical), S7 (5 sources, T2, hierarchical), S8 (5 sources, T2+C, hierarchical + MOA eval). See existing config files. S8 is the hardest existing scenario (370K cells).

**Status**: S8 partial (scenario 8 blockers documented in todo.md). S6-S7 not run.

---

## New scenarios

### wave2 — Wave 2 standalone (TARGET2 + COMPOUND)

| Field | Value |
|-------|-------|
| Sources | source_5, source_9, source_11 (3 sources) |
| Plate types | TARGET2 + COMPOUND |
| Batch key | Metadata_Source |
| Eval key | Metadata_JCP2022 |
| Wave | **Wave 2 only** |
| Microscopes | CV8000 (confocal), Opera Phenix (widefield), Operetta (widefield) — **3 different microscopes, 2 modalities** |
| Plate format | **Mixed: source_5 384-well, source_9 1536-well, source_11 384-well** |
| Compounds | 35,677 COMPOUND + 302 TARGET2 |
| Per-source reps | source_5: median 2, source_9: median 4, source_11: median 2 |
| Combined reps | **Median 7 reps/cpd** (typical: 2 from S5 + 4 from S9 + 2 from S11) |
| Cross-source overlap | **30,244 compounds (84.8%) in ALL 3 sources** |
| Evaluation | 302 TARGET2 + 30,244 cross-source replicated = **30,546 evaluable** |

**Challenges** (all simultaneously):
1. **Cross-microscope**: 3 entirely different microscope systems (CV8000, Opera Phenix, Operetta) from 2 different manufacturers
2. **Cross-modality**: 1 confocal + 2 widefield — tests whether correction can bridge imaging modalities
3. **Cross-plate-format**: source_9 uses 1536-well plates while S5 and S11 use 384-well. Different well sizes → different cell densities, imaging field sizes, edge effects
4. **Small source count**: Only 3 sources — methods that need many batches to learn batch structure may struggle

**Use cases**:
1. **Best evaluation scenario in the entire dataset**: 30K+ compounds replicated across all 3 sources gives extremely well-powered mAP evaluation. No other scenario comes close.
2. **Cross-validated ML features**: Median 7 replicates per compound across 3 independent labs. Per-compound features are the most reliable of any scenario.
3. **Microscope generalization test**: If correction works despite max microscope heterogeneity, the embedding is microscope-agnostic.
4. **Direct comparison with S4/S5**: Same number of microscope types (3), similar data size, but completely different sources and wave. Do rankings match? If not, this proves method performance is dataset-dependent.
5. **Wave 2 reference atlas**: A corrected A9 embedding could serve as a reference for mapping new Wave 2 data. The 84.8% compound sharing means new data from any Wave 2 source has extensive overlap.

**Why this is unique**: A9 has the highest heterogeneity-per-source of any scenario (every source is different) combined with the highest compound replication (84.8% in all sources). This is the toughest correction challenge with the most trustworthy evaluation. If a method scores well here, you know it works.

**Status**: Not started.

---

### wave1 — Wave 1 full (TARGET2 + COMPOUND)

| Field | Value |
|-------|-------|
| Sources | source_1, source_2, source_3, source_6, source_8, source_10, source_15 (7 sources) |
| Plate types | TARGET2 + COMPOUND |
| Batch key | Metadata_Source |
| Eval key | Metadata_JCP2022, Metadata_DRH_MOA |
| Wave | **Wave 1 only** |
| Microscopes | CV8000 ×3 (S2,6,10), Opera Phenix ×2 (S1,3), ImageXpress ×1 (S8) — **3 types, 2 modalities** |
| Plate format | **Mixed: source_1 = 1536-well, rest = 384-well** |
| Compounds | 82,288 COMPOUND + 302 TARGET2 |
| Per-source reps | Median 1 rep per source (Wave 1 design: breadth over depth) |
| Combined reps | **Median 5 reps/cpd** (most compounds exchanged to 4-5 sources, ~1 well each) |
| Cross-source overlap | **61,285 compounds (74.5%) in 5+ sources**, 359 (0.4%) in all 7 |
| Evaluation | 302 TARGET2 + 81,929 cross-source replicated (99.6%) = **82,231 evaluable** |

**Challenges**:
1. **Microscope heterogeneity**: 3 microscope types, 2 imaging modalities — confocal-dominated (4 confocal vs 2 widefield + 1 mixed)
2. **Plate format**: source_1 uses 1536-well plates. Can be dropped for a clean version (→ C4 with 6 sources)
3. **Scale**: 7 sources, likely 400-600K cells
4. **Sparse per-source replication**: Each source has median 1 well per compound (unlike Wave 2's 2-4). The replication comes from cross-source exchange, not within-source repeats.

**Use cases**:
1. **Natural counterpart to wave2**: Direct comparison — same plate types (T2+C), similar compound counts (82K vs 36K), but fundamentally different replication structure (many sources × 1 rep vs few sources × many reps)
2. **Broadest reference atlas**: 82K compounds in a well-replicated design. 74.5% of compounds in 5+ sources means excellent cross-lab validation.
3. **ML features**: Median 5 reps from 5 different labs is excellent for per-compound reliability — even though each lab contributes ~1 rep, having 5 independent measurements from different microscopes/protocols is arguably more informative than 7 reps from 3 labs (wave2).
4. **Community reference**: Wave 1 has the most consortium partners (7 of 10 pharma companies). An embedding here is maximally representative of consortium practices.

**Comparison with wave2**:

| | wave1 | wave2 |
|---|-------|-------|
| Sources | 7 | 3 |
| Microscope types | 3 | 3 |
| Unique compounds | 82,288 | 35,677 |
| Median reps/cpd | 5 (from 5 sources × ~1 well) | 7 (from 3 sources × ~2-4 wells) |
| Cpds in all sources | 359 (0.4%) | 30,244 (84.8%) |
| Cpds in 2+ sources | 81,929 (99.6%) | 30,486 (85.5%) |
| Replication nature | Cross-lab (independent protocols) | Cross-lab (independent protocols) |
| Evaluation power | Very high (82K evaluable) | Very high (30K evaluable) |
| Correction difficulty | Medium-hard (7 sources) | Hard (max heterogeneity per source) |
| ML strength | Breadth (82K cpds) | Depth (7 reps/cpd) |

**Key insight**: wave1 and wave2 are complementary. Wave 1 gives breadth (82K compounds, 5-source cross-validation). Wave 2 gives depth (30K compounds, 7-replicate reliability per compound, 84.8% universal sharing). Together they cover different aspects of ML-readiness.

**Status**: Not started.

---

### A1 — CV8000 cross-wave (TARGET2)

| Field | Value |
|-------|-------|
| Sources | source_2, source_5, source_6, source_10 (4 sources) |
| Plate types | TARGET2 |
| Batch key | Metadata_Source |
| Eval key | Metadata_JCP2022 |
| Wave | **Wave 1 + Wave 2** (S2,6,10 = W1; S5 = W2) |
| Microscopes | CV8000 (confocal) — **all identical** |
| Plate format | 384-well |
| Compounds | 302 TARGET2 (median 61 reps/cpd) |
| Evaluation | 302 TARGET2 |

**Challenges**: Low correction difficulty (same microscope, same plate format). The only new factor vs S2 is cross-wave membership.

**Use cases**:
1. **Isolate wave effect**: Does production wave matter when microscope is controlled? Comparison with S2 (same microscope, Wave 1 only) directly answers this.
2. **Cleanest cross-wave test**: No microscope or plate format confounders.

**Status**: Not started.

---

### A3 — CV8000 cross-wave (TARGET2 + COMPOUND)

| Field | Value |
|-------|-------|
| Sources | source_2, source_5, source_6, source_10 (4 sources) |
| Plate types | TARGET2 + COMPOUND |
| Batch key | Metadata_Source |
| Eval key | Metadata_JCP2022 |
| Wave | Wave 1 + Wave 2 |
| Microscopes | CV8000 (confocal) — all identical |
| Plate format | 384-well |
| Compounds | 112,332 COMPOUND (median 2 reps, 24K at 1x, 60K at 2x) + 302 TARGET2 |
| Cross-source overlap | 843 compounds in 2+ sources (1.8%) |
| Evaluation | 302 TARGET2 + 843 replicated COMPOUND = 1,145 evaluable |

**Challenges**: Low correction difficulty. Largest compound count of any same-microscope scenario (112K). Most compounds have only 1-2 reps.

**Use cases**:
1. **Best same-microscope reference atlas**: Highest expected embedding quality (easiest correction) with broadest compound coverage. Best for downstream users who also have CV8000 data.
2. **ML features for CV8000 labs**: 112K compounds, but median 2 reps limits per-compound reliability.
3. **Cross-wave with compounds**: The 843 cross-wave overlapping compounds enable limited cross-wave evaluation beyond TARGET2.

**Status**: Not started.

---

### C4 — Wave 1 clean (384-well, TARGET2 + COMPOUND, hierarchical)

| Field | Value |
|-------|-------|
| Sources | source_2, source_3, source_6, source_8, source_10, source_15 (6 sources) |
| Plate types | TARGET2 + COMPOUND |
| Batch key | Metadata_Source + Metadata_Batch |
| Eval key | Metadata_JCP2022, Metadata_DRH_MOA |
| Wave | Wave 1 only |
| Microscopes | CV8000 ×3, Opera Phenix ×2, ImageXpress ×1 (3 types, 2 modalities) |
| Plate format | 384-well (all) |
| Compounds | 82,288 COMPOUND (median 4 reps, 62K at 4-5x) + 302 TARGET2 |
| Cross-source overlap | 648 compounds in 2+ sources (0.9%) |
| Evaluation | 302 TARGET2 + 648 replicated COMPOUND = 950 evaluable |

**Challenges**:
1. **Microscope heterogeneity**: 3 microscope types, 2 modalities — same difficulty as S5/S8
2. **Hierarchical batch**: Source + Batch doubles the batch key complexity
3. **Scale**: 6 sources, likely 300-500K cells
4. **MOA evaluation**: Tests whether batch correction preserves mechanism-of-action signal

**Use cases**:
1. **General-purpose reference atlas**: Best balance of microscope diversity (3 types), compound breadth (82K), and clean plate format (all 384-well). If correction succeeds, this is the most useful embedding for the community.
2. **ML features**: Median 4 reps, 82K compounds. Good reliability-breadth tradeoff.
3. **Wave 1 reference for F4 mapping test**: Train reference here, test by mapping Wave 2 (A9) into it.

**Status**: Not started.

---

### C8 — Wave 1 + bridge (TARGET2 + COMPOUND)

| Field | Value |
|-------|-------|
| Sources | source_2, source_3, source_6, source_7, source_10 (5 sources) |
| Plate types | TARGET2 + COMPOUND |
| Batch key | Metadata_Source |
| Eval key | Metadata_JCP2022 |
| Wave | Wave 1 + bridge (source_7) |
| Microscopes | CV8000 ×3, Opera Phenix ×1, CV7000 ×1 (3 types, confocal-dominated) |
| Plate format | 384-well |
| Compounds | 85,565 COMPOUND (median 3 reps) + 302 TARGET2 |
| Cross-source overlap | **3,149 compounds in 2+ sources (6.3%)** — best of any single scenario |
| Key overlap | source_3 ↔ source_7: 2,884 shared compounds |
| Evaluation | 302 TARGET2 + 3,149 replicated COMPOUND = **3,451 evaluable** |

**Challenges**:
1. **Bridge source integration**: source_7 is not part of Wave 1's exchange design — does it integrate cleanly?
2. **CV7000 vs CV8000**: Same manufacturer (Yokogawa) but different generation. Tests within-manufacturer variation.
3. **Unbalanced overlap**: The S3↔S7 overlap dominates. Other pairs have much less sharing.

**Use cases**:
1. **Best within-scenario evaluation**: 3,451 evaluable compounds — 10x more than S5/S8. mAP results are far more trustworthy.
2. **Bridge source validation**: Tests source_7's role as cross-wave connector. If S7 integrates well, it's evidence that bridge compounds help harmonization.
3. **Yokogawa platform comparison**: CV7000 vs CV8000 in the same embedding.

**Status**: Not started.

---

### F4 — Reference mapping: Wave 1 → Wave 2

| Field | Value |
|-------|-------|
| Train | C4 or C3 (Wave 1, 384-well, 5-6 sources) |
| Test | A9 (Wave 2, 3 sources) |
| Shared compounds | **5,086 COMPOUND + 306 TARGET2 = 5,392 evaluable** |
| % test shared | 14.1% of Wave 2 compounds appear in Wave 1 train set |

**Challenges**:
1. **Cross-wave generalization**: Train and test come from completely different production campaigns
2. **Microscope mismatch**: Train has CV8000+OP+IX; test has CV8000+OP+Operetta. Partial overlap.
3. **Plate format mismatch**: Train is all 384-well; test includes 1536-well (source_9)

**Use cases**:
1. **Definitive reference mapping test**: 5,392 evaluable compounds make this the best-powered mapping evaluation. If a reference atlas trained on Wave 1 can accurately place Wave 2 compounds, the reference is truly generalizable.
2. **Practical utility proof**: This is the real-world question — "I have a new experiment from a Wave 2 lab. Can I map it into the existing atlas and get meaningful results?"

**Requires**: C4 (or C3) and A9 to be completed first.

**Status**: Not started.

---

### B7 — Plate format isolation

| Field | Value |
|-------|-------|
| Sources | source_1, source_3, source_9 (3 sources) |
| Plate types | TARGET2 |
| Batch key | Metadata_Source |
| Eval key | Metadata_JCP2022 |
| Wave | Wave 1 (S1,S3) + Wave 2 (S9) |
| Microscopes | Opera Phenix — **all identical** |
| Plate format | **source_3 = 384-well, source_1 + source_9 = 1536-well** |
| Compounds | 302 TARGET2 |
| Cross-source overlap (COMPOUND) | 1,174 in 2+ sources (3.5%) |
| Evaluation | 302 TARGET2 |

**Challenges**: Same microscope, different plate formats. Isolates plate format as the batch variable.

**Use cases**:
1. **Technical question**: Does 384-well vs 1536-well plate format create a correctable batch effect, or a fundamental data incompatibility?
2. **Prerequisite for scale scenarios**: If plate format correction fails, scenarios including S1/S9 (D2, C1, C2) are compromised.
3. **Protocol standardization guidance**: Informs whether labs switching plate formats can still use existing reference atlases.

**Status**: Not started.

---

### E3 — Microscope type as batch key

| Field | Value |
|-------|-------|
| Sources | source_2, source_3, source_6, source_8, source_10 (5 sources, same as S4/S5) |
| Plate types | TARGET2 |
| Batch key | **Metadata_Microscope** (instead of Metadata_Source) |
| Eval key | Metadata_JCP2022 |
| Microscopes | CV8000, Opera Phenix, ImageXpress (3 batch levels instead of 5) |

**Challenges**: Treats microscope type as the batch, not source identity. All CV8000 sources (2, 6, 10) are treated as one "batch." This tests whether microscope is the primary batch axis.

**Use cases**:
1. **Methodological finding**: If correction with 3 microscope batches outperforms 5 source batches, microscope type is the dominant batch factor and the field should rethink batch definition.
2. **Practical**: Fewer batch levels = simpler correction. If this works, labs with the same microscope can be pooled.

**Status**: Not started.

---

## Scenario comparison matrix

| Scenario | #Src | Microscopes | Wave | Plates | Cpds | Med. reps | Evaluable | Difficulty | Best for |
|----------|------|------------|------|--------|------|-----------|-----------|------------|----------|
| S1 | 1 | CV8000 | W1 | T2 | 302 | 21 | 302 | Easy | Baseline |
| S2 | 3 | CV8000 | W1 | T2 | 302 | 37 | 302 | Easy | Same-microscope baseline |
| S3 | 3 | CV8000 | W1 | T2+C | 82K | 2 | 378 | Easy | CV8000 ML features |
| S5 | 5 | 3 types | W1 | T2+C | 82K | 4 | 851 | Medium | Arevalo benchmark |
| S8 | 5 | 3 types | W1 | T2+C | 82K | 4 | 851 | Hard | Hardest existing |
| **wave1** | **7** | 3 types | **W1** | T2+C | **82K** | **5** | **82,231** | Med-Hard | Broadest reference + ML breadth |
| **wave2** | **3** | 3 types | **W2** | T2+C | **36K** | **7** | **30,546** | Hard | Best evaluation + ML depth |
| **A1** | 4 | CV8000 | W1+W2 | T2 | 302 | 61 | 302 | Easy | Isolate wave effect |
| **A3** | 4 | CV8000 | W1+W2 | T2+C | 112K | 2 | 1,145 | Easy | Best same-micro reference |
| **C4** | 6 | 3 types | W1 | T2+C | 82K | 4 | 950 | Med-Hard | General reference atlas (384-only) |
| **C8** | 5 | 3 types | W1+bridge | T2+C | 86K | 3 | 3,451 | Medium | Best within-scenario eval |
| **B7** | 3 | Opera Phenix | W1+W2 | T2 | 302 | - | 302 | Unknown | Plate format test |
| **E3** | 5 | 3 types | W1 | T2 | 302 | 66 | 302 | Unknown | Microscope as batch key |
| **F4** | 6→3 | mixed | W1→W2 | T2+C | - | - | 5,392 | Hard | Reference mapping proof |
