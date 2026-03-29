# Plan: MOA Annotation Mapping for C4 Evaluation

## Background

JUMP-CP compounds are identified by `Metadata_JCP2022` IDs. For the C4 evaluation criterion, we need to assess whether batch correction preserves functional (MOA-level) signal — i.e., do compounds sharing a mechanism of action cluster together after correction?

## What already exists in this repo

### Annotation files (already downloaded)
- `inputs/metadata/repurposinghub_inchi_to_moa_moaexploded.csv` — 7,340 rows, 6,629 unique InChIKeys, 1,132 unique MOAs from Broad Drug Repurposing Hub (DRH, March 2020 version)
- `inputs/metadata/repurposinghub_inchi_to_moa_targetexploded.csv` — 14,427 rows, same compounds exploded by drug target
- `inputs/metadata/opentargets_moa_target2_eval.parquet` — 68 rows, 36 unique InChIKeys, 29 unique MOAs (very sparse)
- `inputs/metadata/compound.csv.gz` — JCP2022 → InChIKey mapping (115,796 entries)
- `inputs/target_annotations.parquet` — ChEMBL bioactivity annotations (590,042 rows, 72,654 unique InChIKeys)

### Code infrastructure (already implemented)
- `metrics/scib.py`: `_merge_with_duplication()`, `_load_repurposinghub_moa_info()`, `_load_opentargets_moa_info()`, `_load_repurposinghub_target_info()`
- `scripts/run_scibmetrics_benchmarker.py`: dispatches on MOA eval_keys via `eval_key_function_mapping`
- Scenario 7 config already includes `"Metadata_DRH_MOA"` as an eval_key (pattern to follow)

### What Arevalo et al. (2024) did
- Did **not** use MOA annotations for primary evaluation — eval_key was always `Metadata_JCP2022` (compound identity)
- Used **copairs mAP** for replicate reproducibility (same compound clusters together)
- Had a **consistency pipeline** (`rules/consistency.smk`) using CLUE drug *target* annotations for median-consensus mAP
- MOA code paths exist but were not activated for scenarios 1-5

## Coverage analysis

### TARGET2-only scenarios (1, 2, 4)
- 302 non-DMSO compounds, **100% have DRH MOA annotations**
- After filtering to MOAs with >= 3 compounds: ~31 MOA classes covering ~110 compounds (36%)
- With >= 5 compounds per MOA: only ~7 classes
- Largest: CDK inhibitor (8), PDGFR inhibitor (7), calcium channel blocker (7)

### COMPOUND+TARGET2 scenarios (3, 5)
- ~82,413 unique compounds, ~3,103 have DRH MOA (3.8%)
- After filtering >= 5: ~185 MOA classes; >= 10: ~72 classes
- Much better class balance for evaluation
- Top: cyclooxygenase inhibitor (59), serotonin receptor antagonist (45), PI3K inhibitor (45)

## Implementation plan

### Step 1: Enable MOA eval_keys in scenario configs (low effort)
Add `"Metadata_DRH_MOA"` to `eval_key` arrays in `inputs/conf/scenario_{1,2,3,4,5}.json`.
This automatically triggers MOA-based scib-metrics evaluation via existing code.

**Files**: `inputs/conf/scenario_{1,2,3,4,5}.json`

### Step 2: Enable MOA-based mAP computation (medium effort)
`metrics/map.py` lines 143-149 explicitly skip MOA eval_keys. This needs modification:
1. Detect when eval_key is an MOA key
2. Load the appropriate annotation table (reuse functions from `metrics/scib.py`)
3. Merge into adata metadata (handle multi-label duplication)
4. Compute mAP

Decision needed: use `copairs.map.multilabel.average_precision` (as in `metrics/consistency.py`) or handle duplication at data level?

**Files**: `metrics/map.py`, possibly `metrics/__init__.py`

### Step 3: Coverage statistics for paper (low effort)
Script to output per-scenario compound counts, annotated counts, MOA class counts at various thresholds, class size distributions.

**Files**: new `scripts/summarize_moa_coverage.py`

### Step 4: Cross-annotation validation (optional, defer to Phase 4)
Compare CLUE DRH vs OpenTargets vs ChEMBL annotations. Not blocking for Phase 1.

## Risks

1. **Class imbalance on TARGET2-only**: Only 7 MOAs have >= 5 compounds. Metrics like silhouette/NMI may be unreliable. Mitigation: use mAP (handles imbalance better), report per-class.
2. **Multi-MOA compounds**: 7.4% of DRH compounds have 2+ MOAs. `_merge_with_duplication` inflates sample counts. Alternative: primary MOA only, or copairs multilabel support.
3. **Annotation version**: Using DRH March 2020 version. Keep for reproducibility, note in Methods.
4. **MOA granularity**: DRH is medium granularity ("CDK inhibitor"), OpenTargets is more specific. Stick with DRH as primary (field standard).
5. **mAP skip logic in `metrics/map.py`**: Was deliberately added — MOA columns aren't in `adata.obs` at mAP time, they require external merge. This is the core technical challenge.

## Key distinction for paper
- **Compound mAP** (eval_key=Metadata_JCP2022): measures replicate reproducibility
- **MOA mAP** (eval_key=Metadata_DRH_MOA): measures functional signal preservation
- These ask fundamentally different questions and both are valuable

## Concrete next steps
1. Add `"Metadata_DRH_MOA"` to eval_keys in scenario configs 1-5
2. Modify `metrics/map.py` to support MOA-based mAP
3. Run a single scenario end-to-end to verify MOA metrics appear
4. Generate coverage statistics for paper Methods section
