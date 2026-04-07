# Plan: Split Seurat into v4 and v5 methods

## Goal
Explicitly version Seurat methods so defaults (Arevalo reproduction) use v4 and HPO runs use v5. Both versions available for all scenarios.

## Renames (v5, existing logic)

| Current | New |
|---------|-----|
| `scripts/correct_with_seurat.R` | `scripts/correct_with_seurat_v5.R` |
| `scripts/run_seurat_trial.R` | `scripts/run_seurat_trial_v5.R` |
| `scripts/optimise_seurat.py` | `scripts/optimise_seurat_v5.py` |
| `containers/r.def` / `r.sif` | `containers/r_v5.def` / `r_v5.sif` |
| `inputs/defaults/optuna_seurat_cca.csv` | `inputs/defaults/optuna_seurat_cca_v5.csv` |
| `inputs/defaults/optuna_seurat_rpca.csv` | `inputs/defaults/optuna_seurat_rpca_v5.csv` |
| Method `seurat_cca` | `seurat_cca_v5` |
| Method `seurat_rpca` | `seurat_rpca_v5` |

## New files (v4)

| File | Based on |
|------|----------|
| `scripts/correct_with_seurat_v4.R` | Original v4 script (commit 1ed542f) + optparse CLI |
| `scripts/run_seurat_trial_v4.R` | v4 API adaptation of run_seurat_trial |
| `scripts/optimise_seurat_v4.py` | Copy of optimise_seurat_v5.py, references v4 R script |
| `containers/r_v4.def` | Pin r-seurat=4.4 from conda-forge |
| `inputs/defaults/optuna_seurat_cca_v4.csv` | Copy current defaults (Arevalo params, tuned for v4) |
| `inputs/defaults/optuna_seurat_rpca_v4.csv` | Copy current defaults |

## Rule changes

### Snakefile METHODS list
```python
METHODS = [
    ...
    "seurat_cca_v4",
    "seurat_rpca_v4",
    "seurat_cca_v5",
    "seurat_rpca_v5",
    ...
]
```

### tune.smk
- Rename `optimize_seurat_cca` → `optimize_seurat_cca_v5`, `optimize_seurat_rpca` → `optimize_seurat_rpca_v5`
- Add `optimize_seurat_cca_v4`, `optimize_seurat_rpca_v4` (same pattern, different script + container)
- Update `use_default_params` ruleorder

### correct.smk
- Rename `methods_seurat_cca` → `methods_seurat_cca_v5`, `methods_seurat_rpca` → `methods_seurat_rpca_v5`
- Add `methods_seurat_cca_v4`, `methods_seurat_rpca_v4` (v4 script + r_v4.sif container)

## Files touched
1. `Snakefile` — METHODS list, auto-skip logic
2. `rules/tune.smk` — rename + add rules, ruleorder
3. `rules/correct.smk` — rename + add rules
4. `scripts/correct_with_seurat.R` → rename to `_v5.R`
5. `scripts/run_seurat_trial.R` → rename to `_v5.R`
6. `scripts/optimise_seurat.py` → rename to `_v5.py`
7. `scripts/correct_with_seurat_v4.R` — new
8. `scripts/run_seurat_trial_v4.R` — new
9. `scripts/optimise_seurat_v4.py` — new
10. `containers/r.def` → rename to `r_v5.def`
11. `containers/r_v4.def` — new
12. `inputs/defaults/optuna_seurat_{cca,rpca}.csv` → rename to `_v5.csv`
13. `inputs/defaults/optuna_seurat_{cca,rpca}_v4.csv` — new (copy of current defaults)
14. `pixi.toml` — build tasks for r_v4.sif

## Risks
- Existing S1-S4 outputs reference `seurat_cca`/`seurat_rpca` (old names) in optuna CSVs and correction parquets. These won't be found by the new method names. Need to either:
  - Rename existing output files, OR
  - Re-run with new names (preferred — HPO is fast for Seurat on small scenarios)
- fastMNN also uses r.sif → needs updating to r_v5.sif (or keep r.sif as symlink to r_v5.sif)

## Build order
1. Create r_v4.def and build r_v4.sif
2. Rename scripts and containers
3. Write v4 scripts
4. Update rules
5. Dry-run to verify DAG
6. Build r_v4.sif on a node with /localscratch
