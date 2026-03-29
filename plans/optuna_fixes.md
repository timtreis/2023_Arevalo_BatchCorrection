# Plan: Optuna Improvements

## Changes

### 1. Add reproducible seed to all Python scripts
- Add `sampler = optuna.samplers.TPESampler(seed=42)` to all `create_study()` calls
- 7 files: optimise_{scvi,scanvi,sysvi,scpoli,harmony,scanorama,desc}.py

### 2. Add `set.seed(42)` to R scripts
- 2 files: optimise_{fastmnn,seurat}.R
- Add `set.seed(42)` before the random search grid generation

### 3. Remove unnecessary `adata.copy()` in Python lambda calls
- The copy was protecting against mutation inside the objective, but:
  - scVI/scANVI/sysVI: mutate adata (X -= min, setup_anndata). Need copy.
  - scPoli: passes adata to constructor, but constructor may mutate. Keep copy for safety.
  - Harmony: does not mutate adata.X (converts to dense array separately). Remove copy.
  - Scanorama: does `adata[...].copy()` inside objective already, but also sorts. The outer adata isn't mutated though — the sorting creates a new view. Actually, line 29 does `adata = adata[...].copy()` which rebinds the local variable. Safe to remove outer copy.
  - DESC: `scale_bygroup(adata, ...)` mutates adata in-place, then `train(adata, ...)` mutates too. Need copy.
- Summary: Remove copy for harmony and scanorama. Keep for scvi, scanvi, sysvi, scpoli, desc.

### 4. Standardize CSV output format
- Python scripts output many extra columns from `study.trials_dataframe()` (datetime_start, datetime_complete, duration, system_attrs_*)
- R scripts output clean format: number, batch, bio, params_*, state, total
- Standardize Python to output same columns as R: number, batch, bio, params_*, state, total
- Filter `df` columns before writing

### 5. R HPO package recommendation
- For now, random search in R is acceptable for the 3 R methods (fastMNN, seurat_cca, seurat_rpca)
- The best pure-R option would be mlr3mbo (ParEGO/SMS-EGO for multi-objective)
- Alternative: Optuna via reticulate (needs Python in r.sif container)
- Document this as a future improvement, not blocking for Phase 1

## Files to modify
- scripts/optimise_scvi.py
- scripts/optimise_scanvi.py
- scripts/optimise_sysvi.py
- scripts/optimise_scpoli.py
- scripts/optimise_harmony.py
- scripts/optimise_scanorama.py
- scripts/optimise_desc.py
- scripts/optimise_fastmnn.R
- scripts/optimise_seurat.R

## Risks
- Changing seed changes HPO results — any existing optimization CSVs will differ on rerun
- Removing adata.copy() for harmony/scanorama: must verify no in-place mutation
- CSV format change: must check if downstream code depends on the extra columns
