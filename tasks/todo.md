# JUMP-CP Integration: Phase 1 Progress

**Last updated**: 2026-03-29
**Target**: Nature Methods, Q3 2026
**Research plan**: `~/lustre/projects/coach/projects/jump-cp-integration/task-list.md`

---

## Milestone 0: Pipeline Validation (smoketest) — DONE

Ran scenario_7 with `smoketest=true`, 2 HPO trials, 5 sources, dual batch keys.

- [x] All 8 containers built (base, scvi, scibmetrics, scpoli, harmony_v1, lightweight, desc, r)
- [x] 13/15 methods: HPO + correction completed
- [x] 2/15 (seurat_cca, seurat_rpca): expected failure on smoketest data size — will work on full data
- [x] Aggregated `all_methods.h5ad` generated
- [x] Apptainer migration merged to main (PR #17)

**Known issues resolved during smoketest** (details in git history):
- sysvi JAX_PLATFORMS=cuda → cpu (later fully removed, see pipeline improvements)
- harmonypy v1/v2 result shape auto-detection
- scib-metrics API version compatibility (bras vs silhouette_batch)
- Lightweight HPO metrics (3 instead of 10)
- annoy SIGILL fix (source compilation with -march=x86-64-v3)

---

## Milestone 1: C1 — Default Benchmark (scenarios 1-5) — IN PROGRESS

Run all 15 methods × 5 original scenarios × 30 HPO trials.

### Scenario status

| Scenario | Sources | Plate types | Batch key | Status |
|----------|---------|-------------|-----------|--------|
| scenario_1 | source_6 | TARGET2 | Metadata_Batch | DONE (2026-03-28). All 15 methods, metrics, plots. PR #29 merged. |
| scenario_2 | 3 sources | TARGET2 | Metadata_Source | DONE (2026-03-28). scpoli re-HPO'd with updated epochs. PR #30 merged. |
| scenario_3 | 3 sources | TARGET2+COMPOUND | Metadata_Source | RELAUNCHED 2026-04-01 on supergpu25 (H100 80GB), PID 4108981. Snakemake picked up sysvi HPO (2 remaining trials) + scpoli HPO+correction + aggregation + metrics + plots. 7/9 methods already done. R methods/DESC/scANVI auto-skipped. |
| scenario_4 | 5 sources | TARGET2 | Metadata_Source | Not started |
| scenario_5 | 5 sources | TARGET2+COMPOUND | Metadata_Source | Not started |

### Tasks

- [x] Launch scenario_1: first run 2026-03-22 (gpusrv30), rerun 2026-03-23 (gpusrv35) with all fixes
- [x] ~~Re-run scenario_1 with all pipeline fixes~~ — completed 2026-03-28 on gpusrv43. All 15 methods done. scpoli fix applied (finite pretraining epochs). DESC HPO+correction ran. Metrics and plots generated.
- [x] ~~Launch scenario_2: `pixi run scenario-2`~~ — completed 2026-03-27 on gpusrv40, re-run scpoli HPO+correction+metrics 2026-03-28 on gpusrv43
- [x] Commit scenario_1 results (PR #29, merged)
- [x] Commit scenario_2 results (PR #30, merged)
- [x] Launch scenario_3: `pixi run scenario-3` — launched 2026-03-28 23:40 on gpusrv43, PID 590014
- [x] ~~Monitor scenario_3~~ — OOM-killed 2026-03-29 18:05 on MIG 20GB partition. ~60% done. Locks cleared, scANVI auto-skip added to Snakefile.
- [x] Relaunch scenario_3 on full A100 80GB (2026-03-29 18:42). Pass 1 running (scvi_multi/scanorama/seurat_cca HPO).
- [x] ~~Pass 2 auto-relaunch via monitor script (PID 3944567)~~ — monitor died. R methods (fastMNN, seurat_cca, seurat_rpca) now auto-skipped for S3 (see below). sysvi reached 28/30 trials before pipeline stopped.
- [x] Investigate R method failures on S3 — root cause: R's 2^31 vector limit in intermediate computations (fastMNN `.average_correction`, Seurat `IntegrateLayers`). Full analysis: `tasks/r_methods_scalability.md`
- [x] Implement auto-skip for R methods (fastMNN, seurat_cca, seurat_rpca) when input >100K cells — added to `Snakefile:45-61`, same pattern as DESC/scANVI skips. Dry-run verified on S1 (keep) and S3 (skip).
- [x] Relaunch S3 — locks cleared, relaunched 2026-04-01 on supergpu25 (H100 80GB), PID 4108981. Snakemake DAG: sysvi HPO (2 trials) → sysvi correction → scpoli HPO → scpoli correction → aggregate → metrics → plots.
- [ ] Verify S3 completion: check `all_methods.h5ad`, metrics parquets, and results_table.pdf
- [x] Fix scANVI broadcast_labels OOM for COMPOUND-plate scenarios — added auto-skip in Snakefile when `plate_types` includes COMPOUND
- [x] Create wave2 config (`inputs/conf/scenario_wave2.json`) — sources 5,9,11, T2+COMPOUND, batch=Metadata_Source
- [ ] Launch scenario_4: `pixi run scenario-4`
- [ ] Launch scenario_5: `pixi run scenario-5`
- [ ] Verify all 15 methods complete per scenario (check `all_methods.h5ad`)
- [ ] Verify results table PDFs generated per scenario
- [ ] Compare Seurat methods on full data (failed on smoketest — should work here)
- [ ] Commit and PR scenario 3-5 results (branch: feature/scenario-3-5-results)

---

## Milestone 2: Reproduce Arevalo et al. — REFERENCE COLLECTED

Published scores documented in `tasks/arevalo_baseline.md` (per-scenario, per-metric tables for scenarios 1-3, 5; aggregate rankings).

Correlation of our default results with published results must be >0.95 for the overlapping methods.

**Overlapping methods** (9 + baseline): Harmony, Seurat CCA, Seurat RPCA, Scanorama, fastMNN, Combat, DESC, scVI, Sphering, Baseline.
**Not in our pipeline**: MNN (replaced by fastMNN — we have it), BBKNN (excluded by Arevalo too).

- [x] Obtain Arevalo et al.'s published metric scores → `tasks/arevalo_baseline.md`
- [ ] After C1 completes: run correlation analysis against published scores
- [ ] If <0.95: check hyperparameter defaults (scVI n_latent=30, Harmony clusters=300/iterations=20, DESC convergence params)
- [ ] If still off: check metric version differences (scib-metrics API changes)

---

## Milestone 3: MOA Annotation Mapping — PLANNED

Link mechanism-of-action annotations to Target2 compounds for C4 evaluation.
Detailed plan in `plans/moa_annotation_mapping.md`.

**Key finding**: Annotation files and merging code already exist in the repo. MOA-based mAP computation has been enabled (commit a841dd9).

- [x] Plan researched and written → `plans/moa_annotation_mapping.md`
- [x] Enable MOA-based mAP computation in `metrics/map.py` (commit a841dd9)
- [ ] Add `"Metadata_DRH_MOA"` to eval_keys in scenario configs 1-5
- [ ] Run end-to-end verification on one scenario
- [ ] Generate coverage statistics for paper Methods section
- [x] Investigate RefChemDB as DRH replacement — POC on TARGET2 (302 cpds). Script: `exploration/refchemdb_mapping.py`, cache: `exploration/cache/`. Bridge: InChIKey → PubChem CID → DTXSID → RefChemDB.
- **Finding**: RefChemDB is not viable as DRH replacement for TARGET2. At evidence count≥2 (expert-curated quality), only 40/302 cpds map, yielding 2 evaluable targets (4 cpds). At count≥1 (unfiltered, includes single-source), 113 cpds map with 77 evaluable. DRH at ≥3 cpds/MOA: 31 MOAs, 110 cpds. Root cause: chemical space mismatch — RefChemDB focuses on EPA tox targets (AR, ESR1, ACHE), JUMP TARGET2 is mostly kinase inhibitors/oncology.
- [ ] Explore alternative annotation sources: ChEMBL target annotations (`inputs/metadata/target_annotations.parquet`, 62% coverage of 116K cpds) or JUMP consortium's own compound metadata. These are more promising for COMPOUND-plate evaluation.

---

## Milestone 4: C2 — Launch Full HPO — NOT STARTED

All 15 methods × 200 trials each (runs in background for 1-2 weeks).

- [ ] Update `shared.json` to `optuna_trials: 200` (or per-scenario override)
- [ ] Launch on HPC with GPU allocation
- [ ] Monitor convergence

---

## Pipeline Improvements (2026-03-23)

All changes on branch `feature/migrate_to_apptainer`, 8 commits:

1. **Optuna improvements** (2c6232e): Reproducible seed (TPESampler seed=42) on all HPO scripts. Standardized CSV output via shared `save_optuna_results`. R-to-Python migration for fastMNN/Seurat HPO orchestration.
2. **GPU resource cleanup** (b761f57): Removed `nvidia_gpu=1` from 13 CPU-only rules (harmony, scanorama, combat, sphering, fastMNN, Seurat, MNN). Allows proper parallelization.
3. **MOA mAP enabled** (a841dd9): Replaced skip logic in `metrics/map.py` with annotation loading + merging, following same pattern as `run_scibmetrics_benchmarker.py`.
4. **sysVI GPU** (02364c0): Removed `JAX_PLATFORMS=cpu` from sysVI — uses same scvi.sif container as scVI which works fine on GPU.
5. **Early stopping standardized** (a3387ac): All scvi-tools models now use `validation_loss` instead of mixed `elbo_validation`/`validation_loss`. Monitors the full training objective.
6. **HPO epochs standardized** (7b00d81 → 29ecdb3): All deep learning HPO capped at 100 epochs (was 50 for scVI/scANVI/scPoli, 999999 for sysVI). Correction scripts still use 999999 + early stopping.
7. **Remaining JAX cleanup** (323f364): Removed unnecessary `JAX_PLATFORMS=cpu` from scanorama and DESC rules. Added `random_state=0` to harmony correction for reproducibility.

---

## Pipeline Fix: scpoli (2026-03-28, PR #29)

- **Root cause**: scPoli's trainer disables early stopping during pretraining. Setting `pretraining_epochs=999999` meant early stopping never activated — training ran forever.
- **Fix**: Use HPO-tuned `pretrain_to_train_ratio` to compute finite pretraining epochs from a 400-epoch budget. Training epochs remain 999999 (early stopping works during training phase).
- **Also fixed**: scpoli.def got cuDNN 9.8 fix + JAX_PLATFORMS="" (GPU enabled). HPO epochs increased 50→75. cLISI disabled in lightweight metrics (degenerate with high-cardinality labels).

---

## Bugs Fixed

1. **batch_key character splitting** (2026-03-22): `','.join(config["batch_key"])` split single-string batch keys into characters. Fixed by normalizing `batch_key` to a list in `rules/common.smk`.
2. **seurat_rpca future.globals.maxSize** (2026-03-23): RPCA creates ~1.26 GiB intermediates, exceeding 500 MiB limit. Fixed with `options(future.globals.maxSize = 4 * 1024^3)` in `run_seurat_trial.R`.
3. **sysVI HPO 999999 epochs** (2026-03-23): sysVI HPO trained to full convergence (~12 min/trial) while all other methods used 50 epochs (~1 min/trial). Fixed by standardizing to 100 epochs for all.
4. **cuDNN 9.1 vs 9.8 in containers** (2026-03-27): `jax[cuda12]` and `torch` pin `nvidia-cudnn-cu12==9.1.0` but JAX 0.7+ needs ≥9.8.0. Fixed by adding `pip install --force-reinstall --no-deps "nvidia-cudnn-cu12==9.8.0.87"` as final step in all GPU container defs.
5. **coarsen_labels ZeroDivisionError** (2026-03-27): `SEMISUP_MIN_BATCHES=5` caused all labels to be marked "Unknown" when scenarios had <5 batches (scenario 2 has 3). scANVI crashed with `n_labels=0`. Fixed with dynamic `min_batches = max(3, min(5, n_batches))` in `scripts/utils.py`.
6. **scANVI HPO redundant scVI training** (2026-03-27): `optimise_scanvi.py` retrained scVI from scratch every trial despite fixed architecture. Refactored to train scVI once and reuse across all 30 trials (~2x speedup).
7. **mAP numpy array indexing** (2026-03-27): `metrics/map.py` called `.loc[]` on numpy arrays from `adata.obsm`. Fixed with `isinstance` check and `get_indexer` for positional indexing.
8. **scpoli infinite pretraining** (2026-03-28): `pretraining_epochs=999999` defeated early stopping (disabled during pretraining by scPoli trainer). Fixed with finite pretraining epochs from HPO ratio × 400 budget. PR #29.

---

## Scenario 8 Blockers (2026-03-25, priority-ranked)

### P1: scvi/scanvi CUDA OOM bug
- **Blocks**: scvi_single, scvi_multi, scanvi_single, scanvi_multi (4 methods × 2 scenarios)
- **Symptom**: `torch.OutOfMemoryError: Tried to allocate 6477.11 GiB` on A100 80GB
- **Root cause**: Not a memory issue — 6.4 TiB allocation is a code/data bug (likely label encoding creating a massive one-hot matrix in `broadcast_labels` during loss computation)
- **Status**: FIXED (2026-03-25)
- [x] Reproduce and diagnose the allocation source
- [x] Fix scVI: removed `labels_key` from `SCVI.setup_anndata()` in `correct_with_scvi.py` and `optimise_scvi.py` (unsupervised model, doesn't need it)
- [x] Fix scANVI: added auto-skip in `Snakefile` when `label_key` has >10K categories (architectural limit: `broadcast_labels` is O(batch × n_labels × latent_dim) per forward pass)
- [ ] Re-run scvi_single, scvi_multi on full GPU (not MIG)

### P2: Clear stale locks + launch missing defaults methods
- **Blocks**: sysvi, scpoli, scvi_multi in scenario_8_defaults (never ran, no errors — just weren't triggered)
- **Action**: Clear `.snakemake/locks/`, re-launch pipeline
- [ ] Clear locks
- [ ] Launch sysvi, scpoli, scvi_multi for scenario_8_defaults

### P3: Re-launch scenario_8 HPO pipeline
- **Blocks**: harmony_v1, harmony_v2 (optimized, correction not run), seurat_cca (never started), all deep learning methods pending P1 fix
- **Action**: Re-launch after clearing locks; harmony/seurat_cca should just work
- [ ] Re-launch scenario_8 pipeline
- [ ] Verify harmony_v1, harmony_v2, seurat_cca corrections complete

### P4: fastMNN R optimization error (scenario_8 only)
- **Blocks**: 1 method in scenario_8
- **Symptom**: `Error: Optimization did not complete successfully`
- **Note**: fastMNN works fine in scenario_8_defaults (using default params)
- [ ] Check R optimization log for root cause
- [ ] Fix and re-run

### P5: DESC 200K cell limit
- **Blocks**: 1 method in both scenarios
- **Symptom**: `ValueError: operands could not be broadcast together with shapes (370231,) (200000,)`
- **Root cause**: DESC library hardcodes 200K cell subsampling limit; scenario_8 has 370K cells
- **Options**: Skip DESC for scenario_8 (it underperforms anyway) or patch subsampling
- [ ] Decide: skip or patch

### Scenario 8 metrics analysis finding (2026-03-25)
- For fast tuning feedback: **graph_conn + PCR_batch** (Adj R²=0.89 vs mean_overall, both <1s per method)
- Best 2-metric predictor overall: **kmeans_nmi + pcr_comparison** (Adj R²=0.92), but kmeans is slow (~7 min on large embeddings)
- Best 4-metric near-perfect proxy: pcr_comparison + iLISI + silhouette_label + kmeans_nmi (Adj R²=0.98)

---

## P0 Blocker: Container Rebuild (2026-03-26) — RESOLVED 2026-03-27

Three containers had stale .sif images. All rebuilt on gpusrv40 using `/localscratch` as TMPDIR.

**Root cause**: cuDNN version mismatch. Apptainer `--nv` injects host cuDNN 8.9.7 (`.so.8`), but JAX inside containers was compiled against cuDNN 9.8.0. Fix: prepend pip-installed cuDNN path in `%environment` so container's own cuDNN takes priority over host's.

**Fix applied to .def files** (not yet built):
- `containers/scvi.def` — added `LD_LIBRARY_PATH=/usr/local/lib/python3.11/dist-packages/nvidia/cudnn/lib:$LD_LIBRARY_PATH`
- `containers/harmony_v1.def` — same fix
- `containers/desc.def` — .def updated Mar 24, .sif stale from Mar 22

- [x] ~~Rebuild scvi.sif, harmony_v1.sif, desc.sif on a node with stable I/O (NOT gpusrv60)~~ — rebuilt 2026-03-27 on gpusrv40 using /localscratch as TMPDIR
- [x] ~~Consider using local `/tmp` if node has enough space + exec permissions~~ — used /localscratch (3.5TB, exec-mounted)
- [x] ~~After rebuild: re-launch scenario_2 GPU methods~~ — completed 2026-03-27, all 15 methods done
- [ ] After rebuild: re-launch scenario_8 GPU methods

---

## Milestone 5: Extended Scenarios (wave1, wave2, reference mapping) — PLANNED

Research completed 2026-03-29. Full documentation in:
- `tasks/jump_metadata_research.md` — JUMP-CP metadata, microscopes, waves, compound replication
- `tasks/scenario_registry.md` — structured scenario cards with challenges, use cases, replication stats

### Key new scenarios identified

| Scenario | Sources | Key property | Status |
|----------|---------|-------------|--------|
| **wave1** | S1,2,3,6,8,10,15 (7 src) | Full Wave 1, 82K cpds, median 5 reps, 74.5% in 5+ sources | Not started |
| **wave2** | S5,9,11 (3 src) | Full Wave 2, 36K cpds, median 7 reps, 84.8% in all 3 sources, 3 diff microscopes + mixed plate format | Not started |
| **A1** | S2,5,6,10 (4 src, CV8000) | Cross-wave, same microscope — isolates wave effect | Not started |
| **A3** | S2,5,6,10 (4 src, CV8000) | Cross-wave + compounds, best same-microscope reference | Not started |
| **C4** | S2,3,6,8,10,15 (6 src) | W1 384-well only, hierarchical batch, general reference atlas | Not started |
| **C8** | S2,3,6,7,10 (5 src) | W1+bridge, best within-scenario evaluation (6.3% replication) | Not started |
| **F4** | C4 → wave2 | Reference mapping: W1 reference, W2 query, 5,392 evaluable cpds | Not started |
| **B7** | S1,3,9 (Opera Phenix) | Plate format isolation (384 vs 1536) | Not started |
| **E3** | S2,3,6,8,10 | Microscope type as batch key (novel methodology) | Not started |

### Tasks

- [x] Research JUMP-CP metadata: microscopes, plate formats, waves, compound overlap
- [x] Analyse compound replication structure (compound_source.csv, well.csv)
- [x] Document Wave 1 vs Wave 2 design differences
- [x] Create scenario registry with structured cards → `tasks/scenario_registry.md`
- [x] Create scenario config JSON for wave2 → `inputs/conf/scenario_wave2.json`
- [ ] Create scenario config JSONs for wave1 and other priority new scenarios
- [ ] Add pixi task for wave2 to `pixi.toml` (or run manually with snakemake --configfile)
- [ ] Run wave2 scenario (3 sources, quick — good first test of new scenarios)
- [ ] Run wave1 scenario (7 sources, large — main reference atlas candidate)
- [ ] Design and run F4 reference mapping experiment (requires wave1 + wave2 embeddings)

---

## Open Questions

- `r_hpo.sif` was never built but turned out unnecessary — fastMNN/Seurat HPO ran via `r.sif`. Should the def file and build task be removed?
- Scenario 7 results table PDFs are stale (Mar 2, before new methods). Regenerate?
- Should scenarios 6-8 and apricot also be included in C1, or only 1-5?
- Pipeline improvements branch needs PR to main after scenario_1 validates successfully.
- wave1 includes source_1 (1536-well). Should we run both wave1 (with S1) and C4 (without S1, 384-only)?
- wave2 source_9 is 1536-well. Should B7 (plate format test) be run first to validate 1536-well compatibility?
