# 12-Hour Reference Mapping Execution Plan

Last updated: 2026-04-20  
Branch: refmap/model-save

---

## Objectives (PI meeting)

1. **Integration quality table** — for each (query_source, method, arm): batch_asw, kBET, iLISI, precision@10
2. **T-arm contrast** — does removing TARGET2 anchors (T+ → T- / T-_matched) hurt integration?
3. **UMAPs** — query projected into reference space; reference + query shuffled before plotting, colored by source and by reference/query split

---

## Experiments

### Exp 1 — LOSO within scenario_5 (pure batch correction test)

- **Reference**: scenario_5 minus source_8 (sources 2, 3, 6, 10)
- **Query**: source_8 (same JUMP compound library, ~100% overlap)
- **Purpose**: isolates batch effect from compound novelty; tests whether each method can re-integrate a held-out source

### Exp 2 — Partial compound overlap (reference mapping proper)

- **Reference**: full scenario_5 (all 5 sources)
- **Query**: source_5 (Wave 2, ~7% compound overlap: 292 T2 + ~2,117 COMPOUND shared)
- **Purpose**: realistic reference mapping — query profiled a largely novel compound set

---

## Execution matrix

| Exp | Query | Atlas | Method | T+ | T- | T-_matched |
|-----|-------|-------|--------|----|----|------------|
| 1 | source_8 | excl_source_8 | Symphony | ✓ mapped | todo | todo |
| 1 | source_8 | excl_source_8 | scVI | todo (atlas missing) | todo | todo |
| 1 | source_8 | excl_source_8 | scPoli | todo (atlas missing) | todo | todo |
| 2 | source_5 | scenario_5_full | Symphony | ✓ mapped | ✓ mapped | ✓ mapped |
| 2 | source_5 | scenario_5_full | scVI | todo (atlas broken*) | todo | todo |
| 2 | source_5 | scenario_5_full | scPoli | todo (atlas missing) | todo | todo |

*scvi_atlas_scenario_5_full/ only has X_min.json — model save failed. Re-run nb10.

---

## Metrics reported per row

| Metric | Higher = | Interpretation |
|--------|----------|----------------|
| `batch_asw` | better mixing | Are query cells interleaved with reference in joint embedding? |
| `kbet` | better mixing | kNN neighborhoods are batch-diverse |
| `ilisi` | better mixing | Inverse LISI; 1 = pure, n_batches = fully mixed |
| `precision@10` | better mapping | Fraction of top-10 reference neighbors with same compound ID |
| `compound_asw` | more coherent | Same-compound query cells cluster; expected near 0 (singletons) |

**Primary metrics for PI presentation**: batch_asw + kbet/ilisi (batch mixing) and precision@10 (compound recovery).

---

## Pre-execution fixes

### Fix 1 — scPoli prototype explosion (CRITICAL, in nb20)

**Problem**: `cell_type_keys=LABEL_KEY` passes all ~32 K unique JCP2022 IDs → ~32 K prototypes → CPU hang (same root cause as the S5 scPoli investigation).

**Fix**: add `T2_ONLY_PROTOTYPES=True` parameter.  
- Create `obs["_cell_type"]`: T2 compound ID for T2 wells, `"BACKGROUND"` for all others.
- Use `cell_type_keys="_cell_type"` → ~293 prototypes (292 T2 + 1 background group).
- Save `training_config.json` to atlas dir so nb21 knows how to set `labeled_indices`.

### Fix 2 — nb21 labeled_indices when T2-only prototypes used

**Problem**: nb21 currently passes ALL compound-matching wells as `labeled_indices`. With T2-only prototypes, only T2-matching wells should be labeled; passing COMPOUND wells confuses the prototype loss.

**Fix**: nb21 reads `training_config.json`; if `t2_only_prototypes=True`, restricts `labeled_indices` to T2-matching wells only.

### Fix 3 — nb41 paradigm order

**Problem**: nb41 still lists `"sysvi"` in `paradigm_order`. Replace with `"scvi"`.  
Also extend to 5-metric layout (add kbet, ilisi, precision_at_10 panels).

### Fix 4 — stale metrics CSV

**Problem**: `results/acute/source_5_metrics.csv` has 3 duplicated/inconsistent rows.  
**Fix**: `run_all.sh` deletes both metric CSVs at startup and recomputes from scratch.

---

## Scripts to create / modify

| File | Action | Purpose |
|------|--------|---------|
| `scripts/run_all.sh` | CREATE | Full 12h orchestrator |
| `scripts/generate_umaps.py` | CREATE | UMAP per (query, method, T_arm=tplus); shuffled z-order |
| `scripts/generate_report.py` | CREATE | Compile CSVs → pivot table + bar charts |
| `notebooks/20_train_atlas_scpoli.py` | PATCH | T2-only prototypes |
| `notebooks/21_map_scpoli.py` | PATCH | Read training_config.json for labeled_indices |
| `notebooks/41_plot_acute.py` | PATCH | scvi paradigm, 5-metric layout |

---

## Execution order and parallelization

```
t=0   Background: Symphony Exp1 remaining + all metrics reset + recompute [CPU]
t=0   GPU track:  scVI atlas excl_source_8 training (~30-45 min)
t=1h  GPU track:  scVI map source_8 × 3 arms (~15 min total)
t=1h  GPU track:  scVI atlas scenario_5_full retraining (~45-60 min)
t=2h  GPU track:  scVI map source_5 × 3 arms (~15 min total)
t=2h  CPU track:  scVI metrics × 6 rows
t=3h  GPU track:  scPoli atlas excl_source_8 training (~2-4h — slowest step)
t=6h  GPU track:  scPoli map source_8 × 3 arms (~30 min total)
t=7h  GPU track:  scPoli atlas scenario_5_full training (~3-5h)
t=11h GPU track:  scPoli map source_5 × 3 arms (~30 min)
t=11h CPU track:  scPoli metrics × 6 rows
t=11h CPU track:  UMAP generation × 6 figures (one per query × method)
t=12h CPU track:  generate_report.py → pivot table + bar charts + LaTeX table
```

Symphony runs fully in background (CPU-only, doesn't compete for GPU).

---

## Expected output files

```
results/acute/source_8_metrics.csv   # 9 rows: 3 methods × 3 arms
results/acute/source_5_metrics.csv   # 9 rows: 3 methods × 3 arms
results/source_8_summary.csv         # pivot table: method × arm × metric
results/source_5_summary.csv
results/source_8_barplot.pdf
results/source_5_barplot.pdf
results/umaps/source_8_symphony_tplus.png
results/umaps/source_8_scvi_tplus.png
results/umaps/source_8_scpoli_tplus.png
results/umaps/source_5_symphony_tplus.png
results/umaps/source_5_scvi_tplus.png
results/umaps/source_5_scpoli_tplus.png
```

---

## Key scientific questions answered

| Question | Metric | T-arm comparison |
|----------|--------|-----------------|
| Does the method integrate the query? | kbet, ilisi | T+ (best case) |
| Do T2 anchors improve integration? | batch_asw, precision@10 | T+ vs T-_matched |
| Is the effect size vs no-anchor? | all | T-_matched vs T- |
| Can methods map novel compound space? | precision@10, batch_asw | Exp2 T+ |
