# Reconcile + commit refmap work (PR refmap/model-save → main)

**Last updated**: 2026-05-09
**Branch**: `refmap/model-save` (2 commits ahead of main: 28c9385 scaffold, 5ecdafa model-save)
**Goal**: Land all reference-mapping experiment work on main as a clean, reviewable PR — without churning unrelated repo state.

## What's actually on disk (audit)

### Modified, tracked files (16)

**Main pipeline (4 files)**
- `scripts/utils.py` — `coarsen_labels` floor lowered `max(3, ...)` → `max(2, ...)` (S5 scPoli investigation; **keep**)
- `scripts/correct_with_scpoli.py` — adds `--min_batches` CLI arg + threading (**keep**)
- `pixi.toml` — adds `scenario-source5` task (preprocessing-only target for the source_5 query parquet) (**keep**)
- `tasks/lessons.md` — 9 new lessons spanning S3 R-scalability through 2026-04-19 sysVI/scPoli/symphony refmap findings (**keep, already current**)
- `tasks/todo.md` — last touched 2026-04-19, **stale** w.r.t. Exp 2 results from 2026-04-21 (**update before commit**)

**Reference mapping (12 files)** — all tracked .py notebooks + `src/metrics.py` + refmap pixi
- `reference_mapping/notebooks/{00,01,10_scvi,11_scvi,20,21,30,31,40,41}.py` — all touched during Apr 19–21 work
- `reference_mapping/src/metrics.py` — refactored metric helpers
- `reference_mapping/pixi.toml` — channel cleanup (drop `pytorch`/`nvidia`, use conda-forge `pytorch-gpu`)

### Untracked, refmap (worth committing)

- `.ipynb` siblings of every `.py` (jupytext-paired) — 14 notebooks
  - **size flag**: `40_metrics_tminus.ipynb` (322K), `40_metrics_tplus.ipynb` (322K), `10_train_atlas_scvi.ipynb` (183K) carry baked-in outputs. Strip before commit.
  - all others <15K — fine as-is
- `01_prepare_query_arms_s5feats.ipynb` (Apr 21) — newer variant that reuses scenario_5 feature space (no separate wave2 preprocessing); pairs with the matching `.py` once jupytext-paired
- `10_train_atlas_sysvi.{py,ipynb}`, `11_map_sysvi.{py,ipynb}` — sysVI scaffolding. Lesson confirms architectural limit (cannot map new batches). **Keep with README header note.**
- `reference_mapping/plans/f4_wave2_mapping.md` (12K), `12h_execution_plan.md` (6K) — planning docs
- `reference_mapping/scripts/` — 9 ad-hoc orchestration scripts (cleanup needed; see below)
- `reference_mapping/results/source_5_summary.csv`, `source_5_barplot.pdf`, `acute/source_5_metrics.csv`, `acute/source_5_plot.pdf`, `precision_at_k.csv`, `umaps/*.png` — the actual Exp 2 deliverables

### Untracked, NOT refmap (handle separately, do not bundle in this PR)

- `inputs/conf/scenario_c4.json` (Apr 19, Wave 1 reference atlas config) — **not refmap**; defer to a separate Wave 1 PR
- `inputs/conf/scenario_source5.json` (Apr 21) — **needed by refmap** (paired with the `scenario-source5` pixi task); include
- `SCENARIO8_LAUNCH_PLAN.md` (Mar 23) — stale top-level doc, **out of scope**; leave for a separate cleanup pass
- `scGPT_human/` (198 MB model weights) — **gitignore**, do not commit
- `exploration/cache/` (90 MB, RefChemDB POC) — already covered by file-extension gitignores; left alone
- `.claude/`, `.keras/`, `.nv/` — local caches, already gitignored or should be
- `reference_mapping/logs/` (run logs) — already covered by `*.log` gitignore but the directory itself is untracked; either gitkeep or skip

## Steps (in dependency order)

### Step 1 — Stop tracking large model weights and stale top-level files

**Files**: `.gitignore`

- Append patterns: `scGPT_human/`, `reference_mapping/papermill_runs/`, `reference_mapping/logs/` (the dir, even though *.log already covers contents).
- Verify with `git status --short` that none of these appear afterwards.

### Step 2 — Consolidate scPoli reference-embedding extractor (one-off cleanup)

**Files**: `reference_mapping/scripts/extract_scpoli_ref_embedding.py`, `reference_mapping/scripts/gen_scpoli_ref_emb.py`

These two scripts do the same thing (extract latent embedding from a trained scPoli model into `reference.h5ad`). The newer `gen_scpoli_ref_emb.py` (Apr 21) supersedes `extract_scpoli_ref_embedding.py` (Apr 19).

- Diff them to confirm gen_scpoli_ref_emb.py is the keeper (better paths handling).
- Delete `extract_scpoli_ref_embedding.py`.
- Rename `gen_scpoli_ref_emb.py` → `extract_scpoli_ref_embedding.py` (more readable name) **OR** keep as-is. Pick one; don't keep both.

### Step 3 — Add sysVI README header note

**File**: new file `reference_mapping/notebooks/SYSVI_NOTE.md` (or top-cell markdown in each sysVI notebook)

Short note (≤15 lines) explaining:
- sysVI is in this repo as scaffolding only.
- It cannot do reference mapping with new (unseen) batches: `SysVI.load_query_data` hard-codes `transfer_batch=False`.
- Use scPoli or Symphony for refmap with new batches; sysVI is appropriate only for in-sample batch correction.
- Cross-link to the lesson in `tasks/lessons.md` for context.

Decision point: separate `.md` file vs. top-cell markdown in `10_train_atlas_sysvi.ipynb` and `11_map_sysvi.ipynb`. **Recommendation**: top-cell markdown in both notebooks — closer to the code, harder to miss when opening the notebook.

### Step 4 — Strip notebook outputs for the three heavy notebooks

**Files**: `reference_mapping/notebooks/40_metrics_tminus.ipynb`, `40_metrics_tplus.ipynb`, `10_train_atlas_scvi.ipynb`

These are 322K/322K/183K vs. <15K for everything else — clear sign of baked-in cell outputs. Strip with `jupyter nbconvert --clear-output --inplace <files>`. Keep code, drop outputs.

Lighter notebooks may also have outputs but the size cost is low; defer unless they balloon a diff.

### Step 5 — Update `tasks/todo.md` to reflect Exp 2 completion

**File**: `tasks/todo.md`

Reconcile the "Reference Mapping Experiment Design" section against the on-disk results:

- Mark **DONE**: nb20 (scPoli atlas), nb21 ×3 (T+/T-/T-_matched), nb31 ×2 (Symphony T-/T-_matched), nb40 ×6 (metrics), Exp 2 source_5 row in execution matrix.
- Mark **abandoned**: S5 scPoli mb=3/mb=4 investigation (no parquet produced; lesson captured).
- Mark **decided**: sysVI dropped from refmap (lesson + sysVI note). Update the open "Decide on nb10/11 (sysVI)" item to `[x]` with note "dropped — see SYSVI_NOTE / lesson".
- Add a new section "Exp 2 results (source_5)" with the headline numbers from `source_5_summary.csv`:
  - scPoli wins compound_asw (0.38) and batch_asw (0.90)
  - scVI wins iLISI (0.36)
  - **T+/T-/T-_matched arms differ by <1% across all metrics** — TARGET2 anchor ablation appears null on this query
- Carry forward the still-current items: Exp 1 (source_8 LOSO) and a confirmation run on source_9 or source_11.

### Step 6 — Stage and commit in three logical commits (don't bundle)

Each commit limited to a coherent scope so review is feasible.

**Commit A — main pipeline plumbing**
```
scripts/utils.py
scripts/correct_with_scpoli.py
pixi.toml
inputs/conf/scenario_source5.json
```
Message: "scPoli min_batches CLI + scenario-source5 preprocessing target"
Body: 1-paragraph rationale referencing the S5 investigation lesson.

**Commit B — refmap notebooks + scripts + plans**
```
reference_mapping/**/*.py
reference_mapping/**/*.ipynb (output-stripped)
reference_mapping/src/metrics.py
reference_mapping/pixi.toml
reference_mapping/plans/*.md
reference_mapping/scripts/* (after step 2 dedup)
reference_mapping/notebooks/SYSVI_NOTE.md (or in-notebook markdown)
```
Message: "Reference mapping Exp 2 (source_5 → scenario_5): scPoli, Symphony, scVI"

**Commit C — docs (lessons + todo)**
```
tasks/lessons.md
tasks/todo.md
```
Message: "Update lessons + todo for refmap Exp 2 + scPoli S5 investigation"

**Commit D — gitignore housekeeping**
```
.gitignore
```
Message: "gitignore scGPT model weights + refmap papermill_runs/logs"

Order matters: D first (so ignored files are filtered cleanly), then A, B, C. Each commit verifiable independently.

### Step 7 — Pre-PR sanity checks

Before pushing:
- `git status --short` — should show only intentional untracked items (`SCENARIO8_LAUNCH_PLAN.md`, `inputs/conf/scenario_c4.json`, `exploration/cache/`, `scGPT_human/`).
- `pixi run -e default python -c "from scripts.correct_with_scpoli import correct_with_scpoli"` — confirm the new arg doesn't break import.
- `git log --oneline main..HEAD` — confirm 4 new commits + the existing 2 (28c9385, 5ecdafa) = 6 total.
- Open the diff in `gh pr diff` (after step 8) to spot unintended changes.

### Step 8 — Open PR

```
gh pr create --base main --head refmap/model-save \
  --title "Reference mapping experiments + scPoli min_batches plumbing" \
  --body "$(cat <<'EOF'
## Summary
- Reference mapping ad-hoc experiment scaffold (`reference_mapping/`)
- Experiment 2 results: source_5 query → scenario_5 reference, three paradigms (scPoli, Symphony, scVI) × three arms (T+, T-, T-_matched)
- scPoli `--min_batches` plumbing in main pipeline (S5 investigation; mb=2 floor in `coarsen_labels`)
- New `scenario-source5` pixi task for preprocessing the wave2 query parquet

## Key findings
- scPoli wins compound_asw (0.38) and batch_asw (0.90) on Exp 2
- scVI wins iLISI (0.36)
- T+/T-/T-_matched arms differ by <1% across all metrics — TARGET2 anchor ablation appears null on this query (needs confirmation on a second query source)

## Caveats
- sysVI scaffolding included but explicitly noted as inapplicable for refmap-with-new-batches (architectural limit of `SysVI.load_query_data`)
- S5 scPoli mb=3/mb=4 investigation abandoned (lesson captured; the CLI plumbing remains for future retries)

## Test plan
- [ ] Verify `pixi run scenario-source5` builds the source_5 query parquet
- [ ] Re-run `01_prepare_query_arms_s5feats.ipynb` against the committed parquet
- [ ] Confirm `tasks/todo.md` reflects on-disk state
EOF
)"
```

## Critique (CLAUDE.md Step 3)

**Public API impact**:
- `scripts/correct_with_scpoli.py` gains a new `--min_batches` CLI flag with a default that matches existing behavior (`SEMISUP_MIN_BATCHES=5`). Pre-existing Snakemake rules call this script without `--min_batches`, so default-path users see no change. Low risk.
- `scripts/utils.py` `coarsen_labels` floor change `3 → 2` is a behavior change for scenarios with ≤3 batches and labels appearing in exactly 2 batches. **In practice**: scenario_2 has 3 batches; the change keeps labels with `min_batches=min(5,3)=3`. The floor of 2 only matters when a caller explicitly passes `min_batches=2`. Effective no-op for the existing pipeline.
- New `scenario-source5` pixi task is additive. No-op for unrelated runs.

**Edge cases**:
- Notebook output stripping: confirm jupytext pairing is preserved after `nbconvert --clear-output --inplace`. If pairing breaks, rerun `jupytext --sync` afterwards.
- The `.gitignore` addition for `scGPT_human/` is post-hoc — verify the dir is currently untracked (no need for `git rm`).

**Simpler path considered**:
- "Just commit everything in one mega-commit" — rejected; review burden makes it likely to be rubber-stamped or deferred indefinitely. The 4-commit split keeps each review unit <300 lines.
- "Rebase the existing 2 commits + new work into 1 squashed PR" — rejected; the 2 prior commits already have clean scope, no need to rewrite history.
- "Skip the sysVI note, just delete the files" — rejected per user decision (keep with README note); the architectural lesson is more valuable as scaffolding-with-warning than as a missing file.

**Risks**:
- The .ipynb files (output-stripped) will still register as binary diffs in PR review. GitHub renders them, but reviewers may struggle. Mitigation: in PR body, point reviewers at the `.py` siblings.
- `reference_mapping/results/*.png` are gitignored by `*png` rule — the UMAP results won't be on the PR. The acute PDFs and CSVs we want will. Confirm the right artifacts land.

**Out of scope (deferred)**:
- Exp 1 (source_8 LOSO) — execution incomplete per `12h_execution_plan.md`; ship Exp 2 results now, follow-up PR for Exp 1.
- Wave 1 reference atlas (`scenario_c4.json`) — separate PR.
- `SCENARIO8_LAUNCH_PLAN.md` cleanup — separate PR.

## Test strategy

- **No new pytest tests proposed**: refmap is ad-hoc experiment code, not part of the pipeline-under-test surface. Existing tests still pass (no main-pipeline behavior change).
- **Manual verification before commit**:
  1. `pixi run -e default python scripts/correct_with_scpoli.py --help` — confirms the new flag is wired
  2. Open one stripped notebook (`40_metrics_tplus.ipynb`) — confirms cells still execute
  3. `pixi run scenario-source5 --dry-run` — confirms the new pixi task resolves
- **Reproducibility check**: with the committed `scenario_source5.json` and the orchestration scripts, a fresh clone should be able to re-run `bash reference_mapping/scripts/run_remap_s5feats.sh` and reproduce the source_5 results.

## Estimated effort

- Step 1 (gitignore): 5 min
- Step 2 (script dedup): 10 min
- Step 3 (sysVI note): 15 min
- Step 4 (strip outputs): 10 min
- Step 5 (todo.md update): 30 min
- Step 6 (4 commits): 15 min
- Step 7 (sanity checks): 10 min
- Step 8 (PR open): 10 min

**Total**: ~1.5 h focused work.
