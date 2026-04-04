# Lessons Learned

## 2026-03-23: Scenario 8 Launch

### DESC 200K cell limit
- DESC's `pretrain()` subsamples to 200K cells for Louvain init, but `fit_on_all()` predicts on full data → shape mismatch crash when N > 200K.
- Arevalo et al. never hit this (their scenarios had <200K cells). Scenario 8 has 370K.
- Fix: Snakefile auto-skips DESC when preprocessed parquet exceeds 200K rows.
- The bug is in the upstream DESC library (`network.py` line ~343), not our code.

### Seurat future.globals.maxSize
- Seurat's `IntegrateLayers` exports ~42 GB of globals for scenario 8, exceeding the default 32 GB limit.
- Fix: set `options(future.globals.maxSize = +Inf)` in `correct_with_seurat.R`. Safe because node RAM (~690 GB) is the actual constraint.

### fastMNN prop_k NULL handling in R
- Empty CSV cells are read as `NA` in R, not `NULL`. The `!is.null()` check passes for `NA`, so the value gets forwarded to `fastMNN()` causing a crash.
- Fix: added `if (is.na(prop_k_val)) prop_k_val = NULL` after reading the CSV.

### MIG vs full GPU nodes
- Some gpu_p nodes have MIG partitions (e.g., `A100 MIG 3g.20gb` = 20 GB), not full 40 GB A100s.
- scvi-tools methods on 370K cells need the full GPU. Always verify with `nvidia-smi` before launching large jobs.

### nvidia_gpu Snakemake resource must match actual hardware
- pixi tasks had `nvidia_gpu=8` but nodes have 1 GPU each → Snakemake would try 8 concurrent GPU jobs.
- Fix: changed to `nvidia_gpu=1` in pixi.toml for scenario-8 tasks.

### Default parameter provenance
- The "defaults" CSVs were a mix of Arevalo explicit choices, library defaults, and incorrect values.
- Correct priority: (1) Arevalo explicit choice → use that, (2) Arevalo used library defaults → use library defaults, (3) method not in Arevalo → use library defaults.
- Key corrections made:
  - scVI: n_layers 1→2 (Arevalo explicit)
  - scANVI: n_epochs_kl_warmup 10→400 (library default; not in Arevalo)
  - sysVI: n_hidden 128→256, n_latent 30→15, prior standard_normal→vamp, z_distance_cycle_weight 0→2.0 (library defaults; not in Arevalo)
  - scPoli: alpha_epoch_anneal 200→100, eta 0.5→1.0, latent_dim 64→10 (library defaults; not in Arevalo)
  - fastMNN: prop_k 0.05→NULL (library default; Arevalo didn't set it)

### cuDNN version mismatch with Apptainer --nv (two layers)
- **Layer 1 (host vs container)**: Apptainer `--nv` injects the host's cuDNN 8.x, shadowing pip-installed 9.x. Fix: prepend pip-installed cuDNN path in `%environment`: `export LD_LIBRARY_PATH=/usr/local/lib/python3.11/dist-packages/nvidia/cudnn/lib:$LD_LIBRARY_PATH`
- **Layer 2 (pip resolver)**: `jax[cuda12]` and `torch` both pin `nvidia-cudnn-cu12==9.1.0` via their deps, but JAX 0.7+ needs cuDNN >=9.8.0. The pip resolver doesn't catch this. Fix: add a **final** `pip install --force-reinstall --no-deps "nvidia-cudnn-cu12==9.8.0.87"` step in the Dockerfile/def, after all other installs, so nothing can downgrade it.
- The `lightweight.sif` container had both fixes — it was the working template.

### Lustre instability kills container builds
- Building Apptainer containers with `APPTAINER_TMPDIR` on Lustre is unreliable under load. Concurrent builds cause metadata storms → `[Errno 5] Input/output error` and client evictions.
- Sequential builds on unstable Lustre nodes (like gpusrv60 on 2026-03-26) also fail during `unsquashfs` or `pip install`.
- `/dev/shm` is mounted `noexec` on this cluster — cannot be used as TMPDIR for container builds.
- Prefer nodes with stable local `/tmp` (SSD-backed, exec-mounted) for container builds.

### nvidia_gpu resource must be 1 for scenarios 1-5
- pixi tasks had `nvidia_gpu=8` → Snakemake scheduled 8 concurrent GPU jobs on nodes with 1 GPU, causing contention/crashes.
- Fixed to `nvidia_gpu=1` in pixi.toml for all scenario tasks.

### Container builds: use /localscratch on gpusrv40
- gpusrv40 has `/localscratch` (3.5 TB, XFS, rw, exec-mounted) — ideal for container builds.
- `/tmp` on gpu nodes is only 10 GB (too small for GPU containers with PyTorch+JAX).
- Successfully built all 3 containers (scvi, harmony_v1, desc) sequentially on gpusrv40.

### coarsen_labels must adapt to scenario batch count
- `SEMISUP_MIN_BATCHES=5` is a static threshold but scenarios have 1-5 batches.
- With 3 batches, no label can appear in ≥5 batches → all labels coarsened → `n_labels=0` → ZeroDivisionError in scANVI.
- Fix: `min_batches = max(3, min(SEMISUP_MIN_BATCHES, n_batches))`.
- The absolute floor of 3 was a user decision to ensure some cross-batch signal.

### scANVI HPO: train scVI once, not per trial
- The scVI base model architecture is fixed from prior HPO. Only scANVI-specific params (classification_ratio, n_epochs_kl_warmup, linear_classifier) are tuned.
- Training scVI 30 times with identical params is wasted compute. Train once, call `SCANVI.from_scvi_model()` per trial.
- Note: `linear_classifier` is an architecture param for `from_scvi_model()`, not a training param, so the scANVI model must still be created per trial.

### Snakemake incomplete markers survive process kills
- Killing snakemake (even with SIGKILL) can leave `.snakemake/incomplete/` markers that block future runs.
- Markers are base64-encoded filenames. Decode with `basename <file> | base64 -d`.
- Removing the marker file directly is simpler than `snakemake --cleanup-metadata` (which has bugs in newer versions).
- Alternative: `--rerun-incomplete` flag, but it re-runs the job from scratch.

### Parallel execution strategy
- Run HPO (scenario-8) and defaults (scenario-8-defaults) on separate GPU nodes via separate SLURM sessions.
- They write to separate output dirs so no conflicts.
- Preprocessing symlinks avoid duplicating 30+ GB of data.

### scPoli disables early stopping during pretraining (2026-03-28)
- scPoli's `scPoliTrainer` explicitly sets `self.use_early_stopping = False` when `self.epoch < self.pretraining_epochs`.
- Setting `pretraining_epochs=999999` means early stopping never activates — the model trains forever until killed.
- Fix: use the HPO-tuned `pretrain_to_train_ratio` to compute finite pretraining epochs from a fixed budget (400 epochs). Only `n_train_epochs` should be 999999 (early stopping works during training).
- This is specific to scarches/scPoli. scVI/scANVI/sysVI don't have a separate pretraining phase, so 999999 epochs + early stopping works fine for them.

### cLISI is degenerate with high-cardinality labels (2026-03-28)
- With ~302 compound labels, random embeddings score ~0.97 on cLISI, masking real signal.
- Replaced `isolated_labels=False` with `clisi_knn=False` in lightweight HPO metrics.
- cLISI remains enabled in the full scibmetrics benchmarker (for reporting), just not in HPO steering.

### scANVI broadcast_labels is incompatible with COMPOUND-plate scenarios (2026-03-29)
- scANVI's loss function calls `broadcast_labels(z1, n_broadcast=self.n_labels)` which creates tensors of O(batch × n_labels × latent_dim).
- With COMPOUND plates, even after `coarsen_labels` (keeping labels in ≥3 sources), the label count remains in the thousands → allocation requests of 555+ GiB.
- This is an architectural limit of scvi-tools' scANVI implementation, not a memory issue. It fails on any GPU.
- Fix: auto-skip `scanvi_single` and `scanvi_multi` in the Snakefile when `plate_types` includes `COMPOUND`.
- scANVI still works for TARGET2-only scenarios (302 labels).
- The previous comment in the Snakefile ("No Snakefile-level skip is needed") was wrong.

### Snakemake lock cleanup (2026-03-28)
- Removing lock files manually causes `FileNotFoundError` on next run — Snakemake expects the locks directory and files to exist.
- Correct approach: `rm -f .snakemake/locks/0.input.lock .snakemake/locks/0.output.lock` then `mkdir -p .snakemake/locks` before relaunching.
- Alternative: use `snakemake --unlock` but this sometimes stalls if the DAG is complex.

### Wave 1 vs Wave 2 compound replication is fundamentally different (2026-03-29)
- **Wave 1** (7 partners): Each partner exchanges with 4 of 6 others. Per-source replication is median 1 well/compound. Cross-source replication comes from exchange — 74.5% of compounds in 5+ sources. Breadth-first design (82K compounds).
- **Wave 2** (3 partners): Each partner exchanges with both others. 84.8% of compounds in ALL 3 sources. Per-source replication is 2-4 wells. Combined median 7 reps. Depth-first design (36K compounds).
- This means Wave 2 has dramatically better per-compound evaluation power (30K evaluable vs ~850 for S5/S8) despite having fewer sources and fewer total compounds.
- The `compound_source.csv.gz` from the JUMP datasets repo only shows source-level presence, not well-level replication. Must use `well.csv.gz` joined with `plate.csv.gz` to get actual replicate counts.

### R methods (fastMNN, Seurat CCA/RPCA) cannot scale beyond ~100K cells (2026-04-01)
- fastMNN's `.average_correction` in the `batchelor` package creates a smoothing kernel matrix of `n_target_cells × n_mnn_pairs`. At 244K cells (scenario_3), this exceeds R's 2^31-1 element vector limit. This is a hard R language constraint, not a RAM issue.
- Seurat CCA/RPCA's `IntegrateLayers` creates dense intermediates during anchor finding that grow superlinearly with cell count. All HPO trials get killed (signal -9) after 1-2 hours on 244K cells, even on nodes with 2 TB RAM.
- When all R trials fail (return None), Optuna's TPE sampler crashes on all-None weights (`TypeError`), causing a cascade where trials 11-29 fail instantly without launching R.
- All COMPOUND-plate scenarios (S3, S5, S8, wave2, wave1) exceed 200K cells and will trigger this. TARGET2-only scenarios (S1, S2, S4, S6, S7) stay under 50K and are safe.
- Fix: auto-skip R methods in Snakefile when `n_cells > 100_000`. Threshold is conservative — the actual failure point varies by HPO parameters but 100K gives safe margin.
- Full analysis: `tasks/r_methods_scalability.md`.

### RefChemDB is not viable as DRH replacement for JUMP-CP (2026-04-01)
- RefChemDB (EPA, Judson 2019) provides expert-curated chemical-target annotations with evidence support scores (`count` = number of independent sources). High-quality core is count ≥ 2.
- Chemical space mismatch: RefChemDB covers ~41K chemicals focused on environmental/industrial/tox targets (AR, ESR1, ACHE, DRD2). JUMP-CP TARGET2 compounds are mostly kinase inhibitors and oncology drugs — minimal overlap at high evidence levels.
- Bridge is lossy: InChIKey → PubChem CID → DTXSID loses 40% of compounds (many JUMP tool compounds aren't in EPA's registry).
- At count ≥ 2: only 40/302 TARGET2 cpds map, yielding 2 evaluable targets. At count ≥ 1 (unfiltered): 113 cpds, 77 evaluable — but single-source annotations have 39.5% precision per the RefChemDB paper.
- For compound annotation alternatives, check `target_annotations.parquet` (ChEMBL, 62% coverage of full 116K cpds) before exploring other sources.
- POC code and cached data: `exploration/refchemdb_mapping.py`, `exploration/cache/`.

### TARGET2 dominates evaluation for most scenarios (2026-03-29)
- For most scenarios, COMPOUND plate cross-source overlap is <1%. The 306 TARGET2 compounds (present in ALL sources) dominate mAP evaluation.
- Exception: Wave 2 (84.8% overlap), C8 with bridge source S7 (6.3%), and cross-wave scenarios.
- This means T2-only vs T2+C scenarios have similar evaluation power. The main value of COMPOUND plates is ML breadth, not evaluation strength.
