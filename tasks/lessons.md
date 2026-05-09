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

### kmeans NMI/ARI is degenerate with high-cardinality labels (2026-04-04)
- scib-metrics sets `n_clusters = n_labels`. With 82K compound labels on 244K cells, that's ~3 cells/cluster — KMeans can't converge meaningfully.
- FAISS warns "please provide at least 3.2M training points" for 82K clusters.
- Results are dominated by initialization noise. Sphering's anomalous 0.45 (vs ~0 for all other methods) was an artifact, not genuine bio conservation.
- Fix: skip `nmi_ari_cluster_labels_kmeans` in `BioConservation()` when `n_labels > 500`. Threshold cleanly separates TARGET2-only (302 labels) from COMPOUND-plate scenarios (thousands+).
- kmeans NMI has the same issue (same n_clusters logic). Both are disabled together.
- `isolated_labels` is NOT degenerate — it works correctly even with high-cardinality labels. In S1/S2 only 0.3% of labels are isolated. In S3, 18.3% are isolated but the metric still has signal (just less informative).

### SLURM time limits kill long-running interactive jobs (2026-04-04)
- Interactive `srun` sessions have time limits (varies by partition). S4 pipeline was killed at 22:57 after ~1.5 hours.
- Snakemake is resilient to this — completed HPO CSVs and correction parquets survive. Just clear locks and relaunch.
- For long pipelines, consider using `sbatch` with explicit `--time` instead of interactive sessions, or request longer time limits.

### FAISS GPU is compiled for specific CUDA compute capabilities (2026-04-05)
- FAISS GPU binaries in pip/conda are compiled for specific architectures (sm_70, sm_80). H100 needs sm_90 which isn't included.
- The failure is a C-level assertion abort (`CUDA error 209: no kernel image available`), not a Python exception — try/except cannot catch it.
- Fix: detect GPU architecture before launching Python (via `nvidia-smi --query-gpu=name`), set `CUDA_VISIBLE_DEVICES=''` to force CPU mode.
- CPU FAISS is fine for <50K cells. For larger datasets, rebuild FAISS from source with sm_90 support or use the scibmetrics-gpu container which has JAX GPU (unaffected).

### Seurat CCA fails on multi-source cross-microscope data (2026-04-05)
- S4 (5 sources, 3 microscope types): all 30 HPO trials produced non-finite values (700K-845K NaN/Inf entries) for every k.weight attempted.
- After ~10 failed trials returning None, Optuna's TPE sampler crashes on all-None weights (`TypeError: NoneType * float`), causing trials 11-29 to fail instantly without launching R.
- RPCA handles the same data fine (uses PCA before integration, more robust to distribution differences).
- Fix: auto-skip methods with 0 COMPLETE HPO trials in Snakefile (reads optuna CSV at DAG time).

### gaushvi/gaushanvi are superseded by mainline scvi-tools (2026-04-05)
- These methods required a custom fork (`git+https://github.com/mBiocoder/scvi-tools.git@scvi_add_normal`) to add `gene_likelihood="normal"`.
- The `gene_likelihood="normal"` option is now in mainline scvi-tools (already in scvi.sif container).
- New `scvi_normal` method replaces both: uses standard scvi-tools, has full HPO (old gaushvi had hardcoded params), skips unnecessary min-value shift, omits labels_key.
- No need to maintain separate containers or conda envs for these methods.

### scVI ZINB likelihood is theoretically wrong for Cell Painting (2026-04-05)
- CellProfiler features are continuous, can be negative, and are not counts. ZINB/NB likelihoods assume non-negative integer data.
- The pipeline's `adata.X -= min_value` shift is a workaround to make ZINB not crash, but it distorts the data.
- `gene_likelihood="normal"` (Gaussian) is the correct likelihood for continuous morphological features. scPoli already uses MSE (equivalent).
- Added `scvi_normal` method to test this hypothesis head-to-head against `scvi_single` (ZINB).

### scPoli early stopping is effectively broken with n_train_epochs=999999 (2026-04-06)
- scPoli's `EarlyStopping` has `patience=20` but also `reduce_lr=True` with `lr_patience=13` and `lr_factor=0.1`.
- When loss stops improving, LR is reduced after 13 epochs. Even tiny fluctuations from reduced LR can reset the patience counter. Then another 13 epochs → another LR reduction → repeat.
- **Observed on S5 (370K cells)**: Loss minimum at epoch ~12 (val_loss=611). Model trained to epoch 705+ before being killed — it would have continued to 999999 epochs. Only 2 LR reductions occurred (LR=0.00001), but micro-fluctuations kept resetting patience indefinitely.
- **Fix applied**: Capped `n_train_epochs` from 999999 → 400 in `correct_with_scpoli.py`. With `reload_best=True`, the saved model uses epoch ~12's weights regardless, so output quality is unaffected. 400 epochs is generous — best model is always found within first 50.
- Alternative mitigation: `early_stopping_kwargs={"reduce_lr": False, "patience": 20}` would fix the root cause but requires passing through the model.train() call.

### Modifying shared scripts triggers Snakemake reruns for existing outputs (2026-04-06)
- Adding `--gene_likelihood` arg to `optimise_scvi.py` caused Snakemake to rerun scvi_single HPO on S5 (script listed as input, content changed).
- S5 had to redo 30 scvi_single trials (~2 hrs) even though the output would be identical (default is still `zinb`).
- This is by design (Snakemake tracks input file hashes), but means script changes should be made BETWEEN scenario runs, not during.
- Alternative: use `--rerun-triggers mtime` to only rerun on mtime changes, not content changes.

### Seurat v4 vs v5 API incompatibility breaks default parameters (2026-04-07)
- Arevalo et al. used Seurat v4.4 (`envs/seurat.yaml` pins `r-seurat=4.4`). Our containers had Seurat 5.4.0.
- v4 uses `FindIntegrationAnchors()` + `IntegrateData()`. v5 uses `IntegrateLayers()` with different anchor-finding internals.
- Arevalo's default params (`k_anchor=5, k_weight=100`) fail on v5 for S1 — too few anchors for the requested k_weight. Our HPO-tuned params (`k_anchor=19, k_weight=51`) work fine on v5.
- Fix: split into `r_v4.sif` (Seurat 4.4.0 for Arevalo reproduction) and `r_v5.sif` (Seurat 5.4.0 for HPO). fastMNN also uses `r_v4.sif` since batchelor hasn't changed.
- When building the v4 container with `remotes::install_version("Seurat", "4.4.0")`, SeuratObject 5.x gets pulled as a dep. The v4 API still works — backwards compatible.

### Snakemake `--touch --forceall` resets code hash metadata (2026-04-07)
- `touch` on files alone doesn't fix Snakemake hash mismatches — it tracks code hashes in `.snakemake/metadata/`, not just mtime.
- `snakemake --touch --forceall` marks ALL outputs as up-to-date with current script hashes without running anything. This is the correct way to "accept" existing outputs after non-behavioral script changes.
- Use this when script changes don't affect output (e.g., adding logging, unused args) but Snakemake wants to rerun everything.

### scvi_normal (Gaussian likelihood scVI) is not viable for Cell Painting (2026-04-07)
- scVI's `gene_likelihood="normal"` uses `exp()` for variance parameterization, which produces NaN when encoder outputs are large.
- All 30 HPO trials on S5 crashed with `ValueError: Expected parameter loc ... to satisfy Real(), but found invalid values: tensor([[nan, ...]])`.
- sysVI already provides a proper Gaussian VAE with softplus + clamping + nan_to_num guards. No need for scvi_normal.
- scvi_normal removed from METHODS list (plumbing kept in scripts for potential future use).

### Always use `pixi run` for snakemake, never call the binary directly (2026-04-12)
- Calling `.pixi/envs/default/bin/snakemake` directly skips pixi's `LD_LIBRARY_PATH` setup.
- System `/lib64/libstdc++.so.6` is too old for pixi env's ICU libs (`CXXABI_1.3.15` not found).
- This manifests as an `ImportError` in `pycytominer → sqlite3 → dbapi2` during Snakefile parsing.
- Fix: always use `pixi run scenario-X` or `pixi run -- bash -c 'PATH=... snakemake ...'`.
- The overnight_pipeline.sh failure (S2/S3 defaults + S1/S2/S4 HPO all failed) was caused by this. resume_pipeline.sh uses `pixi run` and works.

### Scibmetrics on 244K cells takes ~8 hours (2026-04-12)
- `isolated_labels` metric alone is ~55 min per embedding. With 9 embeddings × 10 metrics, total run is ~8 hours on 244K cells.
- Plan SLURM allocations accordingly: scibmetrics on S3-size scenarios needs the full run + buffer. A 24h job with ~15h of method corrections already consumed will NOT fit scibmetrics.
- Mitigation: run corrections and scibmetrics in separate SLURM jobs, or request longer allocations (>24h partitions) upfront.

### Defaults-mode methods can be dramatically slower than HPO'd (2026-04-12)
- harmony_v1 with default `max_iter_harmony=999999` on 244K cells: ~18 min/iteration, needs ~20 iterations = ~6 hours. HPO'd harmony uses 10-50 iterations max.
- scpoli with default params on 244K cells: early stopping never fires (LR reduction resets patience), runs full 400 epoch cap = ~6 hours. HPO'd scpoli converges in ~50 epochs.
- This is actually useful data for the defaults-vs-HPO comparison paper figure — quantifies the compute cost of not tuning.

### source_N scenario number ≠ source in scenario_N (2026-04-19)
- `source_5` is a Wave 2 source (sources 5, 9, 11). It does NOT appear in `scenario_5` (which has sources 2, 3, 6, 8, 10 — a C1 benchmark scenario).
- Scenario numbers and source numbers are independent. Always verify actual source membership by reading the parquet or config JSON before writing QUERY_SOURCE/REFERENCE_SCENARIO params.
- Quick check: `pd.read_parquet(pq, columns=["Metadata_Source"])["Metadata_Source"].value_counts()`.

### sysVI is the correct VAE method for CellProfiler features, not scVI (2026-04-19)
- scVI uses ZINB/NB likelihood, designed for integer count data (UMI counts). CellProfiler features are continuous, can be negative, and are not counts — ZINB is theoretically wrong.
- sysVI uses a Gaussian VAE with proper variance parameterization (softplus + clamping), which is correct for continuous morphological features. scPoli uses MSE (equivalent).
- scvi_normal (Gaussian scVI) was tested and crashed — NaN in encoder due to `exp()` variance parameterization. sysVI avoids this with better numerical guards.
- sysVI has the same scArches surgery API: `SysVI.prepare_query_anndata()` + `SysVI.load_query_data()`, so reference mapping works identically to scVI.
- No X shift needed for sysVI (unlike scVI where `X -= X.min()` is needed for ZINB).

### scPoli prototype count breaks initialization beyond ~5K unique labels (2026-04-19)
- scPoli uses a prototype embedding for every unique label (compound). In S5 with min_batches=5, that's 682 prototypes — works fine, starts GPU training within minutes.
- Lowering min_batches to 3 or 2 adds 67K or 81K prototypes. The process hangs at 98% CPU / 0% GPU for 25+ minutes and never starts training. GPU memory stays at baseline (~447 MiB).
- Root cause is likely scPoli's DataLoader/label-encoder setup being O(n_labels²) or similar. A dense one-hot tensor of shape (370K × 81K) = 30B entries would explain both the memory and time.
- The S5 compound distribution is bimodal: 682 in all 5 sources, then 33K+ in 4 sources, 33K+ in 3 sources — no middle ground. Any min_batches < 5 jumps to 34K+ prototypes.
- scPoli was designed for scRNA-seq with 10–50 cell types. It degrades non-linearly past ~1K unique labels.

### scPoli S5 failure = supervision sparsity × batch effect strength (2026-04-19)
- S5 has 3 different microscope types (CV8000 confocal, Opera Phenix confocal, ImageXpress widefield) — strongest batch effects in the benchmark.
- With min_batches=5, only 27% of cells are labeled (682 compounds in all 5 sources). The 73% unlabeled COMPOUND cells get no prototype loss and cannot be anchored across batches by the VAE alone.
- S4 (same microscopes, T2-only = 99.8% labeled): scPoli ranks 3/15. S3 (T2+COMPOUND, 56% labeled, all-CV8000): scPoli ranks 1/9. Both dimensions must be present to cause failure.
- Harmony v1 (rank 1 in S5, 0.383 overall) is fully unsupervised — iteratively minimizes batch-driven PCA variance for every cell, directly correcting the microscope-type signal without needing labels.
- The labeled fraction and prototype count are in tension: more labeled cells requires fewer prototype anchors to fail (but too many anchors breaks scPoli's initialization).

### SysVI cannot do reference mapping with new (unseen) batch categories (2026-04-19)
- `SysVI.load_query_data` hard-codes `transfer_batch=False`, which triggers a check requiring all query batches to exist in the reference registry. New batches raise `ValueError`.
- This is architectural: SysVI learns system-specific (batch-specific) transformations at train time; there is no surgery path for new batches.
- scVI (SCVI) DOES support new batches via surgery (`transfer_batch=True` default). scPoli is purpose-built for reference mapping with new batches.
- For refmap experiments: use scPoli or Symphony for true new-batch projection. SysVI is only appropriate when all batches are seen at train time (standard batch correction).

### scarches 0.6.1 incompatible with anndata >= 0.10 (2026-04-19)
- `anndata.read` was removed in anndata 0.10. scarches 0.6.1 still does `from anndata import AnnData, read` in 3 files: `models/base/_base.py`, `models/trvae/trvae_model.py`, `models/expimap/expimap_model.py`.
- Fix: patch each file with `try: from anndata import read \n except ImportError: from anndata.io import read_h5ad as read`.
- This is a known upstream issue; the fix must be applied to the installed package until scarches releases a fix.

### harmonypy 0.2.0 changed Z_corr and R to pre-transposed convention (2026-04-19)
- Old harmonypy: `Z_corr` was `(d, N)`, `R` was `(K, N)`. symphonypy 0.2.2 was written for this.
- harmonypy 0.2.0: `Z_corr` property returns `(N, d)` and `R` returns `(N, K)` (already transposed). symphonypy's extra `.T` gives `(d, N)` and `(K, N)` in obsm — wrong shapes.
- Fixed in `.pixi/envs/symphony/lib/.../symphonypy/_utils.py`: removed `.T` from `Z_corr`, added `.T` to `R` in matrix products, changed `sum(axis=1)→sum(axis=0)` for `Nr`, stored `R.T` in uns.
- Check harmonypy version before using symphonypy: 0.1.x uses old convention, 0.2.x uses new convention.

### symphonypy map_embedding requires sc.pp.scale to be run on reference before saving (2026-04-19)
- `sp.tl.map_embedding` reads `ref.var["mean"]` and `ref.var["std"]` to z-score the query before PCA projection. These are written by `sc.pp.scale(ref, zero_center=True)`.
- Must run scale BEFORE pca BEFORE harmony — the order matters. Merge all three into one notebook cell to prevent out-of-order execution in a live kernel.
- If atlas h5ad was saved without scale stats, nb31 can self-heal: load the raw reference from query_arms/, run scale, inject mean/std into `ref.var` before calling `map_embedding`.
- symphonypy `map_embedding` also defaults `use_genes_column="highly_variable"` — pass `use_genes_column=None` for non-genomic data (CellProfiler features have no HVG selection).

### TARGET2 dominates evaluation for most scenarios (2026-03-29)
- For most scenarios, COMPOUND plate cross-source overlap is <1%. The 306 TARGET2 compounds (present in ALL sources) dominate mAP evaluation.
- Exception: Wave 2 (84.8% overlap), C8 with bridge source S7 (6.3%), and cross-wave scenarios.
- This means T2-only vs T2+C scenarios have similar evaluation power. The main value of COMPOUND plates is ML breadth, not evaluation strength.

### Raw metadata is sufficient to compute compound overlap before preprocessing (2026-04-19)
- `inputs/metadata/well.csv.gz` (Source, Plate, Well, JCP2022) joined with `inputs/metadata/plate.csv.gz` (Source, Plate, PlateType) gives compound membership per source per plate type
- Use this to check overlap for any source combination without running the preprocessing pipeline
- Pattern: `meta = wells.merge(plates[["Metadata_Source","Metadata_PlateType","Metadata_Plate"]], on=["Metadata_Source","Metadata_Plate"])`

### Wave2 and scenario_5 share ~2K non-TARGET2 COMPOUND plate compounds per source (2026-04-19)
- From raw metadata: source_5 shares 2,409 compounds with S5 (7.4%), source_9 shares 1,991 (6.2%), source_11 shares 2,410 (7.4%)
- ~292 of those are TARGET2; the remaining ~1,700-2,100 are COMPOUND plate compounds
- This partial overlap is what makes wave2 the correct choice for a reference mapping experiment with genuine compound coverage gaps
- No preprocessed source (apricot extras: TARGET2-only; source_1: 100% overlap, no T2) provides meaningful partial overlap — wave2 preprocessing is required

### Within scenario_5, LOSO gives ~100% compound overlap — this is expected (2026-04-19)
- All 5 scenario_5 sources profiled the same JUMP-CP compound library; LOSO tests pure batch correction, not compound generalisation
- source_10/source_2/source_8: 100% overlap; source_6: 98.2% (1,040 source_6-only); source_3: 99.6%
- All 5 sources have both COMPOUND and TARGET2 plate types → T+/T- ablation is valid for any held-out source

### JUMP-CP data is well-level aggregated profiles, not single-cell (2026-04-19)
- Each row in the preprocessed parquets = one plate well (aggregated CellProfiler features across all cells in that well)
- Confirmed: `observations_per_well == 1.0` for scenario_5 reference
- Consequence: "cells per compound" = wells per compound, typically 1-4 for COMPOUND plates, 54-300+ for TARGET2 compounds across all sources
- This means compound ASW computed on query-alone is usually undefined (98% of COMPOUND plate compounds have 1 well in a single source)
- Correct evaluation metric: cross-source compound precision@k in the joint reference+query embedding

### TARGET2 compounds have vastly more wells than COMPOUND plate compounds (2026-04-19)
- In scenario_5 reference: TARGET2 compounds have median 64 wells each (across 4 sources), DMSO has 37,852 wells
- COMPOUND plate compounds: median 3 wells total, 1 well in a single source is common
- This makes TARGET2 the high-quality prototype anchor tier and COMPOUND plates the noisy tier in scPoli
- The T+/T- ablation tests whether the high-quality TARGET2 anchors are needed — removing them leaves only 1-well COMPOUND plate prototypes

### scPoli with compound identity as cell_type_keys: valid but at unusual scale (2026-04-19)
- Using Metadata_JCP2022 as cell_type_keys creates ~57K prototypes (one per compound) — far beyond the 50-200 designed for scRNA-seq
- Prototype distance computation is O(batch_size × n_prototypes) per forward pass — feasible (256 × 57K × 4B ≈ 58MB per batch)
- Key limitation: prototypes are noisy for 1-well compounds; the loss landscape is flat for those
- Reference prototypes are frozen during query fine-tuning; only the new batch/condition embedding trains
- Passing labeled_indices = shared compound wells focuses prototype loss where signal exists; wave2-only compounds should be unlabeled
