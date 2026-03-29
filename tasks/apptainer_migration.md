# Apptainer Container Migration

## Motivation

Python on Lustre filesystems is slow. Every job startup triggers thousands of
`stat()` / `open()` calls across the module search path for each `import`.
On Lustre, each of those is a network round-trip against the metadata server.
A single `import scvi` chain can take 30–90 s on a cold cache.

Apptainer squashfs images fix this at the kernel level: the entire Python
installation is a single compressed file that the kernel page-caches on first
access. Subsequent imports are local RAM reads.

Secondary benefit: reproducibility — the exact package set is frozen in the
image, not resolved at job time.

---

## Architecture

```
containers/base.def  →  base.sif        (shared foundation)
containers/scvi.def  →  scvi.sif        (built FROM base.sif)
containers/sysvi.def →  sysvi.sif       (future, built FROM base.sif)
containers/scpoli.def → scpoli.sif      (future, may need different base)
...
```

### base.def

- **Docker base**: `nvidia/cuda:12.1.1-cudnn8-runtime-ubuntu22.04`
- **Python**: 3.11 (installed from Ubuntu packages; pip bootstrapped via
  `get-pip.py` because Ubuntu disables `ensurepip` for system Python)
- **Shared packages**: anndata, pandas, scanpy, pyarrow, scib-metrics, scipy,
  scikit-learn, statsmodels, tqdm
- **NumPy pin**: `numpy<2` — scib-metrics uses `np.in1d` which was removed in
  NumPy 2.0; pinning prevents a silent breakage across all derived images

### Method containers

- `Bootstrap: localimage` / `From: containers/base.sif`
- Add only the method-specific packages (torch, JAX, scvi-tools, etc.)
- Keep build time short by maximising what lives in the shared base

### Multiple Python versions (future)

Some methods need older Python (e.g. scpoli → 3.9). Ubuntu 22.04 only ships
3.11 in the default repos. Options:

- **Versioned bases**: `base-py311.def`, `base-py39.def`
  - Python 3.9 requires Ubuntu 20.04 base or the deadsnakes PPA on 22.04
  - Update `pixi.toml` build tasks accordingly
- Method containers then pick their matching versioned base

---

## Environment variables (base %environment)

All are set in `base.def %environment` so every derived container inherits them.
`%environment` is sourced inside the container shell after host variables, so
these reliably override anything the host passes in.

| Variable | Value | Why |
|---|---|---|
| `NUMBA_CACHE_DIR` | `/tmp/numba_cache` | Numba writes JIT cache next to source files by default; the squashfs image is read-only, so Numba raises `RuntimeError: cannot cache function` at import time |
| `MPLCONFIGDIR` | `/tmp/matplotlib` | The host `MPLCONFIGDIR` often points to a Lustre path that is not mounted inside the container |
| `XLA_PYTHON_CLIENT_PREALLOCATE` | `false` | JAX pre-allocates all free GPU memory on first use; this conflicts with PyTorch running in the same process |
| `TORCH_HOME` | `/tmp/torch_home` | Torch kernel / hub cache defaults to `$HOME/.cache/torch`; `$HOME` inside the container may resolve to an unmounted Lustre path |

### Method-specific overrides

**`scvi.def` `%environment`**:

```
JAX_PLATFORMS=cpu
```

scvi trains with PyTorch on GPU; scib-metrics then runs JAX for the
Optuna benchmark step. PyTorch and JAX cannot both initialise cuDNN in the
same process — JAX fails with `FAILED_PRECONDITION: DNN library initialization
failed`. The latent representations are small (n_obs × n_latent ~10–100 dims),
so CPU is sufficient for the benchmark step.

Do **not** set this in `base.def` — future sysVI and other JAX-based method
containers will need JAX on GPU.

---

## Snakemake 8 changes

| Setting | Old (Snakemake 7) | New (Snakemake 8) |
|---|---|---|
| Software deployment | `--use-conda` | `--sdm conda apptainer` |
| Conda plugin | built-in | `snakemake-software-deployment-plugin-conda` (separate pip package) |
| GPU flag | `--singularity-args '--nv'` | `--apptainer-args=--nv` (= form required; argparse rejects `--flag value` when value starts with `--`) |
| Rule: conda env | `conda: "../envs/scvi.yaml"` | `container: "containers/scvi.sif"` |

Bind mounts: Snakemake automatically binds `$PWD` into the container, so all
relative paths in shell commands work without extra `--bind` flags. If scripts
need to access paths outside the project root (e.g. shared data on a different
Lustre mount point), add `--apptainer-args=--bind /path:/path`.

---

## pixi.toml build tasks

```toml
[tasks]
build-base        = "apptainer build containers/base.sif containers/base.def"
build-scvi        = { cmd = "apptainer build containers/scvi.sif containers/scvi.def", depends-on = ["build-base"] }
build-scibmetrics = { cmd = "apptainer build containers/scibmetrics.sif containers/scibmetrics.def", depends-on = ["build-base"] }
build-containers  = { depends-on = ["build-scvi", "build-scibmetrics"] }
```

Always run `build-containers` (not individual tasks) to ensure the dependency
chain is respected. Add new method containers as additional `build-<method>`
tasks with `depends-on = ["build-base"]` (or the appropriate versioned base).

---

## Conda → container checklist (per method)

- [ ] Identify the method's Python version requirement
  - 3.11: use `base.sif` directly
  - Other: create `base-py<ver>.def` first (see Multiple Python versions above)
- [ ] Create `containers/<method>.def` with `Bootstrap: localimage`
- [ ] Install only method-specific packages in `%post`
- [ ] Add method-specific `%environment` overrides if needed (e.g. `JAX_PLATFORMS`)
- [ ] Add `build-<method>` task to `pixi.toml`
- [ ] In `rules/correct.smk`: replace `conda: "../envs/<method>.yaml"` with `container: "containers/<method>.sif"`
- [ ] In `rules/tune.smk`: same substitution for the corresponding optimize rule
- [ ] Remove the GPU resource claim from the metrics rule for CPU-only methods
- [ ] Build and smoke-test: `pixi run build-containers && pixi run pipeline`
- [ ] Check the job log for warnings (unmounted paths, cache dirs, etc.)

---

## Methods to migrate

| Method | Python | Status | Notes |
|---|---|---|---|
| scvi_single | 3.11 | ✅ Done | JAX_PLATFORMS=cpu in scvi.def |
| scibmetrics (metrics rule) | 3.11 | ✅ Done | JAX_PLATFORMS=cpu in scibmetrics.def; jax[cpu] only |
| scvi_multi | 3.11 | ✅ Done | Same container as scvi_single (scvi.sif) |
| scanvi_single | 3.11 | ✅ Done | Same container (scvi.sif) |
| scanvi_multi | 3.11 | ✅ Done | Same container (scvi.sif) |
| sysvi | 3.11 | ✅ Done | Uses scvi.sif; shell overrides JAX_PLATFORMS=cuda for GPU |
| harmony | 3.11 | ✅ Done | lightweight.sif (harmonypy + scanorama + bbknn + optuna) |
| scanorama | 3.11 | ✅ Done | lightweight.sif |
| bbknn | 3.11 | ✅ Done | lightweight.sif |
| combat | 3.11 | ✅ Done | base.sif (only needs scanpy) |
| sphering | 3.11 | ✅ Done | base.sif |
| desc | 3.11 | ✅ Done | desc.sif (standalone, TF 2.14 + CUDA) |
| scpoli | 3.9 | ✅ Done | scpoli.sif (standalone, Python 3.9 via deadsnakes PPA) |
| fastMNN | R 4.4 | ✅ Done | r.sif (rocker/r-ver:4.4.0 + batchelor + arrow) |
| seurat_cca | R 4.4 | ✅ Done | r.sif |
| seurat_rpca | R 4.4 | ✅ Done | r.sif |
| mnn | 3.11 | ⏭️ Skipped | Ancient deps; replaced by fastMNN |
| gaushvi | 3.11 | ⏭️ Deferred | Custom scvi-tools fork needs investigation |
| gaushanvi | 3.11 | ⏭️ Deferred | Custom scvi-tools fork needs investigation |

---

## Known gotchas

1. **`ensurepip` disabled on Ubuntu**: system Python on Ubuntu 22.04/20.04 has
   `ensurepip` disabled. Bootstrap pip with:
   ```
   curl -sS https://bootstrap.pypa.io/get-pip.py | python3.x
   ```

2. **Torch kernel cache**: even with `TORCH_HOME=/tmp/torch_home`, PyTorch's
   CUDA JIT kernel cache (`torch/kernels`) may use a different path derived
   from `XDG_CACHE_HOME`. If the warning persists, also set
   `XDG_CACHE_HOME=/tmp/xdg-cache`.

3. **Optuna in-process GPU memory**: when running multiple Optuna trials in one
   process, PyTorch does not release GPU memory between trials (it stays in the
   allocator cache). If this causes OOM across trials, consider checkpointing
   or using `torch.cuda.empty_cache()` between trials.

4. **Conda env for metrics**: `metrics_run_scibmetrics_benchmarker` still runs
   in the `scibmetrics` conda env (not a container). JAX in that env must be
   forced to CPU via `export JAX_PLATFORMS=cpu` in the rule's shell command,
   because JAX's GPU initialisation in a conda env can hang indefinitely
   (XLA kernel compilation). The benchmark metrics (silhouette, kBET, etc.)
   do not require GPU.

5. **`np.in1d` removed in NumPy 2.0**: pinned `numpy<2` in `base.def`.
   Re-evaluate when scib-metrics releases a NumPy-2-compatible version.

6. **Module-level network fetches**: `preprocessing/metadata.py` originally
   fetched `MICRO_CONFIG` from GitHub at import time. This fails inside
   containers (SSL may not be configured) and adds latency to every job.
   Fixed by making it a lazy proxy that prefers a local file at
   `inputs/metadata/microscope_config.csv`.
