# reference_mapping/

Ad-hoc reference-mapping experiments for the JUMP-CP integration paper. Lives
inside the `2023_Arevalo_BatchCorrection` repo but is intentionally separate
from the Snakemake pipeline. Notebooks here load atlas models produced by the
patched `correct_with_scvi.py` / `correct_with_scpoli.py` rules.

See `projects/jump-cp-integration/research/reference-mapping-implementation-plan.md`
in the coach vault for the full design.

## Quickstart

```bash
# 1. allocate GPU on cluster
srun --gres=gpu:1 --time=8:00:00 --mem=64G --cpus-per-task=8 --partition=gpu_p --pty bash

# 2. install envs (first time only)
cd reference_mapping
pixi install
pixi run install-kernels

# 3. launch jupyter on compute node
pixi run -e symphony jupyter lab --no-browser --ip=0.0.0.0 --port=8888

# 4. tunnel from laptop
ssh -L 8888:<compute-node>:8888 cluster-login
```

## Notebook order (acute first-pass for source_5)

1. `00_derive_target2_manifest`
2. `01_prepare_query_arms`
3. `30_train_atlas_harmony_symphony` → `31_map_symphony` (cheapest smoke test)
4. `40_metrics` (Symphony) → checkpoint 1
5. `10_train_atlas_scvi` → `11_map_scvi` → `40_metrics` (scVI)
6. `20_train_atlas_scpoli` → `21_map_scpoli` → `40_metrics` (scPoli)
7. `41_plot_acute`

## Layout

| Path | Purpose | Tracked? |
|---|---|---|
| `notebooks/` | jupytext-paired `.ipynb` + `.py` | yes |
| `src/` | importable helpers | yes |
| `results/acute/*.csv,*.pdf` | metrics + plots | yes |
| `data/` | query arm + mapped h5ads | no (gitignored) |
| `models/` | atlas snapshots / symlinks | no (gitignored) |
| `papermill_runs/` | papermill execution logs | no (gitignored) |

## Kernels

| Notebook | Kernel |
|---|---|
| 00, 01, 30, 31, 40, 41 | `refmap-symphony` |
| 10, 11 | `refmap-scvi` |
| 20, 21 | `refmap-scpoli` |
