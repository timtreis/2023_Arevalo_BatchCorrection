import logging
import argparse
import subprocess
import tempfile
import os
import optuna

from optuna_utils import save_optuna_results

logger = logging.getLogger(__name__)


def _build_cache(input_data, batch_key, cache_path):
    """Pre-build Seurat object and save as .rds for fast trial loading."""
    logger.info("Building Seurat cache: %s", cache_path)
    result = subprocess.run(
        [
            "Rscript", "scripts/cache_seurat_object.R",
            "--input_data", input_data,
            "--batch_key", batch_key,
            "--max_dims", "50",
            "--output_rds", cache_path,
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        logger.warning("Cache build failed: %s", result.stderr[-500:] if result.stderr else "no stderr")
        return False
    logger.info("Cache built successfully")
    return True


def _get_min_batch_size(input_data, batch_key):
    """Read the input parquet via R (arrow) and return the smallest batch size."""
    r_code = (
        f'suppressPackageStartupMessages(library(arrow)); '
        f'df <- read_parquet("{input_data}", col_select = "{batch_key}"); '
        f'cat(min(table(df[["{batch_key}"]])))'
    )
    result = subprocess.run(
        ["Rscript", "-e", r_code], capture_output=True, text=True,
    )
    if result.returncode != 0:
        logger.warning("Could not read batch sizes: %s", result.stderr[-300:])
        return None
    return int(result.stdout.strip())


def objective(trial, input_data, batch_key, label_key, method, cache_path, k_weight_max):
    dims = trial.suggest_int("dims", 5, 50)
    k_anchor = trial.suggest_int("k_anchor", 3, 30)
    k_weight = trial.suggest_int("k_weight", 5, k_weight_max)

    cmd = [
        "Rscript", "scripts/run_seurat_trial_v5.R",
        "--batch_key", batch_key,
        "--label_key", label_key,
        "--method", method,
        "--dims", str(dims),
        "--k_anchor", str(k_anchor),
        "--k_weight", str(k_weight),
    ]

    if cache_path and os.path.exists(cache_path):
        cmd += ["--cached_rds", cache_path]
    else:
        cmd += ["--input_data", input_data]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.warning("Trial failed: %s", result.stderr[-500:] if result.stderr else "no stderr")
        raise optuna.TrialPruned(f"R trial failed (exit {result.returncode})")

    stdout = result.stdout.strip()

    # Record effective k_weight if R had to retry with a lower value
    for line in stdout.split("\n"):
        if line.startswith("EFFECTIVE_K_WEIGHT:"):
            trial.set_user_attr("effective_k_weight", int(line.split(":")[1]))
            break

    for line in reversed(stdout.split("\n")):
        if line.startswith("RESULT:"):
            batch_score, bio_score = map(float, line[7:].split(","))
            return batch_score, bio_score

    logger.warning("No RESULT line in R output: %s", stdout[-500:])
    raise optuna.TrialPruned("No RESULT line in R output")


def optimize_seurat(input_path, batch_key, label_key, method, n_trials, output_path, smoketest=False):
    if smoketest:
        n_trials = 2

    # Build cached Seurat object (parquet → .rds) once before all trials
    cache_dir = os.path.dirname(output_path) or "."
    cache_path = os.path.join(cache_dir, f"_seurat_{method}_cache.rds")
    if not _build_cache(input_path, batch_key, cache_path):
        logger.warning("Cache build failed, falling back to per-trial parquet loading")
        cache_path = None

    # Adapt k_weight upper bound to the smallest batch.
    # Seurat's FindWeights requires k_weight <= anchor count per pair.
    # Anchor count scales with batch size and cross-batch similarity.
    # With heterogeneous sources (different microscopes), anchor counts can be
    # much lower than batch size, so we use a conservative fraction.
    min_batch = _get_min_batch_size(input_path, batch_key)
    if min_batch is not None:
        k_weight_max = min(200, max(20, min_batch // 10))
        logger.info("min_batch_size=%d → k_weight range [5, %d]", min_batch, k_weight_max)
    else:
        k_weight_max = 200
        logger.info("Could not determine batch sizes, using k_weight range [5, %d]", k_weight_max)

    study = optuna.create_study(
        directions=["maximize", "maximize"],
        sampler=optuna.samplers.TPESampler(seed=42),
    )

    # Seed with a conservative first trial: high k_anchor (more anchors),
    # low k_weight (safe), and Arevalo default dims.
    study.enqueue_trial({"dims": 30, "k_anchor": 30, "k_weight": min(k_weight_max, 10)})

    study.optimize(
        lambda trial: objective(trial, input_path, batch_key, label_key, method, cache_path, k_weight_max),
        n_trials=n_trials,
        catch=(Exception,),
    )
    save_optuna_results(study, output_path)

    # Clean up cache
    if cache_path and os.path.exists(cache_path):
        os.remove(cache_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Use Optuna to tune hyperparameters for Seurat CCA/RPCA.")
    parser.add_argument("--input_data", required=True)
    parser.add_argument("--batch_key", required=True)
    parser.add_argument("--label_key", required=True)
    parser.add_argument("--method", required=True, choices=["cca", "rpca"])
    parser.add_argument("--n_trials", required=True, type=int)
    parser.add_argument("--output_path", required=True)
    parser.add_argument("--smoketest", action="store_true")
    args = parser.parse_args()

    optimize_seurat(
        input_path=args.input_data,
        batch_key=args.batch_key,
        label_key=args.label_key,
        method=args.method,
        n_trials=args.n_trials,
        output_path=args.output_path,
        smoketest=args.smoketest,
    )
