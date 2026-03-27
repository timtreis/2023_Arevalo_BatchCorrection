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


def objective(trial, input_data, batch_key, label_key, method, cache_path):
    dims = trial.suggest_int("dims", 5, 50)
    k_anchor = trial.suggest_int("k_anchor", 3, 30)
    k_weight = trial.suggest_int("k_weight", 50, 200)

    cmd = [
        "Rscript", "scripts/run_seurat_trial.R",
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

    for line in reversed(result.stdout.strip().split("\n")):
        if line.startswith("RESULT:"):
            batch_score, bio_score = map(float, line[7:].split(","))
            return batch_score, bio_score

    logger.warning("No RESULT line in R output: %s", result.stdout[-500:])
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

    study = optuna.create_study(
        directions=["maximize", "maximize"],
        sampler=optuna.samplers.TPESampler(seed=42),
    )
    study.optimize(
        lambda trial: objective(trial, input_path, batch_key, label_key, method, cache_path),
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
