import logging
import argparse
import subprocess
import os
import optuna

from optuna_utils import save_optuna_results

logger = logging.getLogger(__name__)


def _build_cache(input_data, batch_key, label_key, cache_path):
    """Pre-load parquet data and save as .rds for fast trial loading."""
    logger.info("Building fastMNN cache: %s", cache_path)
    result = subprocess.run(
        [
            "Rscript", "scripts/cache_fastmnn_data.R",
            "--input_data", input_data,
            "--batch_key", batch_key,
            "--label_key", label_key,
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


def objective(trial, input_data, batch_key, label_key, cache_path):
    k = trial.suggest_int("k", 5, 50)
    d = trial.suggest_int("d", 5, 50)
    ndist = trial.suggest_int("ndist", 1, 5)
    prop_k = trial.suggest_float("prop_k", 0.01, 0.5)

    cmd = [
        "Rscript", "scripts/run_fastmnn_trial.R",
        "--batch_key", batch_key,
        "--label_key", label_key,
        "--k", str(k),
        "--d", str(d),
        "--ndist", str(ndist),
        "--prop_k", str(prop_k),
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


def optimize_fastmnn(input_path, batch_key, label_key, n_trials, output_path, smoketest=False):
    if smoketest:
        n_trials = 2

    # Build cached data (parquet → .rds) once before all trials
    cache_dir = os.path.dirname(output_path) or "."
    cache_path = os.path.join(cache_dir, "_fastmnn_cache.rds")
    if not _build_cache(input_path, batch_key, label_key, cache_path):
        logger.warning("Cache build failed, falling back to per-trial parquet loading")
        cache_path = None

    study = optuna.create_study(
        directions=["maximize", "maximize"],
        sampler=optuna.samplers.TPESampler(seed=42),
    )
    study.optimize(
        lambda trial: objective(trial, input_path, batch_key, label_key, cache_path),
        n_trials=n_trials,
        catch=(Exception,),
    )
    save_optuna_results(study, output_path)

    # Clean up cache
    if cache_path and os.path.exists(cache_path):
        os.remove(cache_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Use Optuna to tune hyperparameters for fastMNN.")
    parser.add_argument("--input_data", required=True)
    parser.add_argument("--batch_key", required=True)
    parser.add_argument("--label_key", required=True)
    parser.add_argument("--n_trials", required=True, type=int)
    parser.add_argument("--output_path", required=True)
    parser.add_argument("--smoketest", action="store_true")
    args = parser.parse_args()

    optimize_fastmnn(
        input_path=args.input_data,
        batch_key=args.batch_key,
        label_key=args.label_key,
        n_trials=args.n_trials,
        output_path=args.output_path,
        smoketest=args.smoketest,
    )
