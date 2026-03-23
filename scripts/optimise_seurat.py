import logging
import argparse
import subprocess
import optuna

from optuna_utils import save_optuna_results

logger = logging.getLogger(__name__)


def objective(trial, input_data, batch_key, label_key, method):
    dims = trial.suggest_int("dims", 5, 50)
    k_anchor = trial.suggest_int("k_anchor", 3, 30)
    k_weight = trial.suggest_int("k_weight", 50, 200)

    result = subprocess.run(
        [
            "Rscript", "scripts/run_seurat_trial.R",
            "--input_data", input_data,
            "--batch_key", batch_key,
            "--label_key", label_key,
            "--method", method,
            "--dims", str(dims),
            "--k_anchor", str(k_anchor),
            "--k_weight", str(k_weight),
        ],
        capture_output=True,
        text=True,
    )

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

    study = optuna.create_study(
        directions=["maximize", "maximize"],
        sampler=optuna.samplers.TPESampler(seed=42),
    )
    study.optimize(
        lambda trial: objective(trial, input_path, batch_key, label_key, method),
        n_trials=n_trials,
    )
    save_optuna_results(study, output_path)


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
