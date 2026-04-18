"""Pinned paths for the reference_mapping experiments.

All paths are absolute. The notebooks import these symbols rather than
constructing strings inline so that a path change is a single-file edit.
"""
from pathlib import Path

# Pipeline repo root (parent of reference_mapping/)
REPO_ROOT = Path(__file__).resolve().parents[2]

# Pipeline outputs (atlas-trainable parquets + saved models live here)
PIPELINE_OUT = REPO_ROOT / "outputs"

# reference_mapping/ artifact roots
REFMAP_ROOT = REPO_ROOT / "reference_mapping"
DATA_OUT = REFMAP_ROOT / "data"
MODEL_OUT = REFMAP_ROOT / "models"
RESULTS_OUT = REFMAP_ROOT / "results"
PAPERMILL_OUT = REFMAP_ROOT / "papermill_runs"

# Preproc prefix used by the Snakemake rules; matches config["preproc"].
PREPROC = "mad_int_featselect"


def scenario_dir(scenario: str) -> Path:
    """`outputs/<scenario>/`"""
    return PIPELINE_OUT / scenario


def scenario_input_parquet(scenario: str) -> Path:
    """Pre-correction features for a scenario (input to all methods)."""
    return scenario_dir(scenario) / f"{PREPROC}.parquet"


def scenario_method_parquet(scenario: str, method: str) -> Path:
    """Corrected-latent parquet for a method, e.g. method='scvi_single'."""
    return scenario_dir(scenario) / f"{PREPROC}_{method}.parquet"


def scenario_model_dir(scenario: str, method: str) -> Path:
    """Saved-model dir produced by the patched correct_with_*.py rules.

    method is one of: 'scvi_single', 'scpoli'.
    """
    return scenario_dir(scenario) / f"{PREPROC}_{method}_model"


def scenario_optuna_csv(scenario: str, method: str) -> Path:
    """Optuna best-params CSV consumed by the correct_with_* scripts."""
    return scenario_dir(scenario) / "optimization" / f"optuna_{method}.csv"
