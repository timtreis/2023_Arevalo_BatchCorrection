#!/usr/bin/env bash
# Overnight pipeline runner — chains all pending work sequentially.
# Created 2026-04-12. Run from repo root with: bash scripts/overnight_pipeline.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

LOG_DIR="outputs/overnight_20260412"
mkdir -p "$LOG_DIR"
MASTER_LOG="$LOG_DIR/master.log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$MASTER_LOG"; }

# Use pixi run to get correct LD_LIBRARY_PATH (system libstdc++ too old for pixi env's ICU libs)
# Wrap snakemake in pixi run -- bash -c '...' to inherit the pixi environment
PIXI="$REPO_ROOT/.pixi/bin/pixi"
if [ ! -x "$PIXI" ]; then
    PIXI="$(command -v pixi)"
fi

# Set up PATH so snakemake can find python
PYTHON_DIR="$(dirname .pixi/envs/default/bin/python3.13)"
export PATH="$PYTHON_DIR:$PATH"

# Use pixi run for snakemake to get correct library paths
SMK_PREFIX="pixi run --"

smk_run() {
    local name="$1"; shift
    log "START: $name"
    if $SMK_PREFIX bash -c "PATH=$PYTHON_DIR:\$PATH snakemake $* 2>&1" | tee -a "$LOG_DIR/${name}.log"; then
        log "DONE:  $name (success)"
    else
        log "FAIL:  $name (exit $?) — continuing with next task"
    fi
}

smk_touch() {
    local name="$1"; shift
    log "TOUCH: $name — accepting existing outputs"
    $SMK_PREFIX bash -c "PATH=$PYTHON_DIR:\$PATH snakemake --touch --forceall $* 2>&1" | tee -a "$LOG_DIR/${name}_touch.log" || true
}

SHARED="inputs/conf/shared.json"
SMK_COMMON=(-c20 --sdm apptainer --apptainer-args=--nv --resources nvidia_gpu=1 --keep-going)

log "=========================================="
log "OVERNIGHT PIPELINE START"
log "Node: $(hostname), GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'unknown')"
log "=========================================="

# ============================================================
# STEP 0: Wait for currently running pipelines
# ============================================================
log "Waiting for running Snakemake processes to finish..."
while pgrep -u "$(whoami)" -f "snakemake.*scenario_3.json" > /dev/null 2>&1; do
    sleep 30
done
log "S3 pipeline finished."

while pgrep -u "$(whoami)" -f "snakemake.*scenario_1_defaults" > /dev/null 2>&1; do
    sleep 30
done
log "S1 defaults pipeline finished."

# Clear any stale locks
rm -f .snakemake/locks/0.input.lock .snakemake/locks/0.output.lock 2>/dev/null
mkdir -p .snakemake/locks

# ============================================================
# STEP 1: Verify S3 + S1 defaults completed
# ============================================================
log "Verifying S3 completion..."
for f in outputs/scenario_3/metrics/mad_int_featselect_scibmetrics_benchmarker.parquet \
         outputs/scenario_3/plots/results_table.pdf; do
    if [ -f "$f" ]; then
        log "  OK: $(basename $f)"
    else
        log "  MISSING: $f"
    fi
done

log "Verifying S1 defaults completion..."
if [ -f "outputs/scenario_1_defaults/plots/results_table.pdf" ]; then
    log "  OK: S1 defaults complete"
else
    log "  NOTE: S1 defaults may need attention"
fi

# ============================================================
# STEP 2: S2 defaults
# ============================================================
smk_run "s2_defaults" "$SMK" "${SMK_COMMON[@]}" --configfile "$SHARED" inputs/conf/scenario_2_defaults.json

# ============================================================
# STEP 3: S3 defaults
# ============================================================
smk_run "s3_defaults" "$SMK" "${SMK_COMMON[@]}" --configfile "$SHARED" inputs/conf/scenario_3_defaults.json

# ============================================================
# STEP 4: S1 HPO — accept existing outputs, then run missing seurat v4/v5
# ============================================================
log "S1 HPO: touching existing outputs to accept code hash changes..."
smk_touch "s1_hpo" -c20 --configfile "$SHARED" inputs/conf/scenario_1.json --sdm apptainer --apptainer-args=--nv --resources nvidia_gpu=1

log "S1 HPO: running pipeline (should only run missing seurat v4/v5 methods)..."
smk_run "s1_hpo" "$SMK" "${SMK_COMMON[@]}" --configfile "$SHARED" inputs/conf/scenario_1.json

# ============================================================
# STEP 5: S2 HPO — same pattern
# ============================================================
smk_touch "s2_hpo" -c20 --configfile "$SHARED" inputs/conf/scenario_2.json --sdm apptainer --apptainer-args=--nv --resources nvidia_gpu=1

smk_run "s2_hpo" "$SMK" "${SMK_COMMON[@]}" --configfile "$SHARED" inputs/conf/scenario_2.json

# ============================================================
# STEP 6: S4 HPO — already has _v4 corrections, needs _v5 + re-aggregation
# ============================================================
smk_touch "s4_hpo" -c20 --configfile "$SHARED" inputs/conf/scenario_4.json --sdm apptainer --apptainer-args=--nv --resources nvidia_gpu=1

smk_run "s4_hpo" "$SMK" "${SMK_COMMON[@]}" --configfile "$SHARED" inputs/conf/scenario_4.json

# ============================================================
# DONE
# ============================================================
log "=========================================="
log "OVERNIGHT PIPELINE COMPLETE"
log "=========================================="

log ""
log "Summary of output files:"
for scenario in scenario_1 scenario_2 scenario_4 scenario_1_defaults scenario_2_defaults scenario_3_defaults; do
    if [ -f "outputs/$scenario/plots/results_table.pdf" ]; then
        log "  $scenario: plots OK ($(stat -c '%y' outputs/$scenario/plots/results_table.pdf 2>/dev/null | cut -d. -f1))"
    else
        log "  $scenario: plots MISSING"
    fi
done
