#!/usr/bin/env bash
# Resume remaining pipeline work after overnight partial failure.
# Root cause: direct snakemake call missed pixi's LD_LIBRARY_PATH.
# Fix: use pixi run <task-name> which sets up the full environment.
#
# Created 2026-04-12. Run from repo root with: bash scripts/resume_pipeline.sh
set -uo pipefail  # no -e: we handle errors per step

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

LOG_DIR="outputs/overnight_20260412"
mkdir -p "$LOG_DIR"
MASTER_LOG="$LOG_DIR/master_resume.log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$MASTER_LOG"; }

run_step() {
    local name="$1"
    local cmd="$2"
    log "START: $name"
    if eval "$cmd" 2>&1 | tee -a "$LOG_DIR/${name}.log"; then
        log "DONE:  $name (success)"
        return 0
    else
        local rc=$?
        log "FAIL:  $name (exit $rc) — continuing"
        return $rc
    fi
}

# Clear stale locks
rm -f .snakemake/locks/0.input.lock .snakemake/locks/0.output.lock 2>/dev/null
mkdir -p .snakemake/locks

log "=========================================="
log "RESUME PIPELINE START"
log "Node: $(hostname), GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'unknown')"
log "=========================================="

# ============================================================
# STEP 1: S2 defaults
# ============================================================
run_step "s2_defaults" "pixi run scenario-2-defaults"

# ============================================================
# STEP 2: S3 defaults
# ============================================================
run_step "s3_defaults" "pixi run scenario-3-defaults"

# ============================================================
# STEP 3: S1 HPO — touch existing, then run for missing seurat v4/v5
# ============================================================
log "S1 HPO: touching existing outputs..."
pixi run -- bash -c 'PATH=$(dirname $(which python)):$PATH snakemake --touch --forceall -c20 --configfile inputs/conf/shared.json inputs/conf/scenario_1.json --sdm apptainer --apptainer-args=--nv --resources nvidia_gpu=1' 2>&1 | tee -a "$LOG_DIR/s1_hpo_touch.log" || true

run_step "s1_hpo" "pixi run scenario-1"

# ============================================================
# STEP 4: S2 HPO — same pattern
# ============================================================
log "S2 HPO: touching existing outputs..."
pixi run -- bash -c 'PATH=$(dirname $(which python)):$PATH snakemake --touch --forceall -c20 --configfile inputs/conf/shared.json inputs/conf/scenario_2.json --sdm apptainer --apptainer-args=--nv --resources nvidia_gpu=1' 2>&1 | tee -a "$LOG_DIR/s2_hpo_touch.log" || true

run_step "s2_hpo" "pixi run scenario-2"

# ============================================================
# STEP 5: S4 HPO — same pattern
# ============================================================
log "S4 HPO: touching existing outputs..."
pixi run -- bash -c 'PATH=$(dirname $(which python)):$PATH snakemake --touch --forceall -c20 --configfile inputs/conf/shared.json inputs/conf/scenario_4.json --sdm apptainer --apptainer-args=--nv --resources nvidia_gpu=1' 2>&1 | tee -a "$LOG_DIR/s4_hpo_touch.log" || true

run_step "s4_hpo" "pixi run scenario-4"

# ============================================================
# DONE
# ============================================================
log "=========================================="
log "RESUME PIPELINE COMPLETE"
log "=========================================="

log ""
log "Summary:"
for scenario in scenario_1 scenario_2 scenario_4 scenario_1_defaults scenario_2_defaults scenario_3_defaults; do
    if [ -f "outputs/$scenario/plots/results_table.pdf" ]; then
        log "  $scenario: plots OK ($(stat -c '%y' "outputs/$scenario/plots/results_table.pdf" 2>/dev/null | cut -d. -f1))"
    else
        log "  $scenario: plots MISSING"
    fi
done
