#!/usr/bin/env bash
# Resume S1 and S2 HPO to fill in missing Seurat v4/v5 corrections.
# S1: missing seurat_rpca_v4 (HPO + correction). v5 variants skipped (0 COMPLETE trials).
# S2: missing all 4 v4/v5 seurat (HPO + correction).
#
# Uses explicit file targets + --rerun-triggers mtime:
#  - explicit targets force Snakemake to recheck file existence (unlike rule `all`,
#    which short-circuits when the top PDF exists with recent mtime)
#  - --rerun-triggers mtime prevents spurious rebuilds from script-hash changes
#
# NOTE: --touch --forceall (as in resume_pipeline.sh) corrupts metadata when outputs
# are missing — it marks rules "complete" without creating files. Do not use here.
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

LOG_DIR="outputs/overnight_$(date +%Y%m%d_%H%M)"
mkdir -p "$LOG_DIR"
MASTER_LOG="$LOG_DIR/master.log"

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

# Clear stale locks from earlier runs
rm -f .snakemake/locks/0.input.lock .snakemake/locks/0.output.lock 2>/dev/null
mkdir -p .snakemake/locks

log "=========================================="
log "RESUME S1+S2 HPO START"
log "Node: $(hostname), GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo unknown)"
log "SLURM_JOB_ID: ${SLURM_JOB_ID:-none}"
log "=========================================="

# ----- S1: seurat_rpca_v4 (HPO + correction) + downstream refresh -----
S1_TARGETS=(
    "outputs/scenario_1/mad_int_featselect_seurat_rpca_v4.parquet"
    "outputs/scenario_1/plots/results_table.pdf"
)
run_step "s1_hpo" "pixi run -- bash -c 'PATH=\$(dirname \$(which python)):\$PATH snakemake -c20 --rerun-triggers mtime --keep-going --configfile inputs/conf/shared.json inputs/conf/scenario_1.json --sdm apptainer --apptainer-args=--nv --resources nvidia_gpu=1 -- ${S1_TARGETS[*]}'"

# ----- S2: all 4 v4/v5 seurat (HPO + correction) + downstream refresh -----
S2_TARGETS=(
    "outputs/scenario_2/mad_int_featselect_seurat_cca_v4.parquet"
    "outputs/scenario_2/mad_int_featselect_seurat_rpca_v4.parquet"
    "outputs/scenario_2/mad_int_featselect_seurat_cca_v5.parquet"
    "outputs/scenario_2/mad_int_featselect_seurat_rpca_v5.parquet"
    "outputs/scenario_2/plots/results_table.pdf"
)
run_step "s2_hpo" "pixi run -- bash -c 'PATH=\$(dirname \$(which python)):\$PATH snakemake -c20 --rerun-triggers mtime --keep-going --configfile inputs/conf/shared.json inputs/conf/scenario_2.json --sdm apptainer --apptainer-args=--nv --resources nvidia_gpu=1 -- ${S2_TARGETS[*]}'"

log "=========================================="
log "RESUME S1+S2 HPO COMPLETE"
log "=========================================="

for scenario in scenario_1 scenario_2; do
    if [ -f "outputs/$scenario/plots/results_table.pdf" ]; then
        log "  $scenario: plots OK ($(stat -c '%y' "outputs/$scenario/plots/results_table.pdf" 2>/dev/null | cut -d. -f1))"
    else
        log "  $scenario: plots MISSING"
    fi
    for m in seurat_cca_v4 seurat_rpca_v4 seurat_cca_v5 seurat_rpca_v5; do
        if [ -f "outputs/$scenario/mad_int_featselect_${m}.parquet" ]; then
            log "  $scenario: $m parquet OK"
        else
            log "  $scenario: $m parquet MISSING"
        fi
    done
done
