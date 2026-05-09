#!/usr/bin/env bash
# run_remap_s5feats.sh — re-map source_5 using scenario_5-aligned features (1040 real features)
# Atlases are already trained; this only re-runs the mapping + metrics + report stages.
# Usage: bash scripts/run_remap_s5feats.sh [--dry-run]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
LOGS="$ROOT/logs"
RUNS="$ROOT/papermill_runs"
NBS="$ROOT/notebooks"
mkdir -p "$LOGS" "$RUNS" "$ROOT/results/acute" "$ROOT/results/umaps"

DRY_RUN=${1:-""}
LOG="$LOGS/remap_s5feats_$(date +%Y%m%d_%H%M%S).log"
echo "[$(date +%T)] run_remap_s5feats.sh starting — log: $LOG"

log()       { echo "[$(date +%T)] $*" | tee -a "$LOG"; }
milestone() { touch "$LOGS/MILESTONE_${1}"; log "MILESTONE: $1"; }

run_pm() {
    local env="$1" nb="$2" tag="$3"; shift 3
    local nb_ipynb="$RUNS/${tag}_input.ipynb"
    local out="$RUNS/${tag}.ipynb"
    if [[ "$DRY_RUN" == "--dry-run" ]]; then
        log "DRY-RUN: pixi run -e $env papermill $nb $out $*"
        return 0
    fi
    log "START $tag"
    pixi run --manifest-path "$ROOT/pixi.toml" -e "$env" \
        jupytext --to notebook --output "$nb_ipynb" "$NBS/$nb" \
        2>&1 | tee -a "$LOG"
    pixi run --manifest-path "$ROOT/pixi.toml" -e "$env" \
        papermill "$nb_ipynb" "$out" "$@" \
        2>&1 | tee -a "$LOG"
    log "DONE  $tag"
}

# ── Clear stale outputs ────────────────────────────────────────────────────────

log "Clearing stale source_5 mapped h5ads and metrics CSV"
rm -f "$ROOT/data/mapped/source_5_"*.h5ad
rm -f "$ROOT/results/umaps/source_5_"*.png
METRICS="$ROOT/results/acute/source_5_metrics.csv"
rm -f "$METRICS"

# ── Symphony ──────────────────────────────────────────────────────────────────

log "=== Symphony ==="
for arm in tplus tminus tminus_matched; do
    run_pm symphony 31_map_symphony.py "nb31_s5feats_${arm}_symphony" \
        -p QUERY_SOURCE source_5 -p T_ARM "$arm" \
        -p REFERENCE_ATLAS scenario_5_full
done

for arm in tplus tminus tminus_matched; do
    run_pm symphony 40_metrics.py "nb40_s5feats_${arm}_symphony" \
        -p QUERY_SOURCE source_5 -p T_ARM "$arm" \
        -p PARADIGM symphony -p REFERENCE_ATLAS scenario_5_full
done

milestone "SYMPHONY_DONE"

# ── scVI ──────────────────────────────────────────────────────────────────────

log "=== scVI ==="
for arm in tplus tminus tminus_matched; do
    run_pm scvi 11_map_scvi.py "nb11_s5feats_${arm}_scvi" \
        -p QUERY_SOURCE source_5 -p T_ARM "$arm" \
        -p REFERENCE_ATLAS scenario_5_full
done

for arm in tplus tminus tminus_matched; do
    run_pm symphony 40_metrics.py "nb40_s5feats_${arm}_scvi" \
        -p QUERY_SOURCE source_5 -p T_ARM "$arm" \
        -p PARADIGM scvi -p REFERENCE_ATLAS scenario_5_full
done

milestone "SCVI_DONE"

# ── scPoli ────────────────────────────────────────────────────────────────────

log "=== scPoli ==="
for arm in tplus tminus tminus_matched; do
    run_pm scpoli 21_map_scpoli.py "nb21_s5feats_${arm}_scpoli" \
        -p QUERY_SOURCE source_5 -p T_ARM "$arm" \
        -p REFERENCE_ATLAS scenario_5_full \
        -p REFERENCE_SCENARIO scenario_5
done

for arm in tplus tminus tminus_matched; do
    run_pm symphony 40_metrics.py "nb40_s5feats_${arm}_scpoli" \
        -p QUERY_SOURCE source_5 -p T_ARM "$arm" \
        -p PARADIGM scpoli -p REFERENCE_ATLAS scenario_5_full
done

milestone "SCPOLI_DONE"

# ── UMAPs + report ────────────────────────────────────────────────────────────

log "=== UMAPs ==="
pixi run --manifest-path "$ROOT/pixi.toml" -e symphony \
    python "$ROOT/scripts/generate_umaps.py" 2>&1 | tee -a "$LOG"

milestone "UMAPS_DONE"

log "=== Report ==="
pixi run --manifest-path "$ROOT/pixi.toml" -e symphony \
    python "$ROOT/scripts/generate_report.py" 2>&1 | tee -a "$LOG"

run_pm symphony 41_plot_acute.py nb41_s5feats_source5 \
    -p QUERY_SOURCE source_5

milestone "ALL_DONE"

log "=== FINISHED ==="
echo "  Metrics CSV:   $METRICS"
echo "  Summary table: $ROOT/results/source_5_summary.csv"
echo "  Bar chart:     $ROOT/results/source_5_barplot.pdf"
echo "  UMAPs:         $ROOT/results/umaps/"
echo "  Full log:      $LOG"
