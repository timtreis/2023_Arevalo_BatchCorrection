#!/usr/bin/env bash
# run_remaining.sh — resume from scVI mapping onwards (Symphony already done)
# Skips scVI atlas (already trained), skips Symphony entirely.
# Idempotent: skips mapped h5ads that already exist.
# Usage: bash scripts/run_remaining.sh [--dry-run]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
LOGS="$ROOT/logs"
RUNS="$ROOT/papermill_runs"
NBS="$ROOT/notebooks"
mkdir -p "$LOGS" "$RUNS" "$ROOT/results/acute" "$ROOT/results/umaps"

DRY_RUN=${1:-""}
LOG="$LOGS/run_remaining_$(date +%Y%m%d_%H%M%S).log"
echo "[$(date +%T)] run_remaining.sh starting — log: $LOG"

# ── helpers ────────────────────────────────────────────────────────────────────

log() { echo "[$(date +%T)] $*" | tee -a "$LOG"; }

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

mapped_exists() { [[ -f "$ROOT/data/mapped/$1" ]]; }
model_exists()  { [[ -f "$ROOT/models/$1/model.pt" ]]; }
atlas_exists()  { [[ -d "$ROOT/models/$1" ]]; }

milestone() {
    local name="$1"
    touch "$LOGS/MILESTONE_${name}"
    log "MILESTONE: $name"
}

# ── Clean only metrics CSV (preserve mapped h5ads — mapping is idempotent) ─────

log "Trimming metrics CSV to symphony rows only (preserving mapped h5ads)"
# Remove old UMAP outputs so they are regenerated
rm -f "$ROOT/results/umaps/source_5_"*.png
# Trim metrics CSV to only symphony rows (keep header + 3 symphony lines)
METRICS="$ROOT/results/acute/source_5_metrics.csv"
if [[ -f "$METRICS" ]]; then
    head -1 "$METRICS" > "${METRICS}.tmp"
    grep -E ",symphony,|,scvi," "$METRICS" >> "${METRICS}.tmp" || true
    mv "${METRICS}.tmp" "$METRICS"
    log "Trimmed metrics CSV — kept symphony + scvi rows"
fi

# ── WAVE 1: scVI (GPU) ─────────────────────────────────────────────────────────

log "=== WAVE 1: scVI ==="

if ! model_exists "scvi_atlas_scenario_5_full"; then
    run_pm scvi 10_train_atlas_scvi.py nb10_scvi_full \
        -p QUERY_SOURCE "" -p REFERENCE_SCENARIO scenario_5
else
    log "SKIP scVI atlas (already trained)"
fi

for arm in tplus tminus tminus_matched; do
    if ! mapped_exists "source_5_${arm}_scvi.h5ad"; then
        run_pm scvi 11_map_scvi.py "nb11_s5_${arm}_scvi" \
            -p QUERY_SOURCE source_5 -p T_ARM "$arm" \
            -p REFERENCE_ATLAS scenario_5_full
    else
        log "SKIP scVI map $arm (already exists)"
    fi
done

for arm in tplus tminus tminus_matched; do
    if ! mapped_exists "source_5_${arm}_scvi.h5ad"; then
        log "SKIP scVI metrics $arm (mapped h5ad missing)"
    elif grep -q "source_5,${arm},scvi," "$METRICS" 2>/dev/null; then
        log "SKIP scVI metrics $arm (row already in CSV)"
    else
        run_pm symphony 40_metrics.py "nb40_s5_${arm}_scvi" \
            -p QUERY_SOURCE source_5 -p T_ARM "$arm" \
            -p PARADIGM scvi -p REFERENCE_ATLAS scenario_5_full
    fi
done

milestone "SCVI_DONE"

# ── WAVE 2: scPoli (GPU, slowest) ─────────────────────────────────────────────

log "=== WAVE 2: scPoli ==="

if ! atlas_exists "scpoli_atlas_scenario_5_full"; then
    run_pm scpoli 20_train_atlas_scpoli.py nb20_scpoli_full \
        -p QUERY_SOURCE "" -p REFERENCE_SCENARIO scenario_5 \
        -p T2_ONLY_PROTOTYPES true
else
    log "SKIP scPoli atlas (already exists)"
fi

for arm in tplus tminus tminus_matched; do
    if ! mapped_exists "source_5_${arm}_scpoli.h5ad"; then
        run_pm scpoli 21_map_scpoli.py "nb21_s5_${arm}_scpoli" \
            -p QUERY_SOURCE source_5 -p T_ARM "$arm" \
            -p REFERENCE_ATLAS scenario_5_full \
            -p REFERENCE_SCENARIO scenario_5
    else
        log "SKIP scPoli map $arm (already exists)"
    fi
done

for arm in tplus tminus tminus_matched; do
    if mapped_exists "source_5_${arm}_scpoli.h5ad"; then
        run_pm symphony 40_metrics.py "nb40_s5_${arm}_scpoli" \
            -p QUERY_SOURCE source_5 -p T_ARM "$arm" \
            -p PARADIGM scpoli -p REFERENCE_ATLAS scenario_5_full
    else
        log "SKIP scPoli metrics $arm (mapped h5ad missing)"
    fi
done

milestone "SCPOLI_DONE"

# ── WAVE 3: UMAPs + report ────────────────────────────────────────────────────

log "=== WAVE 3: UMAPs + report ==="
pixi run --manifest-path "$ROOT/pixi.toml" -e symphony \
    python "$ROOT/scripts/generate_umaps.py" 2>&1 | tee -a "$LOG"

milestone "UMAPS_DONE"

pixi run --manifest-path "$ROOT/pixi.toml" -e symphony \
    python "$ROOT/scripts/generate_report.py" 2>&1 | tee -a "$LOG"

run_pm symphony 41_plot_acute.py nb41_source_5 \
    -p QUERY_SOURCE source_5

milestone "ALL_DONE"

log "=== FINISHED ==="
echo "  Metrics CSV:   $ROOT/results/acute/source_5_metrics.csv"
echo "  Summary table: $ROOT/results/source_5_summary.csv"
echo "  Bar chart:     $ROOT/results/source_5_barplot.pdf"
echo "  UMAPs:         $ROOT/results/umaps/"
echo "  Full log:      $LOG"
