#!/usr/bin/env bash
# Create a defaults variant of a scenario by symlinking preprocessed data.
# Usage: ./scripts/create_defaults_variant.sh scenario_8
#
# This avoids re-running the expensive preprocessing pipeline (hours, ~30 GB)
# when you just want to run corrections with default parameters.

set -euo pipefail

SCENARIO="${1:?Usage: $0 <scenario_name>}"
DEFAULTS="${SCENARIO}_defaults"
SRC="outputs/${SCENARIO}"
DST="outputs/${DEFAULTS}"

if [ ! -d "$SRC" ]; then
    echo "Error: source directory $SRC does not exist"
    exit 1
fi

mkdir -p "$DST"

# Symlink all preprocessing outputs
for f in raw.parquet neg_stats.parquet variant_feats.parquet \
         mad.parquet mad_int.parquet mad_int_featselect.parquet \
         norm_stats.parquet outliers.parquet; do
    if [ -f "$SRC/$f" ]; then
        ln -sf "$(realpath "$SRC/$f")" "$DST/$f"
        echo "  Linked: $f"
    fi
done

echo "Done. Preprocessing outputs linked from $SRC to $DST."
echo "Run with: pixi run scenario-${SCENARIO##scenario_}-defaults"
