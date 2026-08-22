#!/usr/bin/env bash
# ============================================================================
# rds -> .h5ad pipeline runner (steps 1-3, per analysis/mouse2human.md)
#
# Usage:
#   ./run.sh <tissue_folder> [sample_col] [ct_col] [max_per_group] [assay]
#
#   <tissue_folder>  folder name inside rds_data/ (e.g. PD_LN)
#   sample_col       default: samples2          (the 01 log lists your columns)
#   ct_col           default: final_cell.type
#   max_per_group    default: 200
#   assay            default: RNA  (do NOT pass SCT: its "counts" layer is
#                sctransform-corrected, which silently corrupts tokenization)
#
# Outputs land in $ADPD_ROOT/<tissue>/{logs,export,h5ad}.
# Every step is cached: re-running reuses existing outputs where present.
# Read logs/01_inspect.log after a first run to sanity-check columns,
# gene name format and count integrality.
# ============================================================================
set -euo pipefail

PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TISSUE="${1:?usage: ./run.sh <tissue_folder> [sample_col] [ct_col] [max_per_group] [assay]}"
SAMPLE_COL="${2:-samples2}"
CT_COL="${3:-final_cell.type}"
MAX_PER_GRP="${4:-200}"
ASSAY="${5:-RNA}"

. "$PIPE_DIR/env.sh"
[ -n "${GENEFORMER_DIR:-}" ] || { echo "GENEFORMER_DIR not set (see env.sh)" >&2; exit 1; }

RDS_DIR="$PIPE_DIR/rds_data/$TISSUE"
if [ ! -d "$RDS_DIR" ]; then
  echo "no such tissue folder: $RDS_DIR" >&2
  echo "available: $(ls "$PIPE_DIR/rds_data" 2>/dev/null | tr '\n' ' ')" >&2
  exit 1
fi
RDS="$(find "$RDS_DIR" -maxdepth 2 -name '*.rds' | sort | head -1)"
[ -n "$RDS" ] || { echo "no .rds file in $RDS_DIR" >&2; exit 1; }

ROOT="$ADPD_ROOT/$TISSUE"
mkdir -p "$ROOT"/{logs,export,h5ad}

echo "tissue : $TISSUE"
echo "rds    : $RDS"
echo "cols   : sample=$SAMPLE_COL celltype=$CT_COL max_per_group=$MAX_PER_GRP assay=$ASSAY"
echo "output : $ROOT"

echo
echo "=== [1/3] 01_inspect.R ==="
"$RSCRIPT" "$PIPE_DIR/01_inspect.R" "$RDS" 2>&1 | tee "$ROOT/logs/01_inspect.log"

# cheap pre-check before 02 re-reads the multi-GB object:
# both columns must appear in the 01 log's meta.data listing
for col in "$SAMPLE_COL" "$CT_COL"; do
  grep -qFw "$col" "$ROOT/logs/01_inspect.log" || {
    echo "column '$col' not found in 01_inspect.log — check the meta.data table and re-run:" >&2
    echo "  ./run.sh $TISSUE <sample_col> $col ..." >&2
    exit 1
  }
done
if grep -q "integer-like:  FALSE" "$ROOT/logs/01_inspect.log"; then
  echo "WARNING: 01 log reports non-integer counts (check assay/layer)" >&2
fi

echo
echo "=== [2/3] 02_export.R ==="
"$RSCRIPT" "$PIPE_DIR/02_export.R" "$RDS" "$ROOT/export/data" \
  "$SAMPLE_COL" "$CT_COL" "$MAX_PER_GRP" "$ASSAY" 2>&1 | tee "$ROOT/logs/02_export.log"

echo
echo "=== [3/3] 03_map_and_build_h5ad.py ==="
"$PYTHON" "$PIPE_DIR/03_map_and_build_h5ad.py" \
  "$ROOT/export/data" "$ROOT/h5ad/$TISSUE.h5ad" "$GENEFORMER_DIR/geneformer" \
  2>&1 | tee "$ROOT/logs/03_map.log"

echo
echo "DONE: $ROOT/h5ad/$TISSUE.h5ad"
