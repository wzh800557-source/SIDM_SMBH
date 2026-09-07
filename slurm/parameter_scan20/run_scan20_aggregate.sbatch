#!/bin/bash
#SBATCH --job-name=sidm20_agg
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%x-%A_%a.out

set -euo pipefail
export OMP_NUM_THREADS=1
PYTHON_BIN=${PYTHON_BIN:-python3}
export PYTHON_BIN

: "${SCAN_ROOT:?missing SCAN_ROOT}"
: "${SCAN_PKG:?missing SCAN_PKG}"
export PYTHONPATH="$SCAN_PKG:${PYTHONPATH:-}"

mkdir -p "$SCAN_ROOT/aggregate"
"$PYTHON_BIN" -u "$SCAN_PKG/aggregate_scan20.py" \
  --grid "$SCAN_PKG/scan_grid.json" \
  --root "$SCAN_ROOT" \
  --outdir "$SCAN_ROOT/aggregate"

"$PYTHON_BIN" -u "$SCAN_PKG/plot_scan20_heatmaps.py" \
  --results-csv "$SCAN_ROOT/aggregate/scan20_results.csv" \
  --outdir "$SCAN_ROOT/aggregate"

printf 'SIDM_SMBH_SCAN20_COMPLETE %s ROOT=%s\n' \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$SCAN_ROOT" | tee "$SCAN_ROOT/SCAN20_STATUS.txt"
