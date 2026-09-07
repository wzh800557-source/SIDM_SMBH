#!/bin/bash
#SBATCH --job-name=sidm_scan_agg
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail
export OMP_NUM_THREADS=1
PYTHON_BIN=${PYTHON_BIN:-python3}

: "${SCAN_ROOT:?missing SCAN_ROOT}"
: "${SCAN_PKG:?missing SCAN_PKG}"
export PYTHONPATH="$SCAN_PKG:${PYTHONPATH:-}"

mkdir -p "$SCAN_ROOT/aggregate"
"$PYTHON_BIN" -u "$SCAN_PKG/aggregate_scan.py" \
  --grid "$SCAN_PKG/scan_grid.json" \
  --root "$SCAN_ROOT" \
  --outdir "$SCAN_ROOT/aggregate"

"$PYTHON_BIN" -u "$SCAN_PKG/plot_scan_results.py" \
  --aggregate-json "$SCAN_ROOT/aggregate/parameter_scan_results.json" \
  --radial-csv "$SCAN_ROOT/aggregate/parameter_scan_radial_flux_profiles.csv" \
  --outdir "$SCAN_ROOT/aggregate"

printf 'SIDM_SMBH_PARAMETER_SCAN_COMPLETE %s ROOT=%s\n' \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$SCAN_ROOT" | tee "$SCAN_ROOT/SCAN_STATUS.txt"
