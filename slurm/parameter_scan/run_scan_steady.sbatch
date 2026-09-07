#!/bin/bash
#SBATCH --job-name=sidm_scan_ss
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=192G
#SBATCH --time=16:00:00
#SBATCH --array=0-26%6
#SBATCH --output=slurm-%x-%A_%a.out

set -euo pipefail
export OMP_NUM_THREADS=1
PYTHON_BIN=${PYTHON_BIN:-python3}

: "${SCAN_ROOT:?missing SCAN_ROOT}"
: "${SCAN_PKG:?missing SCAN_PKG}"
: "${PROD_PKG:?missing PROD_PKG}"
export PYTHONPATH="$SCAN_PKG:$PROD_PKG:${PYTHONPATH:-}"

"$PYTHON_BIN" -u "$SCAN_PKG/run_scan_steady.py" \
  --grid "$SCAN_PKG/scan_grid.json" \
  --index "$SLURM_ARRAY_TASK_ID" \
  --root "$SCAN_ROOT" \
  --production-package "$PROD_PKG" \
  --scan-package "$SCAN_PKG"

