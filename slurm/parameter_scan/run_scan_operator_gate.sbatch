#!/bin/bash
#SBATCH --job-name=sidm_scan_op_gate
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail
export OMP_NUM_THREADS=1
PYTHON_BIN=${PYTHON_BIN:-python3}

: "${SCAN_ROOT:?missing SCAN_ROOT}"
: "${SCAN_PKG:?missing SCAN_PKG}"
export PYTHONPATH="$SCAN_PKG:${PYTHONPATH:-}"

"$PYTHON_BIN" -u "$SCAN_PKG/gate_scan_operators.py" \
  --grid "$SCAN_PKG/scan_grid.json" \
  --root "$SCAN_ROOT"
