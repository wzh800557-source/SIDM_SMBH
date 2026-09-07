#!/bin/bash
#SBATCH --job-name=sidm20_prep
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%x-%A_%a.out

set -euo pipefail
export OMP_NUM_THREADS=1
PYTHON_BIN=${PYTHON_BIN:-python3}
export PYTHON_BIN

: "${SCAN_ROOT:?missing SCAN_ROOT}"
: "${SCAN_PKG:?missing SCAN_PKG}"
: "${PROD_PKG:?missing PROD_PKG}"
: "${SOURCE_PROFILE:?missing SOURCE_PROFILE}"
export PYTHONPATH="$SCAN_PKG:$PROD_PKG:${PYTHONPATH:-}"

BATCH_SIZE=${BATCH_SIZE:-24}
START=$((SLURM_ARRAY_TASK_ID * BATCH_SIZE))
END=$((START + BATCH_SIZE - 1))

for INDEX in $(seq "$START" "$END"); do
  "$PYTHON_BIN" -u "$SCAN_PKG/prepare_interface_case.py" \
    --grid "$SCAN_PKG/scan_grid.json" \
    --index "$INDEX" \
    --source-profile "$SOURCE_PROFILE" \
    --root "$SCAN_ROOT" \
    --production-package "$PROD_PKG" \
    --scan-package "$SCAN_PKG"
done
