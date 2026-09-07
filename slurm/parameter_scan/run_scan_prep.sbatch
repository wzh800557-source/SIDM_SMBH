#!/bin/bash
#SBATCH --job-name=sidm_scan_prep
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-26%12
#SBATCH --output=slurm-%x-%A_%a.out

set -euo pipefail
export OMP_NUM_THREADS=1
PYTHON_BIN=${PYTHON_BIN:-python3}

: "${SCAN_ROOT:?missing SCAN_ROOT}"
: "${SCAN_PKG:?missing SCAN_PKG}"
: "${PROD_PKG:?missing PROD_PKG}"
: "${SOURCE_PROFILE:?missing SOURCE_PROFILE}"
export PYTHONPATH="$SCAN_PKG:$PROD_PKG:${PYTHONPATH:-}"

PROFILE_GUARD=()
if [[ -n "${EXPECTED_SOURCE_PROFILE_SHA256:-}" ]]; then
  PROFILE_GUARD=(--expected-source-sha256 "$EXPECTED_SOURCE_PROFILE_SHA256")
fi

"$PYTHON_BIN" -u "$SCAN_PKG/prepare_scan_case.py" \
  --grid "$SCAN_PKG/scan_grid.json" \
  --index "$SLURM_ARRAY_TASK_ID" \
  --source-profile "$SOURCE_PROFILE" \
  --root "$SCAN_ROOT" \
  --production-package "$PROD_PKG" \
  --scan-package "$SCAN_PKG" \
  "${PROFILE_GUARD[@]}"
