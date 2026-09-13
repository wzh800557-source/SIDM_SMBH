#!/bin/bash
#SBATCH --job-name=sidm_e_gate
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
set -euo pipefail

: "${ENERGY_SUMMARY:?capture-energy summary required}"
: "${ENERGY_GATE:?output gate path required}"
: "${PKG:?repository fp_solver directory required}"
PYTHON_BIN=${PYTHON_BIN:-python3}
set +e
"$PYTHON_BIN" -u "$PKG/energy_current_gate.py" \
  --summary "$ENERGY_SUMMARY" --out "$ENERGY_GATE"
RC=$?
set -e
if [[ $RC != 0 && $RC != 3 ]]; then exit "$RC"; fi
# A failed scientific gate is a completed audit, not an infrastructure error.
exit 0

