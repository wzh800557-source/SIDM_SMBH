#!/bin/bash
#SBATCH --job-name=sidm_kin_fit
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
set -euo pipefail

: "${RUN_ROOT:?run directory required}"
: "${PKG:?repository fp_solver directory required}"
PYTHON_BIN=${PYTHON_BIN:-python3}
set +e
"$PYTHON_BIN" -u "$PKG/calibrate_kinetic_hydro_bridge.py" \
  --run-root "$RUN_ROOT" --out "$RUN_ROOT/kinetic_hydro_calibration.json" \
  > "$RUN_ROOT/kinetic_hydro_calibration.log" 2>&1
RC=$?
set -e
if [[ $RC != 0 && $RC != 2 ]]; then exit "$RC"; fi
printf 'CALIBRATION_DIAGNOSTICS_COMPLETE status=%s %s\n' \
  "$RC" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$RUN_ROOT/AGGREGATE_STATUS.txt"

