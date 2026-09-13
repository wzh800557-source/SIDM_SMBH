#!/bin/bash
#SBATCH --job-name=sidm_radial_kin
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=0-25%8
set -euo pipefail

: "${RUN_ROOT:?new output directory required}"
: "${PKG:?repository fp_solver directory required}"
PYTHON_BIN=${PYTHON_BIN:-python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

KNS=(0.01 0.02 0.04 0.08 0.16 0.32 0.64 1.28 2.56 5.12 10.24 20.48 40.96)
INDEX=$((SLURM_ARRAY_TASK_ID % 13))
LEVEL=$((SLURM_ARRAY_TASK_ID / 13))
KN=${KNS[$INDEX]}
if [[ $LEVEL == 0 ]]; then
  LABEL=medium
  NR=40; NVR=56; NJ=32; DT=0.004
else
  LABEL=fine
  NR=56; NVR=72; NJ=44; DT=0.0025
fi
OUT="$RUN_ROOT/${LABEL}_kn${KN}.json"
LOG="$RUN_ROOT/${LABEL}_kn${KN}.log"
STATUS="$RUN_ROOT/${LABEL}_kn${KN}.status"
mkdir -p "$RUN_ROOT"
trap 'echo "FAILED line=$LINENO status=$?" > "$STATUS"' ERR
printf 'RUNNING %s level=%s Kn=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$LABEL" "$KN" > "$STATUS"
set +e
"$PYTHON_BIN" -u "$PKG/radial_kinetic_transition.py" \
  --kn "$KN" --gamma 1.6666666666666667 \
  --r-min 0.03 --r-max 4.0 --velocity-max 20.0 \
  --radial-bins "$NR" --radial-velocity-bins "$NVR" --tangential-bins "$NJ" \
  --dt "$DT" --t-end 16.0 --output-dt 0.1 --wall-seconds 40000 \
  --out "$OUT" > "$LOG" 2>&1
RC=$?
set -e
if [[ $RC != 0 && $RC != 2 ]]; then
  exit "$RC"
fi
printf 'DIAGNOSTICS_COMPLETE %s solver_status=%s\n' \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$RC" > "$STATUS"

