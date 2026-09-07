#!/bin/bash
set -euo pipefail

REPO_ROOT=${REPO_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}
SCAN_PKG=${SCAN_PKG:-$REPO_ROOT/parameter_scan}
PROD_PKG=${PROD_PKG:-$REPO_ROOT/fp_solver}
SCAN_ROOT=${SCAN_ROOT:?set SCAN_ROOT to a new output directory}
SOURCE_PROFILE=${SOURCE_PROFILE:?set SOURCE_PROFILE to the collapse snapshot}
PYTHON_BIN=${PYTHON_BIN:-python3}
export PYTHON_BIN

for path in "$SCAN_PKG/scan_grid.json" "$SOURCE_PROFILE"             "$PROD_PKG/finite_angle_capture.py"; do
  test -s "$path"
done
mkdir -p "$SCAN_ROOT" "$REPO_ROOT/logs"

SBATCH_ROOT="$REPO_ROOT/slurm/parameter_scan"
EXPORTS="ALL,SCAN_ROOT=$SCAN_ROOT,SCAN_PKG=$SCAN_PKG,PROD_PKG=$PROD_PKG,SOURCE_PROFILE=$SOURCE_PROFILE,PYTHON_BIN=$PYTHON_BIN"
cd "$REPO_ROOT/logs"
PREP_JOB=$(sbatch --parsable --export="$EXPORTS" "$SBATCH_ROOT/run_scan_prep.sbatch")
FLUID_JOB=$(sbatch --parsable --dependency="afterok:$PREP_JOB" --export="$EXPORTS" "$SBATCH_ROOT/run_scan_fluid_all.sbatch")
FP_JOB=$(sbatch --parsable --dependency="afterok:$PREP_JOB" --export="$EXPORTS" "$SBATCH_ROOT/run_scan_operator.sbatch")
GATE_JOB=$(sbatch --parsable --dependency="afterok:$FP_JOB" --export="$EXPORTS" "$SBATCH_ROOT/run_scan_operator_gate.sbatch")
STEADY_JOB=$(sbatch --parsable --dependency="afterok:$GATE_JOB" --export="$EXPORTS" "$SBATCH_ROOT/run_scan_steady.sbatch")
AGG_JOB=$(sbatch --parsable --dependency="afterok:$STEADY_JOB:$FLUID_JOB" --export="$EXPORTS" "$SBATCH_ROOT/run_scan_aggregate.sbatch")

cat > "$SCAN_ROOT/SUBMISSION.txt" <<EOF
PREP_JOB=$PREP_JOB
FLUID_JOB=$FLUID_JOB
FP_JOB=$FP_JOB
GATE_JOB=$GATE_JOB
STEADY_JOB=$STEADY_JOB
AGG_JOB=$AGG_JOB
EOF
cat "$SCAN_ROOT/SUBMISSION.txt"
