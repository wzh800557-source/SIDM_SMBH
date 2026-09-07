#!/bin/bash
set -euo pipefail

REPO_ROOT=${REPO_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}
SCAN_PKG=${SCAN_PKG:-$REPO_ROOT/parameter_scan20}
PROD_PKG=${PROD_PKG:-$REPO_ROOT/fp_solver}
SCAN_ROOT=${SCAN_ROOT:?set SCAN_ROOT to a new output directory}
SOURCE_PROFILE=${SOURCE_PROFILE:-$REPO_ROOT/data/production/input/prof_deep.txt}
SBATCH_ROOT=$REPO_ROOT/slurm/parameter_scan20
PYTHON_BIN=${PYTHON_BIN:-python3}
export PYTHON_BIN

for path in \
  "$SCAN_PKG/scan_grid.json" \
  "$SBATCH_ROOT/run_scan20_prep.sbatch" \
  "$SBATCH_ROOT/run_scan20_stage2_submit.sbatch" \
  "$SBATCH_ROOT/run_scan20_fluid.sbatch" \
  "$SBATCH_ROOT/run_scan20_aggregate.sbatch"; do
  test -s "$path"
done
test -s "$SOURCE_PROFILE"
test -s "$PROD_PKG/prepare_bh_remapped_profile.py"
test -s "$PROD_PKG/hydrostatic_bridge.py"
test ! -e "$SCAN_ROOT" || { echo "refusing existing scan root" >&2; exit 2; }
mkdir -p "$SCAN_ROOT" "$REPO_ROOT/logs"
cd "$REPO_ROOT/logs"

BASE_EXPORTS="ALL,SCAN_ROOT=$SCAN_ROOT,SCAN_PKG=$SCAN_PKG,PROD_PKG=$PROD_PKG,SOURCE_PROFILE=$SOURCE_PROFILE,BATCH_SIZE=24,PYTHON_BIN=$PYTHON_BIN,SBATCH_ROOT=$SBATCH_ROOT"

# Batch 24 physical cases into each Slurm task.  This keeps the live job count
# below the ORCD QOS submission limit while retaining parallel execution.
PREP_JOB=$(sbatch --parsable \
  --array=0-49%25 \
  --export="$BASE_EXPORTS" \
  "$SBATCH_ROOT/run_scan20_prep.sbatch")
STAGE2_JOB=$(sbatch --parsable \
  --dependency="afterok:$PREP_JOB" \
  --export="$BASE_EXPORTS" \
  "$SBATCH_ROOT/run_scan20_stage2_submit.sbatch")

cat > "$SCAN_ROOT/SUBMISSION.txt" <<EOF
SIDM_SMBH_SCAN20_SUBMITTED $(date -u +%Y-%m-%dT%H:%M:%SZ)
SCAN_ROOT=$SCAN_ROOT
SCAN_PKG=$SCAN_PKG
PREP_JOB=$PREP_JOB
STAGE2_SUBMIT_JOB=$STAGE2_JOB
GRID=3x20x20
CASES=1200
BATCH_TASKS=50
CASES_PER_TASK=24
SCOPE=exact interface and black-hole-aware fluid luminosity; no dense FP capture-current inference
EOF
cat "$SCAN_ROOT/SUBMISSION.txt"
