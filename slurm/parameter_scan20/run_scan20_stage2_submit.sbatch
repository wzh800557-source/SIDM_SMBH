#!/bin/bash
#SBATCH --job-name=sidm20_stage2
#SBATCH --partition=mit_quicktest
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=slurm-%x-%A_%a.out

set -euo pipefail
PYTHON_BIN=${PYTHON_BIN:-python3}
export PYTHON_BIN
: "${SBATCH_ROOT:?missing SBATCH_ROOT}"
: "${SCAN_ROOT:?missing SCAN_ROOT}"
: "${SCAN_PKG:?missing SCAN_PKG}"
: "${PROD_PKG:?missing PROD_PKG}"
: "${SOURCE_PROFILE:?missing SOURCE_PROFILE}"

EXPORTS="ALL,SCAN_ROOT=$SCAN_ROOT,SCAN_PKG=$SCAN_PKG,PROD_PKG=$PROD_PKG,SOURCE_PROFILE=$SOURCE_PROFILE,BATCH_SIZE=24,PYTHON_BIN=$PYTHON_BIN,SBATCH_ROOT=$SBATCH_ROOT"
FLUID_JOB=$(sbatch --parsable \
  --array=0-49%25 \
  --export="$EXPORTS" \
  "$SBATCH_ROOT/run_scan20_fluid.sbatch")
AGG_JOB=$(sbatch --parsable \
  --dependency="afterok:$FLUID_JOB" \
  --export="$EXPORTS" \
  "$SBATCH_ROOT/run_scan20_aggregate.sbatch")

cat > "$SCAN_ROOT/STAGE2_SUBMISSION.txt" <<EOF
SIDM_SMBH_SCAN20_STAGE2_SUBMITTED $(date -u +%Y-%m-%dT%H:%M:%SZ)
FLUID_JOB=$FLUID_JOB
AGG_JOB=$AGG_JOB
BATCH_TASKS=50
CASES_PER_TASK=24
EOF
cat "$SCAN_ROOT/STAGE2_SUBMISSION.txt"
