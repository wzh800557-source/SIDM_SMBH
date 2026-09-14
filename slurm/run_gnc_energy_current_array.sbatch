#!/bin/bash
#SBATCH --job-name=gnc_e_current
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --array=0-3%4
set -euo pipefail

: "${ENSEMBLE_MANIFEST:?ensemble manifest required}"
: "${GNC_BUILD:?compiled patched GNC tree required}"
: "${PKG:?repository fp_solver directory required}"
PYTHON_BIN=${PYTHON_BIN:-python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

readarray -t MEMBER < <("$PYTHON_BIN" - "$ENSEMBLE_MANIFEST" "$SLURM_ARRAY_TASK_ID" <<'PY'
import json,sys
x=json.load(open(sys.argv[1])); m=x['members'][int(sys.argv[2])]
print(m['run_dir']); print(m['ranks']); print(x['dt_tnr']); print(x['total_tnr'])
print(x['burn_in_tnr_for_gate']); print(x['sigma0_over_m_cm2_g'])
PY
)
RUN_DIR=${MEMBER[0]}
RANKS=${MEMBER[1]}
DT_TNR=${MEMBER[2]}
BURN_TNR=${MEMBER[4]}
SIGMA0=${MEMBER[5]}
STATUS="$RUN_DIR/STATUS.txt"
trap 'rc=$?; echo "FAILED line=$LINENO status=$rc" > "$STATUS"; exit $rc' ERR
if [[ -e "$RUN_DIR/output" ]]; then
  echo "refusing to overwrite an existing GNC output tree" >&2
  exit 2
fi
cp "$GNC_BUILD/main/ini" "$GNC_BUILD/main/main" "$GNC_BUILD/main/pro" "$RUN_DIR/"
cp -a "$GNC_BUILD/examples/1comp_lc/output" "$RUN_DIR/output"
cd "$RUN_DIR"
printf 'INITIALIZING %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$STATUS"
mpirun --bind-to none -np "$RANKS" ./ini > ini.log 2>&1
printf 'EVOLVING %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$STATUS"
mpirun --bind-to none -np "$RANKS" ./main > main.log 2>&1
printf 'POSTPROCESSING %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$STATUS"
./pro > pro.log 2>&1
printf 'ANALYZING %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$STATUS"
set +e
"$PYTHON_BIN" -u "$PKG/analyze_long_run.py" "$RUN_DIR" \
  --dt-tnr "$DT_TNR" --window-tnr 1.0 --burn-in-tnr "$BURN_TNR" \
  --plateau-tolerance 0.25 --inventory-tolerance 0.25 \
  --minimum-boundary-effective-captures 10 \
  --sigma-over-m-cm2-g "$SIGMA0" --kernel-compatible \
  > analysis.log 2>&1
RC=$?
set -e
if [[ $RC != 0 && $RC != 4 ]]; then exit "$RC"; fi
printf 'DIAGNOSTICS_COMPLETE %s analyzer_status=%s\n' \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$RC" > "$STATUS"
