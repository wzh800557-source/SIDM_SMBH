#!/bin/bash
#SBATCH --job-name=sidm_energy_mom
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --array=0-3%2
set -euo pipefail
PYTHON_BIN=${PYTHON_BIN:-python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
: "${ROOT:?production input root required}"
: "${RUN_ROOT:?new output directory required}"
: "${PKG:?new code directory required}"
E=64
if ((SLURM_ARRAY_TASK_ID>=2)); then E=96; fi
SEED=314159
if ((SLURM_ARRAY_TASK_ID%2)); then SEED=271828; fi
TAG=e${E}_j257_s${SEED}
OUT="$RUN_ROOT/$TAG"
mkdir "$OUT"
trap 'echo "FAILED line=$LINENO status=$?" > "$OUT/STATUS.txt"' ERR
printf 'SAMPLING %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$OUT/STATUS.txt"
RB=$("$PYTHON_BIN" - "$ROOT/bridge/boundary_criteria.json" <<'PY'
import json,sys
print('%.16e'%json.load(open(sys.argv[1]))['criteria']['local_circular']['r_boundary_pc'])
PY
)
"$PYTHON_BIN" -u "$PKG/finite_angle_capture.py" \
 --profile "$ROOT/bridge/bridged_profile.txt" \
 --df-table "$ROOT/normalized/f_evolved_normalized.txt" \
 --normalization-json "$ROOT/normalized/df_normalization.json" \
 --out-json "$OUT/run.json" --out-csv "$OUT/run.csv" --out-ej-operator "$OUT/operator.npz" \
 --mbh 4e6 --sigma0-over-m 100 --w-kms 80 --r-outer "$RB" \
 --radial-bins 16 --pairs-per-bin 1600000 --ej-energy-bins "$E" \
 --ej-angular-grid j257 --seed "$SEED" > "$OUT/sampling.log" 2>&1
printf 'SOLVING %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$OUT/STATUS.txt"
# Retain diagnostics on a population or moment failure; never launch a fluid
# evolution from a process exit code or mass-only success marker.
set +e
"$PYTHON_BIN" -u "$PKG/solve_ej_steady_state.py" --operator "$OUT/operator.npz" \
 --out-json "$OUT/steady.json" --out-csv "$OUT/steady.csv" > "$OUT/steady.log" 2>&1
SOLVE_RC=$?
set -e
if [[ $SOLVE_RC != 0 && $SOLVE_RC != 2 ]]; then exit "$SOLVE_RC"; fi
printf 'PROJECTING %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$OUT/STATUS.txt"
for MODEL in direct immediate; do
 set +e
 "$PYTHON_BIN" -u "$PKG/reproject_depleted_df.py" --operator "$OUT/operator.npz" \
  --steady-json "$OUT/steady.json" --steady-csv "$OUT/steady.csv" \
  --df-table "$ROOT/normalized/f_evolved_normalized.txt" \
  --normalization-json "$ROOT/normalized/df_normalization.json" \
  --profile "$ROOT/bridge/bridged_profile.txt" --model "$MODEL" \
  --radii 48 --quadrature-order 8 \
  --out-json "$OUT/projection_${MODEL}.json" --out-csv "$OUT/projection_${MODEL}.csv" \
  > "$OUT/projection_${MODEL}.log" 2>&1
 RC=$?
 set -e
 if [[ $RC != 0 && $RC != 2 ]]; then exit "$RC"; fi
done
"$PYTHON_BIN" - "$OUT" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]);s=json.loads((p/'steady.json').read_text())
x={'status':'DIAGNOSTICS_COMPLETE','population_status':s['status'],
   'energy_moments':{m:s['models'][m]['energy_moment_audit'] for m in ('direct','immediate')},
   'projections':{m:json.loads((p/('projection_'+m+'.json')).read_text()) for m in ('direct','immediate')},
   'fluid_coupling_authorized':False}
(p/'diagnostic_summary.json').write_text(json.dumps(x,indent=2)+'\n')
PY
printf 'DIAGNOSTICS_COMPLETE %s fluid_coupling_authorized=false\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$OUT/STATUS.txt"
