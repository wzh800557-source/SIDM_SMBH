#!/bin/bash
#SBATCH --job-name=sidm_energy_dg1
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --array=0-3%2
set -euo pipefail
PYTHON_BIN=${PYTHON_BIN:-python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
: "${ROOT:?production input root required}"
: "${RUN_ROOT:?new output root required}"
: "${PKG:?new solver directory required}"
E=32; if ((SLURM_ARRAY_TASK_ID>=2)); then E=64; fi
SEED=314159; if ((SLURM_ARRAY_TASK_ID%2)); then SEED=271828; fi
OUT="$RUN_ROOT/e${E}_fine_s${SEED}"
mkdir "$OUT"
trap 'rc=$?; echo "FAILED line=$LINENO status=$rc" > "$OUT/STATUS.txt"; exit $rc' ERR
status() { printf '%s %s\n' "$1" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$OUT/STATUS.txt"; }
RB=$("$PYTHON_BIN" - "$ROOT/bridge/boundary_criteria.json" <<'PY'
import json,sys
print('%.16e'%json.load(open(sys.argv[1]))['criteria']['local_circular']['r_boundary_pc'])
PY
)
status SAMPLING_DG1
"$PYTHON_BIN" -u "$PKG/finite_angle_capture.py" \
 --profile "$ROOT/bridge/bridged_profile.txt" --df-table "$ROOT/normalized/f_evolved_normalized.txt" \
 --normalization-json "$ROOT/normalized/df_normalization.json" \
 --out-json "$OUT/run.json" --out-csv "$OUT/run.csv" --out-ej-operator "$OUT/operator.npz" \
 --mbh 4e6 --sigma0-over-m 100 --w-kms 80 --r-outer "$RB" \
 --radial-bins 16 --pairs-per-bin 200000 --ej-energy-bins "$E" \
 --ej-angular-grid fine --ej-energy-order 1 --seed "$SEED" > "$OUT/sampling.log" 2>&1
status POPULATION_CONTROL
set +e
"$PYTHON_BIN" -u "$PKG/solve_ej_steady_state.py" --operator "$OUT/operator.npz" \
 --out-json "$OUT/population.json" --out-csv "$OUT/population.csv" > "$OUT/population.log" 2>&1
RC=$?
set -e
if [[ $RC != 0 && $RC != 2 ]]; then exit "$RC"; fi
status SOLVING_MASS_AND_ENERGY
set +e
"$PYTHON_BIN" -u "$PKG/solve_ej_energy_dg.py" --operator "$OUT/operator.npz" \
 --initial-population-csv "$OUT/population.csv" --continuation --maximum-evaluations 300 \
 --out-json "$OUT/steady_dg.json" --out-csv "$OUT/steady_dg.csv" > "$OUT/steady_dg.log" 2>&1
RC=$?
set -e
if [[ $RC != 0 && $RC != 2 ]]; then exit "$RC"; fi
"$PYTHON_BIN" - "$OUT" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]);old=json.loads((p/'population.json').read_text());new=json.loads((p/'steady_dg.json').read_text())
s={'status':'DG_PILOT_DIAGNOSTICS_COMPLETE','scope':'production-snapshot discretization validation, not a production-converged flux or fluid response',
   'population_control':old['models'],'dg_models':new['models'],'fluid_coupling_authorized':False}
(p/'diagnostic_summary.json').write_text(json.dumps(s,indent=2)+'\n')
PY
status DIAGNOSTICS_COMPLETE_NO_FLUID_COUPLING
