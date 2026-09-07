#!/bin/bash
#SBATCH --job-name=gnc_res3_agg
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=192G
#SBATCH --time=16:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail
module load deprecated-modules
module load anaconda3/2022.05-x86_64
eval "$(conda shell.bash hook)"
conda activate fesc
export OMP_NUM_THREADS=1
PYTHON_BIN=${PYTHON_BIN:-python3}
: "${ROOT:?submit with ROOT=<production root>}"
PKG=${PKG:?set PKG to the fp_solver directory}
BASE_STEADY_ROOT=${BASE_STEADY_ROOT:-$ROOT/steady_v5}
STEADY_ROOT=${STEADY_ROOT:-$ROOT/steady_v6}
AGGREGATE_ROOT=${AGGREGATE_ROOT:-$ROOT/aggregate_v6}
AGGREGATE_STATUS=${AGGREGATE_STATUS:-$ROOT/AGGREGATE_V6_STATUS.txt}
V5_CONVERGENCE=${V5_CONVERGENCE:-$ROOT/aggregate_v5_reaudit/ej_convergence.json}
mkdir -p "$STEADY_ROOT" "$AGGREGATE_ROOT"

# Continue only from a numerically complete, sentinel-safe V5 result whose
# remaining failures are the declared energy and angular resolution gates.
"$PYTHON_BIN" - "$V5_CONVERGENCE" <<'PY'
import json,sys
d=json.load(open(sys.argv[1]))
assert d['schema']=='finite-angle-ej-convergence-v3'
assert d['status']=='INCOMPLETE_OR_FAILED'
assert d['fractional_tolerance']==0.1
assert set(d['required_dimensions'])=={'energy','angular','boundary'}
path=d['production_resolution_path']
assert path['previous_energy_bins']==32
assert path['energy_ladder_angular_bins']==65
assert path['final_energy_bins']==48
assert path['previous_angular_bins']==65
assert path['final_angular_bins']==129
assert d['dimension_gates']['boundary']['pass']
assert not d['dimension_gates']['energy']['pass']
assert not d['dimension_gates']['angular']['pass']
for group in d['groups']:
 assert group['numerical_solver_gate']
 assert len(group['seed_results']) >= 2
 assert all(
  record['seed_gate']
  for model in group['models'].values()
  for record in model.values()
 )
print('V5_RESOLUTION_ONLY_FAILURE_CONFIRMED')
PY

BASE_GROUPS=(
  local_circular_e8_default
  local_circular_e12_default
  local_circular_e16_default
  local_circular_e16_fine
  orbit_thermal_e_lt_0p99_e16_fine
  orbit_e0p9_e16_fine
  local_circular_e24_fine
  local_circular_e32_fine
  local_circular_e32_ultra
  orbit_thermal_e_lt_0p99_e32_ultra
  orbit_e0p9_e32_ultra
  local_circular_e48_ultra
  local_circular_e48_hyper
  orbit_thermal_e_lt_0p99_e48_hyper
  orbit_e0p9_e48_hyper
)
NEW_GROUPS=(
  local_circular_e64_hyper
  local_circular_e64_j257
  orbit_thermal_e_lt_0p99_e64_j257
  orbit_e0p9_e64_j257
)
RESULTS=()
for GROUP in "${BASE_GROUPS[@]}"; do
  OUT="$BASE_STEADY_ROOT/$GROUP"
  for ITEM in s314159 s271828 combined; do
    test -s "$OUT/steady_${ITEM}.json"
    test -s "$OUT/steady_${ITEM}.csv"
    grep -q '"status": "EJ_STEADY_SMOKE_PASS"' \
      "$OUT/steady_${ITEM}.json"
    RESULTS+=("$OUT/steady_${ITEM}.json")
  done
done

for GROUP in "${NEW_GROUPS[@]}"; do
  A="$ROOT/operators/${GROUP}_s314159"
  B="$ROOT/operators/${GROUP}_s271828"
  test -s "$A/STATUS.txt"; test -s "$B/STATUS.txt"
  grep -q EJ_OPERATOR_DONE "$A/STATUS.txt"
  grep -q EJ_OPERATOR_DONE "$B/STATUS.txt"
  OUT="$STEADY_ROOT/$GROUP"
  mkdir -p "$OUT"
  if test -s "$OUT/operator_combined.npz" && \
     test -s "$OUT/operator_combined.json"; then
    :
  elif test ! -e "$OUT/operator_combined.npz" && \
       test ! -e "$OUT/operator_combined.json"; then
    "$PYTHON_BIN" -u "$PKG/combine_ej_operators.py" \
      "$A/operator.npz" "$B/operator.npz" \
      --out "$OUT/operator_combined.npz" \
      --out-json "$OUT/operator_combined.json"
  else
    echo "refusing a partial combined operator: $OUT" >&2
    exit 2
  fi
  for ITEM in s314159 s271828 combined; do
    if [ "$ITEM" = s314159 ]; then OP="$A/operator.npz"; fi
    if [ "$ITEM" = s271828 ]; then OP="$B/operator.npz"; fi
    if [ "$ITEM" = combined ]; then OP="$OUT/operator_combined.npz"; fi
    if test -s "$OUT/steady_${ITEM}.json" && \
       test -s "$OUT/steady_${ITEM}.csv"; then
      :
    elif test ! -e "$OUT/steady_${ITEM}.json" && \
         test ! -e "$OUT/steady_${ITEM}.csv"; then
      "$PYTHON_BIN" -u "$PKG/solve_ej_steady_state.py" \
        --operator "$OP" --out-json "$OUT/steady_${ITEM}.json" \
        --out-csv "$OUT/steady_${ITEM}.csv" \
        --residual-tolerance 1.0e-7 --maximum-evaluations 10000
    else
      echo "refusing a partial steady-state result: $OUT/steady_${ITEM}" >&2
      exit 2
    fi
    grep -q '"status": "EJ_STEADY_SMOKE_PASS"' \
      "$OUT/steady_${ITEM}.json"
    RESULTS+=("$OUT/steady_${ITEM}.json")
  done
done

PRIMARY=$(python - "$ROOT/bridge/boundary_criteria.json" <<'PY'
import json,sys
print('%.16e'%json.load(open(sys.argv[1]))['criteria']['local_circular']['r_boundary_pc'])
PY
)
"$PYTHON_BIN" -u "$PKG/assess_ej_convergence.py" "${RESULTS[@]}" \
  --out-json "$AGGREGATE_ROOT/ej_convergence.json" \
  --out-csv "$AGGREGATE_ROOT/ej_convergence.csv" \
  --fractional-tolerance 0.10 --minimum-seeds 2 --minimum-boundaries 3 \
  --primary-boundary-radius "$PRIMARY" --require-capture-energy

grep -q '"status": "PRODUCTION_CONVERGENCE_PASS"' \
  "$AGGREGATE_ROOT/ej_convergence.json"
"$PYTHON_BIN" -u "$PKG/assemble_absolute_closure.py" \
  --remap-json "$ROOT/remap/remap.json" \
  --bridge-json "$ROOT/bridge/bridge.json" \
  --normalization-json "$ROOT/normalized/df_normalization.json" \
  --convergence-json "$AGGREGATE_ROOT/ej_convergence.json" \
  --steady-json "$STEADY_ROOT/local_circular_e64_j257/steady_combined.json" \
  --out "$AGGREGATE_ROOT/absolute_closure.json"
grep -q '"status": "ABSOLUTE_CLOSURE_MEASURED"' \
  "$AGGREGATE_ROOT/absolute_closure.json"
printf 'PRODUCTION_EJ_AGGREGATE_DONE %s ROOT=%s AGGREGATE_DIR=%s\n' \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$ROOT" "$AGGREGATE_ROOT" | tee "$AGGREGATE_STATUS"
