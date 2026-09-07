#!/bin/bash
#SBATCH --job-name=gnc_mass_reaudit
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=192G
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail
export OMP_NUM_THREADS=1
PYTHON_BIN=${PYTHON_BIN:-python3}
: "${ROOT:?submit with ROOT=<completed V6 production root>}"
PKG=${PKG:?set PKG to the fp_solver directory}
V6=$ROOT/aggregate_v6/ej_convergence.json
OUT=${OUT:-$ROOT/aggregate_v6_mass_reaudit}
STATUS=${STATUS:-$ROOT/MASS_CLOSURE_V6_STATUS.txt}
test ! -e "$OUT"
test ! -e "$STATUS"
mkdir -p "$OUT"

# The existing two-current result must have failed only because the selected
# capture-binding-energy energy-grid comparisons exceeded ten per cent. This
# guard prevents a mass closure from bypassing any mass, seed, solver, angular,
# or reporting-surface failure.
"$PYTHON_BIN" - "$V6" <<'PY'
import json,sys
d=json.load(open(sys.argv[1]))
assert d['schema']=='finite-angle-ej-convergence-v3'
assert d['status']=='INCOMPLETE_OR_FAILED'
assert d['fractional_tolerance']==0.1
assert set(d['required_metrics'])=={'mass','capture_binding_energy'}
assert all(group['numerical_solver_gate'] for group in d['groups'])
assert all(
 record['seed_gate']
 for group in d['groups']
 for model in group['models'].values()
 for record in model.values()
)
failed=[]
for dimension,items in d['selected_comparisons'].items():
 for item in items:
  if not item['gate']:
   failed.append((dimension,item['model'],item['metric']))
assert sorted(failed)==[
 ('energy','direct','capture_binding_energy'),
 ('energy','immediate','capture_binding_energy'),
]
assert all(
 item['gate']
 for items in d['selected_comparisons'].values()
 for item in items
 if item['metric']=='mass'
)
print('V6_MASS_ONLY_CONVERGENCE_CANDIDATE_CONFIRMED')
PY

BASE_STEADY_ROOT=$ROOT/steady_v5
STEADY_ROOT=$ROOT/steady_v6
GROUPS_BASE=(
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
GROUPS_V6=(
  local_circular_e64_hyper
  local_circular_e64_j257
  orbit_thermal_e_lt_0p99_e64_j257
  orbit_e0p9_e64_j257
)
RESULTS=()
for GROUP_NAME in "${GROUPS_BASE[@]}"; do
  for ITEM in s314159 s271828 combined; do
    FILE=$BASE_STEADY_ROOT/$GROUP_NAME/steady_${ITEM}.json
    test -s "$FILE"
    RESULTS+=("$FILE")
  done
done
for GROUP_NAME in "${GROUPS_V6[@]}"; do
  for ITEM in s314159 s271828 combined; do
    FILE=$STEADY_ROOT/$GROUP_NAME/steady_${ITEM}.json
    test -s "$FILE"
    RESULTS+=("$FILE")
  done
done
PRIMARY=$("$PYTHON_BIN" - "$ROOT/bridge/boundary_criteria.json" <<'PY'
import json,sys
print('%.16e'%json.load(open(sys.argv[1]))['criteria']['local_circular']['r_boundary_pc'])
PY
)
set +e
"$PYTHON_BIN" -u "$PKG/assess_ej_convergence.py" "${RESULTS[@]}" \
  --out-json "$OUT/ej_convergence.json" \
  --out-csv "$OUT/ej_convergence.csv" \
  --fractional-tolerance 0.10 --minimum-seeds 2 --minimum-boundaries 3 \
  --primary-boundary-radius "$PRIMARY" --require-capture-energy
ASSESS_RC=$?
set -e
test "$ASSESS_RC" -eq 2
"$PYTHON_BIN" - "$OUT/ej_convergence.json" <<'PY'
import json,sys
d=json.load(open(sys.argv[1]))
assert d['status']=='INCOMPLETE_OR_FAILED'
assert d['observable_gates']['mass']['status']=='PRODUCTION_OBSERVABLE_CONVERGENCE_PASS'
assert d['observable_gates']['mass']['pass'] is True
assert d['observable_gates']['capture_binding_energy']['pass'] is False
assert d['observable_gates']['capture_binding_energy']['dimensions']['energy']['pass'] is False
print('V6_MASS_OBSERVABLE_CONVERGENCE_PASS')
PY

"$PYTHON_BIN" -u "$PKG/diagnose_ej_capture_tail.py" \
  --case e48_j129 \
    "$BASE_STEADY_ROOT/local_circular_e48_hyper/operator_combined.npz" \
    "$BASE_STEADY_ROOT/local_circular_e48_hyper/steady_combined.json" \
    "$BASE_STEADY_ROOT/local_circular_e48_hyper/steady_combined.csv" \
  --case e64_j129 \
    "$STEADY_ROOT/local_circular_e64_hyper/operator_combined.npz" \
    "$STEADY_ROOT/local_circular_e64_hyper/steady_combined.json" \
    "$STEADY_ROOT/local_circular_e64_hyper/steady_combined.csv" \
  --out "$OUT/capture_tail_diagnostic.json"
grep -q '"status": "CAPTURE_TAIL_DIAGNOSTIC_PASS"' \
  "$OUT/capture_tail_diagnostic.json"

"$PYTHON_BIN" -u "$PKG/assemble_fluid_mass_closure.py" \
  --remap-json "$ROOT/remap/remap.json" \
  --bridge-json "$ROOT/bridge/bridge.json" \
  --normalization-json "$ROOT/normalized/df_normalization.json" \
  --convergence-json "$OUT/ej_convergence.json" \
  --steady-json \
    "$STEADY_ROOT/local_circular_e64_j257/steady_combined.json" \
  --out "$OUT/fluid_mass_closure.json"
grep -q '"status": "FLUID_MASS_CLOSURE_MEASURED"' \
  "$OUT/fluid_mass_closure.json"
printf 'PRODUCTION_MASS_CLOSURE_REAUDIT_DONE %s ROOT=%s OUT=%s\n' \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$ROOT" "$OUT" | tee "$STATUS"
