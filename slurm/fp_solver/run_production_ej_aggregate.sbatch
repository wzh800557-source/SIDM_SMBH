#!/bin/bash
#SBATCH --job-name=gnc_prod_agg
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail
export OMP_NUM_THREADS=1
PYTHON_BIN=${PYTHON_BIN:-python3}
: "${ROOT:?submit with ROOT=<production root>}"
PKG=${PKG:?set PKG to the fp_solver directory}
mkdir -p "$ROOT/steady" "$ROOT/aggregate"

RUN_GROUPS=(
  local_circular_e8_default
  local_circular_e12_default
  local_circular_e16_default
  local_circular_e16_fine
  orbit_thermal_e_lt_0p99_e16_fine
  orbit_e0p9_e16_fine
)
RESULTS=()
for GROUP in "${RUN_GROUPS[@]}"; do
  A="$ROOT/operators/${GROUP}_s314159"
  B="$ROOT/operators/${GROUP}_s271828"
  test -s "$A/STATUS.txt"; test -s "$B/STATUS.txt"
  grep -q EJ_OPERATOR_DONE "$A/STATUS.txt"
  grep -q EJ_OPERATOR_DONE "$B/STATUS.txt"
  OUT="$ROOT/steady/$GROUP"
  mkdir -p "$OUT"
  "$PYTHON_BIN" -u "$PKG/combine_ej_operators.py" \
    "$A/operator.npz" "$B/operator.npz" \
    --out "$OUT/operator_combined.npz" --out-json "$OUT/operator_combined.json"
  for ITEM in s314159 s271828 combined; do
    if [ "$ITEM" = s314159 ]; then OP="$A/operator.npz"; fi
    if [ "$ITEM" = s271828 ]; then OP="$B/operator.npz"; fi
    if [ "$ITEM" = combined ]; then OP="$OUT/operator_combined.npz"; fi
    "$PYTHON_BIN" -u "$PKG/solve_ej_steady_state.py" \
      --operator "$OP" --out-json "$OUT/steady_${ITEM}.json" \
      --out-csv "$OUT/steady_${ITEM}.csv" \
      --residual-tolerance 1.0e-7 --maximum-evaluations 10000
    RESULTS+=("$OUT/steady_${ITEM}.json")
  done
done

PRIMARY=$("$PYTHON_BIN" - "$ROOT/bridge/boundary_criteria.json" <<'PY'
import json,sys
print('%.16e'%json.load(open(sys.argv[1]))['criteria']['local_circular']['r_boundary_pc'])
PY
)
"$PYTHON_BIN" -u "$PKG/assess_ej_convergence.py" "${RESULTS[@]}" \
  --out-json "$ROOT/aggregate/ej_convergence.json" \
  --out-csv "$ROOT/aggregate/ej_convergence.csv" \
  --fractional-tolerance 0.10 --minimum-seeds 2 --minimum-boundaries 3 \
  --primary-boundary-radius "$PRIMARY" \
  --require-capture-energy

grep -q '"status": "PRODUCTION_CONVERGENCE_PASS"' \
  "$ROOT/aggregate/ej_convergence.json"
"$PYTHON_BIN" -u "$PKG/assemble_absolute_closure.py" \
  --remap-json "$ROOT/remap/remap.json" \
  --bridge-json "$ROOT/bridge/bridge.json" \
  --normalization-json "$ROOT/normalized/df_normalization.json" \
  --convergence-json "$ROOT/aggregate/ej_convergence.json" \
  --steady-json "$ROOT/steady/local_circular_e16_fine/steady_combined.json" \
  --out "$ROOT/aggregate/absolute_closure.json"
grep -q '"status": "ABSOLUTE_CLOSURE_MEASURED"' \
  "$ROOT/aggregate/absolute_closure.json"
printf 'PRODUCTION_EJ_AGGREGATE_DONE %s ROOT=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$ROOT" | tee "$ROOT/AGGREGATE_STATUS.txt"
