#!/bin/bash
#SBATCH --job-name=gnc_mass_fluid
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail
export OMP_NUM_THREADS=1
PYTHON_BIN=${PYTHON_BIN:-python3}
: "${ROOT:?submit with ROOT=<production root>}"
PKG=${PKG:?set PKG to the fp_solver directory}
MASS_DIR=${MASS_DIR:-$ROOT/aggregate_v6_mass_reaudit}
MASS_STATUS=${MASS_STATUS:-$ROOT/MASS_CLOSURE_V6_STATUS.txt}
FLUID_DIR=${FLUID_DIR:-$ROOT/fluid_mass_v6}
FLUID_STATUS=${FLUID_STATUS:-$ROOT/FLUID_MASS_V6_STATUS.txt}
test -s "$MASS_STATUS"
grep -q PRODUCTION_MASS_CLOSURE_REAUDIT_DONE "$MASS_STATUS"
test -s "$MASS_DIR/fluid_mass_closure.json"
grep -q '"status": "FLUID_MASS_CLOSURE_MEASURED"' \
  "$MASS_DIR/fluid_mass_closure.json"
mkdir -p "$FLUID_DIR"
for path in records fluid_feedback_trajectories.csv fluid_feedback_summary.json \
            fluid_response_comparison.csv fluid_response_analysis.json; do
  test ! -e "$FLUID_DIR/$path" || {
    echo "refusing to overwrite existing fluid output: $FLUID_DIR/$path" >&2
    exit 2
  }
done
"$PYTHON_BIN" -u "$PKG/measured_fluid_feedback.py" \
  --profile "$ROOT/remap/profile_bh.txt" --profile-state post-bh \
  --diagnostics "$MASS_DIR/fluid_mass_closure.json" \
  --bridge-json "$ROOT/bridge/bridge.json" \
  --outdir "$FLUID_DIR" --mbh 4.0e6 --sigma-over-m 100.0 --w-units 80.0 \
  --branches control,sink --tau-end 0.06 --output-dtau 0.002 \
  --t-epsilon 1.0e-4 --wall-seconds 25000 \
  --max-unmodeled-fluid-shell-mass-fraction 0.02 \
  --max-unmodeled-gnc-mass-fraction 0.20
"$PYTHON_BIN" -u "$PKG/analyze_fluid_response.py" \
  --summary "$FLUID_DIR/fluid_feedback_summary.json" \
  --trajectories "$FLUID_DIR/fluid_feedback_trajectories.csv" \
  --closure-json "$MASS_DIR/fluid_mass_closure.json" \
  --profile "$ROOT/remap/profile_bh.txt" \
  --bridge-json "$ROOT/bridge/bridge.json" \
  --out-json "$FLUID_DIR/fluid_response_analysis.json" \
  --out-csv "$FLUID_DIR/fluid_response_comparison.csv"
grep -q '"schema": "black-hole-aware-fluid-response-v3"' \
  "$FLUID_DIR/fluid_response_analysis.json"
grep -q '"status": "BLACK_HOLE_AWARE_RESPONSE_COMPLETE"' \
  "$FLUID_DIR/fluid_response_analysis.json"
grep -q '"closure_scope": "mass_current_closure_only"' \
  "$FLUID_DIR/fluid_response_analysis.json"
grep -q '"capture_binding_energy_current_used_by_fluid": false' \
  "$FLUID_DIR/fluid_response_analysis.json"
printf 'PRODUCTION_FLUID_MASS_RESPONSE_DONE %s ROOT=%s FLUID_DIR=%s\n' \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$ROOT" "$FLUID_DIR" | tee "$FLUID_STATUS"
