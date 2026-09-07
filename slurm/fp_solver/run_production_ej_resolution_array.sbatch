#!/bin/bash
#SBATCH --job-name=gnc_ej_res
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --array=0-9
#SBATCH --output=slurm-%x-%A_%a.out

set -euo pipefail
export OMP_NUM_THREADS=1
PYTHON_BIN=${PYTHON_BIN:-python3}

: "${ROOT:?submit with ROOT=<completed production root>}"
PKG=${PKG:?set PKG to the fp_solver directory}
PAIRS_PER_BIN=${PAIRS_PER_BIN:-1600000}
if ! [[ "$PAIRS_PER_BIN" =~ ^[1-9][0-9]*$ ]]; then
  echo "PAIRS_PER_BIN must be a positive integer" >&2
  exit 2
fi
test -s "$ROOT/PREP_STATUS.txt"
grep -q PRODUCTION_PREP_DONE "$ROOT/PREP_STATUS.txt"

case "$SLURM_ARRAY_TASK_ID" in
  0) CRIT=local_circular; E=24; GRID=fine;  SEED=314159 ;;
  1) CRIT=local_circular; E=24; GRID=fine;  SEED=271828 ;;
  2) CRIT=local_circular; E=32; GRID=fine;  SEED=314159 ;;
  3) CRIT=local_circular; E=32; GRID=fine;  SEED=271828 ;;
  4) CRIT=local_circular; E=32; GRID=ultra; SEED=314159 ;;
  5) CRIT=local_circular; E=32; GRID=ultra; SEED=271828 ;;
  6) CRIT=orbit_thermal_e_lt_0p99; E=32; GRID=ultra; SEED=314159 ;;
  7) CRIT=orbit_thermal_e_lt_0p99; E=32; GRID=ultra; SEED=271828 ;;
  8) CRIT=orbit_e0p9; E=32; GRID=ultra; SEED=314159 ;;
  9) CRIT=orbit_e0p9; E=32; GRID=ultra; SEED=271828 ;;
  *) exit 2 ;;
esac
RB=$(python - "$ROOT/bridge/boundary_criteria.json" "$CRIT" <<'PY'
import json,sys
print('%.16e'%json.load(open(sys.argv[1]))['criteria'][sys.argv[2]]['r_boundary_pc'])
PY
)
TAG=${CRIT}_e${E}_${GRID}_s${SEED}
OUT=$ROOT/operators/$TAG
mkdir -p "$OUT"
if test -s "$OUT/STATUS.txt" && grep -q EJ_OPERATOR_DONE "$OUT/STATUS.txt"; then
  echo "refusing to overwrite completed operator: $OUT" >&2
  exit 2
fi
printf 'EJ_RESOLUTION_OPERATOR_START %s tag=%s r=%s pairs_per_bin=%s\n' \
  "$(date -Is)" "$TAG" "$RB" "$PAIRS_PER_BIN"
"$PYTHON_BIN" -u "$PKG/finite_angle_capture.py" \
  --profile "$ROOT/bridge/bridged_profile.txt" \
  --df-table "$ROOT/normalized/f_evolved_normalized.txt" \
  --normalization-json "$ROOT/normalized/df_normalization.json" \
  --out-json "$OUT/run.json" --out-csv "$OUT/run.csv" \
  --out-ej-operator "$OUT/operator.npz" \
  --mbh 4.0e6 --sigma0-over-m 100.0 --w-kms 80.0 \
  --r-outer "$RB" --radial-bins 16 --pairs-per-bin "$PAIRS_PER_BIN" \
  --cap-importance-fraction 0.5 --direction-edge-fraction 0.3 \
  --direction-edge-power 4 --survival-anomaly-order 16 \
  --survival-collision-order 48 --ej-energy-bins "$E" \
  --ej-angular-grid "$GRID" --seed "$SEED"
"$PYTHON_BIN" - "$OUT/operator.npz" "$OUT/run.json" <<'PY'
import json,sys,numpy as np
with np.load(sys.argv[1]) as a:
 for key in ('rate_direct_binding_msun_kms2_per_myr','rate_immediate_binding_msun_kms2_per_myr'):
  assert key in a.files and np.all(a[key]>=0)
x=json.load(open(sys.argv[2]))['ej_jump_operator']
assert x['direct_capture_binding_roundtrip_relative_error'] <= 1e-10
assert x['immediate_capture_binding_roundtrip_relative_error'] <= 1e-10
print('MASS_AND_CAPTURE_ENERGY_ROUNDTRIP_PASS')
PY
printf 'EJ_OPERATOR_DONE %s tag=%s\n' "$(date -Is)" "$TAG" | tee "$OUT/STATUS.txt"
