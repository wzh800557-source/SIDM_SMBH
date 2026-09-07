#!/bin/bash
#SBATCH --job-name=gnc_prod_prep
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%x-%j.out

set -euo pipefail
export OMP_NUM_THREADS=1
PYTHON_BIN=${PYTHON_BIN:-python3}

PKG=${PKG:?set PKG to the fp_solver directory}
PROFILE=${PROFILE:?set PROFILE to the input fluid snapshot}
ROOT=${ROOT:?set ROOT to a new output directory}
MBH=${MBH:-4.0e6}
SIGMA0=${SIGMA0:-100.0}
W_KMS=${W_KMS:-80.0}
BH_FORCE_STEP=${BH_FORCE_STEP:-0.02}
EXPECTED_PROFILE_SHA=${EXPECTED_PROFILE_SHA:-}
test ! -e "$ROOT" || { echo "refusing existing production root: $ROOT" >&2; exit 2; }
mkdir -p "$ROOT/input" "$ROOT/remap" "$ROOT/bridge" "$ROOT/normalized" "$ROOT/operators" "$ROOT/steady" "$ROOT/aggregate" "$ROOT/fluid"
printf 'PRODUCTION_PREP_START %s ROOT=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$ROOT"

for f in prepare_bh_remapped_profile.py measured_fluid_feedback.py \
         fluid_bh_scaleheight.py interface_bridge.py hydrostatic_bridge.py \
         boundary_criterion_scan.py normalize_gnc_df.py sidm_born_kernel.py; do
  test -s "$PKG/$f"
done
test -s "$PROFILE"
PROFILE_SHA=$("$PYTHON_BIN" -c 'import hashlib,sys; print(hashlib.sha256(open(sys.argv[1], "rb").read()).hexdigest())' "$PROFILE")
if [[ -n "$EXPECTED_PROFILE_SHA" ]]; then
  test "$PROFILE_SHA" = "$EXPECTED_PROFILE_SHA"
fi
cp "$PROFILE" "$ROOT/input/prof_deep.txt"

"$PYTHON_BIN" -u "$PKG/prepare_bh_remapped_profile.py" \
  --profile "$ROOT/input/prof_deep.txt" \
  --out-profile "$ROOT/remap/profile_bh.txt" \
  --out-json "$ROOT/remap/remap.json" \
  --record-dir "$ROOT/remap/record" \
  --mbh "$MBH" --sigma-over-m "$SIGMA0" --w-kms "$W_KMS" \
  --max-force-fraction "$BH_FORCE_STEP" --r-epsilon 1.0e-9

grep -q '"status": "BH_REMAP_COMPLETE"' "$ROOT/remap/remap.json"
"$PYTHON_BIN" -u "$PKG/hydrostatic_bridge.py" "$ROOT/remap/profile_bh.txt" \
  --out-profile "$ROOT/bridge/bridged_profile.txt" \
  --out-json "$ROOT/bridge/bridge.json" \
  --mbh "$MBH" --sigma-over-m "$SIGMA0" \
  --mode adiabatic-mass-matched --kernel yukawa-tchannel --w-kms "$W_KMS"
grep -q '"status": "OK"' "$ROOT/bridge/bridge.json"
grep -q '"fluid_subgrid_bridge_gate": "PASS"' "$ROOT/bridge/bridge.json"

"$PYTHON_BIN" -u "$PKG/boundary_criterion_scan.py" \
  "$ROOT/bridge/bridged_profile.txt" \
  --bridge-json "$ROOT/bridge/bridge.json" \
  --out "$ROOT/bridge/boundary_criteria.json" \
  --mbh "$MBH" --sigma-over-m "$SIGMA0" \
  --kernel yukawa-tchannel --w-kms "$W_KMS"
"$PYTHON_BIN" - "$ROOT/bridge/boundary_criteria.json" <<'PY'
import json,sys
x=json.load(open(sys.argv[1]))
assert x['fp_loss_cone_angular_kernel_gate']['status'] == 'NONLOCAL_JUMP_OPERATOR_REQUIRED'
for key in ('local_circular','orbit_thermal_e_lt_0p99','orbit_e0p9'):
    assert x['criteria'][key]['r_boundary_pc'] > 0
print('NONLOCAL_OPERATOR_GATE_CONFIRMED')
PY
RB=$("$PYTHON_BIN" - "$ROOT/bridge/boundary_criteria.json" <<'PY'
import json,sys
print('%.16e'%json.load(open(sys.argv[1]))['criteria']['local_circular']['r_boundary_pc'])
PY
)

"$PYTHON_BIN" -u "$PKG/normalize_gnc_df.py" \
  "$ROOT/bridge/bridged_profile.txt" \
  --bridge-json "$ROOT/bridge/bridge.json" \
  --outdir "$ROOT/normalized" --mbh "$MBH" \
  --r-boundary "$RB" --normalization-radius "$RB" \
  --proposal-boundary-radius "$RB" --grid-bins 104 --emax 2.0e7 \
  --df-fit-points 512 --ranks 4 --samples-per-rank 5000
grep -q '"status": "PASS"' "$ROOT/normalized/df_normalization.json"

"$PYTHON_BIN" - "$ROOT" "$PKG" <<'PY'
import hashlib,json,sys
from pathlib import Path
root,pkg=map(Path,sys.argv[1:])
def sha(p):
 h=hashlib.sha256(); h.update(p.read_bytes()); return h.hexdigest()
remap=json.loads((root/'remap/remap.json').read_text())
bridge=json.loads((root/'bridge/bridge.json').read_text())
norm=json.loads((root/'normalized/df_normalization.json').read_text())
out={
 'schema':'gnc-production-preparation-v1',
 'status':'PRODUCTION_PREPARATION_PASS',
 'root':root.name,
 'package':pkg.name,
 'source_profile_sha256':sha(root/'input/prof_deep.txt'),
 'post_bh_profile_sha256':sha(root/'remap/profile_bh.txt'),
 'bridged_profile_sha256':sha(root/'bridge/bridged_profile.txt'),
 'df_table_sha256':sha(root/'normalized/f_evolved_normalized.txt'),
 'r_in_pc':bridge['r_in_pc'], 'rh_pc':bridge['rh_pc'],
 'M_bh_msun':bridge['mbh_msun'],
 'sigma0_over_m_cm2_g':bridge['sigma_over_m_cm2_g'],
 'w_kms':bridge['yukawa_w_kms'],
 'gates':{
  'remap':remap['status'], 'bridge':bridge['status'],
  'subgrid_bridge':bridge['fluid_subgrid_bridge_gate'],
  'df_normalization':norm['status'],
 }
}
(root/'prep_manifest.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))
PY
printf 'PRODUCTION_PREP_DONE %s ROOT=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$ROOT" | tee "$ROOT/PREP_STATUS.txt"
