#!/bin/bash
#SBATCH --job-name=gnc_e_gate
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
set -euo pipefail

: "${RUN_ROOT:?GNC energy-current ensemble root required}"
: "${ENERGY_GATE:?output gate JSON required}"
: "${PKG:?repository fp_solver directory required}"
PYTHON_BIN=${PYTHON_BIN:-python3}
set +e
"$PYTHON_BIN" -u "$PKG/gnc_energy_current_gate.py" \
  --run-root "$RUN_ROOT" --out "$ENERGY_GATE" \
  --minimum-raw-captures 50 --minimum-effective-captures 50 \
  --plateau-tolerance 0.25 --inventory-tolerance 0.25 \
  --seed-tolerance 0.25 --grid-tolerance 0.10 \
  > "$RUN_ROOT/energy_current_gate.log" 2>&1
RC=$?
set -e
if [[ $RC != 0 && $RC != 3 ]]; then exit "$RC"; fi
"$PYTHON_BIN" - "$ENERGY_GATE" "$RUN_ROOT/GATE_STATUS.txt" <<'PY'
import json,sys
x=json.load(open(sys.argv[1]))
open(sys.argv[2],'w').write(x['status']+'\n'+','.join(x['failed_gates'])+'\n')
PY
# A closed scientific gate is a completed audit, not a Slurm failure.
exit 0
