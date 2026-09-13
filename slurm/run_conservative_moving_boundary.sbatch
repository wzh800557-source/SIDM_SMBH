#!/bin/bash
#SBATCH --job-name=sidm_move_bnd
#SBATCH --partition=mit_preemptable
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
set -euo pipefail

: "${ENERGY_GATE:?accepted energy-current gate required}"
: "${BOUNDARY_TRAJECTORY:?live boundary trajectory required}"
: "${MOVING_OUT:?new output directory required}"
: "${PKG:?repository fp_solver directory required}"
PYTHON_BIN=${PYTHON_BIN:-python3}
set +e
"$PYTHON_BIN" -u "$PKG/conservative_moving_boundary.py" \
  --energy-gate "$ENERGY_GATE" --trajectory "$BOUNDARY_TRAJECTORY" \
  --outdir "$MOVING_OUT" \
  --fluid-mass="${FLUID_MASS_MSUN:?required}" \
  --kinetic-mass="${KINETIC_MASS_MSUN:?required}" \
  --black-hole-mass="${BLACK_HOLE_MASS_MSUN:?required}" \
  --fluid-energy="${FLUID_ENERGY_MSUN_KMS2:?required}" \
  --kinetic-energy="${KINETIC_ENERGY_MSUN_KMS2:?required}"
RC=$?
set -e
if [[ $RC == 4 ]]; then
  # Correct, expected result while the energy gate remains closed.
  exit 0
fi
exit "$RC"
