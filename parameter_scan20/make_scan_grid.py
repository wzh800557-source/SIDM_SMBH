#!/usr/bin/env python3
"""Create a 20 by 20 halo-mass and black-hole-mass scan.

Each of the three cross-section slices contains 400 independently remapped
halo snapshots.  The mass axes are logarithmically spaced over the same
factor-of-ten range used by the validated 3 by 3 scan.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


REFERENCE_PROFILE_MASS_MSUN = 1.5283226426126614e10
REFERENCE_BLACK_HOLE_MASS_MSUN = 4.0e6
AXIS_POINTS = 20


def anchored_log_axis() -> np.ndarray:
    """Return 20 increasing values that retain 0.1, 1, and 10 exactly."""
    lower = np.geomspace(0.1, 1.0, 10)
    upper = np.geomspace(1.0, 10.0, 11)[1:]
    return np.concatenate((lower, upper))


HALO_FACTORS = anchored_log_axis()
BLACK_HOLE_FACTORS = anchored_log_axis()
SIGMA_OVER_M = (10.0, 100.0, 1000.0)


def main() -> int:
    cases = []
    index = 0
    for sigma_index, sigma_over_m in enumerate(SIGMA_OVER_M):
        for halo_index, halo_factor_np in enumerate(HALO_FACTORS):
            halo_factor = float(halo_factor_np)
            halo_mass = REFERENCE_PROFILE_MASS_MSUN * halo_factor
            for black_hole_index, black_hole_factor_np in enumerate(
                BLACK_HOLE_FACTORS
            ):
                black_hole_factor = float(black_hole_factor_np)
                black_hole_mass = (
                    REFERENCE_BLACK_HOLE_MASS_MSUN * black_hole_factor
                )
                tag = (
                    f"s{int(sigma_over_m):04d}_"
                    f"h{halo_index:02d}_bh{black_hole_index:02d}"
                )
                cases.append(
                    {
                        "index": index,
                        "tag": tag,
                        "sigma_grid_index": sigma_index,
                        "halo_grid_index": halo_index,
                        "black_hole_grid_index": black_hole_index,
                        "halo_mass_factor": halo_factor,
                        "represented_halo_mass_msun": halo_mass,
                        "black_hole_mass_factor": black_hole_factor,
                        "black_hole_mass_msun": black_hole_mass,
                        "black_hole_to_represented_halo_mass": (
                            black_hole_mass / halo_mass
                        ),
                        "sigma0_over_m_cm2_g": sigma_over_m,
                        "yukawa_w_kms": 80.0,
                    }
                )
                index += 1

    payload = {
        "schema": "sidm-smbh-fixed-collapse-stage-grid-20x20-v1",
        "purpose": (
            "Resolve the collision-defined FP-fluid interface and the "
            "black-hole-aware fluid luminosity on a dense parameter grid."
        ),
        "reference_profile_mass_msun": REFERENCE_PROFILE_MASS_MSUN,
        "reference_black_hole_mass_msun": REFERENCE_BLACK_HOLE_MASS_MSUN,
        "halo_rescaling": {
            "radius_factor": "halo_mass_factor^(1/3)",
            "density_factor": 1.0,
            "one_dimensional_dispersion_factor": "halo_mass_factor^(1/3)",
            "interpretation": (
                "homologous snapshots at the same dimensionless "
                "gravothermal-collapse stage"
            ),
        },
        "axis_spacing": "piecewise logarithmic, anchored at factors 0.1, 1, and 10",
        "grid_shape": [3, AXIS_POINTS, AXIS_POINTS],
        "axis_order": [
            "sigma0_over_m_cm2_g",
            "halo_mass_factor",
            "black_hole_mass_factor",
        ],
        "sigma0_over_m_cm2_g": list(SIGMA_OVER_M),
        "halo_mass_factors": [float(value) for value in HALO_FACTORS],
        "black_hole_mass_factors": [
            float(value) for value in BLACK_HOLE_FACTORS
        ],
        "case_count": len(cases),
        "cases": cases,
    }
    output = Path(__file__).with_name("scan_grid.json")
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(output)
    print(f"CASES={len(cases)} GRID=3x{AXIS_POINTS}x{AXIS_POINTS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
