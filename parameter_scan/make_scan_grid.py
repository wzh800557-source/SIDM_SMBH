#!/usr/bin/env python3
"""Create the controlled SIDM--SMBH parameter grid.

The source profile is rescaled homologously at fixed density and fixed
dimensionless collapse stage.  If the represented halo mass changes by a
factor f, radii and one-dimensional dispersions change by f^(1/3), while the
density profile is unchanged.  This preserves hydrostatic similarity and
changes every enclosed mass by f.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path


REFERENCE_PROFILE_MASS_MSUN = 1.5283226426126614e10
REFERENCE_BLACK_HOLE_MASS_MSUN = 4.0e6


def factor_tag(value: float) -> str:
    return {0.1: "0p1", 1.0: "1", 10.0: "10"}[float(value)]


def main() -> int:
    out = Path(__file__).with_name("scan_grid.json")
    cases = []
    index = 0
    for halo_factor, black_hole_factor, sigma_over_m in itertools.product(
        (0.1, 1.0, 10.0), (0.1, 1.0, 10.0), (10.0, 100.0, 1000.0)
    ):
        mbh = REFERENCE_BLACK_HOLE_MASS_MSUN * black_hole_factor
        tag = (
            f"h{factor_tag(halo_factor)}_"
            f"bh{factor_tag(black_hole_factor)}_"
            f"s{int(sigma_over_m)}"
        )
        cases.append(
            {
                "index": index,
                "tag": tag,
                "halo_mass_factor": halo_factor,
                "represented_halo_mass_msun": (
                    REFERENCE_PROFILE_MASS_MSUN * halo_factor
                ),
                "black_hole_mass_factor": black_hole_factor,
                "black_hole_mass_msun": mbh,
                "black_hole_to_represented_halo_mass": (
                    mbh / (REFERENCE_PROFILE_MASS_MSUN * halo_factor)
                ),
                "sigma0_over_m_cm2_g": sigma_over_m,
                "yukawa_w_kms": 80.0,
            }
        )
        index += 1
    payload = {
        "schema": "sidm-smbh-fixed-collapse-stage-grid-v1",
        "purpose": (
            "Measure interface topology and finite-angle FP currents while "
            "varying halo mass, black-hole mass, and cross-section amplitude."
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
        "grid_shape": [3, 3, 3],
        "case_count": len(cases),
        "cases": cases,
    }
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
