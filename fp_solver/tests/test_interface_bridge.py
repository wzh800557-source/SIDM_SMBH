#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from interface_bridge import validate_interface_bridge  # noqa: E402


def main() -> int:
    bridge = {
        "schema": "absolute-closure-hydrostatic-bridge-v3",
        "status": "OK",
        "input_profile_sha256": "abc123",
        "mbh_msun": 4.0e6,
        "sigma_over_m_cm2_g": 100.0,
        "yukawa_w_kms": 80.0,
        "fluid_subgrid_bridge_gate": "PASS",
        "mass_inside_r_in_fraction_of_first_shell": 0.006,
        "maximum_allowed_subgrid_mass_fraction": 0.02,
        "volume_inside_r_in_fraction_of_first_shell": 8.0e-4,
        "target_first_shell_mass_msun": 2.0e5,
        "reconstructed_first_shell_mass_msun": 2.0e5 * (1.0 + 1.0e-10),
        "pressure_edge_inside": 3000.0,
        "pressure_edge_outside": 3000.0,
        "r_in_pc": 2.3,
        "mass_inside_r_in_msun": 1200.0,
    }
    result = validate_interface_bridge(
        bridge,
        mbh_msun=4.0e6,
        sigma_over_m_cm2_g=100.0,
        w_kms=80.0,
        fluid_profile_sha256="abc123",
        fluid_inner_shell_mass_msun=2.0e5,
        closure_radius_pc=2.3 * (1.0 + 2.0e-4),
    )
    assert result["status"] == "PASS"
    bad = dict(bridge)
    bad["mass_inside_r_in_fraction_of_first_shell"] = 0.03
    try:
        validate_interface_bridge(
            bad,
            mbh_msun=4.0e6,
            sigma_over_m_cm2_g=100.0,
            w_kms=80.0,
            fluid_profile_sha256="abc123",
            fluid_inner_shell_mass_msun=2.0e5,
            closure_radius_pc=2.3,
        )
    except RuntimeError:
        pass
    else:
        raise AssertionError("oversized GNC-owned mass was not rejected")
    try:
        validate_interface_bridge(
            bridge,
            mbh_msun=4.0e6,
            sigma_over_m_cm2_g=100.0,
            w_kms=80.0,
            fluid_profile_sha256="different-profile",
            fluid_inner_shell_mass_msun=2.0e5,
            closure_radius_pc=2.3,
        )
    except RuntimeError:
        pass
    else:
        raise AssertionError("mismatched fluid and bridge profiles were not rejected")
    print("PASS: hydrostatic interface bridge gate")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
