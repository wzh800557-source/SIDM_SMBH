#!/usr/bin/env python3
"""Lightweight local checks for the cluster scan package."""

from __future__ import annotations

import importlib.util
import json
import math
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main() -> int:
    grid = json.loads((HERE / "scan_grid.json").read_text())
    assert grid["case_count"] == 27
    assert len(grid["cases"]) == 27
    assert len({case["tag"] for case in grid["cases"]}) == 27
    assert {case["halo_mass_factor"] for case in grid["cases"]} == {0.1, 1.0, 10.0}
    assert {case["black_hole_mass_msun"] for case in grid["cases"]} == {
        4.0e5,
        4.0e6,
        4.0e7,
    }
    assert {case["sigma0_over_m_cm2_g"] for case in grid["cases"]} == {
        10.0,
        100.0,
        1000.0,
    }

    scale = load_module("scale_physical_profile", HERE / "scale_physical_profile.py")
    r = np.geomspace(1.0, 100.0, 12)
    rho = 3.0 / (1.0 + r) ** 2
    source_mass = scale.enclosed_mass(r, rho)
    factor = 10.0
    length = factor ** (1.0 / 3.0)
    scaled_mass = scale.enclosed_mass(r * length, rho)
    assert math.isclose(scaled_mass / source_mass, factor, rel_tol=2.0e-15)

    summary = load_module("summarize_scan_case", HERE / "summarize_scan_case.py")
    assert summary.symmetric_difference(10.0, 11.0) == 1.0 / 10.5
    assert summary.symmetric_difference(0.0, 0.0) == 0.0

    operator = load_module("run_scan_operator", HERE / "run_scan_operator.py")
    diagnostic = {
        "status": "DIAGNOSTIC_ONLY",
        "rho_df_profile_max_abs_relative_error": 1.0e-5,
        "candidate_relative_standard_error": 0.01,
        "accepted_capture_relative_standard_error": 0.01,
        "candidate_mdot_msun_per_myr": 100.0,
        "direct_no_rescatter_mdot_msun_per_myr": 10.0,
        "unresolved_kepler_domain_fraction_of_candidate": 0.03,
        "nonlocal_jump_balance": {
            "capture_roundtrip_relative_error": 1.0e-15,
            "mass_balance_max_relative_error": 1.0e-15,
            "mass_balance_max_absolute_error_msun_per_myr": 1.0e-15,
        },
        "ej_jump_operator": {
            "direct_capture_roundtrip_relative_error": 1.0e-15,
            "immediate_capture_roundtrip_relative_error": 1.0e-15,
            "direct_capture_binding_roundtrip_relative_error": 1.0e-15,
            "immediate_capture_binding_roundtrip_relative_error": 1.0e-15,
        },
    }
    exclusion = operator.estimator_domain_exclusion(diagnostic)
    assert exclusion is not None
    assert exclusion["reason"] == "CAPTURED_ORBITS_EXCEED_KEPLER_DOMAIN"

    # A redundant J-only projection may miss a fixed relative roundoff gate
    # after many signed additions.  It is acceptable only when the absolute
    # discrepancy is negligible compared with the measured Monte-Carlo error;
    # the authoritative E-J round trips are checked independently.
    diagnostic["unresolved_kepler_domain_fraction_of_candidate"] = 0.001
    diagnostic["nonlocal_jump_balance"].update({
        "mass_balance_max_relative_error": 2.3262938379828222e-10,
        "mass_balance_max_absolute_error_msun_per_myr": 2.0e-6,
    })
    accepted, mode = operator.estimator_acceptance(diagnostic)
    assert accepted
    assert mode == "AUXILIARY_J_BALANCE_BELOW_MC_FLOOR"
    diagnostic["nonlocal_jump_balance"][
        "mass_balance_max_absolute_error_msun_per_myr"
    ] = 2.0e-4
    accepted, mode = operator.estimator_acceptance(diagnostic)
    assert not accepted
    assert "j_mass_balance_below_mc_floor" in mode

    aggregate = load_module("aggregate_scan", HERE / "aggregate_scan.py")
    rows = []
    for case in grid["cases"]:
        mbh = case["black_hole_mass_msun"] / 4.0e6
        halo = case["halo_mass_factor"]
        sigma = case["sigma0_over_m_cm2_g"] / 100.0
        rows.append(
            {
                **case,
                "flux_status": "SCAN_FLUX_CONVERGED",
                "all_halo_flux_status": "ALL_HALO_FLUX_COMPLETE",
                "prep_status": "READY_FOR_FP",
                "synthetic": 7.0 * mbh**2.0 * halo**0.5 * sigma**-0.25,
            }
        )
    fit = aggregate.fit_power_law(rows, "synthetic", "fp")
    assert fit is not None
    assert math.isclose(fit["normalization_at_fiducial"], 7.0, rel_tol=2.0e-14)
    assert math.isclose(fit["black_hole_mass_exponent"], 2.0, abs_tol=2.0e-14)
    assert math.isclose(fit["halo_mass_exponent"], 0.5, abs_tol=2.0e-14)
    assert math.isclose(fit["cross_section_exponent"], -0.25, abs_tol=2.0e-14)

    with tempfile.TemporaryDirectory() as directory:
        assert Path(directory).is_dir()
    print("SCAN_PACKAGE_TESTS_PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
