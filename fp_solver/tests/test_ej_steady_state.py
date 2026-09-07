#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from solve_ej_steady_state import population_rate, solve_model  # noqa: E402


def main() -> int:
    # State 1 is a fixed reservoir. It injects state 0 at 2 Msun/Myr.
    # Capture removes state 0 at 4 D_0 Msun/Myr, so the exact steady
    # occupation is D_0=1/2 and the steady capture current is 2 Msun/Myr.
    arrays = {
        "source1": np.array([1, 0], dtype=np.int64),
        "source2": np.array([1, 1], dtype=np.int64),
        "pre": np.array([1, 0], dtype=np.int64),
        "post": np.array([0, -1], dtype=np.int64),
        "rate_direct_msun_per_myr": np.array([2.0, 4.0]),
        "rate_immediate_msun_per_myr": np.array([2.0, 4.0]),
        "rate_direct_binding_msun_kms2_per_myr": np.array([0.0, 40.0]),
        "rate_immediate_binding_msun_kms2_per_myr": np.array([0.0, 40.0]),
        "sample_count": np.array([100, 100], dtype=np.int64),
        "x_edges": np.array([0.1, 1.0]),
        "j_over_jlc_edges": np.array([1.0, 2.0]),
    }
    metadata = {
        "schema": "finite-angle-ej-jump-operator-v1",
        "state_count": 2,
    }
    result, occupation = solve_model(arrays, metadata, "direct")
    assert result["status"] == "EJ_STEADY_SMOKE_PASS"
    assert abs(occupation[0] - 0.5) < 1.0e-11
    assert occupation[1] == 1.0
    assert abs(result["steady_capture_mdot_msun_per_myr"] - 2.0) < 1.0e-11
    assert abs(
        result["steady_capture_binding_energy_current_msun_kms2_per_myr"]
        - 20.0
    ) < 1.0e-10
    assert abs(result["steady_capture_mean_binding_energy_kms2"] - 10.0) < 1.0e-11
    assert result["mass_identity_relative"] < 1.0e-12
    dot, flux = population_rate(
        occupation,
        arrays["source1"],
        arrays["source2"],
        arrays["pre"],
        arrays["post"],
        arrays["rate_direct_msun_per_myr"],
    )
    assert abs(dot[0]) < 1.0e-11
    assert abs(float(np.sum(dot)) + float(np.sum(flux[arrays["post"] == -1]))) < 1.0e-11

    # A resolved high-energy cell can be depleted far below exp(-15).  The
    # solver bound must not create a false residual failure in that regime.
    depleted = dict(arrays)
    depleted["rate_direct_msun_per_myr"] = np.array([2.0, 4.0e10])
    depleted["rate_immediate_msun_per_myr"] = np.array([2.0, 4.0e10])
    depleted["rate_direct_binding_msun_kms2_per_myr"] = np.array([0.0, 4.0e11])
    depleted["rate_immediate_binding_msun_kms2_per_myr"] = np.array([0.0, 4.0e11])
    depleted_result, depleted_occupation = solve_model(
        depleted, metadata, "direct"
    )
    assert depleted_result["status"] == "EJ_STEADY_SMOKE_PASS"
    assert abs(depleted_occupation[0] / 5.0e-11 - 1.0) < 1.0e-6
    assert depleted_result["states_at_lower_bound"] == 0

    # A drain-only cell has the physical steady occupation D=0.  The
    # non-negative solver must admit that boundary solution while retaining a
    # positive current through the separately supplied cell.
    empty_cell = {
        "source1": np.array([1, 0, 2], dtype=np.int64),
        "source2": np.array([1, 1, 3], dtype=np.int64),
        "pre": np.array([1, 0, 2], dtype=np.int64),
        "post": np.array([0, -1, -1], dtype=np.int64),
        "rate_direct_msun_per_myr": np.array([2.0, 4.0, 3.0]),
        "rate_immediate_msun_per_myr": np.array([2.0, 4.0, 3.0]),
        "rate_direct_binding_msun_kms2_per_myr": np.array([0.0, 40.0, 30.0]),
        "rate_immediate_binding_msun_kms2_per_myr": np.array([0.0, 40.0, 30.0]),
        "sample_count": np.array([100, 100, 100], dtype=np.int64),
        "x_edges": np.array([0.1, 1.0]),
        "j_over_jlc_edges": np.array([1.0, 2.0, 3.0]),
    }
    empty_metadata = {
        "schema": "finite-angle-ej-jump-operator-v1",
        "state_count": 4,
    }
    empty_result, empty_occupation = solve_model(
        empty_cell, empty_metadata, "direct"
    )
    assert empty_result["status"] == "EJ_STEADY_SMOKE_PASS"
    assert abs(empty_occupation[0] - 0.5) < 1.0e-11
    assert empty_occupation[2] <= 1.0e-29
    assert empty_result["states_at_lower_bound"] == 1
    print("PASS: non-local EJ steady-state balance")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
