#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from compare_finite_angle_boundaries import compare_runs  # noqa: E402


def make_run(radius: float, seed: int, current: float) -> dict:
    return {
        "schema": "finite-angle-loss-cone-capture-v6",
        "status": "ISOTROPIC_INJECTION_CEILING_CONVERGED",
        "profile_sha256": "profile",
        "bridge_json_sha256": "bridge",
        "df_table_sha256": "df",
        "df_table_schema": "gnc-normalized-df-v6",
        "mbh_msun": 4.0e6,
        "rh_pc": 17.5,
        "sigma0_over_m_cm2_g": 100.0,
        "w_kms": 80.0,
        "r_inner_pc": 3.0e-6,
        "cap_importance_fraction": 0.9,
        "direction_edge_fraction": 0.3,
        "direction_edge_power": 4.0,
        "survival_anomaly_order": 16,
        "survival_collision_quadrature_order": 48,
        "r_outer_pc": radius,
        "seed": seed,
        "candidate_mdot_msun_per_myr": current,
        "candidate_mdot_standard_error_msun_per_myr": 0.1,
        "direct_no_rescatter_mdot_msun_per_myr": 0.93 * current,
        "direct_no_rescatter_mdot_standard_error_msun_per_myr": 0.1,
        "candidate_binding_energy_current_msun_kms2_per_myr": 4.0e7 * current,
        "candidate_binding_energy_current_standard_error_msun_kms2_per_myr": 4.0e6,
        "rows": [{
            "angular_importance_model": (
                "physical kernel plus projected-disc proposal, with a focused "
                "near-tangent fallback"
            )
        }],
    }


def main() -> int:
    runs = [
        make_run(2.3, 1, 1.80),
        make_run(2.5, 2, 1.84),
        make_run(2.7, 3, 1.79),
    ]
    result = compare_runs(runs)
    assert result["status"] == "FINITE_ANGLE_BOUNDARY_INVARIANCE_PASS"
    bad = [dict(run) for run in runs]
    bad[-1]["candidate_mdot_msun_per_myr"] = 3.0
    assert compare_runs(bad)["status"] == "FAIL"
    print("PASS: finite-angle boundary invariance gates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
