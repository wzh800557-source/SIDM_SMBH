#!/usr/bin/env python3
"""Regression test for finite-angle convergence and weight-tail gates."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "summarize_finite_angle_runs.py"


def mock_run(seed: int, radial_bins: int, maximum_weight: float) -> dict:
    return {
        "schema": "finite-angle-loss-cone-capture-v6",
        "status": "ISOTROPIC_INJECTION_CEILING_CONVERGED",
        "profile": "profile.txt",
        "profile_sha256": "profile-sha",
        "bridge_json_sha256": "bridge-sha",
        "df_table": "df.txt",
        "df_table_sha256": "df-sha",
        "df_table_schema": "gnc-normalized-df-v6",
        "mbh_msun": 4.0e6,
        "rh_pc": 17.5,
        "sigma0_over_m_cm2_g": 100.0,
        "w_kms": 80.0,
        "r_inner_pc": 1.5e-6,
        "r_outer_pc": 2.33,
        "cap_importance_fraction": 0.5,
        "direction_edge_fraction": 0.3,
        "direction_edge_power": 4.0,
        "survival_anomaly_order": 16,
        "survival_collision_quadrature_order": 48,
        "mass_flux_scale_msun_per_myr": 6.0e4,
        "energy_flux_scale_msun_kms2_per_myr": 2.4e8,
        "boundary_sigma1d_kms": 63.0,
        "seed": seed,
        "radial_bins": radial_bins,
        "pairs_per_bin": 10000,
        "candidate_mdot_msun_per_myr": 2.0,
        "candidate_mdot_standard_error_msun_per_myr": 0.10,
        "accepted_capture_mdot_msun_per_myr": 1.86,
        "accepted_capture_mdot_standard_error_msun_per_myr": 0.10,
        "candidate_binding_energy_current_msun_kms2_per_myr": 8.0e7,
        "candidate_binding_energy_current_standard_error_msun_kms2_per_myr": 4.0e6,
        "accepted_capture_binding_energy_current_msun_kms2_per_myr": 7.5e7,
        "accepted_capture_binding_energy_current_standard_error_msun_kms2_per_myr": 4.0e6,
        "rho_df_profile_max_abs_relative_error": 1.0e-5,
        "unresolved_kepler_domain_fraction_of_candidate": 1.0e-4,
        "candidate_current_fraction_from_source_j_below_5_jlc": 0.75,
        "rows": [{
            "angular_importance_model": (
                "physical kernel plus projected-disc proposal, with a focused "
                "near-tangent fallback"
            ),
            "angular_importance_weight_range": [0.01, maximum_weight],
        }],
    }


def summarize(directory: Path, maximum_weight: float) -> dict:
    inputs = []
    for index, (seed, bins) in enumerate(((11, 24), (22, 48), (33, 48))):
        path = directory / f"run{index}.json"
        path.write_text(json.dumps(mock_run(seed, bins, maximum_weight)))
        inputs.append(path)
    output = directory / "summary.json"
    subprocess.run(
        [sys.executable, str(SCRIPT), *map(str, inputs), "--out", str(output)],
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return json.loads(output.read_text())


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)
        accepted = summarize(directory, 2.0)
        assert accepted["status"] == (
            "FINITE_ANGLE_ISOTROPIC_INJECTION_CEILING_CONVERGED"
        )
        assert accepted["gates"]["angular_importance_weight_le_2p01"]

        rejected = summarize(directory, 10.0)
        assert rejected["status"] == "FAIL"
        assert not rejected["gates"]["angular_importance_weight_le_2p01"]
    print("PASS: finite-angle summary rejects heavy importance tails")


if __name__ == "__main__":
    main()
