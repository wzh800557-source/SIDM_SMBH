#!/usr/bin/env python3
"""Test whether finite-angle currents are independent of the reporting surface."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


def weighted_consistency(values: np.ndarray, errors: np.ndarray) -> dict:
    values = np.asarray(values, float)
    errors = np.asarray(errors, float)
    if (
        values.size < 3
        or values.shape != errors.shape
        or np.any(~np.isfinite(values))
        or np.any(~np.isfinite(errors))
        or np.any(errors <= 0.0)
    ):
        raise ValueError("three finite measurements with positive errors are required")
    weight = 1.0 / errors**2
    mean = float(np.sum(weight * values) / np.sum(weight))
    standard_error = float(1.0 / math.sqrt(float(np.sum(weight))))
    chi2 = float(np.sum(((values - mean) / errors) ** 2))
    dof = int(values.size - 1)
    pairwise_z = []
    for i in range(values.size):
        for j in range(i + 1, values.size):
            pairwise_z.append(
                abs(float(values[i] - values[j]))
                / math.hypot(float(errors[i]), float(errors[j]))
            )
    return {
        "weighted_mean": mean,
        "standard_error": standard_error,
        "chi2": chi2,
        "degrees_of_freedom": dof,
        "reduced_chi2": chi2 / dof,
        "maximum_pairwise_difference_sigma": max(pairwise_z),
        "fractional_span": float((np.max(values) - np.min(values)) / abs(mean)),
    }


def compare_runs(runs: list[dict]) -> dict:
    if len(runs) < 3:
        raise ValueError("at least three boundary runs are required")
    for run in runs:
        if run.get("schema") != "finite-angle-loss-cone-capture-v6":
            raise ValueError("boundary run uses a stale finite-angle schema")
    shared_keys = (
        "profile_sha256",
        "bridge_json_sha256",
        "df_table_sha256",
        "df_table_schema",
        "mbh_msun",
        "rh_pc",
        "sigma0_over_m_cm2_g",
        "w_kms",
        "r_inner_pc",
        "cap_importance_fraction",
        "direction_edge_fraction",
        "direction_edge_power",
        "survival_anomaly_order",
        "survival_collision_quadrature_order",
    )
    identity = {key: runs[0][key] for key in shared_keys}
    angular_importance_model = runs[0]["rows"][0][
        "angular_importance_model"
    ]
    for run in runs[1:]:
        for key, expected in identity.items():
            if run[key] != expected:
                raise ValueError(f"boundary runs differ in {key}")
        if any(
            row.get("angular_importance_model") != angular_importance_model
            for row in run["rows"]
        ):
            raise ValueError("boundary runs use different angular proposals")
    identity["angular_importance_model"] = angular_importance_model
    radii = np.asarray([run["r_outer_pc"] for run in runs], float)
    seeds = [int(run["seed"]) for run in runs]
    candidate = weighted_consistency(
        np.asarray([run["candidate_mdot_msun_per_myr"] for run in runs]),
        np.asarray([
            run["candidate_mdot_standard_error_msun_per_myr"] for run in runs
        ]),
    )
    direct = weighted_consistency(
        np.asarray([run["direct_no_rescatter_mdot_msun_per_myr"] for run in runs]),
        np.asarray([
            run["direct_no_rescatter_mdot_standard_error_msun_per_myr"]
            for run in runs
        ]),
    )
    candidate_energy = weighted_consistency(
        np.asarray([
            run["candidate_binding_energy_current_msun_kms2_per_myr"]
            for run in runs
        ]),
        np.asarray([
            run[
                "candidate_binding_energy_current_standard_error_msun_kms2_per_myr"
            ]
            for run in runs
        ]),
    )
    gates = {
        "all_runs_numerically_converged": all(
            run["status"] == "ISOTROPIC_INJECTION_CEILING_CONVERGED"
            for run in runs
        ),
        "three_distinct_boundaries": len(set(radii.tolist())) >= 3,
        "independent_seeds": len(set(seeds)) == len(seeds),
        "candidate_reduced_chi2_le_3": candidate["reduced_chi2"] <= 3.0,
        "candidate_maximum_pairwise_difference_le_3sigma": (
            candidate["maximum_pairwise_difference_sigma"] <= 3.0
        ),
        "direct_reduced_chi2_le_3": direct["reduced_chi2"] <= 3.0,
        "candidate_energy_reduced_chi2_le_3": (
            candidate_energy["reduced_chi2"] <= 3.0
        ),
    }
    return {
        "schema": "finite-angle-boundary-invariance-v1",
        "status": (
            "FINITE_ANGLE_BOUNDARY_INVARIANCE_PASS"
            if all(gates.values()) else "FAIL"
        ),
        "identity": identity,
        "r_outer_pc": radii.tolist(),
        "seeds": seeds,
        "candidate_injection_current": candidate,
        "direct_no_rescatter_current": direct,
        "candidate_binding_energy_current": candidate_energy,
        "gates": gates,
        "note": (
            "All runs use one profile and one normalized DF. Only the outer "
            "reporting surface and the independent Monte Carlo seed change."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", nargs="+", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    result = compare_runs([json.loads(path.read_text()) for path in args.runs])
    result["run_files"] = [str(path) for path in args.runs]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "FINITE_ANGLE_BOUNDARY_INVARIANCE_PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
