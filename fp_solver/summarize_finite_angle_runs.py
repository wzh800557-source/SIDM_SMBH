#!/usr/bin/env python3
"""Combine independent finite-angle capture estimates and enforce convergence gates."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


def inverse_variance_summary(values: np.ndarray, errors: np.ndarray) -> dict:
    if np.any(~np.isfinite(values)) or np.any(~np.isfinite(errors)) or np.any(errors <= 0.0):
        raise ValueError("finite positive uncertainties are required")
    weight = 1.0 / errors**2
    mean = float(np.sum(weight * values) / np.sum(weight))
    standard_error = float(1.0 / math.sqrt(float(np.sum(weight))))
    chi2 = float(np.sum(((values - mean) / errors) ** 2))
    dof = int(values.size - 1)
    return {
        "weighted_mean": mean,
        "standard_error": standard_error,
        "relative_standard_error": standard_error / abs(mean) if mean != 0.0 else math.inf,
        "chi2": chi2,
        "degrees_of_freedom": dof,
        "reduced_chi2": chi2 / dof if dof > 0 else None,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", nargs="+", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if len(args.runs) < 3:
        raise ValueError("at least three independent runs are required")
    runs = [json.loads(path.read_text()) for path in args.runs]
    if any(run.get("schema") != "finite-angle-loss-cone-capture-v6" for run in runs):
        raise ValueError("all inputs must use finite-angle-loss-cone-capture-v6")
    identity_keys = (
        "profile",
        "profile_sha256",
        "bridge_json_sha256",
        "df_table",
        "df_table_sha256",
        "df_table_schema",
        "mbh_msun",
        "rh_pc",
        "sigma0_over_m_cm2_g",
        "w_kms",
        "r_inner_pc",
        "r_outer_pc",
        "cap_importance_fraction",
        "direction_edge_fraction",
        "direction_edge_power",
        "survival_anomaly_order",
        "survival_collision_quadrature_order",
        "mass_flux_scale_msun_per_myr",
        "energy_flux_scale_msun_kms2_per_myr",
        "boundary_sigma1d_kms",
    )
    reference = {key: runs[0][key] for key in identity_keys}
    angular_importance_model = runs[0]["rows"][0][
        "angular_importance_model"
    ]
    for run in runs[1:]:
        for key, expected in reference.items():
            if run[key] != expected:
                raise ValueError(f"incompatible run metadata for {key}")
        if any(
            row.get("angular_importance_model") != angular_importance_model
            for row in run["rows"]
        ):
            raise ValueError("finite-angle runs use different angular proposals")
    reference["angular_importance_model"] = angular_importance_model
    seeds = [int(run["seed"]) for run in runs]
    radial_bins = [int(run["radial_bins"]) for run in runs]
    direct = inverse_variance_summary(
        np.array([run["accepted_capture_mdot_msun_per_myr"] for run in runs]),
        np.array([
            run["accepted_capture_mdot_standard_error_msun_per_myr"]
            for run in runs
        ]),
    )
    candidate = inverse_variance_summary(
        np.array([run["candidate_mdot_msun_per_myr"] for run in runs]),
        np.array([
            run["candidate_mdot_standard_error_msun_per_myr"] for run in runs
        ]),
    )
    candidate_energy = None
    candidate_energy_key = "candidate_binding_energy_current_msun_kms2_per_myr"
    candidate_energy_error_key = (
        "candidate_binding_energy_current_standard_error_msun_kms2_per_myr"
    )
    if all(candidate_energy_error_key in run for run in runs):
        candidate_energy = inverse_variance_summary(
            np.array([run[candidate_energy_key] for run in runs]),
            np.array([run[candidate_energy_error_key] for run in runs]),
        )
    direct_energy = None
    direct_energy_key = "accepted_capture_binding_energy_current_msun_kms2_per_myr"
    direct_energy_error_key = (
        "accepted_capture_binding_energy_current_standard_error_msun_kms2_per_myr"
    )
    if all(direct_energy_error_key in run for run in runs):
        direct_energy = inverse_variance_summary(
            np.array([run[direct_energy_key] for run in runs]),
            np.array([run[direct_energy_error_key] for run in runs]),
        )
    maximum_density_error = max(
        float(run["rho_df_profile_max_abs_relative_error"]) for run in runs
    )
    maximum_unresolved_fraction = max(
        float(run["unresolved_kepler_domain_fraction_of_candidate"])
        for run in runs
    )
    maximum_angular_importance_weight = max(
        float(row["angular_importance_weight_range"][1])
        for run in runs for row in run["rows"]
    )
    depletion_fractions = [
        float(run["candidate_current_fraction_from_source_j_below_5_jlc"])
        for run in runs
    ]
    far_source = {}
    if all("candidate_far_source_currents" in run for run in runs):
        for key in runs[0]["candidate_far_source_currents"]:
            far_source[key] = inverse_variance_summary(
                np.array([
                    run["candidate_far_source_currents"][key]["mdot_msun_per_myr"]
                    for run in runs
                ]),
                np.array([
                    run["candidate_far_source_currents"][key]["standard_error_msun_per_myr"]
                    for run in runs
                ]),
            )
    individual_pass = all(
        run["status"] == "ISOTROPIC_INJECTION_CEILING_CONVERGED"
        for run in runs
    )
    independent_seeds = len(set(seeds)) == len(seeds)
    radial_resolution_check = len(set(radial_bins)) >= 2
    resolution_groups = {}
    for resolution in sorted(set(radial_bins)):
        select = np.asarray([value == resolution for value in radial_bins])
        resolution_groups[str(resolution)] = inverse_variance_summary(
            np.asarray([
                run["candidate_mdot_msun_per_myr"] for run in runs
            ])[select],
            np.asarray([
                run["candidate_mdot_standard_error_msun_per_myr"]
                for run in runs
            ])[select],
        )
    resolution_pairwise_z = []
    group_values = list(resolution_groups.values())
    for i in range(len(group_values)):
        for j in range(i + 1, len(group_values)):
            resolution_pairwise_z.append(
                abs(group_values[i]["weighted_mean"] - group_values[j]["weighted_mean"])
                / math.hypot(
                    group_values[i]["standard_error"],
                    group_values[j]["standard_error"],
                )
            )
    maximum_resolution_difference_sigma = (
        max(resolution_pairwise_z) if resolution_pairwise_z else math.inf
    )
    consistency = (
        candidate["reduced_chi2"] is not None
        and candidate["reduced_chi2"] <= 3.0
    )
    gates = {
        "all_individual_runs_pass": individual_pass,
        "independent_seeds": independent_seeds,
        "at_least_two_radial_resolutions": radial_resolution_check,
        "radial_resolution_difference_le_3sigma": (
            maximum_resolution_difference_sigma <= 3.0
        ),
        "candidate_ceiling_reduced_chi2_le_3": consistency,
        "combined_relative_standard_error_le_0p05": (
            candidate["relative_standard_error"] <= 0.05
        ),
        "direct_component_relative_standard_error_le_0p05": (
            direct["relative_standard_error"] <= 0.05
        ),
        "candidate_energy_reduced_chi2_le_3": (
            candidate_energy is not None
            and candidate_energy["reduced_chi2"] is not None
            and candidate_energy["reduced_chi2"] <= 3.0
        ),
        "candidate_energy_relative_standard_error_le_0p05": (
            candidate_energy is not None
            and candidate_energy["relative_standard_error"] <= 0.05
        ),
        "density_roundtrip_le_0p01": maximum_density_error <= 0.01,
        "unresolved_kepler_fraction_le_0p01": maximum_unresolved_fraction <= 0.01,
        "angular_importance_weight_le_2p01": (
            maximum_angular_importance_weight <= 2.01
        ),
    }
    result = {
        "schema": "finite-angle-loss-cone-convergence-v6",
        "status": (
            "FINITE_ANGLE_ISOTROPIC_INJECTION_CEILING_CONVERGED"
            if all(gates.values()) else "FAIL"
        ),
        "run_files": [str(path) for path in args.runs],
        "run_seeds": seeds,
        "radial_bins": radial_bins,
        "pairs_per_bin": [int(run["pairs_per_bin"]) for run in runs],
        "candidate_current_by_radial_resolution": resolution_groups,
        "maximum_radial_resolution_difference_sigma": (
            maximum_resolution_difference_sigma
        ),
        "identity": reference,
        "candidate_injection_mdot_msun_per_myr": candidate,
        "direct_no_rescatter_mdot_msun_per_myr": direct,
        "candidate_binding_energy_current_msun_kms2_per_myr": candidate_energy,
        "direct_no_rescatter_binding_energy_current_msun_kms2_per_myr": (
            direct_energy
        ),
        "far_source_current_diagnostics": far_source,
        "accepted_to_candidate_current_ratio": (
            direct["weighted_mean"] / candidate["weighted_mean"]
        ),
        "maximum_density_roundtrip_error": maximum_density_error,
        "maximum_unresolved_kepler_domain_fraction": maximum_unresolved_fraction,
        "maximum_angular_importance_weight": maximum_angular_importance_weight,
        "source_j_below_5_jlc_fraction_range": [
            min(depletion_fractions), max(depletion_fractions)
        ],
        "absolute_closure_status": "BLOCKED_BY_NONLOCAL_J_DEPLETION_SOLUTION",
        "gates": gates,
    }
    mass_scale = float(reference["mass_flux_scale_msun_per_myr"])
    sigma_boundary = float(reference["boundary_sigma1d_kms"])
    result.update({
        "isotropic_upper_mdot_msun_per_myr": candidate["weighted_mean"],
        "isotropic_upper_mdot_standard_error_msun_per_myr": candidate["standard_error"],
        "direct_no_rescatter_mdot_value_msun_per_myr": direct["weighted_mean"],
        "direct_no_rescatter_mdot_standard_error_msun_per_myr": direct["standard_error"],
        "C_M_isotropic_injection_ceiling": candidate["weighted_mean"] / mass_scale,
        "C_M_isotropic_injection_ceiling_standard_error": (
            candidate["standard_error"] / mass_scale
        ),
        "C_M_direct_no_rescatter": direct["weighted_mean"] / mass_scale,
        "thermal_sink_current_msun_kms2_per_myr": (
            -candidate["weighted_mean"] * sigma_boundary**2
        ),
        "thermal_sink_upper_msun_kms2_per_myr": (
            -candidate["weighted_mean"] * sigma_boundary**2
        ),
        "thermal_sink_current_standard_error_msun_kms2_per_myr": (
            candidate["standard_error"] * sigma_boundary**2
        ),
    })
    if candidate_energy is not None:
        result.update({
            "returned_energy_source_ceiling_msun_kms2_per_myr": (
                candidate_energy["weighted_mean"]
            ),
            "returned_energy_source_ceiling_standard_error_msun_kms2_per_myr": (
                candidate_energy["standard_error"]
            ),
            "candidate_mean_binding_energy_kms2": (
                candidate_energy["weighted_mean"] / candidate["weighted_mean"]
            ),
            "C_E_candidate_binding_ceiling": (
                candidate_energy["weighted_mean"]
                / float(reference["energy_flux_scale_msun_kms2_per_myr"])
            ),
            "C_E_thermal_sink": -candidate["weighted_mean"] / mass_scale,
        })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return (
        0
        if result["status"] == "FINITE_ANGLE_ISOTROPIC_INJECTION_CEILING_CONVERGED"
        else 2
    )


if __name__ == "__main__":
    raise SystemExit(main())
