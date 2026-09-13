#!/usr/bin/env python3
"""Assemble every measured current and flux scale for one scan case."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fp_solver.cross_regime_closure import (  # noqa: E402
    bondi_lambda,
    bondi_rate_msun_per_myr,
    sound_speed_from_sigma_1d,
)

KMS_TO_PC_PER_MYR = 1.0227121650537077


def symmetric_difference(a: float, b: float) -> float:
    denominator = 0.5 * (abs(a) + abs(b))
    return abs(a - b) / denominator if denominator > 0.0 else 0.0


def two_seed_measurement(rows: list[dict], value_key: str, error_key: str) -> dict:
    values = np.asarray([float(row[value_key]) for row in rows])
    errors = np.asarray([float(row[error_key]) for row in rows])
    return {
        "mean": float(np.mean(values)),
        "independent_seed_values": [float(value) for value in values],
        "independent_seed_sensitivity": symmetric_difference(*values),
        "combined_standard_error": float(np.sqrt(np.sum(errors**2)) / len(errors)),
    }


def model_fluxes(steady: dict, model: str) -> dict:
    row = steady["models"][model]
    return {
        "steady_mdot_msun_per_myr": float(
            row["steady_capture_mdot_msun_per_myr"]
        ),
        "isotropic_mdot_msun_per_myr": float(
            row["isotropic_capture_mdot_msun_per_myr"]
        ),
        "steady_binding_current_msun_kms2_per_myr": float(
            row["steady_capture_binding_energy_current_msun_kms2_per_myr"]
        ),
        "isotropic_binding_current_msun_kms2_per_myr": float(
            row["isotropic_capture_binding_energy_current_msun_kms2_per_myr"]
        ),
        "external_supply_msun_per_myr": float(
            row["external_supply_to_hold_connected_states_msun_per_myr"]
        ),
        "outer_reservoir_return_msun_per_myr": float(
            row["outer_reservoir_return_msun_per_myr"]
        ),
        "minimum_internal_occupation": float(row["minimum_internal_occupation"]),
        "maximum_internal_occupation": float(row["maximum_internal_occupation"]),
        "maximum_scaled_residual": float(row["maximum_scaled_residual"]),
        "solver_status": row["status"],
        "solver_gates": row["gates"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    grid = json.loads(args.grid.read_text())
    case = dict(grid["cases"][args.index])
    case_root = args.root / "cases" / case["tag"]
    prep = json.loads((case_root / "PREP_STATUS.json").read_text())
    bridge = json.loads(Path(prep["bridge_json"]).read_text())
    norm = json.loads(Path(prep["normalization_json"]).read_text())
    fluid = json.loads((case_root / "flux" / "fluid_snapshot_flux.json").read_text())

    resolutions = {}
    for energy_bins in (64, 96):
        steady_dir = case_root / "steady" / f"e{energy_bins}_j257"
        combined = json.loads((steady_dir / "steady_combined.json").read_text())
        seed_rows = []
        for seed in (314159, 271828):
            seed_rows.append(
                json.loads((steady_dir / f"steady_s{seed}.json").read_text())
            )
        models = {}
        for model in ("direct", "immediate"):
            selected = model_fluxes(combined, model)
            seed_flux = [model_fluxes(row, model) for row in seed_rows]
            selected["independent_seed_sensitivity"] = {
                key: symmetric_difference(seed_flux[0][key], seed_flux[1][key])
                for key in (
                    "steady_mdot_msun_per_myr",
                    "steady_binding_current_msun_kms2_per_myr",
                )
            }
            selected["seed_values"] = seed_flux
            models[model] = selected
        resolutions[str(energy_bins)] = {
            "energy_bins": energy_bins,
            "angular_bins": 257,
            "models": models,
            "direct_immediate_mass_model_spread": symmetric_difference(
                models["direct"]["steady_mdot_msun_per_myr"],
                models["immediate"]["steady_mdot_msun_per_myr"],
            ),
            "direct_immediate_binding_model_spread": symmetric_difference(
                models["direct"]["steady_binding_current_msun_kms2_per_myr"],
                models["immediate"]["steady_binding_current_msun_kms2_per_myr"],
            ),
        }

    selected = resolutions["96"]
    energy_sensitivity = {}
    for model in ("direct", "immediate"):
        energy_sensitivity[model] = {
            key: symmetric_difference(
                resolutions["64"]["models"][model][key],
                resolutions["96"]["models"][model][key],
            )
            for key in (
                "steady_mdot_msun_per_myr",
                "steady_binding_current_msun_kms2_per_myr",
            )
        }

    operator_runs = []
    for seed in (314159, 271828):
        path = case_root / "operators" / f"e96_j257_s{seed}" / "run.json"
        operator_runs.append(json.loads(path.read_text()))
    operator_isotropic_fluxes = {
        "candidate_mass_current_msun_per_myr": two_seed_measurement(
            operator_runs,
            "candidate_mdot_msun_per_myr",
            "candidate_mdot_standard_error_msun_per_myr",
        ),
        "direct_no_rescatter_mass_current_msun_per_myr": two_seed_measurement(
            operator_runs,
            "direct_no_rescatter_mdot_msun_per_myr",
            "direct_no_rescatter_mdot_standard_error_msun_per_myr",
        ),
        "candidate_binding_energy_current_msun_kms2_per_myr": two_seed_measurement(
            operator_runs,
            "candidate_binding_energy_current_msun_kms2_per_myr",
            "candidate_binding_energy_current_standard_error_msun_kms2_per_myr",
        ),
        "direct_no_rescatter_binding_energy_current_msun_kms2_per_myr": two_seed_measurement(
            operator_runs,
            "accepted_capture_binding_energy_current_msun_kms2_per_myr",
            "accepted_capture_binding_energy_current_standard_error_msun_kms2_per_myr",
        ),
    }

    rb = float(norm["r_boundary_pc"])
    rho_b = float(norm["rho_boundary_msun_pc3"])
    sigma_b = float(norm["sigma_boundary_kms"])
    mass_scale = (
        4.0 * math.pi * rb**2 * rho_b * sigma_b * KMS_TO_PC_PER_MYR
    )
    energy_scale = mass_scale * sigma_b**2
    selected_fluxes = {}
    for model in ("direct", "immediate"):
        values = selected["models"][model]
        mdot = values["steady_mdot_msun_per_myr"]
        binding = values["steady_binding_current_msun_kms2_per_myr"]
        thermal = -1.5 * mdot * sigma_b**2
        selected_fluxes[model] = {
            **values,
            "thermal_sink_current_msun_kms2_per_myr": thermal,
            "mean_capture_binding_energy_kms2": binding / mdot,
            "C_M": mdot / mass_scale,
            "C_E_binding": binding / energy_scale,
            "C_E_thermal_sink": thermal / energy_scale,
            "steady_over_isotropic_mass_current": (
                mdot / values["isotropic_mdot_msun_per_myr"]
            ),
            "steady_over_isotropic_binding_current": (
                binding / values["isotropic_binding_current_msun_kms2_per_myr"]
            ),
        }

    # Nominal polytropic Bondi benchmark for an ideal monatomic SIDM fluid.
    # The local boundary values are substituted for the asymptotic reservoir
    # values, so this is not a Bondi solution for the saved hydrostatic state.
    bondi_gamma = 5.0 / 3.0
    bondi_sound_speed = sound_speed_from_sigma_1d(sigma_b, bondi_gamma)
    bondi_benchmark = bondi_rate_msun_per_myr(
        case["black_hole_mass_msun"],
        rho_b,
        bondi_sound_speed,
        bondi_gamma,
    )
    max_gate = 0.10
    numerical_sensitivities = []
    solver_gates = []
    for model in ("direct", "immediate"):
        numerical_sensitivities.extend(energy_sensitivity[model].values())
        numerical_sensitivities.extend(
            selected["models"][model]["independent_seed_sensitivity"].values()
        )
        solver_gates.extend(selected["models"][model]["solver_gates"].values())
    numerical_pass = max(numerical_sensitivities) <= max_gate
    solver_pass = all(value is True for value in solver_gates)

    direct_mdot = selected_fluxes["direct"]["steady_mdot_msun_per_myr"]
    direct_thermal = selected_fluxes["direct"][
        "thermal_sink_current_msun_kms2_per_myr"
    ]
    fluid_luminosity = float(
        fluid["first_fluid_shell_luminosity_msun_kms2_per_myr"]
    )
    inner_mass = float(fluid["profile_installation"]["inner_shell_mass_msun"])
    inner_sigma = float(fluid["first_fluid_shell_sigma_1d_kms"])
    inner_kinetic_energy = 1.5 * inner_mass * inner_sigma**2

    output = {
        "schema": "sidm-smbh-parameter-scan-flux-budget-v1",
        "status": (
            "SCAN_FLUX_CONVERGED"
            if numerical_pass and solver_pass
            else "SCAN_FLUX_RESOLUTION_EXTENSION_REQUIRED"
        ),
        "updated_at": datetime.now(timezone.utc).isoformat(),
        "case": case,
        "interface": {
            "implementation_status": prep["interface_implementation_status"],
            "bridge_status": prep["bridge_status"],
            "r_in_pc": prep["r_in_pc"],
            "reporting_radius_pc": rb,
            "r_h_pc": prep["rh_pc"],
            "r_in_over_r_h": prep["r_in_over_rh"],
            "first_fluid_shell_radius_pc": prep["first_shell_radius_pc"],
            "all_local_crossings_pc": prep["all_local_crossings_pc"],
            "N_at_first_fluid_shell": prep["N_at_first_fluid_shell"],
            "fp_owned_mass_fraction_of_first_shell": prep[
                "fp_owned_mass_fraction_of_first_shell"
            ],
        },
        "boundary_state": {
            "rho_msun_pc3": rho_b,
            "sigma_1d_kms": sigma_b,
            "mass_flux_scale_msun_per_myr": mass_scale,
            "energy_flux_scale_msun_kms2_per_myr": energy_scale,
        },
        "selected_resolution": {"energy_bins": 96, "angular_bins": 257},
        "operator_isotropic_fluxes": operator_isotropic_fluxes,
        "selected_fluxes": selected_fluxes,
        "all_resolutions": resolutions,
        "energy_resolution_sensitivity_64_to_96": energy_sensitivity,
        "maximum_numerical_sensitivity": max(numerical_sensitivities),
        "numerical_tolerance": max_gate,
        "all_selected_solver_gates_pass": solver_pass,
        "bondi_dimensional_supply_benchmark": {
            "definition": (
                "4 pi lambda_B (G M_bh)^2 rho_b / c_infinity^3, with "
                "gamma=5/3, lambda_B=1/4, and c_infinity=sqrt(gamma) sigma_b"
            ),
            "gamma": bondi_gamma,
            "lambda_B": bondi_lambda(bondi_gamma),
            "sound_speed_infinity_kms": bondi_sound_speed,
            "mdot_msun_per_myr": bondi_benchmark,
            "bondi_over_direct_fp_mass_current": bondi_benchmark / direct_mdot,
            "interpretation": (
                "nominal steady spherical fluid benchmark obtained by treating "
                "the local boundary as the asymptotic reservoir; use only after "
                "collisionality and flow-boundary gates pass"
            ),
        },
        "black_hole_aware_fluid_flux": fluid,
        "flux_comparisons": {
            "direct_thermal_sink_over_signed_first_shell_gravothermal_luminosity": (
                direct_thermal / fluid_luminosity
                if fluid_luminosity != 0.0
                else None
            ),
            "absolute_direct_thermal_sink_over_absolute_first_shell_gravothermal_luminosity": (
                abs(direct_thermal) / abs(fluid_luminosity)
                if fluid_luminosity != 0.0
                else None
            ),
            "inner_shell_kinetic_energy_msun_kms2": inner_kinetic_energy,
            "direct_thermal_sink_timescale_for_inner_shell_myr": (
                inner_kinetic_energy / abs(direct_thermal)
            ),
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
