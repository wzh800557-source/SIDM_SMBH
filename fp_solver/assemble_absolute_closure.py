#!/usr/bin/env python3
"""Assemble a gated absolute GNC-to-fluid closure from production results."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


KMS_TO_PC_PER_MYR = 1.022712165045695


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def close(a: float, b: float, tolerance: float = 1.0e-10) -> bool:
    return math.isclose(float(a), float(b), rel_tol=tolerance, abs_tol=0.0)


def maximum_resolution_sensitivity(
    convergence: dict, model: str, metric: str
) -> dict:
    by_dimension = {}
    comparison_source = convergence.get(
        "selected_comparisons", convergence.get("comparisons", {})
    )
    for dimension in ("energy", "angular", "boundary"):
        values = [
            float(item.get("fractional_difference", item.get(
                "full_range_fraction_of_mean", math.nan
            )))
            for item in comparison_source.get(dimension, [])
            if item.get("model") == model and item.get("metric") == metric
        ]
        finite = [value for value in values if math.isfinite(value)]
        by_dimension[dimension] = max(finite) if finite else None
    available = [value for value in by_dimension.values() if value is not None]
    return {
        "by_dimension": by_dimension,
        "maximum": max(available) if available else None,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--remap-json", type=Path, required=True)
    parser.add_argument("--bridge-json", type=Path, required=True)
    parser.add_argument("--normalization-json", type=Path, required=True)
    parser.add_argument("--convergence-json", type=Path, required=True)
    parser.add_argument("--steady-json", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--model-spread-tolerance", type=float, default=0.10)
    parser.add_argument("--thermal-specific-energy-factor", type=float, default=1.5)
    args = parser.parse_args()
    if not 0.0 < args.model_spread_tolerance < 1.0:
        raise ValueError("model-spread tolerance must lie in (0,1)")
    if args.thermal_specific_energy_factor <= 0.0:
        raise ValueError("thermal specific-energy factor must be positive")

    remap = json.loads(args.remap_json.read_text())
    bridge = json.loads(args.bridge_json.read_text())
    norm = json.loads(args.normalization_json.read_text())
    convergence = json.loads(args.convergence_json.read_text())
    steady = json.loads(args.steady_json.read_text())
    metadata = steady["operator_metadata"]

    gates = {
        "production_bh_remap": (
            remap.get("schema") == "adiabatic-bh-remapped-fluid-profile-v1"
            and remap.get("status") == "BH_REMAP_COMPLETE"
        ),
        "mass_matched_hydrostatic_bridge": (
            bridge.get("schema") == "absolute-closure-hydrostatic-bridge-v3"
            and bridge.get("status") == "OK"
            and bridge.get("fluid_subgrid_bridge_gate") == "PASS"
        ),
        "physical_df_normalization": (
            norm.get("schema") == "gnc-normalized-df-v6"
            and norm.get("status") == "PASS"
        ),
        "nonlocal_ej_convergence": (
            convergence.get("schema") == "finite-angle-ej-convergence-v3"
            and convergence.get("status") == "PRODUCTION_CONVERGENCE_PASS"
        ),
        "selected_steady_solver": (
            steady.get("schema") == "finite-angle-ej-steady-state-v2"
            and steady.get("status") == "EJ_STEADY_SMOKE_PASS"
            and "seeds" in metadata
        ),
    }

    mbh = float(bridge["mbh_msun"])
    sigma_over_m = float(bridge["sigma_over_m_cm2_g"])
    w_kms = float(bridge["yukawa_w_kms"])
    identity = convergence.get("identity", {})
    gates["single_physical_identity"] = bool(
        remap.get("output_profile_sha256") == bridge.get("input_profile_sha256")
        and bridge.get("output_profile_sha256") == norm.get("profile_sha256")
        and norm.get("df_table_sha256") == identity.get("df_table_sha256")
        and bridge.get("output_profile_sha256") == identity.get("profile_sha256")
        and metadata.get("profile_sha256") == identity.get("profile_sha256")
        and metadata.get("df_table_sha256") == identity.get("df_table_sha256")
        and close(remap["M_bh_msun"], mbh)
        and close(norm["mbh_msun"], mbh)
        and close(metadata["mbh_msun"], mbh)
        and close(identity["mbh_msun"], mbh)
        and close(remap["sigma0_over_m_cm2_g"], sigma_over_m)
        and close(norm["sigma0_over_m_cm2_g"], sigma_over_m)
        and close(metadata["sigma0_over_m_cm2_g"], sigma_over_m)
        and close(identity["sigma0_over_m_cm2_g"], sigma_over_m)
        and close(remap["yukawa_w_kms"], w_kms)
        and close(norm["w_kms"], w_kms)
        and close(metadata["w_kms"], w_kms)
        and close(identity["w_kms"], w_kms)
    )

    r_boundary = float(bridge["r_in_pc"])
    gates["selected_reporting_surface_is_physical_crossing"] = bool(
        close(norm["r_boundary_pc"], r_boundary, 1.0e-8)
        and close(metadata["reservoir_radius_pc"], r_boundary, 1.0e-8)
    )
    matching_groups = [
        group for group in convergence.get("groups", [])
        if close(group["reservoir_radius_pc"], r_boundary, 1.0e-8)
    ]
    maximum_energy = max(int(group["energy_bins"]) for group in matching_groups)
    maximum_angular = max(
        int(group["angular_bins"])
        for group in matching_groups
        if int(group["energy_bins"]) == maximum_energy
    )
    gates["selected_highest_resolution"] = bool(
        int(metadata["energy_bins"]) == maximum_energy
        and int(metadata["angular_bins"]) == maximum_angular
    )
    selected_groups = [
        group for group in matching_groups
        if int(group["energy_bins"]) == int(metadata["energy_bins"])
        and int(group["angular_bins"]) == int(metadata["angular_bins"])
    ]
    gates["unique_selected_convergence_group"] = len(selected_groups) == 1
    selected_group = selected_groups[0] if len(selected_groups) == 1 else None

    direct = steady["models"]["direct"]
    immediate = steady["models"]["immediate"]
    mdot = float(direct["steady_capture_mdot_msun_per_myr"])
    mdot_immediate = float(immediate["steady_capture_mdot_msun_per_myr"])
    energy_current = direct.get(
        "steady_capture_binding_energy_current_msun_kms2_per_myr"
    )
    energy_current_immediate = immediate.get(
        "steady_capture_binding_energy_current_msun_kms2_per_myr"
    )
    gates["positive_direct_capture_current"] = mdot > 0.0
    gates["capture_energy_measured_on_depleted_solution"] = bool(
        energy_current is not None
        and energy_current_immediate is not None
        and float(energy_current) > 0.0
        and float(energy_current_immediate) > 0.0
    )
    model_spread = abs(mdot_immediate - mdot) / max(
        0.5 * (mdot_immediate + mdot), 1.0e-300
    )
    gates["direct_immediate_model_spread"] = (
        model_spread <= args.model_spread_tolerance
    )
    energy_model_spread = (
        abs(float(energy_current_immediate) - float(energy_current))
        / max(
            0.5 * (
                float(energy_current_immediate) + float(energy_current)
            ),
            1.0e-300,
        )
        if energy_current is not None and energy_current_immediate is not None
        else math.inf
    )
    gates["direct_immediate_capture_energy_model_spread"] = (
        energy_model_spread <= args.model_spread_tolerance
    )

    rho_b = float(norm["rho_boundary_msun_pc3"])
    sigma_b = float(norm["sigma_boundary_kms"])
    mass_scale = (
        4.0 * math.pi * r_boundary**2 * rho_b * sigma_b
        * KMS_TO_PC_PER_MYR
    )
    energy_scale = mass_scale * sigma_b**2
    thermal_sink = (
        -args.thermal_specific_energy_factor * mdot * sigma_b**2
    )
    energy_current_value = float(energy_current) if energy_current is not None else 0.0
    seed_mass_by_model = {
        model: float(selected_group["models"][model]["mass"][
            "seed_full_range_fraction_of_mean"
        ])
        if selected_group is not None else math.inf
        for model in ("direct", "immediate")
    }
    seed_energy_by_model = {
        model: float(selected_group["models"][model][
            "capture_binding_energy"
        ]["seed_full_range_fraction_of_mean"])
        if selected_group is not None else math.inf
        for model in ("direct", "immediate")
    }
    seed_mass_range = max(seed_mass_by_model.values())
    seed_energy_range = max(seed_energy_by_model.values())
    mass_resolution_by_model = {
        model: maximum_resolution_sensitivity(convergence, model, "mass")
        for model in ("direct", "immediate")
    }
    energy_resolution_by_model = {
        model: maximum_resolution_sensitivity(
            convergence, model, "capture_binding_energy"
        )
        for model in ("direct", "immediate")
    }
    mass_resolution_max = max(
        value["maximum"] if value["maximum"] is not None else math.inf
        for value in mass_resolution_by_model.values()
    )
    energy_resolution_max = max(
        value["maximum"] if value["maximum"] is not None else math.inf
        for value in energy_resolution_by_model.values()
    )
    mass_validated_envelope = max(
        seed_mass_range,
        model_spread,
        mass_resolution_max if mass_resolution_max is not None else math.inf,
    )
    energy_validated_envelope = max(
        seed_energy_range,
        energy_model_spread,
        energy_resolution_max if energy_resolution_max is not None else math.inf,
    )

    all_pass = all(gates.values())
    result = {
        "schema": "gnc-fluid-absolute-closure-v1",
        "status": "ABSOLUTE_CLOSURE_MEASURED" if all_pass else "FAIL",
        "gates": gates,
        "identity": {
            "source_collapse_profile_sha256": remap["source_profile_sha256"],
            "fluid_profile_sha256": remap["output_profile_sha256"],
            "bridged_profile_sha256": bridge["output_profile_sha256"],
            "df_table_sha256": norm["df_table_sha256"],
            "remap_json_sha256": sha256_file(args.remap_json),
            "bridge_json_sha256": sha256_file(args.bridge_json),
            "normalization_json_sha256": sha256_file(args.normalization_json),
            "convergence_json_sha256": sha256_file(args.convergence_json),
            "steady_json_sha256": sha256_file(args.steady_json),
            "mbh_msun": mbh,
            "sigma0_over_m_cm2_g": sigma_over_m,
            "w_kms": w_kms,
            "r_outer_pc": r_boundary,
            "rh_pc": float(bridge["rh_pc"]),
            "energy_bins": int(metadata["energy_bins"]),
            "angular_bins": int(metadata["angular_bins"]),
            "seeds": metadata["seeds"],
        },
        "measured_mdot_msun_per_myr": mdot,
        "immediate_absorbing_loss_cone_mdot_msun_per_myr": mdot_immediate,
        "direct_immediate_model_spread_fraction_of_mean": model_spread,
        "model_spread_tolerance": args.model_spread_tolerance,
        "measured_capture_binding_energy_current_msun_kms2_per_myr": (
            energy_current_value
        ),
        "immediate_capture_binding_energy_current_msun_kms2_per_myr": (
            float(energy_current_immediate)
            if energy_current_immediate is not None else None
        ),
        "direct_immediate_capture_energy_spread_fraction_of_mean": (
            energy_model_spread
        ),
        "measured_capture_mean_binding_energy_kms2": (
            energy_current_value / mdot if mdot > 0.0 else None
        ),
        "boundary_density_msun_pc3": rho_b,
        "boundary_sigma1d_kms": sigma_b,
        "mass_flux_scale_msun_per_myr": mass_scale,
        "energy_flux_scale_msun_kms2_per_myr": energy_scale,
        "C_M_measured": mdot / mass_scale,
        "C_E_capture_binding_transport": energy_current_value / energy_scale,
        "thermal_specific_energy_factor": args.thermal_specific_energy_factor,
        "thermal_sink_current_msun_kms2_per_myr": thermal_sink,
        "C_E_thermal_sink": thermal_sink / energy_scale,
        "validated_sensitivity_envelope": {
            "interpretation": (
                "These are maximum validated fractional sensitivities, not "
                "Gaussian one-sigma errors. Components are reported separately "
                "and are not combined in quadrature."
            ),
            "mass_current": {
                "independent_seed_full_range_fraction_of_mean": seed_mass_range,
                "independent_seed_by_model": seed_mass_by_model,
                "resolution_and_boundary": {
                    "by_model": mass_resolution_by_model,
                    "maximum": mass_resolution_max,
                },
                "direct_immediate_model_spread_fraction_of_mean": model_spread,
                "maximum_fractional_sensitivity": mass_validated_envelope,
            },
            "capture_binding_energy_current": {
                "independent_seed_full_range_fraction_of_mean": seed_energy_range,
                "independent_seed_by_model": seed_energy_by_model,
                "resolution_and_boundary": {
                    "by_model": energy_resolution_by_model,
                    "maximum": energy_resolution_max,
                },
                "direct_immediate_model_spread_fraction_of_mean": energy_model_spread,
                "maximum_fractional_sensitivity": energy_validated_envelope,
            },
        },
        "returned_energy_source_ceiling_msun_kms2_per_myr": energy_current_value,
        "physical_fluid_branch": "control versus thermal sink",
        "capture_prescription": (
            "direct no-rescattering current from the steady depleted non-local "
            "EJ solution"
        ),
        "energy_interpretation": (
            "The capture binding-energy current is transported into the black "
            "hole and is not deposited in the halo. Its value is retained only "
            "as a nonphysical full-return ceiling. The physical fluid coupling "
            "removes the isotropic reservoir specific kinetic energy."
        ),
        "mass_coupling": (
            "No fluid mass is deleted during the response run. Evolution stops "
            "before omitted captured mass exceeds the declared interface gates."
        ),
        "inputs": {
            "remap_json": args.remap_json.name,
            "bridge_json": args.bridge_json.name,
            "normalization_json": args.normalization_json.name,
            "convergence_json": args.convergence_json.name,
            "selected_steady_json": args.steady_json.name,
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if all_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
