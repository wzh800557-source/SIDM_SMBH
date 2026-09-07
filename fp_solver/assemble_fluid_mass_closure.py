#!/usr/bin/env python3
"""Assemble the converged mass current needed by the physical fluid sink.

The capture-binding-energy current is transported into the black hole and is
not the luminosity applied to the SIDM fluid.  This assembler therefore admits
a mass-current closure only when the mass observable passes every numerical,
seed, energy-grid, angular-grid, and reporting-surface gate.  It records but
does not waive a failed binding-energy-current gate.  The resulting artifact is
not an absolute two-current closure.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from assemble_absolute_closure import (
    KMS_TO_PC_PER_MYR,
    close,
    maximum_resolution_sensitivity,
    sha256_file,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--remap-json", type=Path, required=True)
    parser.add_argument("--bridge-json", type=Path, required=True)
    parser.add_argument("--normalization-json", type=Path, required=True)
    parser.add_argument("--convergence-json", type=Path, required=True)
    parser.add_argument("--steady-json", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--model-spread-tolerance", type=float, default=0.10)
    parser.add_argument(
        "--boundary-radius-tolerance",
        type=float,
        default=1.0e-3,
        help=(
            "maximum relative mismatch between the independently reconstructed "
            "hydrostatic-bridge crossing and the GNC reporting surface"
        ),
    )
    parser.add_argument("--thermal-specific-energy-factor", type=float, default=1.5)
    args = parser.parse_args()
    if not 0.0 < args.model_spread_tolerance < 1.0:
        raise ValueError("model-spread tolerance must lie in (0,1)")
    if not 0.0 < args.boundary_radius_tolerance < 1.0:
        raise ValueError("boundary-radius tolerance must lie in (0,1)")
    if args.thermal_specific_energy_factor <= 0.0:
        raise ValueError("thermal specific-energy factor must be positive")

    remap = json.loads(args.remap_json.read_text())
    bridge = json.loads(args.bridge_json.read_text())
    norm = json.loads(args.normalization_json.read_text())
    convergence = json.loads(args.convergence_json.read_text())
    steady = json.loads(args.steady_json.read_text())
    metadata = steady.get("operator_metadata", {})
    observable_gates = convergence.get("observable_gates", {})
    mass_observable = observable_gates.get("mass", {})
    binding_observable = observable_gates.get("capture_binding_energy", {})

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
        "mass_observable_converged": (
            convergence.get("schema") == "finite-angle-ej-convergence-v3"
            and mass_observable.get("status")
            == "PRODUCTION_OBSERVABLE_CONVERGENCE_PASS"
            and mass_observable.get("pass") is True
        ),
        "binding_energy_status_recorded_separately": (
            isinstance(binding_observable, dict)
            and binding_observable.get("status") in {
                "PRODUCTION_OBSERVABLE_CONVERGENCE_PASS",
                "INCOMPLETE_OR_FAILED",
            }
            and isinstance(binding_observable.get("pass"), bool)
        ),
        "selected_steady_solver": (
            steady.get("schema") == "finite-angle-ej-steady-state-v2"
            and steady.get("status") == "EJ_STEADY_SMOKE_PASS"
            and isinstance(metadata.get("seeds"), list)
            and len(metadata["seeds"]) >= 2
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

    # The GNC operator, DF normalization, and convergence ledger all use the
    # reporting surface chosen by boundary_criterion_scan.py.  The bridge
    # independently reconstructs the same N_orb=1 crossing on its hydrostatic
    # continuation.  Those two root finds need only agree within the declared
    # interface tolerance (also used by interface_bridge.py); forcing bitwise
    # agreement incorrectly rejects harmless interpolation-level differences.
    r_boundary = float(norm["r_boundary_pc"])
    r_bridge = float(bridge["r_in_pc"])
    radius_mismatch = abs(r_boundary / r_bridge - 1.0)
    gates["selected_reporting_surface_is_physical_crossing"] = bool(
        radius_mismatch <= args.boundary_radius_tolerance
        and close(metadata["reservoir_radius_pc"], r_boundary, 1.0e-8)
    )
    matching_groups = [
        group for group in convergence.get("groups", [])
        if close(group["reservoir_radius_pc"], r_boundary, 1.0e-8)
    ]
    maximum_energy = max(
        (int(group["energy_bins"]) for group in matching_groups),
        default=-1,
    )
    maximum_angular = max(
        (
            int(group["angular_bins"])
            for group in matching_groups
            if int(group["energy_bins"]) == maximum_energy
        ),
        default=-1,
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
    model_spread = abs(mdot_immediate - mdot) / max(
        0.5 * (mdot_immediate + mdot), 1.0e-300
    )
    gates["positive_direct_capture_current"] = mdot > 0.0
    gates["direct_immediate_mass_model_spread"] = (
        model_spread <= args.model_spread_tolerance
    )
    gates["capture_binding_energy_not_used_as_fluid_luminosity"] = True

    seed_mass_by_model = {
        model: float(selected_group["models"][model]["mass"][
            "seed_full_range_fraction_of_mean"
        ])
        if selected_group is not None else math.inf
        for model in ("direct", "immediate")
    }
    seed_mass_range = max(seed_mass_by_model.values())
    mass_resolution_by_model = {
        model: maximum_resolution_sensitivity(convergence, model, "mass")
        for model in ("direct", "immediate")
    }
    mass_resolution_max = max(
        record["maximum"] if record["maximum"] is not None else math.inf
        for record in mass_resolution_by_model.values()
    )
    validated_mass_envelope = max(
        seed_mass_range, model_spread, mass_resolution_max
    )
    gates["validated_mass_sensitivity_le_tolerance"] = (
        validated_mass_envelope <= args.model_spread_tolerance
    )

    rho_b = float(norm["rho_boundary_msun_pc3"])
    sigma_b = float(norm["sigma_boundary_kms"])
    mass_scale = (
        4.0 * math.pi * r_boundary**2 * rho_b * sigma_b
        * KMS_TO_PC_PER_MYR
    )
    energy_scale = mass_scale * sigma_b**2
    thermal_sink = -args.thermal_specific_energy_factor * mdot * sigma_b**2
    all_pass = bool(gates and all(gates.values()))
    result = {
        "schema": "gnc-fluid-mass-closure-v1",
        "status": "FLUID_MASS_CLOSURE_MEASURED" if all_pass else "FAIL",
        "scientific_scope": (
            "CONVERGED_MASS_CURRENT_FOR_FIXED_SNAPSHOT_FLUID_RESPONSE"
        ),
        "absolute_two_current_closure_status": convergence.get("status"),
        "capture_binding_energy_observable_status": binding_observable.get(
            "status"
        ),
        "capture_binding_energy_current_used_by_fluid": False,
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
            "bridge_r_in_pc": r_bridge,
            "closure_to_bridge_radius_relative_error": radius_mismatch,
            "boundary_radius_tolerance": args.boundary_radius_tolerance,
            "rh_pc": float(bridge["rh_pc"]),
            "energy_bins": int(metadata["energy_bins"]),
            "angular_bins": int(metadata["angular_bins"]),
            "seeds": metadata["seeds"],
        },
        "measured_mdot_msun_per_myr": mdot,
        "immediate_absorbing_loss_cone_mdot_msun_per_myr": mdot_immediate,
        "direct_immediate_model_spread_fraction_of_mean": model_spread,
        "model_spread_tolerance": args.model_spread_tolerance,
        "boundary_density_msun_pc3": rho_b,
        "boundary_sigma1d_kms": sigma_b,
        "mass_flux_scale_msun_per_myr": mass_scale,
        "energy_flux_scale_msun_kms2_per_myr": energy_scale,
        "C_M_measured": mdot / mass_scale,
        "thermal_specific_energy_factor": args.thermal_specific_energy_factor,
        "thermal_sink_current_msun_kms2_per_myr": thermal_sink,
        "C_E_thermal_sink": thermal_sink / energy_scale,
        "validated_mass_sensitivity_envelope": {
            "interpretation": (
                "This is the largest validated fractional sensitivity across "
                "independent seeds, selected resolution and reporting-surface "
                "tests, and the two capture prescriptions. It is not a "
                "Gaussian one-sigma error."
            ),
            "independent_seed_full_range_fraction_of_mean": seed_mass_range,
            "independent_seed_by_model": seed_mass_by_model,
            "resolution_and_boundary": {
                "by_model": mass_resolution_by_model,
                "maximum": mass_resolution_max,
            },
            "direct_immediate_model_spread_fraction_of_mean": model_spread,
            "maximum_fractional_sensitivity": validated_mass_envelope,
        },
        "physical_fluid_branch": "control versus local thermal sink",
        "capture_prescription": (
            "direct no-rescattering mass current from the steady depleted "
            "non-local EJ solution"
        ),
        "energy_interpretation": (
            "The unresolved capture-binding-energy moment is carried into the "
            "black hole and is excluded from the halo luminosity. The fluid "
            "sink removes the local isotropic kinetic energy associated with "
            "the converged captured mass current."
        ),
        "mass_coupling": (
            "The short response does not delete a fluid shell. It stops before "
            "the omitted captured mass exceeds the declared interface gates."
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
