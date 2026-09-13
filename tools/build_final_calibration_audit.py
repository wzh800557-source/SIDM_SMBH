#!/usr/bin/env python3
"""Assemble the accepted results and open gates into one calibration ledger."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_json(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def build_audit(root: Path = ROOT) -> dict:
    bridge_path = root / "data/production/bridge/bridge.json"
    preparation_path = root / "data/production/prep_manifest.json"
    closure_path = (
        root
        / "data/reference/production/aggregate_v6_mass_reaudit_boundaryfix"
        / "fluid_mass_closure.json"
    )
    convergence_path = closure_path.with_name("ej_convergence.json")
    response_path = (
        root
        / "data/reference/production/fluid_mass_v6_final_20260901_0955"
        / "fluid_response_analysis.json"
    )
    calibration_path = root / "validation/cross_regime_calibration.json"

    bridge = load_json(bridge_path)
    preparation = load_json(preparation_path)
    closure = load_json(closure_path)
    convergence = load_json(convergence_path)
    response = load_json(response_path)
    calibration = load_json(calibration_path)

    require(preparation["status"] == "PRODUCTION_PREPARATION_PASS", "preparation failed")
    require(bridge["status"] == "OK", "hydrostatic bridge failed")
    require(closure["status"] == "FLUID_MASS_CLOSURE_MEASURED", "mass closure failed")
    require(
        convergence["observable_gates"]["mass"]["pass"] is True,
        "mass-current convergence gate failed",
    )
    require(response["status"] == "BLACK_HOLE_AWARE_RESPONSE_COMPLETE", "response failed")
    require(calibration["status"] == "AVAILABLE_BRANCHES_CALIBRATED", "calibration failed")

    fiducial = calibration["fiducial_yukawa_halo"]
    conductive = calibration["conductive_spike"]["digitized_same_form_refit"]
    mass_sensitivity = closure["validated_mass_sensitivity_envelope"]
    energy_pass = convergence["observable_gates"]["capture_binding_energy"]["pass"]
    require(energy_pass is False, "energy-current status changed; review this audit")
    require(
        closure["absolute_two_current_closure_status"] == "INCOMPLETE_OR_FAILED",
        "absolute closure status changed; review this audit",
    )
    require(
        fiducial["radial_collisionality_gate"]["status"]
        == "ORBIT_MEMORY_PRESENT",
        "fiducial solver branch changed; review this audit",
    )

    return {
        "schema": "sidm-smbh-final-calibration-audit-v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "MASS_CURRENT_AND_BRANCH_SELECTION_COMPLETE_ENERGY_CLOSURE_OPEN",
        "scope": (
            "Fixed collapse snapshot, mass-matched black-hole remap, measured "
            "orbit-resolved mass current, archived short fluid sensitivity run, "
            "and calibrated cross-regime branch selection."
        ),
        "accepted_results": {
            "collapse_snapshot_remap": {
                "status": "COMPLETE",
                "preparation_status": preparation["status"],
                "bridge_status": bridge["status"],
                "r_in_pc": bridge["r_in_pc"],
                "r_in_over_r_h": bridge["r_in_over_rh"],
                "fp_owned_fraction_of_first_fluid_shell": bridge[
                    "mass_inside_r_in_fraction_of_first_shell"
                ],
            },
            "nonlocal_depleted_phase_space": {
                "status": "COMPLETE_FOR_CAPTURED_MASS_CURRENT",
                "description": (
                    "The steady finite-angle E-J solution supplies the accepted "
                    "fixed-snapshot mass current. Its capture-energy moment remains "
                    "outside the accepted observable set."
                ),
            },
            "captured_mass_current": {
                "status": "ACCEPTED_FIXED_SNAPSHOT_OBSERVABLE",
                "direct_mdot_msun_per_myr": closure["measured_mdot_msun_per_myr"],
                "immediate_mdot_msun_per_myr": closure[
                    "immediate_absorbing_loss_cone_mdot_msun_per_myr"
                ],
                "C_M": closure["C_M_measured"],
                "maximum_selected_fractional_sensitivity": mass_sensitivity[
                    "maximum_fractional_sensitivity"
                ],
            },
            "black_hole_aware_fluid_response": {
                "status": "NUMERICAL_SENSITIVITY_RUN_COMPLETE",
                "common_duration_in_initial_relaxation_times": response[
                    "common_tau_relax"
                ],
                "sink_over_control_inner_mean_density": response[
                    "inner_mean_density_sink_over_control_at_common_end"
                ],
                "interpretation": (
                    "This run tests the response to a prescribed local kinetic "
                    "sink at the accepted mass rate. It is not an absolute two-way "
                    "evolution because the capture-weighted energy exchange is not "
                    "accepted."
                ),
            },
            "moving_interface_ledger": {
                "status": "IMPLEMENTED_AND_UNIT_TESTED",
                "scope": (
                    "Conservative three-reservoir bookkeeping for supplied mass, "
                    "captured mass, advected energy, luminosity, and boundary work."
                ),
            },
            "cross_regime_calibration": {
                "status": calibration["status"],
                "conductive_spike_refit": {
                    "baseline_msun_per_yr": conductive["baseline_msun_per_yr"],
                    "A": conductive["A"],
                    "B": conductive["B"],
                    "turnover_cm2_g": conductive["sigma_transition_cm2_g"],
                    "rms_msun_per_yr": conductive["rms_msun_per_yr"],
                },
                "fiducial_selected_branch": "orbit-resolved FP",
                "fiducial_N_orb_at_capture": fiducial[
                    "radial_collisionality_gate"
                ]["N_orb_at_capture"],
                "nominal_bondi_benchmark_msun_per_myr": fiducial[
                    "nominal_adiabatic_bondi_benchmark_mdot_msun_per_myr"
                ],
                "nominal_bondi_over_measured_fp": fiducial[
                    "nominal_bondi_benchmark_over_measured_fp"
                ],
                "bondi_authorized_for_fiducial_profile": False,
            },
        },
        "solver_selection": {
            "orbit_resolved": "Use FP where N_orb < 1.",
            "transition": (
                "Use non-orbit-averaged radial kinetics through N_orb near unity."
            ),
            "hydrodynamic_candidate": (
                "Require a continuum path to capture. The released operational "
                "sentinel is N_orb >= 20*pi, equivalent to Kn <= 0.1 only when H=r."
            ),
            "hydrodynamic_endpoint": (
                "Select a conductive spike for hydrostatic boundary conditions or "
                "Bondi for a verified steady, spherical, low-rotation transonic "
                "flow supplied by a maintained reservoir."
            ),
        },
        "open_acceptance_gates": {
            "capture_weighted_energy_current": {
                "status": "NOT_ACCEPTED",
                "reason": "The production energy-resolution gate failed.",
            },
            "transition_regime": {
                "status": "SOLVER_REQUIRED",
                "reason": (
                    "Neither orbit-averaged FP nor hydrodynamics is valid across a "
                    "path whose orbital transport depth is of order unity."
                ),
            },
            "absolute_two_way_evolution": {
                "status": "NOT_AUTHORIZED",
                "reason": (
                    "It requires an accepted energy current and production use of "
                    "the moving-interface mass, energy, luminosity, and work ledger."
                ),
            },
        },
        "scientific_claim_boundary": {
            "supported": [
                "The fiducial fixed-snapshot captured-mass current is accepted.",
                "The fiducial Yukawa nucleus belongs to the orbit-resolved branch.",
                "A hydrostatic conductive spike turns over between LMFP and SMFP limits.",
                "A universal monotonic FP-to-Bondi interpolation is rejected.",
            ],
            "not_supported": [
                "The black hole halts gravothermal collapse.",
                "The black hole heats or cools the halo at the inferred capture-binding energy.",
                "The archived short sink run is a calibrated physical halo trajectory.",
            ],
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    result = build_audit(args.root.resolve())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
