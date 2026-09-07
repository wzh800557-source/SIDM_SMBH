#!/usr/bin/env python3
"""Compare black-hole-aware control and capture-sink fluid trajectories."""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import re
from pathlib import Path


VALID_CONTROL_STOPS = {"tau_end", "rho_factor_end"}
VALID_SINK_STOPS = VALID_CONTROL_STOPS | {
    "unmodeled_fluid_mass_loss_limit",
    "unmodeled_gnc_mass_loss_limit",
}

CANONICAL_RESPONSE_FIELDS = (
    "rho_inner_mean_msun_pc3",
    "sigma_inner_1d_kms",
    "r_inner_pc",
)
LEGACY_RESPONSE_ALIASES = {
    "rho_inner_mean_msun_pc3": "rho_c_msun_pc3",
    "sigma_inner_1d_kms": "sigma_c_kms",
    "r_inner_pc": "r0_pc",
}
INNER_DENSITY_DEFINITION = (
    "mean density inside the innermost resolved Lagrangian shell, "
    "3*M(<r_inner)/(4*pi*r_inner^3)"
)
INNER_DISPERSION_DEFINITION = (
    "one-dimensional velocity dispersion in the innermost resolved "
    "Lagrangian shell"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def is_sha256(value: object) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def close(a: object, b: object, tolerance: float = 1.0e-10) -> bool:
    try:
        left = float(a)
        right = float(b)
    except (TypeError, ValueError):
        return False
    return (
        math.isfinite(left)
        and math.isfinite(right)
        and math.isclose(left, right, rel_tol=tolerance, abs_tol=0.0)
    )


def read_rows(path: Path) -> dict[str, list[dict]]:
    grouped: dict[str, list[dict]] = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            parsed = {
                key: (value if key == "branch" else float(value))
                for key, value in row.items()
            }
            grouped.setdefault(parsed["branch"], []).append(parsed)
    for rows in grouped.values():
        rows.sort(key=lambda item: item["tau_relax"])
    return grouped


def log_interp(x: list[float], y: list[float], target: float) -> float:
    if len(x) != len(y) or len(x) < 2:
        raise ValueError("logarithmic response interpolation needs aligned samples")
    if any(not math.isfinite(value) for value in x + y + [target]):
        raise ValueError("logarithmic response interpolation needs finite values")
    if any(value <= 0.0 for value in y):
        raise ValueError("logarithmic response interpolation needs positive values")
    if any(right <= left for left, right in zip(x, x[1:])):
        raise ValueError("interpolation coordinates must increase strictly")
    if target < x[0] or target > x[-1]:
        raise ValueError("interpolation target lies outside the trajectory")
    if target == x[0]:
        return y[0]
    if target == x[-1]:
        return y[-1]
    lower = bisect.bisect_right(x, target) - 1
    fraction = (target - x[lower]) / (x[lower + 1] - x[lower])
    return math.exp(
        math.log(y[lower])
        + fraction * (math.log(y[lower + 1]) - math.log(y[lower]))
    )


def validate_resolved_diagnostics(
    summary: dict,
    rows: dict[str, list[dict]],
    density_tolerance: float,
) -> float:
    metadata = summary.get("resolved_diagnostics")
    expected_metadata = {
        "density_field": "rho_inner_mean_msun_pc3",
        "density_definition": INNER_DENSITY_DEFINITION,
        "density_is_extrapolated_central_value": False,
        "dispersion_field": "sigma_inner_1d_kms",
        "dispersion_definition": INNER_DISPERSION_DEFINITION,
        "radius_field": "r_inner_pc",
    }
    if not isinstance(metadata, dict) or any(
        metadata.get(key) != value for key, value in expected_metadata.items()
    ):
        raise RuntimeError("fluid feedback has an ambiguous inner-density diagnostic")

    maximum_density_error = 0.0
    for branch, values in rows.items():
        for row in values:
            for canonical, legacy in LEGACY_RESPONSE_ALIASES.items():
                if canonical not in row or legacy not in row:
                    raise RuntimeError(
                        f"{branch} trajectory lacks {canonical} or its legacy alias"
                    )
                if not close(row[canonical], row[legacy], tolerance=5.0e-12):
                    raise RuntimeError(
                        f"{branch} trajectory has inconsistent {canonical} aliases"
                    )
            shell_mass = float(row.get("inner_shell_mass_msun", math.nan))
            radius = float(row["r_inner_pc"])
            expected_density = 3.0 * shell_mass / (4.0 * math.pi * radius**3)
            density_error = abs(
                float(row["rho_inner_mean_msun_pc3"]) / expected_density - 1.0
            ) if expected_density > 0.0 else math.inf
            maximum_density_error = max(maximum_density_error, density_error)
            if (
                not math.isfinite(shell_mass)
                or shell_mass <= 0.0
                or not math.isfinite(density_error)
                or density_error > density_tolerance
            ):
                raise RuntimeError(
                    f"{branch} trajectory does not satisfy the innermost-shell "
                    "mean-density definition"
                )
    return maximum_density_error


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--trajectories", type=Path, required=True)
    parser.add_argument("--closure-json", type=Path, required=True)
    parser.add_argument("--profile", type=Path, required=True)
    parser.add_argument("--bridge-json", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument(
        "--mean-density-consistency-tolerance",
        type=float,
        default=1.0e-5,
        help=(
            "maximum relative mismatch between the stored innermost-shell "
            "density and 3M/(4 pi r^3); this accommodates the accepted "
            "finite hydrostatic-remap tolerance"
        ),
    )
    args = parser.parse_args()
    if not 0.0 < args.mean_density_consistency_tolerance < 1.0e-3:
        raise ValueError("mean-density consistency tolerance must lie in (0,1e-3)")

    summary = json.loads(args.summary.read_text())
    closure = json.loads(args.closure_json.read_text())
    bridge = json.loads(args.bridge_json.read_text())
    rows = read_rows(args.trajectories)
    absolute_closure = (
        closure.get("schema") == "gnc-fluid-absolute-closure-v1"
        and closure.get("status") == "ABSOLUTE_CLOSURE_MEASURED"
    )
    mass_closure = (
        closure.get("schema") == "gnc-fluid-mass-closure-v1"
        and closure.get("status") == "FLUID_MASS_CLOSURE_MEASURED"
        and closure.get("capture_binding_energy_current_used_by_fluid") is False
    )
    expected_summary_status = (
        "MEASURED_FLUID_FEEDBACK_COMPLETE"
        if absolute_closure else "MASS_CLOSURE_FLUID_RESPONSE_COMPLETE"
        if mass_closure else None
    )
    if (
        summary.get("schema") != "gnc-fluid-feedback-v5"
        or summary.get("status") != expected_summary_status
    ):
        raise RuntimeError("fluid feedback did not use an accepted closure")
    if not (absolute_closure or mass_closure):
        raise RuntimeError(
            "fluid response requires an accepted absolute or mass-current closure"
        )
    closure_gates = closure.get("gates")
    if (
        not isinstance(closure_gates, dict)
        or not closure_gates
        or any(value is not True for value in closure_gates.values())
    ):
        raise RuntimeError("fluid response closure has a missing or failed gate")
    if (
        bridge.get("schema") != "absolute-closure-hydrostatic-bridge-v3"
        or bridge.get("status") != "OK"
    ):
        raise RuntimeError("fluid response requires the accepted hydrostatic bridge")
    if set(rows) != {"control", "sink"}:
        raise RuntimeError(
            "the accepted physical response must contain only control and sink"
        )
    resolved_density_max_relative_error = validate_resolved_diagnostics(
        summary,
        rows,
        args.mean_density_consistency_tolerance,
    )
    resolved_diagnostics_verified = True

    control = rows["control"]
    sink = rows["sink"]
    for branch, values in (("control", control), ("sink", sink)):
        tau = [float(row["tau_relax"]) for row in values]
        if (
            len(tau) < 2
            or not all(math.isfinite(value) for value in tau)
            or tau[0] != 0.0
            or any(right <= left for left, right in zip(tau, tau[1:]))
        ):
            raise RuntimeError(f"{branch} trajectory times are not strictly increasing")
    common_tau = min(control[-1]["tau_relax"], sink[-1]["tau_relax"])
    target_tau = [
        float(row["tau_relax"])
        for row in sink if float(row["tau_relax"]) < float(common_tau)
    ]
    if not target_tau or not close(target_tau[-1], common_tau):
        target_tau.append(float(common_tau))
    if len(target_tau) < 2:
        raise RuntimeError("fluid trajectories have no common evolved interval")
    c_tau = [float(row["tau_relax"]) for row in control]
    s_tau = [float(row["tau_relax"]) for row in sink]
    comparison_rows = []
    fields = CANONICAL_RESPONSE_FIELDS
    ratios = {field: [] for field in fields}
    for tau in target_tau:
        out = {"tau_relax": float(tau)}
        for field in fields:
            control_value = float(log_interp(
                c_tau,
                [float(row[field]) for row in control],
                float(tau),
            ))
            sink_value = float(log_interp(
                s_tau,
                [float(row[field]) for row in sink],
                float(tau),
            ))
            ratio = sink_value / control_value
            out[f"control_{field}"] = control_value
            out[f"sink_{field}"] = sink_value
            out[f"sink_over_control_{field}"] = ratio
            ratios[field].append(ratio)
        comparison_rows.append(out)

    branch_summaries = summary["branches"]
    input_identity = summary.get("input_identity", {})
    closure_identity = closure.get("identity", {})
    identity_gates = {
        "summary_hashes_well_formed": all(
            is_sha256(input_identity.get(key))
            for key in (
                "closure_json_sha256",
                "fluid_profile_sha256",
                "bridge_json_sha256",
                "bridged_profile_sha256",
            )
        ),
        "closure_file_identity": (
            input_identity.get("closure_json_sha256") == sha256(args.closure_json)
        ),
        "profile_file_identity": (
            input_identity.get("fluid_profile_sha256") == sha256(args.profile)
            and closure_identity.get("fluid_profile_sha256") == sha256(args.profile)
        ),
        "bridge_file_identity": (
            input_identity.get("bridge_json_sha256") == sha256(args.bridge_json)
            and closure_identity.get("bridge_json_sha256") == sha256(args.bridge_json)
        ),
        "bridged_profile_identity": (
            input_identity.get("bridged_profile_sha256")
            == bridge.get("output_profile_sha256")
            == closure_identity.get("bridged_profile_sha256")
        ),
    }
    applied = summary.get("applied_currents", {})
    fixed_mode = summary.get("closure_update_mode")
    target_mdot = float(applied.get("Mdot_msun_per_myr", math.nan))
    target_sink_luminosity = float(
        applied.get("sink_luminosity_msun_kms2_per_myr", math.nan)
    )
    fixed_trajectory_currents = (
        all(
            close(row.get("L_inner_msun_kms2_per_myr"), 0.0)
            for row in control
        )
        and all(
            close(
                row.get("L_inner_msun_kms2_per_myr"),
                target_sink_luminosity,
            )
            for row in sink
        )
    )
    captured_mass_accounting = (
        all(close(row.get("captured_mass_msun"), 0.0) for row in control)
        and all(
            close(
                row.get("captured_mass_msun"),
                target_mdot * float(row.get("time_myr", math.nan)),
            )
            for row in sink
        )
    )
    black_hole_mass_fixed = all(
        close(row.get("M_bh_msun"), summary.get("M_bh_msun"))
        for values in (control, sink)
        for row in values
    )
    summary_trajectory_endpoints_match = all(
        close(
            branch_summaries[name].get("final_tau_relax"),
            values[-1].get("tau_relax"),
        )
        and close(
            branch_summaries[name].get("final_time_myr"),
            values[-1].get("time_myr"),
        )
        and close(
            branch_summaries[name].get("final_rho_inner_mean_msun_pc3"),
            values[-1].get("rho_inner_mean_msun_pc3"),
        )
        and close(
            branch_summaries[name].get("captured_mass_msun"),
            values[-1].get("captured_mass_msun"),
        )
        and close(
            branch_summaries[name].get("n_conduction"),
            values[-1].get("n_conduction"),
        )
        for name, values in (("control", control), ("sink", sink))
    )
    current_gates = {
        "fixed_snapshot_closure_declared": (
            fixed_mode == "fixed_at_matched_snapshot"
        ),
        "measured_mass_current_applied": close(
            applied.get("Mdot_msun_per_myr"),
            closure.get("measured_mdot_msun_per_myr"),
        ),
        "measured_sink_luminosity_applied": close(
            applied.get("sink_luminosity_msun_kms2_per_myr"),
            closure.get("thermal_sink_current_msun_kms2_per_myr"),
        ),
        "sink_sign_is_physical": (
            float(applied.get("sink_luminosity_msun_kms2_per_myr", math.nan)) < 0.0
        ),
        "control_has_no_current": (
            close(applied.get("control_luminosity_msun_kms2_per_myr"), 0.0)
            and close(branch_summaries["control"].get("Mdot_msun_per_myr"), 0.0)
            and close(branch_summaries["control"].get(
                "L_inner_msun_kms2_per_myr"
            ), 0.0)
        ),
        "sink_branch_matches_applied_current": (
            close(branch_summaries["sink"].get("Mdot_msun_per_myr"),
                  applied.get("Mdot_msun_per_myr"))
            and close(branch_summaries["sink"].get(
                "L_inner_msun_kms2_per_myr"
            ), applied.get("sink_luminosity_msun_kms2_per_myr"))
        ),
        "fixed_boundary_currents_in_trajectories": fixed_trajectory_currents,
        "captured_mass_accounting_in_trajectories": captured_mass_accounting,
        "black_hole_mass_fixed_in_trajectories": black_hole_mass_fixed,
        "summary_trajectory_endpoints_match": summary_trajectory_endpoints_match,
    }
    budget_gates = {}
    for branch in ("control", "sink"):
        value = branch_summaries[branch]["conduction_energy_budget_relative_error"]
        budget_gates[branch] = (
            value is not None
            and math.isfinite(float(value))
            and 0.0 <= float(value) <= 1.0e-10
        )
    initial_match = all(
        close(control[0][field], sink[0][field], tolerance=5.0e-12)
        for field in (
            "rho_inner_mean_msun_pc3", "sigma_inner_1d_kms", "r_inner_pc",
            "inner_shell_mass_msun", "M_bh_msun",
            "lmfp_scaleheight_factor_inner", "lmfp_scaleheight_factor_min",
            "time_myr", "captured_mass_msun", "E_boundary_msun_kms2",
            "E_total_code", "n_conduction",
        )
    )
    bh_conductivity_active = all(
        0.0 < float(row["lmfp_scaleheight_factor_inner"]) < 1.0
        and 0.0 < float(row["lmfp_scaleheight_factor_min"]) < 1.0
        for values in (control, sink)
        for row in values
    )
    gates = {
        "accepted_closure": bool(summary.get("input_closure_accepted")),
        "resolved_inner_density_definition_verified": (
            resolved_diagnostics_verified
        ),
        "closure_scope_matches_summary": (
            bool(summary.get("absolute_two_current_closure")) == absolute_closure
            and bool(summary.get("mass_current_closure")) == mass_closure
        ) if ("absolute_two_current_closure" in summary
              or "mass_current_closure" in summary) else absolute_closure,
        "capture_binding_energy_excluded_for_mass_closure": (
            not mass_closure
            or (
                closure.get("capture_binding_energy_current_used_by_fluid") is False
                and summary.get("applied_currents", {}).get(
                    "source_luminosity_msun_kms2_per_myr"
                ) is None
            )
        ),
        **identity_gates,
        **current_gates,
        "common_evolved_interval": common_tau > 0.0,
        "identical_initial_fluid_state": initial_match,
        "black_hole_limited_conductivity_active": bh_conductivity_active,
        "control_energy_budget": budget_gates["control"],
        "sink_energy_budget": budget_gates["sink"],
        "black_hole_in_both_branches": all(
            math.isclose(
                float(branch_summaries[name]["black_hole_remap"]["M_bh_final_msun"]),
                float(summary["M_bh_msun"]),
                rel_tol=1.0e-10,
            )
            for name in ("control", "sink")
        ),
        "identical_black_hole_remap": (
            branch_summaries["control"].get("black_hole_remap")
            == branch_summaries["sink"].get("black_hole_remap")
        ),
        "capture_mass_validity_respected": (
            float(branch_summaries["sink"][
                "captured_mass_fraction_of_first_fluid_shell"
            ]) <= float(branch_summaries["sink"][
                "maximum_unmodeled_fluid_shell_mass_fraction"
            ]) * 1.01
            and float(branch_summaries["sink"][
                "captured_mass_fraction_of_GNC_owned_domain"
            ]) <= float(branch_summaries["sink"][
                "maximum_unmodeled_GNC_mass_fraction"
            ]) * 1.01
        ),
        "control_completed_physically": (
            branch_summaries["control"].get("stop_reason") in VALID_CONTROL_STOPS
        ),
        "sink_completed_physically": (
            branch_summaries["sink"].get("stop_reason") in VALID_SINK_STOPS
        ),
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(comparison_rows[0]))
        writer.writeheader()
        writer.writerows(comparison_rows)

    result = {
        "schema": "black-hole-aware-fluid-response-v3",
        "status": (
            "BLACK_HOLE_AWARE_RESPONSE_COMPLETE"
            if all(gates.values()) else "FAIL"
        ),
        "gates": gates,
        "common_tau_relax": common_tau,
        "control_stop_reason": branch_summaries["control"]["stop_reason"],
        "sink_stop_reason": branch_summaries["sink"]["stop_reason"],
        "inner_mean_density_sink_over_control_at_common_end": ratios[
            "rho_inner_mean_msun_pc3"
        ][-1],
        "inner_mean_density_max_abs_fractional_difference": max(
            abs(value - 1.0) for value in ratios["rho_inner_mean_msun_pc3"]
        ),
        "inner_dispersion_sink_over_control_at_common_end": ratios[
            "sigma_inner_1d_kms"
        ][-1],
        "inner_radius_sink_over_control_at_common_end": ratios["r_inner_pc"][-1],
        "resolved_density_diagnostic": summary["resolved_diagnostics"],
        "resolved_density_max_relative_consistency_error": (
            resolved_density_max_relative_error
        ),
        "resolved_density_consistency_tolerance": (
            args.mean_density_consistency_tolerance
        ),
        "interpretation": (
            "Both branches include the same central point mass and the same "
            "black-hole-limited conductivity. The closure current is held at "
            "its matched-snapshot value. The branch difference isolates the "
            "thermal energy removed by that current over the interval in "
            "which omitted captured mass remains below the interface validity "
            "thresholds. Density refers to the mean enclosed density of the "
            "innermost resolved Lagrangian shell, not an extrapolated central value."
        ),
        "closure_scope": (
            "absolute_two_current_closure"
            if absolute_closure else "mass_current_closure_only"
        ),
        "capture_binding_energy_current_used_by_fluid": False,
        "fluid_feedback_summary": args.summary.name,
        "fluid_feedback_trajectories": args.trajectories.name,
        "comparison_csv": args.out_csv.name,
        "identity": {
            "closure_json_sha256": sha256(args.closure_json),
            "fluid_profile_sha256": sha256(args.profile),
            "bridge_json_sha256": sha256(args.bridge_json),
            "fluid_feedback_summary_sha256": sha256(args.summary),
            "fluid_feedback_trajectories_sha256": sha256(args.trajectories),
            "fluid_comparison_csv_sha256": sha256(args.out_csv),
        },
    }
    args.out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "BLACK_HOLE_AWARE_RESPONSE_COMPLETE" else 2


if __name__ == "__main__":
    raise SystemExit(main())
