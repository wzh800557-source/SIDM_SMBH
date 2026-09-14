#!/usr/bin/env python3
"""Gate a stationary FP energy current measured from direct plunge records.

The gate accepts only particles whose inherited initial binding energy lies
outside the FP reporting surface.  It therefore measures a boundary-fed
current rather than depletion of the initially populated cusp.  Independent
random seeds, two numerical resolutions, stationary inner inventories, and the
exact binding-current identity are all required before fluid coupling.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np


TINY = np.finfo(float).tiny


def relative_change(a: float, b: float) -> float:
    return abs(float(a) - float(b)) / max(abs(float(a)), abs(float(b)), TINY)


def fractional_range(values: Iterable[float]) -> float:
    values = np.asarray(list(values), dtype=float)
    if values.size < 2 or np.any(~np.isfinite(values)):
        return math.inf
    return float(np.ptp(values) / max(abs(float(np.mean(values))), TINY))


def _finite_float(mapping: dict, key: str) -> float:
    value = float(mapping[key])
    if not math.isfinite(value):
        raise ValueError(f"{key} is not finite")
    return value


def extract(path: Path) -> dict:
    value = json.loads(path.read_text())
    if value.get("schema") != "gnc-direct-energy-current-diagnostic-v3":
        raise ValueError(f"{path}: unsupported diagnostic schema")
    current = value["boundary_supplied_current"]
    inventory = value["inner_reservoir_stationarity"]
    provenance = value["run_provenance"]
    sufficient = current["sufficient_statistics"]
    seed = provenance.get("gnc_seed")
    gx = provenance.get("gx_bins")
    dc = provenance.get("dc_bins")
    if seed is None or gx is None or dc is None:
        raise ValueError(f"{path}: seed or grid provenance is missing")
    return {
        "path": str(path),
        "status": value.get("status"),
        "seed": int(seed),
        "grid": (int(gx), int(dc)),
        "physical_fingerprint": provenance.get(
            "physical_configuration_fingerprint"
        ),
        "reservoir_fingerprint": provenance.get("common_reservoir_fingerprint"),
        "kernel_compatible": value.get("physical_kernel_compatible") is True,
        "mass_consistency": _finite_float(
            value, "event_vs_plunge_mass_relative_difference"
        ),
        "raw_count": int(current["raw_capture_count"]),
        "effective_count": _finite_float(current, "effective_capture_count"),
        "mdot": _finite_float(current, "mdot_msun_per_myr"),
        "capture_current": _finite_float(
            current, "capture_binding_current_msun_kms2_per_myr"
        ),
        "advective_current": _finite_float(
            current, "boundary_advected_binding_current_msun_kms2_per_myr"
        ),
        "outward_current": _finite_float(
            current, "outward_released_binding_current_msun_kms2_per_myr"
        ),
        "identity_relative": _finite_float(
            current, "current_identity_relative_residual"
        ),
        "nonnegative_release": current.get("nonnegative_released_energy_pass") is True,
        "mass_plateau_change": _finite_float(
            current["mass_current_plateau"], "relative_change"
        ),
        "energy_plateau_change": _finite_float(
            current["released_energy_current_plateau"], "relative_change"
        ),
        "inventory_completeness": _finite_float(
            inventory, "completeness_fraction"
        ),
        "inventory_mass_change": _finite_float(
            inventory["mass_inventory_plateau"], "relative_change"
        ),
        "inventory_energy_change": _finite_float(
            inventory["binding_energy_inventory_plateau"], "relative_change"
        ),
        "mass_weight_sum": _finite_float(
            sufficient, "physical_mass_weight_sum_msun"
        ),
        "mass_weight_squared_sum": _finite_float(
            sufficient, "physical_mass_weight_squared_sum_msun2"
        ),
        "weighted_x_capture_sum": _finite_float(
            sufficient, "weighted_x_capture_sum_msun"
        ),
        "weighted_delta_x_sum": _finite_float(
            sufficient, "weighted_delta_x_sum_msun"
        ),
        "time_myr": _finite_float(current, "time_myr"),
        "energy_unit_kms2": _finite_float(value, "energy_unit_kms2"),
        "energy_flux_scale": _finite_float(
            value, "energy_flux_scale_msun_kms2_per_myr"
        ),
        "mass_flux_scale": _finite_float(value, "mass_flux_scale_msun_per_myr"),
        "x_boundary": _finite_float(current, "boundary_x"),
        "tnr_myr": _finite_float(value, "tnr_myr_used_for_physical_flux"),
        "burn_in_tnr": _finite_float(value, "burn_in_tnr"),
        "window_tnr": _finite_float(value, "window_tnr"),
    }


def evaluate(
    diagnostics: Iterable[dict],
    *,
    minimum_raw_captures: int = 50,
    minimum_effective_captures: float = 50.0,
    plateau_tolerance: float = 0.25,
    inventory_tolerance: float = 0.25,
    inventory_completeness: float = 0.95,
    seed_tolerance: float = 0.25,
    grid_tolerance: float = 0.10,
    mass_consistency_tolerance: float = 1.0e-3,
    identity_tolerance: float = 1.0e-10,
) -> dict:
    rows = list(diagnostics)
    if not rows:
        raise ValueError("at least one direct FP diagnostic is required")
    for row in rows:
        if row["mass_weight_sum"] < 0.0 or row["mass_weight_squared_sum"] < 0.0:
            raise ValueError("capture sufficient statistics must be non-negative")

    physical_fingerprints = {row["physical_fingerprint"] for row in rows}
    grids = sorted({row["grid"] for row in rows}, key=lambda grid: (grid[0] * grid[1], grid))
    grouped = {
        grid: [row for row in rows if row["grid"] == grid]
        for grid in grids
    }
    duplicate_seed_grid = len({(row["grid"], row["seed"]) for row in rows}) != len(rows)
    seed_design = bool(grouped) and all(
        len({row["seed"] for row in subset}) >= 2 for subset in grouped.values()
    ) and not duplicate_seed_grid
    resolution_design = len(grids) >= 2
    reservoir_by_grid = {
        f"gx{grid[0]}_dc{grid[1]}": {
            "fingerprints": sorted(
                str(value) for value in {row["reservoir_fingerprint"] for row in subset}
            ),
            "consistent": (
                len({row["reservoir_fingerprint"] for row in subset}) == 1
                and None not in {row["reservoir_fingerprint"] for row in subset}
            ),
        }
        for grid, subset in grouped.items()
    }

    seed_convergence = {}
    for grid, subset in grouped.items():
        key = f"gx{grid[0]}_dc{grid[1]}"
        seed_convergence[key] = {
            "seeds": sorted(row["seed"] for row in subset),
            "mdot_fractional_range": fractional_range(row["mdot"] for row in subset),
            "outward_energy_current_fractional_range": fractional_range(
                row["outward_current"] for row in subset
            ),
            "mean_mdot_msun_per_myr": float(np.mean([row["mdot"] for row in subset])),
            "mean_outward_energy_current_msun_kms2_per_myr": float(
                np.mean([row["outward_current"] for row in subset])
            ),
        }

    resolution_convergence = {
        "available": False,
        "grids": [],
        "mdot_relative_change": math.inf,
        "outward_energy_current_relative_change": math.inf,
    }
    if len(grids) >= 2:
        low, high = grids[-2:]
        low_key = f"gx{low[0]}_dc{low[1]}"
        high_key = f"gx{high[0]}_dc{high[1]}"
        resolution_convergence = {
            "available": True,
            "grids": [list(low), list(high)],
            "mdot_relative_change": relative_change(
                seed_convergence[low_key]["mean_mdot_msun_per_myr"],
                seed_convergence[high_key]["mean_mdot_msun_per_myr"],
            ),
            "outward_energy_current_relative_change": relative_change(
                seed_convergence[low_key][
                    "mean_outward_energy_current_msun_kms2_per_myr"
                ],
                seed_convergence[high_key][
                    "mean_outward_energy_current_msun_kms2_per_myr"
                ],
            ),
        }

    # Only independent seeds at the highest resolution set the accepted value.
    # Lower resolutions are convergence tests, not additional physical exposure.
    highest_grid = grids[-1]
    accepted_rows = grouped[highest_grid]
    sum_w = float(sum(row["mass_weight_sum"] for row in accepted_rows))
    sum_w2 = float(sum(row["mass_weight_squared_sum"] for row in accepted_rows))
    pooled_effective = sum_w * sum_w / sum_w2 if sum_w2 > 0.0 else 0.0
    pooled_raw = int(sum(row["raw_count"] for row in accepted_rows))
    pooled_time = float(sum(row["time_myr"] for row in accepted_rows))
    pooled_mdot = sum_w / pooled_time if pooled_time > 0.0 else math.nan
    mean_x_capture = (
        sum(row["weighted_x_capture_sum"] for row in accepted_rows) / sum_w
        if sum_w > 0.0 else math.nan
    )
    mean_delta_x = (
        sum(row["weighted_delta_x_sum"] for row in accepted_rows) / sum_w
        if sum_w > 0.0 else math.nan
    )
    reference = accepted_rows[0]
    energy_unit = reference["energy_unit_kms2"]
    x_boundary = reference["x_boundary"]
    capture_current = pooled_mdot * mean_x_capture * energy_unit
    advective_current = pooled_mdot * x_boundary * energy_unit
    outward_current = pooled_mdot * mean_delta_x * energy_unit
    identity_residual = capture_current - advective_current - outward_current
    identity_relative = abs(identity_residual) / max(
        abs(capture_current), abs(advective_current), abs(outward_current), TINY
    )

    exact_configuration_keys = (
        "energy_unit_kms2", "energy_flux_scale", "mass_flux_scale", "x_boundary",
        "tnr_myr", "burn_in_tnr", "window_tnr",
    )
    scalar_configuration_consistent = all(
        math.isclose(
            row[key], reference[key], rel_tol=1.0e-12, abs_tol=0.0
        )
        for row in rows for key in exact_configuration_keys
    )
    seed_converged = bool(seed_convergence) and all(
        value["mdot_fractional_range"] <= seed_tolerance
        and value["outward_energy_current_fractional_range"] <= seed_tolerance
        for value in seed_convergence.values()
    )
    grid_converged = (
        resolution_convergence["available"]
        and resolution_convergence["mdot_relative_change"] <= grid_tolerance
        and resolution_convergence["outward_energy_current_relative_change"]
        <= grid_tolerance
    )
    gates = {
        "diagnostics_complete": all(
            row["status"] == "DIRECT_FP_ENERGY_DIAGNOSTIC_COMPLETE" for row in rows
        ),
        "physical_configuration_consistent": (
            len(physical_fingerprints) == 1
            and None not in physical_fingerprints
            and scalar_configuration_consistent
        ),
        "independent_seed_design": seed_design,
        "resolution_design": resolution_design,
        "common_reservoir_within_each_grid": all(
            value["consistent"] for value in reservoir_by_grid.values()
        ),
        "physical_kernel_compatible": all(row["kernel_compatible"] for row in rows),
        "event_plunge_mass_consistent": all(
            row["mass_consistency"] <= mass_consistency_tolerance for row in rows
        ),
        "inner_mass_energy_inventory_stationary": all(
            row["inventory_completeness"] >= inventory_completeness
            and row["inventory_mass_change"] <= inventory_tolerance
            and row["inventory_energy_change"] <= inventory_tolerance
            for row in rows
        ),
        "boundary_fed_currents_stationary": all(
            row["mass_plateau_change"] <= plateau_tolerance
            and row["energy_plateau_change"] <= plateau_tolerance
            for row in rows
        ),
        "pooled_capture_statistics": (
            pooled_raw >= minimum_raw_captures
            and pooled_effective >= minimum_effective_captures
        ),
        "independent_seed_convergence": seed_converged,
        "resolution_convergence": grid_converged,
        "released_binding_energy_nonnegative": (
            all(row["nonnegative_release"] for row in rows)
            and math.isfinite(outward_current)
            and outward_current >= 0.0
        ),
        "binding_current_identity": (
            all(row["identity_relative"] <= identity_tolerance for row in rows)
            and identity_relative <= identity_tolerance
        ),
    }
    authorized = all(gates.values())
    candidate = {
        "selected_grid": list(highest_grid),
        "independent_seeds": sorted(row["seed"] for row in accepted_rows),
        "raw_boundary_supplied_captures": pooled_raw,
        "effective_boundary_supplied_captures": pooled_effective,
        "combined_exposure_myr": pooled_time,
        "capture_mdot_msun_per_myr": pooled_mdot,
        "interface_specific_orbital_energy_kms2": -x_boundary * energy_unit,
        "captured_specific_orbital_energy_kms2": -mean_x_capture * energy_unit,
        "mean_released_binding_energy_kms2": mean_delta_x * energy_unit,
        "capture_binding_current_msun_kms2_per_myr": capture_current,
        "boundary_advected_binding_current_msun_kms2_per_myr": advective_current,
        "outward_luminosity_msun_kms2_per_myr": outward_current,
        "binding_current_identity_relative_residual": identity_relative,
        "C_M": pooled_mdot / reference["mass_flux_scale"],
        "C_E_capture_binding": capture_current / reference["energy_flux_scale"],
        "C_E_boundary_advected_binding": advective_current / reference["energy_flux_scale"],
        "C_E_outward_released_binding": outward_current / reference["energy_flux_scale"],
    }
    return {
        "schema": "gnc-capture-energy-current-gate-v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "ENERGY_CURRENT_GATE_PASS" if authorized else "ENERGY_CURRENT_GATE_FAIL",
        "fluid_coupling_authorized": authorized,
        "measurement": (
            "Direct FP plunge records restricted to boundary-supplied trajectories, "
            "with live inner-reservoir stationarity."
        ),
        "thresholds": {
            "minimum_raw_boundary_supplied_captures": minimum_raw_captures,
            "minimum_effective_boundary_supplied_captures": minimum_effective_captures,
            "mass_and_energy_plateau_relative": plateau_tolerance,
            "inner_inventory_relative": inventory_tolerance,
            "inner_inventory_completeness": inventory_completeness,
            "independent_seed_relative": seed_tolerance,
            "resolution_relative": grid_tolerance,
            "event_plunge_mass_relative": mass_consistency_tolerance,
            "binding_current_identity_relative": identity_tolerance,
        },
        "gates": gates,
        "failed_gates": [key for key, passed in gates.items() if not passed],
        "diagnostics_evaluated": len(rows),
        "grids": [list(grid) for grid in grids],
        "reservoir_consistency": reservoir_by_grid,
        "seed_convergence": seed_convergence,
        "resolution_convergence": resolution_convergence,
        "candidate_closure": candidate,
        "accepted_closure": candidate if authorized else None,
        "decision": (
            "The measured FP energy current is authorized for conservative fluid coupling."
            if authorized else
            "The conservative fluid response remains locked until every direct FP gate passes."
        ),
    }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnostic", type=Path, action="append", default=[])
    parser.add_argument(
        "--run-root",
        type=Path,
        default=None,
        help="also discover */absolute_closure_diagnostics.json below this directory",
    )
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--minimum-raw-captures", type=int, default=50)
    parser.add_argument("--minimum-effective-captures", type=float, default=50.0)
    parser.add_argument("--plateau-tolerance", type=float, default=0.25)
    parser.add_argument("--inventory-tolerance", type=float, default=0.25)
    parser.add_argument("--inventory-completeness", type=float, default=0.95)
    parser.add_argument("--seed-tolerance", type=float, default=0.25)
    parser.add_argument("--grid-tolerance", type=float, default=0.10)
    parser.add_argument("--mass-consistency-tolerance", type=float, default=1.0e-3)
    parser.add_argument("--identity-tolerance", type=float, default=1.0e-10)
    args = parser.parse_args(argv)
    paths = list(args.diagnostic)
    if args.run_root is not None:
        paths.extend(sorted(args.run_root.rglob("absolute_closure_diagnostics.json")))
    paths = list(dict.fromkeys(path.resolve() for path in paths))
    if not paths:
        parser.error("provide --diagnostic or --run-root")
    result = evaluate(
        [extract(path) for path in paths],
        minimum_raw_captures=args.minimum_raw_captures,
        minimum_effective_captures=args.minimum_effective_captures,
        plateau_tolerance=args.plateau_tolerance,
        inventory_tolerance=args.inventory_tolerance,
        inventory_completeness=args.inventory_completeness,
        seed_tolerance=args.seed_tolerance,
        grid_tolerance=args.grid_tolerance,
        mass_consistency_tolerance=args.mass_consistency_tolerance,
        identity_tolerance=args.identity_tolerance,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["fluid_coupling_authorized"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
