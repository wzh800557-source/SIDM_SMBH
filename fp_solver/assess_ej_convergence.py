#!/usr/bin/env python3
"""Assess production convergence of the depleted non-local ``(E,J)`` solve.

The assessment keeps Monte Carlo sampling, energy resolution, angular
resolution, and reporting-surface placement separate. Independent operators
are solved individually and after zero-filled combination. A production pass
requires two independent seeds per represented configuration, converged mass
currents, and, when requested, converged capture-binding-energy currents.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


MODELS = ("direct", "immediate")
METRICS = {
    "mass": ("steady_capture_mdot_msun_per_myr", "Msun/Myr"),
    "capture_binding_energy": (
        "steady_capture_binding_energy_current_msun_kms2_per_myr",
        "Msun (km/s)^2/Myr",
    ),
}


def symmetric_fraction(a: float, b: float) -> float:
    return abs(a - b) / max(0.5 * (abs(a) + abs(b)), 1.0e-300)


def load_case(path: Path) -> dict:
    data = json.loads(path.read_text())
    if data.get("schema") != "finite-angle-ej-steady-state-v2":
        raise RuntimeError(f"{path} is not an EJ steady-state v2 result")
    metadata = data["operator_metadata"]
    return {
        "path": path.name,
        "status": data["status"],
        "energy_bins": int(metadata["energy_bins"]),
        "angular_bins": int(metadata["angular_bins"]),
        "angular_grid": metadata.get("angular_grid", "unspecified"),
        "reservoir_radius_pc": float(metadata["reservoir_radius_pc"]),
        "profile_sha256": metadata["profile_sha256"],
        "df_table_sha256": metadata["df_table_sha256"],
        "mbh_msun": float(metadata["mbh_msun"]),
        "sigma0_over_m_cm2_g": float(metadata["sigma0_over_m_cm2_g"]),
        "w_kms": float(metadata["w_kms"]),
        "combined": "seeds" in metadata,
        "seed": metadata.get("seed"),
        "seeds": metadata.get("seeds"),
        "models": data["models"],
    }


def group_key(case: dict) -> tuple[int, int, float]:
    return (
        case["energy_bins"],
        case["angular_bins"],
        round(case["reservoir_radius_pc"], 12),
    )


def metric_value(case: dict, model: str, metric: str) -> float | None:
    key = METRICS[metric][0]
    value = case["models"][model].get(key)
    return float(value) if value is not None else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", nargs="+", type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    parser.add_argument("--out-csv", required=True, type=Path)
    parser.add_argument("--fractional-tolerance", type=float, default=0.10)
    parser.add_argument("--minimum-seeds", type=int, default=2)
    parser.add_argument("--minimum-boundaries", type=int, default=3)
    parser.add_argument(
        "--primary-boundary-radius",
        type=float,
        default=None,
        help=(
            "physical local-circular handoff radius. The highest-resolution "
            "energy and angular ladders are selected at this surface."
        ),
    )
    parser.add_argument(
        "--require-capture-energy",
        action="store_true",
        help="also require the capture-binding-energy current to converge",
    )
    parser.add_argument(
        "--require-dimensions",
        nargs="+",
        choices=("energy", "angular", "boundary"),
        default=("energy", "angular", "boundary"),
    )
    args = parser.parse_args()
    if not (0.0 < args.fractional_tolerance < 1.0):
        raise ValueError("fractional tolerance must lie in (0,1)")
    if args.minimum_seeds < 2 or args.minimum_boundaries < 2:
        raise ValueError("production convergence requires multiple realizations")

    cases = [load_case(path) for path in args.results]
    invariant_keys = (
        "profile_sha256", "df_table_sha256", "mbh_msun",
        "sigma0_over_m_cm2_g", "w_kms",
    )
    reference = cases[0]
    for case in cases[1:]:
        for key in invariant_keys:
            if case[key] != reference[key]:
                raise RuntimeError(f"steady-state cases differ in {key}")

    required_metrics = ["mass"]
    if args.require_capture_energy:
        required_metrics.append("capture_binding_energy")

    grouped: dict[tuple[int, int, float], list[dict]] = defaultdict(list)
    for case in cases:
        grouped[group_key(case)].append(case)

    rows = []
    group_reports = []
    seed_gates = []
    numerical_gates = []
    combined_map = {}
    for key in sorted(grouped):
        energy_bins, angular_bins, radius = key
        items = grouped[key]
        combined = [item for item in items if item["combined"]]
        seeds = [item for item in items if not item["combined"]]
        if len(combined) != 1:
            raise RuntimeError(
                f"configuration {key} needs exactly one combined result"
            )
        if len(seeds) < args.minimum_seeds:
            raise RuntimeError(
                f"configuration {key} has only {len(seeds)} independent seeds"
            )
        seed_ids = [item["seed"] for item in seeds]
        if any(seed is None for seed in seed_ids):
            raise RuntimeError(f"configuration {key} has an unlabelled realization")
        if len(set(seed_ids)) != len(seed_ids):
            raise RuntimeError(f"configuration {key} repeats an independent seed")
        combined_seed_ids = combined[0].get("seeds")
        if (
            not isinstance(combined_seed_ids, list)
            or sorted(int(seed) for seed in combined_seed_ids)
            != sorted(int(seed) for seed in seed_ids)
        ):
            raise RuntimeError(
                f"configuration {key} combined operator does not contain its seeds"
            )
        angular_grids = {item["angular_grid"] for item in items}
        if len(angular_grids) != 1:
            raise RuntimeError(f"configuration {key} mixes angular-grid definitions")
        combined_map[key] = combined[0]
        report = {
            "energy_bins": energy_bins,
            "angular_bins": angular_bins,
            "angular_grid": combined[0]["angular_grid"],
            "reservoir_radius_pc": radius,
            "combined_result": combined[0]["path"],
            "seed_results": [item["path"] for item in seeds],
            "independent_seeds": sorted(int(seed) for seed in seed_ids),
            "models": {},
        }
        numerical_pass = all(
            item["status"] == "EJ_STEADY_SMOKE_PASS" for item in items
        )
        report["numerical_solver_gate"] = numerical_pass
        numerical_gates.append(numerical_pass)
        for model in MODELS:
            report["models"][model] = {}
            for metric in required_metrics:
                combined_value = metric_value(combined[0], model, metric)
                seed_values = [metric_value(item, model, metric) for item in seeds]
                available = (
                    combined_value is not None
                    and all(value is not None for value in seed_values)
                )
                if not available:
                    metric_report = {
                        "available": False,
                        "combined_value": None,
                        "seed_values": seed_values,
                        "seed_full_range_fraction_of_mean": None,
                        "seed_gate": False,
                    }
                    seed_gate = False
                else:
                    values = [float(value) for value in seed_values]
                    seed_range = (
                        max(values) - min(values)
                    ) / max(sum(values) / len(values), 1.0e-300)
                    seed_gate = seed_range <= args.fractional_tolerance
                    metric_report = {
                        "available": True,
                        "combined_value": combined_value,
                        "seed_values": values,
                        "seed_full_range_fraction_of_mean": seed_range,
                        "seed_gate": seed_gate,
                    }
                seed_gates.append(seed_gate)
                report["models"][model][metric] = metric_report
                rows.append({
                    "energy_bins": energy_bins,
                    "angular_bins": angular_bins,
                    "reservoir_radius_pc": radius,
                    "model": model,
                    "metric": metric,
                    "unit": METRICS[metric][1],
                    "combined_value": combined_value,
                    "seed_count": len(seed_values),
                    "seed_full_range_fraction_of_mean": metric_report[
                        "seed_full_range_fraction_of_mean"
                    ],
                    "seed_gate": seed_gate,
                    "numerical_solver_gate": numerical_pass,
                })
        group_reports.append(report)

    comparisons = {"energy": [], "angular": [], "boundary": []}
    for model in MODELS:
        for metric in required_metrics:
            by_energy = defaultdict(list)
            by_angular = defaultdict(list)
            by_boundary = defaultdict(list)
            for (energy_bins, angular_bins, radius), case in combined_map.items():
                value = metric_value(case, model, metric)
                if value is None:
                    continue
                by_energy[(angular_bins, radius)].append((energy_bins, value))
                by_angular[(energy_bins, radius)].append((angular_bins, value))
                by_boundary[(energy_bins, angular_bins)].append((radius, value))

            for fixed, values in by_energy.items():
                values.sort()
                if len(values) >= 2:
                    lo, hi = values[-2], values[-1]
                    difference = symmetric_fraction(lo[1], hi[1])
                    comparisons["energy"].append({
                        "model": model,
                        "metric": metric,
                        "fixed_angular_bins": fixed[0],
                        "fixed_reservoir_radius_pc": fixed[1],
                        "lower_energy_bins": lo[0],
                        "upper_energy_bins": hi[0],
                        "lower_value": lo[1],
                        "upper_value": hi[1],
                        "fractional_difference": difference,
                        "gate": difference <= args.fractional_tolerance,
                    })
            for fixed, values in by_angular.items():
                values.sort()
                if len(values) >= 2:
                    lo, hi = values[-2], values[-1]
                    difference = symmetric_fraction(lo[1], hi[1])
                    comparisons["angular"].append({
                        "model": model,
                        "metric": metric,
                        "fixed_energy_bins": fixed[0],
                        "fixed_reservoir_radius_pc": fixed[1],
                        "lower_angular_bins": lo[0],
                        "upper_angular_bins": hi[0],
                        "lower_value": lo[1],
                        "upper_value": hi[1],
                        "fractional_difference": difference,
                        "gate": difference <= args.fractional_tolerance,
                    })
            for fixed, values in by_boundary.items():
                values.sort()
                if len(values) >= args.minimum_boundaries:
                    currents = [value for _, value in values]
                    difference = (
                        max(currents) - min(currents)
                    ) / max(sum(currents) / len(currents), 1.0e-300)
                    comparisons["boundary"].append({
                        "model": model,
                        "metric": metric,
                        "fixed_energy_bins": fixed[0],
                        "fixed_angular_bins": fixed[1],
                        "radii_pc": [radius for radius, _ in values],
                        "values": currents,
                        "full_range_fraction_of_mean": difference,
                        "gate": difference <= args.fractional_tolerance,
                    })

    radii = sorted({key[2] for key in combined_map})
    requested_primary = (
        float(args.primary_boundary_radius)
        if args.primary_boundary_radius is not None
        else radii[len(radii) // 2]
    )
    primary_radius = min(radii, key=lambda value: abs(value - requested_primary))
    if symmetric_fraction(primary_radius, requested_primary) > 1.0e-8:
        raise RuntimeError("the requested primary boundary is absent from the scan")
    primary_keys = [key for key in combined_map if key[2] == primary_radius]
    final_energy = max(key[0] for key in primary_keys)
    final_angular = max(
        key[1] for key in primary_keys if key[0] == final_energy
    )
    energy_ladders = defaultdict(list)
    for energy_bins, angular_bins, radius in primary_keys:
        energy_ladders[angular_bins].append(energy_bins)
    energy_ladder_candidates = [
        angular_bins for angular_bins, energies in energy_ladders.items()
        if final_energy in energies and len(set(energies)) >= 2
    ]
    energy_ladder_angular = (
        max(energy_ladder_candidates) if energy_ladder_candidates else None
    )
    angular_ladder = sorted({
        key[1] for key in primary_keys if key[0] == final_energy
    })
    lower_energy = (
        sorted(set(energy_ladders[energy_ladder_angular]))[-2]
        if energy_ladder_angular is not None else None
    )
    lower_angular = angular_ladder[-2] if len(angular_ladder) >= 2 else None

    selected_comparisons = {dimension: [] for dimension in comparisons}
    if energy_ladder_angular is not None:
        selected_comparisons["energy"] = [
            item for item in comparisons["energy"]
            if item["fixed_angular_bins"] == energy_ladder_angular
            and math.isclose(
                item["fixed_reservoir_radius_pc"], primary_radius,
                rel_tol=1.0e-10, abs_tol=0.0,
            )
            and item["upper_energy_bins"] == final_energy
        ]
    if lower_angular is not None:
        selected_comparisons["angular"] = [
            item for item in comparisons["angular"]
            if item["fixed_energy_bins"] == final_energy
            and math.isclose(
                item["fixed_reservoir_radius_pc"], primary_radius,
                rel_tol=1.0e-10, abs_tol=0.0,
            )
            and item["upper_angular_bins"] == final_angular
        ]
    selected_comparisons["boundary"] = [
        item for item in comparisons["boundary"]
        if item["fixed_energy_bins"] == final_energy
        and item["fixed_angular_bins"] == final_angular
    ]
    expected_selected = len(MODELS) * len(required_metrics)
    dimension_gates = {
        dimension: {
            "present": len(selected_comparisons[dimension]) == expected_selected,
            "comparison_count": len(selected_comparisons[dimension]),
            "expected_comparison_count": expected_selected,
            "pass": (
                len(selected_comparisons[dimension]) == expected_selected
                and all(item["gate"] for item in selected_comparisons[dimension])
            ),
        }
        for dimension in comparisons
    }
    required_dimensions_pass = all(
        dimension_gates[dimension]["pass"]
        for dimension in args.require_dimensions
    )
    all_pass = bool(
        numerical_gates and all(numerical_gates)
        and seed_gates and all(seed_gates)
        and required_dimensions_pass
    )
    result = {
        "schema": "finite-angle-ej-convergence-v3",
        "status": (
            "PRODUCTION_CONVERGENCE_PASS" if all_pass
            else "INCOMPLETE_OR_FAILED"
        ),
        "fractional_tolerance": args.fractional_tolerance,
        "minimum_independent_seeds": args.minimum_seeds,
        "minimum_reporting_surfaces": args.minimum_boundaries,
        "required_dimensions": list(args.require_dimensions),
        "required_metrics": required_metrics,
        "identity": {key: reference[key] for key in invariant_keys},
        "dimension_gates": dimension_gates,
        "production_resolution_path": {
            "primary_boundary_radius_pc": primary_radius,
            "primary_boundary_source": (
                "command-line" if args.primary_boundary_radius is not None
                else "median represented reporting radius"
            ),
            "final_energy_bins": final_energy,
            "final_angular_bins": final_angular,
            "energy_ladder_angular_bins": energy_ladder_angular,
            "previous_energy_bins": lower_energy,
            "previous_angular_bins": lower_angular,
        },
        "groups": group_reports,
        "comparisons": comparisons,
        "selected_comparisons": selected_comparisons,
        "note": (
            "The pass applies to the production collapse snapshot and the "
            "sampled finite-angle operator. It is not a statement that the "
            "capture binding energy is returned to the halo."
        ),
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    with args.out_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if all_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
