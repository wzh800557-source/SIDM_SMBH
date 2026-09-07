#!/usr/bin/env python3
"""Recover energy-resolved capture spectra from completed steady EJ solves.

This diagnostic reuses an existing operator and its saved occupation table.  It
does not rerun the nonlinear steady solver and cannot turn a failed convergence
test into a pass.  Its purpose is to locate the source-energy cells responsible
for a changing capture-binding-energy moment before another grid is selected.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

from solve_ej_steady_state import (
    capture_spectra,
    load_operator,
    population_rate,
    sha256_file,
)


MODELS = ("direct", "immediate")


def symmetric_fraction(a: float, b: float) -> float:
    return abs(a - b) / max(0.5 * (abs(a) + abs(b)), 1.0e-300)


def read_occupations(path: Path, state_count: int) -> dict[str, np.ndarray]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != state_count:
        raise RuntimeError(
            f"{path} contains {len(rows)} states; expected {state_count}"
        )
    states = np.asarray([int(row["state"]) for row in rows], dtype=np.int64)
    if not np.array_equal(states, np.arange(state_count, dtype=np.int64)):
        raise RuntimeError(f"{path} does not contain each state exactly once")
    occupations = {}
    for model in MODELS:
        key = f"occupation_{model}"
        if any(key not in row for row in rows):
            raise RuntimeError(f"{path} is missing {key}")
        value = np.asarray([float(row[key]) for row in rows])
        if np.any(~np.isfinite(value)) or np.any(value < 0.0):
            raise RuntimeError(f"{path} contains invalid {key} values")
        occupations[model] = value
    return occupations


def tail_fraction(rows: list[dict], key: str, first_bin: int) -> float | None:
    values = [row[key] for row in rows[first_bin:]]
    if any(value is None for value in values):
        return None
    return float(sum(float(value) for value in values))


def analyze_case(
    label: str,
    operator_path: Path,
    steady_json_path: Path,
    steady_csv_path: Path,
) -> dict:
    arrays, metadata = load_operator(operator_path)
    steady = json.loads(steady_json_path.read_text())
    if steady.get("schema") != "finite-angle-ej-steady-state-v2":
        raise RuntimeError(f"{steady_json_path} is not a steady EJ result")
    if steady.get("status") != "EJ_STEADY_SMOKE_PASS":
        raise RuntimeError(f"{steady_json_path} did not pass its steady solve")
    if steady.get("operator_sha256") != sha256_file(operator_path):
        raise RuntimeError(f"{label} operator hash does not match its steady result")
    if steady.get("operator_metadata") != metadata:
        raise RuntimeError(f"{label} operator metadata changed after the solve")
    state_count = int(metadata["state_count"])
    occupations = read_occupations(steady_csv_path, state_count)
    models = {}
    gates = {}
    for model in MODELS:
        rate = np.asarray(arrays[f"rate_{model}_msun_per_myr"], float)
        _, flux = population_rate(
            occupations[model],
            np.asarray(arrays["source1"], np.int64),
            np.asarray(arrays["source2"], np.int64),
            np.asarray(arrays["pre"], np.int64),
            np.asarray(arrays["post"], np.int64),
            rate,
        )
        spectrum = capture_spectra(arrays, occupations[model], model, flux)
        expected_mass = float(
            steady["models"][model]["steady_capture_mdot_msun_per_myr"]
        )
        expected_energy = float(
            steady["models"][model][
                "steady_capture_binding_energy_current_msun_kms2_per_myr"
            ]
        )
        measured_mass = float(spectrum["mass_current_roundtrip_msun_per_myr"])
        measured_energy = float(
            spectrum["binding_energy_current_roundtrip_msun_kms2_per_myr"]
        )
        mass_error = abs(measured_mass - expected_mass) / max(
            abs(expected_mass), 1.0e-300
        )
        energy_error = abs(measured_energy - expected_energy) / max(
            abs(expected_energy), 1.0e-300
        )
        rows = spectrum["by_source_energy"]
        n_energy = len(rows)
        highest_decile = min(
            max(int(math.floor(0.9 * n_energy)), 0), n_energy - 1
        )
        highest_quarter = min(
            max(int(math.floor(0.75 * n_energy)), 0), n_energy - 1
        )
        model_gates = {
            "mass_roundtrip_le_1e_10": mass_error <= 1.0e-10,
            "binding_energy_roundtrip_le_1e_10": energy_error <= 1.0e-10,
        }
        gates[model] = model_gates
        models[model] = {
            "steady_capture_mdot_msun_per_myr": expected_mass,
            "steady_capture_binding_energy_current_msun_kms2_per_myr": (
                expected_energy
            ),
            "mass_roundtrip_relative_error": mass_error,
            "binding_energy_roundtrip_relative_error": energy_error,
            "highest_energy_bin_mass_fraction": spectrum[
                "highest_energy_bin_mass_fraction"
            ],
            "highest_energy_bin_binding_energy_fraction": spectrum[
                "highest_energy_bin_binding_energy_fraction"
            ],
            "highest_decile_mass_fraction": tail_fraction(
                rows, "capture_mass_fraction", highest_decile
            ),
            "highest_decile_binding_energy_fraction": tail_fraction(
                rows, "capture_binding_energy_fraction", highest_decile
            ),
            "highest_quarter_mass_fraction": tail_fraction(
                rows, "capture_mass_fraction", highest_quarter
            ),
            "highest_quarter_binding_energy_fraction": tail_fraction(
                rows, "capture_binding_energy_fraction", highest_quarter
            ),
            "capture_spectra": spectrum,
        }
    return {
        "label": label,
        "operator": operator_path.name,
        "operator_sha256": sha256_file(operator_path),
        "steady_json": steady_json_path.name,
        "steady_json_sha256": sha256_file(steady_json_path),
        "steady_csv": steady_csv_path.name,
        "steady_csv_sha256": sha256_file(steady_csv_path),
        "energy_bins": int(metadata["energy_bins"]),
        "angular_bins": int(metadata["angular_bins"]),
        "reservoir_radius_pc": float(metadata["reservoir_radius_pc"]),
        "models": models,
        "gates": gates,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        action="append",
        nargs=4,
        metavar=("LABEL", "OPERATOR", "STEADY_JSON", "STEADY_CSV"),
        required=True,
        help="repeat for each completed energy grid",
    )
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    cases = [
        analyze_case(label, Path(operator), Path(steady_json), Path(steady_csv))
        for label, operator, steady_json, steady_csv in args.case
    ]
    labels = [case["label"] for case in cases]
    if len(set(labels)) != len(labels):
        raise RuntimeError("capture-tail case labels must be unique")

    comparisons = []
    grouped = defaultdict(list)
    for case in cases:
        grouped[(
            case["angular_bins"], round(case["reservoir_radius_pc"], 12)
        )].append(case)
    for (angular_bins, radius), group in sorted(grouped.items()):
        group.sort(key=lambda item: item["energy_bins"])
        for lower, upper in zip(group, group[1:]):
            for model in MODELS:
                for metric, key in (
                    ("mass", "steady_capture_mdot_msun_per_myr"),
                    (
                        "capture_binding_energy",
                        "steady_capture_binding_energy_current_msun_kms2_per_myr",
                    ),
                ):
                    lo = float(lower["models"][model][key])
                    hi = float(upper["models"][model][key])
                    comparisons.append({
                        "model": model,
                        "metric": metric,
                        "fixed_angular_bins": angular_bins,
                        "fixed_reservoir_radius_pc": radius,
                        "lower_energy_bins": lower["energy_bins"],
                        "upper_energy_bins": upper["energy_bins"],
                        "lower_value": lo,
                        "upper_value": hi,
                        "fractional_difference": symmetric_fraction(lo, hi),
                    })

    all_gates = [
        gate
        for case in cases
        for model_gates in case["gates"].values()
        for gate in model_gates.values()
    ]
    passed = bool(all_gates and all(all_gates))
    result = {
        "schema": "finite-angle-ej-capture-tail-diagnostic-v1",
        "status": "CAPTURE_TAIL_DIAGNOSTIC_PASS" if passed else "FAIL",
        "scientific_status": (
            "DIAGNOSTIC_ONLY_DOES_NOT_OVERRIDE_CONVERGENCE_GATES"
        ),
        "cases": cases,
        "adjacent_energy_comparisons": comparisons,
        "note": (
            "Binding-energy fractions are grouped by the pre-collision energy "
            "cell of the captured particle. The diagnostic reuses saved "
            "occupations and does not refit or extrapolate the total current."
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
