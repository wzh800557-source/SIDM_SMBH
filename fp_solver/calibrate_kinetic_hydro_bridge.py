#!/usr/bin/env python3
"""Fit a branch-specific BGK-to-hydrodynamic transition formula.

Only the maintained, spherical radial-inflow family is fitted.  The
hydrodynamic and collisionless endpoints are calculated independently.  A fit
is accepted only when every kinetic run is stationary and conservative, both
endpoints are reached, and the inferred parameters agree across resolutions.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


def bridge_rate(knudsen: np.ndarray, hydrodynamic: float, ballistic: float,
                kn_star: float, exponent: float) -> np.ndarray:
    """Hydrodynamic at small Kn and ballistic at large Kn."""

    knudsen = np.asarray(knudsen, dtype=float)
    return ballistic + (hydrodynamic - ballistic) / (
        1.0 + (knudsen / kn_star) ** exponent
    )


def load_runs(root: Path) -> dict[str, list[dict]]:
    groups: dict[str, list[dict]] = {}
    for path in sorted(root.glob("*.json")):
        try:
            item = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if item.get("schema") != "non-orbit-averaged-radial-bgk-v1":
            continue
        grid = item["grid"]
        label = item.get(
            "resolution_label",
            f"r{grid['radial_bins']}_vr{grid['radial_velocity_bins']}_j{grid['tangential_bins']}",
        )
        item["source_file"] = path.name
        groups.setdefault(label, []).append(item)
    for runs in groups.values():
        runs.sort(key=lambda item: float(item["knudsen_infinity"]))
    return groups


def fit_group(runs: list[dict]) -> dict:
    accepted = [
        run for run in runs
        if run.get("status") == "KINETIC_RUN_PASS"
        and run.get("stationarity", {}).get("status") == "PASS"
    ]
    if len(accepted) < 5:
        return {
            "status": "INSUFFICIENT_ACCEPTED_RUNS",
            "available_runs": len(runs),
            "accepted_runs": len(accepted),
        }
    kn = np.asarray([run["knudsen_infinity"] for run in accepted], dtype=float)
    mdot = np.asarray([
        run["stationarity"]["late_mean_mdot"] for run in accepted
    ], dtype=float)
    hydro_values = np.asarray([
        run["hydrodynamic_reference"]["mdot"] for run in accepted
    ], dtype=float)
    ballistic_values = np.asarray([
        run["ballistic_reference"]["mdot"] for run in accepted
    ], dtype=float)
    if np.ptp(hydro_values) > 1.0e-10 * np.mean(hydro_values):
        raise RuntimeError("hydrodynamic endpoints differ within a resolution group")
    if np.ptp(ballistic_values) > 1.0e-10 * max(np.mean(ballistic_values), 1.0e-30):
        raise RuntimeError("ballistic endpoints differ within a resolution group")
    hydro = float(np.mean(hydro_values))
    ballistic = float(np.mean(ballistic_values))

    fit = least_squares(
        lambda par: (
            bridge_rate(kn, hydro, ballistic, math.exp(par[0]), math.exp(par[1]))
            - mdot
        ) / hydro,
        np.log([0.3, 1.0]),
        bounds=(np.log([1.0e-4, 0.1]), np.log([1.0e4, 10.0])),
        xtol=1.0e-13,
        ftol=1.0e-13,
        gtol=1.0e-13,
        max_nfev=2000,
    )
    kn_star, exponent = np.exp(fit.x)
    predicted = bridge_rate(kn, hydro, ballistic, kn_star, exponent)
    fractional = (predicted - mdot) / hydro
    low_index = int(np.argmin(kn))
    high_index = int(np.argmax(kn))
    gates = {
        "optimizer": bool(fit.success),
        "fit_rms_below_five_percent_hydrodynamic_rate": float(
            np.sqrt(np.mean(fractional**2))
        ) < 0.05,
        "small_kn_reaches_hydrodynamic_endpoint": abs(mdot[low_index] / hydro - 1.0) < 0.10,
        "large_kn_reaches_ballistic_endpoint": abs(
            mdot[high_index] / max(ballistic, 1.0e-30) - 1.0
        ) < 0.10,
    }
    return {
        "status": "RESOLUTION_FIT_PASS" if all(gates.values()) else "RESOLUTION_FIT_DIAGNOSTIC",
        "accepted_runs": len(accepted),
        "knudsen": kn.tolist(),
        "measured_mdot": mdot.tolist(),
        "predicted_mdot": predicted.tolist(),
        "hydrodynamic_mdot": hydro,
        "ballistic_mdot": ballistic,
        "Kn_star": float(kn_star),
        "p": float(exponent),
        "rms_over_hydrodynamic_mdot": float(np.sqrt(np.mean(fractional**2))),
        "smallest_kn_endpoint_ratio": float(mdot[low_index] / hydro),
        "largest_kn_endpoint_ratio": float(mdot[high_index] / max(ballistic, 1.0e-30)),
        "gates": gates,
        "source_files": [run["source_file"] for run in accepted],
    }


def aggregate(root: Path) -> dict:
    groups = load_runs(root)
    fits = {label: fit_group(runs) for label, runs in groups.items()}
    passing = [fit for fit in fits.values() if fit.get("status") == "RESOLUTION_FIT_PASS"]
    resolution_gate = False
    comparison = None
    if len(passing) >= 2:
        ordered = sorted(passing, key=lambda fit: fit["accepted_runs"])
        first, second = ordered[-2:]
        kn_difference = abs(first["Kn_star"] / second["Kn_star"] - 1.0)
        p_difference = abs(first["p"] / second["p"] - 1.0)
        resolution_gate = kn_difference < 0.10 and p_difference < 0.10
        comparison = {
            "Kn_star_relative_difference": kn_difference,
            "p_relative_difference": p_difference,
            "tolerance": 0.10,
            "status": "PASS" if resolution_gate else "FAIL",
        }
    accepted = resolution_gate and len(passing) >= 2
    adopted = None
    if accepted:
        adopted = {
            "Kn_star": float(np.mean([fit["Kn_star"] for fit in passing[-2:]])),
            "p": float(np.mean([fit["p"] for fit in passing[-2:]])),
        }
    return {
        "schema": "radial-bgk-hydrodynamic-transition-calibration-v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "CALIBRATION_ACCEPTED" if accepted else "CALIBRATION_NOT_ACCEPTED",
        "scope": (
            "Maintained, steady, spherical radial inflow with a BGK collision model. "
            "The fit does not apply to a hydrostatic conductive SIDM spike."
        ),
        "formula": (
            "Mdot=Mdot_ball+(Mdot_hydro-Mdot_ball)/"
            "[1+(Kn_infinity/Kn_star)^p]"
        ),
        "groups": fits,
        "resolution_comparison": comparison,
        "adopted_parameters": adopted,
        "gates": {
            "two_resolutions_pass": len(passing) >= 2,
            "parameters_resolution_converged": resolution_gate,
        },
        "interpretation": (
            "Kn_star and p are usable only when status is CALIBRATION_ACCEPTED. "
            "A failed gate leaves the branch formula uncalibrated."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = aggregate(args.run_root)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": result["status"],
        "adopted_parameters": result["adopted_parameters"],
        "gates": result["gates"],
    }, indent=2, sort_keys=True))
    return 0 if result["status"] == "CALIBRATION_ACCEPTED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
