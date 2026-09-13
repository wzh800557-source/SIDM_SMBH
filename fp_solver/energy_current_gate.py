#!/usr/bin/env python3
"""Evaluate the production capture-energy gate before fluid coupling."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path

import numpy as np


def evaluate(summary: dict, moment_tolerance: float = 0.01,
             convergence_tolerance: float = 0.05,
             seed_tolerance: float = 0.05) -> dict:
    rows = summary.get("rows", [])
    if not isinstance(rows, list):
        raise ValueError("energy summary rows must be a list")
    required = []
    for row in rows:
        try:
            audit = row["energy_audit"]
            required.append({
                "energy_bins": int(row["energy_bins"]),
                "seed": int(row["seed"]),
                "model": str(row["model"]),
                "mdot": float(row["mdot_msun_per_myr"]),
                "energy_current": float(audit["capture_orbital_energy_rate"]),
                "internal_l1": float(audit["internal_l1_over_boundary_scale"]),
                "moment_status": str(audit["status"]),
                "pair_energy": bool(audit["gates"]["physical_pair_energy_conserved"]),
                "ledger_identity": bool(audit["gates"]["ledger_identity"]),
            })
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"malformed energy row: {error}") from error

    bins = sorted({row["energy_bins"] for row in required})
    models = sorted({row["model"] for row in required})
    highest_two = bins[-2:] if len(bins) >= 2 else bins
    convergence = {}
    seed_checks = {}
    for model in models:
        means = {}
        for energy_bins in bins:
            subset = [
                row for row in required
                if row["model"] == model and row["energy_bins"] == energy_bins
            ]
            if subset:
                energy = np.asarray([row["energy_current"] for row in subset])
                mass = np.asarray([row["mdot"] for row in subset])
                means[energy_bins] = {
                    "energy_current_mean": float(np.mean(energy)),
                    "mdot_mean": float(np.mean(mass)),
                }
                seed_checks[f"{model}_E{energy_bins}"] = {
                    "energy_current_fractional_range": float(
                        np.ptp(energy) / max(abs(float(np.mean(energy))), 1.0e-300)
                    ),
                    "mdot_fractional_range": float(
                        np.ptp(mass) / max(abs(float(np.mean(mass))), 1.0e-300)
                    ),
                }
        if len(highest_two) == 2 and all(value in means for value in highest_two):
            low, high = highest_two
            convergence[model] = {
                "energy_bins": [low, high],
                "capture_energy_relative_change": abs(
                    means[high]["energy_current_mean"]
                    / means[low]["energy_current_mean"] - 1.0
                ),
                "capture_mass_relative_change": abs(
                    means[high]["mdot_mean"] / means[low]["mdot_mean"] - 1.0
                ),
            }

    gates = {
        "summary_complete": (
            summary.get("status") == "FOLLOWUP_DIAGNOSTICS_COMPLETE"
            and not summary.get("missing_runs")
            and len(required) > 0
        ),
        "physical_pair_energy_conserved": bool(required) and all(
            row["pair_energy"] for row in required
        ),
        "ledger_identities": bool(required) and all(
            row["ledger_identity"] for row in required
        ),
        "internal_energy_moments_stationary": bool(required) and all(
            row["moment_status"] == "MOMENT_AUDIT_PASS"
            and row["internal_l1"] <= moment_tolerance
            for row in required
        ),
        "energy_grid_converged": bool(convergence) and all(
            value["capture_energy_relative_change"] <= convergence_tolerance
            for value in convergence.values()
        ),
        "mass_grid_converged": bool(convergence) and all(
            value["capture_mass_relative_change"] <= convergence_tolerance
            for value in convergence.values()
        ),
        "seed_converged": bool(seed_checks) and all(
            value["energy_current_fractional_range"] <= seed_tolerance
            and value["mdot_fractional_range"] <= seed_tolerance
            for value in seed_checks.values()
        ),
        "upstream_authorized": summary.get("fluid_coupling_authorized") is True,
    }
    authorized = all(gates.values())
    maximum_internal = max((row["internal_l1"] for row in required), default=math.inf)
    return {
        "schema": "capture-energy-current-gate-v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "ENERGY_CURRENT_GATE_PASS" if authorized else "ENERGY_CURRENT_GATE_FAIL",
        "fluid_coupling_authorized": authorized,
        "thresholds": {
            "internal_moment_relative": moment_tolerance,
            "grid_convergence_relative": convergence_tolerance,
            "seed_convergence_relative": seed_tolerance,
        },
        "gates": gates,
        "maximum_internal_l1_over_boundary_scale": maximum_internal,
        "resolution_convergence": convergence,
        "seed_convergence": seed_checks,
        "rows_evaluated": len(required),
        "decision": (
            "The conservative moving-boundary evolution may run."
            if authorized else
            "The capture-energy current is not authorized for fluid coupling."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--moment-tolerance", type=float, default=0.01)
    parser.add_argument("--convergence-tolerance", type=float, default=0.05)
    parser.add_argument("--seed-tolerance", type=float, default=0.05)
    args = parser.parse_args()
    result = evaluate(
        json.loads(args.summary.read_text()),
        args.moment_tolerance,
        args.convergence_tolerance,
        args.seed_tolerance,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["fluid_coupling_authorized"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
