#!/usr/bin/env python3
"""Gate-protected conservative evolution of a moving FP--fluid interface.

The driver advances explicit fluid, kinetic, and black-hole reservoirs while
the ``N_orb=1`` surface moves.  Interface motion transfers the swept mass
between the fluid and kinetic domains.  Capture then transfers mass from the
kinetic domain to the black hole.  Every specific energy and work term is an
input from its owning solver; no capture energy is inferred from the mass
current.

This executable refuses to evolve when the production capture-energy gate is
closed.  A blocked invocation writes a machine-readable record and exits
without changing any reservoir.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import json
import math
from pathlib import Path

from moving_interface_ledger import ReservoirState, advance_interface


def swept_mass_rate(radius0: float, radius1: float, density_midpoint: float,
                    dt: float) -> float:
    """Return the signed mass rate swept into the kinetic domain."""

    values = (radius0, radius1, density_midpoint, dt)
    if not all(math.isfinite(float(value)) for value in values):
        raise ValueError("moving-boundary inputs must be finite")
    if min(radius0, radius1, density_midpoint, dt) <= 0.0:
        raise ValueError("moving-boundary inputs must be positive")
    swept_mass = 4.0 * math.pi * density_midpoint * (radius1**3 - radius0**3) / 3.0
    return swept_mass / dt


def load_trajectory(path: Path) -> list[dict]:
    rows = []
    with path.open() as stream:
        for row in csv.DictReader(stream):
            rows.append({key: float(value) for key, value in row.items()})
    required = {
        "time_myr",
        "r_boundary_pc",
        "rho_boundary_msun_pc3",
        "capture_mdot_msun_per_myr",
        "interface_specific_energy_kms2",
        "captured_specific_energy_kms2",
        "outward_luminosity_msun_kms2_per_myr",
        "work_on_system_msun_kms2_per_myr",
    }
    if len(rows) < 2 or any(required - row.keys() for row in rows):
        raise ValueError(f"trajectory must contain at least two rows and {sorted(required)}")
    if any(rows[index + 1]["time_myr"] <= rows[index]["time_myr"] for index in range(len(rows) - 1)):
        raise ValueError("trajectory times must increase")
    return rows


def evolve(initial: ReservoirState, rows: list[dict]) -> tuple[ReservoirState, list[dict], dict]:
    state = initial
    output = []
    maximum_mass_residual = 0.0
    maximum_energy_residual = 0.0
    for left, right in zip(rows[:-1], rows[1:]):
        dt = right["time_myr"] - left["time_myr"]
        density = 0.5 * (
            left["rho_boundary_msun_pc3"] + right["rho_boundary_msun_pc3"]
        )
        interface_rate = swept_mass_rate(
            left["r_boundary_pc"], right["r_boundary_pc"], density, dt
        )
        mdot = 0.5 * (
            left["capture_mdot_msun_per_myr"]
            + right["capture_mdot_msun_per_myr"]
        )
        interface_energy = 0.5 * (
            left["interface_specific_energy_kms2"]
            + right["interface_specific_energy_kms2"]
        )
        captured_energy = 0.5 * (
            left["captured_specific_energy_kms2"]
            + right["captured_specific_energy_kms2"]
        )
        luminosity = 0.5 * (
            left["outward_luminosity_msun_kms2_per_myr"]
            + right["outward_luminosity_msun_kms2_per_myr"]
        )
        work = 0.5 * (
            left["work_on_system_msun_kms2_per_myr"]
            + right["work_on_system_msun_kms2_per_myr"]
        )
        state, ledger = advance_interface(
            state,
            dt=dt,
            fluid_to_kinetic_mdot=interface_rate,
            kinetic_to_black_hole_mdot=mdot,
            fluid_to_kinetic_specific_energy=interface_energy,
            captured_specific_energy=captured_energy,
            outward_luminosity=luminosity,
            work_on_tracked_system=work,
        )
        maximum_mass_residual = max(maximum_mass_residual, abs(ledger["mass_residual"]))
        maximum_energy_residual = max(
            maximum_energy_residual, abs(ledger["energy_residual_after_work"])
        )
        output.append({
            "time_myr": right["time_myr"],
            "r_boundary_pc": right["r_boundary_pc"],
            "fluid_to_kinetic_mdot_msun_per_myr": interface_rate,
            "capture_mdot_msun_per_myr": mdot,
            "fluid_mass_msun": state.fluid_mass,
            "kinetic_mass_msun": state.kinetic_mass,
            "black_hole_mass_msun": state.black_hole_mass,
            "fluid_energy_msun_kms2": state.fluid_energy,
            "kinetic_energy_msun_kms2": state.kinetic_energy,
            "black_hole_advected_energy_msun_kms2": state.black_hole_advected_energy,
            "exported_energy_msun_kms2": state.exported_energy,
            "mass_residual_msun": ledger["mass_residual"],
            "energy_residual_msun_kms2": ledger["energy_residual_after_work"],
        })
    diagnostic = {
        "maximum_step_mass_residual_msun": maximum_mass_residual,
        "maximum_step_energy_residual_msun_kms2": maximum_energy_residual,
        "mass_conserved": maximum_mass_residual < 1.0e-9 * max(initial.tracked_mass, 1.0),
        "energy_conserved_after_work": maximum_energy_residual < 1.0e-9 * max(
            abs(initial.tracked_energy), 1.0
        ),
    }
    return state, output, diagnostic


def write_blocked(path: Path, gate: dict) -> None:
    result = {
        "schema": "conservative-moving-boundary-v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "BLOCKED_BY_ENERGY_CURRENT_GATE",
        "fluid_coupling_authorized": False,
        "failed_gates": [key for key, value in gate.get("gates", {}).items() if value is not True],
        "reservoir_advanced": False,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--energy-gate", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--fluid-mass", type=float, required=True)
    parser.add_argument("--kinetic-mass", type=float, required=True)
    parser.add_argument("--black-hole-mass", type=float, required=True)
    parser.add_argument("--fluid-energy", type=float, required=True)
    parser.add_argument("--kinetic-energy", type=float, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_path = args.outdir / "moving_boundary_summary.json"
    gate = json.loads(args.energy_gate.read_text())
    if (
        gate.get("status") != "ENERGY_CURRENT_GATE_PASS"
        or gate.get("fluid_coupling_authorized") is not True
        or not gate.get("gates")
        or any(value is not True for value in gate["gates"].values())
    ):
        write_blocked(summary_path, gate)
        print("BLOCKED_BY_ENERGY_CURRENT_GATE")
        return 4

    initial = ReservoirState(
        fluid_mass=args.fluid_mass,
        kinetic_mass=args.kinetic_mass,
        black_hole_mass=args.black_hole_mass,
        fluid_energy=args.fluid_energy,
        kinetic_energy=args.kinetic_energy,
    )
    rows = load_trajectory(args.trajectory)
    final, output, diagnostic = evolve(initial, rows)
    trajectory_path = args.outdir / "moving_boundary_trajectory.csv"
    with trajectory_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    gates = {
        "mass_conserved": diagnostic["mass_conserved"],
        "energy_conserved_after_work": diagnostic["energy_conserved_after_work"],
    }
    result = {
        "schema": "conservative-moving-boundary-v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "status": "MOVING_BOUNDARY_PASS" if all(gates.values()) else "MOVING_BOUNDARY_FAIL",
        "fluid_coupling_authorized": True,
        "reservoir_advanced": True,
        "interface_definition": "N_orb=1 radius supplied by the live fluid/kinetic owner",
        "mass_transfer": "4*pi*rho*(r_new^3-r_old^3)/3 per step",
        "initial_state": asdict(initial),
        "final_state": asdict(final),
        "diagnostic": diagnostic,
        "gates": gates,
        "trajectory_csv": trajectory_path.name,
    }
    summary_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(result["status"])
    return 0 if result["status"] == "MOVING_BOUNDARY_PASS" else 5


if __name__ == "__main__":
    raise SystemExit(main())
