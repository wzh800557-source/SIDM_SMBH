#!/usr/bin/env python3
"""Conservative mass and energy ledger for a moving FP--fluid interface.

This is a bookkeeping primitive.  The caller must supply independently
measured specific energies, conductive luminosity, and pressure/gravitational
work.  The module does not infer an energy current from the capture rate.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math


@dataclass(frozen=True)
class ReservoirState:
    fluid_mass: float
    kinetic_mass: float
    black_hole_mass: float
    fluid_energy: float
    kinetic_energy: float
    black_hole_advected_energy: float = 0.0
    exported_energy: float = 0.0

    def validate(self) -> None:
        values = asdict(self)
        if not all(math.isfinite(float(value)) for value in values.values()):
            raise ValueError("all reservoir entries must be finite")
        if min(self.fluid_mass, self.kinetic_mass, self.black_hole_mass) < 0.0:
            raise ValueError("reservoir masses must be non-negative")

    @property
    def tracked_mass(self) -> float:
        return self.fluid_mass + self.kinetic_mass + self.black_hole_mass

    @property
    def tracked_energy(self) -> float:
        return (
            self.fluid_energy
            + self.kinetic_energy
            + self.black_hole_advected_energy
            + self.exported_energy
        )


def advance_interface(
    state: ReservoirState,
    *,
    dt: float,
    fluid_to_kinetic_mdot: float,
    kinetic_to_black_hole_mdot: float,
    fluid_to_kinetic_specific_energy: float,
    captured_specific_energy: float,
    outward_luminosity: float = 0.0,
    work_on_tracked_system: float = 0.0,
) -> tuple[ReservoirState, dict]:
    """Advance a three-reservoir control volume by one explicit step.

    Positive ``fluid_to_kinetic_mdot`` moves material into the FP domain;
    it may be negative when an outward-moving interface returns material to
    the fluid.  Capture is non-negative.  Positive luminosity leaves the
    tracked fluid-plus-kinetic system, and positive work enters it.
    """

    state.validate()
    values = (
        dt,
        fluid_to_kinetic_mdot,
        kinetic_to_black_hole_mdot,
        fluid_to_kinetic_specific_energy,
        captured_specific_energy,
        outward_luminosity,
        work_on_tracked_system,
    )
    if not all(math.isfinite(float(value)) for value in values):
        raise ValueError("all transfer entries must be finite")
    if dt <= 0.0 or kinetic_to_black_hole_mdot < 0.0:
        raise ValueError("dt must be positive and capture rate non-negative")

    dm_fk = dt * fluid_to_kinetic_mdot
    dm_cap = dt * kinetic_to_black_hole_mdot
    de_fk = dm_fk * fluid_to_kinetic_specific_energy
    de_cap = dm_cap * captured_specific_energy
    de_export = dt * outward_luminosity
    de_work = dt * work_on_tracked_system

    updated = ReservoirState(
        fluid_mass=state.fluid_mass - dm_fk,
        kinetic_mass=state.kinetic_mass + dm_fk - dm_cap,
        black_hole_mass=state.black_hole_mass + dm_cap,
        fluid_energy=state.fluid_energy - de_fk - de_export + de_work,
        kinetic_energy=state.kinetic_energy + de_fk - de_cap,
        black_hole_advected_energy=state.black_hole_advected_energy + de_cap,
        exported_energy=state.exported_energy + de_export,
    )
    updated.validate()
    if min(updated.fluid_mass, updated.kinetic_mass) < 0.0:
        raise ValueError("step transfers more mass than a reservoir contains")

    mass_residual = updated.tracked_mass - state.tracked_mass
    energy_residual = updated.tracked_energy - state.tracked_energy - de_work
    return updated, {
        "mass_residual": mass_residual,
        "energy_residual_after_work": energy_residual,
        "fluid_to_kinetic_mass": dm_fk,
        "captured_mass": dm_cap,
        "fluid_to_kinetic_energy": de_fk,
        "captured_advected_energy": de_cap,
        "exported_energy": de_export,
        "work_on_tracked_system": de_work,
    }

