"""Conservation tests for the moving-interface mass and energy ledger."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from moving_interface_ledger import ReservoirState, advance_interface  # noqa: E402


def main():
    initial = ReservoirState(
        fluid_mass=100.0,
        kinetic_mass=10.0,
        black_hole_mass=1.0e6,
        fluid_energy=-300.0,
        kinetic_energy=-80.0,
    )
    updated, ledger = advance_interface(
        initial,
        dt=0.2,
        fluid_to_kinetic_mdot=3.0,
        kinetic_to_black_hole_mdot=1.0,
        fluid_to_kinetic_specific_energy=-2.0,
        captured_specific_energy=-7.0,
        outward_luminosity=0.5,
        work_on_tracked_system=0.25,
    )
    assert abs(ledger["mass_residual"]) < 1.0e-10
    assert abs(ledger["energy_residual_after_work"]) < 1.0e-12
    assert abs(updated.fluid_mass - 99.4) < 1.0e-12
    assert abs(updated.kinetic_mass - 10.4) < 1.0e-12
    assert abs(updated.black_hole_mass - 1000000.2) < 1.0e-10

    # An outward-moving interface is the same conservative transfer with a
    # negative fluid-to-kinetic rate.
    returned, ledger2 = advance_interface(
        updated,
        dt=0.1,
        fluid_to_kinetic_mdot=-1.0,
        kinetic_to_black_hole_mdot=0.0,
        fluid_to_kinetic_specific_energy=-3.0,
        captured_specific_energy=-7.0,
    )
    assert abs(ledger2["mass_residual"]) < 1.0e-10
    assert abs(ledger2["energy_residual_after_work"]) < 1.0e-12
    assert returned.fluid_mass > updated.fluid_mass
    print("PASS: moving-interface mass and energy conservation")


if __name__ == "__main__":
    main()
