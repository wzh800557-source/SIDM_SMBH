"""Conservative open-control-volume primitive, not a calibrated SIDM closure.

Energy and mass inputs must be independently specified. This module does not
infer capture energy from a boundary dispersion or solve hydrostatic readjustment.
It must not be wired into a production fluid run without those measurements.
"""
import math


def transfer_capture(mass, thermal_energy, black_hole_mass, captured_mass,
                     capture_specific_energy, outward_thermal_energy=0.0,
                     work_on_volume=0.0):
    """Apply an integrated transfer in one consistent mass/energy unit system.

Outward thermal energy is positive leaving the fluid. Work is positive entering
it. A full gravitating calculation must separately supply potential and boundary
work terms, including a moving interface. The black-hole mass update here is
Newtonian captured rest mass; relativistic mass-energy conversion is not modeled.
"""
    vals = (mass, thermal_energy, black_hole_mass, captured_mass,
            capture_specific_energy, outward_thermal_energy, work_on_volume)
    if not all(math.isfinite(float(x)) for x in vals):
        raise ValueError("All budget terms must be finite")
    if mass <= 0 or thermal_energy <= 0 or black_hole_mass < 0:
        raise ValueError("Initial mass and thermal energy must be physical")
    if not 0 <= captured_mass < mass or capture_specific_energy < 0:
        raise ValueError("Invalid capture transfer")
    captured_energy = captured_mass * capture_specific_energy
    remaining_energy = thermal_energy - captured_energy - outward_thermal_energy + work_on_volume
    if remaining_energy <= 0:
        raise ValueError("Transfer would exhaust the thermal control volume")
    remaining_mass = mass - captured_mass
    return {
        "mass": remaining_mass,
        "thermal_energy": remaining_energy,
        "specific_energy": remaining_energy / remaining_mass,
        "black_hole_mass": black_hole_mass + captured_mass,
        "captured_energy": captured_energy,
        "mass_residual": remaining_mass + black_hole_mass + captured_mass - mass - black_hole_mass,
        "thermal_budget_residual": remaining_energy + captured_energy + outward_thermal_energy - work_on_volume - thermal_energy,
    }
