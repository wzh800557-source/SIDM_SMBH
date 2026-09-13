#!/usr/bin/env python3
"""Steady spherical hydrodynamic reference for the kinetic transition runs.

The dimensionless normalization is

``G M_bh = rho_infinity = c_infinity = r_B = 1``.

The solver evaluates the critical polytropic Bondi solution from mass and
Bernoulli conservation.  It is deliberately separate from the kinetic code:
the hydrodynamic calculation supplies an endpoint, while the intermediate
Knudsen-number dependence must be measured with a kinetic equation.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq


def bondi_lambda(gamma: float) -> float:
    """Return the critical dimensionless Bondi eigenvalue."""

    gamma = float(gamma)
    if not 1.0 < gamma <= 5.0 / 3.0:
        raise ValueError("gamma must lie in (1, 5/3]")
    if math.isclose(gamma, 5.0 / 3.0, rel_tol=0.0, abs_tol=1.0e-13):
        return 0.25
    exponent = (5.0 - 3.0 * gamma) / (2.0 * (gamma - 1.0))
    return 0.25 * (2.0 / (5.0 - 3.0 * gamma)) ** exponent


def _bernoulli_residual(sound_speed_squared: float, radius: float, gamma: float,
                        eigenvalue: float) -> float:
    density = sound_speed_squared ** (1.0 / (gamma - 1.0))
    speed = eigenvalue / (radius * radius * density)
    return (
        0.5 * speed * speed
        + sound_speed_squared / (gamma - 1.0)
        - 1.0 / radius
        - 1.0 / (gamma - 1.0)
    )


def _positive_roots(radius: float, gamma: float, eigenvalue: float) -> list[float]:
    """Find the two positive roots around the Bernoulli minimum."""

    exponent = (gamma - 1.0) / (gamma + 1.0)
    minimum = (eigenvalue * eigenvalue / radius**4) ** exponent
    f_minimum = _bernoulli_residual(minimum, radius, gamma, eigenvalue)
    scale = 1.0 / (gamma - 1.0) + 1.0 / radius
    if abs(f_minimum) < 2.0e-11 * scale:
        return [float(minimum)]
    if f_minimum > 0.0:
        return []
    low = minimum * 1.0e-12
    high = minimum * 1.0e12
    left = float(brentq(
        _bernoulli_residual, low, minimum,
        args=(radius, gamma, eigenvalue), xtol=1.0e-13, rtol=1.0e-13,
    ))
    right = float(brentq(
        _bernoulli_residual, minimum, high,
        args=(radius, gamma, eigenvalue), xtol=1.0e-13, rtol=1.0e-13,
    ))
    return [left, right]


def critical_bondi_profile(radius: np.ndarray, gamma: float = 5.0 / 3.0) -> dict:
    """Evaluate the critical inflow branch at positive radii.

    ``sound_speed`` is the adiabatic sound speed.  The one-dimensional
    Maxwellian variance used by the BGK solver is ``temperature=c_s^2/gamma``.
    """

    radius = np.asarray(radius, dtype=float)
    if radius.ndim != 1 or radius.size < 2 or np.any(~np.isfinite(radius)):
        raise ValueError("radius must be a finite one-dimensional array")
    if np.any(radius <= 0.0) or np.any(np.diff(radius) <= 0.0):
        raise ValueError("radius must be positive and strictly increasing")
    eigenvalue = bondi_lambda(gamma)
    sound2 = np.empty_like(radius)
    mach = np.empty_like(radius)
    for index, r_value in enumerate(radius):
        roots = _positive_roots(float(r_value), gamma, eigenvalue)
        if not roots:
            raise RuntimeError(f"Bondi root not bracketed at r={r_value:g}")
        candidates = []
        for root in roots:
            density = root ** (1.0 / (gamma - 1.0))
            speed = eigenvalue / (r_value * r_value * density)
            candidates.append((abs(speed) / math.sqrt(root), root))
        sonic_radius = max(0.0, (5.0 - 3.0 * gamma) / 4.0)
        if r_value < sonic_radius * (1.0 - 2.0e-6):
            supersonic = [item for item in candidates if item[0] >= 1.0 - 2.0e-6]
            if not supersonic:
                raise RuntimeError(f"supersonic Bondi root not found at r={r_value:g}")
            chosen_mach, chosen = min(supersonic, key=lambda item: item[1])
        else:
            subsonic = [item for item in candidates if item[0] <= 1.0 + 2.0e-6]
            if subsonic:
                chosen_mach, chosen = max(subsonic, key=lambda item: item[1])
            else:
                chosen_mach, chosen = min(
                    candidates, key=lambda item: abs(item[0] - 1.0)
                )
        sound2[index] = chosen
        mach[index] = chosen_mach

    density = sound2 ** (1.0 / (gamma - 1.0))
    radial_velocity = -eigenvalue / (radius * radius * density)
    temperature = sound2 / gamma
    mdot = -4.0 * math.pi * radius**2 * density * radial_velocity
    bernoulli = (
        0.5 * radial_velocity**2
        + sound2 / (gamma - 1.0)
        - 1.0 / radius
    )
    expected_bernoulli = 1.0 / (gamma - 1.0)
    return {
        "radius": radius,
        "density": density,
        "radial_velocity": radial_velocity,
        "sound_speed": np.sqrt(sound2),
        "temperature": temperature,
        "mach": mach,
        "mdot": mdot,
        "eigenvalue": eigenvalue,
        "mdot_exact": 4.0 * math.pi * eigenvalue,
        "maximum_mass_flux_relative_residual": float(
            np.max(np.abs(mdot / (4.0 * math.pi * eigenvalue) - 1.0))
        ),
        "maximum_bernoulli_absolute_residual": float(
            np.max(np.abs(bernoulli - expected_bernoulli))
        ),
        "normalization": {
            "G_M_bh": 1.0,
            "rho_infinity": 1.0,
            "c_infinity": 1.0,
            "r_B": 1.0,
        },
    }


def serializable_profile(profile: dict) -> dict:
    out = {}
    for key, value in profile.items():
        out[key] = value.tolist() if isinstance(value, np.ndarray) else value
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gamma", type=float, default=5.0 / 3.0)
    parser.add_argument("--r-min", type=float, default=0.1)
    parser.add_argument("--r-max", type=float, default=8.0)
    parser.add_argument("--radial-bins", type=int, default=128)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    radius = np.geomspace(args.r_min, args.r_max, args.radial_bins)
    result = serializable_profile(critical_bondi_profile(radius, args.gamma))
    result.update({
        "schema": "critical-bondi-hydrodynamic-reference-v1",
        "status": "PASS",
        "gamma": args.gamma,
    })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": result["status"],
        "gamma": args.gamma,
        "mdot_exact": result["mdot_exact"],
        "mass_residual": result["maximum_mass_flux_relative_residual"],
        "bernoulli_residual": result["maximum_bernoulli_absolute_residual"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
