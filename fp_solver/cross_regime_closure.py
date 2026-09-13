#!/usr/bin/env python3
"""Mass-current limits and regime checks for an SIDM--SMBH interface.

The orbit-resolved, conductive-spike, and Bondi solutions are different
boundary-value problems.  This module therefore keeps them as separate
branches instead of forcing a Fokker--Planck coefficient to approach Bondi
accretion as the cross-section is increased.

The conductive-spike calibration follows the constant, isotropic
cross-section simulations of Sabarish et al. (2025),
``dotM = C + 1 / (A s + B / s)``.  It is an external benchmark, not a
calibration of the velocity-dependent t-channel model used by the production
halo in this repository.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Iterable

import numpy as np


G_PC_KMS2_MSUN = 4.30091e-3
KMS_TO_PC_PER_MYR = 1.0227121650537077


def _positive_finite(name: str, value: float) -> float:
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return value


def current_scale_msun_per_myr(
    radius_pc: float,
    density_msun_pc3: float,
    sigma_1d_kms: float,
) -> float:
    """Return ``4*pi*r^2*rho*sigma`` in solar masses per Myr."""

    radius_pc = _positive_finite("radius_pc", radius_pc)
    density_msun_pc3 = _positive_finite("density_msun_pc3", density_msun_pc3)
    sigma_1d_kms = _positive_finite("sigma_1d_kms", sigma_1d_kms)
    return (
        4.0
        * math.pi
        * radius_pc**2
        * density_msun_pc3
        * sigma_1d_kms
        * KMS_TO_PC_PER_MYR
    )


def diagnostic_mass_coefficient(
    mdot_msun_per_myr: float,
    radius_pc: float,
    density_msun_pc3: float,
    sigma_1d_kms: float,
) -> float:
    """Return the diagnostic coefficient C_M for a measured mass current."""

    mdot = _positive_finite("mdot_msun_per_myr", mdot_msun_per_myr)
    return mdot / current_scale_msun_per_myr(
        radius_pc, density_msun_pc3, sigma_1d_kms
    )


def bondi_lambda(gamma: float = 5.0 / 3.0) -> float:
    """Dimensionless critical Bondi rate for a polytropic ideal gas.

    The isothermal limit is handled analytically.  For ``gamma=5/3`` the
    transonic solution has ``lambda_B=1/4``.
    """

    gamma = float(gamma)
    if not math.isfinite(gamma) or gamma < 1.0 or gamma > 5.0 / 3.0:
        raise ValueError("gamma must lie in [1, 5/3]")
    if math.isclose(gamma, 1.0, rel_tol=0.0, abs_tol=1.0e-12):
        return math.exp(1.5) / 4.0
    if math.isclose(gamma, 5.0 / 3.0, rel_tol=0.0, abs_tol=1.0e-12):
        return 0.25
    first = 0.5 ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    exponent = -(5.0 - 3.0 * gamma) / (2.0 * (gamma - 1.0))
    second = ((5.0 - 3.0 * gamma) / 4.0) ** exponent
    return first * second


def bondi_radius_pc(
    black_hole_mass_msun: float,
    sound_speed_kms: float,
) -> float:
    """Return ``G M_bh / c_infinity^2`` in pc."""

    mass = _positive_finite("black_hole_mass_msun", black_hole_mass_msun)
    sound_speed = _positive_finite("sound_speed_kms", sound_speed_kms)
    return G_PC_KMS2_MSUN * mass / sound_speed**2


def bondi_rate_msun_per_myr(
    black_hole_mass_msun: float,
    density_infinity_msun_pc3: float,
    sound_speed_infinity_kms: float,
    gamma: float = 5.0 / 3.0,
) -> float:
    """Return the steady spherical Bondi rate in solar masses per Myr."""

    mass = _positive_finite("black_hole_mass_msun", black_hole_mass_msun)
    density = _positive_finite(
        "density_infinity_msun_pc3", density_infinity_msun_pc3
    )
    sound_speed = _positive_finite(
        "sound_speed_infinity_kms", sound_speed_infinity_kms
    )
    return (
        4.0
        * math.pi
        * bondi_lambda(gamma)
        * G_PC_KMS2_MSUN**2
        * mass**2
        * density
        / sound_speed**3
        * KMS_TO_PC_PER_MYR
    )


def sound_speed_from_sigma_1d(sigma_1d_kms: float, gamma: float = 5.0 / 3.0) -> float:
    """Convert one-dimensional random dispersion to ideal-gas sound speed."""

    sigma = _positive_finite("sigma_1d_kms", sigma_1d_kms)
    gamma = _positive_finite("gamma", gamma)
    return math.sqrt(gamma) * sigma


@dataclass(frozen=True)
class ConductiveSpikeCalibration:
    """Constant-cross-section conductive-spike benchmark.

    ``s`` is in cm^2/g and rates are in solar masses per year.  ``baseline``
    is supplied separately because it depends on the initial loss-cone
    population and outer reservoir.  The default A and B are the published
    Sabarish et al. (2025) fit values.
    """

    A: float = 0.056
    B: float = 0.27
    baseline_msun_per_yr: float = 1.7
    source: str = "Sabarish et al. 2025, arXiv:2505.14779"

    def __post_init__(self) -> None:
        _positive_finite("A", self.A)
        _positive_finite("B", self.B)
        if not math.isfinite(float(self.baseline_msun_per_yr)) or self.baseline_msun_per_yr < 0:
            raise ValueError("baseline_msun_per_yr must be non-negative and finite")

    @property
    def sigma_transition_cm2_g(self) -> float:
        return math.sqrt(self.B / self.A)

    @property
    def maximum_excess_msun_per_yr(self) -> float:
        return 1.0 / (2.0 * math.sqrt(self.A * self.B))

    def excess_rate_msun_per_yr(self, sigma_over_m_cm2_g: float | np.ndarray):
        s = np.asarray(sigma_over_m_cm2_g, dtype=float)
        if np.any(~np.isfinite(s)) or np.any(s <= 0.0):
            raise ValueError("sigma_over_m_cm2_g must be positive and finite")
        out = 1.0 / (self.A * s + self.B / s)
        return float(out) if out.ndim == 0 else out

    def total_rate_msun_per_yr(self, sigma_over_m_cm2_g: float | np.ndarray):
        return self.baseline_msun_per_yr + self.excess_rate_msun_per_yr(
            sigma_over_m_cm2_g
        )

    def normalized_excess(self, sigma_over_m_cm2_g: float | np.ndarray):
        """Return ``2*x/(1+x^2)`` with ``x=s/s_transition``."""

        s = np.asarray(sigma_over_m_cm2_g, dtype=float)
        x = s / self.sigma_transition_cm2_g
        out = 2.0 * x / (1.0 + x * x)
        return float(out) if out.ndim == 0 else out

    def as_dict(self) -> dict:
        out = asdict(self)
        out.update(
            {
                "sigma_transition_cm2_g": self.sigma_transition_cm2_g,
                "maximum_excess_msun_per_yr": self.maximum_excess_msun_per_yr,
                "lmfp_scaling": "excess rate proportional to sigma/m",
                "smfp_scaling": "excess rate proportional to (sigma/m)^-1",
            }
        )
        return out


def interpolate_log_profile(radius: Iterable[float], value: Iterable[float], target: float) -> float:
    """Log-log interpolation for a positive radial profile."""

    r = np.asarray(radius, dtype=float)
    y = np.asarray(value, dtype=float)
    target = _positive_finite("target", target)
    if (
        r.ndim != 1
        or y.shape != r.shape
        or r.size < 2
        or np.any(~np.isfinite(r))
        or np.any(~np.isfinite(y))
        or np.any(r <= 0.0)
        or np.any(y <= 0.0)
        or np.any(np.diff(r) <= 0.0)
        or target < r[0]
        or target > r[-1]
    ):
        raise ValueError("invalid positive profile or target outside its domain")
    return float(np.exp(np.interp(np.log(target), np.log(r), np.log(y))))


def radial_collisionality_gate(
    radius_pc: Iterable[float],
    n_orb: Iterable[float],
    capture_radius_pc: float,
    supply_radius_pc: float,
    hydrodynamic_n_orb_threshold: float = 20.0 * math.pi,
) -> dict:
    """Audit whether a hydrodynamic inflow is collisional across its full path.

    A local value at the Bondi radius is insufficient.  Bondi is a candidate
    only when orbital memory is erased from the supply radius down to the
    capture surface.  The physical orbit/kinetic crossing remains
    ``N_orb=1``.  The default continuum sentinel is ``N_orb=20*pi``; when the
    pressure scale height is comparable to radius, this corresponds to
    ``Kn=lambda/H=0.1``.  It is an operational gate, not a fitted transition.
    """

    r = np.asarray(radius_pc, dtype=float)
    n = np.asarray(n_orb, dtype=float)
    capture = _positive_finite("capture_radius_pc", capture_radius_pc)
    supply = _positive_finite("supply_radius_pc", supply_radius_pc)
    threshold = _positive_finite(
        "hydrodynamic_n_orb_threshold", hydrodynamic_n_orb_threshold
    )
    if supply <= capture:
        raise ValueError("supply_radius_pc must exceed capture_radius_pc")
    if (
        r.ndim != 1
        or n.shape != r.shape
        or r.size < 2
        or np.any(~np.isfinite(r))
        or np.any(~np.isfinite(n))
        or np.any(r <= 0.0)
        or np.any(n <= 0.0)
        or np.any(np.diff(r) <= 0.0)
        or capture < r[0] * (1.0 - 1.0e-10)
        or supply > r[-1] * (1.0 + 1.0e-10)
    ):
        raise ValueError("invalid collisionality profile or requested interval")

    sample_r = np.unique(
        np.r_[capture, r[(r > capture) & (r < supply)], supply]
    )
    sample_n = np.exp(np.interp(np.log(sample_r), np.log(r), np.log(n)))
    minimum_index = int(np.argmin(sample_n))
    minimum = float(sample_n[minimum_index])
    if minimum < 1.0:
        status = "ORBIT_MEMORY_PRESENT"
        selected_solver = "FP_INNER_WITH_KINETIC_OR_FLUID_OUTER"
    elif minimum < threshold:
        status = "TRANSITIONAL_COLLISIONALITY"
        selected_solver = "NON_ORBIT_AVERAGED_KINETIC"
    else:
        status = "HYDRODYNAMIC_COLLISIONALITY_PASS"
        selected_solver = "HYDRODYNAMIC_CANDIDATE"
    return {
        "status": status,
        "selected_solver": selected_solver,
        "capture_radius_pc": capture,
        "supply_radius_pc": supply,
        "N_orb_at_capture": float(sample_n[0]),
        "N_orb_at_supply": float(sample_n[-1]),
        "minimum_N_orb": minimum,
        "radius_of_minimum_N_orb_pc": float(sample_r[minimum_index]),
        "hydrodynamic_n_orb_threshold": threshold,
        "maximum_knudsen_if_H_equals_r": 2.0 * math.pi / minimum,
        "continuum_gate_interpretation": (
            "For H approximately r, the default N_orb threshold corresponds "
            "to Kn<=0.1. A hydrodynamic solve must evaluate its own scale height."
        ),
        "physical_crossing": "N_orb=1",
        "hydrodynamic_collisionality_pass": status
        == "HYDRODYNAMIC_COLLISIONALITY_PASS",
        "additional_bondi_requirements": [
            "steady spherical inflow",
            "negligible rotation",
            "outer reservoir fixes density and sound speed",
            "transonic solution connects to the capture surface",
        ],
    }


def select_mass_current_branch(
    gate: dict,
    *,
    radial_flow_family: str,
    flow_conditions_verified: bool = False,
) -> str:
    """Select a solver without interpolating unrelated boundary problems.

    Crossing the collisionality gate only establishes a hydrodynamic candidate.
    The caller must separately verify the hydrostatic-conductive or maintained
    radial-inflow boundary conditions before this routine returns its endpoint.
    """

    allowed = {"hydrostatic-conductive", "maintained-radial-inflow"}
    if radial_flow_family not in allowed:
        raise ValueError(f"radial_flow_family must be one of {sorted(allowed)}")
    status = gate.get("status")
    if status == "ORBIT_MEMORY_PRESENT":
        return "fp"
    if status == "TRANSITIONAL_COLLISIONALITY":
        return "radial-kinetic"
    if status != "HYDRODYNAMIC_COLLISIONALITY_PASS":
        raise ValueError("unrecognized collisionality gate")
    if not flow_conditions_verified:
        return "hydrodynamic-candidate"
    if radial_flow_family == "hydrostatic-conductive":
        return "conductive-spike"
    return "bondi"
