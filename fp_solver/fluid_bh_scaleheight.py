#!/usr/bin/env python3
"""Point-mass correction to the gravothermal LMFP scale height."""

from __future__ import annotations

import numpy as np


def bh_lmfp_scaleheight_factor(
    r: np.ndarray, rho: np.ndarray, mbh: float, r_soft: float = 0.0
) -> np.ndarray:
    """Return H_eff^2/H_J^2 for a central softened point mass.

    GravothermalSIDM uses H_J^2=v^2/rho in its dimensionless units.  Adding
    the point-mass orbital frequency gives

        H_eff^2/H_J^2 = rho/[rho + M_bh/(r^2+r_soft^2)^(3/2)].
    """

    r = np.asarray(r, float)
    rho = np.asarray(rho, float)
    if np.any(r <= 0.0) or np.any(rho <= 0.0) or mbh < 0.0 or r_soft < 0.0:
        raise ValueError("scale-height inputs must be physical")
    omega_bh_sq = float(mbh) / (r**2 + float(r_soft) ** 2) ** 1.5
    return rho / (rho + omega_bh_sq)


def bh_lmfp_scaleheight_factor_physical(r, rho, mbh, r_soft=0.0):
    """Physical-unit form with M_s=4*pi*rho_s*r_s**3 and H_J^2=sigma^2/(4*pi*G*rho).

    Radius and softening share one length unit, mass and density share the
    corresponding mass and volume units. No change to the dimensionless solver.
    """
    return bh_lmfp_scaleheight_factor(r, 4.0*np.pi*np.asarray(rho, float), mbh, r_soft)


def bh_corrected_lmfp_inverse(
    kinv_lmfp: np.ndarray,
    r: np.ndarray,
    rho: np.ndarray,
    mbh: float,
    r_soft: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the BH-corrected LMFP inverse conductivity and H^2 ratio."""

    kinv = np.asarray(kinv_lmfp, float)
    factor = bh_lmfp_scaleheight_factor(r, rho, mbh, r_soft)
    if kinv.shape != factor.shape or np.any(kinv <= 0.0):
        raise ValueError("LMFP inverse conductivity must match the fluid grid")
    return kinv / factor, factor


def harmonic_effective_conductivity(
    kinv_smfp: np.ndarray, kinv_lmfp: np.ndarray
) -> np.ndarray:
    """Harmonic interpolation of the SMFP and corrected LMFP branches."""

    smfp, lmfp = np.broadcast_arrays(
        np.asarray(kinv_smfp, float), np.asarray(kinv_lmfp, float)
    )
    if np.any(smfp <= 0.0) or np.any(lmfp <= 0.0):
        raise ValueError("inverse conductivities must be positive")
    return 1.0 / (smfp + lmfp)


def shell_luminosity(
    r: np.ndarray, specific_energy: np.ndarray, conductivity: np.ndarray
) -> np.ndarray:
    """Return the stock finite-difference luminosity with zero outer flux."""

    r = np.asarray(r, float)
    u = np.asarray(specific_energy, float)
    k = np.asarray(conductivity, float)
    if not (r.ndim == u.ndim == k.ndim == 1 and r.size == u.size == k.size):
        raise ValueError("transport arrays must be one-dimensional and aligned")
    if r.size < 2 or np.any(np.diff(r) <= 0.0) or np.any(k <= 0.0):
        raise ValueError("invalid fluid transport grid")
    luminosity = np.zeros_like(r)
    luminosity[1:-1] = (
        -r[1:-1] ** 2 * (k[1:-1] + k[2:]) / 2.0
        * (u[2:] - u[1:-1])
        / ((r[2:] - r[:-2]) / 2.0)
    )
    luminosity[0] = (
        -r[0] ** 2 * (k[0] + k[1]) / 2.0
        * (u[1] - u[0]) / (r[1] / 2.0)
    )
    luminosity[-1] = 0.0
    return luminosity


def conservative_specific_energy_rate(
    enclosed_mass: np.ndarray,
    luminosity: np.ndarray,
    inner_luminosity: float,
) -> np.ndarray:
    """Finite-volume du/dt whose mass-weighted sum is L_inner-L_outer."""

    mass = np.asarray(enclosed_mass, float)
    flux = np.asarray(luminosity, float)
    if mass.shape != flux.shape or mass.ndim != 1 or mass.size < 2:
        raise ValueError("mass and luminosity arrays must be aligned")
    shell_mass = np.diff(np.r_[0.0, mass])
    if np.any(shell_mass <= 0.0):
        raise ValueError("enclosed mass must increase strictly")
    rate = np.empty_like(mass)
    rate[0] = (float(inner_luminosity) - flux[0]) / shell_mass[0]
    rate[1:] = -(flux[1:] - flux[:-1]) / shell_mass[1:]
    return rate


def unmodeled_mass_loss_stop_reason(
    captured_mass: float,
    first_shell_mass: float,
    gnc_owned_mass: float,
    maximum_shell_fraction: float,
    maximum_gnc_fraction: float,
) -> str | None:
    """Return the first failed mass-conservation gate for a fixed-mass run."""

    values = (
        captured_mass,
        first_shell_mass,
        gnc_owned_mass,
        maximum_shell_fraction,
        maximum_gnc_fraction,
    )
    if any(not np.isfinite(value) for value in values):
        raise ValueError("mass-loss gate inputs must be finite")
    if (
        captured_mass < 0.0
        or first_shell_mass <= 0.0
        or gnc_owned_mass <= 0.0
        or not 0.0 < maximum_shell_fraction < 1.0
        or not 0.0 < maximum_gnc_fraction < 1.0
    ):
        raise ValueError("mass-loss gate inputs must be physical")
    if captured_mass / first_shell_mass >= maximum_shell_fraction:
        return "unmodeled_fluid_mass_loss_limit"
    if captured_mass / gnc_owned_mass >= maximum_gnc_fraction:
        return "unmodeled_gnc_mass_loss_limit"
    return None
