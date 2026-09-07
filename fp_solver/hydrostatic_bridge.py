#!/usr/bin/env python3
"""Build a hydrostatic SIDM-to-SMBH bridge and locate the inner N=1 crossing.

The gravothermal calculation need not be extended through a region in which its
fluid closure is invalid.  For a newly inserted black hole, the primary bridge
resolves the central Lagrangian shell as a hydrostatic, adiabatic atmosphere.  It
preserves that shell's physical mass and the pressure at its outer face.  A
steady weakly collisional conductive cusp remains available as a comparison,

    rho ~ r^-beta,  beta = (3 + a) / 4,

where sigma_DM ~ v^-a.  Pressure is not guessed: it is integrated inward from the
resolved fluid pressure using spherical hydrostatic balance in the combined SMBH
and bridged-DM potential.  The GNC/fluid interface is then the inner root of

    N(r) = 2 pi r rho(r) (sigma/m) = 1.

Input columns are r [pc], rho [Msun/pc^3], sigma_1d [km/s].  The first density
is interpreted correctly as a central-shell average; local edge conditions are
extrapolated from the thin shells outside it.  Extra columns are ignored.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np

from sidm_born_kernel import fluid_tchannel_kp, tchannel_viscosity_over_sigma0

G_PC_KMS2_MSUN = 4.30091e-3
C_KMS = 2.99792458e5
MSUN_G = 1.98847e33
PC_CM = 3.0856775814913673e18
CM2_G_TO_PC2_MSUN = MSUN_G / PC_CM**2


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def provenance_name(path: Path) -> str:
    """Return a root-independent label; the adjacent SHA identifies the file."""
    return path.name


def load_profile(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape[1] < 3:
        raise ValueError("profile must have at least three columns: r rho sigma")
    r, rho, sig = (np.asarray(data[:, k], float) for k in range(3))
    good = np.isfinite(r) & np.isfinite(rho) & np.isfinite(sig)
    good &= (r > 0.0) & (rho > 0.0) & (sig > 0.0)
    r, rho, sig = r[good], rho[good], sig[good]
    order = np.argsort(r)
    r, rho, sig = r[order], rho[order], sig[order]
    unique = np.r_[True, np.diff(r) > 0.0]
    r, rho, sig = r[unique], rho[unique], sig[unique]
    if r.size < 4:
        raise ValueError("profile needs at least four positive, distinct radii")
    return r, rho, sig


def fit_inner_slope(r: np.ndarray, rho: np.ndarray, nfit: int) -> float:
    n = min(max(3, nfit), r.size)
    slope = np.polyfit(np.log(r[:n]), np.log(rho[:n]), 1)[0]
    gamma = float(-slope)
    # A core-to-cusp bridge requires N to decrease inward before the cusp is
    # attached.  Report the raw fit, but prevent noisy shell values from creating
    # a divergent continuation.
    return float(np.clip(gamma, -0.25, 0.95))


def rho_piecewise(
    r: np.ndarray,
    r0: float,
    rho0: float,
    rh: float,
    gamma: float,
    beta: float,
    mode: str,
) -> Tuple[np.ndarray, float, float]:
    if mode == "core-only":
        return rho0 * (r / r0) ** (-gamma), r0, rho0
    if mode == "cusp-from-r0" or rh >= r0:
        rt = r0
        rho_t = rho0
        return rho0 * (r / r0) ** (-beta), rt, rho_t

    rt = rh
    rho_t = rho0 * (rt / r0) ** (-gamma)
    out = np.empty_like(r)
    cusp = r <= rt
    out[cusp] = rho_t * (r[cusp] / rt) ** (-beta)
    out[~cusp] = rho0 * (r[~cusp] / r0) ** (-gamma)
    return out, rt, rho_t


def build_yukawa_cusp(
    r: np.ndarray,
    r0: float,
    rho0: float,
    sig0: float,
    rh: float,
    gamma: float,
    mode: str,
    mbh: float,
    w_kms: float,
    max_iter: int = 200,
    tolerance: float = 2.0e-8,
) -> Tuple[np.ndarray, np.ndarray, float, float, int, float]:
    """Iterate a constant-luminosity Yukawa cusp and hydrostatic pressure.

    In the black-hole dominated LMFP region, the conductive luminosity scales
    as rho^2 K3(sigma/w) sigma^3 r^3.  Relative to the constant-cross-section
    cusp this gives rho proportional to r^-3/4 K3^-1/2.  The pressure and K3
    profiles are iterated together, while the resolved fluid pressure at r0 is
    held fixed.
    """

    if mode == "core-only":
        rho = rho0 * (r / r0) ** (-gamma)
        p, sig, _ = hydrostatic_pressure(r, rho, rho0 * sig0**2, mbh)
        return rho, sig, r0, rho0, 0, 0.0
    rt = r0 if mode == "cusp-from-r0" or rh >= r0 else rh
    rho_t = rho0 * (rt / r0) ** (-gamma)
    outer = r > rt
    rho = np.where(
        outer,
        rho0 * (r / r0) ** (-gamma),
        rho_t * (r / rt) ** (-0.75),
    )
    residual = math.inf
    for iteration in range(1, max_iter + 1):
        _, sig, _ = hydrostatic_pressure(r, rho, rho0 * sig0**2, mbh)
        k3 = np.asarray(fluid_tchannel_kp(sig / w_kms, 3), float)
        k3_t = float(np.exp(np.interp(np.log(rt), np.log(r), np.log(k3))))
        target = rho_t * (r / rt) ** (-0.75) * np.sqrt(k3_t / k3)
        target[outer] = rho0 * (r[outer] / r0) ** (-gamma)
        residual = float(np.max(np.abs(np.log(target / rho))))
        rho *= np.exp(0.45 * np.log(target / rho))
        if residual < tolerance:
            break
    else:
        raise RuntimeError(f"Yukawa cusp iteration failed to converge: {residual:g}")
    _, sig, _ = hydrostatic_pressure(r, rho, rho0 * sig0**2, mbh)
    return rho, sig, rt, rho_t, iteration, residual


def enclosed_mass_numeric(r: np.ndarray, rho: np.ndarray) -> np.ndarray:
    """Mass enclosed by the supplied inner continuation, with M(0)=0."""
    integrand = 4.0 * np.pi * rho * r**2
    m = np.zeros_like(r)
    # The missing interval [0,r[0]] is represented by the local power-law slope.
    if r.size > 1:
        slope = np.log(rho[1] / rho[0]) / np.log(r[1] / r[0])
    else:
        slope = -0.75
    denom = max(3.0 + slope, 1.0e-6)
    m[0] = 4.0 * np.pi * rho[0] * r[0] ** 3 / denom
    dr = np.diff(r)
    m[1:] = m[0] + np.cumsum(0.5 * (integrand[:-1] + integrand[1:]) * dr)
    return m


def hydrostatic_pressure(
    r: np.ndarray,
    rho: np.ndarray,
    p_outer: float,
    mbh: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    mdm = enclosed_mass_numeric(r, rho)
    force = G_PC_KMS2_MSUN * rho * (mbh + mdm) / r**2
    integ = np.zeros_like(r)
    for i in range(r.size - 2, -1, -1):
        integ[i] = integ[i + 1] + 0.5 * (force[i] + force[i + 1]) * (
            r[i + 1] - r[i]
        )
    pressure = p_outer + integ
    sigma = np.sqrt(np.maximum(pressure / rho, 0.0))
    return pressure, sigma, mdm


def extrapolate_inner_edge(
    r: np.ndarray, rho: np.ndarray, sig: np.ndarray, nfit: int = 5
) -> Tuple[float, float]:
    """Estimate local density and pressure just outside the first shell.

    The zeroth fluid density is an average over the central Lagrangian shell,
    not a point value at its outer face.  The thin shells immediately outside
    it provide the appropriate local boundary state for a subgrid bridge.
    """

    if r.size < 3:
        raise ValueError("at least three fluid shells are required")
    stop = min(r.size, 1 + max(2, nfit))
    x = np.log(r[1:stop])
    lrho = np.polyfit(x, np.log(rho[1:stop]), 1)
    lp = np.polyfit(x, np.log(rho[1:stop] * sig[1:stop] ** 2), 1)
    rho_edge = float(np.exp(np.polyval(lrho, np.log(r[0]))))
    pressure_edge = float(np.exp(np.polyval(lp, np.log(r[0]))))
    return rho_edge, pressure_edge


def build_mass_matched_adiabatic_bridge(
    r: np.ndarray,
    r0: float,
    shell_average_rho: float,
    rho_edge_outside: float,
    pressure_edge: float,
    mbh: float,
    gamma_eos: float = 5.0 / 3.0,
    tolerance: float = 2.0e-9,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, int, float, float]:
    """Resolve the central shell while preserving its mass and edge pressure.

    For each trial density on the inner side of the first shell boundary, the
    entropy constant K=P_edge/rho_edge^gamma is fixed and hydrostatic balance
    is iterated in the combined BH and reconstructed-DM potential.  A shooting
    solve then chooses the edge density whose integrated mass equals the
    original Lagrangian shell mass.
    """

    target_mass = 4.0 * np.pi * shell_average_rho * r0**3 / 3.0

    def solve_density(edge_density: float):
        entropy = pressure_edge / edge_density**gamma_eos
        coeff = (gamma_eos - 1.0) * G_PC_KMS2_MSUN * mbh / (
            gamma_eos * entropy
        )
        enthalpy = edge_density ** (gamma_eos - 1.0) + coeff * (
            1.0 / r - 1.0 / r0
        )
        rho_trial = np.maximum(enthalpy, np.finfo(float).tiny) ** (
            1.0 / (gamma_eos - 1.0)
        )
        residual = math.inf
        for iteration in range(1, 301):
            pressure, sigma, mdm = hydrostatic_pressure(
                r, rho_trial, pressure_edge, mbh
            )
            target = np.maximum(pressure / entropy, np.finfo(float).tiny) ** (
                1.0 / gamma_eos
            )
            target[-1] = edge_density
            residual = float(np.max(np.abs(np.log(target / rho_trial))))
            rho_trial *= np.exp(0.45 * np.log(target / rho_trial))
            if residual < tolerance:
                break
        else:
            raise RuntimeError(
                f"adiabatic bridge iteration failed to converge: {residual:g}"
            )
        pressure, sigma, mdm = hydrostatic_pressure(
            r, rho_trial, pressure_edge, mbh
        )
        return rho_trial, pressure, sigma, mdm, iteration, residual

    low = max(rho_edge_outside * 1.0e-4, np.finfo(float).tiny)
    low_state = solve_density(low)
    for _ in range(12):
        if low_state[3][-1] <= target_mass:
            break
        low *= 0.1
        low_state = solve_density(low)
    high = low
    high_state = low_state
    for _ in range(80):
        if high_state[3][-1] >= target_mass:
            break
        high *= 1.8
        high_state = solve_density(high)
    if not (low_state[3][-1] <= target_mass <= high_state[3][-1]):
        raise RuntimeError(
            "central-shell mass cannot be bracketed by the hydrostatic bridge: "
            f"target={target_mass:g}, range=({low_state[3][-1]:g},"
            f"{high_state[3][-1]:g})"
        )
    best = None
    for _ in range(80):
        mid = math.sqrt(low * high)
        state = solve_density(mid)
        best = (mid, state)
        if state[3][-1] < target_mass:
            low = mid
        else:
            high = mid
        if abs(state[3][-1] / target_mass - 1.0) < 2.0e-10:
            break
    edge_density, state = best
    rho_inner, pressure, sigma, mdm, iterations, residual = state
    mass_error = float(mdm[-1] / target_mass - 1.0)
    if abs(mass_error) > 1.0e-8:
        raise RuntimeError(f"adiabatic bridge mass mismatch: {mass_error:g}")
    return (
        rho_inner,
        pressure,
        sigma,
        float(edge_density),
        int(iterations),
        float(residual),
        mass_error,
    )


def crossings_log(r: np.ndarray, y: np.ndarray, level: float = 1.0) -> np.ndarray:
    z = np.log(y / level)
    roots = []
    for i in range(r.size - 1):
        if z[i] == 0.0:
            roots.append(r[i])
        elif z[i] * z[i + 1] < 0.0:
            lr = np.log(r[i]) + (np.log(r[i + 1]) - np.log(r[i])) * (
                -z[i] / (z[i + 1] - z[i])
            )
            roots.append(float(np.exp(lr)))
    if z[-1] == 0.0:
        roots.append(r[-1])
    return np.asarray(roots, float)


def build_one(
    r_res: np.ndarray,
    rho_res: np.ndarray,
    sig_res: np.ndarray,
    mbh: float,
    sigma_over_m: float,
    velocity_exponent: float,
    rh: float,
    gamma: float,
    beta: float,
    mode: str,
    n_inner: int,
    kernel: str,
    w_kms: float,
) -> Dict[str, np.ndarray | float | str]:
    r0, rho0, sig0 = map(float, (r_res[0], rho_res[0], sig_res[0]))
    rho_edge, pressure_edge = extrapolate_inner_edge(r_res, rho_res, sig_res)
    target_shell_mass = 4.0 * np.pi * rho0 * r0**3 / 3.0
    rg = G_PC_KMS2_MSUN * mbh / C_KMS**2
    r_floor = max(6.0 * rg, r0 * 1.0e-10)
    r_inner = np.geomspace(r_floor, r0, n_inner)
    cusp_iterations = 0
    cusp_residual = 0.0
    mass_relative_error = None
    if mode == "adiabatic-mass-matched":
        (
            rho_inner,
            pressure,
            sigma_inner,
            rho_transition,
            cusp_iterations,
            cusp_residual,
            mass_relative_error,
        ) = build_mass_matched_adiabatic_bridge(
            r_inner,
            r0,
            rho0,
            rho_edge,
            pressure_edge,
            mbh,
        )
        mdm = enclosed_mass_numeric(r_inner, rho_inner)
        r_transition = r0
    elif kernel == "yukawa-tchannel":
        (rho_inner, sigma_inner, r_transition, rho_transition,
         cusp_iterations, cusp_residual) = build_yukawa_cusp(
            r_inner, r0, rho0, sig0, rh, gamma, mode, mbh, w_kms
        )
        pressure, sigma_inner, mdm = hydrostatic_pressure(
            r_inner, rho_inner, rho0 * sig0**2, mbh
        )
    else:
        rho_inner, r_transition, rho_transition = rho_piecewise(
            r_inner, r0, rho0, rh, gamma, beta, mode
        )
        pressure, sigma_inner, mdm = hydrostatic_pressure(
            r_inner, rho_inner, rho0 * sig0**2, mbh
        )
    som0 = sigma_over_m * CM2_G_TO_PC2_MSUN
    if kernel == "yukawa-tchannel":
        # The angular-randomization criterion uses the exact viscosity moment.
        # Relative RMS speed between one particle and an isotropic background.
        vc = np.sqrt(G_PC_KMS2_MSUN * (mbh + mdm) / r_inner)
        vrel = np.sqrt(vc**2 + 3.0 * sigma_inner**2)
        som_inner = som0 * tchannel_viscosity_over_sigma0(vrel / w_kms)
        speed_factor_inner = vrel / vc
    else:
        som_inner = som0 * (sigma_inner / sig0) ** (-velocity_exponent)
        speed_factor_inner = np.ones_like(r_inner)
    n_inner_arr = 2.0 * np.pi * r_inner * rho_inner * som_inner * speed_factor_inner
    roots_inner = crossings_log(r_inner, n_inner_arr)

    # The zeroth fluid density is a shell average.  In mass-matched mode replace
    # that one row with the resolved inner profile and preserve all thin shells
    # outside its boundary exactly.
    if mode == "adiabatic-mass-matched":
        r_all = np.r_[r_inner, r_res[1:]]
        rho_all = np.r_[rho_inner, rho_res[1:]]
        sig_all = np.r_[sigma_inner, sig_res[1:]]
    else:
        r_all = np.r_[r_inner[:-1], r_res]
        rho_all = np.r_[rho_inner[:-1], rho_res]
        sig_all = np.r_[sigma_inner[:-1], sig_res]
    p_all = rho_all * sig_all**2
    if kernel == "yukawa-tchannel":
        # Outside the bridge, use the same local relative-speed estimate.
        mdm_all = enclosed_mass_numeric(r_all, rho_all)
        vc_all = np.sqrt(G_PC_KMS2_MSUN * (mbh + mdm_all) / r_all)
        vrel_all = np.sqrt(vc_all**2 + 3.0 * sig_all**2)
        som_all = som0 * tchannel_viscosity_over_sigma0(vrel_all / w_kms)
        speed_factor_all = vrel_all / vc_all
    else:
        som_all = som0 * (sig_all / sig0) ** (-velocity_exponent)
        speed_factor_all = np.ones_like(r_all)
    n_all = 2.0 * np.pi * r_all * rho_all * som_all * speed_factor_all
    roots_all = crossings_log(r_all, n_all)

    if roots_inner.size:
        r_in = float(roots_inner[0])
        rho_in = float(np.exp(np.interp(np.log(r_in), np.log(r_inner), np.log(rho_inner))))
        sig_in = float(np.interp(np.log(r_in), np.log(r_inner), sigma_inner))
    elif n_inner_arr[-1] <= 1.0:
        # The first fluid shell is already orbit dominated; the matching surface
        # lies in the resolved profile rather than in the bridge.
        r_in = float(roots_all[0]) if roots_all.size else float("nan")
        rho_in = float(np.exp(np.interp(np.log(r_in), np.log(r_all), np.log(rho_all)))) if np.isfinite(r_in) else float("nan")
        sig_in = float(np.interp(np.log(r_in), np.log(r_all), sig_all)) if np.isfinite(r_in) else float("nan")
    else:
        r_in = float("nan")
        rho_in = float("nan")
        sig_in = float("nan")
    mass_inside_r_in = (
        float(np.interp(r_in, r_inner, mdm)) if np.isfinite(r_in) else float("nan")
    )

    return {
        "mode": mode,
        "r": r_all,
        "rho": rho_all,
        "sigma": sig_all,
        "pressure": p_all,
        "N": n_all,
        "r_inner_grid": r_inner,
        "rho_inner_grid": rho_inner,
        "sigma_inner_grid": sigma_inner,
        "pressure_inner_grid": pressure,
        "mdm_inner_grid": mdm,
        "r_transition": float(r_transition),
        "rho_transition": float(rho_transition),
        "r_in": r_in,
        "rho_in": rho_in,
        "sigma_in": sig_in,
        "roots": roots_all,
        "kernel": kernel,
        "w_kms": w_kms,
        "cusp_iterations": cusp_iterations,
        "cusp_log_residual": cusp_residual,
        "target_first_shell_mass_msun": target_shell_mass,
        "reconstructed_first_shell_mass_msun": float(mdm[-1]),
        "first_shell_mass_relative_error": mass_relative_error,
        "mass_inside_r_in_msun": mass_inside_r_in,
        "mass_inside_r_in_fraction_of_first_shell": (
            mass_inside_r_in / target_shell_mass
            if np.isfinite(mass_inside_r_in) else float("nan")
        ),
        "volume_inside_r_in_fraction_of_first_shell": (
            (r_in / r0) ** 3 if np.isfinite(r_in) else float("nan")
        ),
        "rho_edge_outside_msun_pc3": rho_edge,
        "pressure_edge_outside": pressure_edge,
        "rho_edge_inside_msun_pc3": float(rho_inner[-1]),
        "pressure_edge_inside": float(pressure[-1]),
    }


def finite_or_none(x: float | None) -> float | None:
    return float(x) if x is not None and np.isfinite(x) else None


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", type=Path)
    parser.add_argument("--out-profile", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--mbh", type=float, default=4.0e6, help="SMBH mass [Msun]")
    parser.add_argument("--sigma-over-m", type=float, default=100.0, help="cm^2/g")
    parser.add_argument("--velocity-exponent", type=float, default=0.0, help="sigma_DM proportional to v^-a")
    parser.add_argument(
        "--kernel", choices=("constant", "powerlaw", "yukawa-tchannel"),
        default="constant", help="collision kernel used by the bridge and N criterion",
    )
    parser.add_argument("--w-kms", type=float, default=80.0, help="Yukawa velocity scale [km/s]")
    parser.add_argument("--rh", type=float, default=None, help="influence radius [pc]; default G Mbh/sigma0^2")
    parser.add_argument("--inner-slope", type=float, default=None, help="resolved core density slope gamma")
    parser.add_argument("--fit-shells", type=int, default=8)
    parser.add_argument("--n-inner", type=int, default=2400)
    parser.add_argument(
        "--max-subgrid-mass-fraction",
        type=float,
        default=0.02,
        help=(
            "maximum fraction of the first Lagrangian shell owned directly by "
            "GNC; below this gate the remainder is a controlled hydrostatic "
            "fluid bridge"
        ),
    )
    parser.add_argument(
        "--mode",
        choices=("adiabatic-mass-matched", "broken-at-rh", "cusp-from-r0", "core-only"),
        default="adiabatic-mass-matched",
    )
    args = parser.parse_args(argv)
    if not 0.0 < args.max_subgrid_mass_fraction < 1.0:
        raise ValueError("max subgrid mass fraction must lie between zero and one")

    r, rho, sig = load_profile(args.profile)
    gamma_fit = fit_inner_slope(r, rho, args.fit_shells)
    gamma = gamma_fit if args.inner_slope is None else float(args.inner_slope)
    beta = (3.0 + args.velocity_exponent) / 4.0
    if args.kernel != "yukawa-tchannel" and not (0.5 < beta < 3.0):
        raise ValueError(
            "this hydrostatic GNC bridge requires 1/2 < beta < 3 so the "
            "isotropic cusp DF is positive and has finite enclosed mass"
        )
    rh = (
        G_PC_KMS2_MSUN * args.mbh / float(sig[0]) ** 2
        if args.rh is None
        else float(args.rh)
    )

    primary = build_one(
        r,
        rho,
        sig,
        args.mbh,
        args.sigma_over_m,
        args.velocity_exponent,
        rh,
        gamma,
        beta,
        args.mode,
        args.n_inner,
        args.kernel,
        args.w_kms,
    )
    alternate_mode = "cusp-from-r0" if args.mode == "broken-at-rh" else "broken-at-rh"
    alternate = build_one(
        r,
        rho,
        sig,
        args.mbh,
        args.sigma_over_m,
        args.velocity_exponent,
        rh,
        gamma,
        beta,
        alternate_mode,
        args.n_inner,
        args.kernel,
        args.w_kms,
    )
    core_only = build_one(
        r,
        rho,
        sig,
        args.mbh,
        args.sigma_over_m,
        args.velocity_exponent,
        rh,
        gamma,
        beta,
        "core-only",
        args.n_inner,
        args.kernel,
        args.w_kms,
    )

    args.out_profile.parent.mkdir(parents=True, exist_ok=True)
    arr = np.column_stack(
        [
            primary["r"],
            primary["rho"],
            primary["sigma"],
            primary["pressure"],
            primary["N"],
        ]
    )
    np.savetxt(
        args.out_profile,
        arr,
        fmt="%.12e",
        header="r_pc rho_Msun_pc3 sigma1d_kms pressure_rho_sigma2 N_scatter_per_orbit",
    )

    rg = G_PC_KMS2_MSUN * args.mbh / C_KMS**2
    subgrid_mass_fraction = float(
        primary["mass_inside_r_in_fraction_of_first_shell"]
    )
    subgrid_volume_fraction = float(
        primary["volume_inside_r_in_fraction_of_first_shell"]
    )
    subgrid_gate = (
        np.isfinite(subgrid_mass_fraction)
        and subgrid_mass_fraction <= args.max_subgrid_mass_fraction
    )
    rin_values = np.asarray([primary["r_in"], alternate["r_in"]], float)
    rin_finite = rin_values[np.isfinite(rin_values)]
    diag = {
        "schema": "absolute-closure-hydrostatic-bridge-v3",
        "input_profile": provenance_name(args.profile),
        "input_profile_sha256": sha256_file(args.profile),
        "output_profile": provenance_name(args.out_profile),
        "output_profile_sha256": sha256_file(args.out_profile),
        "mbh_msun": args.mbh,
        "sigma_over_m_cm2_g": args.sigma_over_m,
        "sigma_over_m_pc2_msun_at_sigma0": args.sigma_over_m * CM2_G_TO_PC2_MSUN,
        "velocity_exponent_a": args.velocity_exponent,
        "collision_kernel": args.kernel,
        "yukawa_w_kms": args.w_kms if args.kernel == "yukawa-tchannel" else None,
        "cusp_slope_beta": beta,
        "core_slope_gamma_fit": gamma_fit,
        "core_slope_gamma_used": gamma,
        "r0_pc": float(r[0]),
        "rho0_msun_pc3": float(rho[0]),
        "sigma0_kms": float(sig[0]),
        "rh_pc": rh,
        "rg_pc": rg,
        "r_isco_pc": 6.0 * rg,
        "primary_mode": args.mode,
        "r_transition_pc": float(primary["r_transition"]),
        "rho_transition_msun_pc3": float(primary["rho_transition"]),
        "r_in_pc": finite_or_none(float(primary["r_in"])),
        "rho_in_msun_pc3": finite_or_none(float(primary["rho_in"])),
        "sigma_in_kms": finite_or_none(float(primary["sigma_in"])),
        "all_N1_crossings_pc": [float(x) for x in primary["roots"]],
        "alternate_mode": alternate_mode,
        "alternate_r_in_pc": finite_or_none(float(alternate["r_in"])),
        "unrelaxed_core_only_r_in_pc": finite_or_none(float(core_only["r_in"])),
        "bridge_rin_envelope_pc": (
            [float(np.min(rin_finite)), float(np.max(rin_finite))]
            if rin_finite.size
            else None
        ),
        "r_in_over_r0": finite_or_none(float(primary["r_in"]) / float(r[0])),
        "r_in_over_rh": finite_or_none(float(primary["r_in"]) / rh),
        "r_in_over_risco": finite_or_none(float(primary["r_in"]) / (6.0 * rg)),
        "N_at_first_fluid_shell": float(primary["N"][np.searchsorted(primary["r"], r[0])]),
        "cusp_iterations": int(primary["cusp_iterations"]),
        "cusp_log_residual": float(primary["cusp_log_residual"]),
        "target_first_shell_mass_msun": float(primary["target_first_shell_mass_msun"]),
        "reconstructed_first_shell_mass_msun": float(primary["reconstructed_first_shell_mass_msun"]),
        "first_shell_mass_relative_error": finite_or_none(primary["first_shell_mass_relative_error"]),
        "mass_inside_r_in_msun": finite_or_none(primary["mass_inside_r_in_msun"]),
        "mass_inside_r_in_fraction_of_first_shell": finite_or_none(subgrid_mass_fraction),
        "volume_inside_r_in_fraction_of_first_shell": finite_or_none(subgrid_volume_fraction),
        "maximum_allowed_subgrid_mass_fraction": args.max_subgrid_mass_fraction,
        "fluid_subgrid_bridge_gate": "PASS" if subgrid_gate else "FAIL",
        "rho_edge_outside_msun_pc3": float(primary["rho_edge_outside_msun_pc3"]),
        "rho_edge_inside_msun_pc3": float(primary["rho_edge_inside_msun_pc3"]),
        "pressure_edge_outside": float(primary["pressure_edge_outside"]),
        "pressure_edge_inside": float(primary["pressure_edge_inside"]),
        "status": (
            "OK"
            if np.isfinite(float(primary["r_in"])) and subgrid_gate
            else (
                "SUBGRID_MASS_FRACTION_TOO_LARGE"
                if np.isfinite(float(primary["r_in"]))
                else "NO_INNER_CROSSING"
            )
        ),
    }
    args.out_json.write_text(json.dumps(diag, indent=2, sort_keys=True) + "\n")
    print(json.dumps(diag, indent=2, sort_keys=True))
    return 0 if diag["status"] == "OK" else 2


if __name__ == "__main__":
    raise SystemExit(main())
