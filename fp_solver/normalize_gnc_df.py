#!/usr/bin/env python3
"""Generate a physically normalized GNC cusp DF and matched Monte-Carlo samples.

GNC evolves a dimensionless phase-space distribution g(x,j), with

    x = E / (G M_bh / R_h),       j = J / J_c(E),

and reconstructs number density as

    n(r) = (2/sqrt(pi)) n0 integral gbar(x) sqrt(R_h/r-x) dx.

The bridge need not be a single power law.  This script therefore solves the
GNC density integral for a non-negative tabulated isotropic DF.  A smooth
outer-reservoir continuation is fixed below the fitted energy range, and a
regularised non-negative least-squares inversion recovers the complete bridged
density profile.  The table is then rescaled at one shared physical radius,
checked against density and enclosed mass, and sampled with the native GNC
measure p(log x,log j) ~ g(x) x^-3/2 j^2.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable, Tuple

import numpy as np
from scipy.integrate import quad
from scipy.optimize import lsq_linear

G_PC_KMS2_MSUN = 4.30091e-3
C_KMS = 2.99792458e5


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonicalize_metadata(value, significant_digits: int = 10):
    """Round diagnostic floats so harmless solver jitter cannot change hashes."""
    if isinstance(value, dict):
        return {
            key: canonicalize_metadata(item, significant_digits)
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [canonicalize_metadata(item, significant_digits) for item in value]
    if isinstance(value, (float, np.floating)):
        numeric = float(value)
        if not math.isfinite(numeric):
            return numeric
        if abs(numeric) < 1.0e-14:
            return 0.0
        return float(f"{numeric:.{significant_digits - 1}e}")
    return value


def load_profile(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    d = np.loadtxt(path, comments="#", ndmin=2)
    if d.shape[1] < 3:
        raise ValueError("bridged profile must contain r, rho, sigma")
    r, rho, sig = (np.asarray(d[:, k], float) for k in range(3))
    good = np.isfinite(r) & np.isfinite(rho) & np.isfinite(sig)
    good &= (r > 0.0) & (rho > 0.0) & (sig > 0.0)
    order = np.argsort(r[good])
    return r[good][order], rho[good][order], sig[good][order]


def log_interp(r: np.ndarray, y: np.ndarray, rq: float) -> float:
    if not (r[0] <= rq <= r[-1]):
        raise ValueError(f"requested radius {rq:g} pc is outside profile")
    return float(np.exp(np.interp(np.log(rq), np.log(r), np.log(y))))


def shape_integral(psi: float, xmin: float, xmax: float, xb: float, power: float) -> float:
    upper = min(psi, xmax)
    if upper <= xmin:
        return 0.0
    val, err = quad(
        lambda x: (x / xb) ** power * math.sqrt(max(psi - x, 0.0)),
        xmin,
        upper,
        epsabs=0.0,
        epsrel=3.0e-10,
        points=[upper] if upper < psi * (1.0 - 1.0e-14) else None,
        limit=300,
    )
    return 2.0 / math.sqrt(math.pi) * val


def reconstructed_density(
    r: np.ndarray,
    rh: float,
    n0: float,
    amplitude: float,
    xmin: float,
    xmax: float,
    xb: float,
    power: float,
    particle_mass: float,
) -> np.ndarray:
    out = np.empty_like(r)
    for i, ri in enumerate(r):
        out[i] = (
            particle_mass
            * n0
            * amplitude
            * shape_integral(rh / ri, xmin, xmax, xb, power)
        )
    return out


def df_forward_matrix(
    psi: np.ndarray, xgrid: np.ndarray, quadrature_order: int = 12
) -> np.ndarray:
    """Map a tabulated isotropic GNC DF to n(psi)/n0.

    The table is interpolated linearly in log(x), matching the logarithmic GNC
    energy grid.  Gauss--Legendre quadrature integrates each intersected grid
    interval, including an interval cut by the local escape energy ``psi``.
    """

    psi = np.atleast_1d(np.asarray(psi, float))
    xgrid = np.asarray(xgrid, float)
    if np.any(psi <= 0.0) or np.any(xgrid <= 0.0):
        raise ValueError("psi and xgrid must be positive")
    if np.any(np.diff(xgrid) <= 0.0):
        raise ValueError("xgrid must increase strictly")
    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
    lx = np.log(xgrid)
    matrix = np.zeros((psi.size, xgrid.size))
    prefactor = 2.0 / math.sqrt(math.pi)
    for i, psii in enumerate(psi):
        for k in range(xgrid.size - 1):
            lo = xgrid[k]
            hi = min(xgrid[k + 1], psii)
            if hi <= lo:
                continue
            xx = 0.5 * (hi - lo) * nodes + 0.5 * (hi + lo)
            ww = 0.5 * (hi - lo) * weights
            t = (np.log(xx) - lx[k]) / (lx[k + 1] - lx[k])
            kernel = prefactor * np.sqrt(np.maximum(psii - xx, 0.0)) * ww
            matrix[i, k] += float(np.sum(kernel * (1.0 - t)))
            matrix[i, k + 1] += float(np.sum(kernel * t))
    return matrix


def table_reconstructed_density(
    radius: np.ndarray,
    rh: float,
    n0: float,
    xgrid: np.ndarray,
    gx: np.ndarray,
    particle_mass: float,
) -> np.ndarray:
    matrix = df_forward_matrix(rh / np.asarray(radius, float), xgrid)
    return particle_mass * n0 * (matrix @ np.asarray(gx, float))


def local_log_slope(r: np.ndarray, rho: np.ndarray, rq: float) -> float:
    """Return -d log(rho)/d log(r) from a local polynomial fit."""

    lr = np.log(r)
    i = int(np.searchsorted(r, rq))
    lo = max(0, i - 4)
    hi = min(r.size, i + 4)
    if hi - lo < 4:
        raise ValueError("too few profile points for a local density slope")
    return float(-np.polyfit(lr[lo:hi], np.log(rho[lo:hi]), 1)[0])


def construct_tabulated_isotropic_df(
    r: np.ndarray,
    rho: np.ndarray,
    rh: float,
    n0: float,
    particle_mass: float,
    xgrid: np.ndarray,
    normalization_radius: float,
    fit_outer_radius: float,
    rg: float,
    regularization: float,
    fit_points: int,
) -> tuple[np.ndarray, dict]:
    """Invert the bridged density into a non-negative isotropic Kepler DF."""

    if not (r[0] <= fit_outer_radius <= r[-1]):
        raise ValueError("DF fit outer radius lies outside the bridged profile")
    if normalization_radius > fit_outer_radius:
        raise ValueError("normalization radius lies outside the fitted DF domain")
    if regularization <= 0.0 or fit_points < 4 * xgrid.size:
        raise ValueError("DF inversion needs positive regularization and adequate fit points")

    gamma_outer = local_log_slope(r, rho, fit_outer_radius)
    if not (0.5 < gamma_outer < 1.5):
        raise ValueError(
            "the Kepler reservoir continuation requires 1/2 < local density "
            f"slope < 3/2, found {gamma_outer:g}"
        )
    xfit_requested = rh / (2.0 * fit_outer_radius)
    kfit = int(np.searchsorted(xgrid, xfit_requested))
    if not (1 <= kfit < xgrid.size - 3):
        raise ValueError("DF outer-reservoir join is too close to an energy-grid edge")
    xfit = float(xgrid[kfit])

    # Lock the underconstrained weak-binding reservoir to the local positive
    # Kepler continuation.  Above xfit every tabulated value is solved for.
    gbase = np.ones_like(xgrid)
    outer_power = gamma_outer - 1.5
    gbase[: kfit + 1] = (xgrid[: kfit + 1] / xfit) ** outer_power
    transform = np.zeros((xgrid.size, 1 + xgrid.size - kfit - 1))
    transform[: kfit + 1, 0] = 1.0
    for k in range(kfit + 1, xgrid.size):
        transform[k, 1 + k - (kfit + 1)] = 1.0

    fit_inner_radius = max(float(r[0]) * (1.0 + 1.0e-10),
                           6.0 * rg * (1.0 + 1.0e-8),
                           rh / (0.8 * float(xgrid[-1])))
    if fit_inner_radius >= fit_outer_radius:
        raise ValueError("empty radial interval for DF inversion")
    rfit = np.geomspace(fit_inner_radius, fit_outer_radius, fit_points)
    target = np.array([log_interp(r, rho, ri) for ri in rfit])
    target_dimensionless = target / (particle_mass * n0)
    forward = df_forward_matrix(rh / rfit, xgrid)
    design = forward @ (gbase[:, None] * transform)

    # Relative residuals give every decade in density comparable influence.
    row_weight = 1.0 / target_dimensionless
    d1 = np.zeros((xgrid.size - 1, xgrid.size))
    d2 = np.zeros((xgrid.size - 2, xgrid.size))
    for k in range(xgrid.size - 1):
        d1[k, k:k + 2] = (-1.0, 1.0)
    for k in range(xgrid.size - 2):
        d2[k, k:k + 3] = (1.0, -2.0, 1.0)
    system = np.vstack([
        design * row_weight[:, None],
        regularization * (d2 @ transform),
        0.1 * regularization * (d1 @ transform),
    ])
    rhs = np.concatenate([
        target_dimensionless * row_weight,
        np.zeros(d2.shape[0] + d1.shape[0]),
    ])
    solution = lsq_linear(
        system,
        rhs,
        bounds=(0.0, np.inf),
        tol=1.0e-11,
        lsmr_tol=1.0e-11,
        max_iter=2000,
    )
    if solution.status <= 0:
        raise RuntimeError(f"non-negative DF inversion failed: {solution.message}")
    gx = gbase * (transform @ solution.x)
    if np.any(~np.isfinite(gx)) or np.any(gx <= 0.0):
        raise RuntimeError("DF inversion produced a non-positive table")

    # The inversion fixes the shape.  One shared physical anchor fixes its
    # amplitude exactly and therefore preserves a common reservoir across a
    # boundary-position scan.
    norm_forward = df_forward_matrix(
        np.array([rh / normalization_radius]), xgrid
    )[0]
    target_norm = log_interp(r, rho, normalization_radius) / (
        particle_mass * n0
    )
    normalization_factor = target_norm / float(norm_forward @ gx)
    gx *= normalization_factor

    # Dense least-squares backends can differ at the last few bits even for
    # identical inputs.  Canonicalize the solved shape before the final physical
    # normalization.  Eight significant digits change the fitted DF by far less
    # than the declared density-roundtrip tolerance, while making the common
    # reservoir byte-reproducible across reporting-surface runs.
    gx_before_canonicalization = gx.copy()
    gx = np.asarray([float(f"{value:.7e}") for value in gx], dtype=float)
    canonicalization_change = float(np.max(np.abs(
        gx / gx_before_canonicalization - 1.0
    )))
    canonical_normalization = target_norm / float(norm_forward @ gx)
    gx *= canonical_normalization
    normalization_factor *= canonical_normalization
    # Use the exact values that will be emitted by ``write_df_table`` for all
    # subsequent density checks and Monte-Carlo sampling.  Otherwise two arrays
    # that serialize to the same table can retain sub-ulp in-memory differences
    # and generate different importance weights.
    gx = np.asarray([float(f"{value:.12e}") for value in gx], dtype=float)

    fitted = particle_mass * n0 * (forward @ gx)
    relative = fitted / target - 1.0
    diagnostics = {
        "df_shape_model": "regularized non-negative Abel inversion",
        "df_fit_outer_radius_pc": fit_outer_radius,
        "df_fit_inner_radius_pc": fit_inner_radius,
        "df_fit_points": fit_points,
        "df_fit_regularization": regularization,
        "df_outer_join_x": xfit,
        "df_outer_join_grid_index": kfit,
        "df_outer_density_slope": gamma_outer,
        "df_outer_power": outer_power,
        "df_solver_status": int(solution.status),
        "df_solver_message": str(solution.message),
        "df_solver_cost": float(solution.cost),
        "df_amplitude_normalization_factor": float(normalization_factor),
        "df_reproducibility_significant_digits": 8,
        "df_reproducibility_max_abs_relative_change": float(
            f"{canonicalization_change:.2e}"
        ),
        "df_fit_max_abs_relative_error": float(np.max(np.abs(relative))),
        "df_fit_median_abs_relative_error": float(np.median(np.abs(relative))),
        "df_table_min": float(np.min(gx)),
        "df_table_max": float(np.max(gx)),
    }
    if diagnostics["df_fit_max_abs_relative_error"] > 0.01:
        raise RuntimeError(
            "tabulated DF misses the bridged density by more than one per cent: "
            f"{diagnostics['df_fit_max_abs_relative_error']:g}"
        )
    return gx, diagnostics


def tabulated_log_pdf(
    xgrid: np.ndarray, gx: np.ndarray, exponent: float = -1.5,
    ngrid: int = 16384,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a normalized PDF and CDF in z=log(x)."""

    z = np.linspace(math.log(float(xgrid[0])), math.log(float(xgrid[-1])), ngrid)
    # Match ``df_forward_matrix`` and the native table lookup: g is linear on
    # the logarithmic energy coordinate.
    g = np.interp(z, np.log(xgrid), gx)
    raw = g * np.exp(exponent * z)
    dz = np.diff(z)
    cdf = np.zeros_like(z)
    cdf[1:] = np.cumsum(0.5 * (raw[:-1] + raw[1:]) * dz)
    norm = float(cdf[-1])
    if not np.isfinite(norm) or norm <= 0.0:
        raise RuntimeError("tabulated x PDF has invalid normalization")
    return z, raw / norm, cdf / norm


def stratified_tabulated_samples(
    n: int, z: np.ndarray, cdf: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    u = (np.arange(n, dtype=float) + rng.random(n)) / n
    rng.shuffle(u)
    return np.exp(np.interp(u, cdf, z))


def integrate_mass(r: np.ndarray, rho: np.ndarray) -> float:
    integrand = 4.0 * np.pi * rho * r**2
    return float(np.sum(0.5 * (integrand[:-1] + integrand[1:]) * np.diff(r)))


def stratified_power_samples(
    n: int,
    lo: float,
    hi: float,
    exponent_in_log_measure: float,
    rng: np.random.Generator,
) -> np.ndarray:
    u = (np.arange(n, dtype=float) + rng.random(n)) / n
    rng.shuffle(u)
    k = exponent_in_log_measure
    if abs(k) < 1.0e-12:
        return np.exp(np.log(lo) + u * np.log(hi / lo))
    return (lo**k + u * (hi**k - lo**k)) ** (1.0 / k)


def write_df_table(path: Path, x: np.ndarray, j: np.ndarray, g: np.ndarray, metadata: dict) -> None:
    metadata = canonicalize_metadata(metadata)
    with path.open("w") as f:
        # The cluster's load_df_from_table routine consumes exactly one header
        # record before reading ``nx nj``.  Keep the metadata on that same record.
        f.write("# physically normalized GNC DF table; metadata="
                + json.dumps(metadata, sort_keys=True) + "\n")
        f.write(f"{x.size:d} {j.size:d}\n")
        f.write(" ".join(f"{v:.12e}" for v in x) + "\n")
        f.write(" ".join(f"{v:.12e}" for v in j) + "\n")
        for row in g:
            f.write(" ".join(f"{v:.12e}" for v in row) + "\n")


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", type=Path, help="hydrostatically bridged physical profile")
    parser.add_argument("--bridge-json", type=Path, default=None)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--mbh", type=float, default=4.0e6)
    parser.add_argument("--rh", type=float, default=None)
    parser.add_argument("--r-boundary", type=float, default=None)
    parser.add_argument(
        "--normalization-radius",
        type=float,
        default=None,
        help=(
            "radius [pc] at which the physical DF amplitude is anchored; "
            "defaults to r-boundary.  Use one shared value when comparing "
            "several admissible handoff radii."
        ),
    )
    parser.add_argument(
        "--proposal-boundary-radius",
        type=float,
        default=None,
        help=(
            "boundary radius [pc] used only to focus the importance proposal; "
            "defaults to r-boundary.  A shared value gives identical samples "
            "across a handoff-invariance scan."
        ),
    )
    parser.add_argument("--beta", type=float, default=0.75)
    parser.add_argument(
        "--df-fit-outer-radius",
        type=float,
        default=None,
        help=(
            "outer radius [pc] used to invert the bridged density; defaults to "
            "R_h/2 and must enclose every compared handoff radius"
        ),
    )
    parser.add_argument("--df-fit-points", type=int, default=256)
    parser.add_argument("--df-regularization", type=float, default=1.0e-2)
    parser.add_argument("--particle-mass", type=float, default=1.0, help="GNC mass packet unit [Msun]")
    parser.add_argument("--guard-factor", type=float, default=100.0)
    parser.add_argument("--emin-floor", type=float, default=0.03)
    parser.add_argument("--emax", type=float, default=1.0e5)
    parser.add_argument("--jmin", type=float, default=5.0e-4)
    parser.add_argument("--jmax", type=float, default=0.99999)
    parser.add_argument("--grid-bins", type=int, default=48)
    parser.add_argument("--ranks", type=int, default=12)
    parser.add_argument("--samples-per-rank", type=int, default=5000)
    parser.add_argument(
        "--importance-uniform-fraction",
        type=float,
        default=0.2,
        help="fraction sampled uniformly in log x; weights correct back to the physical DF",
    )
    parser.add_argument(
        "--importance-capture-fraction",
        type=float,
        default=0.65,
        help="fraction sampled uniformly over the capture-dominant low-binding x interval",
    )
    parser.add_argument(
        "--importance-capture-xmax-factor",
        type=float,
        default=5.0,
        help="upper edge of focused x proposal in units of x_boundary",
    )
    parser.add_argument(
        "--importance-uniform-j-fraction",
        type=float,
        default=0.2,
        help="fraction sampled uniformly in log j to resolve radial and loss-cone orbits",
    )
    parser.add_argument(
        "--importance-capture-j-fraction",
        type=float,
        default=0.65,
        help="fraction sampled over the angular-momentum interval dominating late captures",
    )
    parser.add_argument("--importance-capture-jmin", type=float, default=0.015)
    parser.add_argument("--importance-capture-jmax", type=float, default=0.30)
    parser.add_argument("--weight-n", type=float, default=800.0)
    parser.add_argument("--clone-factor", type=int, default=30)
    parser.add_argument("--seed", type=int, default=73421)
    parser.add_argument(
        "--boundary-roundtrip-tolerance",
        type=float,
        default=None,
        help=(
            "maximum fractional density mismatch at r-boundary.  The default "
            "is 2e-8 when the DF is anchored there and 0.05 for a distinct "
            "shared normalization anchor."
        ),
    )
    args = parser.parse_args(argv)

    bridge = {}
    if args.bridge_json is not None:
        bridge = json.loads(args.bridge_json.read_text())
        if bridge.get("status") != "OK":
            raise RuntimeError("hydrostatic bridge did not pass its acceptance gate")
        if bridge.get("schema") != "absolute-closure-hydrostatic-bridge-v3":
            raise RuntimeError("hydrostatic bridge metadata schema is stale")
        if bridge.get("output_profile_sha256") != sha256_file(args.profile):
            raise RuntimeError("DF profile does not match the hydrostatic bridge output")
    rh = float(args.rh if args.rh is not None else bridge.get("rh_pc", math.nan))
    rb = float(
        args.r_boundary
        if args.r_boundary is not None
        else bridge.get("r_in_pc", math.nan)
    )
    rnorm = float(args.normalization_radius if args.normalization_radius is not None else rb)
    rproposal = float(
        args.proposal_boundary_radius
        if args.proposal_boundary_radius is not None
        else rb
    )
    if not (
        np.isfinite(rh)
        and rh > 0.0
        and np.isfinite(rb)
        and rb > 0.0
        and np.isfinite(rnorm)
        and rnorm > 0.0
        and np.isfinite(rproposal)
        and rproposal > 0.0
    ):
        raise ValueError("positive --rh and --r-boundary (or --bridge-json) are required")
    if bridge:
        for key, measured, expected in (
            ("mbh_msun", float(bridge["mbh_msun"]), args.mbh),
            ("rh_pc", float(bridge["rh_pc"]), rh),
        ):
            if not math.isclose(measured, expected, rel_tol=1.0e-10, abs_tol=0.0):
                raise RuntimeError(f"DF normalization metadata mismatch for {key}")
    rg = G_PC_KMS2_MSUN * args.mbh / C_KMS**2
    if max(rb, rnorm, rproposal) >= rh:
        raise ValueError(
            "the boundary, normalization anchor, and proposal anchor must all "
            "lie inside R_h for this Kepler-potential GNC implementation"
        )
    if min(rb, rnorm, rproposal) <= 6.0 * rg:
        raise ValueError("a requested matching radius is at or inside the Schwarzschild ISCO")
    if args.guard_factor <= 1.0:
        raise ValueError("guard-factor must exceed one")

    r, rho, sig = load_profile(args.profile)
    rho_b = log_interp(r, rho, rb)
    sig_b = log_interp(r, sig, rb)
    rho_norm = log_interp(r, rho, rnorm)
    rho_h = log_interp(r, rho, rh)
    n0 = rho_h / args.particle_mass

    # This factor of two is exact in GNC: a = R_h/(2x).
    xb = rh / (2.0 * rb)
    xnorm = rh / (2.0 * rnorm)
    xproposal = rh / (2.0 * rproposal)
    # Weakly bound eccentric orbits cross the spatial matching surface and
    # contribute to its density.  They must remain in the numerical reservoir even
    # when the formal GNC boundary is at much larger x.  Keeping the standard low-x
    # floor also gives the constructor many populated guard cells below x_b.
    xmin = min(args.emin_floor, xb / args.guard_factor)
    xmax = max(args.emax, 20.0 * xb)
    if not xmin < xb < xmax:
        raise ValueError(f"bad energy ordering: xmin={xmin:g}, xb={xb:g}, xmax={xmax:g}")
    xgrid = np.geomspace(xmin, xmax, args.grid_bins)
    jgrid = np.geomspace(args.jmin, args.jmax, args.grid_bins)
    fit_outer_radius = float(
        args.df_fit_outer_radius if args.df_fit_outer_radius is not None else 0.5 * rh
    )
    if fit_outer_radius <= max(rb, rnorm, rproposal):
        raise ValueError(
            "DF fit outer radius must enclose the boundary, normalization, and "
            "proposal radii; pass one shared larger value for a boundary scan"
        )
    gx, df_fit = construct_tabulated_isotropic_df(
        r,
        rho,
        rh,
        n0,
        args.particle_mass,
        xgrid,
        rnorm,
        fit_outer_radius,
        rg,
        args.df_regularization,
        args.df_fit_points,
    )
    asymp_boundary = float(np.interp(
        math.log(xb), np.log(xgrid), gx
    ))
    asymp_norm = float(np.interp(
        math.log(xnorm), np.log(xgrid), gx
    ))
    beta = local_log_slope(r, rho, rnorm)

    # The table stores the physical dimensionless GNC DF.  Its angular
    # dependence is constant for an isotropic cusp; j enters through the
    # phase-space measure used for sampling and for bar-g integration.
    gxy = np.repeat(gx[:, None], args.grid_bins, axis=1)

    args.outdir.mkdir(parents=True, exist_ok=True)
    table_meta = {
        "schema": "gnc-normalized-df-v6",
        "profile_sha256": sha256_file(args.profile),
        "bridge_json_sha256": (
            sha256_file(args.bridge_json) if args.bridge_json is not None else None
        ),
        "mbh_msun": args.mbh,
        "collision_kernel": bridge.get("collision_kernel"),
        "sigma0_over_m_cm2_g": bridge.get("sigma_over_m_cm2_g"),
        "w_kms": bridge.get("yukawa_w_kms"),
        "rh_pc": rh,
        "r_boundary_pc": rb,
        "x_boundary": xb,
        "normalization_radius_pc": rnorm,
        "normalization_x": xnorm,
        "normalization_density_msun_pc3": rho_norm,
        "proposal_boundary_radius_pc": rproposal,
        "proposal_x_boundary": xproposal,
        "xmin": xmin,
        "xmax": xmax,
        "jmin": args.jmin,
        "jmax": args.jmax,
        "beta": beta,
        "power": df_fit["df_outer_power"],
        "n0_pc3": n0,
        "rho_h_msun_pc3": rho_h,
        "rho_boundary_msun_pc3": rho_b,
        "sigma_boundary_kms": sig_b,
        "particle_mass_msun": args.particle_mass,
        "capture_species": "SBH point-particle proxy",
        "capture_radius_model": "GNC compact-object plunge, r_lc=16 Rg/(1+e)",
        "asymp_g_at_normalization_x": asymp_norm,
        "asymp_g_at_x_boundary": asymp_boundary,
        **df_fit,
    }
    table_meta = canonicalize_metadata(table_meta)
    write_df_table(args.outdir / "f_evolved_normalized.txt", xgrid, jgrid, gxy, table_meta)
    df_table_sha256 = sha256_file(args.outdir / "f_evolved_normalized.txt")
    asymp_boundary_emitted = float(table_meta["asymp_g_at_x_boundary"])

    # GNC estimates g from weighted samples.  In logarithmic coordinates the
    # physical number measure is g x^-3/2 j^2.  Direct sampling would put almost
    # every orbit in the weakly bound reservoir when x_b/xmin is large, leaving no
    # loss-cone statistics.  We therefore mix the physical law with broad and
    # capture-focused proposals, and store p/q as a per-particle weight_N
    # multiplier.  This is importance sampling, not a change to the DF.
    target_z, target_pdf_z, target_cdf_z = tabulated_log_pdf(xgrid, gx)
    mix = float(args.importance_uniform_fraction)
    mix_cap = float(args.importance_capture_fraction)
    if not (0.0 <= mix < 1.0 and 0.0 <= mix_cap < 1.0 and mix + mix_cap < 1.0):
        raise ValueError("x proposal fractions must be non-negative and sum to less than one")
    capture_x_lo = xproposal
    capture_x_hi = min(xmax, xproposal * float(args.importance_capture_xmax_factor))
    if not capture_x_lo < capture_x_hi:
        raise ValueError("focused capture x interval is empty")

    def target_pdf_logx(xv: np.ndarray) -> np.ndarray:
        return np.interp(np.log(xv), target_z, target_pdf_z)

    uniform_pdf_logx = 1.0 / math.log(xmax / xmin)
    # Generate all MPI shards before normalizing the weights.  Normalizing each
    # rank independently conditions the physical DF on a separate finite-sample
    # fluctuation on every rank.  GNC later divides the combined particle weights
    # by the MPI task count, so only one normalization over the full ensemble is
    # appropriate.
    rank_samples: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for rank in range(args.ranks):
        rng = np.random.default_rng(args.seed + rank)
        component_x = rng.choice(
            3, size=args.samples_per_rank,
            p=[1.0 - mix - mix_cap, mix, mix_cap],
        )
        xs = np.empty(args.samples_per_rank)
        select_target = component_x == 0
        select_uniform = component_x == 1
        select_capture = component_x == 2
        xs[select_target] = stratified_tabulated_samples(
            int(np.count_nonzero(select_target)), target_z, target_cdf_z, rng
        )
        xs[select_uniform] = stratified_power_samples(
            int(np.count_nonzero(select_uniform)), xmin, xmax, 0.0, rng
        )
        xs[select_capture] = stratified_power_samples(
            int(np.count_nonzero(select_capture)), capture_x_lo, capture_x_hi, 0.0, rng
        )
        pt = target_pdf_logx(xs)
        focused_pdf_logx = np.where(
            (xs >= capture_x_lo) & (xs <= capture_x_hi),
            1.0 / math.log(capture_x_hi / capture_x_lo),
            0.0,
        )
        proposal = (
            (1.0 - mix - mix_cap) * pt
            + mix * uniform_pdf_logx
            + mix_cap * focused_pdf_logx
        )
        importance_x = pt / proposal

        # The isotropic physical measure is dN/dln(j) proportional to j^2.
        # Drawing it directly leaves too few eccentric orbits to represent the
        # density at r_in or the loss cone.  A log-uniform mixture resolves those
        # orbits, while p/q preserves the original isotropic DF exactly.
        mix_j = float(args.importance_uniform_j_fraction)
        mix_cap_j = float(args.importance_capture_j_fraction)
        if not (0.0 <= mix_j < 1.0 and 0.0 <= mix_cap_j < 1.0
                and mix_j + mix_cap_j < 1.0):
            raise ValueError("j proposal fractions must be non-negative and sum to less than one")
        capture_j_lo = max(args.jmin, float(args.importance_capture_jmin))
        capture_j_hi = min(args.jmax, float(args.importance_capture_jmax))
        if not capture_j_lo < capture_j_hi:
            raise ValueError("focused capture j interval is empty")
        component_j = rng.choice(
            3, size=args.samples_per_rank,
            p=[1.0 - mix_j - mix_cap_j, mix_j, mix_cap_j],
        )
        js = np.empty(args.samples_per_rank)
        select_target_j = component_j == 0
        select_uniform_j = component_j == 1
        select_capture_j = component_j == 2
        js[select_target_j] = stratified_power_samples(
            int(np.count_nonzero(select_target_j)), args.jmin, args.jmax, 2.0, rng
        )
        js[select_uniform_j] = stratified_power_samples(
            int(np.count_nonzero(select_uniform_j)), args.jmin, args.jmax, 0.0, rng
        )
        js[select_capture_j] = stratified_power_samples(
            int(np.count_nonzero(select_capture_j)), capture_j_lo, capture_j_hi, 0.0, rng
        )
        target_pdf_logj = 2.0 * js**2 / (args.jmax**2 - args.jmin**2)
        uniform_pdf_logj = 1.0 / math.log(args.jmax / args.jmin)
        focused_pdf_logj = np.where(
            (js >= capture_j_lo) & (js <= capture_j_hi),
            1.0 / math.log(capture_j_hi / capture_j_lo),
            0.0,
        )
        proposal_j = (
            (1.0 - mix_j - mix_cap_j) * target_pdf_logj
            + mix_j * uniform_pdf_logj
            + mix_cap_j * focused_pdf_logj
        )
        importance = importance_x * target_pdf_logj / proposal_j
        rank_samples.append((xs, js, importance))

    raw_weight_sum = math.fsum(
        float(wi) for _, _, weights in rank_samples for wi in weights
    )
    total_samples = args.ranks * args.samples_per_rank
    global_raw_mean = raw_weight_sum / total_samples
    if not np.isfinite(global_raw_mean) or global_raw_mean <= 0.0:
        raise RuntimeError("global importance-weight mean is not positive and finite")

    # Normalize once, then canonicalize the finite-sample weights before they
    # are written.  BLAS reductions can differ at the final bit across repeated
    # inversions even when the canonical DF table is identical.  Rounding the
    # already-normalized weights to twelve significant digits is far below the
    # Monte-Carlo precision, and a second global normalization restores unit
    # mean without rank-local conditioning.
    provisional_rank_samples: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for rank, (xs, js, raw_importance) in enumerate(rank_samples):
        importance = raw_importance / global_raw_mean
        importance = np.asarray(
            [float(f"{value:.11e}") for value in importance], dtype=float
        )
        provisional_rank_samples.append((xs, js, importance))
    canonical_weight_mean = math.fsum(
        float(weight)
        for _, _, weights in provisional_rank_samples
        for weight in weights
    ) / total_samples
    normalized_rank_samples: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for rank, (xs, js, importance) in enumerate(provisional_rank_samples):
        importance = importance / canonical_weight_mean
        normalized_rank_samples.append((xs, js, importance))
        p = args.outdir / f"initial_xj_rank{rank:04d}.dat"
        with p.open("w") as f:
            f.write(f"{args.samples_per_rank:d}\n")
            for xv, jv, wi in zip(xs, js, importance):
                f.write(f"{xv:.16e} {jv:.16e} {wi:.16e}\n")

    # Physical normalization enters native GNC through dso%asymp and
    # dms%weight_asym.  The patched sampler count is deliberately independent of
    # this amplitude, so normalization changes physical weights, not resolution.
    (args.outdir / "mfrac.normalized.in").write_text(
        "#mass bin(given) = number GIVEN\n"
        "1 GIVEN\n"
        "# m1 mc m2 asymptotic_g(xb) Weight_n Clone_factor alpha_ini\n"
        f"1.0 1.0 1.0 {asymp_boundary_emitted:.16e} {args.weight_n:.8g} {args.clone_factor:d} 0.25\n"
        # SIDM particles are point-like at the loss cone.  The MS component would
        # be tidally disrupted at a stellar tidal radius and is therefore not a
        # physical proxy.  GNC's SBH component uses its relativistic plunge rule.
        "0.0 1.0 0.0 0.0 0.0\n"
    )
    (args.outdir / "profile_norm.in").write_text(
        f"{float(table_meta['rh_pc']):.16e} {float(table_meta['n0_pc3']):.16e} "
        f"{float(table_meta['r_boundary_pc']):.16e} "
        f"{float(table_meta['particle_mass_msun']):.16e} "
        f"{float(table_meta['beta']):.16e} "
        f"{float(table_meta['rho_boundary_msun_pc3']):.16e} "
        f"{float(table_meta['sigma_boundary_kms']):.16e} "
        f"{asymp_boundary_emitted:.16e}\n"
        "# rh_pc n0_pc^-3 r_boundary_pc particle_mass_Msun beta rho_b sigma_b asymp\n"
    )

    # Round-trip checks.  The boundary is an exact amplitude constraint.  The
    # comparison over the outermost factor 20 in radius avoids interpreting the
    # finite emax cutoff as a physical inner density truncation.
    r_eval_lo = max(6.0 * rg, rb / 20.0)
    r_eval = np.geomspace(r_eval_lo, rb, 300)
    rho_target = np.array([log_interp(r, rho, ri) for ri in r_eval])
    rho_back = table_reconstructed_density(
        r_eval, rh, n0, xgrid, gx, args.particle_mass
    )
    rel = rho_back / rho_target - 1.0
    m_target = integrate_mass(r_eval, rho_target)
    m_back = integrate_mass(r_eval, rho_back)
    rho_boundary_back = float(rho_back[-1])

    # Test the emitted law over the complete MPI ensemble.  Rank-local summaries
    # are retained only to expose an accidental shard imbalance; they do not alter
    # the physical normalization.
    x_sample = np.concatenate([x for x, _, _ in normalized_rank_samples])
    j_sample = np.concatenate([j for _, j, _ in normalized_rank_samples])
    w_sample = np.concatenate([w for _, _, w in normalized_rank_samples])
    rank_weight_means = [float(np.mean(w)) for _, _, w in normalized_rank_samples]
    global_weight_mean = float(math.fsum(float(wi) for wi in w_sample) / w_sample.size)
    active_fraction = float(np.mean(x_sample >= xb))
    active_weight_fraction = float(np.sum(w_sample[x_sample >= xb]) / np.sum(w_sample))
    effective_sample_size = float(np.sum(w_sample) ** 2 / np.sum(w_sample**2))
    radius_ratio = x_sample / xb  # r_b/a for a=R_h/(2x)
    reaches_boundary = (
        (radius_ratio <= 2.0)
        & (j_sample**2 <= 2.0 * radius_ratio - radius_ratio**2)
    )

    diag = {
        **table_meta,
        "df_table_sha256": df_table_sha256,
        "boundary_density_integral_I": float(
            df_forward_matrix(np.array([rh / rb]), xgrid)[0] @ gx
        ),
        "rho_boundary_roundtrip_msun_pc3": rho_boundary_back,
        "boundary_density_relative_error": rho_boundary_back / rho_b - 1.0,
        "roundtrip_radius_range_pc": [float(r_eval[0]), float(r_eval[-1])],
        "roundtrip_density_max_abs_relative_error": float(np.max(np.abs(rel))),
        "roundtrip_density_median_abs_relative_error": float(np.median(np.abs(rel))),
        "roundtrip_mass_target_msun": m_target,
        "roundtrip_mass_gnc_msun": m_back,
        "roundtrip_mass_relative_error": m_back / m_target - 1.0,
        "sample_x_range": [float(x_sample.min()), float(x_sample.max())],
        "sample_j_range": [float(j_sample.min()), float(j_sample.max())],
        "sample_active_fraction_x_ge_xb": active_fraction,
        "sample_active_physical_weight_fraction_x_ge_xb": active_weight_fraction,
        "sample_importance_weight_range": [float(w_sample.min()), float(w_sample.max())],
        "sample_effective_size_global": effective_sample_size,
        "importance_weight_raw_global_mean": global_raw_mean,
        "importance_weight_canonical_pre_normalization_mean": canonical_weight_mean,
        "importance_weight_reproducibility_significant_digits": 12,
        "importance_weight_global_mean": global_weight_mean,
        "importance_weight_rank_means": rank_weight_means,
        "importance_weight_rank_mean_range": [
            float(min(rank_weight_means)),
            float(max(rank_weight_means)),
        ],
        "importance_weight_normalization": "one global mean over all MPI ranks",
        "importance_uniform_fraction": mix,
        "importance_capture_fraction": mix_cap,
        "importance_capture_x_interval": [capture_x_lo, capture_x_hi],
        "importance_uniform_j_fraction": float(args.importance_uniform_j_fraction),
        "importance_capture_j_fraction": float(args.importance_capture_j_fraction),
        "importance_capture_j_interval": [capture_j_lo, capture_j_hi],
        "sample_reaches_boundary_fraction": float(np.mean(reaches_boundary)),
        "sample_reaches_boundary_physical_weight_fraction": float(
            np.sum(w_sample[reaches_boundary]) / np.sum(w_sample)
        ),
        "r_boundary_over_rh": rb / rh,
        "normalization_radius_over_rh": rnorm / rh,
        "proposal_boundary_radius_over_rh": rproposal / rh,
        "r_boundary_over_risco": rb / (6.0 * rg),
        "n_ranks": args.ranks,
        "samples_per_rank": args.samples_per_rank,
        "total_initial_samples": total_samples,
        "normalization_path": "bridged density -> non-negative Abel inversion -> shared physical DF anchor -> g(x_boundary) -> dms weight_asym -> sample weight_real",
        "status": "PASS",
    }
    tol_boundary = (
        float(args.boundary_roundtrip_tolerance)
        if args.boundary_roundtrip_tolerance is not None
        else (2.0e-8 if math.isclose(rnorm, rb, rel_tol=1.0e-14) else 0.05)
    )
    if not (np.isfinite(tol_boundary) and tol_boundary > 0.0):
        raise ValueError("boundary-roundtrip-tolerance must be positive and finite")
    diag["boundary_roundtrip_tolerance"] = tol_boundary
    tol_mass = 0.08
    if abs(diag["boundary_density_relative_error"]) > tol_boundary:
        diag["status"] = "FAIL_BOUNDARY_NORMALIZATION"
    elif abs(diag["roundtrip_mass_relative_error"]) > tol_mass:
        diag["status"] = "FAIL_MASS_ROUNDTRIP"
    elif abs(global_weight_mean - 1.0) > 5.0e-13:
        diag["status"] = "FAIL_GLOBAL_IMPORTANCE_NORMALIZATION"
    elif active_fraction < 0.01:
        diag["status"] = "FAIL_ACTIVE_SAMPLE_COUNT"

    diag = canonicalize_metadata(diag)
    (args.outdir / "df_normalization.json").write_text(
        json.dumps(diag, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(diag, indent=2, sort_keys=True))
    return 0 if diag["status"] == "PASS" else 3


if __name__ == "__main__":
    raise SystemExit(main())
