#!/usr/bin/env python3
"""Apply an accepted current or an injection ceiling to the SIDM fluid.

The input profile is a physical fluid snapshot (r, shell-averaged density, and
one-dimensional velocity dispersion).  It is reconstructed exactly in the native
Lagrangian variables of GravothermalSIDM.  The snapshot is remapped adiabatically
into hydrostatic equilibrium in the central point-mass potential.  The production
calculation evolves two otherwise identical branches: no inner current and the
thermal sink carried into the unresolved capture domain.  A positive full-return
ceiling remains available only as an explicitly provisional sensitivity test.

GNC owns the unresolved nuclear region. An accepted closure, or an explicitly
provisional ceiling, enters the first fluid control volume through the conservative
inner-face luminosity L_in. No local fluid mass is deleted. In the code's sign
convention, positive L_in deposits energy and negative L_in removes it.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import time
from pathlib import Path

import astropy.units as u
import numpy as np

from fluid_engine import load_fluid_engine

evolve, record = load_fluid_engine()

from fluid_bh_scaleheight import (
    bh_corrected_lmfp_inverse,
    conservative_specific_energy_rate,
    harmonic_effective_conductivity,
    shell_luminosity,
    unmodeled_mass_loss_stop_reason,
)
from interface_bridge import validate_interface_bridge


INNER_DENSITY_DEFINITION = (
    "mean density inside the innermost resolved Lagrangian shell, "
    "3*M(<r_inner)/(4*pi*r_inner^3)"
)
INNER_DISPERSION_DEFINITION = (
    "one-dimensional velocity dispersion in the innermost resolved "
    "Lagrangian shell"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class MeasuredBoundaryHalo(evolve.Halo):
    """Gravothermal halo with a conservative luminosity at its inner face."""

    L_inner = 0.0
    Mdot_msun_per_myr = 0.0

    def update_derived_parameters(self) -> None:
        """Update the fluid state using the BH-limited LMFP scale height.

        The stock gravothermal closure takes the long-mean-free-path scale
        height to be the self-gravitating Jeans length.  In a point-mass
        potential the available orbital scale is smaller.  Starting from the
        stock conductivity, multiply its LMFP branch by

            H_eff^2/H_J^2 = rho/(rho + M_bh/(r^2+r_soft^2)^(3/2)).

        This tends to unity in the halo and to the orbital-scale closure in the
        black-hole dominated region.  The SMFP branch is local and unchanged.
        """

        super().update_derived_parameters()
        mbh = float(getattr(self, "M_bh", 0.0))
        if mbh > 0.0:
            corrected_lmfp, factor = bh_corrected_lmfp_inverse(
                self.Kinv_lmfp,
                self.r,
                self.rho,
                mbh,
                self.r_soft,
            )
        else:
            factor = np.ones_like(self.rho)
            corrected_lmfp = self.Kinv_lmfp.copy()
        if np.any(~np.isfinite(factor)) or np.any(factor <= 0.0) or np.any(factor > 1.0):
            raise FloatingPointError("invalid black-hole LMFP scale-height factor")
        self.lmfp_scaleheight_factor = factor
        self.Kinv_lmfp = corrected_lmfp
        self.Kn = self.Kn / np.sqrt(factor)

        keff = harmonic_effective_conductivity(
            self.Kinv_smfp, self.Kinv_lmfp
        )
        self.L[:] = shell_luminosity(self.r, self.u, keff)

    def set_measured_boundary(self, luminosity_code: float,
                              mdot_msun_per_myr: float) -> None:
        self.L_inner = float(luminosity_code)
        self.Mdot_msun_per_myr = float(mdot_msun_per_myr)
        self.E_boundary_code = 0.0
        self.E_conduction_budget_error_code = 0.0
        self.E_conduction_budget_scale_code = 0.0

    def _specific_energy_rate(self) -> np.ndarray:
        """Return du/dt for all shells, including the inner-face current."""
        return conservative_specific_energy_rate(
            self.m, self.L, self.L_inner
        )

    def get_timestep(self) -> float:
        rate = self._specific_energy_rate()
        good = np.isfinite(rate) & (np.abs(rate) > 0.0)
        delta_t1 = float(np.min(np.abs(self.u[good] / rate[good]))) if np.any(good) else math.inf
        delta_t2 = float(np.min(1.0 / (self.rho * self.v)))
        if self.flag_timestep_use_relaxation and self.flag_timestep_use_energy:
            delta_t = min(delta_t1, delta_t2)
        elif self.flag_timestep_use_relaxation:
            delta_t = delta_t2
        elif self.flag_timestep_use_energy:
            delta_t = delta_t1
        else:
            raise RuntimeError("at least one timestep criterion is required")
        return self.t_epsilon * delta_t

    def conduct_heat(self) -> None:
        self.n_conduction += 1
        delta_t = self.get_timestep()
        rate = self._specific_energy_rate()
        self.delta_uc[:] = delta_t * rate

        # Finite-volume identity: sum dm du = (L_inner-L_outer) dt, and
        # L_outer is fixed to zero by update_derived_parameters().
        dm = np.diff(np.r_[0.0, self.m])
        actual = float(np.dot(dm, self.delta_uc))
        expected = self.L_inner * delta_t
        self.E_conduction_budget_error_code += actual - expected
        self.E_conduction_budget_scale_code += (
            float(np.dot(dm, np.abs(self.delta_uc))) + abs(expected)
        )

        delta_pc = self.p * self.delta_uc / self.u
        self.p += delta_pc
        self.u += self.delta_uc
        if np.any(self.p <= 0.0) or np.any(self.u <= 0.0):
            raise FloatingPointError("non-positive pressure or specific energy")

        self.E_boundary_code += expected
        self.t_before = self.t
        self.t += delta_t


def load_profile(path: Path) -> np.ndarray:
    a = np.loadtxt(path)
    if a.ndim != 2 or a.shape[1] < 3:
        raise ValueError(f"{path} must contain r_pc, rho_msun_pc3, sigma_kms")
    if np.any(~np.isfinite(a[:, :3])) or np.any(a[:, :3] <= 0.0):
        raise ValueError("profile contains non-finite or non-positive values")
    if np.any(np.diff(a[:, 0]) <= 0.0):
        raise ValueError("profile radii are not strictly increasing")
    return a[:, :3]


def install_physical_profile(h: MeasuredBoundaryHalo, profile: np.ndarray) -> dict:
    if profile.shape[0] != h.n_shells:
        raise ValueError(f"profile has {profile.shape[0]} shells; solver expects {h.n_shells}")
    rscale = h.scale_r.to_value(u.pc)
    rhoscale = h.scale_rho.to_value(u.Msun / u.pc**3)
    vscale = h.scale_v.to_value(u.km / u.s)
    mscale = h.scale_m.to_value(u.Msun)

    r = profile[:, 0] / rscale
    rho = profile[:, 1] / rhoscale
    sig = profile[:, 2] / vscale
    p = rho * sig**2
    dm = rho * np.diff(np.r_[0.0, r**3]) / 3.0
    m = np.cumsum(dm)

    h.r[:] = r
    h.rho[:] = rho
    h.p[:] = p
    h.m[:] = m
    h.t = 0.0
    h.t_before = 0.0
    h.n_adjustment = 0
    h.n_conduction = 0
    h.n_save = 0
    h.rho_center = float(rho[0])
    h.update_derived_parameters()

    roundtrip = np.column_stack([
        h.r * rscale,
        h.rho * rhoscale,
        np.sqrt(h.p / h.rho) * vscale,
    ])
    rel = np.max(np.abs(roundtrip / profile - 1.0))
    if rel > 5.0e-12:
        raise RuntimeError(f"physical profile round trip failed: {rel:g}")
    return {
        "profile_roundtrip_max_relative_error": float(rel),
        "inner_shell_mass_msun": float(h.m[0] * mscale),
        "total_enclosed_mass_msun": float(h.m[-1] * mscale),
        "r0_pc": float(profile[0, 0]),
        "rho0_msun_pc3": float(profile[0, 1]),
        "sigma0_kms": float(profile[0, 2]),
    }


def adiabatically_grow_black_hole(
    h: MeasuredBoundaryHalo,
    target_msun: float,
    max_force_fraction: float,
    r_epsilon: float,
    max_growth_steps: int,
) -> dict:
    """Grow the point mass without an impulsive change in the resolved force."""

    if target_msun < 0.0 or not 0.0 < max_force_fraction <= 0.25:
        raise ValueError("invalid black-hole growth controls")
    target_code = (target_msun * u.Msun / h.scale_m).decompose().value
    h.M_bh = 0.0
    h.r_epsilon = r_epsilon
    h.update_derived_parameters()
    shell_mass0 = h.m.copy()
    r0_before = float(h.r[0] * h.scale_r.to_value(u.pc))
    rho0_before = float(h.rho[0] * h.scale_rho.to_value(u.Msun / u.pc**3))
    entropy_raw = h.p / h.rho**h.gamma
    trajectory = []

    # Remove interpolation-level imbalance at zero BH mass before changing the
    # potential.  Hydrostatic adjustment is adiabatic and preserves shell mass.
    h.hydrostatic_adjustment()
    h.update_derived_parameters()
    entropy0 = h.p / h.rho**h.gamma
    initial_entropy_drift = float(
        np.max(np.abs(entropy0 / entropy_raw - 1.0))
    )
    r0_equilibrated = float(h.r[0] * h.scale_r.to_value(u.pc))
    rho0_equilibrated = float(
        h.rho[0] * h.scale_rho.to_value(u.Msun / u.pc**3)
    )
    for step in range(max_growth_steps):
        if h.M_bh >= target_code * (1.0 - 1.0e-14):
            break
        soften = h.r[0]**3 / (h.r[0]**2 + h.r_soft**2) ** 1.5
        enclosed0 = h.m[0] + h.M_bh * soften
        dmbh = max_force_fraction * enclosed0 / max(soften, 1.0e-300)
        old_mbh = h.M_bh
        h.M_bh = min(target_code, h.M_bh + dmbh)
        h.hydrostatic_adjustment()
        h.update_derived_parameters()
        arrays = np.r_[h.r, h.rho, h.p, h.u]
        if np.any(~np.isfinite(arrays)) or np.any(arrays <= 0.0):
            raise FloatingPointError("non-physical state during adiabatic BH growth")
        trajectory.append({
            "step": step + 1,
            "M_bh_msun": float(h.M_bh * h.scale_m.to_value(u.Msun)),
            "fractional_force_increment_at_r0": float(
                (h.M_bh - old_mbh) * soften / enclosed0
            ),
            "r0_pc": float(h.r[0] * h.scale_r.to_value(u.pc)),
            "rho0_msun_pc3": float(
                h.rho[0] * h.scale_rho.to_value(u.Msun / u.pc**3)
            ),
            "hydrostatic_last_max_dr_over_r": float(
                np.max(np.abs(h.delta_r / h.r))
            ),
            "hydrostatic_adjustments": int(h.n_adjustment),
        })
    else:
        raise RuntimeError("black-hole growth exceeded max_growth_steps")

    entropy = h.p / h.rho**h.gamma
    entropy_drift = float(np.max(np.abs(entropy / entropy0 - 1.0)))
    mass_drift = float(np.max(np.abs(h.m - shell_mass0)))
    final_residual = float(np.max(np.abs(h.delta_r / h.r)))
    if mass_drift > 1.0e-13 * max(float(h.m[-1]), 1.0):
        raise RuntimeError("Lagrangian shell masses changed during BH remap")
    if final_residual > max(10.0 * r_epsilon, 1.0e-8):
        raise RuntimeError(f"BH hydrostatic remap did not converge: {final_residual:g}")
    if initial_entropy_drift > 1.0e-2:
        raise RuntimeError(
            "zero-BH profile equilibration changes entropy by too much: "
            f"{initial_entropy_drift:g}"
        )
    if entropy_drift > 5.0e-3:
        raise RuntimeError(f"BH remap entropy drift is too large: {entropy_drift:g}")
    return {
        "method": "adaptive adiabatic point-mass growth",
        "M_bh_target_msun": target_msun,
        "M_bh_final_msun": float(h.M_bh * h.scale_m.to_value(u.Msun)),
        "growth_steps": len(trajectory),
        "max_force_fraction": max_force_fraction,
        "initial_equilibration_entropy_max_relative_drift": initial_entropy_drift,
        "entropy_max_relative_drift": entropy_drift,
        "shell_mass_max_absolute_drift_code": mass_drift,
        "hydrostatic_final_max_dr_over_r": final_residual,
        "r0_initial_pc": r0_before,
        "rho0_initial_msun_pc3": rho0_before,
        "r0_pre_bh_equilibrated_pc": r0_equilibrated,
        "rho0_pre_bh_equilibrated_msun_pc3": rho0_equilibrated,
        "r0_final_pc": float(h.r[0] * h.scale_r.to_value(u.pc)),
        "rho0_final_msun_pc3": float(
            h.rho[0] * h.scale_rho.to_value(u.Msun / u.pc**3)
        ),
        "trajectory": trajectory,
    }


def validate_precomputed_black_hole_state(
    h: MeasuredBoundaryHalo,
    target_msun: float,
    r_epsilon: float,
) -> dict:
    """Verify a profile already remapped into the target BH potential."""

    target_code = (target_msun * u.Msun / h.scale_m).decompose().value
    h.M_bh = target_code
    h.r_epsilon = r_epsilon
    h.update_derived_parameters()
    entropy0 = h.p / h.rho**h.gamma
    mass0 = h.m.copy()
    r0_before = float(h.r[0] * h.scale_r.to_value(u.pc))
    rho0_before = float(h.rho[0] * h.scale_rho.to_value(u.Msun / u.pc**3))
    h.hydrostatic_adjustment()
    h.update_derived_parameters()
    entropy = h.p / h.rho**h.gamma
    entropy_drift = float(np.max(np.abs(entropy / entropy0 - 1.0)))
    mass_drift = float(np.max(np.abs(h.m - mass0)))
    final_residual = float(np.max(np.abs(h.delta_r / h.r)))
    radius_drift = float(abs(
        h.r[0] * h.scale_r.to_value(u.pc) / r0_before - 1.0
    ))
    if mass_drift > 1.0e-13 * max(float(h.m[-1]), 1.0):
        raise RuntimeError("precomputed BH profile changed Lagrangian shell masses")
    if entropy_drift > 5.0e-4:
        raise RuntimeError(
            f"precomputed BH profile entropy drift is too large: {entropy_drift:g}"
        )
    if radius_drift > 5.0e-4:
        raise RuntimeError(
            f"precomputed BH profile is not hydrostatic: radius drift {radius_drift:g}"
        )
    if final_residual > max(10.0 * r_epsilon, 1.0e-8):
        raise RuntimeError(
            f"precomputed BH profile residual is too large: {final_residual:g}"
        )
    return {
        "method": "validated precomputed adiabatic point-mass remap",
        "M_bh_target_msun": target_msun,
        "M_bh_final_msun": float(h.M_bh * h.scale_m.to_value(u.Msun)),
        "growth_steps": 0,
        "max_force_fraction": None,
        "initial_equilibration_entropy_max_relative_drift": 0.0,
        "entropy_max_relative_drift": entropy_drift,
        "shell_mass_max_absolute_drift_code": mass_drift,
        "hydrostatic_final_max_dr_over_r": final_residual,
        "validation_inner_radius_relative_drift": radius_drift,
        "r0_initial_pc": r0_before,
        "rho0_initial_msun_pc3": rho0_before,
        "r0_pre_bh_equilibrated_pc": r0_before,
        "rho0_pre_bh_equilibrated_msun_pc3": rho0_before,
        "r0_final_pc": float(h.r[0] * h.scale_r.to_value(u.pc)),
        "rho0_final_msun_pc3": float(
            h.rho[0] * h.scale_rho.to_value(u.Msun / u.pc**3)
        ),
        "trajectory": [],
    }


def total_energy_code(h: MeasuredBoundaryHalo) -> float:
    dm = np.diff(np.r_[0.0, h.m])
    rin = np.r_[0.0, h.r[:-1]]
    rmid = np.cbrt(0.5 * (rin**3 + h.r**3))
    mmid = np.r_[0.0, h.m[:-1]] + 0.5 * dm
    kinetic = float(np.dot(dm, h.u))
    potential_self = float(-np.dot(dm, mmid / rmid))
    potential_bh = float(
        -np.dot(dm, h.M_bh / np.sqrt(rmid**2 + h.r_soft**2))
    )
    return kinetic + potential_self + potential_bh


def snapshot(h: MeasuredBoundaryHalo, branch: str, lphys: float,
             scale_t_myr: float, scale_e: float, mdot: float) -> dict:
    t_myr = h.t * scale_t_myr
    rho_inner_mean = (
        h.rho[0] * h.scale_rho.to_value(u.Msun / u.pc**3)
    )
    sigma_inner_1d = h.v[0] * h.scale_v.to_value(u.km / u.s)
    r_inner = h.r[0] * h.scale_r.to_value(u.pc)
    return {
        "branch": branch,
        "tau_relax": h.t / h.t_relax,
        "time_myr": t_myr,
        "rho_inner_mean_msun_pc3": rho_inner_mean,
        "sigma_inner_1d_kms": sigma_inner_1d,
        "r_inner_pc": r_inner,
        "inner_shell_mass_msun": h.m[0] * h.scale_m.to_value(u.Msun),
        # Retain the old column names for readers of development trajectories.
        # New analysis products require the canonical names above and verify
        # that these aliases are numerically identical.
        "rho_c_msun_pc3": rho_inner_mean,
        "sigma_c_kms": sigma_inner_1d,
        "r0_pc": r_inner,
        "M_bh_msun": h.M_bh * h.scale_m.to_value(u.Msun),
        "L_inner_msun_kms2_per_myr": lphys,
        "captured_mass_msun": mdot * t_myr,
        "E_boundary_msun_kms2": h.E_boundary_code * scale_e,
        "E_conduction_budget_error_msun_kms2": h.E_conduction_budget_error_code * scale_e,
        "E_total_code": total_energy_code(h),
        "n_conduction": h.n_conduction,
        "lmfp_scaleheight_factor_inner": float(h.lmfp_scaleheight_factor[0]),
        "lmfp_scaleheight_factor_min": float(np.min(h.lmfp_scaleheight_factor)),
    }


def run_branch(branch: str, lphys: float, mdot: float, profile: np.ndarray,
               args: argparse.Namespace, diag: dict,
               bridge: dict) -> tuple[list[dict], dict]:
    recdir = args.outdir / "records" / branch
    h = MeasuredBoundaryHalo(
        record.HaloRecord(str(recdir)),
        sigma_m_with_units=args.sigma_over_m,
        w_units=args.w_units,
        model_elastic_scattering_lmfp=(
            "YukawaBornViscosityApproxTchannel_K3_order1"
        ),
        model_elastic_scattering_smfp=(
            "YukawaBornViscosityApproxTchannel_K5_order1"
        ),
        M_bh_with_units=args.mbh,
        r_soft=0.0,
        n_shells=profile.shape[0],
        n_adjustment_max=args.bh_adjustment_max,
    )
    # Accept either the original fluid snapshot or the deterministic post-remap
    # snapshot used to construct the GNC reservoir.
    if args.profile_state == "pre-bh":
        h.M_bh = 0.0
    profile_check = install_physical_profile(h, profile)
    closure_radius = (
        diag.get("identity", {}).get("r_outer_pc")
        if isinstance(diag.get("identity"), dict)
        else None
    )
    if closure_radius is None:
        closure_radius = diag.get("r_boundary_pc")
    bridge_check = validate_interface_bridge(
        bridge,
        mbh_msun=args.mbh,
        sigma_over_m_cm2_g=args.sigma_over_m,
        w_kms=args.w_units,
        fluid_profile_sha256=sha256_file(args.profile),
        fluid_inner_shell_mass_msun=profile_check["inner_shell_mass_msun"],
        closure_radius_pc=(
            float(closure_radius) if closure_radius is not None else None
        ),
    )
    if args.profile_state == "pre-bh":
        bh_check = adiabatically_grow_black_hole(
            h,
            args.mbh,
            args.bh_max_force_fraction,
            args.bh_r_epsilon,
            args.bh_max_growth_steps,
        )
    else:
        bh_check = validate_precomputed_black_hole_state(
            h, args.mbh, args.bh_r_epsilon
        )
    l_quantity = lphys * u.Msun * (u.km / u.s)**2 / u.Myr
    lcode = (l_quantity / h.scale_L).decompose().value
    h.set_measured_boundary(lcode, mdot)
    h.t_epsilon = args.t_epsilon
    h.r_epsilon = args.r_epsilon
    scale_t_myr = h.scale_t.to_value(u.Myr)
    scale_e = (h.scale_L * h.scale_t).to_value(u.Msun * (u.km / u.s)**2)

    rows = [snapshot(h, branch, lphys, scale_t_myr, scale_e, mdot)]
    next_output = args.output_dtau
    rho_inner_initial = float(h.rho[0])
    wall0 = time.time()
    stop_reason = "tau_end"
    first_shell_mass = float(profile_check["inner_shell_mass_msun"])
    gnc_owned_mass = float(bridge_check["GNC_owned_mass_msun"])
    while h.t / h.t_relax < args.tau_end:
        if h.n_conduction >= args.max_steps:
            stop_reason = "max_steps"
            break
        if time.time() - wall0 > args.wall_seconds:
            stop_reason = "wall_limit"
            break
        h.conduct_heat()
        h.hydrostatic_adjustment()
        arrays = np.r_[h.r, h.rho, h.p, h.u, h.L]
        if np.any(~np.isfinite(arrays)) or np.any(h.r <= 0.0) or np.any(h.rho <= 0.0):
            raise FloatingPointError(f"{branch}: non-physical fluid state")
        tau = h.t / h.t_relax
        if tau >= next_output or tau >= args.tau_end:
            rows.append(snapshot(h, branch, lphys, scale_t_myr, scale_e, mdot))
            next_output += args.output_dtau
        if h.rho[0] / rho_inner_initial >= args.rho_factor_end:
            stop_reason = "rho_factor_end"
            if rows[-1]["n_conduction"] != h.n_conduction:
                rows.append(snapshot(h, branch, lphys, scale_t_myr, scale_e, mdot))
            break
        captured_mass = mdot * h.t * scale_t_myr
        mass_stop = unmodeled_mass_loss_stop_reason(
            captured_mass,
            first_shell_mass,
            gnc_owned_mass,
            args.max_unmodeled_fluid_shell_mass_fraction,
            args.max_unmodeled_gnc_mass_fraction,
        )
        if mass_stop is not None:
            stop_reason = mass_stop
            if rows[-1]["n_conduction"] != h.n_conduction:
                rows.append(snapshot(h, branch, lphys, scale_t_myr, scale_e, mdot))
            break

    # Always preserve the state at which the branch actually stopped.  This is
    # needed not only for physical stopping conditions, but also for a step or
    # wall-clock guard reached between scheduled output times.  The branch
    # summary and trajectory must describe the same endpoint.
    if rows[-1]["n_conduction"] != h.n_conduction:
        rows.append(snapshot(h, branch, lphys, scale_t_myr, scale_e, mdot))

    budget_rel = (
        abs(h.E_conduction_budget_error_code)
        / h.E_conduction_budget_scale_code
        if h.E_conduction_budget_scale_code > 1.0e-300 else 0.0
    )
    summary = {
        "branch": branch,
        "stop_reason": stop_reason,
        "final_tau_relax": h.t / h.t_relax,
        "final_time_myr": h.t * scale_t_myr,
        "final_rho_inner_mean_msun_pc3": rows[-1][
            "rho_inner_mean_msun_pc3"
        ],
        "rho_inner_mean_factor": float(h.rho[0] / rho_inner_initial),
        # Deprecated aliases retained in the raw branch summary only.
        "final_rho_c_msun_pc3": rows[-1]["rho_inner_mean_msun_pc3"],
        "rho_c_factor": float(h.rho[0] / rho_inner_initial),
        "n_conduction": h.n_conduction,
        "L_inner_code": lcode,
        "L_inner_msun_kms2_per_myr": lphys,
        "Mdot_msun_per_myr": mdot,
        "captured_mass_msun": mdot * h.t * scale_t_myr,
        "captured_mass_fraction_of_first_fluid_shell": (
            mdot * h.t * scale_t_myr / first_shell_mass
        ),
        "captured_mass_fraction_of_GNC_owned_domain": (
            mdot * h.t * scale_t_myr / gnc_owned_mass
        ),
        "captured_mass_fraction_of_black_hole": (
            mdot * h.t * scale_t_myr / args.mbh
        ),
        "maximum_unmodeled_fluid_shell_mass_fraction": (
            args.max_unmodeled_fluid_shell_mass_fraction
        ),
        "maximum_unmodeled_GNC_mass_fraction": (
            args.max_unmodeled_gnc_mass_fraction
        ),
        "conduction_energy_budget_error_code": h.E_conduction_budget_error_code,
        "conduction_energy_budget_scale_code": h.E_conduction_budget_scale_code,
        "conduction_energy_budget_relative_error": budget_rel,
        "collapse_proxy_time_myr": (
            h.t * scale_t_myr if stop_reason == "rho_factor_end" else None
        ),
        "wall_seconds": time.time() - wall0,
        **profile_check,
        "black_hole_remap": bh_check,
        "interface_bridge": bridge_check,
    }
    print(json.dumps(summary, sort_keys=True), flush=True)
    return rows, summary


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--profile", type=Path, required=True)
    p.add_argument("--diagnostics", type=Path, required=True)
    p.add_argument("--bridge-json", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--sigma-over-m", type=float, default=100.0)
    p.add_argument("--w-units", type=float, default=80.0)
    p.add_argument("--mbh", type=float, default=4.0e6)
    p.add_argument(
        "--profile-state", choices=("pre-bh", "post-bh"), default="pre-bh",
        help="whether the supplied snapshot precedes or follows the adiabatic BH remap",
    )
    p.add_argument("--bh-max-force-fraction", type=float, default=0.05)
    p.add_argument("--bh-r-epsilon", type=float, default=1.0e-10)
    p.add_argument("--bh-adjustment-max", type=int, default=500)
    p.add_argument("--bh-max-growth-steps", type=int, default=512)
    p.add_argument("--tau-end", type=float, default=3.0,
                   help="elapsed halo relaxation times after the matched snapshot")
    p.add_argument("--output-dtau", type=float, default=0.02)
    p.add_argument("--rho-factor-end", type=float, default=1.0e3)
    p.add_argument("--t-epsilon", type=float, default=1.0e-4)
    p.add_argument("--r-epsilon", type=float, default=1.0e-12)
    p.add_argument("--max-steps", type=int, default=2_000_000)
    p.add_argument("--wall-seconds", type=float, default=10_000.0)
    p.add_argument(
        "--max-unmodeled-fluid-shell-mass-fraction",
        type=float,
        default=0.02,
        help="stop before omitted capture changes the first fluid shell by two per cent",
    )
    p.add_argument(
        "--max-unmodeled-gnc-mass-fraction",
        type=float,
        default=0.20,
        help="stop before the fixed GNC reservoir loses twenty per cent of its mass",
    )
    p.add_argument(
        "--allow-provisional", action="store_true",
        help=("run a sensitivity calculation even when the input GNC current "
              "has not yet passed its stationarity/counting gates")
    )
    p.add_argument(
        "--branches", default="control,sink,source",
        help="comma-separated subset of control,sink,source",
    )
    args = p.parse_args()
    if not (
        0.0 < args.max_unmodeled_fluid_shell_mass_fraction < 1.0
        and 0.0 < args.max_unmodeled_gnc_mass_fraction < 1.0
    ):
        raise ValueError("unmodeled mass-loss fractions must lie between zero and one")
    args.outdir.mkdir(parents=True, exist_ok=True)

    profile = load_profile(args.profile)
    diag = json.loads(args.diagnostics.read_text())
    bridge = json.loads(args.bridge_json.read_text())
    input_status = diag.get("status")
    absolute_input = input_status == "ABSOLUTE_CLOSURE_MEASURED"
    mass_closure_input = input_status == "FLUID_MASS_CLOSURE_MEASURED"
    accepted_input = absolute_input or mass_closure_input
    provisional_input = (
        input_status == "FINITE_ANGLE_ISOTROPIC_INJECTION_CEILING_CONVERGED"
    )
    if accepted_input:
        expected_schema = (
            "gnc-fluid-absolute-closure-v1"
            if absolute_input else "gnc-fluid-mass-closure-v1"
        )
        if diag.get("schema") != expected_schema:
            raise RuntimeError("accepted closure metadata schema is stale")
        if mass_closure_input and diag.get(
            "capture_binding_energy_current_used_by_fluid"
        ) is not False:
            raise RuntimeError(
                "mass closure must exclude the capture binding-energy current"
            )
        closure_gates = diag.get("gates")
        if (
            not isinstance(closure_gates, dict)
            or not closure_gates
            or any(value is not True for value in closure_gates.values())
        ):
            raise RuntimeError("accepted closure contains a missing or failed gate")
        identity = diag.get("identity", {})
        if identity.get("fluid_profile_sha256") != sha256_file(args.profile):
            raise RuntimeError(
                "accepted closure and fluid response use different post-BH profiles"
            )
        if identity.get("bridge_json_sha256") != sha256_file(args.bridge_json):
            raise RuntimeError(
                "accepted closure and fluid response use different bridge records"
            )
        if identity.get("bridged_profile_sha256") != bridge.get(
            "output_profile_sha256"
        ):
            raise RuntimeError(
                "accepted closure does not use this hydrostatic bridge profile"
            )
        for key, expected in (
            ("mbh_msun", args.mbh),
            ("sigma0_over_m_cm2_g", args.sigma_over_m),
            ("w_kms", args.w_units),
        ):
            measured = float(identity.get(key, math.nan))
            if not math.isclose(
                measured, expected, rel_tol=1.0e-10, abs_tol=0.0
            ):
                raise RuntimeError(
                    f"accepted closure and fluid response differ in {key}"
                )
    else:
        if not args.allow_provisional:
            raise RuntimeError(
                "the input current did not pass the absolute-closure gates"
            )
        if not provisional_input:
            raise RuntimeError(
                "--allow-provisional accepts only a converged finite-angle "
                "isotropic-reservoir injection ceiling"
            )
        identity = diag.get("identity", {})
        for key, expected in (
            ("mbh_msun", args.mbh),
            ("sigma0_over_m_cm2_g", args.sigma_over_m),
            ("w_kms", args.w_units),
        ):
            measured = float(identity.get(key, math.nan))
            if not math.isclose(measured, expected, rel_tol=1.0e-10, abs_tol=0.0):
                raise RuntimeError(
                    f"finite-angle and fluid metadata differ for {key}"
                )
        if identity.get("bridge_json_sha256") != sha256_file(args.bridge_json):
            raise RuntimeError(
                "finite-angle current and fluid response use different bridge records"
            )
        if identity.get("profile_sha256") != bridge.get("output_profile_sha256"):
            raise RuntimeError(
                "finite-angle current does not use the bridged profile supplied to the fluid"
            )
    mdot = float(
        diag["measured_mdot_msun_per_myr"]
        if accepted_input else diag["isotropic_upper_mdot_msun_per_myr"]
    )
    currents = {
        "control": 0.0,
        "sink": float(diag["thermal_sink_current_msun_kms2_per_myr"]),
    }
    if not mass_closure_input:
        currents["source"] = float(
            diag["returned_energy_source_ceiling_msun_kms2_per_myr"]
        )
    if (
        not math.isfinite(mdot) or mdot <= 0.0
        or not math.isfinite(currents["sink"]) or currents["sink"] >= 0.0
        or (
            "source" in currents
            and (
                not math.isfinite(currents["source"])
                or currents["source"] <= 0.0
            )
        )
    ):
        raise RuntimeError("accepted closure currents have an unphysical sign")

    all_rows: list[dict] = []
    summaries = {}
    branches = tuple(x.strip() for x in args.branches.split(",") if x.strip())
    if not branches or any(x not in currents for x in branches):
        raise ValueError(
            "--branches requests a current that is absent from this closure"
        )
    for branch in branches:
        branch_mdot = 0.0 if branch == "control" else mdot
        rows, summary = run_branch(
            branch, currents[branch], branch_mdot, profile, args, diag, bridge
        )
        all_rows.extend(rows)
        summaries[branch] = summary

    fields = list(all_rows[0].keys())
    with (args.outdir / "fluid_feedback_trajectories.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(all_rows)

    out = {
        "schema": "gnc-fluid-feedback-v5",
        "status": (
            "MEASURED_FLUID_FEEDBACK_COMPLETE"
            if absolute_input
            else "MASS_CLOSURE_FLUID_RESPONSE_COMPLETE"
            if mass_closure_input
            else "FLUID_FEEDBACK_SENSITIVITY_COMPLETE"
        ),
        "input_closure_status": input_status,
        "input_closure_accepted": accepted_input,
        "absolute_two_current_closure": absolute_input,
        "mass_current_closure": mass_closure_input,
        "provisional_sensitivity_only": provisional_input,
        "method": (
            "shared adiabatic central point-mass remap; BH-limited LMFP scale height; "
            "matched t-channel Yukawa K3/K5 conductivity; fixed orbit-domain "
            "current at the matched deep snapshot; conservative inner-face "
            "luminosity; no bare fluid mass sink"
        ),
        "current_role": (
            "accepted absolute closure"
            if absolute_input
            else "accepted mass closure with unresolved binding-energy moment"
            if mass_closure_input
            else "provisional isotropic-reservoir injection-ceiling sensitivity"
        ),
        "closure_update_mode": "fixed_at_matched_snapshot",
        "resolved_diagnostics": {
            "density_field": "rho_inner_mean_msun_pc3",
            "density_definition": INNER_DENSITY_DEFINITION,
            "density_is_extrapolated_central_value": False,
            "dispersion_field": "sigma_inner_1d_kms",
            "dispersion_definition": INNER_DISPERSION_DEFINITION,
            "radius_field": "r_inner_pc",
        },
        "profile": args.profile.name,
        "gnc_diagnostics": args.diagnostics.name,
        "interface_bridge_json": args.bridge_json.name,
        "input_identity": {
            "closure_json_sha256": sha256_file(args.diagnostics),
            "fluid_profile_sha256": sha256_file(args.profile),
            "bridge_json_sha256": sha256_file(args.bridge_json),
            "bridged_profile_sha256": bridge.get("output_profile_sha256"),
        },
        "applied_currents": {
            "Mdot_msun_per_myr": mdot,
            "control_luminosity_msun_kms2_per_myr": currents["control"],
            "sink_luminosity_msun_kms2_per_myr": currents["sink"],
            "source_luminosity_msun_kms2_per_myr": currents.get("source"),
        },
        "tau_end_requested": args.tau_end,
        "sigma_over_m_cm2_g": args.sigma_over_m,
        "yukawa_w_kms": args.w_units,
        "M_bh_msun": args.mbh,
        "input_profile_state": args.profile_state,
        "branches": summaries,
        "trajectory_csv": "fluid_feedback_trajectories.csv",
    }
    (args.outdir / "fluid_feedback_summary.json").write_text(
        json.dumps(out, indent=2, sort_keys=True) + "\n"
    )
    print(out["status"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
