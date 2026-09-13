#!/usr/bin/env python3
"""Non-orbit-averaged spherical BGK transition solver.

The distribution is represented as ``f(r, v_r, J)`` with
``J=r v_t``.  This is exactly equivalent to ``f(r, v_r, v_t)`` at a fixed
radius, but angular momentum is the natural conservative coordinate for the
collisionless drift.  Keeping the sign of ``v_r`` allows the calculation to
represent a coherent inflow, which an orbit-averaged ``f(E,J)`` solver cannot.

The code advances

``d f/dt + v_r d f/dr + (J^2/r^3-GM/r^2) d f/dv_r = C_BGK[f]``

between an absorbing inner surface and a maintained outer Maxwellian.  The
BGK target is a discrete maximum-entropy Maxwellian with the same mass, radial
momentum, and kinetic energy as the local distribution.  Thus the collision
step conserves all three represented invariants to roundoff.  This model is a
controlled kinetic bridge, not a replacement for the finite-angle Yukawa
collision operator used by the production loss-cone calculation.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import time

import numpy as np

from bondi_hydrodynamic import critical_bondi_profile


PHASE_FACTOR = 8.0 * math.pi**2


@dataclass(frozen=True)
class KineticGrid:
    r_edges: np.ndarray
    vr_edges: np.ndarray
    j_edges: np.ndarray

    @property
    def r(self) -> np.ndarray:
        return np.sqrt(self.r_edges[:-1] * self.r_edges[1:])

    @property
    def vr(self) -> np.ndarray:
        return 0.5 * (self.vr_edges[:-1] + self.vr_edges[1:])

    @property
    def angular_momentum(self) -> np.ndarray:
        return 0.5 * (self.j_edges[:-1] + self.j_edges[1:])

    @property
    def dr(self) -> np.ndarray:
        return np.diff(self.r_edges)

    @property
    def dvr(self) -> np.ndarray:
        return np.diff(self.vr_edges)

    @property
    def dj(self) -> np.ndarray:
        return np.diff(self.j_edges)


def make_grid(r_min: float, r_max: float, radial_bins: int, radial_velocity_bins: int,
              tangential_bins: int, velocity_max: float) -> KineticGrid:
    if not (0.0 < r_min < r_max and radial_bins >= 8):
        raise ValueError("invalid radial grid")
    if radial_velocity_bins < 12 or tangential_bins < 8 or velocity_max <= 0.0:
        raise ValueError("invalid velocity grid")
    r_edges = np.geomspace(r_min, r_max, radial_bins + 1)
    vr_edges = np.linspace(-velocity_max, velocity_max, radial_velocity_bins + 1)
    # Quartic spacing resolves the loss-cone end of J while retaining the
    # thermal support J~r_max*v_t at the reservoir.  A uniform or merely
    # quadratic grid leaves the innermost capture cone in one oversized cell.
    unit = np.linspace(0.0, 1.0, tangential_bins + 1)
    j_edges = r_max * velocity_max * unit**4
    return KineticGrid(r_edges, vr_edges, j_edges)


def _moments_without_geometry(h: np.ndarray, grid: KineticGrid) -> tuple[np.ndarray, ...]:
    """Return integral moments of ``h=J f`` at each radius."""

    vr = grid.vr[None, :, None]
    vt = grid.angular_momentum[None, None, :] / grid.r[:, None, None]
    weight = grid.dvr[None, :, None] * grid.dj[None, None, :]
    m0 = np.sum(h * weight, axis=(1, 2))
    m1 = np.sum(h * vr * weight, axis=(1, 2))
    m2 = np.sum(h * (vr * vr + vt * vt) * weight, axis=(1, 2))
    return m0, m1, m2


def macroscopic_moments(h: np.ndarray, grid: KineticGrid) -> dict:
    m0, m1, m2 = _moments_without_geometry(h, grid)
    if np.any(m0 <= 0.0):
        raise FloatingPointError("non-positive kinetic density")
    radial_velocity = m1 / m0
    mean_v2 = m2 / m0
    temperature = (mean_v2 - radial_velocity**2) / 3.0
    if np.any(temperature <= 0.0):
        raise FloatingPointError("non-positive kinetic temperature")
    density = 2.0 * math.pi * m0 / grid.r**2
    return {
        "density": density,
        "radial_velocity": radial_velocity,
        "temperature": temperature,
        "mean_v2": mean_v2,
        "m0": m0,
        "m1": m1,
        "m2": m2,
    }


def discrete_maxwellian_from_moments(
    grid: KineticGrid,
    target_m0: np.ndarray,
    target_m1: np.ndarray,
    target_m2: np.ndarray,
    *,
    maximum_iterations: int = 24,
    tolerance: float = 2.0e-11,
) -> tuple[np.ndarray, dict]:
    """Return the maximum-entropy ``h=J f`` matching three discrete moments."""

    m0 = np.asarray(target_m0, dtype=float)
    m1 = np.asarray(target_m1, dtype=float)
    m2 = np.asarray(target_m2, dtype=float)
    if m0.shape != grid.r.shape or m1.shape != m0.shape or m2.shape != m0.shape:
        raise ValueError("target moments do not match the radial grid")
    if np.any(m0 <= 0.0) or np.any(~np.isfinite(np.r_[m0, m1, m2])):
        raise ValueError("invalid target moments")
    mean1 = m1 / m0
    mean2 = m2 / m0
    temperature = (mean2 - mean1**2) / 3.0
    if np.any(temperature <= 0.0):
        raise ValueError("target kinetic temperature is not positive")

    vr = grid.vr[None, :, None]
    vt = grid.angular_momentum[None, None, :] / grid.r[:, None, None]
    v2 = vr * vr + vt * vt
    j = grid.angular_momentum[None, None, :]
    weight = grid.dvr[None, :, None] * grid.dj[None, None, :]
    b = mean1 / temperature
    c = -0.5 / temperature
    converged = np.zeros(grid.r.size, dtype=bool)
    residual = np.full(grid.r.size, np.inf)

    for iteration in range(maximum_iterations):
        exponent = b[:, None, None] * vr + c[:, None, None] * v2
        offset = np.max(exponent, axis=(1, 2))
        raw = j * np.exp(np.clip(exponent - offset[:, None, None], -745.0, 0.0))
        z = np.sum(raw * weight, axis=(1, 2))
        p = raw * weight / z[:, None, None]
        e_vr = np.sum(p * vr, axis=(1, 2))
        e_v2 = np.sum(p * v2, axis=(1, 2))
        d1 = mean1 - e_vr
        d2 = mean2 - e_v2
        scale1 = np.sqrt(np.maximum(temperature, 1.0e-30))
        scale2 = np.maximum(mean2, 1.0e-30)
        residual = np.maximum(np.abs(d1) / scale1, np.abs(d2) / scale2)
        converged = residual < tolerance
        if np.all(converged):
            break

        var_vr = np.sum(p * (vr - e_vr[:, None, None]) ** 2, axis=(1, 2))
        var_v2 = np.sum(p * (v2 - e_v2[:, None, None]) ** 2, axis=(1, 2))
        cov = np.sum(
            p
            * (vr - e_vr[:, None, None])
            * (v2 - e_v2[:, None, None]),
            axis=(1, 2),
        )
        matrix = np.empty((grid.r.size, 2, 2))
        matrix[:, 0, 0] = var_vr
        matrix[:, 0, 1] = cov
        matrix[:, 1, 0] = cov
        matrix[:, 1, 1] = var_v2
        rhs = np.column_stack([d1, d2])
        determinant = var_vr * var_v2 - cov * cov
        if np.any(determinant <= 1.0e-24):
            raise FloatingPointError("discrete Maxwellian moment matrix is singular")
        delta = np.linalg.solve(matrix, rhs[..., None])[..., 0]
        delta[converged] = 0.0
        # On the finite velocity grid the maximum-entropy solution can have a
        # positive quadratic multiplier when a cell is dominated by the edge
        # of the represented support.  We allow that discrete solution so the
        # run finishes, then reject the velocity domain through its separate
        # boundary-loss gate instead of crashing in the collision step.
        step_norm = np.maximum(np.abs(delta[:, 0]), np.abs(delta[:, 1]))
        damping = np.minimum(1.0, 20.0 / np.maximum(step_norm, 1.0e-30))
        b += damping * delta[:, 0]
        c += damping * delta[:, 1]
    else:
        raise RuntimeError(
            "discrete Maxwellian did not converge; maximum scaled residual "
            f"{float(np.max(residual)):g}"
        )

    exponent = b[:, None, None] * vr + c[:, None, None] * v2
    offset = np.max(exponent, axis=(1, 2))
    raw = j * np.exp(np.clip(exponent - offset[:, None, None], -745.0, 0.0))
    z = np.sum(raw * weight, axis=(1, 2))
    equilibrium = raw * (m0 / z)[:, None, None]
    check0, check1, check2 = _moments_without_geometry(equilibrium, grid)
    errors = np.column_stack([
        np.abs(check0 / m0 - 1.0),
        np.abs(check1 - m1) / np.maximum(np.abs(m1), m0 * np.sqrt(temperature) * 1.0e-8),
        np.abs(check2 / m2 - 1.0),
    ])
    return equilibrium, {
        "iterations": iteration + 1,
        "maximum_scaled_newton_residual": float(np.max(residual)),
        "maximum_discrete_moment_relative_error": float(np.max(errors)),
    }


def bondi_maxwellian(grid: KineticGrid, gamma: float) -> tuple[np.ndarray, dict]:
    profile = critical_bondi_profile(grid.r, gamma)
    m0 = grid.r**2 * profile["density"] / (2.0 * math.pi)
    m1 = m0 * profile["radial_velocity"]
    m2 = m0 * (
        profile["radial_velocity"] ** 2 + 3.0 * profile["temperature"]
    )
    h, diagnostic = discrete_maxwellian_from_moments(grid, m0, m1, m2)
    return h, {"hydrodynamic_profile": profile, "maxwellian": diagnostic}


def _cumulative_at(q: np.ndarray, edges: np.ndarray, query: np.ndarray,
                   left_value: np.ndarray, right_value: np.ndarray) -> np.ndarray:
    """Piecewise-constant antiderivative for many columns."""

    ncell, ncolumn = q.shape
    width = np.diff(edges)
    cumulative = np.vstack([
        np.zeros((1, ncolumn)), np.cumsum(q * width[:, None], axis=0)
    ])
    index = np.searchsorted(edges, query, side="right") - 1
    columns = np.broadcast_to(np.arange(ncolumn)[None, :], query.shape)
    clipped = np.clip(index, 0, ncell - 1)
    values = (
        cumulative[clipped, columns]
        + q[clipped, columns] * (query - edges[clipped])
    )
    low = query <= edges[0]
    high = query >= edges[-1]
    values[low] = ((query - edges[0]) * left_value[None, :])[low]
    values[high] = (
        cumulative[-1][None, :] + (query - edges[-1]) * right_value[None, :]
    )[high]
    return values


def conservative_shift(q: np.ndarray, edges: np.ndarray, shift: np.ndarray,
                       left_value: np.ndarray, right_value: np.ndarray) -> tuple[np.ndarray, dict]:
    """Translate columns by arbitrary distances with a conservative remap."""

    q = np.asarray(q, dtype=float)
    trailing = q.shape[1:]
    flat = q.reshape(q.shape[0], -1)
    shifts = np.broadcast_to(np.asarray(shift, dtype=float), trailing).reshape(-1)
    left = np.broadcast_to(np.asarray(left_value, dtype=float), trailing).reshape(-1)
    right = np.broadcast_to(np.asarray(right_value, dtype=float), trailing).reshape(-1)
    query = edges[:, None] - shifts[None, :]
    cumulative = _cumulative_at(flat, edges, query, left, right)
    shifted = np.diff(cumulative, axis=0) / np.diff(edges)[:, None]
    shifted = shifted.reshape(q.shape)
    if np.any(shifted < -1.0e-13 * max(float(np.max(q)), 1.0)):
        raise FloatingPointError("conservative remap produced a negative distribution")
    shifted = np.maximum(shifted, 0.0)

    # Exact phase-density crossing ledgers for shifts smaller than the domain.
    domain = edges[-1] - edges[0]
    old_edges = np.broadcast_to(edges[:, None], query.shape)
    c_old = _cumulative_at(
        flat, edges, old_edges, np.zeros_like(left), np.zeros_like(right)
    )
    lower_sample = np.clip(edges[0] - shifts, edges[0], edges[-1])[None, :]
    upper_sample = np.clip(edges[-1] - shifts, edges[0], edges[-1])[None, :]
    c_lower = _cumulative_at(flat, edges, lower_sample, np.zeros_like(left), np.zeros_like(right))[0]
    c_upper = _cumulative_at(flat, edges, upper_sample, np.zeros_like(left), np.zeros_like(right))[0]
    left_out = np.where(shifts < 0.0, c_lower, 0.0)
    right_out = np.where(shifts > 0.0, c_old[-1] - c_upper, 0.0)
    left_in = np.where(shifts > 0.0, np.minimum(shifts, domain) * left, 0.0)
    right_in = np.where(shifts < 0.0, np.minimum(-shifts, domain) * right, 0.0)
    return shifted, {
        "left_out": left_out.reshape(trailing),
        "right_out": right_out.reshape(trailing),
        "left_in": left_in.reshape(trailing),
        "right_in": right_in.reshape(trailing),
    }


def total_mass(h: np.ndarray, grid: KineticGrid) -> float:
    return float(PHASE_FACTOR * np.sum(
        h * grid.dr[:, None, None] * grid.dvr[None, :, None] * grid.dj[None, None, :]
    ))


def outer_reservoir(grid: KineticGrid, gamma: float) -> tuple[np.ndarray, dict]:
    # The finite outer surface is held at the critical Bondi solution evaluated
    # there, rather than at an inconsistent finite-radius copy of infinity.
    profile = critical_bondi_profile(np.array([grid.r[-1], grid.r_edges[-1]]), gamma)
    rho = float(profile["density"][-1])
    u = float(profile["radial_velocity"][-1])
    temperature = float(profile["temperature"][-1])
    r = grid.r_edges[-1]
    vr = grid.vr[:, None]
    vt = grid.angular_momentum[None, :] / r
    f = rho / (2.0 * math.pi * temperature) ** 1.5 * np.exp(
        -((vr - u) ** 2 + vt**2) / (2.0 * temperature)
    )
    h = grid.angular_momentum[None, :] * f
    return h, {"density": rho, "radial_velocity": u, "temperature": temperature}


def ballistic_capture_rate(grid: KineticGrid, reservoir: np.ndarray) -> float:
    """Discrete collisionless capture current from the maintained reservoir."""

    vr = grid.vr[:, None]
    j = grid.angular_momentum[None, :]
    r_outer = grid.r_edges[-1]
    r_inner = grid.r_edges[0]
    vt_outer = j / r_outer
    radial_speed_squared_at_inner = (
        vr * vr
        + vt_outer * vt_outer
        + 2.0 * (1.0 / r_inner - 1.0 / r_outer)
        - (j / r_inner) ** 2
    )
    capture = (vr < 0.0) & (radial_speed_squared_at_inner >= 0.0)
    weight = grid.dvr[:, None] * grid.dj[None, :]
    vr_full = np.broadcast_to(vr, capture.shape)
    weight_full = np.broadcast_to(weight, capture.shape)
    return float(-PHASE_FACTOR * np.sum(
        vr_full[capture] * reservoir[capture] * weight_full[capture]
    ))


def collision_step(h: np.ndarray, grid: KineticGrid, dt: float, knudsen_infinity: float,
                   temperature_infinity: float) -> tuple[np.ndarray, dict]:
    moments = macroscopic_moments(h, grid)
    equilibrium, diagnostic = discrete_maxwellian_from_moments(
        grid, moments["m0"], moments["m1"], moments["m2"]
    )
    tau = knudsen_infinity / (
        moments["density"]
        * np.sqrt(np.maximum(moments["temperature"] / temperature_infinity, 1.0e-30))
    )
    relaxation = -np.expm1(-dt / np.maximum(tau, 1.0e-14))
    updated = h + relaxation[:, None, None] * (equilibrium - h)
    before = np.column_stack(_moments_without_geometry(h, grid))
    after = np.column_stack(_moments_without_geometry(updated, grid))
    scale = np.maximum(np.abs(before), np.column_stack([
        before[:, 0], before[:, 0], before[:, 0]
    ]) * 1.0e-12)
    diagnostic.update({
        "minimum_tau": float(np.min(tau)),
        "maximum_tau": float(np.max(tau)),
        "maximum_collision_invariant_relative_error": float(
            np.max(np.abs(after - before) / scale)
        ),
    })
    return updated, diagnostic


def drift_step(h: np.ndarray, grid: KineticGrid, dt: float, reservoir: np.ndarray) -> tuple[np.ndarray, dict]:
    shifted, ledger = conservative_shift(
        h,
        grid.r_edges,
        grid.vr[:, None] * dt,
        np.zeros((grid.vr.size, grid.angular_momentum.size)),
        reservoir,
    )
    weight = grid.dvr[:, None] * grid.dj[None, :]
    vr = grid.vr[:, None]
    j = grid.angular_momentum[None, :]
    capture_weight = ledger["left_out"] * weight
    specific_energy = 0.5 * (vr * vr + (j / grid.r_edges[0]) ** 2) - 1.0 / grid.r_edges[0]
    return shifted, {
        "captured_mass": float(PHASE_FACTOR * np.sum(capture_weight)),
        "captured_mechanical_energy": float(
            PHASE_FACTOR * np.sum(capture_weight * specific_energy)
        ),
        "outer_injected_mass": float(PHASE_FACTOR * np.sum(ledger["right_in"] * weight)),
        "outer_exported_mass": float(PHASE_FACTOR * np.sum(ledger["right_out"] * weight)),
    }


def acceleration_step(h: np.ndarray, grid: KineticGrid, dt: float) -> tuple[np.ndarray, dict]:
    r = grid.r[:, None]
    j = grid.angular_momentum[None, :]
    acceleration = j * j / r**3 - 1.0 / r**2
    q = np.moveaxis(h, 1, 0)
    shifted, ledger = conservative_shift(
        q,
        grid.vr_edges,
        acceleration * dt,
        np.zeros_like(acceleration),
        np.zeros_like(acceleration),
    )
    radial_weight = grid.dr[:, None] * grid.dj[None, :]
    lost = ledger["left_out"] + ledger["right_out"]
    return np.moveaxis(shifted, 0, 1), {
        "velocity_boundary_mass_loss": float(PHASE_FACTOR * np.sum(lost * radial_weight)),
    }


def transport_step(h: np.ndarray, grid: KineticGrid, dt: float,
                   reservoir: np.ndarray) -> tuple[np.ndarray, dict]:
    h, first = drift_step(h, grid, 0.5 * dt, reservoir)
    h, kick = acceleration_step(h, grid, dt)
    h, second = drift_step(h, grid, 0.5 * dt, reservoir)
    return h, {
        "captured_mass": first["captured_mass"] + second["captured_mass"],
        "captured_mechanical_energy": (
            first["captured_mechanical_energy"] + second["captured_mechanical_energy"]
        ),
        "outer_injected_mass": first["outer_injected_mass"] + second["outer_injected_mass"],
        "outer_exported_mass": first["outer_exported_mass"] + second["outer_exported_mass"],
        "velocity_boundary_mass_loss": kick["velocity_boundary_mass_loss"],
    }


def _stationarity(time_values: np.ndarray, current_values: np.ndarray) -> dict:
    keep = time_values >= time_values[0] + 0.7 * (time_values[-1] - time_values[0])
    t = time_values[keep]
    current = current_values[keep]
    mean = float(np.mean(current))
    if t.size < 4 or mean <= 0.0:
        return {"status": "FAIL", "reason": "too few positive late samples"}
    slope = float(np.polyfit(t, current, 1)[0])
    fractional_change = abs(slope) * (t[-1] - t[0]) / mean
    coefficient_of_variation = float(np.std(current) / mean)
    return {
        "late_mean_mdot": mean,
        "late_fractional_linear_change": fractional_change,
        "late_coefficient_of_variation": coefficient_of_variation,
        "samples": int(t.size),
        "status": (
            "PASS" if fractional_change < 0.03 and coefficient_of_variation < 0.03 else "FAIL"
        ),
    }


def run_transition(
    *,
    knudsen_infinity: float,
    gamma: float = 5.0 / 3.0,
    r_min: float = 0.1,
    r_max: float = 8.0,
    radial_bins: int = 48,
    radial_velocity_bins: int = 56,
    tangential_bins: int = 32,
    velocity_max: float = 8.0,
    dt: float = 0.005,
    t_end: float = 24.0,
    output_dt: float = 0.2,
    wall_seconds: float = math.inf,
) -> dict:
    if knudsen_infinity <= 0.0 or dt <= 0.0 or t_end <= 0.0:
        raise ValueError("Kn, dt, and t_end must be positive")
    grid = make_grid(
        r_min, r_max, radial_bins, radial_velocity_bins, tangential_bins, velocity_max
    )
    h, initial = bondi_maxwellian(grid, gamma)
    reservoir, reservoir_state = outer_reservoir(grid, gamma)
    ballistic_mdot = ballistic_capture_rate(grid, reservoir)
    initial_mass = total_mass(h, grid)
    cumulative = {
        "captured_mass": 0.0,
        "captured_mechanical_energy": 0.0,
        "outer_injected_mass": 0.0,
        "outer_exported_mass": 0.0,
        "velocity_boundary_mass_loss": 0.0,
    }
    maximum_collision_error = 0.0
    maximum_newton_residual = 0.0
    rows = []
    current_time = 0.0
    next_output = 0.0
    started = time.monotonic()
    status = "COMPLETE"
    steps = 0
    while current_time < t_end * (1.0 - 1.0e-13):
        if time.monotonic() - started > wall_seconds:
            status = "WALLTIME_STOP"
            break
        step_dt = min(dt, t_end - current_time)
        h, transport = transport_step(h, grid, step_dt, reservoir)
        for key, value in transport.items():
            cumulative[key] += value
        h, collision = collision_step(
            h, grid, step_dt, knudsen_infinity, 1.0 / gamma
        )
        maximum_collision_error = max(
            maximum_collision_error,
            collision["maximum_collision_invariant_relative_error"],
        )
        maximum_newton_residual = max(
            maximum_newton_residual,
            collision["maximum_scaled_newton_residual"],
        )
        current_time += step_dt
        steps += 1
        if current_time + 1.0e-12 >= next_output or current_time >= t_end:
            moments = macroscopic_moments(h, grid)
            first_h = h[0]
            inward = grid.vr < 0.0
            weight = grid.dvr[:, None] * grid.dj[None, :]
            mdot_face = float(-PHASE_FACTOR * np.sum(
                grid.vr[inward, None] * first_h[inward] * weight[inward]
            ))
            energy = 0.5 * (
                grid.vr[:, None] ** 2
                + (grid.angular_momentum[None, :] / grid.r_edges[0]) ** 2
            ) - 1.0 / grid.r_edges[0]
            energy_face = float(-PHASE_FACTOR * np.sum(
                grid.vr[inward, None]
                * first_h[inward]
                * energy[inward]
                * weight[inward]
            ))
            rows.append({
                "time": current_time,
                "mdot_inner_face": mdot_face,
                "mechanical_energy_current_inner_face": energy_face,
                "mass_on_grid": total_mass(h, grid),
                "rho_inner": float(moments["density"][0]),
                "u_inner": float(moments["radial_velocity"][0]),
                "temperature_inner": float(moments["temperature"][0]),
                "rho_outer": float(moments["density"][-1]),
                "u_outer": float(moments["radial_velocity"][-1]),
            })
            next_output += output_dt

    final_mass = total_mass(h, grid)
    expected_final_mass = (
        initial_mass
        + cumulative["outer_injected_mass"]
        - cumulative["outer_exported_mass"]
        - cumulative["captured_mass"]
        - cumulative["velocity_boundary_mass_loss"]
    )
    mass_residual = final_mass - expected_final_mass
    times = np.asarray([row["time"] for row in rows])
    mdot = np.asarray([row["mdot_inner_face"] for row in rows])
    stationarity = _stationarity(times, mdot)
    hydro_mdot = float(initial["hydrodynamic_profile"]["mdot_exact"])
    gates = {
        "completed_requested_time": status == "COMPLETE",
        "positive_distribution": bool(np.all(np.isfinite(h)) and np.all(h >= 0.0)),
        "collision_invariants": maximum_collision_error < 2.0e-9,
        "phase_mass_ledger": abs(mass_residual) / max(initial_mass, 1.0e-300) < 2.0e-8,
        "velocity_domain": (
            cumulative["velocity_boundary_mass_loss"] / max(initial_mass, 1.0e-300) < 1.0e-4
        ),
        "late_current_stationary": stationarity.get("status") == "PASS",
    }
    result = {
        "schema": "non-orbit-averaged-radial-bgk-v1",
        "status": "KINETIC_RUN_PASS" if all(gates.values()) else "KINETIC_RUN_DIAGNOSTIC",
        "physical_scope": (
            "Spherical BGK transition model with a maintained Bondi-family reservoir. "
            "It resolves odd radial-velocity moments but does not yet use the Yukawa "
            "finite-angle collision kernel."
        ),
        "coordinates": {
            "stored": "f(r,v_r,J) represented through h=J f",
            "equivalent_requested_form": "f(r,v_r,v_t), with v_t=J/r",
            "orbit_averaged": False,
        },
        "normalization": initial["hydrodynamic_profile"]["normalization"],
        "gamma": gamma,
        "knudsen_infinity": knudsen_infinity,
        "grid": {
            "r_min": r_min,
            "r_max": r_max,
            "radial_bins": radial_bins,
            "radial_velocity_bins": radial_velocity_bins,
            "tangential_bins": tangential_bins,
            "velocity_max": velocity_max,
        },
        "integration": {
            "dt": dt,
            "t_end_requested": t_end,
            "t_end_reached": current_time,
            "steps": steps,
            "wall_seconds": time.monotonic() - started,
        },
        "outer_reservoir": reservoir_state,
        "hydrodynamic_reference": {
            "mdot": hydro_mdot,
            "bondi_lambda": float(initial["hydrodynamic_profile"]["eigenvalue"]),
            "maximum_mass_flux_relative_residual": float(
                initial["hydrodynamic_profile"]["maximum_mass_flux_relative_residual"]
            ),
            "maximum_bernoulli_absolute_residual": float(
                initial["hydrodynamic_profile"]["maximum_bernoulli_absolute_residual"]
            ),
        },
        "ballistic_reference": {
            "mdot": ballistic_mdot,
            "definition": (
                "Outer-reservoir inward current whose Kepler periapse reaches "
                "the absorbing inner surface, integrated on the same velocity grid."
            ),
        },
        "stationarity": stationarity,
        "late_mdot_over_hydrodynamic": (
            stationarity.get("late_mean_mdot", math.nan) / hydro_mdot
        ),
        "cumulative_ledger": cumulative,
        "mass_ledger": {
            "initial_mass": initial_mass,
            "final_mass": final_mass,
            "expected_final_mass": expected_final_mass,
            "residual": mass_residual,
            "relative_residual": mass_residual / max(initial_mass, 1.0e-300),
        },
        "maximum_collision_invariant_relative_error": maximum_collision_error,
        "maximum_maxwellian_scaled_residual": maximum_newton_residual,
        "gates": gates,
        "trajectory": rows,
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kn", type=float, required=True)
    parser.add_argument("--gamma", type=float, default=5.0 / 3.0)
    parser.add_argument("--r-min", type=float, default=0.1)
    parser.add_argument("--r-max", type=float, default=8.0)
    parser.add_argument("--radial-bins", type=int, default=48)
    parser.add_argument("--radial-velocity-bins", type=int, default=56)
    parser.add_argument("--tangential-bins", type=int, default=32)
    parser.add_argument("--velocity-max", type=float, default=8.0)
    parser.add_argument("--dt", type=float, default=0.005)
    parser.add_argument("--t-end", type=float, default=24.0)
    parser.add_argument("--output-dt", type=float, default=0.2)
    parser.add_argument("--wall-seconds", type=float, default=math.inf)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = run_transition(
        knudsen_infinity=args.kn,
        gamma=args.gamma,
        r_min=args.r_min,
        r_max=args.r_max,
        radial_bins=args.radial_bins,
        radial_velocity_bins=args.radial_velocity_bins,
        tangential_bins=args.tangential_bins,
        velocity_max=args.velocity_max,
        dt=args.dt,
        t_end=args.t_end,
        output_dt=args.output_dt,
        wall_seconds=args.wall_seconds,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": result["status"],
        "knudsen_infinity": args.kn,
        "late_mdot": result["stationarity"].get("late_mean_mdot"),
        "late_mdot_over_hydrodynamic": result["late_mdot_over_hydrodynamic"],
        "gates": result["gates"],
    }, indent=2, sort_keys=True))
    return 0 if result["status"] == "KINETIC_RUN_PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
