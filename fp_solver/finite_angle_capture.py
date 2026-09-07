#!/usr/bin/env python3
"""Estimate finite-angle SIDM transitions into the black-hole loss cone.

The GNC diffusion sink is invalid when fewer than many physical collisions occur
during one loss-cone diffusion time.  This diagnostic evaluates the corresponding
single-collision Boltzmann source directly.  It samples two particles from the
tabulated isotropic DF at each radius, excludes particles already inside the loss
cone, draws the exact t-channel scattering angle, and counts bound post-collision
states with J < J_lc.

For every post-collision loss-cone state, the code then integrates the exact
single-particle total collision rate along its Kepler orbit. The no-collision
probability to the relativistic plunge surface isolates the direct,
no-rescattering part of the current. The total injection rate is an isotropic
ceiling. Neither quantity replaces the steady capture current until the
depleted, non-local distribution has been solved.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from functools import lru_cache
from pathlib import Path
from typing import Optional

import numpy as np

from normalize_gnc_df import df_forward_matrix, load_profile, log_interp
import ej_energy_ledger
import ej_energy_dg

G_PC_KMS2_MSUN = 4.30091e-3
C_KMS = 2.99792458e5
MSUN_G = 1.98847e33
PC_CM = 3.0856775814913673e18
CM2_G_TO_PC2_MSUN = MSUN_G / PC_CM**2
KMS_TO_PC_PER_MYR = 1.022712165045695
SOURCE_J_OVER_JLC_EDGES = np.array([1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0, np.inf])
SOURCE_J_FLOOR_THRESHOLDS = (3.0, 5.0, 10.0, 30.0)
JUMP_J_OVER_JLC_EDGES = np.array([
    1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0,
    10.0, 15.0, 20.0, 30.0, 50.0, 100.0, 300.0, 1000.0, np.inf,
])
EJ_J_OVER_JLC_EDGES_COARSE = np.array([
    1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0, 300.0, 1000.0, np.inf,
])
_ej_finite_default = JUMP_J_OVER_JLC_EDGES[np.isfinite(JUMP_J_OVER_JLC_EDGES)]
EJ_J_OVER_JLC_EDGES_FINE = np.concatenate([
    np.column_stack([
        _ej_finite_default[:-1],
        np.sqrt(_ej_finite_default[:-1] * _ej_finite_default[1:]),
    ]).reshape(-1),
    _ej_finite_default[-1:],
    [np.inf],
])
_ej_finite_fine = EJ_J_OVER_JLC_EDGES_FINE[
    np.isfinite(EJ_J_OVER_JLC_EDGES_FINE)
]
EJ_J_OVER_JLC_EDGES_ULTRA = np.concatenate([
    np.column_stack([
        _ej_finite_fine[:-1],
        np.sqrt(_ej_finite_fine[:-1] * _ej_finite_fine[1:]),
    ]).reshape(-1),
    _ej_finite_fine[-1:],
    [np.inf],
])
_ej_finite_ultra = EJ_J_OVER_JLC_EDGES_ULTRA[
    np.isfinite(EJ_J_OVER_JLC_EDGES_ULTRA)
]
EJ_J_OVER_JLC_EDGES_HYPER = np.concatenate([
    np.column_stack([
        _ej_finite_ultra[:-1],
        np.sqrt(_ej_finite_ultra[:-1] * _ej_finite_ultra[1:]),
    ]).reshape(-1),
    _ej_finite_ultra[-1:],
    [np.inf],
])
_ej_finite_hyper = EJ_J_OVER_JLC_EDGES_HYPER[
    np.isfinite(EJ_J_OVER_JLC_EDGES_HYPER)
]
EJ_J_OVER_JLC_EDGES_J257 = np.concatenate([
    np.column_stack([
        _ej_finite_hyper[:-1],
        np.sqrt(_ej_finite_hyper[:-1] * _ej_finite_hyper[1:]),
    ]).reshape(-1),
    _ej_finite_hyper[-1:],
    [np.inf],
])
EJ_ANGULAR_GRIDS = {
    "coarse": EJ_J_OVER_JLC_EDGES_COARSE,
    "default": JUMP_J_OVER_JLC_EDGES,
    "fine": EJ_J_OVER_JLC_EDGES_FINE,
    "ultra": EJ_J_OVER_JLC_EDGES_ULTRA,
    "hyper": EJ_J_OVER_JLC_EDGES_HYPER,
    "j257": EJ_J_OVER_JLC_EDGES_J257,
}
EJ_POST_RESERVOIR = -2
EJ_POST_CAPTURE = -1


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@lru_cache(maxsize=None)
def _leggauss(order: int) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return nodes, weights


def read_df_table(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    lines = path.read_text().splitlines()
    if len(lines) < 5 or not lines[0].startswith("#"):
        raise ValueError("malformed GNC DF table")
    metadata = {}
    marker = "metadata="
    if marker in lines[0]:
        metadata = json.loads(lines[0].split(marker, 1)[1])
    nx, nj = map(int, lines[1].split())
    xgrid = np.fromstring(lines[2], sep=" ")
    jgrid = np.fromstring(lines[3], sep=" ")
    table = np.loadtxt(lines[4:])
    if table.shape != (nx, nj) or xgrid.size != nx or jgrid.size != nj:
        raise ValueError("DF table dimensions do not match its header")
    if np.any(table <= 0.0) or np.any(~np.isfinite(table)):
        raise ValueError("DF table must be positive and finite")
    return xgrid, jgrid, table, metadata


def sample_tchannel_u(a: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Sample u=sin^2(theta/2) from the exact conditional kernel."""

    a = np.asarray(a, float)
    if np.any(a < 0.0) or np.any(~np.isfinite(a)):
        raise ValueError("t-channel angular parameter must be finite and non-negative")
    uniform = rng.random(a.shape)
    return uniform / (1.0 + a * (1.0 - uniform))


def scatter_equal_mass(
    v1: np.ndarray, v2: np.ndarray, w_kms: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Draw exact equal-mass post-collision velocities and sigma_tot/sigma0."""

    relative = v1 - v2
    speed = np.linalg.norm(relative, axis=1)
    if np.any(speed <= 0.0):
        # Exact duplicates have zero collision weight.  Give them a harmless
        # direction so the vector algebra remains finite.
        speed = np.maximum(speed, np.finfo(float).tiny)
    direction = relative / speed[:, None]
    a = (speed / w_kms) ** 2
    u = sample_tchannel_u(a, rng)
    cos_theta = 1.0 - 2.0 * u
    sin_theta = 2.0 * np.sqrt(np.maximum(u * (1.0 - u), 0.0))
    azimuth = 2.0 * math.pi * rng.random(speed.size)

    reference = np.zeros_like(direction)
    use_x = np.abs(direction[:, 2]) > 0.9
    reference[~use_x, 2] = 1.0
    reference[use_x, 0] = 1.0
    e1 = np.cross(direction, reference)
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(direction, e1)
    rotated = (
        cos_theta[:, None] * direction
        + sin_theta[:, None]
        * (np.cos(azimuth)[:, None] * e1 + np.sin(azimuth)[:, None] * e2)
    )
    relative_after = speed[:, None] * rotated
    centre = 0.5 * (v1 + v2)
    return (
        centre + 0.5 * relative_after,
        centre - 0.5 * relative_after,
        1.0 / (1.0 + a),
    )


def circle_intersection_area(
    separation: np.ndarray,
    radius: np.ndarray,
    unit_radius: float = 1.0,
) -> np.ndarray:
    """Area shared by a unit disc and a displaced disc of given radius."""

    d, r = np.broadcast_arrays(
        np.asarray(separation, float), np.asarray(radius, float)
    )
    if np.any(d < 0.0) or np.any(r < 0.0) or unit_radius <= 0.0:
        raise ValueError("circle separations and radii must be non-negative")
    out = np.zeros_like(d)
    disjoint = d >= unit_radius + r
    contained = d <= np.abs(unit_radius - r)
    out[contained] = math.pi * np.minimum(unit_radius, r[contained]) ** 2
    partial = ~(disjoint | contained)
    if np.any(partial):
        dp = d[partial]
        rp = r[partial]
        cosine_unit = np.clip(
            (dp**2 + unit_radius**2 - rp**2) / (2.0 * dp * unit_radius),
            -1.0,
            1.0,
        )
        cosine_r = np.clip(
            (dp**2 + rp**2 - unit_radius**2) / (2.0 * dp * rp),
            -1.0,
            1.0,
        )
        radical = np.maximum(
            (-dp + unit_radius + rp)
            * (dp + unit_radius - rp)
            * (dp - unit_radius + rp)
            * (dp + unit_radius + rp),
            0.0,
        )
        out[partial] = (
            unit_radius**2 * np.arccos(cosine_unit)
            + rp**2 * np.arccos(cosine_r)
            - 0.5 * np.sqrt(radical)
        )
    return out


def scatter_equal_mass_disc_importance(
    v1: np.ndarray,
    v2: np.ndarray,
    w_kms: float,
    radius_pc: float,
    j_lc_pc_kms: float,
    rng: np.random.Generator,
    cap_fraction: float = 0.5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Scatter with an exactly weighted loss-cone-cap proposal.

    The proposal mixes the physical t-channel law with directions whose
    projected post-collision velocity lies inside either particle's loss-cone
    disc.  The returned importance weight restores the exact differential
    cross section, while replacing rare zero/one candidate counts by a smooth
    weighted estimate.
    """

    if not (0.0 <= cap_fraction < 1.0):
        raise ValueError("cap_fraction must lie in [0,1)")
    relative = v1 - v2
    speed = np.linalg.norm(relative, axis=1)
    safe_speed = np.maximum(speed, np.finfo(float).tiny)
    direction = relative / safe_speed[:, None]
    centre = 0.5 * (v1 + v2)
    half_speed = 0.5 * safe_speed
    a = (safe_speed / w_kms) ** 2

    cap_radius = j_lc_pc_kms / (radius_pc * half_speed)
    cap_centre_1 = -centre[:, :2] / half_speed[:, None]
    cap_centre_2 = -cap_centre_1
    cap_separation = np.linalg.norm(cap_centre_1, axis=1)
    cap_area = circle_intersection_area(cap_separation, cap_radius)
    parent_area = math.pi * np.minimum(1.0, cap_radius**2)
    rejection_efficiency = np.divide(
        cap_area,
        parent_area,
        out=np.zeros_like(cap_area),
        where=parent_area > 0.0,
    )
    # Uniform projected-disc sampling is efficient for ordinary intersections.
    # A focused von Mises-Fisher component covers thin near-tangent lenses,
    # which would otherwise rely on rare baseline events with large variance.
    efficient_disc = rejection_efficiency >= 0.02
    alpha = np.where(cap_area > 0.0, cap_fraction, 0.0)
    use_cap = rng.random(speed.size) < alpha
    cap_for_particle_1 = rng.random(speed.size) < 0.5
    chosen_centre = np.where(
        cap_for_particle_1[:, None], cap_centre_1, cap_centre_2
    )

    post_direction = np.empty_like(direction)
    baseline = ~use_cap
    if np.any(baseline):
        ids = np.flatnonzero(baseline)
        u = sample_tchannel_u(a[ids], rng)
        cos_theta = 1.0 - 2.0 * u
        sin_theta = 2.0 * np.sqrt(np.maximum(u * (1.0 - u), 0.0))
        azimuth = 2.0 * math.pi * rng.random(ids.size)
        reference = np.zeros((ids.size, 3))
        use_x = np.abs(direction[ids, 2]) > 0.9
        reference[~use_x, 2] = 1.0
        reference[use_x, 0] = 1.0
        e1 = np.cross(direction[ids], reference)
        e1 /= np.linalg.norm(e1, axis=1)[:, None]
        e2 = np.cross(direction[ids], e1)
        post_direction[ids] = (
            cos_theta[:, None] * direction[ids]
            + sin_theta[:, None]
            * (np.cos(azimuth)[:, None] * e1 + np.sin(azimuth)[:, None] * e2)
        )

    # Draw uniformly from each chosen loss-cone disc intersected with the unit
    # projected sphere.  Rejection from the smaller parent disc is efficient
    # both for narrow caps and for caps covering most of the sphere.
    unresolved = np.flatnonzero(use_cap & efficient_disc)
    while unresolved.size:
        delta = cap_radius[unresolved]
        centres = chosen_centre[unresolved]
        sample_target_disc = delta <= 1.0
        angle = 2.0 * math.pi * rng.random(unresolved.size)
        radial = np.sqrt(rng.random(unresolved.size))
        xy = np.empty((unresolved.size, 2))
        xy[sample_target_disc] = (
            centres[sample_target_disc]
            + (delta[sample_target_disc] * radial[sample_target_disc])[:, None]
            * np.column_stack([
                np.cos(angle[sample_target_disc]),
                np.sin(angle[sample_target_disc]),
            ])
        )
        xy[~sample_target_disc] = radial[~sample_target_disc, None] * np.column_stack([
            np.cos(angle[~sample_target_disc]),
            np.sin(angle[~sample_target_disc]),
        ])
        inside_unit = np.sum(xy * xy, axis=1) <= 1.0
        inside_target = (
            np.sum((xy - centres) ** 2, axis=1) <= delta**2
        )
        accepted = inside_unit & inside_target
        accepted_ids = unresolved[accepted]
        accepted_xy = xy[accepted]
        z = np.sqrt(np.maximum(1.0 - np.sum(accepted_xy**2, axis=1), 0.0))
        z *= np.where(rng.random(accepted_ids.size) < 0.5, -1.0, 1.0)
        post_direction[accepted_ids] = np.column_stack([accepted_xy, z])
        unresolved = unresolved[~accepted]

    thin = use_cap & ~efficient_disc
    if np.any(thin):
        focus1 = _nearest_loss_cone_direction(
            direction[thin], cap_centre_1[thin], cap_radius[thin]
        )
        focus2 = _nearest_loss_cone_direction(
            direction[thin], cap_centre_2[thin], cap_radius[thin]
        )
        concentration = np.clip(
            np.maximum(
                0.5 * a[thin],
                2.0 / np.maximum(cap_radius[thin], 1.0e-8) ** 2,
            ),
            1.0,
            1.0e8,
        )
        chosen_focus = np.where(
            cap_for_particle_1[thin, None], focus1, focus2
        )
        post_direction[thin] = _sample_von_mises_fisher(
            chosen_focus, concentration, rng
        )

    dot = np.sum(direction * post_direction, axis=1)
    u = np.clip(0.5 * (1.0 - dot), 0.0, 1.0)
    qphysical = (1.0 + a) / (4.0 * math.pi * (1.0 + a * u) ** 2)
    xy = post_direction[:, :2]
    in_cap_1 = np.sum((xy - cap_centre_1) ** 2, axis=1) <= cap_radius**2
    in_cap_2 = np.sum((xy - cap_centre_2) ** 2, axis=1) <= cap_radius**2
    qcap_prefactor = np.divide(
        np.abs(post_direction[:, 2]),
        2.0 * cap_area,
        out=np.zeros_like(cap_area),
        where=cap_area > 0.0,
    )
    qdisc_mixture = 0.5 * qcap_prefactor * (
        in_cap_1.astype(float) + in_cap_2.astype(float)
    )
    focus1_all = _nearest_loss_cone_direction(
        direction, cap_centre_1, cap_radius
    )
    focus2_all = _nearest_loss_cone_direction(
        direction, cap_centre_2, cap_radius
    )
    concentration_all = np.clip(
        np.maximum(0.5 * a, 2.0 / np.maximum(cap_radius, 1.0e-8) ** 2),
        1.0,
        1.0e8,
    )
    qfocus_mixture = 0.5 * (
        _von_mises_fisher_density(
            post_direction, focus1_all, concentration_all
        )
        + _von_mises_fisher_density(
            post_direction, focus2_all, concentration_all
        )
    )
    qcap_mixture = np.where(efficient_disc, qdisc_mixture, qfocus_mixture)
    qproposal = (1.0 - alpha) * qphysical + alpha * qcap_mixture
    importance_weight = qphysical / qproposal
    relative_after = safe_speed[:, None] * post_direction
    return (
        centre + 0.5 * relative_after,
        centre - 0.5 * relative_after,
        1.0 / (1.0 + a),
        importance_weight,
    )


def _nearest_loss_cone_direction(
    physical_direction: np.ndarray,
    cap_centre: np.ndarray,
    cap_radius: np.ndarray,
) -> np.ndarray:
    """Direction in a projected loss-cone disc nearest the physical kernel."""

    kxy = physical_direction[:, :2]
    offset = kxy - cap_centre
    distance = np.linalg.norm(offset, axis=1)
    safe_distance = np.maximum(distance, np.finfo(float).tiny)
    xy = cap_centre + (
        cap_radius / safe_distance
    )[:, None] * offset
    degenerate = distance <= np.finfo(float).tiny
    xy[degenerate] = cap_centre[degenerate] + np.column_stack([
        cap_radius[degenerate], np.zeros(np.count_nonzero(degenerate))
    ])
    inside = distance <= cap_radius
    xy[inside] = kxy[inside]
    normxy = np.linalg.norm(xy, axis=1)
    outside_unit = normxy > 1.0
    xy[outside_unit] /= normxy[outside_unit, None]

    # Only a narrow circle-circle lens needs scalar endpoint handling.
    bad = np.linalg.norm(xy - cap_centre, axis=1) > cap_radius * (1.0 + 1.0e-10)
    for i in np.flatnonzero(bad):
        centre = cap_centre[i]
        delta = cap_radius[i]
        d = float(np.linalg.norm(centre))
        if d <= 0.0:
            xy[i] = kxy[i] / max(float(np.linalg.norm(kxy[i])), 1.0)
            continue
        along = (1.0 - delta**2 + d**2) / (2.0 * d)
        height = math.sqrt(max(1.0 - along**2, 0.0))
        axis = centre / d
        perpendicular = np.array([-axis[1], axis[0]])
        options = np.vstack([
            along * axis + height * perpendicular,
            along * axis - height * perpendicular,
        ])
        xy[i] = options[int(np.argmax(options @ kxy[i]))]
    z = np.sqrt(np.maximum(1.0 - np.sum(xy * xy, axis=1), 0.0))
    z *= np.where(physical_direction[:, 2] < 0.0, -1.0, 1.0)
    out = np.column_stack([xy, z])
    out /= np.linalg.norm(out, axis=1)[:, None]
    return out


def _sample_von_mises_fisher(
    mean_direction: np.ndarray,
    concentration: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Sample the three-dimensional von Mises-Fisher distribution."""

    n = mean_direction.shape[0]
    uniform = rng.random(n)
    exp_tail = np.exp(-2.0 * concentration)
    cosine = 1.0 + np.log(uniform + (1.0 - uniform) * exp_tail) / concentration
    cosine = np.clip(cosine, -1.0, 1.0)
    sine = np.sqrt(np.maximum(1.0 - cosine**2, 0.0))
    azimuth = 2.0 * math.pi * rng.random(n)
    reference = np.zeros_like(mean_direction)
    use_x = np.abs(mean_direction[:, 2]) > 0.9
    reference[~use_x, 2] = 1.0
    reference[use_x, 0] = 1.0
    e1 = np.cross(mean_direction, reference)
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(mean_direction, e1)
    return (
        cosine[:, None] * mean_direction
        + sine[:, None]
        * (np.cos(azimuth)[:, None] * e1 + np.sin(azimuth)[:, None] * e2)
    )


def _von_mises_fisher_density(
    direction: np.ndarray,
    mean_direction: np.ndarray,
    concentration: np.ndarray,
) -> np.ndarray:
    dot = np.clip(np.sum(direction * mean_direction, axis=1), -1.0, 1.0)
    log_density = np.empty_like(concentration)
    large = concentration > 50.0
    log_density[large] = (
        np.log(concentration[large] / (2.0 * math.pi))
        + concentration[large] * (dot[large] - 1.0)
    )
    small = ~large
    log_density[small] = (
        np.log(concentration[small])
        - math.log(4.0 * math.pi)
        - np.log(np.sinh(concentration[small]))
        + concentration[small] * dot[small]
    )
    return np.exp(np.clip(log_density, -745.0, 700.0))


def scatter_equal_mass_importance(
    v1: np.ndarray,
    v2: np.ndarray,
    w_kms: float,
    radius_pc: float,
    j_lc_pc_kms: float,
    rng: np.random.Generator,
    cap_fraction: float = 0.9,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Scatter with a forward-focused, exactly weighted loss-cone proposal."""

    if not (0.0 <= cap_fraction < 1.0):
        raise ValueError("cap_fraction must lie in [0,1)")
    relative = v1 - v2
    speed = np.linalg.norm(relative, axis=1)
    safe_speed = np.maximum(speed, np.finfo(float).tiny)
    direction = relative / safe_speed[:, None]
    centre = 0.5 * (v1 + v2)
    half_speed = 0.5 * safe_speed
    a = (safe_speed / w_kms) ** 2
    cap_radius = j_lc_pc_kms / (radius_pc * half_speed)
    cap_centre_1 = -centre[:, :2] / half_speed[:, None]
    cap_centre_2 = -cap_centre_1
    cap_area = circle_intersection_area(
        np.linalg.norm(cap_centre_1, axis=1), cap_radius
    )
    alpha = np.where(cap_area > 0.0, cap_fraction, 0.0)
    focus1 = _nearest_loss_cone_direction(direction, cap_centre_1, cap_radius)
    focus2 = _nearest_loss_cone_direction(direction, cap_centre_2, cap_radius)
    concentration = np.clip(
        np.maximum(0.5 * a, 2.0 / np.maximum(cap_radius, 1.0e-8) ** 2),
        1.0,
        1.0e8,
    )

    use_focus = rng.random(speed.size) < alpha
    focus_on_1 = rng.random(speed.size) < 0.5
    chosen_focus = np.where(focus_on_1[:, None], focus1, focus2)
    post_direction = np.empty_like(direction)
    focused_ids = np.flatnonzero(use_focus)
    if focused_ids.size:
        post_direction[focused_ids] = _sample_von_mises_fisher(
            chosen_focus[focused_ids], concentration[focused_ids], rng
        )
    baseline_ids = np.flatnonzero(~use_focus)
    if baseline_ids.size:
        u = sample_tchannel_u(a[baseline_ids], rng)
        cos_theta = 1.0 - 2.0 * u
        sin_theta = 2.0 * np.sqrt(np.maximum(u * (1.0 - u), 0.0))
        azimuth = 2.0 * math.pi * rng.random(baseline_ids.size)
        reference = np.zeros((baseline_ids.size, 3))
        use_x = np.abs(direction[baseline_ids, 2]) > 0.9
        reference[~use_x, 2] = 1.0
        reference[use_x, 0] = 1.0
        e1 = np.cross(direction[baseline_ids], reference)
        e1 /= np.linalg.norm(e1, axis=1)[:, None]
        e2 = np.cross(direction[baseline_ids], e1)
        post_direction[baseline_ids] = (
            cos_theta[:, None] * direction[baseline_ids]
            + sin_theta[:, None]
            * (
                np.cos(azimuth)[:, None] * e1
                + np.sin(azimuth)[:, None] * e2
            )
        )

    dot = np.sum(direction * post_direction, axis=1)
    u = np.clip(0.5 * (1.0 - dot), 0.0, 1.0)
    qphysical = (1.0 + a) / (4.0 * math.pi * (1.0 + a * u) ** 2)
    qfocus = 0.5 * (
        _von_mises_fisher_density(post_direction, focus1, concentration)
        + _von_mises_fisher_density(post_direction, focus2, concentration)
    )
    qproposal = (1.0 - alpha) * qphysical + alpha * qfocus
    importance_weight = qphysical / qproposal
    relative_after = safe_speed[:, None] * post_direction
    return (
        centre + 0.5 * relative_after,
        centre - 0.5 * relative_after,
        1.0 / (1.0 + a),
        importance_weight,
    )


def angular_mean_v_sigma_total_ratio(
    test_speed_kms: np.ndarray | float,
    background_speed_kms: np.ndarray | float,
    w_kms: float,
) -> np.ndarray | float:
    """Return the isotropic angular mean of ``V sigma_tot(V)/sigma0``.

    For the t-channel kernel ``sigma_tot/sigma0 = 1/[1+(V/w)^2]``.
    The angular integral is analytic, which removes Monte-Carlo noise from the
    post-jump survival probability.
    """

    u, v = np.broadcast_arrays(
        np.asarray(test_speed_kms, float),
        np.asarray(background_speed_kms, float),
    )
    if w_kms <= 0.0 or np.any(u < 0.0) or np.any(v < 0.0):
        raise ValueError("speeds and w_kms must be non-negative, with w_kms > 0")
    product = u * v
    slow = product <= 1.0e-14 * max(w_kms * w_kms, 1.0)
    out = np.empty_like(product)
    relative = u + v
    out[slow] = relative[slow] / (1.0 + (relative[slow] / w_kms) ** 2)
    if np.any(~slow):
        lo = np.abs(u[~slow] - v[~slow])
        hi = u[~slow] + v[~slow]

        def primitive(s: np.ndarray) -> np.ndarray:
            return s - w_kms * np.arctan(s / w_kms)

        out[~slow] = (
            w_kms**2
            * (primitive(hi) - primitive(lo))
            / (2.0 * product[~slow])
        )
    return float(out) if out.ndim == 0 else out


def local_total_collision_rates_per_myr(
    radius_pc: np.ndarray | float,
    binding_energy_kms2: np.ndarray | float,
    rh_pc: float,
    v0_kms: float,
    n0_pc3: float,
    particle_mass_msun: float,
    sigma0_over_m_cm2_g: float,
    w_kms: float,
    xgrid: np.ndarray,
    gx: np.ndarray,
    quadrature_order: int = 96,
) -> np.ndarray | float:
    """Total collision rates for test particles in the isotropic DF."""

    radius, binding = np.broadcast_arrays(
        np.asarray(radius_pc, float), np.asarray(binding_energy_kms2, float)
    )
    shape = radius.shape
    radius = radius.reshape(-1)
    binding = binding.reshape(-1)
    if np.any(radius <= 0.0) or np.any(binding <= 0.0):
        raise ValueError("radius and binding energy must be positive")
    psi = rh_pc / radius
    test_speed2 = 2.0 * (v0_kms**2 * psi - binding)
    if np.any(test_speed2 < -1.0e-8 * v0_kms**2):
        raise ValueError("test orbit cannot reach the requested radius")
    test_speed = np.sqrt(np.maximum(test_speed2, 0.0))
    upper_x = np.minimum(float(xgrid[-1]), psi)
    if np.any(upper_x <= xgrid[0]):
        raise RuntimeError("post-jump orbit leaves the represented DF support")
    # x=psi-y^2 removes the square-root endpoint from the local density
    # integral.  One high-order rule then covers the whole tabulated DF and is
    # substantially faster than revisiting every energy interval for every
    # point along every candidate orbit.
    ylo = np.sqrt(np.maximum(psi - upper_x, 0.0))
    yhi = np.sqrt(np.maximum(psi - float(xgrid[0]), 0.0))
    nodes, weights = _leggauss(quadrature_order)
    yy = (
        0.5 * (yhi - ylo)[:, None] * nodes[None, :]
        + 0.5 * (yhi + ylo)[:, None]
    )
    ww = 0.5 * (yhi - ylo)[:, None] * weights[None, :]
    xx = np.clip(
        psi[:, None] - yy * yy,
        float(xgrid[0]),
        upper_x[:, None],
    )
    g = np.interp(
        np.log(xx).reshape(-1), np.log(xgrid), gx
    ).reshape(xx.shape)
    # sqrt(psi-x) dx = 2 y^2 dy.
    phase_weight = 2.0 * ww * g * yy**2
    background_speed = v0_kms * np.sqrt(2.0) * yy
    rate_integral = np.sum(
        phase_weight
        * angular_mean_v_sigma_total_ratio(
            test_speed[:, None], background_speed, w_kms
        ),
        axis=1,
    )
    density_integral = np.sum(phase_weight, axis=1)
    if np.any(density_integral <= 0.0):
        raise RuntimeError("post-jump orbit leaves the represented DF support")
    density_factor = 2.0 / math.sqrt(math.pi) * density_integral
    rho = particle_mass_msun * n0_pc3 * density_factor
    mean_v_sigma_ratio = rate_integral / density_integral
    som = sigma0_over_m_cm2_g * CM2_G_TO_PC2_MSUN
    result = rho * som * mean_v_sigma_ratio * KMS_TO_PC_PER_MYR
    result = result.reshape(shape)
    return float(result) if result.ndim == 0 else result


def local_total_collision_rate_per_myr(
    radius_pc: float,
    binding_energy_kms2: float,
    rh_pc: float,
    v0_kms: float,
    n0_pc3: float,
    particle_mass_msun: float,
    sigma0_over_m_cm2_g: float,
    w_kms: float,
    xgrid: np.ndarray,
    gx: np.ndarray,
    quadrature_order: int = 96,
) -> float:
    """Scalar wrapper for ``local_total_collision_rates_per_myr``."""

    return float(local_total_collision_rates_per_myr(
        radius_pc,
        binding_energy_kms2,
        rh_pc,
        v0_kms,
        n0_pc3,
        particle_mass_msun,
        sigma0_over_m_cm2_g,
        w_kms,
        xgrid,
        gx,
        quadrature_order,
    ))


def orbit_survival_probability(
    radius_pc: float,
    velocity_kms: np.ndarray,
    rh_pc: float,
    mbh_msun: float,
    n0_pc3: float,
    particle_mass_msun: float,
    sigma0_over_m_cm2_g: float,
    w_kms: float,
    xgrid: np.ndarray,
    gx: np.ndarray,
    anomaly_order: int = 24,
) -> tuple[float, float, float, float]:
    """Return no-rescattering survival to the orbit-dependent plunge surface."""

    velocity = np.asarray(velocity_kms, float)
    if velocity.shape != (3,):
        raise ValueError("velocity_kms must be a three-vector")
    gm = G_PC_KMS2_MSUN * mbh_msun
    binding = gm / radius_pc - 0.5 * float(np.dot(velocity, velocity))
    angular_momentum = radius_pc * float(np.linalg.norm(velocity[:2]))
    jlc = 4.0 * gm / C_KMS
    if binding <= 0.0 or angular_momentum >= jlc:
        raise ValueError("survival is defined only for bound loss-cone states")
    semimajor = gm / (2.0 * binding)
    eccentricity2 = 1.0 - angular_momentum**2 / (gm * semimajor)
    if eccentricity2 < -1.0e-10:
        raise ValueError("post-jump state does not define a bound Kepler orbit")
    eccentricity = math.sqrt(max(eccentricity2, 0.0))
    apocentre = semimajor * (1.0 + eccentricity)
    rg = gm / C_KMS**2
    capture_radius = 16.0 * rg / (1.0 + eccentricity)
    if radius_pc <= capture_radius * (1.0 + 1.0e-12):
        return 1.0, 0.0, apocentre, capture_radius

    def anomaly_at_radius(r: float) -> float:
        if eccentricity <= 1.0e-14:
            return 0.0
        cosine = (1.0 - r / semimajor) / eccentricity
        return math.acos(float(np.clip(cosine, -1.0, 1.0)))

    current_anomaly = anomaly_at_radius(radius_pc)
    capture_anomaly = anomaly_at_radius(capture_radius)
    if capture_anomaly > current_anomaly + 1.0e-10:
        raise ValueError("collision occurred inside the orbit-specific plunge surface")
    v0 = math.sqrt(gm / rh_pc)
    time_scale_myr = math.sqrt(semimajor**3 / gm) / KMS_TO_PC_PER_MYR
    qnodes, qweights = _leggauss(anomaly_order)

    def optical_depth_between(lo: float, hi: float) -> float:
        if hi <= lo:
            return 0.0
        anomaly = 0.5 * (hi - lo) * qnodes + 0.5 * (hi + lo)
        weight = 0.5 * (hi - lo) * qweights
        radii = semimajor * (1.0 - eccentricity * np.cos(anomaly))
        dt_danomaly = time_scale_myr * (
            1.0 - eccentricity * np.cos(anomaly)
        )
        rates = np.array([
            local_total_collision_rate_per_myr(
                float(r), binding, rh_pc, v0, n0_pc3,
                particle_mass_msun, sigma0_over_m_cm2_g, w_kms,
                xgrid, gx,
            )
            for r in radii
        ])
        return float(np.sum(weight * dt_danomaly * rates))

    inner_leg = optical_depth_between(capture_anomaly, current_anomaly)
    if velocity[2] < 0.0:
        optical_depth = inner_leg
    else:
        outward_leg = optical_depth_between(current_anomaly, math.pi)
        return_leg = optical_depth_between(capture_anomaly, math.pi)
        optical_depth = outward_leg + return_leg
    return (
        math.exp(-min(optical_depth, 745.0)),
        optical_depth,
        apocentre,
        capture_radius,
    )


def orbit_survival_probabilities(
    radius_pc: float,
    velocities_kms: np.ndarray,
    rh_pc: float,
    mbh_msun: float,
    n0_pc3: float,
    particle_mass_msun: float,
    sigma0_over_m_cm2_g: float,
    w_kms: float,
    xgrid: np.ndarray,
    gx: np.ndarray,
    anomaly_order: int = 16,
    collision_quadrature_order: int = 48,
    max_apocentre_pc: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized no-rescattering probabilities for loss-cone states."""

    velocity = np.asarray(velocities_kms, float)
    if velocity.ndim != 2 or velocity.shape[1] != 3:
        raise ValueError("velocities_kms must have shape (N,3)")
    if velocity.shape[0] == 0:
        empty = np.empty(0)
        return empty, empty, empty, empty
    gm = G_PC_KMS2_MSUN * mbh_msun
    binding = gm / radius_pc - 0.5 * np.sum(velocity * velocity, axis=1)
    angular_momentum = radius_pc * np.linalg.norm(velocity[:, :2], axis=1)
    jlc = 4.0 * gm / C_KMS
    if np.any(binding <= 0.0) or np.any(angular_momentum >= jlc):
        raise ValueError("survival is defined only for bound loss-cone states")
    semimajor = gm / (2.0 * binding)
    eccentricity2 = 1.0 - angular_momentum**2 / (gm * semimajor)
    if np.any(eccentricity2 < -1.0e-10):
        raise ValueError("post-jump state does not define a bound Kepler orbit")
    eccentricity = np.sqrt(np.maximum(eccentricity2, 0.0))
    apocentre = semimajor * (1.0 + eccentricity)
    rg = gm / C_KMS**2
    capture_radius = 16.0 * rg / (1.0 + eccentricity)
    supported = (
        np.ones(velocity.shape[0], dtype=bool)
        if max_apocentre_pc is None
        else apocentre <= max_apocentre_pc
    )
    immediate = supported & (
        radius_pc <= capture_radius * (1.0 + 1.0e-12)
    )
    active = supported & ~immediate
    if not np.all(active):
        survival = np.zeros(velocity.shape[0])
        optical_depth = np.full(velocity.shape[0], math.inf)
        survival[immediate] = 1.0
        optical_depth[immediate] = 0.0
        if np.any(active):
            sub = orbit_survival_probabilities(
                radius_pc,
                velocity[active],
                rh_pc,
                mbh_msun,
                n0_pc3,
                particle_mass_msun,
                sigma0_over_m_cm2_g,
                w_kms,
                xgrid,
                gx,
                anomaly_order=anomaly_order,
                collision_quadrature_order=collision_quadrature_order,
                max_apocentre_pc=None,
            )
            survival[active] = sub[0]
            optical_depth[active] = sub[1]
        return survival, optical_depth, apocentre, capture_radius

    def anomaly_at_radius(r: np.ndarray | float) -> np.ndarray:
        radius = np.broadcast_to(np.asarray(r, float), eccentricity.shape)
        cosine = np.divide(
            1.0 - radius / semimajor,
            eccentricity,
            out=np.ones_like(eccentricity),
            where=eccentricity > 1.0e-14,
        )
        return np.arccos(np.clip(cosine, -1.0, 1.0))

    current_anomaly = anomaly_at_radius(radius_pc)
    capture_anomaly = anomaly_at_radius(capture_radius)
    if np.any(capture_anomaly > current_anomaly + 1.0e-10):
        raise ValueError("collision occurred inside the orbit-specific plunge surface")
    v0 = math.sqrt(gm / rh_pc)
    time_scale_myr = np.sqrt(semimajor**3 / gm) / KMS_TO_PC_PER_MYR
    qnodes, qweights = _leggauss(anomaly_order)

    def optical_depth_between(lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
        width = np.maximum(hi - lo, 0.0)
        anomaly = (
            0.5 * width[:, None] * qnodes[None, :]
            + 0.5 * (hi + lo)[:, None]
        )
        weight = 0.5 * width[:, None] * qweights[None, :]
        factor = 1.0 - eccentricity[:, None] * np.cos(anomaly)
        radii = semimajor[:, None] * factor
        rates = local_total_collision_rates_per_myr(
            radii.reshape(-1),
            np.repeat(binding, anomaly_order),
            rh_pc,
            v0,
            n0_pc3,
            particle_mass_msun,
            sigma0_over_m_cm2_g,
            w_kms,
            xgrid,
            gx,
            quadrature_order=collision_quadrature_order,
        ).reshape(radii.shape)
        return np.sum(
            weight * time_scale_myr[:, None] * factor * rates,
            axis=1,
        )

    inner_leg = optical_depth_between(capture_anomaly, current_anomaly)
    outward = velocity[:, 2] >= 0.0
    optical_depth = inner_leg.copy()
    if np.any(outward):
        pi = np.full_like(current_anomaly, math.pi)
        outward_leg = optical_depth_between(current_anomaly, pi)
        return_leg = optical_depth_between(capture_anomaly, pi)
        optical_depth[outward] = (
            outward_leg[outward] + return_leg[outward]
        )
    return (
        np.exp(-np.minimum(optical_depth, 745.0)),
        optical_depth,
        apocentre,
        capture_radius,
    )


def _cdf_from_raw(z: np.ndarray, raw: np.ndarray) -> tuple[np.ndarray, float]:
    cdf = np.zeros_like(z)
    cdf[1:] = np.cumsum(
        0.5 * (raw[:-1] + raw[1:]) * np.diff(z)
    )
    norm = float(cdf[-1])
    if not np.isfinite(norm) or norm <= 0.0:
        raise RuntimeError("empty local velocity distribution outside the loss cone")
    return cdf / norm, norm


def jump_balance_tensor(
    pre_y1: np.ndarray,
    pre_y2: np.ndarray,
    post_y1: np.ndarray,
    post_y2: np.ndarray,
    post_bound1: np.ndarray,
    post_bound2: np.ndarray,
    capture1: np.ndarray,
    capture2: np.ndarray,
    survival1: np.ndarray,
    survival2: np.ndarray,
    event_weight: np.ndarray,
    edges: np.ndarray = JUMP_J_OVER_JLC_EDGES,
) -> tuple[np.ndarray, np.ndarray]:
    """Accumulate the finite-jump population balance in J/J_lc bins.

    A captured post-collision state removes its pre-collision population with
    the no-rescattering probability.  The complementary probability is
    returned to the source bin, a conservative closure for the unmodelled
    transient inside the loss cone.  Unbound outcomes are likewise returned,
    so the tensor isolates black-hole capture rather than evaporation.
    """

    arrays = np.broadcast_arrays(
        pre_y1, pre_y2, post_y1, post_y2, post_bound1, post_bound2,
        capture1, capture2, survival1, survival2, event_weight,
    )
    (
        pre_y1, pre_y2, post_y1, post_y2, post_bound1, post_bound2,
        capture1, capture2, survival1, survival2, event_weight,
    ) = [np.asarray(x) for x in arrays]
    nbins = edges.size - 1

    def bins(y: np.ndarray) -> np.ndarray:
        index = np.searchsorted(edges, y, side="right") - 1
        return np.clip(index, 0, nbins - 1)

    pre1 = bins(pre_y1)
    pre2 = bins(pre_y2)
    post1 = bins(post_y1)
    post2 = bins(post_y2)
    pair = pre1 * nbins + pre2
    tensor = np.zeros((nbins * nbins, nbins))
    capture = np.zeros(nbins * nbins)

    def add_particle(
        pre: np.ndarray,
        post: np.ndarray,
        bound: np.ndarray,
        captured: np.ndarray,
        survival: np.ndarray,
    ) -> None:
        outside = bound & ~captured
        np.add.at(tensor, (pair[outside], pre[outside]), -event_weight[outside])
        np.add.at(tensor, (pair[outside], post[outside]), event_weight[outside])
        captured_weight = event_weight * survival
        np.add.at(tensor, (pair[captured], pre[captured]), -captured_weight[captured])
        np.add.at(capture, pair[captured], captured_weight[captured])

    add_particle(pre1, post1, post_bound1, capture1, survival1)
    add_particle(pre2, post2, post_bound2, capture2, survival2)
    return tensor.reshape(nbins, nbins, nbins), capture.reshape(nbins, nbins)


def kepler_apocentre_pc(
    binding_energy_kms2: np.ndarray,
    angular_momentum_pc_kms: np.ndarray,
    gm_pc_kms2: float,
) -> np.ndarray:
    """Return the Kepler apocentre for bound phase-space points."""

    binding, angular_momentum = np.broadcast_arrays(
        np.asarray(binding_energy_kms2, float),
        np.asarray(angular_momentum_pc_kms, float),
    )
    out = np.full(binding.shape, np.inf)
    bound = binding > 0.0
    if not np.any(bound):
        return out
    semimajor = gm_pc_kms2 / (2.0 * binding[bound])
    eccentricity2 = 1.0 - angular_momentum[bound] ** 2 / (
        gm_pc_kms2 * semimajor
    )
    if np.any(eccentricity2 < -1.0e-8):
        raise ValueError("phase-space point exceeds the circular angular momentum")
    eccentricity = np.sqrt(np.maximum(eccentricity2, 0.0))
    out[bound] = semimajor * (1.0 + eccentricity)
    return out


def ej_state_indices(
    binding_energy_kms2: np.ndarray,
    angular_momentum_pc_kms: np.ndarray,
    gm_pc_kms2: float,
    v0_kms: float,
    j_lc_pc_kms: float,
    reservoir_radius_pc: float,
    x_edges: np.ndarray,
    y_edges: np.ndarray = JUMP_J_OVER_JLC_EDGES,
) -> np.ndarray:
    """Map bound states outside the loss cone to ``(E,J,domain)`` bins.

    Every energy--angular-momentum cell is split into an internal state and a
    reservoir-connected state.  The split is evaluated from the exact Kepler
    apocentre of each sampled state, so no coarse bin can straddle the handoff.
    States below the represented weak-binding range are assigned directly to
    the outer reservoir.
    """

    binding, angular_momentum = np.broadcast_arrays(
        np.asarray(binding_energy_kms2, float),
        np.asarray(angular_momentum_pc_kms, float),
    )
    x_edges = np.asarray(x_edges, float)
    y_edges = np.asarray(y_edges, float)
    if (
        x_edges.ndim != 1
        or x_edges.size < 3
        or np.any(np.diff(x_edges) <= 0.0)
        or y_edges.ndim != 1
        or y_edges.size < 3
        or np.any(np.diff(y_edges) <= 0.0)
    ):
        raise ValueError("EJ bin edges must be strictly increasing")
    n_energy = x_edges.size - 1
    n_angular = y_edges.size - 1
    state = np.full(binding.shape, EJ_POST_RESERVOIR, dtype=np.int64)
    x = binding / v0_kms**2
    y = angular_momentum / j_lc_pc_kms
    represented = (
        (binding > 0.0)
        & (y >= y_edges[0] * (1.0 - 2.0e-10))
        & (x >= x_edges[0])
        & (x <= x_edges[-1] * (1.0 + 1.0e-12))
    )
    if not np.any(represented):
        return state
    energy_bin = np.searchsorted(x_edges, x[represented], side="right") - 1
    energy_bin = np.clip(energy_bin, 0, n_energy - 1)
    angular_bin = np.searchsorted(
        y_edges, np.maximum(y[represented], y_edges[0]), side="right"
    ) - 1
    angular_bin = np.clip(angular_bin, 0, n_angular - 1)
    apocentre = kepler_apocentre_pc(
        binding[represented], angular_momentum[represented], gm_pc_kms2
    )
    connected = (apocentre >= reservoir_radius_pc).astype(np.int64)
    state[represented] = 2 * (
        energy_bin * n_angular + angular_bin
    ) + connected
    return state


def aggregate_ej_jump_records(
    source1: np.ndarray,
    source2: np.ndarray,
    pre_state: np.ndarray,
    post_state: np.ndarray,
    direct_rate: np.ndarray,
    immediate_rate: np.ndarray,
    sample_count: Optional[np.ndarray] = None,
    direct_binding_energy_rate: Optional[np.ndarray] = None,
    immediate_binding_energy_rate: Optional[np.ndarray] = None,
    nstates: int = 0,
) -> dict[str, np.ndarray]:
    """Aggregate sparse two-particle jump records without losing currents.

    The optional binding-energy columns carry capture-rate times the positive
    post-collision binding energy.  They are zero for ordinary phase-space
    transitions and reservoir exchange.  Keeping them on the same sparse keys
    ensures that the steady mass and capture-energy currents are evaluated
    with exactly the same depleted occupations.
    """

    arrays = np.broadcast_arrays(
        source1, source2, pre_state, post_state, direct_rate, immediate_rate
    )
    source1, source2, pre_state, post_state = [
        np.asarray(value, np.int64).reshape(-1) for value in arrays[:4]
    ]
    direct_rate, immediate_rate = [
        np.asarray(value, float).reshape(-1) for value in arrays[4:]
    ]
    if nstates <= 0:
        raise ValueError("the EJ operator needs a positive state count")
    if sample_count is None:
        sample_count = np.ones(source1.size, dtype=np.int64)
    else:
        sample_count = np.asarray(sample_count, np.int64).reshape(-1)
    if direct_binding_energy_rate is None:
        direct_binding_energy_rate = np.zeros(source1.size, dtype=float)
    else:
        direct_binding_energy_rate = np.asarray(
            direct_binding_energy_rate, float
        ).reshape(-1)
    if immediate_binding_energy_rate is None:
        immediate_binding_energy_rate = np.zeros(source1.size, dtype=float)
    else:
        immediate_binding_energy_rate = np.asarray(
            immediate_binding_energy_rate, float
        ).reshape(-1)
    if not all(array.size == source1.size for array in (
        source2, pre_state, post_state, direct_rate, immediate_rate, sample_count,
        direct_binding_energy_rate, immediate_binding_energy_rate,
    )):
        raise ValueError("sparse EJ record arrays have inconsistent lengths")
    valid = (
        (source1 >= 0) & (source1 < nstates)
        & (source2 >= 0) & (source2 < nstates)
        & (pre_state >= 0) & (pre_state < nstates)
        & (post_state >= EJ_POST_RESERVOIR) & (post_state < nstates)
        & np.isfinite(direct_rate) & np.isfinite(immediate_rate)
        & (direct_rate >= 0.0) & (immediate_rate >= 0.0)
        & np.isfinite(direct_binding_energy_rate)
        & np.isfinite(immediate_binding_energy_rate)
        & (direct_binding_energy_rate >= 0.0)
        & (immediate_binding_energy_rate >= 0.0)
        & (sample_count > 0)
    )
    if not np.all(valid):
        raise ValueError("invalid sparse EJ jump record")
    post_span = nstates + 2
    key = (
        ((source1 * nstates + source2) * nstates + pre_state)
        * post_span + (post_state + 2)
    )
    unique, inverse = np.unique(key, return_inverse=True)
    direct_sum = np.bincount(inverse, weights=direct_rate)
    immediate_sum = np.bincount(inverse, weights=immediate_rate)
    direct_energy_sum = np.bincount(
        inverse, weights=direct_binding_energy_rate
    )
    immediate_energy_sum = np.bincount(
        inverse, weights=immediate_binding_energy_rate
    )
    count_sum = np.bincount(inverse, weights=sample_count).astype(np.int64)
    post_shifted = unique % post_span
    reduced = unique // post_span
    pre = reduced % nstates
    reduced //= nstates
    source_b = reduced % nstates
    source_a = reduced // nstates
    return {
        "source1": source_a.astype(np.int64),
        "source2": source_b.astype(np.int64),
        "pre": pre.astype(np.int64),
        "post": (post_shifted - 2).astype(np.int64),
        "rate_direct_msun_per_myr": direct_sum,
        "rate_immediate_msun_per_myr": immediate_sum,
        "rate_direct_binding_msun_kms2_per_myr": direct_energy_sum,
        "rate_immediate_binding_msun_kms2_per_myr": immediate_energy_sum,
        "sample_count": count_sum,
    }


def ej_jump_records(
    pre_binding1: np.ndarray,
    pre_binding2: np.ndarray,
    pre_j1: np.ndarray,
    pre_j2: np.ndarray,
    post_binding1: np.ndarray,
    post_binding2: np.ndarray,
    post_j1: np.ndarray,
    post_j2: np.ndarray,
    capture1: np.ndarray,
    capture2: np.ndarray,
    survival1: np.ndarray,
    survival2: np.ndarray,
    event_weight: np.ndarray,
    gm_pc_kms2: float,
    v0_kms: float,
    j_lc_pc_kms: float,
    reservoir_radius_pc: float,
    x_edges: np.ndarray,
    include_nonkepler_exchange: bool = True,
    y_edges: np.ndarray = JUMP_J_OVER_JLC_EDGES,
    energy_order: int = 0,
) -> dict[str, np.ndarray]:
    """Construct a sparse non-local collision operator in energy and angular momentum.

    ``rate_direct`` removes only the fraction that reaches the plunge surface
    before another collision. ``rate_immediate`` applies an absorbing loss cone
    to every injected state. Ordinary bound transitions are identical in both
    models. Unbound post-collision states can be transferred to the outer fluid
    reservoir instead of being silently returned to their source bin.
    """

    arrays = np.broadcast_arrays(
        pre_binding1, pre_binding2, pre_j1, pre_j2,
        post_binding1, post_binding2, post_j1, post_j2,
        capture1, capture2, survival1, survival2, event_weight,
    )
    (
        pre_binding1, pre_binding2, pre_j1, pre_j2,
        post_binding1, post_binding2, post_j1, post_j2,
        capture1, capture2, survival1, survival2, event_weight,
    ) = [np.asarray(value) for value in arrays]
    nstates = 2 * (len(x_edges) - 1) * (len(y_edges) - 1)
    source1 = ej_state_indices(
        pre_binding1, pre_j1, gm_pc_kms2, v0_kms, j_lc_pc_kms,
        reservoir_radius_pc, x_edges, y_edges,
    )
    source2 = ej_state_indices(
        pre_binding2, pre_j2, gm_pc_kms2, v0_kms, j_lc_pc_kms,
        reservoir_radius_pc, x_edges, y_edges,
    )
    if np.any(source1 < 0) or np.any(source2 < 0):
        bad_binding = np.concatenate([
            np.asarray(pre_binding1)[source1 < 0],
            np.asarray(pre_binding2)[source2 < 0],
        ]) / v0_kms**2
        bad_j = np.concatenate([
            np.asarray(pre_j1)[source1 < 0],
            np.asarray(pre_j2)[source2 < 0],
        ]) / j_lc_pc_kms
        raise RuntimeError(
            "sampled source states lie outside the EJ operator grid: "
            f"count={bad_binding.size}, x_range="
            f"[{np.min(bad_binding):.6e},{np.max(bad_binding):.6e}], "
            f"J/J_lc_range=[{np.min(bad_j):.6e},{np.max(bad_j):.6e}], "
            f"grid=[{x_edges[0]:.6e},{x_edges[-1]:.6e}]"
        )

    # Separate sidecar preserves self-cell energy exchange without changing
    # the population operator or its convergence normalization.
    post_states = [ej_state_indices(
        binding, j, gm_pc_kms2, v0_kms, j_lc_pc_kms,
        reservoir_radius_pc, x_edges, y_edges,
    ) for binding, j in ((post_binding1, post_j1), (post_binding2, post_j2))]
    energy_records = ej_energy_ledger.build_records(
        source1, source2, (pre_binding1, pre_binding2),
        (post_binding1, post_binding2), post_states, (capture1, capture2),
        (survival1, survival2), event_weight, nstates,
        nonkepler_exchange=include_nonkepler_exchange,
    )
    if energy_order not in (0, 1):
        raise ValueError("energy_order must be zero or one")
    if energy_order == 1:
        lo, hi = ej_energy_dg.state_edges(x_edges, len(y_edges)-1, v0_kms**2)
        energy_records.update(ej_energy_dg.build_records(
            source1, source2, (pre_binding1, pre_binding2),
            (post_binding1, post_binding2), post_states, (capture1, capture2),
            (survival1, survival2), event_weight, lo, hi,
            nonkepler_exchange=include_nonkepler_exchange))
    chunks: list[tuple[np.ndarray, ...]] = []
    for pre, post_binding, post_j, captured, survival in (
        (source1, post_binding1, post_j1, capture1, survival1),
        (source2, post_binding2, post_j2, capture2, survival2),
    ):
        bound_outside = (post_binding > 0.0) & ~captured
        post = ej_state_indices(
            post_binding, post_j, gm_pc_kms2, v0_kms, j_lc_pc_kms,
            reservoir_radius_pc, x_edges, y_edges,
        )
        ordinary = bound_outside & (post != pre)
        if np.any(ordinary):
            weight = event_weight[ordinary]
            zero = np.zeros(np.count_nonzero(ordinary), dtype=float)
            chunks.append((
                source1[ordinary], source2[ordinary], pre[ordinary], post[ordinary],
                weight, weight, np.ones(np.count_nonzero(ordinary), dtype=np.int64),
                zero, zero,
            ))
        captured = np.asarray(captured, bool)
        if np.any(captured):
            direct = event_weight[captured] * np.asarray(survival, float)[captured]
            immediate = event_weight[captured]
            keep = (direct > 0.0) | (immediate > 0.0)
            ids = np.flatnonzero(captured)[keep]
            binding = np.asarray(post_binding, float)[ids]
            chunks.append((
                source1[ids], source2[ids], pre[ids],
                np.full(ids.size, EJ_POST_CAPTURE, dtype=np.int64),
                direct[keep], immediate[keep], np.ones(ids.size, dtype=np.int64),
                direct[keep] * binding, immediate[keep] * binding,
            ))
        if include_nonkepler_exchange:
            nonkepler = post_binding <= 0.0
            if np.any(nonkepler):
                weight = event_weight[nonkepler]
                zero = np.zeros(np.count_nonzero(nonkepler), dtype=float)
                chunks.append((
                    source1[nonkepler], source2[nonkepler], pre[nonkepler],
                    np.full(np.count_nonzero(nonkepler), EJ_POST_RESERVOIR, dtype=np.int64),
                    weight, weight,
                    np.ones(np.count_nonzero(nonkepler), dtype=np.int64),
                    zero, zero,
                ))
    if not chunks:
        empty_i = np.empty(0, dtype=np.int64)
        empty_f = np.empty(0, dtype=float)
        return {
            "source1": empty_i, "source2": empty_i, "pre": empty_i,
            "post": empty_i, "rate_direct_msun_per_myr": empty_f,
            "rate_immediate_msun_per_myr": empty_f, "sample_count": empty_i,
            "rate_direct_binding_msun_kms2_per_myr": empty_f,
            "rate_immediate_binding_msun_kms2_per_myr": empty_f,
            **energy_records,
        }
    concatenated = [np.concatenate(items) for items in zip(*chunks)]
    return {**aggregate_ej_jump_records(*concatenated, nstates=nstates), **energy_records}


def sample_outside_radial_cosine(
    cmax: np.ndarray,
    rng: np.random.Generator,
    edge_fraction: float,
    edge_power: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample an outside-loss-cone direction with exact edge importance."""

    cmax = np.asarray(cmax, float)
    if (
        np.any(cmax < 0.0)
        or np.any(cmax > 1.0)
        or not (0.0 <= edge_fraction < 1.0)
        or edge_power < 1.0
    ):
        raise ValueError("invalid outside-loss-cone direction proposal")
    use_edge = rng.random(cmax.size) < edge_fraction
    t = rng.random(cmax.size)
    t[use_edge] **= edge_power
    qedge = (
        np.power(np.maximum(t, np.finfo(float).tiny), 1.0 / edge_power - 1.0)
        / edge_power
    )
    proposal_over_target = (1.0 - edge_fraction) + edge_fraction * qedge
    importance = 1.0 / proposal_over_target
    sign = np.where(rng.random(cmax.size) < 0.5, -1.0, 1.0)
    return sign * cmax * (1.0 - t), importance


def sample_outside_loss_cone(
    n: int,
    radius_pc: float,
    rh_pc: float,
    v0_kms: float,
    j_lc_pc_kms: float,
    xgrid: np.ndarray,
    gx: np.ndarray,
    rng: np.random.Generator,
    integration_points: int = 4096,
    direction_edge_fraction: float = 0.0,
    direction_edge_power: float = 4.0,
) -> tuple[np.ndarray, float, np.ndarray]:
    """Sample the local isotropic DF conditional on J > J_lc."""

    psi = rh_pc / radius_pc
    minimum_speed = j_lc_pc_kms / radius_pc
    loss_cone_energy_ceiling = psi - 0.5 * (minimum_speed / v0_kms) ** 2
    upper = min(
        float(xgrid[-1]),
        loss_cone_energy_ceiling
        - 1.0e-12 * max(abs(loss_cone_energy_ceiling), 1.0),
    )
    if upper <= xgrid[0]:
        raise RuntimeError("the DF has no bound particles at this radius")
    z = np.linspace(math.log(float(xgrid[0])), math.log(upper), integration_points)
    x = np.exp(z)
    # Match the table lookup and density forward model: g is linear in the
    # logarithmic energy coordinate, not log-linear in amplitude.
    g = np.interp(z, np.log(xgrid), gx)
    v = v0_kms * np.sqrt(np.maximum(2.0 * (psi - x), 0.0))
    ratio = np.divide(
        j_lc_pc_kms,
        radius_pc * v,
        out=np.full_like(v, np.inf),
        where=v > 0.0,
    )
    outside_direction_fraction = np.sqrt(np.maximum(1.0 - ratio**2, 0.0))
    raw_outside = (
        g * np.sqrt(np.maximum(psi - x, 0.0)) * x
        * outside_direction_fraction
    )
    cdf, outside_norm = _cdf_from_raw(z, raw_outside)
    # The conditional sampler ends at the largest energy for which an
    # outside-loss-cone direction exists.  The denominator must instead be the
    # complete local density integral.  States above the loss-cone energy
    # ceiling have too little speed to satisfy J>J_lc, but they still belong to
    # the isotropic density used to normalize the outside population.
    full_upper = min(float(xgrid[-1]), psi)
    full_z = np.linspace(
        math.log(float(xgrid[0])), math.log(full_upper), integration_points
    )
    full_x = np.exp(full_z)
    full_g = np.interp(full_z, np.log(xgrid), gx)
    full_raw = (
        full_g * np.sqrt(np.maximum(psi - full_x, 0.0)) * full_x
    )
    _, full_norm = _cdf_from_raw(full_z, full_raw)
    outside_fraction = outside_norm / full_norm
    if not (0.0 < outside_fraction <= 1.0 + 1.0e-8):
        raise RuntimeError("outside-loss-cone density fraction is not physical")
    outside_fraction = min(outside_fraction, 1.0)

    uniform = (np.arange(n, dtype=float) + rng.random(n)) / n
    rng.shuffle(uniform)
    # Remove zero-density plateaus before inverting the CDF.  Passing repeated
    # CDF values to np.interp can place a small but heavily weighted subset in
    # the forbidden J<J_lc speed interval.
    keep = np.r_[True, np.diff(cdf) > 1.0e-15]
    sampled_z = np.interp(uniform, cdf[keep], z[keep])
    sampled_x = np.exp(sampled_z)
    sampled_v = v0_kms * np.sqrt(2.0 * (psi - sampled_x))
    sampled_ratio = j_lc_pc_kms / (radius_pc * sampled_v)
    if np.any(sampled_ratio > 1.0 + 2.0e-10):
        raise RuntimeError("inverse local CDF sampled a state inside the loss cone")
    sampled_ratio = np.minimum(sampled_ratio, 1.0)
    cmax = np.sqrt(np.maximum(1.0 - sampled_ratio**2, 0.0))
    radial_cosine, direction_importance = sample_outside_radial_cosine(
        cmax,
        rng,
        direction_edge_fraction,
        direction_edge_power,
    )
    azimuth = 2.0 * math.pi * rng.random(n)
    tangential_sine = np.sqrt(np.maximum(1.0 - radial_cosine**2, 0.0))
    velocity = sampled_v[:, None] * np.column_stack([
        tangential_sine * np.cos(azimuth),
        tangential_sine * np.sin(azimuth),
        radial_cosine,
    ])
    return velocity, float(outside_fraction), direction_importance


def local_candidate_current(
    radius_pc: float,
    rho_msun_pc3: float,
    rh_pc: float,
    mbh_msun: float,
    sigma0_over_m_cm2_g: float,
    w_kms: float,
    n0_pc3: float,
    particle_mass_msun: float,
    xgrid: np.ndarray,
    gx: np.ndarray,
    pairs: int,
    rng: np.random.Generator,
    cap_importance_fraction: float = 0.5,
    direction_edge_fraction: float = 0.3,
    direction_edge_power: float = 4.0,
    survival_anomaly_order: int = 16,
    survival_collision_order: int = 48,
    ej_x_edges: Optional[np.ndarray] = None,
    ej_y_edges: Optional[np.ndarray] = None,
    ej_reservoir_radius_pc: Optional[float] = None,
    ej_include_nonkepler_exchange: bool = True,
    ej_energy_order: int = 0,
) -> dict:
    v0 = math.sqrt(G_PC_KMS2_MSUN * mbh_msun / rh_pc)
    jlc = 4.0 * G_PC_KMS2_MSUN * mbh_msun / C_KMS
    v1, fraction1, direction_weight1 = sample_outside_loss_cone(
        pairs, radius_pc, rh_pc, v0, jlc, xgrid, gx, rng,
        direction_edge_fraction=direction_edge_fraction,
        direction_edge_power=direction_edge_power,
    )
    v2, fraction2, direction_weight2 = sample_outside_loss_cone(
        pairs, radius_pc, rh_pc, v0, jlc, xgrid, gx, rng,
        direction_edge_fraction=direction_edge_fraction,
        direction_edge_power=direction_edge_power,
    )
    v1_after, v2_after, total_ratio, angular_importance = (
        scatter_equal_mass_disc_importance(
            v1, v2, w_kms, radius_pc, jlc, rng,
            cap_fraction=cap_importance_fraction,
        )
    )
    relative_speed = np.linalg.norm(v1 - v2, axis=1)
    psi_physical = G_PC_KMS2_MSUN * mbh_msun / radius_pc

    def capture_state(
        v: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        binding = psi_physical - 0.5 * np.sum(v * v, axis=1)
        angular_momentum = radius_pc * np.linalg.norm(v[:, :2], axis=1)
        captured = (binding > 0.0) & (angular_momentum < jlc)
        inward = captured & (v[:, 2] < 0.0)
        return captured, inward, binding, angular_momentum

    cap1, inward1, bind1, jafter1 = capture_state(v1_after)
    cap2, inward2, bind2, jafter2 = capture_state(v2_after)
    candidate_count = cap1.astype(float) + cap2.astype(float)
    inward_count = inward1.astype(float) + inward2.astype(float)
    survival1 = np.zeros(pairs)
    survival2 = np.zeros(pairs)
    apocentre1 = np.zeros(pairs)
    apocentre2 = np.zeros(pairs)
    optical_depths: list[float] = []
    apocentres: list[float] = []
    capture_radii: list[float] = []
    for mask, velocities, survival, apocentre_store in (
        (cap1, v1_after, survival1, apocentre1),
        (cap2, v2_after, survival2, apocentre2),
    ):
        indices = np.flatnonzero(mask)
        for start in range(0, indices.size, 256):
            batch = indices[start:start + 256]
            probability, tau, apo, rcapture = orbit_survival_probabilities(
                radius_pc,
                velocities[batch],
                rh_pc,
                mbh_msun,
                n0_pc3,
                particle_mass_msun,
                sigma0_over_m_cm2_g,
                w_kms,
                xgrid,
                gx,
                anomaly_order=survival_anomaly_order,
                collision_quadrature_order=survival_collision_order,
                max_apocentre_pc=rh_pc,
            )
            survival[batch] = probability
            apocentre_store[batch] = apo
            optical_depths.extend(tau.tolist())
            apocentres.extend(apo.tolist())
            capture_radii.extend(rcapture.tolist())
    accepted_count = survival1 + survival2
    unsupported_count = (
        (cap1 & (apocentre1 > rh_pc)).astype(float)
        + (cap2 & (apocentre2 > rh_pc)).astype(float)
    )
    collision_weight = (
        relative_speed
        * total_ratio
        * angular_importance
        * direction_weight1
        * direction_weight2
    )
    weighted_candidates = collision_weight * candidate_count
    weighted_accepted = collision_weight * accepted_count
    weighted_unsupported = collision_weight * unsupported_count
    weighted_inward = collision_weight * inward_count
    weighted_candidate_energy = collision_weight * (
        cap1 * bind1 + cap2 * bind2
    )
    weighted_accepted_energy = collision_weight * (
        survival1 * bind1 + survival2 * bind2
    )
    source_j1 = radius_pc * np.linalg.norm(v1[:, :2], axis=1) / jlc
    source_j2 = radius_pc * np.linalg.norm(v2[:, :2], axis=1) / jlc
    pre_bind1 = psi_physical - 0.5 * np.sum(v1 * v1, axis=1)
    pre_bind2 = psi_physical - 0.5 * np.sum(v2 * v2, axis=1)
    candidate_source_j_contributions = []
    accepted_source_j_contributions = []
    for lo, hi in zip(
        SOURCE_J_OVER_JLC_EDGES[:-1], SOURCE_J_OVER_JLC_EDGES[1:]
    ):
        candidate_contribution = collision_weight * (
            cap1 * ((source_j1 >= lo) & (source_j1 < hi))
            + cap2 * ((source_j2 >= lo) & (source_j2 < hi))
        )
        accepted_contribution = collision_weight * (
            survival1 * ((source_j1 >= lo) & (source_j1 < hi))
            + survival2 * ((source_j2 >= lo) & (source_j2 < hi))
        )
        candidate_source_j_contributions.append(candidate_contribution)
        accepted_source_j_contributions.append(accepted_contribution)
    candidate_source_j_contributions = np.asarray(candidate_source_j_contributions)
    accepted_source_j_contributions = np.asarray(accepted_source_j_contributions)
    candidate_source_j_hist = np.sum(candidate_source_j_contributions, axis=1)
    candidate_source_j_hist_se = np.std(
        candidate_source_j_contributions, axis=1, ddof=1
    ) / math.sqrt(pairs)
    accepted_source_j_hist = np.sum(accepted_source_j_contributions, axis=1)
    accepted_source_j_hist_se = np.std(
        accepted_source_j_contributions, axis=1, ddof=1
    ) / math.sqrt(pairs)
    candidate_source_j_floor = {}
    accepted_source_j_floor = {}
    for threshold in SOURCE_J_FLOOR_THRESHOLDS:
        candidate_contribution = collision_weight * (
            cap1 * (source_j1 >= threshold)
            + cap2 * (source_j2 >= threshold)
        )
        accepted_contribution = collision_weight * (
            survival1 * (source_j1 >= threshold)
            + survival2 * (source_j2 >= threshold)
        )
        candidate_source_j_floor[str(int(threshold))] = {
            "mean": float(np.mean(candidate_contribution)),
            "standard_error": float(
                np.std(candidate_contribution, ddof=1) / math.sqrt(pairs)
            ),
        }
        accepted_source_j_floor[str(int(threshold))] = {
            "mean": float(np.mean(accepted_contribution)),
            "standard_error": float(
                np.std(accepted_contribution, ddof=1) / math.sqrt(pairs)
            ),
        }
    outside_density = rho_msun_pc3 * math.sqrt(fraction1 * fraction2)
    som = sigma0_over_m_cm2_g * CM2_G_TO_PC2_MSUN
    prefactor = (
        0.5 * outside_density**2 * som * KMS_TO_PC_PER_MYR
    )
    jump_tensor, jump_capture = jump_balance_tensor(
        source_j1,
        source_j2,
        jafter1 / jlc,
        jafter2 / jlc,
        bind1 > 0.0,
        bind2 > 0.0,
        cap1,
        cap2,
        survival1,
        survival2,
        prefactor * collision_weight / pairs,
    )
    ej_operator = None
    if ej_x_edges is not None:
        if ej_reservoir_radius_pc is None or ej_reservoir_radius_pc <= 0.0:
            raise ValueError("an EJ operator requires a positive reservoir radius")
        ej_operator = ej_jump_records(
            pre_bind1,
            pre_bind2,
            source_j1 * jlc,
            source_j2 * jlc,
            bind1,
            bind2,
            jafter1,
            jafter2,
            cap1,
            cap2,
            survival1,
            survival2,
            prefactor * collision_weight / pairs,
            G_PC_KMS2_MSUN * mbh_msun,
            v0,
            jlc,
            ej_reservoir_radius_pc,
            np.asarray(ej_x_edges, float),
            include_nonkepler_exchange=ej_include_nonkepler_exchange,
            energy_order=ej_energy_order,
            y_edges=(
                np.asarray(ej_y_edges, float)
                if ej_y_edges is not None else JUMP_J_OVER_JLC_EDGES
            ),
        )
    mean = float(np.mean(weighted_candidates))
    standard_error = float(np.std(weighted_candidates, ddof=1) / math.sqrt(pairs))
    accepted_mean = float(np.mean(weighted_accepted))
    accepted_standard_error = float(
        np.std(weighted_accepted, ddof=1) / math.sqrt(pairs)
    )
    inward_mean = float(np.mean(weighted_inward))
    unsupported_mean = float(np.mean(weighted_unsupported))
    candidate_energy_mean = float(np.mean(weighted_candidate_energy))
    candidate_energy_standard_error = float(
        np.std(weighted_candidate_energy, ddof=1) / math.sqrt(pairs)
    )
    accepted_energy_mean = float(np.mean(weighted_accepted_energy))
    accepted_energy_standard_error = float(
        np.std(weighted_accepted_energy, ddof=1) / math.sqrt(pairs)
    )
    result = {
        "_jump_balance_tensor_rate_density": jump_tensor,
        "_jump_capture_matrix_rate_density": jump_capture,
        "radius_pc": radius_pc,
        "rho_df_msun_pc3": rho_msun_pc3,
        "outside_loss_cone_density_fraction": math.sqrt(fraction1 * fraction2),
        "raw_pairs": pairs,
        "cap_importance_fraction": cap_importance_fraction,
        "angular_importance_model": (
            "physical kernel plus projected-disc proposal, with a focused "
            "near-tangent fallback"
        ),
        "direction_edge_fraction": direction_edge_fraction,
        "direction_edge_power": direction_edge_power,
        "direction_importance_weight_mean": float(np.mean(
            direction_weight1 * direction_weight2
        )),
        "direction_importance_weight_range": [
            float(np.min(direction_weight1 * direction_weight2)),
            float(np.max(direction_weight1 * direction_weight2)),
        ],
        "survival_anomaly_order": survival_anomaly_order,
        "survival_collision_quadrature_order": survival_collision_order,
        "angular_importance_weight_mean": float(np.mean(angular_importance)),
        "angular_importance_weight_range": [
            float(np.min(angular_importance)), float(np.max(angular_importance))
        ],
        "raw_candidate_particles": int(np.sum(candidate_count)),
        "raw_inward_candidate_particles": int(np.sum(inward_count)),
        "raw_candidates_outside_kepler_domain": int(np.sum(unsupported_count)),
        "mean_v_sigma_total_over_sigma0_times_candidates_kms": mean,
        "mean_v_sigma_total_over_sigma0_times_accepted_kms": accepted_mean,
        "mean_v_sigma_total_over_sigma0_times_inward_candidates_kms": inward_mean,
        "candidate_rate_density_msun_per_pc3_per_myr": prefactor * mean,
        "candidate_rate_density_standard_error_msun_per_pc3_per_myr": prefactor * standard_error,
        "accepted_capture_rate_density_msun_per_pc3_per_myr": prefactor * accepted_mean,
        "accepted_capture_rate_density_standard_error_msun_per_pc3_per_myr": (
            prefactor * accepted_standard_error
        ),
        "inward_candidate_rate_density_msun_per_pc3_per_myr": prefactor * inward_mean,
        "unresolved_kepler_domain_rate_density_msun_per_pc3_per_myr": (
            prefactor * unsupported_mean
        ),
        "candidate_binding_energy_rate_density_msun_kms2_per_pc3_per_myr": (
            prefactor * candidate_energy_mean
        ),
        "candidate_binding_energy_rate_density_standard_error_msun_kms2_per_pc3_per_myr": (
            prefactor * candidate_energy_standard_error
        ),
        "accepted_binding_energy_rate_density_msun_kms2_per_pc3_per_myr": (
            prefactor * accepted_energy_mean
        ),
        "accepted_binding_energy_rate_density_standard_error_msun_kms2_per_pc3_per_myr": (
            prefactor * accepted_energy_standard_error
        ),
        "candidate_source_j_over_jlc_histogram_rate_density_msun_per_pc3_per_myr": (
            prefactor * candidate_source_j_hist / pairs
        ).tolist(),
        "candidate_source_j_over_jlc_histogram_rate_density_standard_error_msun_per_pc3_per_myr": (
            prefactor * candidate_source_j_hist_se
        ).tolist(),
        "accepted_capture_source_j_over_jlc_histogram_rate_density_msun_per_pc3_per_myr": (
            prefactor * accepted_source_j_hist / pairs
        ).tolist(),
        "accepted_capture_source_j_over_jlc_histogram_rate_density_standard_error_msun_per_pc3_per_myr": (
            prefactor * accepted_source_j_hist_se
        ).tolist(),
        "candidate_far_source_rate_density_msun_per_pc3_per_myr": {
            key: prefactor * value["mean"]
            for key, value in candidate_source_j_floor.items()
        },
        "candidate_far_source_rate_density_standard_error_msun_per_pc3_per_myr": {
            key: prefactor * value["standard_error"]
            for key, value in candidate_source_j_floor.items()
        },
        "accepted_capture_far_source_rate_density_msun_per_pc3_per_myr": {
            key: prefactor * value["mean"]
            for key, value in accepted_source_j_floor.items()
        },
        "accepted_capture_far_source_rate_density_standard_error_msun_per_pc3_per_myr": {
            key: prefactor * value["standard_error"]
            for key, value in accepted_source_j_floor.items()
        },
        "candidate_mean_survival_probability": (
            float(np.mean(np.concatenate([survival1[cap1], survival2[cap2]])))
            if optical_depths else None
        ),
        "candidate_optical_depth_range": (
            [float(min(optical_depths)), float(max(optical_depths))]
            if optical_depths else None
        ),
        "candidate_apocentre_range_pc": (
            [float(min(apocentres)), float(max(apocentres))]
            if apocentres else None
        ),
        "candidate_capture_radius_range_pc": (
            [float(min(capture_radii)), float(max(capture_radii))]
            if capture_radii else None
        ),
    }
    if ej_operator is not None:
        result["_ej_sparse_operator_rate_density"] = ej_operator
    return result


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--profile", type=Path, required=True)
    p.add_argument("--df-table", type=Path, required=True)
    p.add_argument("--normalization-json", type=Path, required=True)
    p.add_argument("--out-json", type=Path, required=True)
    p.add_argument("--out-csv", type=Path, required=True)
    p.add_argument("--mbh", type=float, default=4.0e6)
    p.add_argument("--sigma0-over-m", type=float, default=100.0)
    p.add_argument("--w-kms", type=float, default=80.0)
    p.add_argument("--r-outer", type=float, default=None)
    p.add_argument("--r-inner", type=float, default=None)
    p.add_argument("--radial-bins", type=int, default=48)
    p.add_argument("--pairs-per-bin", type=int, default=50000)
    p.add_argument("--cap-importance-fraction", type=float, default=0.5)
    p.add_argument("--direction-edge-fraction", type=float, default=0.3)
    p.add_argument("--direction-edge-power", type=float, default=4.0)
    p.add_argument("--survival-anomaly-order", type=int, default=16)
    p.add_argument("--survival-collision-order", type=int, default=48)
    p.add_argument("--out-ej-operator", type=Path, default=None)
    p.add_argument("--ej-energy-bins", type=int, default=12)
    p.add_argument("--ej-energy-order", type=int, choices=(0,1), default=0)
    p.add_argument(
        "--ej-angular-grid",
        choices=tuple(EJ_ANGULAR_GRIDS),
        default="default",
        help="angular resolution of the non-local EJ operator",
    )
    p.add_argument(
        "--ej-return-nonkepler-to-source",
        action="store_true",
        help=(
            "regression mode: cancel unbound post-collision outcomes instead of "
            "transferring them to the outer fluid reservoir"
        ),
    )
    p.add_argument("--seed", type=int, default=91327)
    args = p.parse_args()
    if args.radial_bins < 8 or args.pairs_per_bin < 1000:
        raise ValueError("finite-angle estimate needs at least 8 bins and 1000 pairs per bin")
    if not (0.0 <= args.cap_importance_fraction < 1.0):
        raise ValueError("cap importance fraction must lie in [0,1)")
    if not (0.0 <= args.direction_edge_fraction < 1.0):
        raise ValueError("direction edge fraction must lie in [0,1)")
    if args.direction_edge_power < 1.0:
        raise ValueError("direction edge power must be at least one")
    if args.survival_anomaly_order < 8 or args.survival_collision_order < 24:
        raise ValueError("survival quadrature orders are too small")
    if args.out_ej_operator is not None:
        if args.out_ej_operator.suffix != ".npz":
            raise ValueError("the EJ operator output must have suffix .npz")
        if args.ej_energy_bins < 4:
            raise ValueError("the EJ operator needs at least four energy bins")

    r_profile, rho_profile, sigma_profile = load_profile(args.profile)
    xgrid, _, table, table_metadata = read_df_table(args.df_table)
    gx = np.mean(table, axis=1)
    norm = json.loads(args.normalization_json.read_text())
    if norm.get("status") != "PASS" or norm.get("schema") != "gnc-normalized-df-v6":
        raise RuntimeError("DF normalization did not pass the v6 acceptance gates")
    if table_metadata.get("schema") != "gnc-normalized-df-v6":
        raise RuntimeError("DF table metadata schema is stale")
    profile_sha256 = sha256_file(args.profile)
    df_table_sha256 = sha256_file(args.df_table)
    if norm.get("profile_sha256") != profile_sha256:
        raise RuntimeError("finite-angle profile does not match the normalized DF")
    if table_metadata.get("profile_sha256") != profile_sha256:
        raise RuntimeError("DF table and finite-angle profile fingerprints differ")
    if norm.get("df_table_sha256") != df_table_sha256:
        raise RuntimeError("DF table fingerprint differs from its normalization record")
    for key in (
        "profile_sha256",
        "bridge_json_sha256",
        "mbh_msun",
        "rh_pc",
        "n0_pc3",
        "particle_mass_msun",
        "collision_kernel",
        "sigma0_over_m_cm2_g",
        "w_kms",
        "xmin",
        "xmax",
        "df_fit_inner_radius_pc",
    ):
        if table_metadata.get(key) != norm.get(key):
            raise RuntimeError(f"DF table and normalization metadata differ for {key}")
    for key, measured, expected in (
        ("mbh_msun", float(norm["mbh_msun"]), args.mbh),
        ("sigma0_over_m_cm2_g", float(norm["sigma0_over_m_cm2_g"]), args.sigma0_over_m),
        ("w_kms", float(norm["w_kms"]), args.w_kms),
    ):
        if not math.isclose(measured, expected, rel_tol=1.0e-10, abs_tol=0.0):
            raise RuntimeError(f"finite-angle metadata mismatch for {key}")
    if norm.get("collision_kernel") != "yukawa-tchannel":
        raise RuntimeError("finite-angle estimator requires the matched t-channel kernel")
    rh = float(norm["rh_pc"])
    n0 = float(norm["n0_pc3"])
    particle_mass = float(norm["particle_mass_msun"])
    rg = G_PC_KMS2_MSUN * args.mbh / C_KMS**2
    represented_inner_radius = float(norm["df_fit_inner_radius_pc"])
    r_inner = float(
        args.r_inner
        if args.r_inner is not None
        else max(8.0 * rg * (1.0 + 1.0e-10), r_profile[0], represented_inner_radius)
    )
    r_outer = float(args.r_outer if args.r_outer is not None else norm["r_boundary_pc"])
    if not (r_profile[0] <= r_inner < r_outer <= r_profile[-1]):
        raise ValueError("finite-angle radial interval lies outside the profile")

    ej_x_edges = None
    ej_y_edges = None
    ej_nstates = None
    ej_chunks: list[dict[str, np.ndarray]] = []
    if args.out_ej_operator is not None:
        v0 = math.sqrt(G_PC_KMS2_MSUN * args.mbh / rh)
        # J_c(E)=J_lc at this binding energy. More tightly bound states cannot
        # lie outside the loss cone and therefore are not source states of the
        # non-local operator.
        x_lc_limit = (C_KMS / (4.0 * math.sqrt(2.0) * v0)) ** 2
        x_upper = min(float(xgrid[-1]), x_lc_limit * (1.0 - 1.0e-12))
        if x_upper <= float(xgrid[0]):
            raise RuntimeError("the DF has no represented EJ source domain")
        ej_x_edges = np.asarray([
            float(f"{value:.13e}") for value in np.geomspace(
                float(xgrid[0]), x_upper, args.ej_energy_bins + 1
            )
        ])
        ej_y_edges = np.asarray(EJ_ANGULAR_GRIDS[args.ej_angular_grid], float)
        ej_nstates = (
            2 * args.ej_energy_bins * (ej_y_edges.size - 1)
        )

    edges = np.geomspace(r_inner, r_outer, args.radial_bins + 1)
    radii = np.sqrt(edges[:-1] * edges[1:])
    volumes = 4.0 * math.pi * (edges[1:]**3 - edges[:-1]**3) / 3.0
    rho_df = particle_mass * n0 * (
        df_forward_matrix(rh / radii, xgrid) @ gx
    )
    rng = np.random.default_rng(args.seed)
    rows = []
    jump_nbins = JUMP_J_OVER_JLC_EDGES.size - 1
    jump_tensor_total = np.zeros((jump_nbins, jump_nbins, jump_nbins))
    jump_capture_total = np.zeros((jump_nbins, jump_nbins))
    for radius, density, volume in zip(radii, rho_df, volumes):
        row = local_candidate_current(
            float(radius), float(density), rh, args.mbh,
            args.sigma0_over_m, args.w_kms,
            n0, particle_mass,
            xgrid, gx, args.pairs_per_bin, rng,
            cap_importance_fraction=args.cap_importance_fraction,
            direction_edge_fraction=args.direction_edge_fraction,
            direction_edge_power=args.direction_edge_power,
            survival_anomaly_order=args.survival_anomaly_order,
            survival_collision_order=args.survival_collision_order,
            ej_x_edges=ej_x_edges,
            ej_y_edges=ej_y_edges,
            ej_reservoir_radius_pc=r_outer,
            ej_energy_order=args.ej_energy_order,
            ej_include_nonkepler_exchange=(
                not args.ej_return_nonkepler_to_source
            ),
        )
        jump_tensor_total += row.pop(
            "_jump_balance_tensor_rate_density"
        ) * volume
        jump_capture_total += row.pop(
            "_jump_capture_matrix_rate_density"
        ) * volume
        ej_shell = row.pop("_ej_sparse_operator_rate_density", None)
        if ej_shell is not None:
            ej_shell = dict(ej_shell)
            for key in (
                "rate_direct_msun_per_myr",
                "rate_immediate_msun_per_myr",
                "rate_direct_binding_msun_kms2_per_myr",
                "rate_immediate_binding_msun_kms2_per_myr",
            ) + ej_energy_ledger.RATE_KEYS + (ej_energy_dg.RATE_KEYS if args.ej_energy_order else ()):
                ej_shell[key] = ej_shell[key] * volume
            ej_chunks.append(ej_shell)
        row["shell_volume_pc3"] = float(volume)
        row["shell_candidate_mdot_msun_myr"] = (
            row["candidate_rate_density_msun_per_pc3_per_myr"] * volume
        )
        row["shell_candidate_mdot_standard_error_msun_myr"] = (
            row["candidate_rate_density_standard_error_msun_per_pc3_per_myr"] * volume
        )
        row["shell_accepted_capture_mdot_msun_myr"] = (
            row["accepted_capture_rate_density_msun_per_pc3_per_myr"] * volume
        )
        row["shell_accepted_capture_mdot_standard_error_msun_myr"] = (
            row["accepted_capture_rate_density_standard_error_msun_per_pc3_per_myr"]
            * volume
        )
        row["shell_inward_candidate_mdot_msun_myr"] = (
            row["inward_candidate_rate_density_msun_per_pc3_per_myr"] * volume
        )
        row["shell_unresolved_kepler_domain_mdot_msun_myr"] = (
            row["unresolved_kepler_domain_rate_density_msun_per_pc3_per_myr"]
            * volume
        )
        row["shell_candidate_binding_energy_current_msun_kms2_myr"] = (
            row["candidate_binding_energy_rate_density_msun_kms2_per_pc3_per_myr"]
            * volume
        )
        row["shell_candidate_binding_energy_current_standard_error_msun_kms2_myr"] = (
            row[
                "candidate_binding_energy_rate_density_standard_error_msun_kms2_per_pc3_per_myr"
            ] * volume
        )
        row["shell_accepted_capture_energy_current_msun_kms2_myr"] = (
            row["accepted_binding_energy_rate_density_msun_kms2_per_pc3_per_myr"]
            * volume
        )
        row["shell_accepted_capture_energy_current_standard_error_msun_kms2_myr"] = (
            row["accepted_binding_energy_rate_density_standard_error_msun_kms2_per_pc3_per_myr"]
            * volume
        )
        row["shell_accepted_capture_source_j_over_jlc_histogram_msun_myr"] = [
            float(value * volume)
            for value in row[
                "accepted_capture_source_j_over_jlc_histogram_rate_density_msun_per_pc3_per_myr"
            ]
        ]
        row["shell_accepted_capture_source_j_over_jlc_histogram_standard_error_msun_myr"] = [
            float(value * volume)
            for value in row[
                "accepted_capture_source_j_over_jlc_histogram_rate_density_standard_error_msun_per_pc3_per_myr"
            ]
        ]
        row["shell_candidate_source_j_over_jlc_histogram_msun_myr"] = [
            float(value * volume)
            for value in row[
                "candidate_source_j_over_jlc_histogram_rate_density_msun_per_pc3_per_myr"
            ]
        ]
        row["shell_candidate_source_j_over_jlc_histogram_standard_error_msun_myr"] = [
            float(value * volume)
            for value in row[
                "candidate_source_j_over_jlc_histogram_rate_density_standard_error_msun_per_pc3_per_myr"
            ]
        ]
        row["shell_candidate_far_source_mdot_msun_myr"] = {
            key: float(value * volume)
            for key, value in row[
                "candidate_far_source_rate_density_msun_per_pc3_per_myr"
            ].items()
        }
        row["shell_candidate_far_source_mdot_standard_error_msun_myr"] = {
            key: float(value * volume)
            for key, value in row[
                "candidate_far_source_rate_density_standard_error_msun_per_pc3_per_myr"
            ].items()
        }
        row["shell_accepted_capture_far_source_mdot_msun_myr"] = {
            key: float(value * volume)
            for key, value in row[
                "accepted_capture_far_source_rate_density_msun_per_pc3_per_myr"
            ].items()
        }
        row["shell_accepted_capture_far_source_mdot_standard_error_msun_myr"] = {
            key: float(value * volume)
            for key, value in row[
                "accepted_capture_far_source_rate_density_standard_error_msun_per_pc3_per_myr"
            ].items()
        }
        row["rho_profile_msun_pc3"] = log_interp(
            r_profile, rho_profile, float(radius)
        )
        row["rho_df_over_profile"] = density / row["rho_profile_msun_pc3"]
        rows.append(row)

    mdot = float(sum(row["shell_candidate_mdot_msun_myr"] for row in rows))
    mdot_se = float(math.sqrt(sum(
        row["shell_candidate_mdot_standard_error_msun_myr"] ** 2
        for row in rows
    )))
    mdot_inward = float(sum(
        row["shell_inward_candidate_mdot_msun_myr"] for row in rows
    ))
    accepted_mdot = float(sum(
        row["shell_accepted_capture_mdot_msun_myr"] for row in rows
    ))
    accepted_mdot_se = float(math.sqrt(sum(
        row["shell_accepted_capture_mdot_standard_error_msun_myr"] ** 2
        for row in rows
    )))
    unresolved_mdot = float(sum(
        row["shell_unresolved_kepler_domain_mdot_msun_myr"] for row in rows
    ))
    candidate_energy_current = float(sum(
        row["shell_candidate_binding_energy_current_msun_kms2_myr"]
        for row in rows
    ))
    candidate_energy_current_se = float(math.sqrt(sum(
        row["shell_candidate_binding_energy_current_standard_error_msun_kms2_myr"] ** 2
        for row in rows
    )))
    accepted_energy_current = float(sum(
        row["shell_accepted_capture_energy_current_msun_kms2_myr"]
        for row in rows
    ))
    accepted_energy_current_se = float(math.sqrt(sum(
        row["shell_accepted_capture_energy_current_standard_error_msun_kms2_myr"] ** 2
        for row in rows
    )))
    candidate_source_j_hist = np.sum(np.array([
        row["shell_candidate_source_j_over_jlc_histogram_msun_myr"]
        for row in rows
    ]), axis=0)
    candidate_source_j_hist_se = np.sqrt(np.sum(np.array([
        row["shell_candidate_source_j_over_jlc_histogram_standard_error_msun_myr"]
        for row in rows
    ]) ** 2, axis=0))
    accepted_source_j_hist = np.sum(np.array([
        row["shell_accepted_capture_source_j_over_jlc_histogram_msun_myr"]
        for row in rows
    ]), axis=0)
    accepted_source_j_hist_se = np.sqrt(np.sum(np.array([
        row["shell_accepted_capture_source_j_over_jlc_histogram_standard_error_msun_myr"]
        for row in rows
    ]) ** 2, axis=0))
    candidate_far_source = {
        key: float(sum(
            row["shell_candidate_far_source_mdot_msun_myr"][key]
            for row in rows
        ))
        for key in (str(int(x)) for x in SOURCE_J_FLOOR_THRESHOLDS)
    }
    candidate_far_source_se = {
        key: float(math.sqrt(sum(
            row["shell_candidate_far_source_mdot_standard_error_msun_myr"][key] ** 2
            for row in rows
        )))
        for key in candidate_far_source
    }
    accepted_far_source = {
        key: float(sum(
            row["shell_accepted_capture_far_source_mdot_msun_myr"][key]
            for row in rows
        ))
        for key in (str(int(x)) for x in SOURCE_J_FLOOR_THRESHOLDS)
    }
    accepted_far_source_se = {
        key: float(math.sqrt(sum(
            row["shell_accepted_capture_far_source_mdot_standard_error_msun_myr"][key] ** 2
            for row in rows
        )))
        for key in accepted_far_source
    }
    rho_boundary = log_interp(r_profile, rho_profile, r_outer)
    sigma_boundary = log_interp(r_profile, sigma_profile, r_outer)
    mass_flux_scale = (
        4.0 * math.pi * r_outer**2 * rho_boundary * sigma_boundary
        * KMS_TO_PC_PER_MYR
    )
    energy_flux_scale = mass_flux_scale * sigma_boundary**2
    density_error = float(max(
        abs(row["rho_df_over_profile"] - 1.0) for row in rows
    ))
    candidate_relative_se = (
        mdot_se / mdot if mdot > 0.0 else math.inf
    )
    accepted_relative_se = (
        accepted_mdot_se / accepted_mdot if accepted_mdot > 0.0 else math.inf
    )
    unresolved_fraction = (
        unresolved_mdot / mdot if mdot > 0.0 else math.inf
    )
    jump_capture_isotropic = float(np.sum(jump_capture_total))
    jump_mass_identity_error = float(np.max(np.abs(
        np.sum(jump_tensor_total, axis=2) + jump_capture_total
    )))
    jump_scale = max(float(np.max(jump_capture_total)), accepted_mdot, 1.0e-300)
    jump_capture_relative_error = abs(
        jump_capture_isotropic / accepted_mdot - 1.0
    ) if accepted_mdot > 0.0 else math.inf
    numerical_pass = (
        True
        if (
            density_error <= 0.01
            and candidate_relative_se <= 0.10
            and accepted_relative_se <= 0.10
            and unresolved_fraction <= 0.01
            and jump_capture_relative_error <= 1.0e-10
            and jump_mass_identity_error / jump_scale <= 1.0e-10
        )
        else False
    )
    candidate_near_loss_cone_fraction = (
        float(np.sum(candidate_source_j_hist[:4]) / mdot)
        if mdot > 0.0 else math.inf
    )
    accepted_near_loss_cone_fraction = (
        float(np.sum(accepted_source_j_hist[:4]) / accepted_mdot)
        if accepted_mdot > 0.0 else math.inf
    )
    status = (
        "ISOTROPIC_INJECTION_CEILING_CONVERGED"
        if numerical_pass else "DIAGNOSTIC_ONLY"
    )
    ej_operator_summary = None
    if args.out_ej_operator is not None:
        if (
            not ej_chunks or ej_x_edges is None or ej_y_edges is None
            or ej_nstates is None
        ):
            raise RuntimeError("the requested EJ operator has no transition records")
        keys = (
            "source1", "source2", "pre", "post",
            "rate_direct_msun_per_myr", "rate_immediate_msun_per_myr",
            "sample_count",
            "rate_direct_binding_msun_kms2_per_myr",
            "rate_immediate_binding_msun_kms2_per_myr",
        )
        combined = aggregate_ej_jump_records(
            *(np.concatenate([chunk[key] for chunk in ej_chunks]) for key in keys),
            nstates=ej_nstates,
        )
        combined.update(ej_energy_ledger.combine(ej_chunks, ej_nstates))
        if args.ej_energy_order:
            combined.update(ej_energy_dg.aggregate(ej_chunks, ej_nstates))
        capture_mask = combined["post"] == EJ_POST_CAPTURE
        reservoir_mask = combined["post"] == EJ_POST_RESERVOIR
        direct_roundtrip = float(np.sum(
            combined["rate_direct_msun_per_myr"][capture_mask]
        ))
        immediate_roundtrip = float(np.sum(
            combined["rate_immediate_msun_per_myr"][capture_mask]
        ))
        direct_energy_roundtrip = float(np.sum(
            combined["rate_direct_binding_msun_kms2_per_myr"][capture_mask]
        ))
        immediate_energy_roundtrip = float(np.sum(
            combined["rate_immediate_binding_msun_kms2_per_myr"][capture_mask]
        ))
        direct_error = abs(direct_roundtrip / accepted_mdot - 1.0)
        immediate_error = abs(immediate_roundtrip / mdot - 1.0)
        direct_energy_error = abs(
            direct_energy_roundtrip / accepted_energy_current - 1.0
        )
        immediate_energy_error = abs(
            immediate_energy_roundtrip / candidate_energy_current - 1.0
        )
        if max(
            direct_error, immediate_error,
            direct_energy_error, immediate_energy_error,
        ) > 1.0e-10:
            raise RuntimeError(
                "the EJ operator does not reproduce the isotropic mass and "
                "capture-energy currents"
            )
        metadata = {
            "schema": "finite-angle-ej-jump-operator-v1",
            "energy_ledger_schema": ej_energy_ledger.SCHEMA,
            "energy_occupation_order": args.ej_energy_order,
            "energy_dg_schema": ej_energy_dg.SCHEMA if args.ej_energy_order else None,
            "energy_ledger_includes_self_cell_collisions": True,
            "profile_sha256": profile_sha256,
            "df_table_sha256": df_table_sha256,
            "mbh_msun": args.mbh,
            "rh_pc": rh,
            "r_inner_pc": r_inner,
            "reservoir_radius_pc": r_outer,
            "sigma0_over_m_cm2_g": args.sigma0_over_m,
            "w_kms": args.w_kms,
            "radial_bins": args.radial_bins,
            "pairs_per_bin": args.pairs_per_bin,
            "seed": args.seed,
            "energy_bins": args.ej_energy_bins,
            "angular_grid": args.ej_angular_grid,
            "angular_bins": int(ej_y_edges.size - 1),
            "domain_states_per_cell": 2,
            "state_count": ej_nstates,
            "post_capture_sentinel": EJ_POST_CAPTURE,
            "post_reservoir_sentinel": EJ_POST_RESERVOIR,
            "reservoir_state_rule": "exact Kepler apocentre >= reporting radius",
            "nonkepler_exchange": (
                "outer-fluid reservoir"
                if not args.ej_return_nonkepler_to_source
                else "returned to source, regression only"
            ),
            "direct_capture_roundtrip_msun_per_myr": direct_roundtrip,
            "immediate_capture_roundtrip_msun_per_myr": immediate_roundtrip,
            "direct_capture_roundtrip_relative_error": direct_error,
            "immediate_capture_roundtrip_relative_error": immediate_error,
            "direct_capture_binding_roundtrip_msun_kms2_per_myr": (
                direct_energy_roundtrip
            ),
            "immediate_capture_binding_roundtrip_msun_kms2_per_myr": (
                immediate_energy_roundtrip
            ),
            "direct_capture_binding_roundtrip_relative_error": (
                direct_energy_error
            ),
            "immediate_capture_binding_roundtrip_relative_error": (
                immediate_energy_error
            ),
            "record_count": int(combined["source1"].size),
            "represented_event_records": int(np.sum(combined["sample_count"])),
            "reservoir_exchange_rate_isotropic_msun_per_myr": float(np.sum(
                combined["rate_direct_msun_per_myr"][reservoir_mask]
            )),
            "note": (
                "Each rate is multiplied by D[source1] D[source2]. Internal "
                "states have even indices and reservoir-connected states have "
                "odd indices. The direct and immediate capture columns differ "
                "only for post-collision loss-cone states."
            ),
        }
        args.out_ej_operator.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            args.out_ej_operator,
            **combined,
            x_edges=np.asarray(ej_x_edges),
            j_over_jlc_edges=np.asarray(ej_y_edges),
            metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)),
        )
        ej_operator_summary = dict(metadata)
        ej_operator_summary.update({
            "path": args.out_ej_operator.name,
            "sha256": sha256_file(args.out_ej_operator),
        })
    result = {
        "schema": "finite-angle-loss-cone-capture-v6",
        "status": status,
        "profile": args.profile.name,
        "profile_sha256": profile_sha256,
        "bridge_json_sha256": norm["bridge_json_sha256"],
        "df_table": args.df_table.name,
        "df_table_sha256": df_table_sha256,
        "df_table_schema": table_metadata.get("schema"),
        "mbh_msun": args.mbh,
        "rh_pc": rh,
        "sigma0_over_m_cm2_g": args.sigma0_over_m,
        "w_kms": args.w_kms,
        "r_inner_pc": r_inner,
        "r_inner_over_rg": r_inner / rg,
        "radial_inner_limit_model": (
            "max(8 Rg outside-loss-cone kinematic limit, profile support, "
            "DF fit support)"
        ),
        "represented_df_inner_radius_pc": represented_inner_radius,
        "r_outer_pc": r_outer,
        "radial_bins": args.radial_bins,
        "pairs_per_bin": args.pairs_per_bin,
        "cap_importance_fraction": args.cap_importance_fraction,
        "direction_edge_fraction": args.direction_edge_fraction,
        "direction_edge_power": args.direction_edge_power,
        "survival_anomaly_order": args.survival_anomaly_order,
        "survival_collision_quadrature_order": args.survival_collision_order,
        "seed": args.seed,
        "ej_jump_operator": ej_operator_summary,
        "candidate_mdot_msun_per_myr": mdot,
        "candidate_mdot_standard_error_msun_per_myr": mdot_se,
        "candidate_relative_standard_error": candidate_relative_se,
        "inward_candidate_mdot_msun_per_myr": mdot_inward,
        "accepted_capture_mdot_msun_per_myr": accepted_mdot,
        "accepted_capture_mdot_standard_error_msun_per_myr": accepted_mdot_se,
        "accepted_capture_relative_standard_error": accepted_relative_se,
        "direct_no_rescatter_mdot_msun_per_myr": accepted_mdot,
        "direct_no_rescatter_mdot_standard_error_msun_per_myr": accepted_mdot_se,
        "unresolved_kepler_domain_mdot_upper_bound_msun_per_myr": unresolved_mdot,
        "unresolved_kepler_domain_fraction_of_candidate": unresolved_fraction,
        "loss_cone_depletion_status": "NONLOCAL_STEADY_STATE_REQUIRED",
        "candidate_current_fraction_from_source_j_below_5_jlc": (
            candidate_near_loss_cone_fraction
        ),
        "accepted_current_fraction_from_source_j_below_5_jlc": (
            accepted_near_loss_cone_fraction
        ),
        "candidate_binding_energy_current_msun_kms2_per_myr": (
            candidate_energy_current
        ),
        "candidate_binding_energy_current_standard_error_msun_kms2_per_myr": (
            candidate_energy_current_se
        ),
        "candidate_mean_binding_energy_kms2": (
            candidate_energy_current / mdot if mdot > 0.0 else None
        ),
        "accepted_capture_binding_energy_current_msun_kms2_per_myr": (
            accepted_energy_current
        ),
        "accepted_capture_binding_energy_current_standard_error_msun_kms2_per_myr": (
            accepted_energy_current_se
        ),
        "accepted_capture_mean_binding_energy_kms2": (
            accepted_energy_current / accepted_mdot
            if accepted_mdot > 0.0 else None
        ),
        "boundary_density_msun_pc3": rho_boundary,
        "boundary_sigma1d_kms": sigma_boundary,
        "mass_flux_scale_msun_per_myr": mass_flux_scale,
        "energy_flux_scale_msun_kms2_per_myr": energy_flux_scale,
        "C_M_candidate_ceiling": mdot / mass_flux_scale,
        "C_M_candidate_ceiling_standard_error": mdot_se / mass_flux_scale,
        "C_M_direct_no_rescatter": accepted_mdot / mass_flux_scale,
        "C_M_direct_no_rescatter_standard_error": (
            accepted_mdot_se / mass_flux_scale
        ),
        "C_E_candidate_binding_ceiling": (
            candidate_energy_current / energy_flux_scale
        ),
        "C_E_candidate_binding_ceiling_standard_error": (
            candidate_energy_current_se / energy_flux_scale
        ),
        "C_E_direct_no_rescatter_binding": (
            accepted_energy_current / energy_flux_scale
        ),
        "C_E_direct_no_rescatter_binding_standard_error": (
            accepted_energy_current_se / energy_flux_scale
        ),
        "thermal_sink_current_msun_kms2_per_myr": None,
        "C_E_thermal_sink": None,
        "returned_capture_energy_source_ceiling_msun_kms2_per_myr": None,
        "thermal_closure_status": "NOT_MEASURED",
        "thermal_closure_note": "Capture orbital energy is not a thermal boundary current. No source or sink is assigned from the boundary dispersion.",
        "candidate_source_j_over_jlc_histogram": {
            "edges": [1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0, None],
            "mdot_msun_per_myr": candidate_source_j_hist.tolist(),
            "standard_error_msun_per_myr": candidate_source_j_hist_se.tolist(),
            "fraction": (
                candidate_source_j_hist / mdot
            ).tolist() if mdot > 0.0 else None,
            "note": (
                "Pre-collision angular momentum of the particle injected into "
                "the loss cone, before applying the no-rescattering factor."
            ),
        },
        "accepted_capture_source_j_over_jlc_histogram": {
            "edges": [1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0, None],
            "mdot_msun_per_myr": accepted_source_j_hist.tolist(),
            "standard_error_msun_per_myr": accepted_source_j_hist_se.tolist(),
            "fraction": (
                accepted_source_j_hist / accepted_mdot
            ).tolist() if accepted_mdot > 0.0 else None,
            "note": (
                "Pre-collision angular momentum of the particle that enters "
                "the loss cone. Concentration near unity diagnoses sensitivity "
                "to a depleted non-local boundary layer in J."
            ),
        },
        "candidate_far_source_currents": {
            key: {
                "minimum_source_j_over_jlc": float(key),
                "mdot_msun_per_myr": candidate_far_source[key],
                "standard_error_msun_per_myr": candidate_far_source_se[key],
                "fraction_of_isotropic_candidate_ceiling": (
                    candidate_far_source[key] / mdot if mdot > 0.0 else None
                ),
            }
            for key in candidate_far_source
        },
        "accepted_capture_far_source_currents": {
            key: {
                "minimum_source_j_over_jlc": float(key),
                "mdot_msun_per_myr": accepted_far_source[key],
                "standard_error_msun_per_myr": accepted_far_source_se[key],
                "fraction_of_isotropic_direct_component": (
                    accepted_far_source[key] / accepted_mdot
                    if accepted_mdot > 0.0 else None
                ),
            }
            for key in accepted_far_source
        },
        "nonlocal_jump_balance": {
            "j_over_jlc_edges": [
                float(x) if np.isfinite(x) else None
                for x in JUMP_J_OVER_JLC_EDGES
            ],
            "population_rate_tensor_msun_per_myr": jump_tensor_total.tolist(),
            "capture_rate_matrix_msun_per_myr": jump_capture_total.tolist(),
            "isotropic_direct_no_rescatter_mdot_msun_per_myr": (
                jump_capture_isotropic
            ),
            "capture_roundtrip_relative_error": jump_capture_relative_error,
            "mass_balance_max_absolute_error_msun_per_myr": jump_mass_identity_error,
            "mass_balance_max_relative_error": jump_mass_identity_error / jump_scale,
            "definition": (
                "T[i,j,k] is the net population rate in output bin k from "
                "collisions whose two source particles occupy bins i and j. "
                "For depletion factors D, dot N_k=sum_ij D_i D_j T_ijk."
            ),
        },
        "rho_df_profile_max_abs_relative_error": density_error,
        "raw_candidate_particles": int(sum(
            row["raw_candidate_particles"] for row in rows
        )),
        "note": (
            "The candidate current counts finite-angle transitions from outside "
            "to bound J<J_lc states and is the isotropic-reservoir injection "
            "ceiling. The direct no-rescattering current multiplies each "
            "transition by the probability of reaching the plunge surface before "
            "another collision. Neither quantity is an absolute steady capture "
            "rate because the non-local depleted (E,J) distribution has not been "
            "solved."
        ),
        "rows": rows,
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    with args.out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
