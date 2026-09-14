#!/usr/bin/env python3
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from finite_angle_capture import (  # noqa: E402
    EJ_ANGULAR_GRIDS,
    angular_mean_v_sigma_total_ratio,
    circle_intersection_area,
    ej_jump_records,
    ej_state_indices,
    jump_balance_tensor,
    orbit_survival_probability,
    orbit_survival_probabilities,
    sample_outside_radial_cosine,
    sample_outside_loss_cone,
    sample_tchannel_u,
    scatter_equal_mass,
    scatter_equal_mass_disc_importance,
    scatter_equal_mass_importance,
)


def main() -> int:
    rng = np.random.default_rng(8472)
    assert set(EJ_ANGULAR_GRIDS) == {
        "coarse", "default", "fine", "ultra", "hyper", "j257"
    }
    assert (
        len(EJ_ANGULAR_GRIDS["coarse"])
        < len(EJ_ANGULAR_GRIDS["default"])
        < len(EJ_ANGULAR_GRIDS["fine"])
        < len(EJ_ANGULAR_GRIDS["ultra"])
        < len(EJ_ANGULAR_GRIDS["hyper"])
        < len(EJ_ANGULAR_GRIDS["j257"])
    )
    assert len(EJ_ANGULAR_GRIDS["hyper"]) - 1 == 129
    assert len(EJ_ANGULAR_GRIDS["j257"]) - 1 == 257
    for edges in EJ_ANGULAR_GRIDS.values():
        assert edges[0] == 1.0 and np.isinf(edges[-1])
        assert np.all(np.diff(edges) > 0.0)
    a = np.full(300000, 9.0)
    u = sample_tchannel_u(a, rng)
    # Exact conditional mean from numerical quadrature.
    grid = np.linspace(0.0, 1.0, 200001)
    pdf = 10.0 / (1.0 + 9.0 * grid) ** 2
    mean_exact = np.trapezoid(grid * pdf, grid)
    assert abs(float(np.mean(u)) / mean_exact - 1.0) < 8.0e-3

    v1 = rng.normal(size=(10000, 3)) * 300.0
    v2 = rng.normal(size=(10000, 3)) * 300.0
    p0 = v1 + v2
    e0 = np.sum(v1 * v1 + v2 * v2, axis=1)
    out1, out2, total_ratio = scatter_equal_mass(v1, v2, 80.0, rng)
    assert np.max(np.abs(out1 + out2 - p0)) < 2.0e-12
    assert np.max(np.abs(np.sum(out1 * out1 + out2 * out2, axis=1) - e0)) < 2.0e-9
    assert np.all((total_ratio > 0.0) & (total_ratio <= 1.0))

    area = circle_intersection_area(
        np.array([0.0, 3.0, 0.0]), np.array([0.25, 0.5, 2.0])
    )
    assert np.allclose(area, [math.pi * 0.25**2, 0.0, math.pi])

    # The cap-biased scattering proposal must integrate back to the exact
    # t-channel law.  A repeated fixed pair isolates the angular estimator.
    nimportance = 180000
    fixed1 = np.repeat([[100.0, 0.0, 50.0]], nimportance, axis=0)
    fixed2 = np.repeat([[-100.0, 0.0, -20.0]], nimportance, axis=0)
    gm = 4.30091e-3 * 4.0e6
    jlc = 4.0 * gm / 2.99792458e5
    imp1, imp2, _, angular_weight = scatter_equal_mass_importance(
        fixed1, fixed2, 80.0, 0.01, jlc, rng, cap_fraction=0.9
    )
    pre_direction = (fixed1[0] - fixed2[0])
    pre_direction /= np.linalg.norm(pre_direction)
    post_direction = imp1 - imp2
    post_direction /= np.linalg.norm(post_direction, axis=1)[:, None]
    sampled_u = 0.5 * (1.0 - post_direction @ pre_direction)
    afixed = (np.linalg.norm(fixed1[0] - fixed2[0]) / 80.0) ** 2
    grid_fixed = np.linspace(0.0, 1.0, 300001)
    pdf_fixed = (1.0 + afixed) / (1.0 + afixed * grid_fixed) ** 2
    exact_u = np.trapezoid(grid_fixed * pdf_fixed, grid_fixed)
    assert abs(float(np.mean(angular_weight)) - 1.0) < 1.2e-2
    assert abs(float(np.mean(angular_weight * sampled_u)) / exact_u - 1.0) < 1.5e-2
    direct1, direct2, _ = scatter_equal_mass(fixed1, fixed2, 80.0, rng)
    direct_count = (
        (0.01 * np.linalg.norm(direct1[:, :2], axis=1) < jlc).astype(float)
        + (0.01 * np.linalg.norm(direct2[:, :2], axis=1) < jlc).astype(float)
    )
    importance_count = angular_weight * (
        (0.01 * np.linalg.norm(imp1[:, :2], axis=1) < jlc).astype(float)
        + (0.01 * np.linalg.norm(imp2[:, :2], axis=1) < jlc).astype(float)
    )
    direct_mean = float(np.mean(direct_count))
    importance_mean = float(np.mean(importance_count))
    combined_se = math.sqrt(
        float(np.var(direct_count, ddof=1) + np.var(importance_count, ddof=1))
        / nimportance
    )
    assert abs(importance_mean - direct_mean) < 5.0 * combined_se

    disc1, disc2, _, disc_weight = scatter_equal_mass_disc_importance(
        fixed1, fixed2, 80.0, 0.01, jlc, rng, cap_fraction=0.7
    )
    disc_count = disc_weight * (
        (0.01 * np.linalg.norm(disc1[:, :2], axis=1) < jlc).astype(float)
        + (0.01 * np.linalg.norm(disc2[:, :2], axis=1) < jlc).astype(float)
    )
    disc_mean = float(np.mean(disc_count))
    disc_se = float(np.std(disc_count, ddof=1) / math.sqrt(nimportance))
    direct_se = float(np.std(direct_count, ddof=1) / math.sqrt(nimportance))
    assert abs(disc_mean - direct_mean) < 5.0 * math.hypot(disc_se, direct_se)

    cmax = np.full(300000, 0.8)
    radial_cosine, direction_weight = sample_outside_radial_cosine(
        cmax, rng, edge_fraction=0.8, edge_power=5.0
    )
    t = 1.0 - np.abs(radial_cosine) / cmax
    assert abs(float(np.mean(direction_weight)) - 1.0) < 8.0e-3
    assert abs(float(np.mean(direction_weight * t)) - 0.5) < 8.0e-3
    local_x = np.geomspace(0.03, 1.0e5, 64)
    local_g = np.ones_like(local_x)
    local_v, _, _ = sample_outside_loss_cone(
        40000,
        0.003,
        17.5,
        np.sqrt(gm / 17.5),
        jlc,
        local_x,
        local_g,
        rng,
        direction_edge_fraction=0.8,
        direction_edge_power=5.0,
    )
    assert np.min(0.003 * np.linalg.norm(local_v[:, :2], axis=1) / jlc) >= 1.0 - 1.0e-10

    # The outside-density fraction is normalized by the complete local DF,
    # including high-binding-energy states whose speed is too small to remain
    # outside the loss cone.  For g=constant, psi=1, and J_lc/r=1, the exact
    # fraction is 1/(2 sqrt(2)).
    fraction_grid = np.geomspace(1.0e-8, 2.0, 128)
    _, outside_fraction, _ = sample_outside_loss_cone(
        2000,
        1.0,
        1.0,
        1.0,
        1.0,
        fraction_grid,
        np.ones_like(fraction_grid),
        rng,
    )
    assert abs(outside_fraction - 1.0 / (2.0 * math.sqrt(2.0))) < 3.0e-4

    tensor, capture = jump_balance_tensor(
        np.array([1.2, 4.0]),
        np.array([3.0, 8.0]),
        np.array([0.5, 5.0]),
        np.array([4.0, 9.0]),
        np.array([True, True]),
        np.array([True, True]),
        np.array([True, False]),
        np.array([False, False]),
        np.array([0.75, 0.0]),
        np.array([0.0, 0.0]),
        np.array([2.0, 3.0]),
    )
    assert np.max(np.abs(np.sum(tensor, axis=2) + capture)) < 1.0e-14
    assert abs(float(np.sum(capture)) - 1.5) < 1.0e-14

    # The two-dimensional jump operator separates internal orbits from states
    # whose exact apocentre connects to the fixed outer reservoir. It must also
    # retain both the direct no-rescattering and immediate-capture weights.
    ej_x_edges = np.array([0.01, 0.1, 1.0])
    source_states = ej_state_indices(
        np.array([0.5, 0.02]),
        np.array([0.2, 0.3]),
        1.0,
        1.0,
        0.1,
        5.0,
        ej_x_edges,
    )
    assert source_states[0] % 2 == 0
    assert source_states[1] % 2 == 1
    ej = ej_jump_records(
        np.array([0.5]),
        np.array([0.02]),
        np.array([0.2]),
        np.array([0.3]),
        np.array([0.5]),
        np.array([0.4]),
        np.array([0.05]),
        np.array([0.3]),
        np.array([True]),
        np.array([False]),
        np.array([0.75]),
        np.array([0.0]),
        np.array([2.0]),
        1.0,
        1.0,
        0.1,
        5.0,
        ej_x_edges,
    )
    capture_record = ej["post"] == -1
    assert np.count_nonzero(capture_record) == 1
    assert abs(float(np.sum(ej["rate_direct_msun_per_myr"][capture_record])) - 1.5) < 1.0e-14
    assert abs(float(np.sum(ej["rate_immediate_msun_per_myr"][capture_record])) - 2.0) < 1.0e-14
    assert np.any(ej["post"] >= 0)

    mu = np.linspace(-1.0, 1.0, 500001)
    test_speed = 180.0
    background_speed = 95.0
    relative = np.sqrt(
        test_speed**2 + background_speed**2
        - 2.0 * test_speed * background_speed * mu
    )
    brute = np.trapezoid(relative / (1.0 + (relative / 80.0) ** 2), mu) / 2.0
    exact = angular_mean_v_sigma_total_ratio(
        test_speed, background_speed, 80.0
    )
    assert abs(exact / brute - 1.0) < 2.0e-10

    # A zero cross section gives unit survival for both directions.  The
    # outward orbit exercises the apocentre and return-leg integration.
    mbh = 4.0e6
    rh = 20.0
    gm = 4.30091e-3 * mbh
    radius = 0.01
    speed = np.sqrt(gm / radius)
    jlc = 4.0 * gm / 2.99792458e5
    vt = 0.5 * jlc / radius
    vr = np.sqrt(speed**2 - vt**2)
    xgrid = np.geomspace(0.03, 1.0e6, 80)
    gx = np.ones_like(xgrid)
    inward = np.array([vt, 0.0, -vr])
    outward = np.array([vt, 0.0, vr])
    pin, tauin, _, _ = orbit_survival_probability(
        radius, inward, rh, mbh, 1.0, 1.0, 0.0, 80.0, xgrid, gx
    )
    pout, tauout, _, _ = orbit_survival_probability(
        radius, outward, rh, mbh, 1.0, 1.0, 0.0, 80.0, xgrid, gx
    )
    assert pin == 1.0 and pout == 1.0 and tauin == 0.0 and tauout == 0.0
    batch_p, batch_tau, _, _ = orbit_survival_probabilities(
        radius,
        np.vstack([inward, outward]),
        rh,
        mbh,
        3.0,
        1.0,
        100.0,
        80.0,
        xgrid,
        gx,
        anomaly_order=24,
        collision_quadrature_order=96,
    )
    for index, velocity in enumerate((inward, outward)):
        scalar_p, scalar_tau, _, _ = orbit_survival_probability(
            radius, velocity, rh, mbh, 3.0, 1.0, 100.0, 80.0, xgrid, gx
        )
        assert abs(batch_p[index] - scalar_p) < 2.0e-12
        assert abs(batch_tau[index] - scalar_tau) < 2.0e-12
    fast_p, fast_tau, _, _ = orbit_survival_probabilities(
        radius,
        np.vstack([inward, outward]),
        rh,
        mbh,
        3.0,
        1.0,
        100.0,
        80.0,
        xgrid,
        gx,
    )
    assert np.max(np.abs(fast_p / batch_p - 1.0)) < 1.0e-3
    assert np.max(np.abs(fast_tau / batch_tau - 1.0)) < 2.0e-3

    # A jump occurring inside the orbit-dependent plunge surface is captured
    # immediately; it should not be rejected by the orbit integral.
    rg = gm / 2.99792458e5**2
    immediate_radius = 12.0 * rg
    immediate_speed = np.sqrt(gm / immediate_radius)
    immediate_velocity = np.array([immediate_speed, 0.0, 0.0])
    immediate_p, immediate_tau, _, immediate_capture_radius = (
        orbit_survival_probability(
            immediate_radius,
            immediate_velocity,
            rh,
            mbh,
            3.0,
            1.0,
            100.0,
            80.0,
            xgrid,
            gx,
        )
    )
    assert immediate_capture_radius > immediate_radius
    assert immediate_p == 1.0 and immediate_tau == 0.0
    print("PASS: finite-angle t-channel kinematics")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
