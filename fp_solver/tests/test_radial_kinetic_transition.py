"""Invariant and endpoint tests for the non-orbit-averaged transition solver."""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bondi_hydrodynamic import critical_bondi_profile  # noqa: E402
from radial_kinetic_transition import (  # noqa: E402
    ballistic_capture_rate,
    collision_step,
    conservative_shift,
    discrete_maxwellian_from_moments,
    macroscopic_moments,
    make_grid,
    outer_reservoir,
)


def main():
    for gamma in (1.4, 5.0 / 3.0):
        profile = critical_bondi_profile(np.geomspace(0.05, 8.0, 48), gamma)
        assert profile["maximum_mass_flux_relative_residual"] < 1.0e-12
        assert profile["maximum_bernoulli_absolute_residual"] < 2.0e-10
        assert np.all(profile["radial_velocity"] < 0.0)
    assert critical_bondi_profile(np.geomspace(0.03, 8.0, 48), 1.4)["mach"][0] > 1.0

    grid = make_grid(0.1, 4.0, 12, 24, 16, 9.0)
    rho = 1.0 + 0.2 / grid.r
    velocity = -0.15 / np.sqrt(grid.r)
    temperature = 0.7 + 0.1 / grid.r
    m0 = grid.r**2 * rho / (2.0 * math.pi)
    m1 = m0 * velocity
    m2 = m0 * (velocity**2 + 3.0 * temperature)
    h, diagnostic = discrete_maxwellian_from_moments(grid, m0, m1, m2)
    moments = macroscopic_moments(h, grid)
    assert diagnostic["maximum_discrete_moment_relative_error"] < 2.0e-9
    assert np.max(np.abs(moments["radial_velocity"] - velocity)) < 2.0e-9
    assert np.max(np.abs(moments["temperature"] / temperature - 1.0)) < 2.0e-9
    # The odd radial moment is explicit and non-zero.
    assert np.all(moments["radial_velocity"] < -0.01)

    collided, collision = collision_step(h * (1.0 + 0.03 * np.sin(
        grid.vr[None, :, None]
    )), grid, 0.1, 0.5, 1.0 / (5.0 / 3.0))
    assert np.all(collided >= 0.0)
    assert collision["maximum_collision_invariant_relative_error"] < 2.0e-9

    q = np.ones((8, 3))
    edges = np.linspace(0.0, 1.0, 9)
    shifted, ledger = conservative_shift(
        q, edges, np.array([0.03, -0.04, 0.0]), np.zeros(3), np.zeros(3)
    )
    before = np.sum(q * np.diff(edges)[:, None], axis=0)
    after = np.sum(shifted * np.diff(edges)[:, None], axis=0)
    expected = before - ledger["left_out"] - ledger["right_out"]
    assert np.max(np.abs(after - expected)) < 1.0e-13

    reservoir, _ = outer_reservoir(grid, 5.0 / 3.0)
    ballistic = ballistic_capture_rate(grid, reservoir)
    assert math.isfinite(ballistic) and ballistic > 0.0
    print("PASS: radial kinetic invariants and Bondi endpoint")


if __name__ == "__main__":
    main()
