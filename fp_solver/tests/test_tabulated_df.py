#!/usr/bin/env python3
"""Regression test for the non-power-law physical DF inversion."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from normalize_gnc_df import (  # noqa: E402
    construct_tabulated_isotropic_df,
    df_forward_matrix,
    table_reconstructed_density,
)


def main() -> int:
    rh = 10.0
    n0 = 2.0
    particle_mass = 1.0
    xgrid = np.geomspace(0.03, 1.0e4, 64)

    # A positive DF that changes continuously from a weak-binding power law to
    # an approximately flat high-binding distribution.  Its density is not a
    # single radial power law.
    gx_true = 0.2 * xgrid**-0.75 * (
        1.0 + (xgrid / 2.0) ** 4
    ) ** (0.75 / 4.0)
    radius = np.geomspace(2.0e-3, 50.0, 900)
    rho = particle_mass * n0 * (
        df_forward_matrix(rh / radius, xgrid) @ gx_true
    )
    sigma = np.sqrt(1.0 / radius)
    good = rho > 0.0
    radius, rho, sigma = radius[good], rho[good], sigma[good]

    gx_fit, diag = construct_tabulated_isotropic_df(
        radius,
        rho,
        rh,
        n0,
        particle_mass,
        xgrid,
        normalization_radius=1.0,
        fit_outer_radius=5.0,
        rg=1.0e-8,
        regularization=1.0e-2,
        fit_points=256,
    )
    evaluate = np.geomspace(0.02, 1.0, 160)
    target = np.exp(np.interp(
        np.log(evaluate), np.log(radius), np.log(rho)
    ))
    recovered = table_reconstructed_density(
        evaluate, rh, n0, xgrid, gx_fit, particle_mass
    )
    relative = recovered / target - 1.0
    assert diag["df_shape_model"] == "regularized non-negative Abel inversion"
    assert np.all(np.isfinite(gx_fit)) and np.all(gx_fit > 0.0)
    assert np.max(np.abs(relative)) < 3.0e-3
    assert abs(relative[-1]) < 5.0e-10
    print("PASS: non-power-law physical DF inversion")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
