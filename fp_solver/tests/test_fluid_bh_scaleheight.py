#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fluid_bh_scaleheight import (  # noqa: E402
    bh_corrected_lmfp_inverse,
    bh_lmfp_scaleheight_factor,
    conservative_specific_energy_rate,
    harmonic_effective_conductivity,
    shell_luminosity,
    unmodeled_mass_loss_stop_reason,
)


def main() -> int:
    r = np.geomspace(1.0e-4, 1.0e4, 400)
    rho = 3.0 * np.ones_like(r)
    factor = bh_lmfp_scaleheight_factor(r, rho, 2.0)
    assert np.all((factor > 0.0) & (factor <= 1.0))
    assert np.all(np.diff(factor) > 0.0)
    assert factor[0] < 1.0e-10
    assert abs(factor[-1] - 1.0) < 1.0e-10
    assert np.allclose(bh_lmfp_scaleheight_factor(r, rho, 0.0), 1.0)
    corrected, factor2 = bh_corrected_lmfp_inverse(
        np.ones_like(r), r, rho, 2.0
    )
    assert np.allclose(factor2, factor)
    assert np.allclose(corrected * factor, 1.0)
    stock = harmonic_effective_conductivity(
        np.full_like(r, 1.0e-6), np.ones_like(r)
    )
    with_bh = harmonic_effective_conductivity(
        np.full_like(r, 1.0e-6), corrected
    )
    assert np.all(with_bh <= stock)

    grid = np.geomspace(0.1, 10.0, 32)
    conductivity = np.linspace(0.5, 1.5, grid.size)
    assert np.all(shell_luminosity(
        grid, np.ones_like(grid), conductivity
    ) == 0.0)
    luminosity = shell_luminosity(
        grid, np.log1p(grid), conductivity
    )
    assert np.all(luminosity[:-1] < 0.0) and luminosity[-1] == 0.0
    shell_mass = np.linspace(1.0, 2.0, grid.size)
    enclosed = np.cumsum(shell_mass)
    inner = -0.37
    rate = conservative_specific_energy_rate(enclosed, luminosity, inner)
    assert abs(float(np.dot(shell_mass, rate)) - inner) < 2.0e-14
    dt = 0.013
    delta_u = dt * rate
    actual = float(np.dot(shell_mass, delta_u))
    expected = inner * dt
    scale = float(np.dot(shell_mass, np.abs(delta_u))) + abs(expected)
    assert abs(actual - expected) / scale < 2.0e-14
    assert unmodeled_mass_loss_stop_reason(10.0, 1000.0, 100.0, 0.02, 0.2) is None
    assert (
        unmodeled_mass_loss_stop_reason(21.0, 1000.0, 1000.0, 0.02, 0.2)
        == "unmodeled_fluid_mass_loss_limit"
    )
    assert (
        unmodeled_mass_loss_stop_reason(21.0, 10000.0, 100.0, 0.02, 0.2)
        == "unmodeled_gnc_mass_loss_limit"
    )
    print("PASS: BH-limited LMFP scale height")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
