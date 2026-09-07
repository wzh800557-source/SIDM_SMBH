#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hydrostatic_bridge import (  # noqa: E402
    C_KMS,
    G_PC_KMS2_MSUN,
    build_mass_matched_adiabatic_bridge,
    build_one,
    extrapolate_inner_edge,
    provenance_name,
)


def main() -> int:
    assert provenance_name(Path("/one/root/profile.txt")) == "profile.txt"
    assert provenance_name(Path("/another/root/profile.txt")) == "profile.txt"
    r0 = 25.7
    r = np.geomspace(r0, 300.0, 128)
    rho = 2.9 * (r / r0) ** -0.12
    sig = 31.3 * (r / r0) ** -0.02
    rho_edge, pressure_edge = extrapolate_inner_edge(r, rho, sig)
    mbh = 4.0e6
    rg = G_PC_KMS2_MSUN * mbh / C_KMS**2
    rinner = np.geomspace(6.0 * rg, r0, 600)
    result = build_mass_matched_adiabatic_bridge(
        rinner, r0, rho[0], rho_edge, pressure_edge, mbh
    )
    rhoi, pressure, sigma, rho_edge_inside, iterations, residual, mass_error = result
    assert np.all(np.isfinite(rhoi)) and np.all(rhoi > 0.0)
    assert np.all(np.isfinite(sigma)) and np.all(sigma > 0.0)
    assert abs(mass_error) < 1.0e-8
    assert abs(pressure[-1] / pressure_edge - 1.0) < 1.0e-13
    assert abs(rhoi[-1] / rho_edge_inside - 1.0) < 1.0e-13
    assert iterations <= 300 and residual < 2.0e-9
    bridged = build_one(
        r,
        rho,
        sig,
        mbh,
        100.0,
        0.0,
        G_PC_KMS2_MSUN * mbh / sig[0] ** 2,
        0.12,
        0.75,
        "adiabatic-mass-matched",
        600,
        "yukawa-tchannel",
        80.0,
    )
    assert 0.0 < bridged["mass_inside_r_in_fraction_of_first_shell"] < 0.02
    assert 0.0 < bridged["volume_inside_r_in_fraction_of_first_shell"] < 0.02
    print("PASS: mass-matched adiabatic BH bridge")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
