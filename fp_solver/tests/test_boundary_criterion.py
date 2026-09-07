#!/usr/bin/env python3
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from boundary_criterion_scan import (  # noqa: E402
    C_KMS,
    G,
    Norb,
    Norb_all_components,
    Norb_components,
    loss_cone_R,
    loss_cone_orbit_components,
    loss_cone_jump_components,
    total_cross_section,
)


def main() -> int:
    r = np.geomspace(1.0, 100.0, 64)
    rho = np.ones_like(r)
    sigma = np.ones_like(r)

    # An eccentric orbit that leaves the tabulated radial domain must return a
    # two-component sentinel so that both the scalar root finder and the
    # large-angle diagnostic remain well defined.
    components = Norb_components(
        r, rho, sigma, 1.0, 1.0e6, 2.0, 0.9, "constant", 80.0, 128
    )
    assert len(components) == 2
    assert all(math.isnan(value) for value in components)
    assert math.isnan(Norb(
        r, rho, sigma, 1.0, 1.0e6, 2.0, 0.9, "constant", 80.0, 128
    ))
    all_components = Norb_all_components(
        r, rho, sigma, 1.0, 1.0e6, 2.0, 0.9, "constant", 80.0, 128
    )
    assert len(all_components) == 3
    assert all(math.isnan(value) for value in all_components)

    mbh = 4.0e6
    rg = G * mbh / C_KMS**2
    assert math.isclose(loss_cone_R(mbh, 16.0 * rg), 1.0)
    assert math.isclose(loss_cone_R(mbh, 160.0 * rg), 0.1)

    # The critical-orbit helper must return q=N_orb/R_lc and an orbit on the
    # loss-cone edge, independently of the collision-kernel choice.
    rw = np.geomspace(6.0 * rg, 1.0e3, 2048)
    rhod = 1.0e4 * (rw / rw[100]) ** -0.75
    sigd = np.sqrt(G * mbh / (2.0 * rw))
    q, nv, nl, rlc, elc = loss_cone_orbit_components(
        rw, rhod, sigd, 1.0e-4, mbh, 1.0e-2,
        "yukawa-tchannel", 80.0, 4096,
    )
    assert all(math.isfinite(value) for value in (q, nv, nl, rlc, elc))
    assert math.isclose(q, nv / rlc, rel_tol=1.0e-14)
    assert math.isclose(elc * elc, 1.0 - rlc, rel_tol=1.0e-14)
    assert 0.0 < nl < nv

    jump = loss_cone_jump_components(
        rw, rhod, sigd, 1.0e-4, mbh, 1.0e-2,
        "yukawa-tchannel", 80.0, 4096,
    )
    q2, nv2, nl2, nt2, rlc2, elc2, nevent, kick = jump
    assert math.isclose(q2, q, rel_tol=1.0e-14)
    assert math.isclose(nv2, nv, rel_tol=1.0e-14)
    assert math.isclose(nl2, nl, rel_tol=1.0e-14)
    assert nt2 > nv2 > nl2 > 0.0
    assert math.isclose(nevent, nt2 / q2, rel_tol=1.0e-14)
    assert math.isclose(kick**2, (nv2 / nt2) / rlc2, rel_tol=1.0e-14)
    assert math.isclose(elc2, elc, rel_tol=1.0e-14)

    vrel = np.array([0.0, 80.0, 800.0])
    total = total_cross_section(2.0, vrel, "yukawa-tchannel", 80.0)
    assert np.allclose(total, 2.0 / (1.0 + (vrel / 80.0) ** 2))

    print("PASS: boundary orbit-domain sentinel")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
