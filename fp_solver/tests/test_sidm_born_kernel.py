#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sidm_born_kernel import (  # noqa: E402
    fluid_tchannel_kp,
    gnc_effective_log_corrected,
    tchannel_viscosity_ratio,
)


def main() -> int:
    xlo = np.array([1.0e-6, 3.0e-5, 1.0e-4])
    slo = tchannel_viscosity_ratio(xlo)
    assert np.max(np.abs(slo - 1.0)) < 1.1e-8

    xhi = np.array([30.0, 100.0, 300.0])
    shi = tchannel_viscosity_ratio(xhi)
    asym = 12.0 * (np.log(xhi) - 1.0) / xhi**4
    assert np.max(np.abs(shi / asym - 1.0)) < 0.08

    grid = np.geomspace(1.0e-5, 1.0e3, 2000)
    for p in (3, 5):
        kp = fluid_tchannel_kp(grid, p)
        assert np.all(np.isfinite(kp)) and np.all(kp > 0.0)
        # The published order-one fit has a sub-0.3 per cent regularization
        # ripple below s=0.01, but must decrease thereafter.
        assert np.max(np.diff(kp[grid >= 1.0e-2])) < 0.0
        assert abs(kp[0] - 1.0) < 3.0e-3

    # The mapped GNC multiplier must reproduce the exact raw viscosity moment:
    # sigma_V/sigma0 = 4 Lambda_eff/y^4 = (2/3) S(y).
    v0, w = 0.7, 2.686
    xe = np.geomspace(1.0e-12, 1.0e-7, 20)
    y = np.sqrt(2.0 * xe) * v0 / w
    le = gnc_effective_log_corrected(xe, v0, w)
    assert np.max(np.abs(le / y**4 - 1.0 / 6.0)) < 1.0e-8
    yall = np.geomspace(1.0e-5, 1.0e3, 500)
    xall = 0.5 * (yall * w / v0) ** 2
    mapped = 4.0 * gnc_effective_log_corrected(xall, v0, w) / yall**4
    exact = (2.0 / 3.0) * tchannel_viscosity_ratio(yall)
    assert np.max(np.abs(mapped / exact - 1.0)) < 3.0e-9
    assert np.all(np.diff(gnc_effective_log_corrected(np.geomspace(1e-8, 1e8, 200), v0, w)) > 0)
    print("PASS: Born/Yukawa angular, K3/K5, and GNC asymptotic tests")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
