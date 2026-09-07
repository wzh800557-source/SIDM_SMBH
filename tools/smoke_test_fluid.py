#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Exercise native fluid evolution, BH gravity, and inner-face energy accounting."""
from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "fp_solver"))
from measured_fluid_feedback import MeasuredBoundaryHalo, evolve, record


def main() -> int:
    results = []
    with tempfile.TemporaryDirectory(prefix="sidm-fluid-smoke-") as temporary:
        for name, mbh, sink in [("pure_halo", 0.0, False),
                                ("bh_control", 4.0e6, False),
                                ("bh_test_current", 4.0e6, True)]:
            h = MeasuredBoundaryHalo(record.HaloRecord(str(Path(temporary) / name)),
                                     n_shells=48, r_min=0.1, n_adjustment_max=500,
                                     M_bh_with_units=mbh, sigma_m_with_units=100.0)
            h.t_epsilon = 1.0e-4
            h.r_epsilon = 1.0e-10
            mass0 = h.m.copy()
            dm = np.diff(np.r_[0.0, h.m])
            luminosity = -1.0e-6 * float(np.dot(dm, h.u)) / h.t_relax if sink else 0.0
            h.set_measured_boundary(luminosity, 0.0)
            for _ in range(8):
                before = float(h.t)
                h.conduct_heat()
                h.hydrostatic_adjustment()
                assert np.max(np.abs(h.delta_r / h.r)) <= h.r_epsilon
                assert np.isfinite(h.t) and h.t > before
                assert np.all(np.isfinite(h.rho)) and np.all(h.rho > 0)
                assert np.all(np.isfinite(h.p)) and np.all(h.p > 0)
            np.testing.assert_array_equal(h.m, mass0)
            h.set_hydrostatic_coefficients()
            expected_mass = h.m + h.M_bh * h.r**3 / (h.r**2 + h.r_soft**2)**1.5
            np.testing.assert_allclose(h.m_ext[1:], expected_mass, rtol=1.0e-13)
            error = abs(h.E_conduction_budget_error_code) / max(h.E_conduction_budget_scale_code, 1.0e-300)
            assert error < 1.0e-10
            if sink:
                assert h.E_boundary_code < 0.0
            results.append({"branch": name, "steps": 8, "time_code": float(h.t),
                            "conduction_budget_relative_error": float(error)})
    print(json.dumps({"status": "PASS", "native_engine": "GravothermalSIDM",
                      "branches": results,
                      "scope": "Installation and integration smoke test. The imposed test current is synthetic, not an FP measurement."}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
