#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Generate and advance a pure SIDM halo with the native fluid engine."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np
from astropy import units as u

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, os.environ.get("GRAVOTHERMAL_ROOT") or str(ROOT / "fluid_solver" / "GravothermalSIDM"))
from SourcePy import evolve, record


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--shells", type=int, default=64)
    p.add_argument("--steps", type=int, default=20)
    p.add_argument("--r-s-kpc", type=float, default=2.586)
    p.add_argument("--rho-s", type=float, default=0.0194, help="Msun/pc^3")
    p.add_argument("--sigma-over-m", type=float, default=100.0, help="cm^2/g")
    p.add_argument("--t-epsilon", type=float, default=1.0e-4)
    args = p.parse_args()
    if args.shells < 8 or args.steps < 0:
        p.error("use at least eight shells and a nonnegative step count")
    if min(args.r_s_kpc, args.rho_s, args.sigma_over_m, args.t_epsilon) <= 0:
        p.error("physical scales and timestep tolerance must be positive")
    args.outdir.mkdir(parents=True, exist_ok=False)
    h = evolve.Halo(record.HaloRecord(str(args.outdir / "record")),
                    n_shells=args.shells, r_s=args.r_s_kpc, rho_s=args.rho_s,
                    sigma_m_with_units=args.sigma_over_m, M_bh_with_units=0.0,
                    n_adjustment_max=500)
    h.t_epsilon = args.t_epsilon
    h.r_epsilon = 1.0e-10
    initial_density = float(h.rho[0])
    for _ in range(args.steps):
        before = float(h.t)
        h.conduct_heat()
        h.hydrostatic_adjustment()
        if np.max(np.abs(h.delta_r / h.r)) > h.r_epsilon:
            raise RuntimeError("hydrostatic adjustment did not reach its tolerance")
        if not np.isfinite(h.t) or h.t <= before:
            raise RuntimeError("fluid timestep did not advance")
        if np.any(~np.isfinite(h.rho)) or np.any(h.rho <= 0) or np.any(h.p <= 0):
            raise RuntimeError("invalid fluid state")
    profile = np.column_stack((h.r * h.scale_r.to_value(u.pc),
                               h.rho * h.scale_rho.to_value(u.Msun / u.pc**3),
                               np.sqrt(h.p / h.rho) * h.scale_v.to_value(u.km / u.s)))
    np.savetxt(args.outdir / "profile.txt", profile,
               header="r_pc rho_shell_average_msun_pc3 sigma_1d_kms")
    metadata = {"schema": "generated-fluid-halo-v1", "status": "COMPLETE",
                "steps": args.steps, "shells": args.shells, "M_bh_msun": 0.0,
                "r_s_kpc": args.r_s_kpc, "rho_s_msun_pc3": args.rho_s,
                "sigma_over_m_cm2_g": args.sigma_over_m,
                "scattering_model": "constant",
                "time_code": float(h.t), "time_relax": float(h.t / h.t_relax),
                "central_density_ratio_to_initial": float(h.rho[0] / initial_density),
                "scope": "Generated example halo. No claim of deep collapse or production equivalence."}
    (args.outdir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
