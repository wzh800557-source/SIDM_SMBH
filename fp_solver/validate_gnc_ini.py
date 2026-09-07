#!/usr/bin/env python3
"""Validate GNC's physical density over the spatial domain owned by GNC."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import h5py
import numpy as np

PC_AU = 206264.80624709636


def table_columns(ds: h5py.Dataset) -> tuple[np.ndarray, np.ndarray]:
    a = np.asarray(ds)
    if a.dtype.names:
        names = {n.strip().upper(): n for n in a.dtype.names}
        return np.asarray(a[names["X"]], float), np.asarray(a[names["FX"]], float)
    if a.ndim == 2 and 2 in a.shape:
        b = a if a.shape[1] == 2 else a.T
        return np.asarray(b[:, 0], float), np.asarray(b[:, 1], float)
    raise ValueError(f"cannot decode HDF5 table {ds.name}: shape={a.shape}, dtype={a.dtype}")


def interp_positive(x: np.ndarray, y: np.ndarray, xq: float) -> float:
    o = np.argsort(x)
    x, y = x[o], y[o]
    good = np.isfinite(x) & np.isfinite(y) & (y > 0.0)
    if np.count_nonzero(good) < 2 or not (x[good][0] <= xq <= x[good][-1]):
        raise ValueError(f"query {xq:g} outside positive density table")
    return float(np.exp(np.interp(xq, x[good], np.log(y[good]))))


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("hdf5", type=Path)
    p.add_argument("--run-manifest", type=Path, required=True)
    p.add_argument("--ini-log", type=Path, required=True)
    p.add_argument("--out-json", type=Path, required=True)
    p.add_argument("--tolerance", type=float, default=0.15)
    p.add_argument("--mc-tolerance", type=float, default=0.25)
    args = p.parse_args()

    m = json.loads(args.run_manifest.read_text())
    rb = float(m["r_in_pc"])
    rh = float(m["rh_pc"])
    rho_b = float(m["rho_boundary_msun_pc3"])
    mpart = float(m["particle_mass_msun"])
    expected_n0 = float(m["n0_pc3"])
    beta = float(m["cusp_slope_beta"])

    with h5py.File(args.hdf5, "r") as f:
        species = "sbh" if "sbh" in f else "star"
        if species not in f or "fden" not in f[species]:
            raise KeyError(f"no aggregate compact-particle fden in {args.hdf5}")
        x, fden = table_columns(f[species]["fden"])
        xs, fden_simu = table_columns(f[species]["fden_simu"])

    # GNC's radial coordinate is log10(r/AU), while fden is number/AU^3.
    log_rb_au = np.log10(rb * PC_AU)
    log_rh_au = np.log10(rh * PC_AU)
    conv = mpart * PC_AU**3
    rho_b_gnc = interp_positive(x, fden, log_rb_au) * conv
    rho_h_gnc = interp_positive(x, fden, log_rh_au) * conv
    rho_b_mc = interp_positive(xs, fden_simu, log_rb_au) * conv
    rho_h_target = expected_n0 * mpart
    rel_b = rho_b_gnc / rho_b - 1.0
    rel_h = rho_h_gnc / rho_h_target - 1.0
    # GNC writes fden_simu before update_weights applies the final physical
    # dms%weight_asym, while all subsequently evolved particles and events use
    # the updated weights.  Recover that explicit factor from the initialization
    # log and validate the same weights that the production run will evolve.
    matches = re.findall(
        r"weight_asym\s*=\s*([+\-0-9.EeDd]+)", args.ini_log.read_text(errors="replace")
    )
    if not matches:
        raise ValueError(f"no final weight_asym in {args.ini_log}")
    weight_asym = float(matches[-1].replace("D", "E").replace("d", "e"))
    rho_b_mc_evolved = rho_b_mc * weight_asym
    rel_mc_raw = rho_b_mc / rho_b - 1.0
    rel_mc = rho_b_mc_evolved / rho_b - 1.0

    # R_h is a useful diagnostic but is outside the inner GNC-owned domain by
    # orders of magnitude.  Gate instead on the full density profile from
    # r_in/20 to r_in, where the hydrostatic bridge and GNC overlap.
    r_inner = np.geomspace(rb / 20.0, rb, 16)
    target_inner = rho_b * (r_inner / rb) ** (-beta)
    gnc_inner = np.array(
        [interp_positive(x, fden, np.log10(rq * PC_AU)) * conv for rq in r_inner]
    )
    rel_inner = gnc_inner / target_inner - 1.0
    inner_max = float(np.max(np.abs(rel_inner)))
    inner_median = float(np.median(np.abs(rel_inner)))

    status = (
        "PASS"
        if max(abs(rel_b), inner_max) <= args.tolerance
        and abs(rel_mc) <= args.mc_tolerance
        else "FAIL"
    )
    diag = {
        "schema": "gnc-ini-density-roundtrip-v1",
        "status": status,
        "hdf5": args.hdf5.name,
        "species_group": species,
        "r_in_pc": rb,
        "rh_pc": rh,
        "rho_boundary_target_msun_pc3": rho_b,
        "rho_boundary_gnc_msun_pc3": rho_b_gnc,
        "rho_boundary_relative_error": rel_b,
        "rho_boundary_direct_mc_msun_pc3": rho_b_mc,
        "rho_boundary_direct_mc_relative_error_before_final_weight_update": rel_mc_raw,
        "final_event_weight_asym": weight_asym,
        "rho_boundary_evolved_particles_msun_pc3": rho_b_mc_evolved,
        "rho_boundary_evolved_particles_relative_error": rel_mc,
        "rho_rh_target_msun_pc3": rho_h_target,
        "rho_rh_gnc_msun_pc3": rho_h_gnc,
        "rho_rh_relative_error": rel_h,
        "rho_rh_role": "diagnostic only; R_h is outside the inner GNC-owned spatial domain",
        "inner_validation_radius_pc": [float(r_inner[0]), float(r_inner[-1])],
        "inner_profile_max_abs_relative_error": inner_max,
        "inner_profile_median_abs_relative_error": inner_median,
        "tolerance": args.tolerance,
        "mc_tolerance": args.mc_tolerance,
    }
    args.out_json.write_text(json.dumps(diag, indent=2, sort_keys=True) + "\n")
    print(json.dumps(diag, indent=2, sort_keys=True))
    return 0 if status == "PASS" else 3


if __name__ == "__main__":
    raise SystemExit(main())
