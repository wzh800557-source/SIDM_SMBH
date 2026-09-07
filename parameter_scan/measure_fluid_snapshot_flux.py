#!/usr/bin/env python3
"""Measure the black-hole-aware gravothermal luminosity of one snapshot."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from pathlib import Path

import astropy.units as u
import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", type=Path, required=True)
    parser.add_argument("--production-package", type=Path, required=True)
    parser.add_argument("--record-dir", type=Path, required=True)
    parser.add_argument("--mbh", type=float, required=True)
    parser.add_argument("--sigma-over-m", type=float, required=True)
    parser.add_argument("--w-kms", type=float, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--out-csv", type=Path, required=True)
    args = parser.parse_args()

    sys.path.insert(0, str(args.production_package))
    # The response adapter finds the bundled engine unless an override is set.
    from measured_fluid_feedback import (  # pylint: disable=import-error,import-outside-toplevel
        MeasuredBoundaryHalo,
        install_physical_profile,
        load_profile,
        record,
    )

    profile = load_profile(args.profile)
    args.record_dir.mkdir(parents=True, exist_ok=True)
    h = MeasuredBoundaryHalo(
        record.HaloRecord(str(args.record_dir)),
        sigma_m_with_units=args.sigma_over_m,
        w_units=args.w_kms,
        model_elastic_scattering_lmfp=(
            "YukawaBornViscosityApproxTchannel_K3_order1"
        ),
        model_elastic_scattering_smfp=(
            "YukawaBornViscosityApproxTchannel_K5_order1"
        ),
        M_bh_with_units=args.mbh,
        r_soft=0.0,
        n_shells=profile.shape[0],
        n_adjustment_max=500,
    )
    installation = install_physical_profile(h, profile)
    h.update_derived_parameters()

    r_pc = h.r * h.scale_r.to_value(u.pc)
    rho = h.rho * h.scale_rho.to_value(u.Msun / u.pc**3)
    sigma = h.v * h.scale_v.to_value(u.km / u.s)
    enclosed_mass = h.m * h.scale_m.to_value(u.Msun)
    luminosity = h.L * h.scale_L.to_value(
        u.Msun * (u.km / u.s) ** 2 / u.Myr
    )
    fields = (
        "radius_pc",
        "rho_shell_msun_pc3",
        "sigma_1d_kms",
        "enclosed_mass_msun",
        "gravothermal_luminosity_msun_kms2_per_myr",
        "lmfp_scaleheight_factor",
    )
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        writer.writerows(
            zip(
                r_pc,
                rho,
                sigma,
                enclosed_mass,
                luminosity,
                h.lmfp_scaleheight_factor,
            )
        )

    finite = np.isfinite(luminosity)
    result = {
        "schema": "black-hole-aware-fluid-snapshot-flux-v1",
        "status": "PASS",
        "profile": args.profile.name,
        "profile_sha256": sha256(args.profile),
        "black_hole_mass_msun": args.mbh,
        "sigma0_over_m_cm2_g": args.sigma_over_m,
        "yukawa_w_kms": args.w_kms,
        "profile_installation": installation,
        "first_fluid_shell_radius_pc": float(r_pc[0]),
        "first_fluid_shell_density_msun_pc3": float(rho[0]),
        "first_fluid_shell_sigma_1d_kms": float(sigma[0]),
        "first_fluid_shell_luminosity_msun_kms2_per_myr": float(luminosity[0]),
        "maximum_outward_luminosity_msun_kms2_per_myr": float(
            np.max(luminosity[finite])
        ),
        "minimum_signed_luminosity_msun_kms2_per_myr": float(
            np.min(luminosity[finite])
        ),
        "maximum_absolute_luminosity_msun_kms2_per_myr": float(
            np.max(np.abs(luminosity[finite]))
        ),
        "radius_of_maximum_absolute_luminosity_pc": float(
            r_pc[np.nanargmax(np.abs(luminosity))]
        ),
        "lmfp_scaleheight_factor_inner": float(h.lmfp_scaleheight_factor[0]),
        "lmfp_scaleheight_factor_minimum": float(
            np.min(h.lmfp_scaleheight_factor)
        ),
        "radial_flux_csv": args.out_csv.name,
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
