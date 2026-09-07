#!/usr/bin/env python3
"""Create the single post-BH fluid snapshot shared by GNC and fluid runs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import astropy.units as u
import numpy as np

from measured_fluid_feedback import (
    MeasuredBoundaryHalo,
    adiabatically_grow_black_hole,
    install_physical_profile,
    load_profile,
    record,
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--profile", type=Path, required=True)
    p.add_argument("--out-profile", type=Path, required=True)
    p.add_argument("--out-json", type=Path, required=True)
    p.add_argument("--record-dir", type=Path, required=True)
    p.add_argument("--mbh", type=float, default=4.0e6)
    p.add_argument("--sigma-over-m", type=float, default=100.0)
    p.add_argument("--w-kms", type=float, default=80.0)
    p.add_argument("--max-force-fraction", type=float, default=0.05)
    p.add_argument("--r-epsilon", type=float, default=1.0e-9)
    p.add_argument("--adjustment-max", type=int, default=500)
    p.add_argument("--max-growth-steps", type=int, default=512)
    args = p.parse_args()

    profile = load_profile(args.profile)
    args.out_profile.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.record_dir.parent.mkdir(parents=True, exist_ok=True)

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
        n_adjustment_max=args.adjustment_max,
    )
    h.M_bh = 0.0
    profile_check = install_physical_profile(h, profile)
    remap = adiabatically_grow_black_hole(
        h,
        args.mbh,
        args.max_force_fraction,
        args.r_epsilon,
        args.max_growth_steps,
    )

    post = np.column_stack((
        h.r * h.scale_r.to_value(u.pc),
        h.rho * h.scale_rho.to_value(u.Msun / u.pc**3),
        h.v * h.scale_v.to_value(u.km / u.s),
    ))
    np.savetxt(
        args.out_profile,
        post,
        fmt="%.16e",
        header=(
            "r_pc rho_shell_Msun_pc3 sigma1d_kms; adiabatic post-BH "
            "Lagrangian fluid snapshot"
        ),
    )
    result = {
        "schema": "adiabatic-bh-remapped-fluid-profile-v1",
        "status": "BH_REMAP_COMPLETE",
        "source_profile": args.profile.name,
        "source_profile_sha256": sha256(args.profile),
        "output_profile": args.out_profile.name,
        "output_profile_sha256": sha256(args.out_profile),
        "M_bh_msun": args.mbh,
        "sigma0_over_m_cm2_g": args.sigma_over_m,
        "yukawa_w_kms": args.w_kms,
        "profile_installation": profile_check,
        "black_hole_remap": remap,
        "lmfp_scaleheight_factor_inner": float(h.lmfp_scaleheight_factor[0]),
        "lmfp_scaleheight_factor_min": float(np.min(h.lmfp_scaleheight_factor)),
        "post_profile_rows": int(post.shape[0]),
    }
    args.out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)
    print("BH_REMAP_COMPLETE", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
