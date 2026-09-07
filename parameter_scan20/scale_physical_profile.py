#!/usr/bin/env python3
"""Homologously rescale a physical gravothermal snapshot by halo mass."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def enclosed_mass(r_pc: np.ndarray, rho_msun_pc3: np.ndarray) -> float:
    return float(
        np.sum(
            4.0
            * np.pi
            * rho_msun_pc3
            * np.diff(np.r_[0.0, r_pc**3])
            / 3.0
        )
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", type=Path)
    parser.add_argument("--halo-mass-factor", type=float, required=True)
    parser.add_argument("--out-profile", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    args = parser.parse_args()
    if not np.isfinite(args.halo_mass_factor) or args.halo_mass_factor <= 0.0:
        raise ValueError("halo mass factor must be positive")

    data = np.loadtxt(args.profile, comments="#", ndmin=2)
    if data.shape[1] < 3 or data.shape[0] < 4:
        raise ValueError("profile needs at least four rows and three columns")
    source = np.asarray(data[:, :3], float)
    if np.any(~np.isfinite(source)) or np.any(source <= 0.0):
        raise ValueError("profile values must be positive and finite")
    if np.any(np.diff(source[:, 0]) <= 0.0):
        raise ValueError("profile radii must increase")

    length_factor = args.halo_mass_factor ** (1.0 / 3.0)
    scaled = source.copy()
    scaled[:, 0] *= length_factor
    scaled[:, 2] *= length_factor

    args.out_profile.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        args.out_profile,
        scaled,
        fmt="%.16e",
        header=(
            "r_pc rho_shell_Msun_pc3 sigma1d_kms; homologous fixed-collapse-"
            "stage halo rescaling"
        ),
    )
    source_mass = enclosed_mass(source[:, 0], source[:, 1])
    scaled_mass = enclosed_mass(scaled[:, 0], scaled[:, 1])
    measured_factor = scaled_mass / source_mass
    if abs(measured_factor / args.halo_mass_factor - 1.0) > 2.0e-12:
        raise RuntimeError("enclosed-mass scaling failed")
    virial_residual = np.max(
        np.abs(
            (scaled[:, 2] / source[:, 2]) ** 2
            / (args.halo_mass_factor / length_factor)
            - 1.0
        )
    )
    if virial_residual > 2.0e-12:
        raise RuntimeError("velocity scaling is not virial-homologous")

    result = {
        "schema": "homologous-physical-profile-v1",
        "status": "PASS",
        "source_profile": args.profile.name,
        "source_profile_sha256": sha256(args.profile),
        "output_profile": args.out_profile.name,
        "output_profile_sha256": sha256(args.out_profile),
        "halo_mass_factor": args.halo_mass_factor,
        "radius_factor": length_factor,
        "density_factor": 1.0,
        "one_dimensional_dispersion_factor": length_factor,
        "source_represented_mass_msun": source_mass,
        "scaled_represented_mass_msun": scaled_mass,
        "measured_mass_factor": measured_factor,
        "virial_scaling_max_relative_residual": float(virial_residual),
        "rows": int(scaled.shape[0]),
    }
    args.out_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
