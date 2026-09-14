#!/usr/bin/env python3
"""Prepare a two-grid, two-seed GNC ensemble for the energy-current gate."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
FP = ROOT / "fp_solver"


def run(arguments: list[object]) -> None:
    subprocess.run([sys.executable, *map(str, arguments)], check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", type=Path, required=True)
    parser.add_argument("--bridge-json", type=Path, required=True)
    parser.add_argument("--base-model", type=Path, required=True)
    parser.add_argument("--gnc-build", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--grid", type=int, action="append", default=[])
    parser.add_argument("--seed", type=int, action="append", default=[])
    parser.add_argument("--ranks", type=int, default=8)
    parser.add_argument("--samples-per-rank", type=int, default=8000)
    parser.add_argument("--normalization-seed", type=int, default=73421)
    parser.add_argument("--dt-tnr", type=float, default=0.1)
    parser.add_argument("--total-tnr", type=float, default=40.0)
    parser.add_argument("--burn-in-tnr", type=float, default=10.0)
    parser.add_argument("--updates-per-snapshot", type=int, default=10)
    parser.add_argument("--sigma-over-m", type=float, default=100.0)
    parser.add_argument("--w-kms", type=float, default=80.0)
    parser.add_argument("--importance-uniform-fraction", type=float, default=0.05)
    parser.add_argument("--importance-capture-fraction", type=float, default=0.15)
    parser.add_argument("--importance-boundary-fed-fraction", type=float, default=0.70)
    parser.add_argument("--importance-boundary-fed-xmin-factor", type=float, default=0.5)
    args = parser.parse_args()

    grids = args.grid or [104, 128]
    seeds = args.seed or [271828, 314159]
    if len(set(grids)) < 2 or len(set(seeds)) < 2:
        parser.error("the production gate requires at least two grids and two seeds")
    if any(grid % args.ranks for grid in grids):
        parser.error("each grid is also used for dc-bins and must divide by ranks")
    if not 0.0 < args.burn_in_tnr < args.total_tnr:
        parser.error("burn-in-tnr must lie between zero and total-tnr")
    cfs = args.gnc_build / "common_data" / "cfuns_34.bin"
    executables = [args.gnc_build / "main" / name for name in ("ini", "main", "pro")]
    for path in [args.profile, args.bridge_json, args.base_model, cfs, *executables]:
        if not path.is_file():
            raise FileNotFoundError(path)
    bridge = json.loads(args.bridge_json.read_text())
    if bridge.get("status") != "OK":
        raise ValueError("hydrostatic bridge did not pass")
    boundary = float(bridge["r_in_pc"])
    mbh_msun = float(bridge["mbh_msun"])
    # The patched GNC initializer constructs mbh / Weight_n particles on each
    # MPI rank.  Match that native count to the imported phase-space catalogue
    # so increasing --samples-per-rank changes the Monte-Carlo resolution
    # without changing the represented physical mass.
    particle_weight_msun = mbh_msun / args.samples_per_rank
    args.run_root.mkdir(parents=True, exist_ok=False)

    members = []
    for grid in sorted(set(grids)):
        normalized = args.run_root / f"normalized_gx{grid}"
        run([
            FP / "normalize_gnc_df.py",
            args.profile,
            "--bridge-json", args.bridge_json,
            "--outdir", normalized,
            "--mbh", mbh_msun,
            "--r-boundary", boundary,
            "--normalization-radius", boundary,
            "--proposal-boundary-radius", boundary,
            "--grid-bins", grid,
            "--emax", 2.0e7,
            "--df-fit-points", 512,
            "--ranks", args.ranks,
            "--samples-per-rank", args.samples_per_rank,
            "--weight-n", particle_weight_msun,
            "--seed", args.normalization_seed,
            "--importance-uniform-fraction", args.importance_uniform_fraction,
            "--importance-capture-fraction", args.importance_capture_fraction,
            "--importance-boundary-fed-fraction", (
                args.importance_boundary_fed_fraction
            ),
            "--importance-boundary-fed-xmin-factor", (
                args.importance_boundary_fed_xmin_factor
            ),
        ])
        mfrac_lines = (normalized / "mfrac.normalized.in").read_text().splitlines()
        component_fields = next(
            line.split()
            for line in mfrac_lines
            if line.strip() and not line.lstrip().startswith("#")
            and len(line.split()) >= 7 and float(line.split()[0]) > 0.0
        )
        written_particle_weight = float(component_fields[4])
        implied_samples = mbh_msun / written_particle_weight
        if not math.isclose(
            implied_samples, args.samples_per_rank, rel_tol=1.0e-7, abs_tol=0.0
        ):
            raise RuntimeError(
                "GNC particle count does not match imported catalogue: "
                f"{implied_samples:g} != {args.samples_per_rank:d}"
            )
        for seed in sorted(set(seeds)):
            tag = f"gx{grid}_dc{grid}_s{seed}"
            run_dir = args.run_root / tag
            run([
                FP / "prepare_absolute_run.py",
                "--base-model", args.base_model,
                "--normalized-dir", normalized,
                "--run-dir", run_dir,
                "--ranks", args.ranks,
                "--gx-bins", grid,
                "--dc-bins", grid,
                "--dt-tnr", args.dt_tnr,
                "--total-tnr", args.total_tnr,
                "--updates-per-snapshot", args.updates_per_snapshot,
                "--gnc-seed", seed,
                "--sidm-kernel", "yukawa-tchannel",
                "--sigma-over-m", args.sigma_over_m,
                "--w-kms", args.w_kms,
                "--cfs-file", cfs,
            ])
            members.append({
                "array_index": len(members),
                "tag": tag,
                "run_dir": str(run_dir),
                "grid": [grid, grid],
                "seed": seed,
                "ranks": args.ranks,
            })

    manifest = {
        "schema": "gnc-energy-current-ensemble-v1",
        "status": "PREFLIGHT_PASS",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "profile": args.profile.name,
        "bridge": args.bridge_json.name,
        "r_in_pc": boundary,
        "mbh_msun": mbh_msun,
        "samples_per_rank": args.samples_per_rank,
        "particle_weight_msun": particle_weight_msun,
        "grids": sorted(set(grids)),
        "seeds": sorted(set(seeds)),
        "dt_tnr": args.dt_tnr,
        "total_tnr": args.total_tnr,
        "burn_in_tnr_for_gate": args.burn_in_tnr,
        "importance_proposal": {
            "uniform_fraction": args.importance_uniform_fraction,
            "initial_cusp_fraction": args.importance_capture_fraction,
            "boundary_fed_fraction": args.importance_boundary_fed_fraction,
            "boundary_fed_xmin_factor": args.importance_boundary_fed_xmin_factor,
        },
        "sigma0_over_m_cm2_g": args.sigma_over_m,
        "w_kms": args.w_kms,
        "members": members,
    }
    path = args.run_root / "ensemble_manifest.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
