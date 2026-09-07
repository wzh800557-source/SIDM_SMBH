#!/usr/bin/env python3
"""Measure the instantaneous gravothermal flux profile for every scan halo."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def write_status(path: Path, payload: dict) -> None:
    payload = dict(payload)
    payload["updated_at"] = datetime.now(timezone.utc).isoformat()
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)


def locate_remapped_profile(case_root: Path, prep: dict) -> Path | None:
    candidates = []
    if prep.get("remapped_profile"):
        candidates.append(Path(prep["remapped_profile"]))
    candidates.extend(
        (
            case_root / "remap_refined" / "profile_bh.txt",
            case_root / "remap" / "profile_bh.txt",
        )
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--production-package", type=Path, required=True)
    parser.add_argument("--scan-package", type=Path, required=True)
    args = parser.parse_args()

    grid = json.loads(args.grid.read_text())
    case = dict(grid["cases"][args.index])
    case_root = args.root / "cases" / case["tag"]
    outdir = case_root / "all_halo_flux"
    outdir.mkdir(parents=True, exist_ok=True)
    status_path = outdir / "STATUS.json"

    prep_path = case_root / "PREP_STATUS.json"
    if not prep_path.is_file():
        write_status(
            status_path,
            {
                "schema": "sidm-smbh-all-halo-fluid-flux-v1",
                "status": "ALL_HALO_FLUX_FAILED",
                "case_index": args.index,
                "case_tag": case["tag"],
                "reason": "missing preparation status",
            },
        )
        return 0
    prep = json.loads(prep_path.read_text())
    profile = locate_remapped_profile(case_root, prep)
    if profile is None:
        write_status(
            status_path,
            {
                "schema": "sidm-smbh-all-halo-fluid-flux-v1",
                "status": "ALL_HALO_FLUX_FAILED",
                "case_index": args.index,
                "case_tag": case["tag"],
                "preparation_status": prep.get("status"),
                "reason": "no completed black-hole-remapped profile",
            },
        )
        return 0

    expected_sha = sha256(profile)
    if status_path.is_file():
        old = json.loads(status_path.read_text())
        if (
            old.get("status") == "ALL_HALO_FLUX_COMPLETE"
            and old.get("profile_sha256") == expected_sha
        ):
            print(json.dumps(old, indent=2, sort_keys=True))
            return 0

    out_json = outdir / "fluid_snapshot_flux.json"
    out_csv = outdir / "fluid_snapshot_flux_profile.csv"
    command = [
        sys.executable,
        str(args.scan_package / "measure_fluid_snapshot_flux.py"),
        "--profile",
        str(profile),
        "--production-package",
        str(args.production_package),
        "--record-dir",
        str(outdir / "record"),
        "--mbh",
        str(case["black_hole_mass_msun"]),
        "--sigma-over-m",
        str(case["sigma0_over_m_cm2_g"]),
        "--w-kms",
        str(case["yukawa_w_kms"]),
        "--out-json",
        str(out_json),
        "--out-csv",
        str(out_csv),
    ]
    with (outdir / "measure.log").open("w") as stream:
        stream.write("COMMAND " + " ".join(command) + "\n")
        stream.flush()
        completed = subprocess.run(
            command,
            stdout=stream,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        stream.write(f"RETURN_CODE {completed.returncode}\n")

    if completed.returncode != 0 or not out_json.is_file():
        write_status(
            status_path,
            {
                "schema": "sidm-smbh-all-halo-fluid-flux-v1",
                "status": "ALL_HALO_FLUX_FAILED",
                "case_index": args.index,
                "case_tag": case["tag"],
                "preparation_status": prep.get("status"),
                "profile": str(profile),
                "profile_sha256": expected_sha,
                "reason": f"measurement returned {completed.returncode}",
            },
        )
        return 0

    measured = json.loads(out_json.read_text())
    if measured.get("status") != "PASS":
        write_status(
            status_path,
            {
                "schema": "sidm-smbh-all-halo-fluid-flux-v1",
                "status": "ALL_HALO_FLUX_FAILED",
                "case_index": args.index,
                "case_tag": case["tag"],
                "preparation_status": prep.get("status"),
                "profile": str(profile),
                "profile_sha256": expected_sha,
                "reason": f"unexpected measurement status {measured.get('status')}",
            },
        )
        return 0

    write_status(
        status_path,
        {
            "schema": "sidm-smbh-all-halo-fluid-flux-v1",
            "status": "ALL_HALO_FLUX_COMPLETE",
            "case_index": args.index,
            "case_tag": case["tag"],
            "preparation_status": prep.get("status"),
            "profile": str(profile),
            "profile_sha256": expected_sha,
            "flux_json": str(out_json),
            "radial_flux_csv": str(out_csv),
            "first_shell_luminosity_msun_kms2_per_myr": measured[
                "first_fluid_shell_luminosity_msun_kms2_per_myr"
            ],
            "maximum_absolute_luminosity_msun_kms2_per_myr": measured[
                "maximum_absolute_luminosity_msun_kms2_per_myr"
            ],
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
