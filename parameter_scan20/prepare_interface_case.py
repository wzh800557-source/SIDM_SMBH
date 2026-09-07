#!/usr/bin/env python3
"""Remap one halo and measure its collision-defined FP-fluid interface.

This high-resolution stage deliberately stops before distribution-function
normalization or Monte Carlo operator construction.  It measures the exact
interface geometry needed for the dense heat maps without pretending that a
1200-case FP flux calculation has been performed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


REFERENCE_SOURCE_SHA256 = (
    "a8199391def81777c58e6d2b1d80c0be4ef242d9730c276be95e654fa95c6103"
)
REFERENCE_BH_TO_HALO_MASS = 4.0e6 / 1.5283226426126614e10
G_PC_KMS2_MSUN = 4.30091e-3
C_KMS = 2.99792458e5


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def run_stage(name: str, command: list[str], case_root: Path) -> int:
    log = case_root / f"{name}.log"
    with log.open("w") as stream:
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
    return int(completed.returncode)


def write_status(path: Path, payload: dict) -> None:
    payload = dict(payload)
    payload["updated_at"] = datetime.now(timezone.utc).isoformat()
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)


def finite(value: object) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--source-profile", type=Path, required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--production-package", type=Path, required=True)
    parser.add_argument("--scan-package", type=Path, required=True)
    args = parser.parse_args()

    grid = json.loads(args.grid.read_text())
    case = dict(grid["cases"][args.index])
    case_root = args.root / "cases" / case["tag"]
    case_root.mkdir(parents=True, exist_ok=True)
    status_path = case_root / "PREP_STATUS.json"
    if status_path.is_file():
        old = json.loads(status_path.read_text())
        old_profile = Path(old.get("remapped_profile", ""))
        if old.get("status") == "INTERFACE_COMPLETE" and old_profile.is_file():
            print(json.dumps(old, indent=2, sort_keys=True))
            return 0

    base = {
        "schema": "sidm-smbh-interface-scan-20x20-v1",
        "case": case,
        "status": "INTERFACE_PREPARATION_FAILED",
        "stage": "initialization",
    }
    try:
        if sha256(args.source_profile) != REFERENCE_SOURCE_SHA256:
            raise RuntimeError("source profile fingerprint differs from production")
        for name in ("prepare_bh_remapped_profile.py", "hydrostatic_bridge.py"):
            if not (args.production_package / name).is_file():
                raise FileNotFoundError(args.production_package / name)

        profile_dir = case_root / "profile"
        remap_dir = case_root / "remap"
        bridge_dir = case_root / "bridge"
        for directory in (profile_dir, remap_dir, bridge_dir):
            directory.mkdir(parents=True, exist_ok=True)

        scaled_profile = profile_dir / "prof_deep_scaled.txt"
        scale_json = profile_dir / "scaling.json"
        relative_bh_fraction = (
            case["black_hole_to_represented_halo_mass"]
            / REFERENCE_BH_TO_HALO_MASS
        )
        force_step = min(0.02, 0.01 / max(relative_bh_fraction, 1.0))
        maximum_growth_steps = max(768, int(math.ceil(10.0 / force_step)))

        command = [
            sys.executable,
            str(args.scan_package / "scale_physical_profile.py"),
            str(args.source_profile),
            "--halo-mass-factor",
            str(case["halo_mass_factor"]),
            "--out-profile",
            str(scaled_profile),
            "--out-json",
            str(scale_json),
        ]
        if run_stage("01_scale", command, case_root) != 0:
            raise RuntimeError("homologous profile scaling failed")

        remapped_profile = remap_dir / "profile_bh.txt"
        remap_json = remap_dir / "remap.json"
        command = [
            sys.executable,
            str(args.production_package / "prepare_bh_remapped_profile.py"),
            "--profile",
            str(scaled_profile),
            "--out-profile",
            str(remapped_profile),
            "--out-json",
            str(remap_json),
            "--record-dir",
            str(remap_dir / "record"),
            "--mbh",
            str(case["black_hole_mass_msun"]),
            "--sigma-over-m",
            str(case["sigma0_over_m_cm2_g"]),
            "--w-kms",
            str(case["yukawa_w_kms"]),
            "--max-force-fraction",
            str(force_step),
            "--r-epsilon",
            "1e-9",
            "--max-growth-steps",
            str(maximum_growth_steps),
        ]
        if run_stage("02_bh_remap", command, case_root) != 0:
            raise RuntimeError("adiabatic black-hole remap failed")

        bridged_profile = bridge_dir / "bridged_profile.txt"
        bridge_json = bridge_dir / "bridge.json"
        command = [
            sys.executable,
            str(args.production_package / "hydrostatic_bridge.py"),
            str(remapped_profile),
            "--out-profile",
            str(bridged_profile),
            "--out-json",
            str(bridge_json),
            "--mbh",
            str(case["black_hole_mass_msun"]),
            "--sigma-over-m",
            str(case["sigma0_over_m_cm2_g"]),
            "--mode",
            "adiabatic-mass-matched",
            "--kernel",
            "yukawa-tchannel",
            "--w-kms",
            str(case["yukawa_w_kms"]),
        ]
        bridge_return_code = run_stage("03_bridge", command, case_root)
        if not bridge_json.is_file() or not bridged_profile.is_file():
            raise RuntimeError(
                f"bridge produced no usable output, return code {bridge_return_code}"
            )
        bridge = json.loads(bridge_json.read_text())

        boundary_json = bridge_dir / "boundary_criteria.json"
        command = [
            sys.executable,
            str(args.scan_package / "boundary_criterion_scan_scan.py"),
            str(bridged_profile),
            "--bridge-json",
            str(bridge_json),
            "--out",
            str(boundary_json),
            "--mbh",
            str(case["black_hole_mass_msun"]),
            "--sigma-over-m",
            str(case["sigma0_over_m_cm2_g"]),
            "--kernel",
            "yukawa-tchannel",
            "--w-kms",
            str(case["yukawa_w_kms"]),
        ]
        boundary_return_code = run_stage("04_orbit_criteria", command, case_root)
        boundary = (
            json.loads(boundary_json.read_text())
            if boundary_return_code == 0 and boundary_json.is_file()
            else {}
        )

        local = boundary.get("criteria", {}).get("local_circular", {})
        bridge_r_in = bridge.get("r_in_pc")
        reporting_radius = local.get("r_boundary_pc")
        r_in = reporting_radius if finite(reporting_radius) else bridge_r_in
        rh = boundary.get("rh_pc", bridge.get("rh_pc"))
        rg = G_PC_KMS2_MSUN * case["black_hole_mass_msun"] / C_KMS**2
        bridge_status = bridge.get("status")

        if bridge_status == "OK":
            implementation = "SUBGRID_COUPLING_READY"
        elif bridge_status == "SUBGRID_MASS_FRACTION_TOO_LARGE":
            implementation = "RESOLVED_FP_REPLACEMENT_REQUIRED"
        else:
            implementation = "NO_USABLE_INNER_CROSSING"

        admissible = True
        reason = "GEOMETRIC_FP_DOMAIN_PASS"
        if not finite(r_in) or not finite(rh) or float(rh) <= 0.0:
            admissible = False
            reason = "NO_USABLE_INNER_CROSSING"
        elif float(r_in) >= 0.5 * float(rh):
            admissible = False
            reason = "MATCHING_RADIUS_OUTSIDE_CONSERVATIVE_KEPLER_DOMAIN"
        elif float(r_in) <= 8.0 * rg:
            admissible = False
            reason = "MATCHING_RADIUS_REACHES_PLUNGE_DOMAIN"
        elif bridge_status not in ("OK", "SUBGRID_MASS_FRACTION_TOO_LARGE"):
            admissible = False
            reason = f"BRIDGE_STATUS_{bridge_status}"
        elif not local:
            admissible = False
            reason = "ORBIT_CRITERION_UNAVAILABLE"

        all_roots = local.get("all_roots_pc", [])
        outer_crossing = None
        if isinstance(all_roots, list):
            finite_roots = sorted(float(value) for value in all_roots if finite(value))
            if finite_roots:
                r_in = finite_roots[0]
                outer_crossing = finite_roots[-1] if len(finite_roots) > 1 else None
                if finite(rh) and float(r_in) >= 0.5 * float(rh):
                    admissible = False
                    reason = "MATCHING_RADIUS_OUTSIDE_CONSERVATIVE_KEPLER_DOMAIN"

        write_status(
            status_path,
            {
                **base,
                "status": "INTERFACE_COMPLETE",
                "stage": "complete",
                "fp_geometric_admissible": admissible,
                "fp_geometric_reason": reason,
                "interface_implementation_status": implementation,
                "bridge_status": bridge_status,
                "bridge_return_code": bridge_return_code,
                "boundary_return_code": boundary_return_code,
                "scaled_profile": str(scaled_profile),
                "remapped_profile": str(remapped_profile),
                "bridged_profile": str(bridged_profile),
                "bridge_json": str(bridge_json),
                "boundary_json": str(boundary_json) if boundary_json.is_file() else None,
                "r_in_pc": r_in,
                "r_out_pc": outer_crossing,
                "rh_pc": rh,
                "r_in_over_rh": (
                    float(r_in) / float(rh)
                    if finite(r_in) and finite(rh) and float(rh) > 0.0
                    else None
                ),
                "r_in_over_rg": float(r_in) / rg if finite(r_in) else None,
                "N_at_first_fluid_shell": bridge.get("N_at_first_fluid_shell"),
                "first_shell_radius_pc": bridge.get("r0_pc"),
                "rho_in_msun_pc3": bridge.get("rho_in_msun_pc3"),
                "sigma_in_kms": bridge.get("sigma_in_kms"),
                "fp_boundary_angular_kernel_status": local.get(
                    "fp_angular_kernel_status"
                ),
                "fp_loss_cone_kernel_status": boundary.get(
                    "fp_loss_cone_angular_kernel_gate", {}
                ).get("status"),
                "remap_force_step": force_step,
                "remap_maximum_growth_steps": maximum_growth_steps,
            },
        )
    except Exception as exc:
        write_status(
            status_path,
            {
                **base,
                "status": "INTERFACE_PREPARATION_FAILED",
                "reason": f"{type(exc).__name__}: {exc}",
            },
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
