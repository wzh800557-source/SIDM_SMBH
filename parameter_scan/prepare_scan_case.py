#!/usr/bin/env python3
"""Prepare one independently remapped halo for the FP flux scan."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


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


def write_status(case_root: Path, payload: dict) -> None:
    payload = dict(payload)
    payload["updated_at"] = datetime.now(timezone.utc).isoformat()
    (case_root / "PREP_STATUS.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--source-profile", type=Path, required=True)
    parser.add_argument("--expected-source-sha256")
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--production-package", type=Path, required=True)
    parser.add_argument("--scan-package", type=Path, required=True)
    args = parser.parse_args()

    grid = json.loads(args.grid.read_text())
    case = dict(grid["cases"][args.index])
    case_root = args.root / "cases" / case["tag"]
    case_root.mkdir(parents=True, exist_ok=True)
    existing = case_root / "PREP_STATUS.json"
    if existing.exists():
        old = json.loads(existing.read_text())
        if old.get("status") in {
            "READY_FOR_FP",
            "FP_DOMAIN_NOT_REPRESENTABLE",
        }:
            print(json.dumps(old, indent=2, sort_keys=True))
            return 0

    base = {
        "schema": "sidm-smbh-scan-preparation-v1",
        "case": case,
        "status": "PREPARATION_FAILED",
        "stage": "initialization",
    }
    try:
        if (
            args.expected_source_sha256 is not None
            and sha256(args.source_profile) != args.expected_source_sha256
        ):
            raise RuntimeError("source profile fingerprint differs from production")
        for name in (
            "prepare_bh_remapped_profile.py",
            "hydrostatic_bridge.py",
            "finite_angle_capture.py",
            "combine_ej_operators.py",
            "solve_ej_steady_state.py",
        ):
            if not (args.production_package / name).is_file():
                raise FileNotFoundError(args.production_package / name)

        profile_dir = case_root / "profile"
        # Use a versioned directory so a refined remap never overwrites a
        # failed coarse-step attempt.
        remap_dir = case_root / "remap_refined"
        bridge_dir = case_root / "bridge"
        normalized_root = case_root / "normalized_refined"
        for directory in (profile_dir, remap_dir, bridge_dir, normalized_root):
            directory.mkdir(parents=True, exist_ok=True)

        scaled_profile = profile_dir / "prof_deep_scaled.txt"
        scale_json = profile_dir / "scaling.json"
        relative_bh_fraction = (
            case["black_hole_to_represented_halo_mass"]
            / REFERENCE_BH_TO_HALO_MASS
        )
        force_step = min(0.02, 0.01 / max(relative_bh_fraction, 1.0))
        maximum_growth_steps = max(768, int(math.ceil(10.0 / force_step)))
        base["remap_force_step"] = force_step
        base["remap_maximum_growth_steps"] = maximum_growth_steps
        base["black_hole_fraction_relative_to_fiducial"] = relative_bh_fraction
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
        bridge_rc = run_stage("03_bridge", command, case_root)
        if not bridge_json.is_file():
            raise RuntimeError(f"bridge failed without metadata, return code {bridge_rc}")
        bridge = json.loads(bridge_json.read_text())
        bridge_status = bridge.get("status")
        r_in = bridge.get("r_in_pc")
        rh = bridge.get("rh_pc")
        rg = G_PC_KMS2_MSUN * case["black_hole_mass_msun"] / C_KMS**2
        if bridge_status == "OK":
            interface_status = "SUBGRID_COUPLING_READY"
        elif bridge_status == "SUBGRID_MASS_FRACTION_TOO_LARGE":
            interface_status = "RESOLVED_FP_REPLACEMENT_REQUIRED"
        else:
            interface_status = "NO_USABLE_INNER_CROSSING"

        unavailable_reason = None
        if not isinstance(r_in, (int, float)) or not math.isfinite(r_in):
            unavailable_reason = "NO_USABLE_INNER_CROSSING"
        elif not isinstance(rh, (int, float)) or r_in >= 0.5 * rh:
            unavailable_reason = "MATCHING_RADIUS_OUTSIDE_CONSERVATIVE_KEPLER_DOMAIN"
        elif r_in <= 8.0 * rg:
            unavailable_reason = "MATCHING_RADIUS_REACHES_PLUNGE_DOMAIN"
        elif bridge_status not in ("OK", "SUBGRID_MASS_FRACTION_TOO_LARGE"):
            unavailable_reason = f"BRIDGE_STATUS_{bridge_status}"

        if unavailable_reason is not None:
            base.update(
                {
                    "status": "FP_DOMAIN_NOT_REPRESENTABLE",
                    "stage": "interface_classification",
                    "reason": unavailable_reason,
                    "interface_implementation_status": interface_status,
                    "bridge_status": bridge_status,
                    "bridge_return_code": bridge_rc,
                    "r_in_pc": r_in,
                    "rh_pc": rh,
                    "r_in_over_rh": (
                        r_in / rh
                        if isinstance(r_in, (int, float))
                        and isinstance(rh, (int, float))
                        and rh > 0
                        else None
                    ),
                    "bridge_json": str(bridge_json),
                }
            )
            write_status(case_root, base)
            return 0

        boundary_json = bridge_dir / "boundary_criteria.json"
        command = [
            sys.executable,
            str(args.scan_package / "boundary_criterion_scan.py"),
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
        if run_stage("04_orbit_criteria", command, case_root) != 0:
            raise RuntimeError("orbit-integrated interface calculation failed")
        boundary = json.loads(boundary_json.read_text())
        rb = boundary["criteria"]["local_circular"]["r_boundary_pc"]
        if rb is None or not math.isfinite(rb) or not 8.0 * rg < rb < 0.5 * rh:
            base.update(
                {
                    "status": "FP_DOMAIN_NOT_REPRESENTABLE",
                    "stage": "orbit_interface_classification",
                    "reason": "ORBIT_CRITERION_OUTSIDE_FP_KEPLER_WINDOW",
                    "interface_implementation_status": interface_status,
                    "bridge_status": bridge_status,
                    "r_in_pc": r_in,
                    "orbit_reporting_radius_pc": rb,
                    "rh_pc": rh,
                }
            )
            write_status(case_root, base)
            return 0

        x_at_six_rg = rh / (12.0 * rg)
        emax = max(2.0e7, 1.25 * x_at_six_rg, 20.0 * rh / (2.0 * rb))
        bins = max(104, int(math.ceil(104.0 * math.log10(emax / 0.03) /
                                      math.log10(2.0e7 / 0.03))))
        # The dynamic range and the slope at the reservoir join vary across
        # the mass scan. Preserve the one-per-cent density gate, but refine the
        # energy grid and move the join inward when the first inversion cannot
        # meet it. Every attempt has its own directory and log.
        refined_bins = max(bins, 160)
        highly_refined_bins = max(bins, 208)
        attempts = (
            (bins, max(512, 5 * bins), 0.75 * rh),
            (refined_bins, max(1024, 6 * refined_bins), 0.75 * rh),
            (refined_bins, max(1024, 6 * refined_bins), 0.50 * rh),
            (highly_refined_bins, max(1536, 7 * highly_refined_bins), 0.50 * rh),
        )
        df_attempts = []
        normalization = None
        normalization_json = None
        normalized_dir = None
        selected_grid_bins = None
        for attempt_index, (attempt_bins, fit_points, fit_outer) in enumerate(
            attempts, start=1
        ):
            if fit_outer <= 1.05 * rb:
                df_attempts.append(
                    {
                        "attempt": attempt_index,
                        "grid_bins": attempt_bins,
                        "fit_points": fit_points,
                        "fit_outer_radius_pc": fit_outer,
                        "status": "SKIPPED_JOIN_TOO_CLOSE_TO_BOUNDARY",
                    }
                )
                continue
            attempt_dir = normalized_root / f"attempt_{attempt_index}"
            attempt_dir.mkdir(parents=True, exist_ok=True)
            command = [
                sys.executable,
                str(args.scan_package / "normalize_gnc_df.py"),
                str(bridged_profile),
                "--bridge-json",
                str(bridge_json),
                "--outdir",
                str(attempt_dir),
                "--mbh",
                str(case["black_hole_mass_msun"]),
                "--r-boundary",
                str(rb),
                "--normalization-radius",
                str(rb),
                "--proposal-boundary-radius",
                str(rb),
                "--grid-bins",
                str(attempt_bins),
                "--emax",
                str(emax),
                "--df-fit-points",
                str(fit_points),
                "--df-fit-outer-radius",
                str(fit_outer),
                "--ranks",
                "4",
                "--samples-per-rank",
                "5000",
            ]
            return_code = run_stage(
                f"05_df_normalization_attempt_{attempt_index}",
                command,
                case_root,
            )
            metadata_path = attempt_dir / "df_normalization.json"
            metadata = (
                json.loads(metadata_path.read_text())
                if metadata_path.is_file()
                else None
            )
            passed = return_code == 0 and metadata is not None and metadata.get(
                "status"
            ) == "PASS"
            df_attempts.append(
                {
                    "attempt": attempt_index,
                    "grid_bins": attempt_bins,
                    "fit_points": fit_points,
                    "fit_outer_radius_pc": fit_outer,
                    "return_code": return_code,
                    "status": "PASS" if passed else "FAILED",
                }
            )
            if passed:
                normalized_dir = attempt_dir
                normalization_json = metadata_path
                normalization = metadata
                selected_grid_bins = attempt_bins
                break
        base["df_normalization_attempts"] = df_attempts
        if normalization is None or normalized_dir is None:
            raise RuntimeError(
                "physical DF construction failed every fixed-tolerance refinement"
            )

        base.update(
            {
                "status": "READY_FOR_FP",
                "stage": "complete",
                "interface_implementation_status": interface_status,
                "bridge_status": bridge_status,
                "r_in_pc": r_in,
                "orbit_reporting_radius_pc": rb,
                "rh_pc": rh,
                "r_in_over_rh": r_in / rh,
                "r_in_over_rg": r_in / rg,
                "all_local_crossings_pc": bridge.get("all_N1_crossings_pc", []),
                "first_shell_radius_pc": bridge.get("r0_pc"),
                "N_at_first_fluid_shell": bridge.get("N_at_first_fluid_shell"),
                "fp_owned_mass_fraction_of_first_shell": bridge.get(
                    "mass_inside_r_in_fraction_of_first_shell"
                ),
                "scaled_profile": str(scaled_profile),
                "remapped_profile": str(remapped_profile),
                "bridged_profile": str(bridged_profile),
                "bridge_json": str(bridge_json),
                "boundary_json": str(boundary_json),
                "normalization_json": str(normalization_json),
                "df_table": str(normalized_dir / "f_evolved_normalized.txt"),
                "energy_table_max": normalization.get("xmax"),
                "energy_table_bins": selected_grid_bins,
                "profile_sha256": sha256(bridged_profile),
                "df_table_sha256": sha256(
                    normalized_dir / "f_evolved_normalized.txt"
                ),
            }
        )
        write_status(case_root, base)
        return 0
    except Exception as exc:
        base.update(
            {
                "status": "PREPARATION_FAILED",
                "reason": f"{type(exc).__name__}: {exc}",
            }
        )
        write_status(case_root, base)
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
