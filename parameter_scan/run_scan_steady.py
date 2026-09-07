#!/usr/bin/env python3
"""Combine independent operators, solve steady f(E,J), and build flux budgets."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def run_logged(command: list[str], log: Path) -> int:
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
    result_dir = case_root / "flux"
    result_dir.mkdir(parents=True, exist_ok=True)
    status_path = result_dir / "STATUS.json"
    if status_path.exists():
        old = json.loads(status_path.read_text())
        if old.get("status") in {
            "SCAN_CASE_COMPLETE",
            "SCAN_CASE_SKIPPED",
            "SCAN_CASE_FAILED",
        }:
            print(json.dumps(old, indent=2, sort_keys=True))
            return 0

    base = {
        "schema": "sidm-smbh-scan-steady-driver-v1",
        "case_index": args.index,
        "case_tag": case["tag"],
    }
    try:
        prep = json.loads((case_root / "PREP_STATUS.json").read_text())
        admissibility_path = case_root / "FP_ADMISSIBILITY.json"
        if not admissibility_path.is_file():
            raise RuntimeError("missing final FP admissibility classification")
        admissibility = json.loads(admissibility_path.read_text())
        if admissibility.get("status") != "FP_ADMISSIBLE":
            write_status(
                status_path,
                {
                    **base,
                    "status": "SCAN_CASE_SKIPPED",
                    "reason": admissibility.get("status"),
                    "admissibility_reason": admissibility.get("reason"),
                },
            )
            return 0
        if prep.get("status") != "READY_FOR_FP":
            write_status(
                status_path,
                {
                    **base,
                    "status": "SCAN_CASE_SKIPPED",
                    "reason": prep.get("status"),
                    "preparation_reason": prep.get("reason"),
                },
            )
            return 0

        for energy_bins in (64, 96):
            operator_dirs = [
                case_root / "operators" / f"e{energy_bins}_j257_s{seed}"
                for seed in (314159, 271828)
            ]
            for directory in operator_dirs:
                status = json.loads((directory / "STATUS.json").read_text())
                if status.get("status") != "FP_OPERATOR_COMPLETE":
                    raise RuntimeError(
                        f"operator unavailable at {energy_bins} bins: {status}"
                    )

            steady_dir = case_root / "steady" / f"e{energy_bins}_j257"
            steady_dir.mkdir(parents=True, exist_ok=True)
            combined_npz = steady_dir / "operator_combined.npz"
            combined_json = steady_dir / "operator_combined.json"
            command = [
                sys.executable,
                str(args.production_package / "combine_ej_operators.py"),
                str(operator_dirs[0] / "operator.npz"),
                str(operator_dirs[1] / "operator.npz"),
                "--out",
                str(combined_npz),
                "--out-json",
                str(combined_json),
            ]
            if run_logged(command, steady_dir / "01_combine.log") != 0:
                raise RuntimeError(f"operator combination failed for E={energy_bins}")

            operator_map = {
                "s314159": operator_dirs[0] / "operator.npz",
                "s271828": operator_dirs[1] / "operator.npz",
                "combined": combined_npz,
            }
            for label, operator in operator_map.items():
                command = [
                    sys.executable,
                    str(args.production_package / "solve_ej_steady_state.py"),
                    "--operator",
                    str(operator),
                    "--out-json",
                    str(steady_dir / f"steady_{label}.json"),
                    "--out-csv",
                    str(steady_dir / f"steady_{label}.csv"),
                    "--residual-tolerance",
                    "1e-7",
                    "--maximum-evaluations",
                    "10000",
                ]
                if run_logged(command, steady_dir / f"02_solve_{label}.log") != 0:
                    raise RuntimeError(
                        f"steady solver command failed for E={energy_bins}, {label}"
                    )
                solved = json.loads(
                    (steady_dir / f"steady_{label}.json").read_text()
                )
                if solved.get("status") != "EJ_STEADY_SMOKE_PASS":
                    raise RuntimeError(
                        f"steady solver gates failed for E={energy_bins}, {label}"
                    )

        command = [
            sys.executable,
            str(args.scan_package / "measure_fluid_snapshot_flux.py"),
            "--profile",
            prep["remapped_profile"],
            "--production-package",
            str(args.production_package),
            "--record-dir",
            str(result_dir / "fluid_flux_record"),
            "--mbh",
            str(case["black_hole_mass_msun"]),
            "--sigma-over-m",
            str(case["sigma0_over_m_cm2_g"]),
            "--w-kms",
            str(case["yukawa_w_kms"]),
            "--out-json",
            str(result_dir / "fluid_snapshot_flux.json"),
            "--out-csv",
            str(result_dir / "fluid_snapshot_flux_profile.csv"),
        ]
        if run_logged(command, result_dir / "03_fluid_flux.log") != 0:
            raise RuntimeError("black-hole-aware fluid flux measurement failed")

        summary_path = result_dir / "flux_budget.json"
        command = [
            sys.executable,
            str(args.scan_package / "summarize_scan_case.py"),
            "--grid",
            str(args.grid),
            "--index",
            str(args.index),
            "--root",
            str(args.root),
            "--out",
            str(summary_path),
        ]
        if run_logged(command, result_dir / "04_flux_summary.log") != 0:
            raise RuntimeError("flux-budget construction failed")
        summary = json.loads(summary_path.read_text())
        write_status(
            status_path,
            {
                **base,
                "status": "SCAN_CASE_COMPLETE",
                "flux_status": summary["status"],
                "flux_budget": str(summary_path),
                "maximum_numerical_sensitivity": summary[
                    "maximum_numerical_sensitivity"
                ],
            },
        )
    except Exception as exc:
        write_status(
            status_path,
            {
                **base,
                "status": "SCAN_CASE_FAILED",
                "reason": f"{type(exc).__name__}: {exc}",
            },
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
