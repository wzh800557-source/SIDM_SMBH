#!/usr/bin/env python3
"""Aggregate the dense interface and fluid-luminosity scan."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def read_json(path: Path) -> dict | None:
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    grid = json.loads(args.grid.read_text())
    rows = []
    prep_counts: Counter[str] = Counter()
    fluid_counts: Counter[str] = Counter()

    for case in grid["cases"]:
        case_root = args.root / "cases" / case["tag"]
        prep = read_json(case_root / "PREP_STATUS.json") or {"status": "MISSING"}
        fluid_status = read_json(case_root / "all_halo_flux" / "STATUS.json") or {
            "status": "MISSING"
        }
        fluid = read_json(
            case_root / "all_halo_flux" / "fluid_snapshot_flux.json"
        )
        prep_counts[str(prep.get("status", "MISSING"))] += 1
        fluid_counts[str(fluid_status.get("status", "MISSING"))] += 1

        row = {
            **case,
            "prep_status": prep.get("status"),
            "prep_reason": prep.get("reason"),
            "fp_geometric_admissible": prep.get("fp_geometric_admissible"),
            "fp_geometric_reason": prep.get("fp_geometric_reason"),
            "interface_implementation_status": prep.get(
                "interface_implementation_status"
            ),
            "bridge_status": prep.get("bridge_status"),
            "r_in_pc": prep.get("r_in_pc"),
            "r_out_pc": prep.get("r_out_pc"),
            "r_h_pc": prep.get("rh_pc"),
            "r_in_over_r_h": prep.get("r_in_over_rh"),
            "r_in_over_r_g": prep.get("r_in_over_rg"),
            "N_at_first_fluid_shell": prep.get("N_at_first_fluid_shell"),
            "rho_boundary_msun_pc3": prep.get("rho_in_msun_pc3"),
            "sigma_boundary_kms": prep.get("sigma_in_kms"),
            "fp_boundary_angular_kernel_status": prep.get(
                "fp_boundary_angular_kernel_status"
            ),
            "fp_loss_cone_kernel_status": prep.get(
                "fp_loss_cone_kernel_status"
            ),
            "all_halo_flux_status": fluid_status.get("status"),
        }
        if fluid is not None and fluid.get("status") == "PASS":
            row.update(
                {
                    "first_shell_gravothermal_luminosity_msun_kms2_per_myr": fluid.get(
                        "first_fluid_shell_luminosity_msun_kms2_per_myr"
                    ),
                    "maximum_outward_gravothermal_luminosity_msun_kms2_per_myr": fluid.get(
                        "maximum_outward_luminosity_msun_kms2_per_myr"
                    ),
                    "maximum_absolute_gravothermal_luminosity_msun_kms2_per_myr": fluid.get(
                        "maximum_absolute_luminosity_msun_kms2_per_myr"
                    ),
                    "radius_of_maximum_absolute_luminosity_pc": fluid.get(
                        "radius_of_maximum_absolute_luminosity_pc"
                    ),
                }
            )
        rows.append(row)

    args.outdir.mkdir(parents=True, exist_ok=True)
    columns = sorted({key for row in rows for key in row})
    csv_path = args.outdir / "scan20_results.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)

    prep_complete = prep_counts.get("INTERFACE_COMPLETE", 0)
    fluid_complete = fluid_counts.get("ALL_HALO_FLUX_COMPLETE", 0)
    geometric_admissible = sum(
        row.get("fp_geometric_admissible") is True for row in rows
    )
    status = {
        "schema": "sidm-smbh-interface-fluid-scan-20x20-aggregate-v1",
        "status": (
            "SCAN20_COMPLETE"
            if prep_complete == grid["case_count"]
            and fluid_complete == grid["case_count"]
            else "SCAN20_PARTIAL"
        ),
        "updated_at": datetime.now(timezone.utc).isoformat(),
        "grid_shape": grid["grid_shape"],
        "case_count": grid["case_count"],
        "interface_complete_count": prep_complete,
        "fluid_complete_count": fluid_complete,
        "fp_geometric_admissible_count": geometric_admissible,
        "preparation_status_counts": dict(prep_counts),
        "fluid_status_counts": dict(fluid_counts),
        "results_csv": str(csv_path),
        "scope": (
            "Exact black-hole remap, interface, and fluid luminosity. "
            "No high-resolution FP capture current is inferred."
        ),
    }
    status_path = args.outdir / "scan20_summary.json"
    status_path.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")
    print(json.dumps(status, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
