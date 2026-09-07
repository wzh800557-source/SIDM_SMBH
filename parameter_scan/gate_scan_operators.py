#!/usr/bin/env python3
"""Audit all scan operators and repair roundoff-only status failures."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from run_scan_operator import (
    SLOTS,
    complete_payload,
    domain_exclusion_payload,
    estimator_domain_exclusion,
    validate_products,
    write_status,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()

    grid = json.loads(args.grid.read_text())
    failures = []
    repaired = []
    statuses = []
    admissibility_statuses = []
    for case in grid["cases"]:
        case_index = int(case["index"])
        case_root = args.root / "cases" / case["tag"]
        prep = json.loads((case_root / "PREP_STATUS.json").read_text())
        case_failures = []
        domain_exclusions = []
        case_statuses = []
        for energy_bins, seed in SLOTS:
            outdir = case_root / "operators" / f"e{energy_bins}_j257_s{seed}"
            status_path = outdir / "STATUS.json"
            if not status_path.is_file():
                case_failures.append(
                    {
                        "case_tag": case["tag"],
                        "energy_bins": energy_bins,
                        "seed": seed,
                        "reason": "MISSING_STATUS",
                    }
                )
                continue
            status = json.loads(status_path.read_text())
            if status.get("status") == "FP_OPERATOR_SKIPPED":
                if prep.get("status") == "READY_FOR_FP":
                    case_failures.append(
                        {
                            "case_tag": case["tag"],
                            "energy_bins": energy_bins,
                            "seed": seed,
                            "reason": "READY_CASE_WAS_SKIPPED",
                        }
                    )
                else:
                    statuses.append("FP_OPERATOR_SKIPPED")
                    case_statuses.append("FP_OPERATOR_SKIPPED")
                continue
            if status.get("status") == "FP_OPERATOR_COMPLETE":
                statuses.append("FP_OPERATOR_COMPLETE")
                case_statuses.append("FP_OPERATOR_COMPLETE")
                continue

            run_path = outdir / "run.json"
            operator_path = outdir / "operator.npz"
            if run_path.is_file() and operator_path.is_file():
                run = json.loads(run_path.read_text())
                exclusion = estimator_domain_exclusion(run)
                if exclusion is not None:
                    write_status(
                        status_path,
                        domain_exclusion_payload(
                            {
                                "schema": "sidm-smbh-scan-operator-v1",
                                "case_index": case_index,
                                "case_tag": case["tag"],
                                "energy_bins": energy_bins,
                                "angular_bins": 257,
                                "seed": seed,
                                "pairs_per_bin": run.get("pairs_per_bin"),
                            },
                            outdir,
                            run,
                            exclusion,
                        ),
                    )
                    statuses.append("FP_OPERATOR_DOMAIN_UNRESOLVED")
                    case_statuses.append("FP_OPERATOR_DOMAIN_UNRESOLVED")
                    domain_exclusions.append(
                        {
                            "energy_bins": energy_bins,
                            "seed": seed,
                            **exclusion,
                        }
                    )
                    continue

            if not run_path.is_file() or not operator_path.is_file():
                case_failures.append(
                    {
                        "case_tag": case["tag"],
                        "energy_bins": energy_bins,
                        "seed": seed,
                        "reason": status.get("reason", "MISSING_PRODUCTS"),
                    }
                )
                continue
            try:
                run = json.loads(run_path.read_text())
                mode = validate_products(outdir, energy_bins, run)
                base = {
                    "schema": "sidm-smbh-scan-operator-v1",
                    "case_index": case_index,
                    "case_tag": case["tag"],
                    "energy_bins": energy_bins,
                    "angular_bins": 257,
                    "seed": seed,
                    "pairs_per_bin": run.get("pairs_per_bin"),
                }
                write_status(
                    status_path,
                    complete_payload(base, outdir, run, mode),
                )
                statuses.append("FP_OPERATOR_COMPLETE")
                case_statuses.append("FP_OPERATOR_COMPLETE")
                repaired.append(
                    {
                        "case_tag": case["tag"],
                        "energy_bins": energy_bins,
                        "seed": seed,
                        "acceptance_mode": mode,
                    }
                )
            except Exception as exc:
                case_failures.append(
                    {
                        "case_tag": case["tag"],
                        "energy_bins": energy_bins,
                        "seed": seed,
                        "reason": f"{type(exc).__name__}: {exc}",
                    }
                )

        if case_failures:
            failures.extend(case_failures)
            admissibility = {
                "schema": "sidm-smbh-scan-fp-admissibility-v1",
                "status": "FP_ADMISSIBILITY_FAILED",
                "case_index": case_index,
                "case_tag": case["tag"],
                "reason": "UNCLASSIFIED_OPERATOR_FAILURE",
                "operator_statuses": case_statuses,
                "failures": case_failures,
            }
        elif prep.get("status") != "READY_FOR_FP":
            admissibility = {
                "schema": "sidm-smbh-scan-fp-admissibility-v1",
                "status": "FP_NOT_APPLICABLE",
                "case_index": case_index,
                "case_tag": case["tag"],
                "reason": prep.get("status"),
                "operator_statuses": case_statuses,
            }
        elif domain_exclusions:
            admissibility = {
                "schema": "sidm-smbh-scan-fp-admissibility-v1",
                "status": "FP_DOMAIN_NOT_REPRESENTABLE_AFTER_OPERATOR",
                "case_index": case_index,
                "case_tag": case["tag"],
                "reason": "CAPTURED_ORBITS_EXCEED_KEPLER_DOMAIN",
                "gate_limit": 0.01,
                "maximum_unresolved_kepler_domain_fraction": max(
                    item["unresolved_kepler_domain_fraction_of_candidate"]
                    for item in domain_exclusions
                ),
                "operator_statuses": case_statuses,
                "domain_exclusions": domain_exclusions,
            }
        else:
            admissibility = {
                "schema": "sidm-smbh-scan-fp-admissibility-v1",
                "status": "FP_ADMISSIBLE",
                "case_index": case_index,
                "case_tag": case["tag"],
                "reason": "ALL_OPERATOR_GATES_PASS",
                "operator_statuses": case_statuses,
            }
        admissibility["updated_at"] = datetime.now(timezone.utc).isoformat()
        (case_root / "FP_ADMISSIBILITY.json").write_text(
            json.dumps(admissibility, indent=2, sort_keys=True) + "\n"
        )
        admissibility_statuses.append(admissibility["status"])

    expected = len(grid["cases"]) * len(SLOTS)
    passed = len(statuses) == expected and not failures
    payload = {
        "schema": "sidm-smbh-scan-operator-gate-v1",
        "status": "OPERATOR_GATE_PASS" if passed else "OPERATOR_GATE_FAIL",
        "updated_at": datetime.now(timezone.utc).isoformat(),
        "expected_operator_count": expected,
        "classified_operator_count": len(statuses),
        "usable_operator_count": statuses.count("FP_OPERATOR_COMPLETE"),
        "status_counts": dict(Counter(statuses)),
        "fp_admissibility_status_counts": dict(Counter(admissibility_statuses)),
        "roundoff_only_repairs": repaired,
        "failures": failures,
    }
    output = args.root / "OPERATOR_GATE.json"
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
