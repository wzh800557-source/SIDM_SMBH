#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from diagnose_ej_capture_tail import analyze_case  # noqa: E402
from solve_ej_steady_state import sha256_file  # noqa: E402


def main() -> int:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        operator = root / "operator.npz"
        steady_json = root / "steady.json"
        steady_csv = root / "steady.csv"
        metadata = {
            "schema": "finite-angle-ej-jump-operator-v1",
            "state_count": 4,
            "energy_bins": 2,
            "angular_bins": 1,
            "angular_grid": "synthetic",
            "reservoir_radius_pc": 1.0,
        }
        arrays = {
            "source1": np.array([1, 0, 3, 2], dtype=np.int64),
            "source2": np.array([1, 1, 3, 3], dtype=np.int64),
            "pre": np.array([1, 0, 3, 2], dtype=np.int64),
            "post": np.array([0, -1, 2, -1], dtype=np.int64),
            "rate_direct_msun_per_myr": np.array([2.0, 4.0, 3.0, 6.0]),
            "rate_immediate_msun_per_myr": np.array([2.0, 4.0, 3.0, 6.0]),
            "rate_direct_binding_msun_kms2_per_myr": np.array([
                0.0, 40.0, 0.0, 120.0,
            ]),
            "rate_immediate_binding_msun_kms2_per_myr": np.array([
                0.0, 40.0, 0.0, 120.0,
            ]),
            "sample_count": np.array([100, 100, 100, 100], dtype=np.int64),
            "x_edges": np.array([0.1, 1.0, 10.0]),
            "j_over_jlc_edges": np.array([1.0, 2.0]),
        }
        np.savez_compressed(
            operator,
            **arrays,
            metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)),
        )
        steady = {
            "schema": "finite-angle-ej-steady-state-v2",
            "status": "EJ_STEADY_SMOKE_PASS",
            "operator_sha256": sha256_file(operator),
            "operator_metadata": metadata,
            "models": {
                model: {
                    "steady_capture_mdot_msun_per_myr": 5.0,
                    "steady_capture_binding_energy_current_msun_kms2_per_myr": 80.0,
                }
                for model in ("direct", "immediate")
            },
        }
        steady_json.write_text(json.dumps(steady, indent=2, sort_keys=True) + "\n")
        with steady_csv.open("w", newline="") as stream:
            writer = csv.DictWriter(
                stream,
                fieldnames=["state", "occupation_direct", "occupation_immediate"],
            )
            writer.writeheader()
            for state, occupation in enumerate([0.5, 1.0, 0.5, 1.0]):
                writer.writerow({
                    "state": state,
                    "occupation_direct": occupation,
                    "occupation_immediate": occupation,
                })

        result = analyze_case("synthetic", operator, steady_json, steady_csv)
        for model in ("direct", "immediate"):
            diagnostics = result["models"][model]
            assert diagnostics["mass_roundtrip_relative_error"] < 1.0e-12
            assert diagnostics["binding_energy_roundtrip_relative_error"] < 1.0e-12
            assert abs(diagnostics["highest_energy_bin_mass_fraction"] - 0.6) < 1.0e-12
            assert abs(
                diagnostics["highest_energy_bin_binding_energy_fraction"] - 0.75
            ) < 1.0e-12
            assert all(result["gates"][model].values())
    print("PASS: saved-solution capture-tail diagnostic")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
