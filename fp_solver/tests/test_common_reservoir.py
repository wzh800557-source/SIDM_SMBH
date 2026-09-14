#!/usr/bin/env python3
"""Regression test for distributed DF normalization and handoff scans."""

from __future__ import annotations

import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path


PACKAGE = Path(__file__).resolve().parents[1]
PROFILE = Path(__file__).with_name("deep_profile_proxy.txt")
MODEL = Path(__file__).with_name("model_template.in")


def run(*args: object) -> None:
    subprocess.run([sys.executable, *map(str, args)], check=True, stdout=subprocess.DEVNULL)


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="gnc-common-reservoir-") as tmp:
        root = Path(tmp)
        bridged = root / "bridged_profile.txt"
        bridge_json = root / "bridge.json"
        run(
            PACKAGE / "hydrostatic_bridge.py",
            PROFILE,
            "--out-profile",
            bridged,
            "--out-json",
            bridge_json,
            "--mbh",
            4.0e6,
            "--sigma-over-m",
            100.0,
            "--mode",
            "broken-at-rh",
        )
        bridge = json.loads(bridge_json.read_text())
        anchor = float(bridge["r_in_pc"])
        boundaries = {
            "outer": anchor,
            "inner": 0.8 * anchor,
            "outer_repeat": anchor,
        }

        manifests = {}
        diagnostics = {}
        for name, boundary in boundaries.items():
            norm = root / f"norm_{name}"
            rundir = root / f"run_{name}"
            run(
                PACKAGE / "normalize_gnc_df.py",
                bridged,
                "--bridge-json",
                bridge_json,
                "--outdir",
                norm,
                "--mbh",
                4.0e6,
                "--r-boundary",
                boundary,
                "--normalization-radius",
                anchor,
                "--proposal-boundary-radius",
                anchor,
                "--ranks",
                3,
                "--samples-per-rank",
                300,
                "--grid-bins",
                48,
            )
            run(
                PACKAGE / "prepare_absolute_run.py",
                "--base-model",
                MODEL,
                "--normalized-dir",
                norm,
                "--run-dir",
                rundir,
                "--ranks",
                3,
                "--gx-bins",
                48,
                "--dc-bins",
                72,
                "--dt-tnr",
                0.1,
                "--total-tnr",
                0.2,
                "--gnc-seed",
                7300 + len(manifests),
            )
            diagnostics[name] = json.loads((norm / "df_normalization.json").read_text())
            manifests[name] = json.loads((rundir / "run_manifest.json").read_text())

        for diag in diagnostics.values():
            assert diag["status"] == "PASS"
            assert math.isclose(
                float(diag["importance_weight_global_mean"]),
                1.0,
                rel_tol=0.0,
                abs_tol=5.0e-13,
            )
            rank_means = diag["importance_weight_rank_means"]
            assert max(abs(float(x) - 1.0) for x in rank_means) > 1.0e-4

        assert (
            manifests["outer"]["common_reservoir_fingerprint"]
            == manifests["inner"]["common_reservoir_fingerprint"]
            == manifests["outer_repeat"]["common_reservoir_fingerprint"]
        )
        assert manifests["outer"]["gnc_seed"] != manifests["inner"]["gnc_seed"]
        for name, manifest in manifests.items():
            model = (root / f"run_{name}" / "model.in").read_text()
            assert f"SEED_VALUE\t\t\t\t\t\t= {manifest['gnc_seed']}" in model
            assert "SAME_INI_SEED\t\t\t\t\t= 1" in model
            assert "SAME_EVL_SEED\t\t\t\t\t= 1" in model
        assert (
            diagnostics["outer"]["df_table_sha256"]
            == diagnostics["outer_repeat"]["df_table_sha256"]
        )
        assert manifests["outer"]["r_in_pc"] != manifests["inner"]["r_in_pc"]
        assert math.isclose(
            manifests["outer"]["normalization_radius_pc"],
            manifests["inner"]["normalization_radius_pc"],
        )
    print("PASS: global MPI normalization and common-reservoir fingerprint")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
