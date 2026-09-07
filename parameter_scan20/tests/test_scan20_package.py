#!/usr/bin/env python3
"""Local structural tests for the 20 by 20 cluster scan package."""

from __future__ import annotations

import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parents[1]


def main() -> int:
    grid = json.loads((HERE / "scan_grid.json").read_text())
    assert grid["grid_shape"] == [3, 20, 20]
    assert grid["case_count"] == 1200
    assert len(grid["cases"]) == 1200
    assert len({case["tag"] for case in grid["cases"]}) == 1200
    assert len({case["halo_mass_factor"] for case in grid["cases"]}) == 20
    assert len({case["black_hole_mass_factor"] for case in grid["cases"]}) == 20
    assert {case["sigma0_over_m_cm2_g"] for case in grid["cases"]} == {
        10.0,
        100.0,
        1000.0,
    }

    halo = np.array(grid["halo_mass_factors"])
    black_hole = np.array(grid["black_hole_mass_factors"])
    assert np.all(np.diff(halo) > 0.0)
    assert np.all(np.diff(black_hole) > 0.0)
    assert np.isclose(halo[0], 0.1) and np.isclose(halo[-1], 10.0)
    assert np.isclose(black_hole[0], 0.1) and np.isclose(black_hole[-1], 10.0)

    for case in grid["cases"]:
        expected = case["black_hole_mass_msun"] / case[
            "represented_halo_mass_msun"
        ]
        assert np.isclose(
            expected, case["black_hole_to_represented_halo_mass"]
        )

    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        csv_path = root / "synthetic.csv"
        fields = (
            "represented_halo_mass_msun",
            "black_hole_mass_msun",
            "sigma0_over_m_cm2_g",
            "r_in_over_r_h",
            "fp_geometric_admissible",
            "first_shell_gravothermal_luminosity_msun_kms2_per_myr",
            "maximum_absolute_gravothermal_luminosity_msun_kms2_per_myr",
        )
        with csv_path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for case in grid["cases"]:
                hf = case["halo_mass_factor"]
                bf = case["black_hole_mass_factor"]
                sf = case["sigma0_over_m_cm2_g"] / 100.0
                ratio = 0.12 * hf**0.68 * bf**0.05 * sf**-0.96
                writer.writerow(
                    {
                        "represented_halo_mass_msun": case[
                            "represented_halo_mass_msun"
                        ],
                        "black_hole_mass_msun": case["black_hole_mass_msun"],
                        "sigma0_over_m_cm2_g": case["sigma0_over_m_cm2_g"],
                        "r_in_over_r_h": ratio,
                        "fp_geometric_admissible": ratio < 0.5,
                        "first_shell_gravothermal_luminosity_msun_kms2_per_myr": (
                            5.0e8 * hf**1.2 * bf**0.4 * sf**0.5
                        ),
                        "maximum_absolute_gravothermal_luminosity_msun_kms2_per_myr": (
                            1.4e9 * hf**1.1 * bf**0.6 * sf**0.3
                        ),
                    }
                )
        if importlib.util.find_spec("matplotlib") is not None:
            subprocess.run(
                [
                    sys.executable,
                    str(HERE / "plot_scan20_heatmaps.py"),
                    "--results-csv",
                    str(csv_path),
                    "--outdir",
                    str(root),
                ],
                check=True,
            )
            for name in (
                "fig_scan20_interface_heatmap.png",
                "fig_scan20_interface_heatmap.pdf",
                "fig_scan20_fluid_luminosity_heatmap.png",
                "fig_scan20_fluid_luminosity_heatmap.pdf",
            ):
                assert (root / name).stat().st_size > 1000
        else:
            print("SCAN20_PLOT_TEST_SKIPPED_NO_LOCAL_MATPLOTLIB")

    print("SCAN20_PACKAGE_TESTS_PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
