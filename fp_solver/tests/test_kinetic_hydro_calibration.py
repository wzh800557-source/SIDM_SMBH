"""Synthetic recovery test for the branch-specific Kn transition fit."""

import json
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from calibrate_kinetic_hydro_bridge import aggregate, bridge_rate  # noqa: E402


def main():
    hydro = 3.141592653589793
    ballistic = 0.12
    kn_star = 0.42
    exponent = 1.35
    kn = np.geomspace(0.01, 40.0, 11)
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        for label, grid in (
            ("medium", (40, 56, 32)),
            ("fine", (56, 72, 44)),
        ):
            for index, value in enumerate(kn):
                mdot = float(bridge_rate(value, hydro, ballistic, kn_star, exponent))
                item = {
                    "schema": "non-orbit-averaged-radial-bgk-v1",
                    "status": "KINETIC_RUN_PASS",
                    "resolution_label": label,
                    "grid": {
                        "radial_bins": grid[0],
                        "radial_velocity_bins": grid[1],
                        "tangential_bins": grid[2],
                    },
                    "knudsen_infinity": float(value),
                    "stationarity": {"status": "PASS", "late_mean_mdot": mdot},
                    "hydrodynamic_reference": {"mdot": hydro},
                    "ballistic_reference": {"mdot": ballistic},
                }
                (root / f"{label}_{index}.json").write_text(json.dumps(item))
        result = aggregate(root)
        json.dumps(result)
    assert result["status"] == "CALIBRATION_ACCEPTED"
    assert abs(result["adopted_parameters"]["Kn_star"] / kn_star - 1.0) < 2.0e-5
    assert abs(result["adopted_parameters"]["p"] / exponent - 1.0) < 2.0e-5
    print("PASS: synthetic kinetic-hydrodynamic transition calibration")


if __name__ == "__main__":
    main()
