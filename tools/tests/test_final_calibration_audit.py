"""Regression test for the integrated final-calibration ledger."""

import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from build_final_calibration_audit import build_audit  # noqa: E402


def main():
    result = build_audit(ROOT)
    assert result["status"] == (
        "MASS_CURRENT_AND_BRANCH_SELECTION_COMPLETE_ENERGY_CLOSURE_OPEN"
    )
    accepted = result["accepted_results"]
    assert math.isclose(
        accepted["captured_mass_current"]["direct_mdot_msun_per_myr"],
        14.06080600205032,
        rel_tol=1.0e-12,
    )
    assert accepted["cross_regime_calibration"]["fiducial_selected_branch"] == (
        "orbit-resolved FP"
    )
    assert not accepted["cross_regime_calibration"][
        "bondi_authorized_for_fiducial_profile"
    ]
    assert result["open_acceptance_gates"]["capture_weighted_energy_current"][
        "status"
    ] == "NOT_ACCEPTED"
    assert result["open_acceptance_gates"]["absolute_two_way_evolution"][
        "status"
    ] == "NOT_AUTHORIZED"
    print("PASS: integrated final-calibration acceptance ledger")


if __name__ == "__main__":
    main()
