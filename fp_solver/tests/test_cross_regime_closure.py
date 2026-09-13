"""Analytic and regression tests for the cross-regime closure branches."""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cross_regime_closure import (  # noqa: E402
    ConductiveSpikeCalibration,
    bondi_lambda,
    bondi_rate_msun_per_myr,
    current_scale_msun_per_myr,
    diagnostic_mass_coefficient,
    radial_collisionality_gate,
    select_mass_current_branch,
)
from calibrate_cross_regime import (  # noqa: E402
    fiducial_audit,
    fit_generalized_shape,
    fit_published_form,
    leave_one_out_published_form,
    read_digitized,
)


def main():
    assert math.isclose(bondi_lambda(1.0), math.exp(1.5) / 4.0, rel_tol=1e-14)
    assert math.isclose(bondi_lambda(5.0 / 3.0), 0.25, rel_tol=1e-14)
    assert bondi_lambda(1.4) > 0.25

    mdot = bondi_rate_msun_per_myr(4.0e6, 46.855, 102.77, 5.0 / 3.0)
    assert 4.0e4 < mdot < 4.2e4
    scale = current_scale_msun_per_myr(1.295, 46.855, 79.605)
    coefficient = diagnostic_mass_coefficient(14.060806, 1.295, 46.855, 79.605)
    assert math.isclose(coefficient, 14.060806 / scale, rel_tol=1e-14)

    calibration = ConductiveSpikeCalibration()
    transition = calibration.sigma_transition_cm2_g
    assert math.isclose(transition, math.sqrt(0.27 / 0.056), rel_tol=1e-14)
    assert math.isclose(calibration.normalized_excess(transition), 1.0, rel_tol=1e-14)
    low = calibration.excess_rate_msun_per_yr(1.0e-6)
    high = calibration.excess_rate_msun_per_yr(1.0e6)
    assert math.isclose(low / 1.0e-6, 1.0 / calibration.B, rel_tol=1e-5)
    assert math.isclose(high * 1.0e6, 1.0 / calibration.A, rel_tol=1e-5)
    assert calibration.excess_rate_msun_per_yr(transition / 2.0) < calibration.maximum_excess_msun_per_yr
    assert calibration.excess_rate_msun_per_yr(transition * 2.0) < calibration.maximum_excess_msun_per_yr

    radius = np.geomspace(1.0e-6, 10.0, 200)
    orbit_memory = np.geomspace(1.0e-8, 20.0, 200)
    gate = radial_collisionality_gate(radius, orbit_memory, radius[0], radius[-1])
    assert gate["status"] == "ORBIT_MEMORY_PRESENT"
    assert select_mass_current_branch(gate, radial_flow_family="maintained-radial-inflow") == "fp"

    transition_gate = radial_collisionality_gate(
        radius, np.full_like(radius, 2.0), radius[0], radius[-1]
    )
    assert transition_gate["status"] == "TRANSITIONAL_COLLISIONALITY"
    assert select_mass_current_branch(
        transition_gate, radial_flow_family="hydrostatic-conductive"
    ) == "radial-kinetic"

    hydro_gate = radial_collisionality_gate(
        radius, np.full_like(radius, 100.0), radius[0], radius[-1]
    )
    assert hydro_gate["hydrodynamic_collisionality_pass"]
    assert select_mass_current_branch(
        hydro_gate, radial_flow_family="hydrostatic-conductive"
    ) == "hydrodynamic-candidate"
    assert select_mass_current_branch(
        hydro_gate,
        radial_flow_family="hydrostatic-conductive",
        flow_conditions_verified=True,
    ) == "conductive-spike"
    assert select_mass_current_branch(
        hydro_gate,
        radial_flow_family="maintained-radial-inflow",
        flow_conditions_verified=True,
    ) == "bondi"

    root = Path(__file__).resolve().parents[2]
    s, measured = read_digitized(
        root / "data/calibration/sabarish_2025_imfp_digitized.csv"
    )
    refit = fit_published_form(s, measured)
    generalized = fit_generalized_shape(s, measured)
    leave_one_out = leave_one_out_published_form(s, measured)
    assert refit["optimizer_success"] and refit["rms_msun_per_yr"] < 0.13
    assert abs(refit["A"] - 0.056) < 0.003
    assert abs(refit["B"] - 0.27) < 0.02
    assert 2.18 < refit["sigma_transition_cm2_g"] < 2.21
    assert 0.8 < generalized["shape_exponent_p"] < 1.0
    assert leave_one_out["number_of_refits"] == measured.size
    assert leave_one_out["ranges"]["sigma_transition_cm2_g"]["minimum"] > 2.0
    assert leave_one_out["ranges"]["sigma_transition_cm2_g"]["maximum"] < 2.4

    fiducial = fiducial_audit(
        root / "data/production/bridge/bridged_profile.txt",
        root / "data/production/bridge/bridge.json",
        root
        / "data/reference/production/aggregate_v6_mass_reaudit_boundaryfix/fluid_mass_closure.json",
    )
    assert fiducial["radial_collisionality_gate"]["status"] == "ORBIT_MEMORY_PRESENT"
    assert fiducial["radial_collisionality_gate"]["N_orb_at_capture"] < 1.0e-7
    assert 2800.0 < fiducial["nominal_bondi_benchmark_over_measured_fp"] < 3000.0
    assert not fiducial["radial_collisionality_gate"]["hydrodynamic_collisionality_pass"]
    print(
        "PASS: Bondi limits, conductive-spike calibration, fiducial Yukawa "
        "regime gate, and state-dependent solver selection"
    )


if __name__ == "__main__":
    main()
