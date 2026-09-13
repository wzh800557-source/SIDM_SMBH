#!/usr/bin/env python3
"""Calibrate the available cross-regime benchmarks and audit the fiducial halo."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from cross_regime_closure import (
    ConductiveSpikeCalibration,
    bondi_rate_msun_per_myr,
    bondi_radius_pc,
    current_scale_msun_per_myr,
    diagnostic_mass_coefficient,
    interpolate_log_profile,
    radial_collisionality_gate,
    sound_speed_from_sigma_1d,
)


def read_digitized(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows = list(csv.DictReader(path.open()))
    if len(rows) < 5:
        raise ValueError("calibration table is too short")
    s = np.asarray([float(row["sigma_over_m_cm2_g"]) for row in rows])
    rate = np.asarray([float(row["accretion_rate_msun_per_yr"]) for row in rows])
    if np.any(s <= 0.0) or np.any(rate <= 0.0):
        raise ValueError("calibration values must be positive")
    return s, rate


def fit_published_form(s: np.ndarray, rate: np.ndarray) -> dict:
    def model(par):
        baseline, a, b = par
        return baseline + 1.0 / (a * s + b / s)

    fit = least_squares(
        lambda par: model(par) - rate,
        np.array([1.7, 0.056, 0.27]),
        bounds=(0.0, np.inf),
        xtol=1.0e-14,
        ftol=1.0e-14,
        gtol=1.0e-14,
    )
    baseline, a, b = map(float, fit.x)
    residual = model(fit.x) - rate
    rms = float(np.sqrt(np.mean(residual**2)))
    return {
        "baseline_msun_per_yr": baseline,
        "A": a,
        "B": b,
        "sigma_transition_cm2_g": math.sqrt(b / a),
        "maximum_excess_msun_per_yr": 1.0 / (2.0 * math.sqrt(a * b)),
        "rms_msun_per_yr": rms,
        "maximum_absolute_residual_msun_per_yr": float(np.max(np.abs(residual))),
        "optimizer_success": bool(fit.success),
    }


def fit_generalized_shape(s: np.ndarray, rate: np.ndarray) -> dict:
    def model(par):
        baseline, amplitude, transition, p = par
        x = s / transition
        return baseline + amplitude * 2.0 * x**p / (1.0 + x ** (2.0 * p))

    fit = least_squares(
        lambda par: model(par) - rate,
        np.array([1.7, 4.0, 2.2, 1.0]),
        bounds=(np.array([0.0, 0.0, 1.0e-4, 0.1]), np.array([np.inf, np.inf, np.inf, 4.0])),
        xtol=1.0e-14,
        ftol=1.0e-14,
        gtol=1.0e-14,
    )
    baseline, amplitude, transition, p = map(float, fit.x)
    residual = model(fit.x) - rate
    return {
        "baseline_msun_per_yr": baseline,
        "maximum_excess_msun_per_yr": amplitude,
        "sigma_transition_cm2_g": transition,
        "shape_exponent_p": p,
        "rms_msun_per_yr": float(np.sqrt(np.mean(residual**2))),
        "optimizer_success": bool(fit.success),
    }


def leave_one_out_published_form(s: np.ndarray, rate: np.ndarray) -> dict:
    """Return deterministic leave-one-out ranges for the harmonic refit.

    These ranges measure sensitivity to individual plotted points.  They are
    not statistical confidence intervals because the source figure does not
    provide a covariance model for the digitized rates.
    """

    fits = []
    for omitted in range(s.size):
        keep = np.arange(s.size) != omitted
        fits.append(fit_published_form(s[keep], rate[keep]))
    keys = (
        "baseline_msun_per_yr",
        "A",
        "B",
        "sigma_transition_cm2_g",
        "maximum_excess_msun_per_yr",
    )
    return {
        "interpretation": (
            "Sensitivity to omission of one vector point; not a statistical "
            "confidence interval."
        ),
        "number_of_refits": len(fits),
        "ranges": {
            key: {
                "minimum": float(min(fit[key] for fit in fits)),
                "maximum": float(max(fit[key] for fit in fits)),
            }
            for key in keys
        },
    }


def fiducial_audit(
    bridge_profile: Path,
    bridge_json: Path,
    mass_closure_json: Path,
) -> dict:
    profile = np.loadtxt(bridge_profile, comments="#", ndmin=2)
    radius, density, sigma, n_orb = (
        np.asarray(profile[:, index], dtype=float) for index in (0, 1, 2, 4)
    )
    bridge = json.loads(bridge_json.read_text())
    closure = json.loads(mass_closure_json.read_text())
    black_hole_mass = float(bridge["mbh_msun"])
    capture_radius = float(bridge["r_isco_pc"])
    boundary_radius = float(bridge["r_in_pc"])
    boundary_density = float(bridge["rho_in_msun_pc3"])
    boundary_sigma = float(bridge["sigma_in_kms"])
    gamma = 5.0 / 3.0
    sound_speed = sound_speed_from_sigma_1d(boundary_sigma, gamma)
    nominal_r_bondi = bondi_radius_pc(black_hole_mass, sound_speed)
    nominal_bondi_rate = bondi_rate_msun_per_myr(
        black_hole_mass, boundary_density, sound_speed, gamma
    )
    measured_rate = float(closure["measured_mdot_msun_per_myr"])
    gate = radial_collisionality_gate(
        radius,
        n_orb,
        capture_radius,
        nominal_r_bondi,
    )
    n_boundary = interpolate_log_profile(radius, n_orb, boundary_radius)
    n_bondi = interpolate_log_profile(radius, n_orb, nominal_r_bondi)
    sigma0 = float(bridge["sigma_over_m_cm2_g"])
    sigma_critical_fixed_profile = sigma0 / gate["N_orb_at_capture"]
    sigma_continuum_fixed_profile = (
        sigma0
        * gate["hydrodynamic_n_orb_threshold"]
        / gate["N_orb_at_capture"]
    )
    scale = current_scale_msun_per_myr(
        boundary_radius, boundary_density, boundary_sigma
    )
    return {
        "black_hole_mass_msun": black_hole_mass,
        "gamma": gamma,
        "boundary_radius_pc": boundary_radius,
        "boundary_density_msun_pc3": boundary_density,
        "boundary_sigma_1d_kms": boundary_sigma,
        "sound_speed_infinity_kms": sound_speed,
        "nominal_bondi_radius_pc": nominal_r_bondi,
        "bondi_reservoir_identification": (
            "Diagnostic substitution rho_infinity=rho(r_in) and "
            "sigma_infinity=sigma_1d(r_in); the saved hydrostatic snapshot does "
            "not provide an asymptotic radial-flow reservoir."
        ),
        "measured_fp_mdot_msun_per_myr": measured_rate,
        "measured_fp_C_M": diagnostic_mass_coefficient(
            measured_rate, boundary_radius, boundary_density, boundary_sigma
        ),
        "nominal_adiabatic_bondi_benchmark_mdot_msun_per_myr": nominal_bondi_rate,
        "nominal_bondi_benchmark_over_measured_fp": nominal_bondi_rate
        / measured_rate,
        "nominal_bondi_C_M_in_boundary_normalization": nominal_bondi_rate / scale,
        "N_orb_at_r_in": n_boundary,
        "N_orb_at_r_B": n_bondi,
        "radial_collisionality_gate": gate,
        "sigma0_over_m_for_N_ISCO_equal_one_fixed_profile_cm2_g": sigma_critical_fixed_profile,
        "sigma0_over_m_for_continuum_sentinel_at_ISCO_fixed_profile_cm2_g": (
            sigma_continuum_fixed_profile
        ),
        "fixed_profile_scaling_warning": (
            "This threshold scales the saved profile linearly in sigma0/m. "
            "It is not a self-consistent high-cross-section halo solution."
        ),
    }


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--digitized-csv",
        type=Path,
        default=root / "data/calibration/sabarish_2025_imfp_digitized.csv",
    )
    p.add_argument(
        "--digitized-metadata",
        type=Path,
        default=root / "data/calibration/sabarish_2025_imfp_digitized.json",
    )
    p.add_argument(
        "--bridge-profile",
        type=Path,
        default=root / "data/production/bridge/bridged_profile.txt",
    )
    p.add_argument(
        "--bridge-json",
        type=Path,
        default=root / "data/production/bridge/bridge.json",
    )
    p.add_argument(
        "--mass-closure-json",
        type=Path,
        default=(
            root
            / "data/reference/production/aggregate_v6_mass_reaudit_boundaryfix/fluid_mass_closure.json"
        ),
    )
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()

    s, rate = read_digitized(args.digitized_csv)
    metadata = json.loads(args.digitized_metadata.read_text())
    published = ConductiveSpikeCalibration(
        A=float(metadata["published_A"]),
        B=float(metadata["published_B"]),
        baseline_msun_per_yr=float(metadata["collisionless_baseline_msun_per_yr"]),
    )
    same_form = fit_published_form(s, rate)
    generalized = fit_generalized_shape(s, rate)
    leave_one_out = leave_one_out_published_form(s, rate)
    result = {
        "schema": "sidm-smbh-cross-regime-calibration-v2",
        "status": "AVAILABLE_BRANCHES_CALIBRATED",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "conductive_spike": {
            "metadata": metadata,
            "reference_parameters": published.as_dict(),
            "digitized_same_form_refit": same_form,
            "leave_one_out_refit_range": leave_one_out,
            "digitized_generalized_shape_check": generalized,
            "adopted_shape_exponent_p": 1.0,
            "adoption_reason": (
                "The harmonic-conductivity asymptotes fix p=1. The free-p fit is "
                "retained only as a diagnostic of the measured turnover shape."
            ),
        },
        "fiducial_yukawa_halo": fiducial_audit(
            args.bridge_profile, args.bridge_json, args.mass_closure_json
        ),
        "closure_decision": {
            "universal_monotonic_fp_to_bondi_bridge": "REJECTED",
            "reason": (
                "A hydrostatic conductive spike and a maintained radial Bondi "
                "flow have different high-collisionality limits. Collisionality "
                "alone does not specify the outer reservoir or radial flow state."
            ),
            "implemented_branches": [
                "measured orbit-resolved FP current",
                "published constant-cross-section conductive-spike benchmark",
                "exact polytropic Bondi endpoint",
            ],
            "transition_solver": (
                "Use non-orbit-averaged radial kinetics when N_orb crosses unity "
                "inside the inflow path."
            ),
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
