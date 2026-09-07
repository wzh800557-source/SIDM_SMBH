#!/usr/bin/env python3
"""Aggregate the 3x3x3 interface and FP flux scan and fit power laws."""

from __future__ import annotations

import argparse
import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


def flatten_case(
    case: dict,
    prep: dict,
    admissibility: dict | None,
    steady_status: dict | None,
    flux: dict | None,
    all_halo_status: dict | None,
    all_halo_fluid: dict | None,
) -> dict:
    row = {
        "case_index": case["index"],
        "case_tag": case["tag"],
        "halo_mass_factor": case["halo_mass_factor"],
        "represented_halo_mass_msun": case["represented_halo_mass_msun"],
        "black_hole_mass_msun": case["black_hole_mass_msun"],
        "black_hole_to_halo_mass": case[
            "black_hole_to_represented_halo_mass"
        ],
        "sigma0_over_m_cm2_g": case["sigma0_over_m_cm2_g"],
        "prep_status": prep.get("status"),
        "prep_reason": prep.get("reason"),
        "interface_implementation_status": prep.get(
            "interface_implementation_status"
        ),
        "bridge_status": prep.get("bridge_status"),
        "r_in_pc": prep.get("r_in_pc"),
        "r_h_pc": prep.get("rh_pc"),
        "r_in_over_r_h": prep.get("r_in_over_rh"),
        "first_fluid_shell_radius_pc": prep.get("first_shell_radius_pc"),
        "N_at_first_fluid_shell": prep.get("N_at_first_fluid_shell"),
        "fp_owned_mass_fraction_of_first_shell": prep.get(
            "fp_owned_mass_fraction_of_first_shell"
        ),
        "flux_status": None,
        "fp_admissibility_status": (
            admissibility.get("status") if admissibility else None
        ),
        "fp_admissibility_reason": (
            admissibility.get("reason") if admissibility else None
        ),
        "maximum_unresolved_kepler_domain_fraction": (
            admissibility.get("maximum_unresolved_kepler_domain_fraction")
            if admissibility else None
        ),
        "steady_driver_status": (
            steady_status.get("status") if steady_status else None
        ),
        "all_halo_flux_status": (
            all_halo_status.get("status") if all_halo_status else None
        ),
    }
    if all_halo_fluid is not None:
        row.update(
            {
                "first_shell_gravothermal_luminosity_msun_kms2_per_myr": all_halo_fluid[
                    "first_fluid_shell_luminosity_msun_kms2_per_myr"
                ],
                "maximum_outward_gravothermal_luminosity_msun_kms2_per_myr": all_halo_fluid[
                    "maximum_outward_luminosity_msun_kms2_per_myr"
                ],
                "minimum_signed_gravothermal_luminosity_msun_kms2_per_myr": all_halo_fluid[
                    "minimum_signed_luminosity_msun_kms2_per_myr"
                ],
                "maximum_absolute_gravothermal_luminosity_msun_kms2_per_myr": all_halo_fluid[
                    "maximum_absolute_luminosity_msun_kms2_per_myr"
                ],
                "radius_of_maximum_absolute_gravothermal_luminosity_pc": all_halo_fluid[
                    "radius_of_maximum_absolute_luminosity_pc"
                ],
            }
        )
    if flux is None:
        return row
    direct = flux["selected_fluxes"]["direct"]
    immediate = flux["selected_fluxes"]["immediate"]
    boundary = flux["boundary_state"]
    fluid = flux["black_hole_aware_fluid_flux"]
    comparisons = flux["flux_comparisons"]
    operator = flux.get("operator_isotropic_fluxes", {})
    row.update(
        {
            "flux_status": flux["status"],
            "maximum_numerical_sensitivity": flux[
                "maximum_numerical_sensitivity"
            ],
            "rho_boundary_msun_pc3": boundary["rho_msun_pc3"],
            "sigma_boundary_kms": boundary["sigma_1d_kms"],
            "mass_flux_scale_msun_per_myr": boundary[
                "mass_flux_scale_msun_per_myr"
            ],
            "direct_steady_mdot_msun_per_myr": direct[
                "steady_mdot_msun_per_myr"
            ],
            "immediate_steady_mdot_msun_per_myr": immediate[
                "steady_mdot_msun_per_myr"
            ],
            "direct_isotropic_mdot_msun_per_myr": direct[
                "isotropic_mdot_msun_per_myr"
            ],
            "immediate_isotropic_mdot_msun_per_myr": immediate[
                "isotropic_mdot_msun_per_myr"
            ],
            "direct_binding_current_msun_kms2_per_myr": direct[
                "steady_binding_current_msun_kms2_per_myr"
            ],
            "immediate_binding_current_msun_kms2_per_myr": immediate[
                "steady_binding_current_msun_kms2_per_myr"
            ],
            "direct_thermal_sink_msun_kms2_per_myr": direct[
                "thermal_sink_current_msun_kms2_per_myr"
            ],
            "immediate_thermal_sink_msun_kms2_per_myr": immediate[
                "thermal_sink_current_msun_kms2_per_myr"
            ],
            "direct_mean_binding_energy_kms2": direct[
                "mean_capture_binding_energy_kms2"
            ],
            "direct_C_M": direct["C_M"],
            "direct_C_E_binding": direct["C_E_binding"],
            "direct_C_E_thermal_sink": direct["C_E_thermal_sink"],
            "steady_over_isotropic_direct_mass": direct[
                "steady_over_isotropic_mass_current"
            ],
            "candidate_isotropic_mdot_msun_per_myr": operator.get(
                "candidate_mass_current_msun_per_myr", {}
            ).get("mean"),
            "direct_no_rescatter_isotropic_mdot_msun_per_myr": operator.get(
                "direct_no_rescatter_mass_current_msun_per_myr", {}
            ).get("mean"),
            "candidate_isotropic_binding_current_msun_kms2_per_myr": operator.get(
                "candidate_binding_energy_current_msun_kms2_per_myr", {}
            ).get("mean"),
            "direct_no_rescatter_isotropic_binding_current_msun_kms2_per_myr": operator.get(
                "direct_no_rescatter_binding_energy_current_msun_kms2_per_myr",
                {},
            ).get("mean"),
            "bondi_benchmark_mdot_msun_per_myr": flux[
                "bondi_dimensional_supply_benchmark"
            ]["mdot_msun_per_myr"],
            "bondi_over_direct_fp": flux["bondi_dimensional_supply_benchmark"][
                "bondi_over_direct_fp_mass_current"
            ],
            "first_shell_gravothermal_luminosity_msun_kms2_per_myr": (
                all_halo_fluid or fluid
            )["first_fluid_shell_luminosity_msun_kms2_per_myr"],
            "maximum_absolute_gravothermal_luminosity_msun_kms2_per_myr": (
                all_halo_fluid or fluid
            )["maximum_absolute_luminosity_msun_kms2_per_myr"],
            "absolute_thermal_sink_over_first_shell_gravothermal_luminosity": comparisons[
                "absolute_direct_thermal_sink_over_absolute_first_shell_gravothermal_luminosity"
            ],
            "thermal_sink_inner_shell_timescale_myr": comparisons[
                "direct_thermal_sink_timescale_for_inner_shell_myr"
            ],
        }
    )
    return row


def fit_power_law(
    rows: list[dict], field: str, eligibility: str
) -> dict | None:
    selected = []
    for row in rows:
        value = row.get(field)
        if eligibility == "fp":
            eligible = row.get("flux_status") == "SCAN_FLUX_CONVERGED"
        elif eligibility == "fluid":
            eligible = row.get("all_halo_flux_status") == "ALL_HALO_FLUX_COMPLETE"
        elif eligibility == "interface":
            eligible = row.get("prep_status") in {
                "READY_FOR_FP",
                "FP_DOMAIN_NOT_REPRESENTABLE",
            }
        else:
            raise ValueError(f"unknown fit eligibility {eligibility}")
        if (
            eligible
            and isinstance(value, (int, float))
            and math.isfinite(value)
            and value != 0.0
        ):
            selected.append(row)
    if len(selected) < 6:
        return None
    reference_halo = 1.5283226426126614e10
    x = np.array(
        [
            [
                1.0,
                math.log10(row["black_hole_mass_msun"] / 4.0e6),
                math.log10(row["represented_halo_mass_msun"] / reference_halo),
                math.log10(row["sigma0_over_m_cm2_g"] / 100.0),
            ]
            for row in selected
        ]
    )
    y = np.log10(np.abs([row[field] for row in selected]))
    coefficients, _, rank, singular = np.linalg.lstsq(x, y, rcond=None)
    prediction = x @ coefficients
    residual = y - prediction
    rms_residual = float(np.sqrt(np.mean(residual**2)))
    if rms_residual <= 0.15:
        fit_quality = "GOOD"
    elif rms_residual <= 0.30:
        fit_quality = "DESCRIPTIVE_ONLY"
    else:
        fit_quality = "POOR_NON_POWER_LAW_BEHAVIOR"
    return {
        "field": field,
        "eligibility": eligibility,
        "case_count": len(selected),
        "rank": int(rank),
        "normalization_at_fiducial": float(10.0 ** coefficients[0]),
        "black_hole_mass_exponent": float(coefficients[1]),
        "halo_mass_exponent": float(coefficients[2]),
        "cross_section_exponent": float(coefficients[3]),
        "rms_residual_dex": rms_residual,
        "maximum_absolute_residual_dex": float(np.max(np.abs(residual))),
        "fit_quality": fit_quality,
        "recommended_for_interpretation": fit_quality != "POOR_NON_POWER_LAW_BEHAVIOR",
        "singular_values": [float(value) for value in singular],
        "model": (
            "|Y| = A (M_bh/4e6 Msun)^alpha "
            "(M_halo,profile/1.5283226e10 Msun)^beta "
            "[(sigma0/m)/100 cm2 g-1]^gamma"
        ),
    }


def append_radial_rows(
    radial_rows: list[dict],
    case: dict,
    prep: dict,
    flux: dict | None,
    profile_csv: Path | None,
) -> None:
    common = {
        "case_index": case["index"],
        "case_tag": case["tag"],
        "represented_halo_mass_msun": case["represented_halo_mass_msun"],
        "black_hole_mass_msun": case["black_hole_mass_msun"],
        "sigma0_over_m_cm2_g": case["sigma0_over_m_cm2_g"],
        "prep_status": prep.get("status"),
        "r_in_pc": prep.get("r_in_pc"),
        "r_h_pc": prep.get("rh_pc"),
    }
    if flux is not None:
        boundary = flux["boundary_state"]
        direct = flux["selected_fluxes"]["direct"]
        radial_rows.append(
            {
                **common,
                "radial_source": "fp_boundary",
                "radius_pc": flux["interface"]["reporting_radius_pc"],
                "rho_msun_pc3": boundary["rho_msun_pc3"],
                "sigma_1d_kms": boundary["sigma_1d_kms"],
                "enclosed_mass_msun": None,
                "gravothermal_luminosity_msun_kms2_per_myr": None,
                "fp_direct_mass_current_msun_per_myr": direct[
                    "steady_mdot_msun_per_myr"
                ],
                "fp_direct_binding_current_msun_kms2_per_myr": direct[
                    "steady_binding_current_msun_kms2_per_myr"
                ],
                "fp_direct_thermal_sink_msun_kms2_per_myr": direct[
                    "thermal_sink_current_msun_kms2_per_myr"
                ],
            }
        )
    if profile_csv is None or not profile_csv.is_file():
        return
    with profile_csv.open(newline="") as stream:
        for item in csv.DictReader(stream):
            radial_rows.append(
                {
                    **common,
                    "radial_source": "fluid_shell",
                    "radius_pc": float(item["radius_pc"]),
                    "rho_msun_pc3": float(item["rho_shell_msun_pc3"]),
                    "sigma_1d_kms": float(item["sigma_1d_kms"]),
                    "enclosed_mass_msun": float(item["enclosed_mass_msun"]),
                    "gravothermal_luminosity_msun_kms2_per_myr": float(
                        item["gravothermal_luminosity_msun_kms2_per_myr"]
                    ),
                    "fp_direct_mass_current_msun_per_myr": None,
                    "fp_direct_binding_current_msun_kms2_per_myr": None,
                    "fp_direct_thermal_sink_msun_kms2_per_myr": None,
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    grid = json.loads(args.grid.read_text())
    rows = []
    radial_rows = []
    preparation_counts = {}
    flux_counts = {}
    all_halo_flux_counts = {}
    admissibility_counts = {}
    for case in grid["cases"]:
        case_root = args.root / "cases" / case["tag"]
        prep_path = case_root / "PREP_STATUS.json"
        if not prep_path.exists():
            prep = {"status": "MISSING"}
        else:
            prep = json.loads(prep_path.read_text())
        flux_path = case_root / "flux" / "flux_budget.json"
        flux = json.loads(flux_path.read_text()) if flux_path.exists() else None
        admissibility_path = case_root / "FP_ADMISSIBILITY.json"
        admissibility = (
            json.loads(admissibility_path.read_text())
            if admissibility_path.exists()
            else None
        )
        steady_status_path = case_root / "flux" / "STATUS.json"
        steady_status = (
            json.loads(steady_status_path.read_text())
            if steady_status_path.exists()
            else None
        )
        all_halo_status_path = case_root / "all_halo_flux" / "STATUS.json"
        all_halo_status = (
            json.loads(all_halo_status_path.read_text())
            if all_halo_status_path.exists()
            else None
        )
        all_halo_flux_path = case_root / "all_halo_flux" / "fluid_snapshot_flux.json"
        all_halo_fluid = (
            json.loads(all_halo_flux_path.read_text())
            if all_halo_flux_path.exists()
            else None
        )
        row = flatten_case(
            case,
            prep,
            admissibility,
            steady_status,
            flux,
            all_halo_status,
            all_halo_fluid,
        )
        rows.append(row)
        prep_key = row["prep_status"] or "MISSING"
        flux_key = row["flux_status"] or "NOT_APPLICABLE_OR_MISSING"
        all_halo_key = row["all_halo_flux_status"] or "MISSING"
        preparation_counts[prep_key] = preparation_counts.get(prep_key, 0) + 1
        flux_counts[flux_key] = flux_counts.get(flux_key, 0) + 1
        all_halo_flux_counts[all_halo_key] = (
            all_halo_flux_counts.get(all_halo_key, 0) + 1
        )
        admissibility_key = row["fp_admissibility_status"] or "MISSING"
        admissibility_counts[admissibility_key] = (
            admissibility_counts.get(admissibility_key, 0) + 1
        )
        append_radial_rows(
            radial_rows,
            case,
            prep,
            flux,
            (
                case_root
                / "all_halo_flux"
                / "fluid_snapshot_flux_profile.csv"
            ),
        )

    args.outdir.mkdir(parents=True, exist_ok=True)
    columns = sorted({key for row in rows for key in row})
    with (args.outdir / "parameter_scan_results.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)

    radial_columns = (
        "case_index",
        "case_tag",
        "represented_halo_mass_msun",
        "black_hole_mass_msun",
        "sigma0_over_m_cm2_g",
        "prep_status",
        "r_in_pc",
        "r_h_pc",
        "radial_source",
        "radius_pc",
        "rho_msun_pc3",
        "sigma_1d_kms",
        "enclosed_mass_msun",
        "gravothermal_luminosity_msun_kms2_per_myr",
        "fp_direct_mass_current_msun_per_myr",
        "fp_direct_binding_current_msun_kms2_per_myr",
        "fp_direct_thermal_sink_msun_kms2_per_myr",
    )
    with (args.outdir / "parameter_scan_radial_flux_profiles.csv").open(
        "w", newline=""
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=radial_columns)
        writer.writeheader()
        writer.writerows(radial_rows)

    fit_fields = {
        "r_in_pc": "interface",
        "direct_steady_mdot_msun_per_myr": "fp",
        "direct_binding_current_msun_kms2_per_myr": "fp",
        "direct_thermal_sink_msun_kms2_per_myr": "fp",
        "direct_C_M": "fp",
        "bondi_over_direct_fp": "fp",
        "first_shell_gravothermal_luminosity_msun_kms2_per_myr": "fluid",
        "maximum_absolute_gravothermal_luminosity_msun_kms2_per_myr": "fluid",
        "absolute_thermal_sink_over_first_shell_gravothermal_luminosity": "fp",
    }
    fits = {
        field: fit
        for field, eligibility in fit_fields.items()
        if (fit := fit_power_law(rows, field, eligibility)) is not None
    }
    prep_complete = all(
        row["prep_status"] in {"READY_FOR_FP", "FP_DOMAIN_NOT_REPRESENTABLE"}
        for row in rows
    )
    all_halo_complete = all(
        row["all_halo_flux_status"] == "ALL_HALO_FLUX_COMPLETE"
        for row in rows
    )
    fp_complete = all(
        (
            row["flux_status"] is not None
            if row["fp_admissibility_status"] == "FP_ADMISSIBLE"
            else row["steady_driver_status"] == "SCAN_CASE_SKIPPED"
        )
        for row in rows
    )
    complete = prep_complete and all_halo_complete and fp_complete
    output = {
        "schema": "sidm-smbh-parameter-scan-aggregate-v1",
        "status": "SCAN_AGGREGATE_COMPLETE" if complete else "SCAN_AGGREGATE_PARTIAL",
        "updated_at": datetime.now(timezone.utc).isoformat(),
        "grid_schema": grid["schema"],
        "case_count": grid["case_count"],
        "preparation_status_counts": preparation_counts,
        "fp_flux_status_counts": flux_counts,
        "fp_admissibility_status_counts": admissibility_counts,
        "all_halo_flux_status_counts": all_halo_flux_counts,
        "completion_gates": {
            "all_preparations_classified": prep_complete,
            "all_halo_fluid_fluxes_measured": all_halo_complete,
            "all_fp_admissible_cases_have_flux_budgets": fp_complete,
        },
        "power_law_fits": fits,
        "rows": rows,
        "csv": "parameter_scan_results.csv",
        "radial_flux_csv": "parameter_scan_radial_flux_profiles.csv",
    }
    (args.outdir / "parameter_scan_results.json").write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0 if complete else 2


if __name__ == "__main__":
    raise SystemExit(main())
