#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def run_analyzer(package: Path, root: Path, suffix: str) -> subprocess.CompletedProcess:
    return subprocess.run([
        sys.executable, str(package / "analyze_fluid_response.py"),
        "--summary", str(root / "summary.json"),
        "--trajectories", str(root / "trajectories.csv"),
        "--closure-json", str(root / "closure.json"),
        "--profile", str(root / "profile.txt"),
        "--bridge-json", str(root / "bridge.json"),
        "--out-json", str(root / f"analysis_{suffix}.json"),
        "--out-csv", str(root / f"comparison_{suffix}.csv"),
    ], capture_output=True, text=True)


def main() -> int:
    package = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        profile = root / "profile.txt"
        profile.write_text("1 2 3\n")
        bridged_sha = hashlib.sha256(b"bridged profile").hexdigest()
        bridge = root / "bridge.json"
        write_json(bridge, {
            "schema": "absolute-closure-hydrostatic-bridge-v3",
            "status": "OK",
            "output_profile_sha256": bridged_sha,
        })
        closure = root / "closure.json"
        write_json(closure, {
            "schema": "gnc-fluid-absolute-closure-v1",
            "status": "ABSOLUTE_CLOSURE_MEASURED",
            "gates": {"all_scientific_gates": True},
            "measured_mdot_msun_per_myr": 10.0,
            "thermal_sink_current_msun_kms2_per_myr": -13500.0,
            "identity": {
                "fluid_profile_sha256": sha(profile),
                "bridge_json_sha256": sha(bridge),
                "bridged_profile_sha256": bridged_sha,
            },
        })

        inner_shell_mass = 8.0 * math.pi / 3.0

        def mean_density(radius: float) -> float:
            return 3.0 * inner_shell_mass / (4.0 * math.pi * radius**3)

        branch = {
            "stop_reason": "tau_end",
            "conduction_energy_budget_relative_error": 1.0e-14,
            "black_hole_remap": {"M_bh_final_msun": 4.0e6},
            "captured_mass_fraction_of_first_fluid_shell": 0.0,
            "captured_mass_fraction_of_GNC_owned_domain": 0.0,
            "maximum_unmodeled_fluid_shell_mass_fraction": 0.02,
            "maximum_unmodeled_GNC_mass_fraction": 0.20,
            "Mdot_msun_per_myr": 0.0,
            "L_inner_msun_kms2_per_myr": 0.0,
            "final_tau_relax": 1.0,
            "final_time_myr": 2.0,
            "final_rho_inner_mean_msun_pc3": mean_density(0.8),
            "final_rho_c_msun_pc3": mean_density(0.8),
            "inner_shell_mass_msun": inner_shell_mass,
            "captured_mass_msun": 0.0,
            "n_conduction": 1,
        }
        sink_branch = dict(branch)
        sink_branch.update({
            "captured_mass_fraction_of_first_fluid_shell": 0.01,
            "captured_mass_fraction_of_GNC_owned_domain": 0.10,
            "Mdot_msun_per_myr": 10.0,
            "L_inner_msun_kms2_per_myr": -13500.0,
            "final_rho_inner_mean_msun_pc3": mean_density(0.79),
            "final_rho_c_msun_pc3": mean_density(0.79),
            "captured_mass_msun": 20.0,
        })
        summary = {
            "schema": "gnc-fluid-feedback-v5",
            "status": "MEASURED_FLUID_FEEDBACK_COMPLETE",
            "input_closure_accepted": True,
            "closure_update_mode": "fixed_at_matched_snapshot",
            "M_bh_msun": 4.0e6,
            "resolved_diagnostics": {
                "density_field": "rho_inner_mean_msun_pc3",
                "density_definition": (
                    "mean density inside the innermost resolved Lagrangian shell, "
                    "3*M(<r_inner)/(4*pi*r_inner^3)"
                ),
                "density_is_extrapolated_central_value": False,
                "dispersion_field": "sigma_inner_1d_kms",
                "dispersion_definition": (
                    "one-dimensional velocity dispersion in the innermost "
                    "resolved Lagrangian shell"
                ),
                "radius_field": "r_inner_pc",
            },
            "input_identity": {
                "closure_json_sha256": sha(closure),
                "fluid_profile_sha256": sha(profile),
                "bridge_json_sha256": sha(bridge),
                "bridged_profile_sha256": bridged_sha,
            },
            "applied_currents": {
                "Mdot_msun_per_myr": 10.0,
                "control_luminosity_msun_kms2_per_myr": 0.0,
                "sink_luminosity_msun_kms2_per_myr": -13500.0,
                "source_luminosity_msun_kms2_per_myr": 0.0,
            },
            "branches": {"control": branch, "sink": sink_branch},
        }
        write_json(root / "summary.json", summary)

        fields = [
            "branch", "tau_relax", "rho_inner_mean_msun_pc3",
            "sigma_inner_1d_kms", "r_inner_pc", "inner_shell_mass_msun",
            "rho_c_msun_pc3", "sigma_c_kms", "r0_pc",
            "M_bh_msun", "lmfp_scaleheight_factor_inner",
            "lmfp_scaleheight_factor_min", "time_myr",
            "captured_mass_msun", "L_inner_msun_kms2_per_myr",
            "E_boundary_msun_kms2", "E_total_code", "n_conduction",
        ]

        def trajectory_row(
            branch_name: str, tau: float, radius: float, sigma: float,
            scaleheight_inner: float, scaleheight_min: float, time_myr: float,
            captured_mass: float, luminosity: float, boundary_energy: float,
            total_energy: float, n_conduction: int,
        ) -> list:
            density = mean_density(radius)
            return [
                branch_name, tau, density, sigma, radius, inner_shell_mass,
                density, sigma, radius, 4.0e6, scaleheight_inner,
                scaleheight_min, time_myr, captured_mass, luminosity,
                boundary_energy, total_energy, n_conduction,
            ]

        rows = [
            trajectory_row(
                "control", 0.0, 1.0, 30.0, 0.2, 0.1,
                0.0, 0.0, 0.0, 0.0, -10.0, 0,
            ),
            trajectory_row(
                "control", 1.0, 0.8, 31.0, 0.18, 0.09,
                2.0, 0.0, 0.0, 0.0, -9.0, 1,
            ),
            trajectory_row(
                "sink", 0.0, 1.0, 30.0, 0.2, 0.1,
                0.0, 0.0, -13500.0, 0.0, -10.0, 0,
            ),
            trajectory_row(
                "sink", 1.0, 0.79, 31.1, 0.18, 0.09,
                2.0, 20.0, -13500.0, -27000.0, -9.1, 1,
            ),
        ]
        with (root / "trajectories.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(fields)
            writer.writerows(rows)

        passed = run_analyzer(package, root, "pass")
        assert passed.returncode == 0, passed.stdout + passed.stderr
        accepted = json.loads((root / "analysis_pass.json").read_text())
        assert accepted["schema"] == "black-hole-aware-fluid-response-v3"
        assert accepted["status"] == "BLACK_HOLE_AWARE_RESPONSE_COMPLETE"
        assert all(accepted["gates"].values())
        assert accepted["identity"]["closure_json_sha256"] == sha(closure)
        assert accepted["identity"]["fluid_comparison_csv_sha256"] == sha(
            root / "comparison_pass.csv"
        )

        rows[-1][6] *= 1.01
        with (root / "trajectories.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(fields)
            writer.writerows(rows)
        ambiguous_density = run_analyzer(package, root, "ambiguous_density")
        assert ambiguous_density.returncode != 0
        assert "inconsistent rho_inner_mean_msun_pc3 aliases" in (
            ambiguous_density.stderr
        )
        rows[-1][6] /= 1.01
        with (root / "trajectories.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(fields)
            writer.writerows(rows)

        rows[-1][2] *= 1.01
        rows[-1][6] *= 1.01
        with (root / "trajectories.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(fields)
            writer.writerows(rows)
        wrong_mean_density = run_analyzer(package, root, "wrong_mean_density")
        assert wrong_mean_density.returncode != 0
        assert "does not satisfy the innermost-shell mean-density definition" in (
            wrong_mean_density.stderr
        )
        rows[-1][2] /= 1.01
        rows[-1][6] /= 1.01
        with (root / "trajectories.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(fields)
            writer.writerows(rows)

        # A converged mass-current closure may drive the physical thermal sink
        # even when the capture-binding-energy moment remains unresolved. The
        # latter must be absent from the applied fluid currents.
        absolute_closure_value = json.loads(closure.read_text())
        mass_closure_value = dict(absolute_closure_value)
        mass_closure_value.update({
            "schema": "gnc-fluid-mass-closure-v1",
            "status": "FLUID_MASS_CLOSURE_MEASURED",
            "capture_binding_energy_current_used_by_fluid": False,
        })
        write_json(closure, mass_closure_value)
        summary.update({
            "status": "MASS_CLOSURE_FLUID_RESPONSE_COMPLETE",
            "absolute_two_current_closure": False,
            "mass_current_closure": True,
        })
        summary["applied_currents"][
            "source_luminosity_msun_kms2_per_myr"
        ] = None
        summary["input_identity"]["closure_json_sha256"] = sha(closure)
        write_json(root / "summary.json", summary)
        mass_passed = run_analyzer(package, root, "mass_pass")
        assert mass_passed.returncode == 0, mass_passed.stdout + mass_passed.stderr
        mass_analysis = json.loads((root / "analysis_mass_pass.json").read_text())
        assert mass_analysis["status"] == "BLACK_HOLE_AWARE_RESPONSE_COMPLETE"
        assert mass_analysis["closure_scope"] == "mass_current_closure_only"
        assert mass_analysis["capture_binding_energy_current_used_by_fluid"] is False
        assert all(mass_analysis["gates"].values())

        write_json(closure, absolute_closure_value)
        summary.update({
            "status": "MEASURED_FLUID_FEEDBACK_COMPLETE",
            "absolute_two_current_closure": True,
            "mass_current_closure": False,
        })
        summary["applied_currents"][
            "source_luminosity_msun_kms2_per_myr"
        ] = 0.0
        summary["input_identity"]["closure_json_sha256"] = sha(closure)
        write_json(root / "summary.json", summary)

        closure_value = json.loads(closure.read_text())
        closure_value["gates"]["all_scientific_gates"] = False
        write_json(closure, closure_value)
        failed_gate = run_analyzer(package, root, "failed_gate")
        assert failed_gate.returncode != 0
        assert "missing or failed gate" in failed_gate.stderr
        closure_value["gates"]["all_scientific_gates"] = True
        write_json(closure, closure_value)

        # If the control branch ends between two sink outputs, the comparison
        # must interpolate both branches at the actual common endpoint rather
        # than report the preceding sink output as that endpoint.
        staggered_rows = [
            trajectory_row(
                "control", 0.0, 1.0, 30.0, 0.2, 0.1,
                0.0, 0.0, 0.0, 0.0, -10.0, 0,
            ),
            trajectory_row(
                "control", 0.75, 0.82, 31.0, 0.18, 0.09,
                1.5, 0.0, 0.0, 0.0, -8.5, 1,
            ),
            trajectory_row(
                "sink", 0.0, 1.0, 30.0, 0.2, 0.1,
                0.0, 0.0, -13500.0, 0.0, -10.0, 0,
            ),
            trajectory_row(
                "sink", 0.5, 0.90, 30.5, 0.19, 0.095,
                1.0, 10.0, -13500.0, -13500.0, -9.5, 1,
            ),
            trajectory_row(
                "sink", 1.0, 0.80, 31.5, 0.17, 0.085,
                2.0, 20.0, -13500.0, -27000.0, -8.8, 2,
            ),
        ]
        summary["branches"]["control"].update({
            "final_tau_relax": 0.75,
            "final_time_myr": 1.5,
            "final_rho_inner_mean_msun_pc3": mean_density(0.82),
            "final_rho_c_msun_pc3": mean_density(0.82),
            "captured_mass_msun": 0.0,
            "n_conduction": 1,
        })
        summary["branches"]["sink"].update({
            "final_tau_relax": 1.0,
            "final_time_myr": 2.0,
            "final_rho_inner_mean_msun_pc3": mean_density(0.80),
            "final_rho_c_msun_pc3": mean_density(0.80),
            "captured_mass_msun": 20.0,
            "n_conduction": 2,
        })
        write_json(root / "summary.json", summary)
        with (root / "trajectories.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(fields)
            writer.writerows(staggered_rows)
        staggered = run_analyzer(package, root, "staggered")
        assert staggered.returncode == 0, staggered.stdout + staggered.stderr
        staggered_analysis = json.loads(
            (root / "analysis_staggered.json").read_text()
        )
        with (root / "comparison_staggered.csv").open(newline="") as stream:
            staggered_comparison = list(csv.DictReader(stream))
        assert float(staggered_comparison[-1]["tau_relax"]) == 0.75
        assert abs(
            float(staggered_comparison[-1][
                "sink_over_control_rho_inner_mean_msun_pc3"
            ])
            - staggered_analysis[
                "inner_mean_density_sink_over_control_at_common_end"
            ]
        ) < 1.0e-12

        summary["branches"]["control"].update({
            "final_tau_relax": 1.0,
            "final_time_myr": 2.0,
            "final_rho_inner_mean_msun_pc3": mean_density(0.8),
            "final_rho_c_msun_pc3": mean_density(0.8),
            "captured_mass_msun": 0.0,
            "n_conduction": 1,
        })
        summary["branches"]["sink"].update({
            "final_tau_relax": 1.0,
            "final_time_myr": 2.0,
            "final_rho_inner_mean_msun_pc3": mean_density(0.79),
            "final_rho_c_msun_pc3": mean_density(0.79),
            "captured_mass_msun": 20.0,
            "n_conduction": 1,
        })
        write_json(root / "summary.json", summary)

        with (root / "trajectories.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(fields)
            writer.writerows(rows)

        rows[-1][14] = -13000.0
        with (root / "trajectories.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(fields)
            writer.writerows(rows)
        varying_current = run_analyzer(package, root, "varying_current")
        assert varying_current.returncode == 2, varying_current.stderr
        varying = json.loads((root / "analysis_varying_current.json").read_text())
        assert not varying["gates"]["fixed_boundary_currents_in_trajectories"]
        rows[-1][14] = -13500.0
        with (root / "trajectories.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(fields)
            writer.writerows(rows)

        summary["branches"]["control"]["stop_reason"] = "wall_limit"
        write_json(root / "summary.json", summary)
        wall_limited = run_analyzer(package, root, "wall")
        assert wall_limited.returncode == 2, wall_limited.stderr
        wall = json.loads((root / "analysis_wall.json").read_text())
        assert not wall["gates"]["control_completed_physically"]

        summary["branches"]["control"]["stop_reason"] = "tau_end"
        write_json(root / "summary.json", summary)
        closure_value = json.loads(closure.read_text())
        closure_value["measured_mdot_msun_per_myr"] = 11.0
        write_json(closure, closure_value)
        wrong_closure = run_analyzer(package, root, "wrong_closure")
        assert wrong_closure.returncode == 2, wrong_closure.stderr
        rejected = json.loads((root / "analysis_wrong_closure.json").read_text())
        assert not rejected["gates"]["closure_file_identity"]
        assert not rejected["gates"]["measured_mass_current_applied"]

    print("PASS: fluid-response provenance and completion gates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
