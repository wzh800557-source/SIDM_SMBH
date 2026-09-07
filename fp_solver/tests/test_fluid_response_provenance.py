#!/usr/bin/env python3
from __future__ import annotations

import csv
import hashlib
import json
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
            "final_rho_c_msun_pc3": 4.0,
            "captured_mass_msun": 0.0,
            "n_conduction": 1,
        }
        sink_branch = dict(branch)
        sink_branch.update({
            "captured_mass_fraction_of_first_fluid_shell": 0.01,
            "captured_mass_fraction_of_GNC_owned_domain": 0.10,
            "Mdot_msun_per_myr": 10.0,
            "L_inner_msun_kms2_per_myr": -13500.0,
            "final_rho_c_msun_pc3": 4.04,
            "captured_mass_msun": 20.0,
        })
        summary = {
            "schema": "gnc-fluid-feedback-v4",
            "status": "MEASURED_FLUID_FEEDBACK_COMPLETE",
            "input_closure_accepted": True,
            "closure_update_mode": "fixed_at_matched_snapshot",
            "M_bh_msun": 4.0e6,
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
            "branch", "tau_relax", "rho_c_msun_pc3", "sigma_c_kms",
            "r0_pc", "M_bh_msun", "lmfp_scaleheight_factor_inner",
            "lmfp_scaleheight_factor_min", "time_myr",
            "captured_mass_msun", "L_inner_msun_kms2_per_myr",
            "E_boundary_msun_kms2", "E_total_code", "n_conduction",
        ]
        rows = [
            ["control", 0.0, 2.0, 30.0, 1.0, 4.0e6, 0.2, 0.1,
             0.0, 0.0, 0.0, 0.0, -10.0, 0],
            ["control", 1.0, 4.0, 31.0, 0.8, 4.0e6, 0.18, 0.09,
             2.0, 0.0, 0.0, 0.0, -9.0, 1],
            ["sink", 0.0, 2.0, 30.0, 1.0, 4.0e6, 0.2, 0.1,
             0.0, 0.0, -13500.0, 0.0, -10.0, 0],
            ["sink", 1.0, 4.04, 31.1, 0.79, 4.0e6, 0.18, 0.09,
             2.0, 20.0, -13500.0, -27000.0, -9.1, 1],
        ]
        with (root / "trajectories.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(fields)
            writer.writerows(rows)

        passed = run_analyzer(package, root, "pass")
        assert passed.returncode == 0, passed.stdout + passed.stderr
        accepted = json.loads((root / "analysis_pass.json").read_text())
        assert accepted["schema"] == "black-hole-aware-fluid-response-v2"
        assert accepted["status"] == "BLACK_HOLE_AWARE_RESPONSE_COMPLETE"
        assert all(accepted["gates"].values())
        assert accepted["identity"]["closure_json_sha256"] == sha(closure)
        assert accepted["identity"]["fluid_comparison_csv_sha256"] == sha(
            root / "comparison_pass.csv"
        )

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
            ["control", 0.0, 2.0, 30.0, 1.0, 4.0e6, 0.2, 0.1,
             0.0, 0.0, 0.0, 0.0, -10.0, 0],
            ["control", 0.75, 7.0, 31.0, 0.82, 4.0e6, 0.18, 0.09,
             1.5, 0.0, 0.0, 0.0, -8.5, 1],
            ["sink", 0.0, 2.0, 30.0, 1.0, 4.0e6, 0.2, 0.1,
             0.0, 0.0, -13500.0, 0.0, -10.0, 0],
            ["sink", 0.5, 4.0, 30.5, 0.90, 4.0e6, 0.19, 0.095,
             1.0, 10.0, -13500.0, -13500.0, -9.5, 1],
            ["sink", 1.0, 8.0, 31.5, 0.80, 4.0e6, 0.17, 0.085,
             2.0, 20.0, -13500.0, -27000.0, -8.8, 2],
        ]
        summary["branches"]["control"].update({
            "final_tau_relax": 0.75,
            "final_time_myr": 1.5,
            "final_rho_c_msun_pc3": 7.0,
            "captured_mass_msun": 0.0,
            "n_conduction": 1,
        })
        summary["branches"]["sink"].update({
            "final_tau_relax": 1.0,
            "final_time_myr": 2.0,
            "final_rho_c_msun_pc3": 8.0,
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
                "sink_over_control_rho_c_msun_pc3"
            ])
            - staggered_analysis["density_sink_over_control_at_common_end"]
        ) < 1.0e-12

        summary["branches"]["control"].update({
            "final_tau_relax": 1.0,
            "final_time_myr": 2.0,
            "final_rho_c_msun_pc3": 4.0,
            "captured_mass_msun": 0.0,
            "n_conduction": 1,
        })
        summary["branches"]["sink"].update({
            "final_tau_relax": 1.0,
            "final_time_myr": 2.0,
            "final_rho_c_msun_pc3": 4.04,
            "captured_mass_msun": 20.0,
            "n_conduction": 1,
        })
        write_json(root / "summary.json", summary)

        with (root / "trajectories.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(fields)
            writer.writerows(rows)

        rows[-1][10] = -13000.0
        with (root / "trajectories.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(fields)
            writer.writerows(rows)
        varying_current = run_analyzer(package, root, "varying_current")
        assert varying_current.returncode == 2, varying_current.stderr
        varying = json.loads((root / "analysis_varying_current.json").read_text())
        assert not varying["gates"]["fixed_boundary_currents_in_trajectories"]
        rows[-1][10] = -13500.0
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
