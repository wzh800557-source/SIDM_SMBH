#!/usr/bin/env python3
from __future__ import annotations

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


def model(current: float, binding: float) -> dict:
    return {
        "mass": {
            "combined_value": current,
            # Deliberately put the combined estimator just outside the seed
            # interval for one-sided-error-bar coverage.
            "seed_values": [current * 1.002, current * 1.018],
        },
        "capture_binding_energy": {
            "combined_value": binding,
            "seed_values": [binding * 0.985, binding * 1.015],
        },
    }


def group(energy: int, angular: int, radius: float, current: float) -> dict:
    return {
        "energy_bins": energy,
        "angular_bins": angular,
        "reservoir_radius_pc": radius,
        "models": {
            "direct": model(current, 100.0 * current),
            "immediate": model(1.04 * current, 102.0 * current),
        },
    }


def selected_records(dimension: str) -> list[dict]:
    records = []
    for model_name in ("direct", "immediate"):
        for metric_name in ("mass", "capture_binding_energy"):
            item = {
                "model": model_name,
                "metric": metric_name,
                "fractional_difference": 0.03,
                "gate": True,
            }
            if dimension == "boundary":
                item["full_range_fraction_of_mean"] = item.pop("fractional_difference")
                scale = 1.0 if metric_name == "mass" else 100.0
                model_scale = 1.0 if model_name == "direct" else 1.04
                item.update({
                    "radii_pc": [0.9, 1.0, 1.1],
                    "values": [
                        10.4 * scale * model_scale,
                        10.5 * scale * model_scale,
                        10.6 * scale * model_scale,
                    ],
                })
            records.append(item)
    return records


def main() -> int:
    package = Path(__file__).resolve().parents[1]
    plotter = package / "plot_production_results.py"
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        convergence_path = root / "convergence.json"
        convergence = {
            "schema": "finite-angle-ej-convergence-v3",
            "status": "PRODUCTION_CONVERGENCE_PASS",
            "fractional_tolerance": 0.1,
            "required_dimensions": ["energy", "angular", "boundary"],
            "required_metrics": ["mass", "capture_binding_energy"],
            "dimension_gates": {
                name: {"pass": True}
                for name in ("energy", "angular", "boundary")
            },
            "production_resolution_path": {
                "primary_boundary_radius_pc": 1.0,
                "previous_energy_bins": 48,
                "energy_ladder_angular_bins": 129,
                "final_energy_bins": 64,
                "previous_angular_bins": 129,
                "final_angular_bins": 257,
            },
            "groups": [
                # These lower-resolution records must remain in the ledger but
                # must not crowd the production figure.
                group(32, 65, 1.0, 8.0),
                group(48, 65, 1.0, 9.0),
                group(48, 129, 1.0, 10.0),
                group(64, 129, 1.0, 10.3),
                group(64, 257, 1.0, 10.5),
                group(64, 257, 0.9, 10.4),
                group(64, 257, 1.1, 10.6),
            ],
            "selected_comparisons": {
                dimension: selected_records(dimension)
                for dimension in ("energy", "angular", "boundary")
            },
        }
        write_json(convergence_path, convergence)

        closure_path = root / "closure.json"
        closure = {
            "schema": "gnc-fluid-absolute-closure-v1",
            "status": "ABSOLUTE_CLOSURE_MEASURED",
            "gates": {
                "mass_normalization": True,
                "energy_convergence": True,
                "angular_convergence": True,
                "boundary_convergence": True,
            },
            "identity": {
                "convergence_json_sha256": sha(convergence_path),
                "energy_bins": 64,
                "angular_bins": 257,
            },
            "validated_sensitivity_envelope": {
                "mass_current": {"maximum_fractional_sensitivity": 0.08},
                "capture_binding_energy_current": {
                    "maximum_fractional_sensitivity": 0.09
                },
            },
            "measured_mdot_msun_per_myr": 10.5,
            "C_M_measured": 1.1e-4,
            "C_E_thermal_sink": -1.65e-4,
            "thermal_specific_energy_factor": 1.5,
        }
        write_json(closure_path, closure)

        comparison_path = root / "fluid_comparison.csv"
        comparison_path.write_text(
            "tau_relax,control_rho_inner_mean_msun_pc3,"
            "sink_rho_inner_mean_msun_pc3,control_sigma_inner_1d_kms,"
            "sink_sigma_inner_1d_kms,control_r_inner_pc,sink_r_inner_pc,"
            "sink_over_control_rho_inner_mean_msun_pc3,"
            "sink_over_control_sigma_inner_1d_kms,"
            "sink_over_control_r_inner_pc\n"
            "0,2,2,30,30,25,25,1,1,1\n"
            "0.5,4,4.02,31,31.031,20,19.98,1.005,1.001,0.999\n"
            "0.75,5,5.04,32,32.064,18,17.964,1.008,1.002,0.998\n"
        )
        analysis_path = root / "fluid_analysis.json"
        analysis = {
            "schema": "black-hole-aware-fluid-response-v3",
            "status": "BLACK_HOLE_AWARE_RESPONSE_COMPLETE",
            "closure_scope": "absolute_two_current_closure",
            "capture_binding_energy_current_used_by_fluid": False,
            "gates": {
                "control_completed_physically": True,
                "sink_completed_physically": True,
                "captured_mass_below_limit": True,
            },
            "identity": {
                "closure_json_sha256": sha(closure_path),
                "fluid_comparison_csv_sha256": sha(comparison_path),
            },
            "common_tau_relax": 0.75,
            "inner_mean_density_sink_over_control_at_common_end": 1.008,
            "inner_dispersion_sink_over_control_at_common_end": 1.002,
            "inner_radius_sink_over_control_at_common_end": 0.998,
            "inner_mean_density_max_abs_fractional_difference": 0.008,
            "resolved_density_diagnostic": {
                "density_field": "rho_inner_mean_msun_pc3",
                "density_is_extrapolated_central_value": False,
                "dispersion_field": "sigma_inner_1d_kms",
                "radius_field": "r_inner_pc",
            },
        }
        write_json(analysis_path, analysis)

        outdir = root / "figures"
        command = [
            sys.executable, str(plotter),
            "--convergence-json", str(convergence_path),
            "--closure-json", str(closure_path),
            "--fluid-analysis-json", str(analysis_path),
            "--fluid-comparison-csv", str(comparison_path),
            "--outdir", str(outdir),
        ]
        completed = subprocess.run(
            command, check=True, capture_output=True, text=True,
        )
        assert "PRODUCTION_FIGURES_COMPLETE" in completed.stdout
        ledger_path = outdir / "production_figure_ledger.json"
        ledger = json.loads(ledger_path.read_text())
        assert ledger["schema"] == "gnc-production-figure-ledger-v3"
        assert ledger["status"] == "PRODUCTION_FIGURES_COMPLETE"
        assert ledger["plotted_resolution_path"] == [
            "48/129", "64/129", "64/257"
        ]
        assert len(ledger["outputs"]) == 4
        for name in ledger["outputs"]:
            path = outdir / name
            assert path.is_file() and path.stat().st_size > 1000
            assert ledger["output_sha256"][name] == sha(path)
            if path.suffix == ".png":
                assert path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
            else:
                assert path.read_bytes().startswith(b"%PDF-")
        assert ledger["sources"] == {
            "convergence_sha256": sha(convergence_path),
            "closure_sha256": sha(closure_path),
            "fluid_analysis_sha256": sha(analysis_path),
            "fluid_comparison_sha256": sha(comparison_path),
        }

        # A physical mass-only closure retains the failed binding-energy
        # comparison in the figure while plotting the accepted mass current
        # and fluid response.
        mass_convergence = json.loads(convergence_path.read_text())
        mass_convergence["status"] = "INCOMPLETE_OR_FAILED"
        mass_convergence["dimension_gates"]["energy"]["pass"] = False
        for item in mass_convergence["selected_comparisons"]["energy"]:
            if item["metric"] == "capture_binding_energy":
                item["fractional_difference"] = 0.141
                item["gate"] = False
        mass_convergence["observable_gates"] = {
            "mass": {
                "status": "PRODUCTION_OBSERVABLE_CONVERGENCE_PASS",
                "pass": True,
            },
            "capture_binding_energy": {
                "status": "INCOMPLETE_OR_FAILED",
                "pass": False,
            },
        }
        write_json(convergence_path, mass_convergence)
        mass_closure = {
            "schema": "gnc-fluid-mass-closure-v1",
            "status": "FLUID_MASS_CLOSURE_MEASURED",
            "capture_binding_energy_current_used_by_fluid": False,
            "gates": {"all_mass_gates": True},
            "identity": {
                "convergence_json_sha256": sha(convergence_path),
                "energy_bins": 64,
                "angular_bins": 257,
            },
            "validated_mass_sensitivity_envelope": {
                "maximum_fractional_sensitivity": 0.08,
            },
            "measured_mdot_msun_per_myr": 10.5,
            "C_M_measured": 1.1e-4,
            "C_E_thermal_sink": -1.65e-4,
            "thermal_specific_energy_factor": 1.5,
        }
        write_json(closure_path, mass_closure)
        analysis["closure_scope"] = "mass_current_closure_only"
        analysis["identity"]["closure_json_sha256"] = sha(closure_path)
        write_json(analysis_path, analysis)
        mass_outdir = root / "mass_figures"
        mass_command = command[:-1] + [str(mass_outdir)]
        mass_completed = subprocess.run(
            mass_command, check=True, capture_output=True, text=True,
        )
        assert "PRODUCTION_FIGURES_COMPLETE" in mass_completed.stdout
        mass_ledger = json.loads(
            (mass_outdir / "production_figure_ledger.json").read_text()
        )
        assert mass_ledger["closure_scope"] == "mass_current_closure_only"
        assert mass_ledger["capture_binding_energy_current_accepted"] is False

        analysis_uses_binding = json.loads(analysis_path.read_text())
        analysis_uses_binding["capture_binding_energy_current_used_by_fluid"] = True
        write_json(analysis_path, analysis_uses_binding)
        rejected_binding = subprocess.run(
            mass_command[:-1] + [str(root / "rejected_binding_figures")],
            check=False, capture_output=True, text=True,
        )
        assert rejected_binding.returncode != 0
        assert "uses the unresolved binding-energy current" in rejected_binding.stderr
        analysis_uses_binding["capture_binding_energy_current_used_by_fluid"] = False
        write_json(analysis_path, analysis_uses_binding)

        # Provenance is fail-closed: a table changed after analysis must not be
        # rendered into an apparently accepted production figure.
        comparison_path.write_text(comparison_path.read_text() + "\n")
        rejected = subprocess.run(
            command[:-1] + [str(root / "rejected_figures")],
            check=False, capture_output=True, text=True,
        )
        assert rejected.returncode != 0
        assert "does not match the comparison table" in rejected.stderr

    print("PASS: production plots are sparse, validated, and provenance-locked")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
