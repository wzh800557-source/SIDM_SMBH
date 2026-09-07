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


def write_steady(
    path: Path, profile_sha: str, df_sha: str, energy: int, angular: int,
    radius: float, current: float, seed: int | None,
) -> None:
    metadata = {
        "schema": "finite-angle-ej-jump-operator-v1",
        "profile_sha256": profile_sha,
        "df_table_sha256": df_sha,
        "mbh_msun": 4.0e6,
        "sigma0_over_m_cm2_g": 100.0,
        "w_kms": 80.0,
        "energy_bins": energy,
        "angular_bins": angular,
        "angular_grid": {
            17: "default", 33: "fine", 65: "ultra", 129: "hyper",
            257: "j257",
        }[angular],
        "reservoir_radius_pc": radius,
        "state_count": 2,
    }
    if seed is None:
        metadata["seeds"] = [314159, 271828]
    else:
        metadata["seed"] = seed
    models = {}
    for model, factor in (("direct", 1.0), ("immediate", 1.05)):
        value = current * factor
        models[model] = {
            "status": "EJ_STEADY_SMOKE_PASS",
            "steady_capture_mdot_msun_per_myr": value,
            "steady_capture_binding_energy_current_msun_kms2_per_myr": 100.0 * value,
        }
    path.write_text(json.dumps({
        "schema": "finite-angle-ej-steady-state-v2",
        "status": "EJ_STEADY_SMOKE_PASS",
        "operator_metadata": metadata,
        "models": models,
    }, indent=2) + "\n")


def main() -> int:
    package = Path(__file__).resolve().parents[1]
    slurm = package.parent / "slurm" / "fp_solver"
    aggregate_script = slurm / "run_production_ej_aggregate.sbatch"
    aggregate_text = aggregate_script.read_text()
    assert "\nGROUPS=(" not in aggregate_text
    assert "\nRUN_GROUPS=(" in aggregate_text
    assert '${RUN_GROUPS[@]}' in aggregate_text
    subprocess.run(
        ["bash", "-n", str(aggregate_script)],
        check=True,
        capture_output=True,
        text=True,
    )
    array_script = slurm / "run_production_ej_array.sbatch"
    array_text = array_script.read_text()
    assert "PAIRS_PER_BIN=${PAIRS_PER_BIN:-40000}" in array_text
    assert '--pairs-per-bin "$PAIRS_PER_BIN"' in array_text
    assert "--pairs-per-bin 40000" not in array_text
    subprocess.run(
        ["bash", "-n", str(array_script)],
        check=True,
        capture_output=True,
        text=True,
    )
    resolution_array = slurm / "run_production_ej_resolution_array.sbatch"
    resolution_array_text = resolution_array.read_text()
    assert "#SBATCH --array=0-9" in resolution_array_text
    assert "E=24; GRID=fine" in resolution_array_text
    assert "E=32; GRID=fine" in resolution_array_text
    assert "E=32; GRID=ultra" in resolution_array_text
    assert "PAIRS_PER_BIN=${PAIRS_PER_BIN:-1600000}" in resolution_array_text
    subprocess.run(
        ["bash", "-n", str(resolution_array)],
        check=True,
        capture_output=True,
        text=True,
    )
    resolution_aggregate = (
        slurm / "run_production_ej_resolution_aggregate.sbatch"
    )
    resolution_aggregate_text = resolution_aggregate.read_text()
    assert "local_circular_e32_ultra" in resolution_aggregate_text
    assert '--primary-boundary-radius "$PRIMARY"' in resolution_aggregate_text
    assert "$ROOT/aggregate_v3" in resolution_aggregate_text
    subprocess.run(
        ["bash", "-n", str(resolution_aggregate)],
        check=True,
        capture_output=True,
        text=True,
    )
    resolution2_array = slurm / "run_production_ej_resolution2_array.sbatch"
    resolution2_array_text = resolution2_array.read_text()
    assert "#SBATCH --array=0-7" in resolution2_array_text
    assert "E=48; GRID=ultra" in resolution2_array_text
    assert "E=48; GRID=hyper" in resolution2_array_text
    assert "PAIRS_PER_BIN=${PAIRS_PER_BIN:-1600000}" in resolution2_array_text
    subprocess.run(
        ["bash", "-n", str(resolution2_array)],
        check=True, capture_output=True, text=True,
    )
    resolution2_aggregate = (
        slurm / "run_production_ej_resolution2_aggregate.sbatch"
    )
    resolution2_aggregate_text = resolution2_aggregate.read_text()
    assert "local_circular_e48_hyper" in resolution2_aggregate_text
    assert "$ROOT/aggregate_v4" in resolution2_aggregate_text
    assert "$ROOT/AGGREGATE_V4_STATUS.txt" in resolution2_aggregate_text
    assert "STEADY_ROOT=${STEADY_ROOT:-$ROOT/steady}" in resolution2_aggregate_text
    assert 'tee "$AGGREGATE_STATUS"' in resolution2_aggregate_text
    subprocess.run(
        ["bash", "-n", str(resolution2_aggregate)],
        check=True, capture_output=True, text=True,
    )
    resolution3_array = slurm / "run_production_ej_resolution3_array.sbatch"
    resolution3_array_text = resolution3_array.read_text()
    assert "#SBATCH --array=0-7" in resolution3_array_text
    assert "E=64; GRID=hyper" in resolution3_array_text
    assert "E=64; GRID=j257" in resolution3_array_text
    assert "PAIRS_PER_BIN=${PAIRS_PER_BIN:-1600000}" in resolution3_array_text
    subprocess.run(
        ["bash", "-n", str(resolution3_array)],
        check=True, capture_output=True, text=True,
    )
    resolution3_aggregate = (
        slurm / "run_production_ej_resolution3_aggregate.sbatch"
    )
    resolution3_aggregate_text = resolution3_aggregate.read_text()
    assert "local_circular_e64_hyper" in resolution3_aggregate_text
    assert "local_circular_e64_j257" in resolution3_aggregate_text
    assert "$ROOT/steady_v5" in resolution3_aggregate_text
    assert "$ROOT/aggregate_v6" in resolution3_aggregate_text
    assert "$ROOT/AGGREGATE_V6_STATUS.txt" in resolution3_aggregate_text
    assert "V5_RESOLUTION_ONLY_FAILURE_CONFIRMED" in resolution3_aggregate_text
    subprocess.run(
        ["bash", "-n", str(resolution3_aggregate)],
        check=True, capture_output=True, text=True,
    )
    bridge_text = (package / "hydrostatic_bridge.py").read_text()
    assert '"input_profile": provenance_name(args.profile)' in bridge_text
    assert '"output_profile": provenance_name(args.out_profile)' in bridge_text
    assert 'str(args.profile.resolve())' not in bridge_text
    assert 'str(args.out_profile.resolve())' not in bridge_text
    feedback_text = (package / "measured_fluid_feedback.py").read_text()
    assert "E_conduction_budget_scale_code" in feedback_text
    assert "np.dot(dm, np.abs(self.delta_uc))" in feedback_text
    assert '"schema": "gnc-fluid-feedback-v4"' in feedback_text
    assert '"closure_json_sha256": sha256_file(args.diagnostics)' in feedback_text
    response_text = (package / "analyze_fluid_response.py").read_text()
    assert "value is not None" in response_text
    assert "math.isfinite(float(value))" in response_text
    assert '"schema": "black-hole-aware-fluid-response-v2"' in response_text
    assert '"control_completed_physically"' in response_text
    assert 'set(rows) != {"control", "sink"}' in response_text
    assert '"summary_trajectory_endpoints_match"' in response_text
    measured_text = (package / "measured_fluid_feedback.py").read_text()
    assert 'closure_gates = diag.get("gates")' in measured_text
    assert 'accepted closure contains a missing or failed gate' in measured_text
    fluid_script = (slurm / "run_production_fluid_response.sbatch").read_text()
    assert '--closure-json "$AGGREGATE_DIR/absolute_closure.json"' in fluid_script
    assert '--profile "$ROOT/remap/profile_bh.txt"' in fluid_script
    assert '--bridge-json "$ROOT/bridge/bridge.json"' in fluid_script
    assert '--branches control,sink' in fluid_script
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        source = root / "source.txt"; source.write_text("source\n")
        post = root / "post.txt"; post.write_text("post\n")
        bridged = root / "bridged.txt"; bridged.write_text("bridge\n")
        df = root / "df.txt"; df.write_text("df\n")
        configurations = [
            # The deliberately non-converged lower ladder must remain in the
            # report without vetoing a converged highest-resolution path.
            (8, 17, 1.0, 6.0),
            (12, 17, 1.0, 7.0),
            (16, 17, 1.0, 8.0),
            (16, 33, 1.0, 8.0),
            (24, 33, 1.0, 10.0),
            (32, 33, 1.0, 10.2),
            (32, 65, 1.0, 10.3),
            (32, 65, 0.9, 10.2),
            (32, 65, 1.1, 10.4),
            (48, 65, 1.0, 10.35),
            (48, 129, 1.0, 10.4),
            (48, 129, 0.9, 10.3),
            (48, 129, 1.1, 10.5),
            (64, 129, 1.0, 10.45),
            (64, 257, 1.0, 10.5),
            (64, 257, 0.9, 10.4),
            (64, 257, 1.1, 10.6),
        ]
        results = []
        selected = None
        for index, (energy, angular, radius, current) in enumerate(configurations):
            for seed, perturbation in ((314159, 0.98), (271828, 1.02)):
                path = root / f"steady_{index}_{seed}.json"
                write_steady(
                    path, sha(bridged), sha(df), energy, angular, radius,
                    current * perturbation, seed,
                )
                results.append(path)
            combined = root / f"steady_{index}_combined.json"
            write_steady(
                combined, sha(bridged), sha(df), energy, angular, radius,
                current, None,
            )
            results.append(combined)
            if (energy, angular, radius) == (64, 257, 1.0):
                selected = combined
        convergence = root / "convergence.json"
        command = [
            sys.executable, str(package / "assess_ej_convergence.py"),
            *(str(path) for path in results),
            "--out-json", str(convergence),
            "--out-csv", str(root / "convergence.csv"),
            "--primary-boundary-radius", "1.0",
            "--require-capture-energy",
        ]
        subprocess.run(command, check=True, capture_output=True, text=True)
        convergence_result = json.loads(convergence.read_text())
        assert convergence_result["status"] == "PRODUCTION_CONVERGENCE_PASS"
        assert all(
            group["numerical_solver_gate"]
            for group in convergence_result["groups"]
        )
        path = convergence_result["production_resolution_path"]
        assert path["final_energy_bins"] == 64
        assert path["final_angular_bins"] == 257
        assert path["energy_ladder_angular_bins"] == 129
        assert path["previous_energy_bins"] == 48
        assert path["previous_angular_bins"] == 129
        assert any(
            not item["gate"]
            for item in convergence_result["comparisons"]["energy"]
            if item["upper_energy_bins"] < 32
        )
        assert all(
            item["gate"]
            for values in convergence_result["selected_comparisons"].values()
            for item in values
        )
        victim = results[0]
        original_victim = victim.read_text()
        failed_value = json.loads(original_victim)
        failed_value["status"] = "FAIL"
        victim.write_text(json.dumps(failed_value, indent=2) + "\n")
        failed_convergence = root / "failed_numerical.json"
        failed_numerical = subprocess.run([
            sys.executable, str(package / "assess_ej_convergence.py"),
            *(str(path) for path in results),
            "--out-json", str(failed_convergence),
            "--out-csv", str(root / "failed_numerical.csv"),
            "--primary-boundary-radius", "1.0",
            "--require-capture-energy",
        ], capture_output=True, text=True)
        victim.write_text(original_victim)
        assert failed_numerical.returncode == 2
        failed_report = json.loads(failed_convergence.read_text())
        assert failed_report["status"] == "INCOMPLETE_OR_FAILED"
        assert any(
            not group["numerical_solver_gate"]
            for group in failed_report["groups"]
        )
        duplicate = subprocess.run([
            sys.executable, str(package / "assess_ej_convergence.py"),
            str(results[0]), str(results[0]), str(results[2]),
            "--out-json", str(root / "duplicate.json"),
            "--out-csv", str(root / "duplicate.csv"),
        ], capture_output=True, text=True)
        assert duplicate.returncode != 0
        assert "repeats an independent seed" in duplicate.stderr

        remap = root / "remap.json"
        remap.write_text(json.dumps({
            "schema": "adiabatic-bh-remapped-fluid-profile-v1",
            "status": "BH_REMAP_COMPLETE",
            "source_profile_sha256": sha(source),
            "output_profile_sha256": sha(post),
            "M_bh_msun": 4.0e6,
            "sigma0_over_m_cm2_g": 100.0,
            "yukawa_w_kms": 80.0,
        }) + "\n")
        bridge = root / "bridge.json"
        bridge.write_text(json.dumps({
            "schema": "absolute-closure-hydrostatic-bridge-v3",
            "status": "OK",
            "fluid_subgrid_bridge_gate": "PASS",
            "input_profile_sha256": sha(post),
            "output_profile_sha256": sha(bridged),
            "mbh_msun": 4.0e6,
            "sigma_over_m_cm2_g": 100.0,
            "yukawa_w_kms": 80.0,
            "r_in_pc": 1.0,
            "rh_pc": 20.0,
        }) + "\n")
        norm = root / "norm.json"
        norm.write_text(json.dumps({
            "schema": "gnc-normalized-df-v6",
            "status": "PASS",
            "profile_sha256": sha(bridged),
            "df_table_sha256": sha(df),
            "mbh_msun": 4.0e6,
            "sigma0_over_m_cm2_g": 100.0,
            "w_kms": 80.0,
            "r_boundary_pc": 1.0,
            "rho_boundary_msun_pc3": 5.0,
            "sigma_boundary_kms": 30.0,
        }) + "\n")
        closure = root / "closure.json"
        subprocess.run([
            sys.executable, str(package / "assemble_absolute_closure.py"),
            "--remap-json", str(remap),
            "--bridge-json", str(bridge),
            "--normalization-json", str(norm),
            "--convergence-json", str(convergence),
            "--steady-json", str(selected),
            "--out", str(closure),
        ], check=True, capture_output=True, text=True)
        result = json.loads(closure.read_text())
        assert result["status"] == "ABSOLUTE_CLOSURE_MEASURED"
        assert result["gates"]["direct_immediate_capture_energy_model_spread"]
        assert result["direct_immediate_capture_energy_spread_fraction_of_mean"] < 0.10
        envelope = result["validated_sensitivity_envelope"]
        assert envelope["mass_current"]["maximum_fractional_sensitivity"] < 0.10
        assert envelope["capture_binding_energy_current"][
            "maximum_fractional_sensitivity"
        ] < 0.10
        for metric in ("mass_current", "capture_binding_energy_current"):
            record = envelope[metric]
            assert record[
                "independent_seed_full_range_fraction_of_mean"
            ] == max(record["independent_seed_by_model"].values())
            by_model = record["resolution_and_boundary"]["by_model"]
            assert record["resolution_and_boundary"]["maximum"] == max(
                value["maximum"] for value in by_model.values()
            )
        assert abs(result["thermal_sink_current_msun_kms2_per_myr"] / (
            -1.5 * result["measured_mdot_msun_per_myr"] * 30.0**2
        ) - 1.0) < 1.0e-14
    print("PASS: production EJ convergence and absolute-closure assembly")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
