#!/usr/bin/env python3
"""Regression tests for the direct boundary-fed FP energy-current gate."""

from __future__ import annotations

import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


PACKAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE))

from gnc_energy_current_gate import evaluate  # noqa: E402


def gate_row(grid, seed, scale=1.0):
    mdot = scale
    return {
        "path": f"grid{grid[0]}_seed{seed}",
        "status": "DIRECT_FP_ENERGY_DIAGNOSTIC_COMPLETE",
        "seed": seed,
        "grid": grid,
        "physical_fingerprint": "same-physical-halo",
        "reservoir_fingerprint": f"reservoir-{grid}",
        "solver_fingerprint": "same-patched-binary",
        "cfs_sha256": "same-cfs-table",
        "validated_build": True,
        "kernel_compatible": True,
        "mass_consistency": 1.0e-6,
        "raw_count": 30,
        "effective_count": 30.0,
        "mdot": mdot,
        "capture_current": 100.0 * scale,
        "advective_current": 40.0 * scale,
        "outward_current": 60.0 * scale,
        "identity_relative": 0.0,
        "nonnegative_release": True,
        "mass_plateau_change": 0.04,
        "energy_plateau_change": 0.05,
        "inventory_completeness": 1.0,
        "inventory_mass_change": 0.03,
        "inventory_energy_change": 0.04,
        "mass_weight_sum": 30.0 * scale,
        "mass_weight_squared_sum": 30.0 * scale * scale,
        "weighted_x_capture_sum": 300.0 * scale,
        "weighted_delta_x_sum": 180.0 * scale,
        "time_myr": 30.0,
        "energy_unit_kms2": 10.0,
        "energy_flux_scale": 1000.0,
        "mass_flux_scale": 100.0,
        "x_boundary": 4.0,
        "tnr_myr": 1.0,
        "burn_in_tnr": 2.0,
        "window_tnr": 1.0,
    }


def write_event_table(path: Path, weighted: bool) -> None:
    with path.open("w") as stream:
        stream.write("Tsnap N_all N_norm N_norm_bd N_emax N_plunge\n")
        for snapshot in range(1, 61):
            capture = 1.0 if weighted and snapshot > 20 else 0.0
            if not weighted:
                capture = 0.5 if snapshot > 20 else 0.0
            stream.write(
                f"{snapshot * 0.1:.6f} 100 100 100 0 {capture:.6f}\n"
            )


def test_synthetic_long_run(root: Path) -> None:
    run = root / "run"
    output = run / "output" / "pro" / "BH"
    output.mkdir(parents=True)
    manifest = {
        "status": "PREFLIGHT_PASS",
        "ranks": 2,
        "particle_mass_msun": 1.0,
        "rh_pc": 10.0,
        "r_in_pc": 1.0,
        "x_boundary": 4.0,
        "rho_boundary_msun_pc3": 2.0,
        "sigma_boundary_kms": 20.0,
        "n0_pc3": 3.0,
        "gnc_seed": 11,
        "same_initialization_seed": True,
        "same_evolution_seed": True,
        "gx_bins": 48,
        "dc_bins": 72,
        "common_reservoir_fingerprint": "reservoir",
        "physical_configuration_fingerprint": "physical",
        "sidm_kernel": "yukawa-tchannel",
        "sigma0_over_m_cm2_g": 100.0,
        "yukawa_w_kms": 80.0,
    }
    (run / "run_manifest.json").write_text(json.dumps(manifest))
    (run / "df_normalization.json").write_text(json.dumps({"mbh_msun": 4.0e6}))
    for executable in ("ini", "main", "pro"):
        (run / executable).write_bytes(("synthetic-" + executable).encode())
    (run / "build_status.json").write_text(json.dumps({
        "status": "COMPILED",
        "compiled": True,
        "absolute_normalization_patch": True,
        "plunge_records_patch": True,
        "inner_inventory_patch": True,
        "snapshot_terminal_coefficients_only": True,
        "weighted_xj_loader": True,
        "born_kernel": True,
    }))
    write_event_table(output / "BH_event_Nweight.txt", True)
    write_event_table(output / "BH_event_N.txt", False)
    for snapshot in range(1, 61):
        inventory = output / f"inner_inventory_{snapshot:04d}.txt"
        inventory.write_text(
            "# isnap weighted_number weighted_mass_msun weighted_mass_x "
            "weighted_mass_x2 raw_count\n"
            f"{snapshot} 100 100 500 3000 200\n"
        )
        plunge = output / f"plunge_records_{snapshot:04d}.txt"
        if snapshot <= 20:
            plunge.write_text("# no post-burn capture\n")
        else:
            # weight_real=2 and ranks=2 gives one physical solar mass.
            plunge.write_text(
                "# isnap x_final x_initial j_final j_initial weight_real "
                "m_msun exit_time rp create_time\n"
                f"{snapshot} 10 2 0.01 0.2 2 1 {snapshot * 0.1} 1e-6 0\n"
            )
    diagnostic = run / "diagnostic.json"
    subprocess.run(
        [
            sys.executable,
            str(PACKAGE / "analyze_long_run.py"),
            str(run),
            "--physical-tnr-myr",
            "1",
            "--kernel-compatible",
            "--out-json",
            str(diagnostic),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    result = json.loads(diagnostic.read_text())
    current = result["boundary_supplied_current"]
    assert result["status"] == "DIRECT_FP_ENERGY_DIAGNOSTIC_COMPLETE"
    assert current["raw_capture_count"] == 40
    assert math.isclose(current["mdot_msun_per_myr"], 10.0)
    assert math.isclose(current["mean_x_capture"], 10.0)
    assert math.isclose(current["mean_delta_x_released"], 6.0)
    assert current["current_identity_relative_residual"] < 1.0e-14
    assert result["inner_reservoir_stationarity"]["stationary"] is True


def main() -> int:
    rows = [
        gate_row((48, 72), 11, 0.98),
        gate_row((48, 72), 22, 1.00),
        gate_row((64, 96), 11, 1.00),
        gate_row((64, 96), 22, 1.02),
    ]
    passed = evaluate(rows)
    assert passed["status"] == "ENERGY_CURRENT_GATE_PASS"
    assert passed["fluid_coupling_authorized"] is True
    assert passed["accepted_closure"] is not None
    assert passed["candidate_closure"]["effective_boundary_supplied_captures"] > 50
    assert passed["candidate_closure"]["binding_current_identity_relative_residual"] < 1e-14

    transient = [dict(row) for row in rows]
    transient[0]["inventory_energy_change"] = 0.6
    failed = evaluate(transient)
    assert failed["status"] == "ENERGY_CURRENT_GATE_FAIL"
    assert failed["gates"]["inner_mass_energy_inventory_stationary"] is False
    assert failed["accepted_closure"] is None

    with tempfile.TemporaryDirectory(prefix="direct-fp-energy-") as tmp:
        test_synthetic_long_run(Path(tmp))
    print("PASS: boundary-fed FP current, live inventory, seed, and grid gates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
