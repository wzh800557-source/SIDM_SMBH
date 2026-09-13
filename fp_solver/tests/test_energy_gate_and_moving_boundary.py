"""Gate and conservation tests for the production moving-boundary path."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from conservative_moving_boundary import evolve  # noqa: E402
from energy_current_gate import evaluate  # noqa: E402
from moving_interface_ledger import ReservoirState  # noqa: E402


def energy_row(energy_bins, seed, energy, mdot, passed=True):
    return {
        "energy_bins": energy_bins,
        "seed": seed,
        "model": "direct",
        "mdot_msun_per_myr": mdot,
        "energy_audit": {
            "capture_orbital_energy_rate": energy,
            "internal_l1_over_boundary_scale": 0.005 if passed else 2.0,
            "status": "MOMENT_AUDIT_PASS" if passed else "MOMENT_AUDIT_FAIL",
            "gates": {
                "physical_pair_energy_conserved": True,
                "ledger_identity": True,
            },
        },
    }


def main():
    summary = {
        "status": "FOLLOWUP_DIAGNOSTICS_COMPLETE",
        "missing_runs": [],
        "fluid_coupling_authorized": True,
        "rows": [
            energy_row(64, 1, -100.0, 10.0),
            energy_row(64, 2, -101.0, 10.1),
            energy_row(96, 1, -102.0, 10.2),
            energy_row(96, 2, -103.0, 10.3),
        ],
    }
    passed = evaluate(summary)
    assert passed["fluid_coupling_authorized"] is True
    summary["rows"][0]["energy_audit"]["internal_l1_over_boundary_scale"] = 2.7
    summary["rows"][0]["energy_audit"]["status"] = "MOMENT_AUDIT_FAIL"
    failed = evaluate(summary)
    assert failed["fluid_coupling_authorized"] is False
    assert failed["gates"]["internal_energy_moments_stationary"] is False

    initial = ReservoirState(100.0, 10.0, 1.0e6, -300.0, -80.0)
    base = {
        "rho_boundary_msun_pc3": 2.0,
        "capture_mdot_msun_per_myr": 0.1,
        "interface_specific_energy_kms2": -2.0,
        "captured_specific_energy_kms2": -7.0,
        "outward_luminosity_msun_kms2_per_myr": 0.2,
        "work_on_system_msun_kms2_per_myr": 0.05,
    }
    rows = [
        dict(base, time_myr=0.0, r_boundary_pc=1.0),
        dict(base, time_myr=0.1, r_boundary_pc=1.01),
        dict(base, time_myr=0.2, r_boundary_pc=0.99),
    ]
    final, trajectory, diagnostic = evolve(initial, rows)
    assert len(trajectory) == 2
    assert diagnostic["mass_conserved"]
    assert diagnostic["energy_conserved_after_work"]
    assert abs(final.tracked_mass - initial.tracked_mass) < 1.0e-9
    print("PASS: capture-energy gate and moving-boundary conservation")


if __name__ == "__main__":
    main()
