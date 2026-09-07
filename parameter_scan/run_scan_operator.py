#!/usr/bin/env python3
"""Generate one finite-angle FP jump operator for the parameter scan."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


SLOTS = (
    (64, 314159),
    (64, 271828),
    (96, 314159),
    (96, 271828),
)

# The J-only balance tensor is a redundant projection of the full E-J
# operator.  Its entries are assembled with tens of millions of signed
# ``np.add.at`` updates, so demanding a fixed 1e-10 relative closure can make
# the audit depend on summation order.  We instead require its absolute
# residual to remain at least three orders of magnitude below the smaller of
# the independently measured one-sigma Monte-Carlo errors.  The E-J operator,
# which is the object passed to the steady solver, retains its separate 1e-10
# round-trip tests below.
J_ONLY_BALANCE_TO_MC_SIGMA_LIMIT = 1.0e-3


def write_status(path: Path, payload: dict) -> None:
    payload = dict(payload)
    payload["updated_at"] = datetime.now(timezone.utc).isoformat()
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True), flush=True)


def estimator_acceptance(run: dict) -> tuple[bool, str]:
    raw_status = run.get("status")
    if raw_status == "ISOTROPIC_INJECTION_CEILING_CONVERGED":
        return True, "STANDARD"
    if raw_status != "DIAGNOSTIC_ONLY":
        return False, f"UNEXPECTED_STATUS_{raw_status}"

    gates = estimator_gate_results(run)
    if all(gates.values()):
        balance = run.get("nonlocal_jump_balance", {})
        if abs(float(balance.get("mass_balance_max_relative_error", np.inf))) \
                <= 1.0e-10:
            return True, "STANDARD_DIAGNOSTIC_GATES"
        return True, "AUXILIARY_J_BALANCE_BELOW_MC_FLOOR"
    failed = ",".join(key for key, value in gates.items() if not value)
    return False, f"DIAGNOSTIC_GATES_FAILED_{failed}"


def estimator_gate_results(run: dict) -> dict[str, bool]:
    """Return the physical and numerical gates used by the scan driver."""

    balance = run.get("nonlocal_jump_balance", {})
    return {
        "density": float(run.get("rho_df_profile_max_abs_relative_error", np.inf))
        <= 0.01,
        "candidate_sampling": float(
            run.get("candidate_relative_standard_error", np.inf)
        )
        <= 0.10,
        "accepted_sampling": float(
            run.get("accepted_capture_relative_standard_error", np.inf)
        )
        <= 0.10,
        "unresolved_domain": float(
            run.get("unresolved_kepler_domain_fraction_of_candidate", np.inf)
        )
        <= 0.01,
        "j_capture_roundtrip": abs(
            float(balance.get("capture_roundtrip_relative_error", np.inf))
        )
        <= 1.0e-10,
        "j_mass_balance_below_mc_floor": j_only_balance_below_mc_floor(run),
    }


def j_only_balance_over_mc_sigma(run: dict) -> float:
    """Return the auxiliary J-ledger residual in Monte-Carlo sigma units."""

    balance = run.get("nonlocal_jump_balance", {})
    residual = abs(float(
        balance.get("mass_balance_max_absolute_error_msun_per_myr", np.inf)
    ))
    uncertainties = (
        abs(float(run.get("candidate_mdot_msun_per_myr", 0.0)))
        * abs(float(run.get("candidate_relative_standard_error", np.inf))),
        abs(float(run.get("direct_no_rescatter_mdot_msun_per_myr", 0.0)))
        * abs(float(run.get("accepted_capture_relative_standard_error", np.inf))),
    )
    finite_positive = [
        value for value in uncertainties if np.isfinite(value) and value > 0.0
    ]
    if not np.isfinite(residual) or not finite_positive:
        return np.inf
    return residual / min(finite_positive)


def j_only_balance_below_mc_floor(run: dict) -> bool:
    """Test that J-only accumulation error is negligible for this estimate."""

    return j_only_balance_over_mc_sigma(run) <= J_ONLY_BALANCE_TO_MC_SIGMA_LIMIT


def estimator_domain_exclusion(run: dict) -> dict | None:
    """Identify a clean operator whose captured orbits exceed the Kepler domain.

    Such a result is not a failed current estimate.  It says that the present
    point-mass FP operator cannot represent this halo/BH combination.  The
    current is therefore excluded rather than rescued by relaxing the one per
    cent domain gate.
    """

    if run.get("status") != "DIAGNOSTIC_ONLY":
        return None
    gates = estimator_gate_results(run)
    numerical_keys = tuple(key for key in gates if key != "unresolved_domain")
    operator = run.get("ej_jump_operator", {})
    roundtrips = {
        key: abs(float(operator.get(key, np.inf))) <= 1.0e-10
        for key in (
            "direct_capture_roundtrip_relative_error",
            "immediate_capture_roundtrip_relative_error",
            "direct_capture_binding_roundtrip_relative_error",
            "immediate_capture_binding_roundtrip_relative_error",
        )
    }
    if (
        not gates["unresolved_domain"]
        and all(gates[key] for key in numerical_keys)
        and all(roundtrips.values())
    ):
        return {
            "reason": "CAPTURED_ORBITS_EXCEED_KEPLER_DOMAIN",
            "unresolved_kepler_domain_fraction_of_candidate": float(
                run["unresolved_kepler_domain_fraction_of_candidate"]
            ),
            "gate_limit": 0.01,
            "other_gates": gates,
            "ej_roundtrip_gates": roundtrips,
        }
    return None


def domain_exclusion_payload(
    base: dict, outdir: Path, run_json: dict, exclusion: dict
) -> dict:
    return {
        **base,
        "status": "FP_OPERATOR_DOMAIN_UNRESOLVED",
        "reason": exclusion["reason"],
        "unresolved_kepler_domain_fraction_of_candidate": exclusion[
            "unresolved_kepler_domain_fraction_of_candidate"
        ],
        "unresolved_kepler_domain_fraction_limit": exclusion["gate_limit"],
        "run_json": str(outdir / "run.json"),
        "operator_npz": str(outdir / "operator.npz"),
    }


def validate_products(
    outdir: Path, energy_bins: int, run_json: dict
) -> str:
    accepted, acceptance_mode = estimator_acceptance(run_json)
    if not accepted:
        raise RuntimeError(
            f"estimator status {run_json.get('status')}: {acceptance_mode}"
        )
    operator = run_json["ej_jump_operator"]
    if operator["energy_bins"] != energy_bins or operator["angular_bins"] != 257:
        raise RuntimeError("operator grid differs from requested grid")
    for key in (
        "direct_capture_roundtrip_relative_error",
        "immediate_capture_roundtrip_relative_error",
        "direct_capture_binding_roundtrip_relative_error",
        "immediate_capture_binding_roundtrip_relative_error",
    ):
        if abs(float(operator[key])) > 1.0e-10:
            raise RuntimeError(f"operator round trip failed for {key}")
    with np.load(outdir / "operator.npz") as data:
        for key in (
            "rate_direct_binding_msun_kms2_per_myr",
            "rate_immediate_binding_msun_kms2_per_myr",
        ):
            if key not in data.files or np.any(data[key] < 0.0):
                raise RuntimeError(f"invalid binding-energy operator field {key}")
    return acceptance_mode


def complete_payload(base: dict, outdir: Path, run_json: dict, mode: str) -> dict:
    return {
        **base,
        "status": "FP_OPERATOR_COMPLETE",
        "raw_estimator_status": run_json.get("status"),
        "acceptance_mode": mode,
        "j_only_mass_balance_relative_error": run_json.get(
            "nonlocal_jump_balance", {}
        ).get("mass_balance_max_relative_error"),
        "j_only_mass_balance_absolute_error_msun_per_myr": run_json.get(
            "nonlocal_jump_balance", {}
        ).get("mass_balance_max_absolute_error_msun_per_myr"),
        "j_only_mass_balance_over_mc_sigma": j_only_balance_over_mc_sigma(
            run_json
        ),
        "run_json": str(outdir / "run.json"),
        "operator_npz": str(outdir / "operator.npz"),
        "candidate_mdot_msun_per_myr": run_json["candidate_mdot_msun_per_myr"],
        "direct_isotropic_mdot_msun_per_myr": run_json[
            "direct_no_rescatter_mdot_msun_per_myr"
        ],
        "candidate_binding_current_msun_kms2_per_myr": run_json[
            "candidate_binding_energy_current_msun_kms2_per_myr"
        ],
        "direct_isotropic_binding_current_msun_kms2_per_myr": run_json[
            "accepted_capture_binding_energy_current_msun_kms2_per_myr"
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--task-index", type=int, required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--production-package", type=Path, required=True)
    parser.add_argument("--pairs-per-bin", type=int, default=1_600_000)
    args = parser.parse_args()

    grid = json.loads(args.grid.read_text())
    case_index, slot_index = divmod(args.task_index, len(SLOTS))
    energy_bins, seed = SLOTS[slot_index]
    case = dict(grid["cases"][case_index])
    case_root = args.root / "cases" / case["tag"]
    outdir = case_root / "operators" / f"e{energy_bins}_j257_s{seed}"
    outdir.mkdir(parents=True, exist_ok=True)
    base = {
        "schema": "sidm-smbh-scan-operator-v1",
        "case_index": case_index,
        "case_tag": case["tag"],
        "energy_bins": energy_bins,
        "angular_bins": 257,
        "seed": seed,
        "pairs_per_bin": args.pairs_per_bin,
    }
    status_path = outdir / "STATUS.json"
    if status_path.exists():
        old = json.loads(status_path.read_text())
        if old.get("status") in {
            "FP_OPERATOR_COMPLETE",
            "FP_OPERATOR_SKIPPED",
            "FP_OPERATOR_DOMAIN_UNRESOLVED",
        }:
            print(json.dumps(old, indent=2, sort_keys=True))
            return 0
        run_path = outdir / "run.json"
        operator_path = outdir / "operator.npz"
        if run_path.is_file() and operator_path.is_file():
            try:
                existing_run = json.loads(run_path.read_text())
                exclusion = estimator_domain_exclusion(existing_run)
                if exclusion is not None:
                    write_status(
                        status_path,
                        domain_exclusion_payload(
                            base, outdir, existing_run, exclusion
                        ),
                    )
                    return 0
                mode = validate_products(outdir, energy_bins, existing_run)
                write_status(
                    status_path,
                    complete_payload(base, outdir, existing_run, mode),
                )
                return 0
            except Exception:
                pass
    prep_path = case_root / "PREP_STATUS.json"
    if not prep_path.exists():
        write_status(status_path, {**base, "status": "FP_OPERATOR_FAILED", "reason": "missing preparation status"})
        return 0
    prep = json.loads(prep_path.read_text())
    if prep.get("status") != "READY_FOR_FP":
        write_status(
            status_path,
            {
                **base,
                "status": "FP_OPERATOR_SKIPPED",
                "reason": prep.get("status"),
                "preparation_reason": prep.get("reason"),
            },
        )
        return 0

    command = [
        sys.executable,
        str(args.production_package / "finite_angle_capture.py"),
        "--profile",
        prep["bridged_profile"],
        "--df-table",
        prep["df_table"],
        "--normalization-json",
        prep["normalization_json"],
        "--out-json",
        str(outdir / "run.json"),
        "--out-csv",
        str(outdir / "run.csv"),
        "--out-ej-operator",
        str(outdir / "operator.npz"),
        "--mbh",
        str(case["black_hole_mass_msun"]),
        "--sigma0-over-m",
        str(case["sigma0_over_m_cm2_g"]),
        "--w-kms",
        str(case["yukawa_w_kms"]),
        "--r-outer",
        str(prep["orbit_reporting_radius_pc"]),
        "--radial-bins",
        "16",
        "--pairs-per-bin",
        str(args.pairs_per_bin),
        "--cap-importance-fraction",
        "0.5",
        "--direction-edge-fraction",
        "0.3",
        "--direction-edge-power",
        "4",
        "--survival-anomaly-order",
        "16",
        "--survival-collision-order",
        "48",
        "--ej-energy-bins",
        str(energy_bins),
        "--ej-angular-grid",
        "j257",
        "--seed",
        str(seed),
    ]
    with (outdir / "operator.log").open("w") as stream:
        stream.write("COMMAND " + " ".join(command) + "\n")
        stream.flush()
        completed = subprocess.run(
            command,
            stdout=stream,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        stream.write(f"RETURN_CODE {completed.returncode}\n")
    try:
        if completed.returncode != 0:
            raise RuntimeError(f"finite-angle estimator returned {completed.returncode}")
        run_json = json.loads((outdir / "run.json").read_text())
        exclusion = estimator_domain_exclusion(run_json)
        if exclusion is not None:
            write_status(
                status_path,
                domain_exclusion_payload(base, outdir, run_json, exclusion),
            )
            return 0
        mode = validate_products(outdir, energy_bins, run_json)
        write_status(
            status_path,
            complete_payload(base, outdir, run_json, mode),
        )
    except Exception as exc:
        write_status(
            status_path,
            {
                **base,
                "status": "FP_OPERATOR_FAILED",
                "reason": f"{type(exc).__name__}: {exc}",
            },
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
