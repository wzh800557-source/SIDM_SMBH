"""Fail-closed loading of a validated fluid comparison for publication.

This verifies the numerical response record, not the physical energy closure.
Test fixtures and historical *_old files must never be used as production data.
"""
from pathlib import Path
import csv
import hashlib
import json
import math


def load_validated_response(directory):
    directory = Path(directory)
    analysis_path = directory / "fluid_response_analysis.json"
    comparison_path = directory / "fluid_response_comparison.csv"
    for path in (analysis_path, comparison_path):
        if not path.is_file():
            raise RuntimeError(f"Missing validated response input: {path}. No historical fallback is allowed.")
    analysis = json.loads(analysis_path.read_text())
    gates = analysis.get("gates", {})
    if analysis.get("status") != "PASS" or not gates or any(v is not True for v in gates.values()):
        raise RuntimeError("Fluid-response validation is not PASS with all required gates satisfied")
    if analysis.get("is_fixture") or analysis.get("synthetic"):
        raise RuntimeError("Synthetic inputs cannot be published as a production response")
    with comparison_path.open(newline="") as stream:
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]
    if len(rows) < 2 or any(not math.isfinite(v) for row in rows for v in row.values()):
        raise RuntimeError("Fluid comparison is incomplete or non-finite")
    tau = [row["tau_relax"] for row in rows]
    if any(b <= a for a, b in zip(tau, tau[1:])):
        raise RuntimeError("Fluid comparison times are not strictly increasing")
    comparisons = {
        "common_tau_relax": "tau_relax",
        "inner_mean_density_sink_over_control_at_common_end": "sink_over_control_rho_inner_mean_msun_pc3",
        "inner_radius_sink_over_control_at_common_end": "sink_over_control_r_inner_pc",
        "inner_dispersion_sink_over_control_at_common_end": "sink_over_control_sigma_inner_1d_kms",
    }
    for ledger_key, csv_key in comparisons.items():
        if not math.isclose(float(analysis[ledger_key]), rows[-1][csv_key], rel_tol=1e-10, abs_tol=1e-12):
            raise RuntimeError(f"Ledger and comparison endpoints disagree: {ledger_key}")
    expected = analysis.get("identity", {}).get("fluid_comparison_csv_sha256")
    if not expected or hashlib.sha256(comparison_path.read_bytes()).hexdigest() != expected:
        raise RuntimeError("Comparison table does not match the validated input identity")
    return analysis, rows


def measured_response_endpoint(analysis):
    """Require actual common-endpoint values, never infer them from a stop limit."""
    names = ("common_time_myr", "captured_mass_at_common_end_msun")
    values = []
    for name in names:
        value = analysis.get(name)
        if value is None or not math.isfinite(float(value)) or float(value) <= 0:
            raise RuntimeError(f"Missing measured common-endpoint quantity: {name}")
        values.append(float(value))
    return tuple(values)
