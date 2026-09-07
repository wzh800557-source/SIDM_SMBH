#!/usr/bin/env python3
"""Plot accepted production E-J convergence and black-hole fluid response."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from pathlib import Path
from typing import Optional

os.environ.setdefault("MPLBACKEND", "Agg")
_platform = sys.platform
if _platform == "darwin":
    sys.platform = "linux"
import matplotlib.pyplot as plt
sys.platform = _platform
import numpy as np


MODELS = (("direct", "direct", "#275d95", "o"),
          ("immediate", "immediate", "#cc6b32", "s"))
METRICS = {
    "mass": ("mass current /\nfinal direct mass current", "mass"),
    "capture_binding_energy": (
        "binding current /\nfinal direct binding current", "binding energy"
    ),
}


def load_json(path: Path, status: str, schema: Optional[str] = None) -> dict:
    value = json.loads(path.read_text())
    if value.get("status") != status:
        raise RuntimeError(f"{path} has status {value.get('status')!r}, not {status}")
    if schema is not None and value.get("schema") != schema:
        raise RuntimeError(f"{path} has schema {value.get('schema')!r}, not {schema}")
    return value


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(a: object, b: object, tolerance: float = 1.0e-10) -> bool:
    return math.isclose(
        float(a), float(b), rel_tol=tolerance, abs_tol=1.0e-12
    )


def require_all_true(mapping: dict, label: str) -> None:
    if not isinstance(mapping, dict) or not mapping:
        raise RuntimeError(f"{label} has no recorded gates")
    failed = [key for key, value in mapping.items() if value is not True]
    if failed:
        raise RuntimeError(f"{label} has failed gates: {', '.join(failed)}")


def validate_convergence(convergence: dict, closure_scope: str) -> None:
    tolerance = float(convergence.get("fractional_tolerance", math.nan))
    if not math.isfinite(tolerance) or not 0.0 < tolerance <= 0.1:
        raise RuntimeError("convergence ledger has an invalid production tolerance")
    expected = {
        (model, metric)
        for model, _label, _colour, _marker in MODELS
        for metric in METRICS
    }
    required = set(convergence.get("required_dimensions", ()))
    if required != {"energy", "angular", "boundary"}:
        raise RuntimeError("convergence ledger did not require all three dimensions")
    if set(convergence.get("required_metrics", ())) != set(METRICS):
        raise RuntimeError("convergence ledger did not require both currents")
    if closure_scope == "mass_current_closure_only":
        mass_gate = convergence.get("observable_gates", {}).get("mass", {})
        binding_gate = convergence.get("observable_gates", {}).get(
            "capture_binding_energy", {}
        )
        if (
            mass_gate.get("status")
            != "PRODUCTION_OBSERVABLE_CONVERGENCE_PASS"
            or mass_gate.get("pass") is not True
        ):
            raise RuntimeError("mass current did not pass its observable gate")
        if binding_gate.get("pass") not in (True, False):
            raise RuntimeError("binding-energy observable status is absent")
    for dimension in ("energy", "angular", "boundary"):
        gate = convergence.get("dimension_gates", {}).get(dimension, {})
        if closure_scope == "absolute_two_current_closure" and gate.get(
            "pass"
        ) is not True:
            raise RuntimeError(f"convergence dimension {dimension} did not pass")
        records = convergence.get("selected_comparisons", {}).get(dimension, [])
        observed = {(item.get("model"), item.get("metric")) for item in records}
        if len(records) != len(expected) or observed != expected:
            raise RuntimeError(
                f"convergence dimension {dimension} lacks the four selected tests"
            )
        required_records = (
            records
            if closure_scope == "absolute_two_current_closure"
            else [item for item in records if item.get("metric") == "mass"]
        )
        if any(item.get("gate") is not True for item in required_records):
            raise RuntimeError(
                f"required {dimension} comparison did not pass"
            )


def validate_closure(closure: dict, tolerance: float) -> str:
    absolute = (
        closure.get("schema") == "gnc-fluid-absolute-closure-v1"
        and closure.get("status") == "ABSOLUTE_CLOSURE_MEASURED"
    )
    mass_only = (
        closure.get("schema") == "gnc-fluid-mass-closure-v1"
        and closure.get("status") == "FLUID_MASS_CLOSURE_MEASURED"
        and closure.get("capture_binding_energy_current_used_by_fluid") is False
    )
    if not (absolute or mass_only):
        raise RuntimeError("closure is neither an accepted absolute nor mass closure")
    require_all_true(closure.get("gates", {}), "accepted closure")
    mdot = float(closure.get("measured_mdot_msun_per_myr", math.nan))
    cm = float(closure.get("C_M_measured", math.nan))
    ce = float(closure.get("C_E_thermal_sink", math.nan))
    factor = float(closure.get("thermal_specific_energy_factor", math.nan))
    if not all(math.isfinite(value) for value in (mdot, cm, ce, factor)):
        raise RuntimeError("closure currents contain a non-finite value")
    if mdot <= 0.0 or cm <= 0.0 or ce >= 0.0 or factor <= 0.0:
        raise RuntimeError("closure currents have an unphysical sign")
    if not close(ce, -factor * cm):
        raise RuntimeError("thermal sink coefficient is inconsistent with mass current")
    if absolute:
        envelopes = closure.get("validated_sensitivity_envelope", {})
        records = (
            envelopes.get("mass_current", {}),
            envelopes.get("capture_binding_energy_current", {}),
        )
    else:
        records = (closure.get("validated_mass_sensitivity_envelope", {}),)
    for record in records:
        value = float(
            record.get("maximum_fractional_sensitivity", math.nan)
        )
        if not math.isfinite(value) or value < 0.0 or value > tolerance:
            raise RuntimeError("closure current exceeds the accepted envelope")
    return (
        "absolute_two_current_closure" if absolute else "mass_current_closure_only"
    )


def primary_groups(convergence: dict) -> list[dict]:
    radius = float(
        convergence["production_resolution_path"]["primary_boundary_radius_pc"]
    )
    groups = [
        group for group in convergence["groups"]
        if np.isclose(float(group["reservoir_radius_pc"]), radius, rtol=1.0e-10)
    ]
    groups.sort(key=lambda group: (int(group["energy_bins"]), int(group["angular_bins"])))
    if not groups:
        raise RuntimeError("convergence ledger has no primary-boundary groups")
    return groups


def production_path_groups(convergence: dict) -> list[dict]:
    """Return only the three adjacent grids used by the accepted gates."""

    path = convergence["production_resolution_path"]
    keys = [
        (int(path["previous_energy_bins"]),
         int(path["energy_ladder_angular_bins"])),
        (int(path["final_energy_bins"]),
         int(path["energy_ladder_angular_bins"])),
        (int(path["final_energy_bins"]), int(path["final_angular_bins"])),
    ]
    if len(set(keys)) != 3:
        raise RuntimeError("accepted convergence path does not contain three grids")
    groups = primary_groups(convergence)
    selected = []
    for energy, angular in keys:
        matches = [
            group for group in groups
            if int(group["energy_bins"]) == energy
            and int(group["angular_bins"]) == angular
        ]
        if len(matches) != 1:
            raise RuntimeError(
                f"cannot identify accepted path grid {energy}/{angular}"
            )
        selected.append(matches[0])
    return selected


def metric(group: dict, model: str, name: str) -> tuple[float, float, float]:
    record = group["models"][model][name]
    combined = float(record["combined_value"])
    seeds = np.asarray(record["seed_values"], dtype=float)
    if (
        not math.isfinite(combined) or combined <= 0.0
        or seeds.size < 2 or np.any(~np.isfinite(seeds)) or np.any(seeds <= 0.0)
    ):
        raise RuntimeError("production current estimator is non-positive or incomplete")
    return combined, float(np.min(seeds)), float(np.max(seeds))


def save(fig: plt.Figure, stem: Path) -> list[str]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    png = stem.with_suffix(".png")
    pdf = stem.with_suffix(".pdf")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return [png.name, pdf.name]


def plot_convergence(convergence: dict, closure: dict, outdir: Path) -> list[str]:
    groups = production_path_groups(convergence)
    final_energy = int(
        convergence["production_resolution_path"]["final_energy_bins"]
    )
    final_angular = int(
        convergence["production_resolution_path"]["final_angular_bins"]
    )
    selected = [
        group for group in groups
        if int(group["energy_bins"]) == final_energy
        and int(group["angular_bins"]) == final_angular
    ]
    if len(selected) != 1:
        raise RuntimeError("cannot identify the unique final production group")

    plt.rcParams.update({
        "font.family": "serif", "font.size": 9, "axes.labelsize": 9,
        "axes.titlesize": 9, "legend.fontsize": 7.5,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.linewidth": 0.8,
    })
    fig, axes = plt.subplots(1, 3, figsize=(7.15, 3.05), constrained_layout=True)
    x = np.arange(len(groups), dtype=float)
    labels = [
        f"{int(group['energy_bins'])}/{int(group['angular_bins'])}"
        for group in groups
    ]
    for ax, (name, (ylabel, _title)) in zip(axes[:2], METRICS.items()):
        normalizer = metric(selected[0], "direct", name)[0]
        for model, label, colour, marker in MODELS:
            centres, lower, upper = [], [], []
            for group in groups:
                centre, lo, hi = metric(group, model, name)
                centres.append(centre / normalizer)
                lower.append((centre - min(centre, lo)) / normalizer)
                upper.append((max(centre, hi) - centre) / normalizer)
            ax.errorbar(
                x, centres, yerr=np.asarray([lower, upper]), marker=marker,
                ms=3.8, lw=1.0, capsize=2.0, color=colour, label=label,
            )
        ax.axhspan(0.9, 1.1, color="#7aa974", alpha=0.15, lw=0)
        ax.axhline(1.0, color="0.2", lw=0.8, ls=":")
        ax.set_xticks(x, labels, rotation=42, ha="right")
        ax.set_xlabel(r"energy/angular bins, $N_E/N_J$")
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.20)
    mass_only = closure.get("schema") == "gnc-fluid-mass-closure-v1"
    axes[0].set_title("(a) Mass current\n(accepted)", loc="left", fontsize=9.0)
    axes[1].set_title(
        "(b) Binding current\n"
        + ("(unresolved)" if mass_only else "(accepted)"),
        loc="left", fontsize=9.0,
    )
    axes[0].legend(frameon=False, loc="best")

    ax = axes[2]
    colours = {"mass": "#275d95", "capture_binding_energy": "#4c956c"}
    markers = {"mass": "o", "capture_binding_energy": "s"}
    labels_by_metric = {"mass": "mass", "capture_binding_energy": "binding energy"}
    primary = float(
        convergence["production_resolution_path"]["primary_boundary_radius_pc"]
    )
    for name in ("mass", "capture_binding_energy"):
        matches = [
            item for item in convergence["selected_comparisons"]["boundary"]
            if item["model"] == "direct" and item["metric"] == name
        ]
        if len(matches) != 1:
            raise RuntimeError(f"missing selected direct boundary comparison for {name}")
        item = matches[0]
        radii = np.asarray(item["radii_pc"], dtype=float)
        values = np.asarray(item["values"], dtype=float)
        if (
            radii.size < 3 or values.size != radii.size
            or np.any(~np.isfinite(radii)) or np.any(~np.isfinite(values))
            or np.any(radii <= 0.0) or np.any(values <= 0.0)
            or np.any(np.diff(radii) <= 0.0)
        ):
            raise RuntimeError(f"invalid boundary comparison for {name}")
        at_primary = np.flatnonzero(
            np.isclose(radii, primary, rtol=1.0e-10, atol=0.0)
        )
        if at_primary.size != 1:
            raise RuntimeError(f"boundary comparison for {name} lacks r_in")
        reference = values[int(at_primary[0])]
        ax.plot(
            radii / primary, values / reference, marker=markers[name],
            ms=4.0, lw=1.0,
            color=colours[name], label=labels_by_metric[name],
        )
    ax.axhspan(0.9, 1.1, color="#7aa974", alpha=0.15, lw=0)
    ax.axhline(1.0, color="0.2", lw=0.8, ls=":")
    ax.axvline(1.0, color="0.5", lw=0.7, ls="--")
    ax.set_xlabel(r"reporting radius / $r_{\rm in}$")
    ax.set_ylabel(r"current / current at $r_{\rm in}$")
    ax.set_title("(c) Reporting surface\n(final grid)", loc="left", fontsize=9.0)
    ax.legend(frameon=False, loc="upper left")
    ax.grid(axis="y", alpha=0.20)

    if mass_only:
        mass_env = 100.0 * float(
            closure["validated_mass_sensitivity_envelope"][
                "maximum_fractional_sensitivity"
            ]
        )
        energy_changes = [
            float(item["full_range_fraction_of_mean"] if dimension == "boundary"
                  else item["fractional_difference"])
            for dimension in ("energy", "angular", "boundary")
            for item in convergence["selected_comparisons"][dimension]
            if item["metric"] == "capture_binding_energy"
        ]
        binding_pass = convergence["observable_gates"]["capture_binding_energy"]["pass"]
        energy_text = (f"binding energy: {100.0 * max(energy_changes):.1f}% "
                       f"({'passes' if binding_pass else 'fails'})")
    else:
        mass_env = 100.0 * float(
            closure["validated_sensitivity_envelope"]["mass_current"][
                "maximum_fractional_sensitivity"
            ]
        )
        energy_env = 100.0 * float(
            closure["validated_sensitivity_envelope"][
                "capture_binding_energy_current"
            ]["maximum_fractional_sensitivity"]
        )
        energy_text = f"binding: {energy_env:.1f}%"
    axes[2].text(
        0.04, 0.05,
        f"mass envelope: {mass_env:.1f}%\n{energy_text}",
        transform=axes[2].transAxes, fontsize=7.5, va="bottom",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82,
              "pad": 1.5},
    )
    return save(fig, outdir / "fig_production_ej_convergence")


def read_csv(path: Path) -> list[dict[str, float]]:
    with path.open(newline="") as stream:
        return [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(stream)
        ]


def validate_response_rows(analysis: dict, rows: list[dict[str, float]]) -> None:
    if len(rows) < 2:
        raise RuntimeError("fluid comparison needs at least two output times")
    required = {
        "tau_relax",
        "control_rho_inner_mean_msun_pc3", "sink_rho_inner_mean_msun_pc3",
        "control_sigma_inner_1d_kms", "sink_sigma_inner_1d_kms",
        "control_r_inner_pc", "sink_r_inner_pc",
        "sink_over_control_rho_inner_mean_msun_pc3",
        "sink_over_control_sigma_inner_1d_kms",
        "sink_over_control_r_inner_pc",
    }
    for row in rows:
        if not required <= row.keys():
            raise RuntimeError("fluid comparison is missing a required column")
        if any(not math.isfinite(float(row[key])) for key in required):
            raise RuntimeError("fluid comparison contains a non-finite value")
        for field in (
            "rho_inner_mean_msun_pc3", "sigma_inner_1d_kms", "r_inner_pc"
        ):
            control = float(row[f"control_{field}"])
            sink = float(row[f"sink_{field}"])
            ratio = float(row[f"sink_over_control_{field}"])
            if control <= 0.0 or sink <= 0.0 or not close(ratio, sink / control):
                raise RuntimeError("fluid comparison ratio is inconsistent")
    tau = [float(row["tau_relax"]) for row in rows]
    if tau[0] != 0.0 or any(right <= left for left, right in zip(tau, tau[1:])):
        raise RuntimeError("fluid comparison times are not strictly increasing")
    if float(analysis["common_tau_relax"]) <= 0.0:
        raise RuntimeError("fluid comparison has no evolved common interval")
    if not close(tau[-1], analysis["common_tau_relax"]):
        raise RuntimeError("fluid comparison does not end at the common time")
    final = rows[-1]
    expected_endpoints = {
        "inner_mean_density_sink_over_control_at_common_end":
            final["sink_over_control_rho_inner_mean_msun_pc3"],
        "inner_dispersion_sink_over_control_at_common_end":
            final["sink_over_control_sigma_inner_1d_kms"],
        "inner_radius_sink_over_control_at_common_end":
            final["sink_over_control_r_inner_pc"],
    }
    for key, value in expected_endpoints.items():
        if not close(analysis[key], value):
            raise RuntimeError(f"fluid analysis does not match {key}")
    density_max = max(
        abs(float(row["sink_over_control_rho_inner_mean_msun_pc3"]) - 1.0)
        for row in rows
    )
    if not close(
        analysis["inner_mean_density_max_abs_fractional_difference"], density_max
    ):
        raise RuntimeError("fluid analysis does not match the density history")


def plot_response(analysis: dict, rows: list[dict[str, float]], outdir: Path) -> list[str]:
    validate_response_rows(analysis, rows)
    tau = np.asarray([row["tau_relax"] for row in rows])
    control_rho = np.asarray([
        row["control_rho_inner_mean_msun_pc3"] for row in rows
    ])
    sink_rho = np.asarray([
        row["sink_rho_inner_mean_msun_pc3"] for row in rows
    ])
    fig, axes = plt.subplots(1, 2, figsize=(7.15, 2.75), constrained_layout=True)
    axes[0].semilogy(tau, control_rho, color="#444444", lw=1.5, label="control")
    axes[0].semilogy(tau, sink_rho, color="#b44b3c", lw=1.5, label="capture sink")
    axes[0].set_xlabel(r"time [$t_{\rm relax}$]")
    axes[0].set_ylabel(r"mean density inside $r_0$ [$M_\odot\,{\rm pc}^{-3}$]")
    axes[0].set_title("(a) Innermost resolved density", loc="left")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.20)

    fields = (
        (
            "sink_over_control_rho_inner_mean_msun_pc3",
            r"$\bar{\rho}(<r_0)$", "#b44b3c",
        ),
        (
            "sink_over_control_sigma_inner_1d_kms",
            r"$\sigma_{1{\rm D},0}$", "#275d95",
        ),
        ("sink_over_control_r_inner_pc", r"$r_0$", "#4c956c"),
    )
    for field, label, colour in fields:
        axes[1].plot(
            tau, [row[field] for row in rows], lw=1.4, color=colour, label=label,
        )
    axes[1].axhline(1.0, color="0.25", lw=0.8, ls=":")
    axes[1].set_xlabel(r"time [$t_{\rm relax}$]")
    axes[1].set_ylabel("capture sink / control")
    axes[1].set_title("(b) Black-hole-aware response", loc="left")
    axes[1].legend(frameon=False)
    axes[1].grid(alpha=0.20)
    density_max = 100.0 * float(
        analysis["inner_mean_density_max_abs_fractional_difference"]
    )
    axes[1].text(
        0.04, 0.06, f"maximum inner-density difference: {density_max:.2f}%",
        transform=axes[1].transAxes, fontsize=7.7,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82,
              "pad": 1.5},
    )
    return save(fig, outdir / "fig_production_fluid_response")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--convergence-json", type=Path, required=True)
    parser.add_argument("--closure-json", type=Path, required=True)
    parser.add_argument("--fluid-analysis-json", type=Path, required=True)
    parser.add_argument("--fluid-comparison-csv", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    convergence = json.loads(args.convergence_json.read_text())
    if convergence.get("schema") != "finite-angle-ej-convergence-v3":
        raise RuntimeError("convergence ledger schema is stale")
    closure = json.loads(args.closure_json.read_text())
    closure_scope = validate_closure(
        closure, float(convergence.get("fractional_tolerance", math.nan))
    )
    allowed_convergence_status = (
        {"PRODUCTION_CONVERGENCE_PASS"}
        if closure_scope == "absolute_two_current_closure"
        else {"PRODUCTION_CONVERGENCE_PASS", "INCOMPLETE_OR_FAILED"}
    )
    if convergence.get("status") not in allowed_convergence_status:
        raise RuntimeError(
            "convergence status is inconsistent with the closure scope"
        )
    analysis = load_json(
        args.fluid_analysis_json,
        "BLACK_HOLE_AWARE_RESPONSE_COMPLETE",
        "black-hole-aware-fluid-response-v3",
    )
    rows = read_csv(args.fluid_comparison_csv)
    require_all_true(analysis.get("gates", {}), "fluid response")
    validate_convergence(convergence, closure_scope)
    if analysis.get("closure_scope") != closure_scope:
        raise RuntimeError("fluid analysis and closure report different scopes")
    if (
        closure_scope == "mass_current_closure_only"
        and analysis.get("capture_binding_energy_current_used_by_fluid") is not False
    ):
        raise RuntimeError(
            "mass-current figure input uses the unresolved binding-energy current"
        )
    closure_identity = closure.get("identity", {})
    path = convergence["production_resolution_path"]
    if closure_identity.get("convergence_json_sha256") != sha256(
        args.convergence_json
    ):
        raise RuntimeError("closure does not match the convergence ledger")
    if (
        int(closure_identity.get("energy_bins", -1))
        != int(path["final_energy_bins"])
        or int(closure_identity.get("angular_bins", -1))
        != int(path["final_angular_bins"])
    ):
        raise RuntimeError("closure does not use the final phase-space grid")
    analysis_identity = analysis.get("identity", {})
    if analysis_identity.get("closure_json_sha256") != sha256(args.closure_json):
        raise RuntimeError("fluid analysis does not match the closure ledger")
    if analysis_identity.get("fluid_comparison_csv_sha256") != sha256(
        args.fluid_comparison_csv
    ):
        raise RuntimeError("fluid analysis does not match the comparison table")
    validate_response_rows(analysis, rows)
    outputs = plot_convergence(convergence, closure, args.outdir)
    outputs.extend(plot_response(analysis, rows, args.outdir))
    report = {
        "schema": "gnc-production-figure-ledger-v3",
        "status": "PRODUCTION_FIGURES_COMPLETE",
        "closure_scope": closure_scope,
        "capture_binding_energy_current_accepted": (
            closure_scope == "absolute_two_current_closure"
        ),
        "sources": {
            "convergence_sha256": sha256(args.convergence_json),
            "closure_sha256": sha256(args.closure_json),
            "fluid_analysis_sha256": sha256(args.fluid_analysis_json),
            "fluid_comparison_sha256": sha256(args.fluid_comparison_csv),
        },
        "outputs": outputs,
        "output_sha256": {
            name: sha256(args.outdir / name) for name in outputs
        },
        "measured_mdot_msun_per_myr": closure["measured_mdot_msun_per_myr"],
        "C_M_measured": closure["C_M_measured"],
        "C_E_thermal_sink": closure["C_E_thermal_sink"],
        "inner_mean_density_max_abs_fractional_difference": analysis[
            "inner_mean_density_max_abs_fractional_difference"
        ],
        "plotted_resolution_path": [
            f"{int(group['energy_bins'])}/{int(group['angular_bins'])}"
            for group in production_path_groups(convergence)
        ],
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "production_figure_ledger.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
