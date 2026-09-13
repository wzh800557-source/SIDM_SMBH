#!/usr/bin/env python3
"""Plot the calibrated conductive benchmark and the fiducial Yukawa regime audit."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "fp_solver"))

from cross_regime_closure import ConductiveSpikeCalibration  # noqa: E402


def read_points(path: Path):
    rows = list(csv.DictReader(path.open()))
    return (
        np.asarray([float(row["sigma_over_m_cm2_g"]) for row in rows]),
        np.asarray([float(row["accretion_rate_msun_per_yr"]) for row in rows]),
    )


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--calibration-json",
        type=Path,
        default=ROOT / "validation/cross_regime_calibration.json",
    )
    p.add_argument(
        "--digitized-csv",
        type=Path,
        default=ROOT / "data/calibration/sabarish_2025_imfp_digitized.csv",
    )
    p.add_argument(
        "--bridge-profile",
        type=Path,
        default=ROOT / "data/production/bridge/bridged_profile.txt",
    )
    p.add_argument("--outdir", type=Path, required=True)
    args = p.parse_args()

    result = json.loads(args.calibration_json.read_text())
    points_s, points_rate = read_points(args.digitized_csv)
    refit = result["conductive_spike"]["digitized_same_form_refit"]
    calibration = ConductiveSpikeCalibration(
        A=refit["A"],
        B=refit["B"],
        baseline_msun_per_yr=refit["baseline_msun_per_yr"],
    )
    profile = np.loadtxt(args.bridge_profile, comments="#", ndmin=2)
    radius = profile[:, 0]
    n_fiducial = profile[:, 4]
    fiducial = result["fiducial_yukawa_halo"]

    plt.rcParams.update(
        {
            "font.size": 9.0,
            "axes.labelsize": 9.5,
            "legend.fontsize": 7.6,
            "xtick.labelsize": 8.2,
            "ytick.labelsize": 8.2,
            "axes.linewidth": 0.8,
            "savefig.bbox": "tight",
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(7.15, 3.05))

    ax = axes[0]
    s_grid = np.geomspace(0.015, 12.0, 500)
    model_rate = calibration.total_rate_msun_per_yr(s_grid)
    ax.plot(s_grid, model_rate, color="#1f4e79", lw=2.0, label="harmonic-form refit")
    ax.scatter(
        points_s,
        points_rate,
        s=23,
        color="#d95f02",
        edgecolor="white",
        linewidth=0.5,
        zorder=3,
        label="digitized simulation points",
    )
    transition = calibration.sigma_transition_cm2_g
    ax.axvline(transition, color="0.35", lw=1.0, ls="--")
    ax.text(
        transition * 1.08,
        5.82,
        rf"$s_\ast={transition:.2f}$",
        va="top",
        color="0.25",
    )
    ax.text(
        0.038,
        2.72,
        r"LMFP: $\Delta\dot M\propto s$",
        color="#1f4e79",
        fontsize=8.2,
    )
    ax.text(
        4.6,
        2.72,
        r"SMFP: $\Delta\dot M\propto s^{-1}$",
        color="#1f4e79",
        fontsize=8.0,
        ha="center",
    )
    ax.set_xscale("log")
    ax.set_xlim(0.015, 12.0)
    ax.set_ylim(1.45, 6.05)
    ax.set_xlabel(r"constant cross-section $s=\sigma_\chi/m$ [$\mathrm{cm^2\,g^{-1}}$]")
    ax.set_ylabel(r"capture rate [$\mathrm{M_\odot\,yr^{-1}}$]")
    ax.legend(loc="upper left", frameon=False, bbox_to_anchor=(0.035, 1.0))
    ax.text(-0.055, 1.035, "(a)", transform=ax.transAxes, va="top", fontweight="bold")

    ax = axes[1]
    colors = {10.0: "#2ca25f", 100.0: "#2b8cbe", 1000.0: "#88419d"}
    mask = radius <= max(3.0, fiducial["nominal_bondi_radius_pc"] * 1.2)
    for sigma0 in (10.0, 100.0, 1000.0):
        n_curve = n_fiducial * sigma0 / 100.0
        ax.plot(
            radius[mask],
            n_curve[mask],
            color=colors[sigma0],
            lw=1.7,
            label=rf"$\sigma_{{\chi,0}}/m={sigma0:g}$",
        )
    ax.axhspan(1.0e-10, 1.0, color="#edf8fb", zorder=-5)
    ax.axhline(1.0, color="0.15", lw=1.0, ls="--")
    ax.text(2.0e-6, 2.5e-7, "orbit memory retained", color="#24566b")
    ax.axvline(fiducial["boundary_radius_pc"], color="0.45", lw=0.9, ls=":")
    ax.axvline(fiducial["nominal_bondi_radius_pc"], color="0.15", lw=0.9, ls="-.")
    ax.text(
        fiducial["boundary_radius_pc"] * 0.85,
        1.8e2,
        r"$r_{\rm in}$",
        ha="right",
        va="center",
        color="0.35",
    )
    ax.text(
        fiducial["nominal_bondi_radius_pc"] * 1.08,
        4.0e1,
        r"$r_{\rm B}$",
        ha="left",
        va="center",
        color="0.15",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(radius[0], 3.0)
    ax.set_ylim(1.0e-10, 5.0e2)
    ax.set_xlabel("radius [pc]")
    ax.set_ylabel(r"orbital transport depth $N_{\rm orb}$")
    legend = ax.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.92)
    legend.get_frame().set_linewidth(0.0)
    ax.text(-0.055, 1.035, "(b)", transform=ax.transAxes, va="top", fontweight="bold")
    ax.text(
        0.98,
        0.035,
        r"fiducial Yukawa profile, $\sigma_{\chi,0}/m=100$"
        + "\n"
        + rf"$N_{{\rm orb}}(r_{{\rm ISCO}})={fiducial['radial_collisionality_gate']['N_orb_at_capture']:.2g}$",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="0.25",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.86, "pad": 1.6},
    )

    for ax in axes:
        ax.tick_params(direction="in", which="both", top=True, right=True)
    fig.subplots_adjust(wspace=0.30)
    args.outdir.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "png"):
        fig.savefig(args.outdir / f"fig_cross_regime_calibration.{suffix}", dpi=300)
    plt.close(fig)
    print("CROSS_REGIME_CALIBRATION_FIGURE_COMPLETE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
