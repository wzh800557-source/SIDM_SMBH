#!/usr/bin/env python3
"""Create publication-style summary figures for the SIDM--SMBH scan."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize


HALO_FACTORS = (0.1, 1.0, 10.0)
BH_MASSES = (4.0e5, 4.0e6, 4.0e7)
SIGMAS = (10.0, 100.0, 1000.0)
SIGMA_COLORS = {10.0: "#3b82f6", 100.0: "#d97706", 1000.0: "#16803c"}


def finite(value) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(value)


def save(fig: plt.Figure, outdir: Path, stem: str) -> None:
    fig.savefig(outdir / f"{stem}.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def interface_figure(rows: list[dict], outdir: Path) -> None:
    values = [
        row["r_in_over_r_h"]
        for row in rows
        if finite(row.get("r_in_over_r_h")) and row["r_in_over_r_h"] > 0
    ]
    norm = Normalize(
        vmin=math.log10(min(values)),
        vmax=math.log10(max(values)),
    )
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.8), sharex=True, sharey=True)
    image = None
    for ax, sigma in zip(axes, SIGMAS):
        matrix = np.full((3, 3), np.nan)
        labels = [["" for _ in range(3)] for _ in range(3)]
        for row in rows:
            if row["sigma0_over_m_cm2_g"] != sigma:
                continue
            ix = HALO_FACTORS.index(row["halo_mass_factor"])
            iy = BH_MASSES.index(row["black_hole_mass_msun"])
            ratio = row.get("r_in_over_r_h")
            if finite(ratio) and ratio > 0:
                matrix[iy, ix] = math.log10(ratio)
                admitted = row.get("fp_admissibility_status") == "FP_ADMISSIBLE"
                labels[iy][ix] = f"{ratio:.2g}" + ("\nFP" if admitted else "")
            else:
                labels[iy][ix] = "no crossing"
        image = ax.imshow(matrix, origin="lower", cmap="viridis", norm=norm)
        for iy in range(3):
            for ix in range(3):
                value = matrix[iy, ix]
                color = "white" if finite(value) and norm(value) < 0.55 else "black"
                ax.text(ix, iy, labels[iy][ix], ha="center", va="center", color=color, fontsize=8)
        ax.set_title(rf"$\sigma_0/m={sigma:g}\;{{\rm cm^2\,g^{{-1}}}}$")
        ax.set_xticks(range(3), ["0.1", "1", "10"])
        ax.set_yticks(range(3), [r"$4\times10^5$", r"$4\times10^6$", r"$4\times10^7$"])
        ax.set_xlabel(r"represented halo mass / fiducial")
    axes[0].set_ylabel(r"$M_\bullet\ [{\rm M_\odot}]$")
    colorbar = fig.colorbar(image, ax=axes, pad=0.02, fraction=0.035)
    colorbar.set_label(r"$\log_{10}(r_{\rm in}/R_h)$")
    fig.suptitle(
        "Collision-defined matching radius after the black-hole remap\n"
        "Cells marked FP pass the full phase-space domain audit",
        y=1.08,
    )
    save(fig, outdir, "fig_scan_interface_map")


def mass_flux_figure(rows: list[dict], outdir: Path) -> None:
    ready = [row for row in rows if finite(row.get("direct_steady_mdot_msun_per_myr"))]
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.2))
    x = np.arange(len(ready))
    if ready:
        series = (
            ("candidate_isotropic_mdot_msun_per_myr", "isotropic injection", "#a3a3a3", "o"),
            ("direct_no_rescatter_isotropic_mdot_msun_per_myr", "surviving direct", "#60a5fa", "s"),
            ("direct_steady_mdot_msun_per_myr", "steady depleted FP", "#1d4ed8", "D"),
            ("bondi_benchmark_mdot_msun_per_myr", "Bondi dimensional scale", "#dc2626", "^"),
        )
        for field, label, color, marker in series:
            axes[0].plot(
                x,
                [row[field] for row in ready],
                marker=marker,
                ms=4.2,
                lw=1.1,
                color=color,
                label=label,
            )
        axes[0].set_yscale("log")
        axes[0].set_xlabel("FP-ready scan case")
        axes[0].set_ylabel(r"mass current $[\mathrm{M_\odot\,Myr^{-1}}]$")
        axes[0].legend(frameon=False, fontsize=8)

        axes[1].plot(
            x,
            [abs(row["direct_binding_current_msun_kms2_per_myr"]) for row in ready],
            "o-",
            ms=4.2,
            lw=1.1,
            color="#7c3aed",
            label="capture binding-energy current",
        )
        axes[1].plot(
            x,
            [abs(row["direct_thermal_sink_msun_kms2_per_myr"]) for row in ready],
            "s-",
            ms=4.2,
            lw=1.1,
            color="#111827",
            label="local thermal sink",
        )
        axes[1].plot(
            x,
            [abs(row["maximum_absolute_gravothermal_luminosity_msun_kms2_per_myr"]) for row in ready],
            "^-",
            ms=4.2,
            lw=1.1,
            color="#d97706",
            label="maximum fluid conductive flux",
        )
        axes[1].set_yscale("log")
        axes[1].set_xlabel("FP-ready scan case")
        axes[1].set_ylabel(r"energy current $[\mathrm{M_\odot\,(km\,s^{-1})^2\,Myr^{-1}}]$")
        axes[1].legend(frameon=False, fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.2)
    fig.suptitle("Mass and energy-current hierarchy")
    save(fig, outdir, "fig_scan_flux_hierarchy")


def radial_figure(radial_csv: Path, outdir: Path) -> None:
    groups = defaultdict(list)
    fp_points = {}
    with radial_csv.open(newline="") as stream:
        for row in csv.DictReader(stream):
            key = (
                float(row["represented_halo_mass_msun"]),
                float(row["black_hole_mass_msun"]),
                float(row["sigma0_over_m_cm2_g"]),
            )
            if row["radial_source"] == "fluid_shell":
                groups[key].append(
                    (
                        float(row["radius_pc"]),
                        float(row["gravothermal_luminosity_msun_kms2_per_myr"]),
                        float(row["r_h_pc"]) if row["r_h_pc"] else np.nan,
                    )
                )
            elif row["fp_direct_thermal_sink_msun_kms2_per_myr"]:
                fp_points[key] = (
                    float(row["radius_pc"]),
                    float(row["fp_direct_thermal_sink_msun_kms2_per_myr"]),
                    float(row["r_h_pc"]),
                )

    fig, axes = plt.subplots(3, 3, figsize=(11.2, 9.0), sharex=True, sharey=True)
    halo_masses = sorted({key[0] for key in groups})
    for row_index, bh_mass in enumerate(BH_MASSES):
        for col_index, halo_mass in enumerate(halo_masses):
            ax = axes[row_index, col_index]
            for sigma in SIGMAS:
                values = sorted(groups.get((halo_mass, bh_mass, sigma), []))
                if not values:
                    continue
                radius = np.asarray([item[0] for item in values])
                luminosity = np.asarray([item[1] for item in values])
                rh = values[0][2]
                x = radius / rh if finite(rh) and rh > 0 else radius
                ax.plot(x, luminosity, lw=1.15, color=SIGMA_COLORS[sigma], label=f"{sigma:g}")
                point = fp_points.get((halo_mass, bh_mass, sigma))
                if point is not None:
                    ax.scatter(
                        point[0] / point[2],
                        point[1],
                        marker="v",
                        s=26,
                        color=SIGMA_COLORS[sigma],
                        edgecolor="black",
                        linewidth=0.35,
                        zorder=4,
                    )
            ax.axhline(0.0, color="0.5", lw=0.6)
            ax.set_xscale("log")
            ax.set_yscale("symlog", linthresh=1.0e3, linscale=0.5)
            ax.grid(alpha=0.15)
            if row_index == 0:
                factor = halo_mass / 1.5283226426126614e10
                ax.set_title(rf"$M_h/M_{{h,0}}={factor:g}$")
            if col_index == 0:
                ax.set_ylabel(
                    rf"$M_\bullet={bh_mass:.0e}\,\mathrm{{M_\odot}}$" + "\n" + r"$L\ [\mathrm{M_\odot(km\,s^{-1})^2\,Myr^{-1}}]$"
                )
            if row_index == 2:
                ax.set_xlabel(r"$r/R_h$")
    handles = [
        plt.Line2D([0], [0], color=SIGMA_COLORS[sigma], lw=1.5, label=rf"$\sigma_0/m={sigma:g}$")
        for sigma in SIGMAS
    ]
    handles.append(
        plt.Line2D([0], [0], marker="v", color="none", markerfacecolor="0.6", markeredgecolor="black", label="FP thermal sink at boundary")
    )
    fig.suptitle("Gravothermal luminosity after the black-hole remap", y=0.995)
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.965),
        ncol=4,
        frameon=False,
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    save(fig, outdir, "fig_scan_radial_flux_atlas")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aggregate-json", type=Path, required=True)
    parser.add_argument("--radial-csv", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    aggregate = json.loads(args.aggregate_json.read_text())
    args.outdir.mkdir(parents=True, exist_ok=True)
    interface_figure(aggregate["rows"], args.outdir)
    mass_flux_figure(aggregate["rows"], args.outdir)
    radial_figure(args.radial_csv, args.outdir)
    print("SCAN_FIGURES_COMPLETE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
