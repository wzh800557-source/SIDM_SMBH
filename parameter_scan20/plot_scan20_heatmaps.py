#!/usr/bin/env python3
"""Draw publication-style heat maps from the dense scan."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D


SIGMAS = (10.0, 100.0, 1000.0)


def maybe_float(value: str | None) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return math.nan
    return result if math.isfinite(result) else math.nan


def geometric_edges(values: np.ndarray) -> np.ndarray:
    log_values = np.log10(values)
    edges = np.empty(values.size + 1)
    edges[1:-1] = 0.5 * (log_values[:-1] + log_values[1:])
    edges[0] = log_values[0] - 0.5 * (log_values[1] - log_values[0])
    edges[-1] = log_values[-1] + 0.5 * (log_values[-1] - log_values[-2])
    return 10.0**edges


def read_rows(path: Path) -> list[dict]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def grid_axes(rows: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    halo = np.array(sorted({maybe_float(row["represented_halo_mass_msun"]) for row in rows}))
    black_hole = np.array(sorted({maybe_float(row["black_hole_mass_msun"]) for row in rows}))
    return halo[np.isfinite(halo)], black_hole[np.isfinite(black_hole)]


def matrix(
    rows: list[dict],
    sigma: float,
    field: str,
    halo_values: np.ndarray,
    black_hole_values: np.ndarray,
    absolute: bool = False,
) -> np.ndarray:
    out = np.full((black_hole_values.size, halo_values.size), np.nan)
    hx = {round(math.log10(value), 12): i for i, value in enumerate(halo_values)}
    by = {
        round(math.log10(value), 12): i
        for i, value in enumerate(black_hole_values)
    }
    for row in rows:
        if not math.isclose(maybe_float(row["sigma0_over_m_cm2_g"]), sigma):
            continue
        value = maybe_float(row.get(field))
        if absolute and math.isfinite(value):
            value = abs(value)
        x = hx.get(round(math.log10(maybe_float(row["represented_halo_mass_msun"])), 12))
        y = by.get(round(math.log10(maybe_float(row["black_hole_mass_msun"])), 12))
        if x is not None and y is not None:
            out[y, x] = value
    return out


def admissibility_matrix(
    rows: list[dict],
    sigma: float,
    halo_values: np.ndarray,
    black_hole_values: np.ndarray,
) -> np.ndarray:
    out = np.full((black_hole_values.size, halo_values.size), np.nan)
    hx = {round(math.log10(value), 12): i for i, value in enumerate(halo_values)}
    by = {
        round(math.log10(value), 12): i
        for i, value in enumerate(black_hole_values)
    }
    for row in rows:
        if not math.isclose(maybe_float(row["sigma0_over_m_cm2_g"]), sigma):
            continue
        x = hx.get(round(math.log10(maybe_float(row["represented_halo_mass_msun"])), 12))
        y = by.get(round(math.log10(maybe_float(row["black_hole_mass_msun"])), 12))
        if x is None or y is None:
            continue
        value = str(row.get("fp_geometric_admissible", "")).lower()
        out[y, x] = 1.0 if value == "true" else 0.0
    return out


def positive_limits(matrices: list[np.ndarray]) -> tuple[float, float]:
    values = np.concatenate(
        [item[np.isfinite(item) & (item > 0.0)] for item in matrices]
    )
    if values.size == 0:
        return 1.0, 10.0
    return float(np.nanmin(values)), float(np.nanmax(values))


def style_axis(ax: plt.Axes) -> None:
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.tick_params(direction="out", length=3.5, width=0.8)
    ax.set_xlabel(r"represented halo mass $M_h\ [\mathrm{M_\odot}]$")


def save(fig: plt.Figure, outdir: Path, stem: str) -> None:
    fig.savefig(outdir / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def interface_heatmap(rows: list[dict], outdir: Path) -> None:
    halo_values, black_hole_values = grid_axes(rows)
    x_edges = geometric_edges(halo_values)
    y_edges = geometric_edges(black_hole_values)
    data = [
        matrix(rows, sigma, "r_in_over_r_h", halo_values, black_hole_values)
        for sigma in SIGMAS
    ]
    admissible = [
        admissibility_matrix(rows, sigma, halo_values, black_hole_values)
        for sigma in SIGMAS
    ]
    vmin, vmax = positive_limits(data)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#e5e7eb")
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.1), sharex=True, sharey=True)
    image = None
    xx, yy = np.meshgrid(halo_values, black_hole_values)
    for ax, sigma, values, gate in zip(axes, SIGMAS, data, admissible):
        image = ax.pcolormesh(
            x_edges,
            y_edges,
            values,
            cmap=cmap,
            norm=LogNorm(vmin=vmin, vmax=vmax),
            shading="flat",
            rasterized=True,
        )
        valid_gate = np.where(np.isfinite(gate), gate, np.nan)
        if np.any(np.isfinite(valid_gate)):
            ax.contour(
                xx,
                yy,
                valid_gate,
                levels=[0.5],
                colors="white",
                linewidths=1.4,
            )
            invalid = np.ma.masked_where(valid_gate != 0.0, valid_gate)
            ax.contourf(
                xx,
                yy,
                invalid,
                levels=[-0.5, 0.5],
                colors="none",
                hatches=["////"],
            )
        style_axis(ax)
        ax.set_title(
            rf"$\sigma_0/m={sigma:g}\ \mathrm{{cm^2\,g^{{-1}}}}$",
            fontsize=11,
        )
    axes[0].set_ylabel(r"black-hole mass $M_\bullet\ [\mathrm{M_\odot}]$")
    colorbar = fig.colorbar(image, ax=axes, pad=0.02, fraction=0.035)
    colorbar.set_label(r"$r_{\rm in}/R_h$")
    legend = [
        Line2D([0], [0], color="white", lw=1.5, label="geometric FP-domain boundary"),
        Line2D([0], [0], color="0.35", lw=0, marker="s", markersize=9,
               markerfacecolor="white", markeredgecolor="0.35",
               label="hatched: outside conservative FP domain"),
    ]
    fig.legend(handles=legend, loc="lower center", ncol=2, frameon=False, fontsize=8)
    fig.suptitle(
        "Collision-defined inner interface after the black-hole remap",
        fontsize=14,
        y=1.02,
    )
    fig.subplots_adjust(bottom=0.20, wspace=0.08)
    save(fig, outdir, "fig_scan20_interface_heatmap")


def fluid_heatmap(rows: list[dict], outdir: Path) -> None:
    halo_values, black_hole_values = grid_axes(rows)
    x_edges = geometric_edges(halo_values)
    y_edges = geometric_edges(black_hole_values)
    fields = (
        (
            "first_shell_gravothermal_luminosity_msun_kms2_per_myr",
            r"$|L_{\rm fluid}(r_0)|$",
        ),
        (
            "maximum_absolute_gravothermal_luminosity_msun_kms2_per_myr",
            r"$\max_r |L_{\rm fluid}(r)|$",
        ),
    )
    matrices = {
        field: [
            matrix(rows, sigma, field, halo_values, black_hole_values, absolute=True)
            for sigma in SIGMAS
        ]
        for field, _ in fields
    }
    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("#e5e7eb")
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 7.4), sharex=True, sharey=True)
    for row_index, (field, label) in enumerate(fields):
        vmin, vmax = positive_limits(matrices[field])
        image = None
        for col_index, (sigma, values) in enumerate(zip(SIGMAS, matrices[field])):
            ax = axes[row_index, col_index]
            image = ax.pcolormesh(
                x_edges,
                y_edges,
                values,
                cmap=cmap,
                norm=LogNorm(vmin=vmin, vmax=vmax),
                shading="flat",
                rasterized=True,
            )
            style_axis(ax)
            if row_index == 0:
                ax.set_title(
                    rf"$\sigma_0/m={sigma:g}\ \mathrm{{cm^2\,g^{{-1}}}}$",
                    fontsize=11,
                )
            if col_index == 0:
                ax.set_ylabel(
                    r"black-hole mass $M_\bullet\ [\mathrm{M_\odot}]$"
                )
        colorbar = fig.colorbar(image, ax=axes[row_index, :], pad=0.02, fraction=0.02)
        colorbar.set_label(
            label + r" $[\mathrm{M_\odot\,(km\,s^{-1})^2\,Myr^{-1}}]$"
        )
        axes[row_index, 0].text(
            -0.27,
            0.5,
            "first fluid shell" if row_index == 0 else "maximum over radius",
            rotation=90,
            transform=axes[row_index, 0].transAxes,
            va="center",
            ha="center",
            fontsize=10,
            color="0.35",
        )
    fig.suptitle(
        "Black-hole-aware gravothermal luminosity on the dense scan",
        fontsize=14,
        y=0.995,
    )
    fig.subplots_adjust(left=0.10, right=0.89, bottom=0.10, top=0.93, hspace=0.16, wspace=0.08)
    save(fig, outdir, "fig_scan20_fluid_luminosity_heatmap")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-csv", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = read_rows(args.results_csv)
    interface_heatmap(rows, args.outdir)
    fluid_heatmap(rows, args.outdir)
    print("SCAN20_HEATMAPS_COMPLETE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
