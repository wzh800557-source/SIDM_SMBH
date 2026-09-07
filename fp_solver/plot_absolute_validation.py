#!/usr/bin/env python3
"""Plot the hydrostatic bridge and physically normalized GNC round trip."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np

import hydrostatic_bridge as hb
import normalize_gnc_df as ng


def main(argv: Iterable[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--resolved-profile", type=Path, required=True)
    p.add_argument("--bridged-profile", type=Path, required=True)
    p.add_argument("--bridge-json", type=Path, required=True)
    p.add_argument("--df-json", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args(argv)

    dbridge = json.loads(args.bridge_json.read_text())
    ddf = json.loads(args.df_json.read_text())
    rr, rhor, sigr = hb.load_profile(args.resolved_profile)
    rp, rhop, sigp = hb.load_profile(args.bridged_profile)
    beta = float(dbridge["cusp_slope_beta"])
    alt = hb.build_one(
        rr,
        rhor,
        sigr,
        float(dbridge["mbh_msun"]),
        float(dbridge["sigma_over_m_cm2_g"]),
        float(dbridge["velocity_exponent_a"]),
        float(dbridge["rh_pc"]),
        float(dbridge["core_slope_gamma_used"]),
        beta,
        str(dbridge["alternate_mode"]),
        2400,
    )
    core_only = hb.build_one(
        rr,
        rhor,
        sigr,
        float(dbridge["mbh_msun"]),
        float(dbridge["sigma_over_m_cm2_g"]),
        float(dbridge["velocity_exponent_a"]),
        float(dbridge["rh_pc"]),
        float(dbridge["core_slope_gamma_used"]),
        beta,
        "core-only",
        2400,
    )

    fig, ax = plt.subplots(1, 3, figsize=(12.4, 3.85))
    blue, orange, gray = "#3569a8", "#d9792b", "#555b61"

    inner = rp <= rr[0]
    ax[0].loglog(rr, rhor, color=gray, lw=2.1, label="resolved fluid")
    ax[0].loglog(rp[inner], rhop[inner], color=blue, lw=2.3, label=r"bridge at $R_h$")
    alt_inner = np.asarray(alt["r"]) <= rr[0]
    ax[0].loglog(
        np.asarray(alt["r"])[alt_inner],
        np.asarray(alt["rho"])[alt_inner],
        color=orange,
        lw=1.7,
        ls="--",
        label=r"cusp from $r_0$",
    )
    core_inner = np.asarray(core_only["r"]) <= rr[0]
    ax[0].loglog(
        np.asarray(core_only["r"])[core_inner],
        np.asarray(core_only["rho"])[core_inner],
        color="#8a8a8a",
        lw=1.2,
        ls=":",
        label="unrelaxed core",
    )
    ax[0].set_xlabel(r"radius [pc]")
    ax[0].set_ylabel(r"density [$M_\odot\,{\rm pc}^{-3}$]")
    ax[0].legend(frameon=False, fontsize=8.5, loc="lower left")

    som = float(dbridge["sigma_over_m_cm2_g"]) * hb.CM2_G_TO_PC2_MSUN
    n_alt = 2.0 * np.pi * np.asarray(alt["r"]) * np.asarray(alt["rho"]) * som
    n_core = 2.0 * np.pi * np.asarray(core_only["r"]) * np.asarray(core_only["rho"]) * som
    n_primary = 2.0 * np.pi * rp * rhop * som
    ax[1].loglog(rp, n_primary, color=blue, lw=2.3)
    ax[1].loglog(np.asarray(alt["r"]), n_alt, color=orange, lw=1.7, ls="--")
    ax[1].loglog(
        np.asarray(core_only["r"]), n_core, color="#8a8a8a", lw=1.2, ls=":"
    )
    ax[1].axhline(1.0, color="black", lw=1.0)
    envelope = dbridge["bridge_rin_envelope_pc"]
    ax[1].axvspan(envelope[0], envelope[1], color="#78b879", alpha=0.24, lw=0)
    ax[1].text(
        0.05,
        0.08,
        r"GNC: $N<1$" + "\n" + r"fluid: $N>1$",
        transform=ax[1].transAxes,
        fontsize=9,
    )
    ax[1].set_xlabel(r"radius [pc]")
    ax[1].set_ylabel(r"scatterings per orbit, $N$")

    rb = float(ddf["r_boundary_pc"])
    reval = np.geomspace(float(ddf["roundtrip_radius_range_pc"][0]), rb, 300)
    rho_target = np.array([ng.log_interp(rp, rhop, x) for x in reval])
    rho_back = ng.reconstructed_density(
        reval,
        float(ddf["rh_pc"]),
        float(ddf["n0_pc3"]),
        float(ddf["asymp_g_at_x_boundary"]),
        float(ddf["xmin"]),
        float(ddf["xmax"]),
        float(ddf["x_boundary"]),
        float(ddf["power"]),
        float(ddf["particle_mass_msun"]),
    )
    ax[2].semilogx(reval / rb, rho_back / rho_target, color=blue, lw=2.3)
    ax[2].axhline(1.0, color="black", lw=1.0)
    ax[2].scatter([1.0], [1.0], s=32, color=orange, zorder=3)
    ax[2].text(
        0.06,
        0.08,
        "boundary error: 0\n"
        + f"mass error: {100*float(ddf['roundtrip_mass_relative_error']):.2f}%",
        transform=ax[2].transAxes,
        fontsize=9,
    )
    ax[2].set_xlabel(r"$r/r_{\rm in}$")
    ax[2].set_ylabel(r"$\rho_{\rm GNC}/\rho_{\rm bridge}$")
    ax[2].set_ylim(0.94, 1.055)

    rh = float(dbridge["rh_pc"])
    r0 = float(dbridge["r0_pc"])
    for a in ax[:2]:
        a.axvline(rh, color="#6f6f6f", lw=1.0, ls=":")
        a.axvline(r0, color="#6f6f6f", lw=1.0, ls="-.")
        a.grid(which="major", color="#d8d8d8", lw=0.55, alpha=0.7)
    ax[2].grid(which="major", color="#d8d8d8", lw=0.55, alpha=0.7)
    ax[0].text(
        rh,
        0.95,
        r"$R_h$",
        transform=ax[0].get_xaxis_transform(),
        ha="right",
        va="top",
        fontsize=8,
    )
    ax[0].text(
        r0,
        0.84,
        r"$r_0$",
        transform=ax[0].get_xaxis_transform(),
        ha="left",
        va="top",
        fontsize=8,
    )
    for label, a in zip(("a", "b", "c"), ax):
        a.text(0.02, 0.96, f"({label})", transform=a.transAxes, va="top", fontweight="bold")

    fig.suptitle(
        "Hydrostatic bridge and absolute GNC normalization (local validation)",
        fontsize=12.2,
        y=1.01,
    )
    fig.tight_layout(w_pad=2.0)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
