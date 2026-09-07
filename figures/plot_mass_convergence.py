from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt



import argparse
p=argparse.ArgumentParser(description="Plot the archived mass-current convergence test.")
p.add_argument('--convergence-json',type=Path,required=True)
p.add_argument('--outdir',type=Path,required=True)
args=p.parse_args()
DATA=args.convergence_json
OUT=args.outdir
OUT.mkdir(parents=True,exist_ok=True)


def make_mass_convergence() -> Path:
    data = json.loads(DATA.read_text())
    mass_gate = data.get("observable_gates", {}).get("mass", {})
    if mass_gate.get("pass") is not True:
        raise RuntimeError("Mass convergence did not pass")
    radius = float(data["production_resolution_path"]["primary_boundary_radius_pc"])
    keys = [(48, 129), (64, 129), (64, 257)]
    groups = []
    for ne, nj in keys:
        match = [
            g for g in data["groups"]
            if int(g["energy_bins"]) == ne
            and int(g["angular_bins"]) == nj
            and abs(float(g["reservoir_radius_pc"]) - radius) < 1e-9
        ]
        if len(match) != 1:
            raise RuntimeError(f"missing mass-current group {ne}/{nj}")
        groups.append(match[0])

    final = float(groups[-1]["models"]["direct"]["mass"]["combined_value"])
    labels = ["48/129", "64/129", "64/257"]
    colours = {"direct": "#2b6ca3", "immediate": "#d0742c"}
    markers = {"direct": "o", "immediate": "s"}

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.linewidth": 0.7,
        "xtick.direction": "out",
        "ytick.direction": "out",
    })
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9)
    )
    fig.subplots_adjust(left=.09,right=.985,top=.9,bottom=.23,wspace=.32)
    #)

    ax = axes[0]
    x = list(range(3))
    ax.axhspan(0.9, 1.1, color="#dcebd8", alpha=0.72, zorder=0)
    ax.axhline(1.0, color="#666666", lw=1.0, ls="--")
    for model in ("direct", "immediate"):
        centres, lo, hi = [], [], []
        for group in groups:
            record = group["models"][model]["mass"]
            centre = float(record["combined_value"]) / final
            seeds = [float(v) / final for v in record["seed_values"]]
            centres.append(centre)
            lo.append(centre - min(seeds))
            hi.append(max(seeds) - centre)
        ax.errorbar(
            x, centres, yerr=[lo, hi], color=colours[model], marker=markers[model],
            ms=4, lw=1.2, capsize=2, label=model,
        )
    ax.set_xticks(x, labels)
    ax.set_ylim(0.88, 1.13)
    ax.set_ylabel("current / final direct current")
    ax.set_xlabel("energy bins / angular bins")
    ax.set_title("(a) Captured-mass convergence", loc="left", fontweight="bold")
    ax.text(0.03, 0.93, "shading: 10% tolerance", transform=ax.transAxes,
            color="#2f6f2f")
    ax.grid(axis="y", color="#d9d9d9", lw=0.7)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    row = [
        r for r in data["selected_comparisons"]["boundary"]
        if r["metric"] == "mass" and r["model"] == "direct"
    ][0]
    radii = [float(v) / radius for v in row["radii_pc"]]
    values = [float(v) for v in row["values"]]
    values = [v / values[0] for v in values]
    ax.axhspan(0.9, 1.1, color="#dcebd8", alpha=0.72, zorder=0)
    ax.axhline(1.0, color="#666666", lw=1.0, ls="--")
    ax.plot(radii, values, color=colours["direct"], marker="o", ms=4, lw=1.2)
    ax.set_xticks(radii, [f"{r:.3f}" for r in radii])
    ax.set_ylim(0.90, 1.10)
    ax.set_ylabel("current / current at inner crossing")
    ax.set_xlabel("reporting radius / inner crossing")
    ax.set_title("(b) Reporting surface", loc="left", fontweight="bold")
    ax.text(0.03, 0.93, f"mass span {100*(max(values)-min(values))/(sum(values)/len(values)):.2f}%", transform=ax.transAxes,
            color="#2f6f2f")
    ax.grid(axis="y", color="#d9d9d9", lw=0.7)
    ax.spines[["top", "right"]].set_visible(False)

    handles, names = axes[0].get_legend_handles_labels()
    fig.legend(handles, names, loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, -0.01))
    path = OUT / "fig_production_mass_convergence.png"
    fig.savefig(path.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    fig.savefig(path, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return path


if __name__ == "__main__":
    print(make_mass_convergence())
