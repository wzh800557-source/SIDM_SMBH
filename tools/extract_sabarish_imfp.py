#!/usr/bin/env python3
"""Extract the ten black simulation markers from the vector IMFP figure.

The input is ``images/imfp.pdf`` from the arXiv source package for
Sabarish et al. (2025). The source PDF is not redistributed here.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

try:
    import pdfplumber
except ImportError as exc:  # pragma: no cover - exercised only without the optional tool
    raise SystemExit(
        "pdfplumber is required for vector extraction; install "
        "requirements-calibration-extraction.txt"
    ) from exc


# Affine coordinates of the plotted axes. The vertical coordinate follows
# pdfplumber's distance from the top of the page.
X_AT_ZERO = 63.123352
X_PER_UNIT = 40.5818181818
TOP_AT_RATE_TWO = 254.1278125
TOP_PER_RATE_UNIT = 40.32
EXPECTED_CROSS_SECTIONS = np.asarray(
    [0.02, 0.2, 0.4, 0.6, 0.8, 1.4, 2.0, 4.0, 6.0, 8.0]
)


def extract_points(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with pdfplumber.open(path) as document:
        if len(document.pages) != 1:
            raise ValueError("expected a one-page vector figure")
        page = document.pages[0]
        markers = []
        for curve in page.curves:
            black = curve.get("stroking_color") == 0
            five_point_marker = (
                abs(float(curve.get("width", 0.0)) - 5.0) < 1.0e-6
                and abs(float(curve.get("height", 0.0)) - 5.0) < 1.0e-6
            )
            inside_axes = (
                46.8 < float(curve.get("x0", -1.0)) < 404.1
                and 12.0 < float(curve.get("top", -1.0)) < 278.5
            )
            if black and five_point_marker and inside_axes:
                x = 0.5 * (float(curve["x0"]) + float(curve["x1"]))
                top = 0.5 * (float(curve["top"]) + float(curve["bottom"]))
                markers.append((x, top))
    if len(markers) != EXPECTED_CROSS_SECTIONS.size:
        raise ValueError(f"expected 10 simulation markers, found {len(markers)}")

    markers.sort()
    x = np.asarray([point[0] for point in markers])
    top = np.asarray([point[1] for point in markers])
    extracted_cross_section = (x - X_AT_ZERO) / X_PER_UNIT
    rate = 2.0 + (TOP_AT_RATE_TWO - top) / TOP_PER_RATE_UNIT
    if not np.allclose(
        extracted_cross_section, EXPECTED_CROSS_SECTIONS, rtol=0.0, atol=2.0e-8
    ):
        raise ValueError("marker abscissae do not match the expected source figure")
    return EXPECTED_CROSS_SECTIONS.copy(), rate


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("figure", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--reference", type=Path)
    args = parser.parse_args()

    cross_section, rate = extract_points(args.figure)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sigma_over_m_cm2_g", "accretion_rate_msun_per_yr"])
        writer.writerows(zip(cross_section, rate))

    if args.reference is not None:
        reference = np.genfromtxt(args.reference, delimiter=",", names=True)
        if not np.allclose(
            reference["sigma_over_m_cm2_g"],
            cross_section,
            rtol=0.0,
            atol=2.0e-8,
        ) or not np.allclose(
            reference["accretion_rate_msun_per_yr"],
            rate,
            rtol=0.0,
            atol=2.0e-8,
        ):
            raise RuntimeError("extracted vector markers differ from the reference CSV")
    print(f"EXTRACTED_SABARISH_IMFP_MARKERS count={rate.size}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
