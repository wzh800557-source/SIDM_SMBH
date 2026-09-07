#!/usr/bin/env python3
"""Install the validated t-channel Born kernel into a copied GNC source tree."""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path


def replace_once(text: str, pattern: str, replacement: str, label: str) -> str:
    out, n = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if n != 1:
        raise RuntimeError(f"expected one {label} insertion point, found {n}")
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("build_root", type=Path)
    p.add_argument(
        "--module-template",
        type=Path,
        default=Path(__file__).with_name("sidm_params_born_fixed.f90"),
    )
    args = p.parse_args()
    source = args.build_root / "source"
    md = source / "md_main_gw.f90"
    sidm = source / "sidm_params.f90"
    if not md.is_file() or not sidm.is_file() or not args.module_template.is_file():
        raise FileNotFoundError("GNC source tree or SIDM module template is missing")

    text = md.read_text()
    if "initialize_sidm_kernel" not in text:
        text = replace_once(
            text,
            r"(subroutine\s+init_model_ctl\(\)\s*\n\s*use\s+com_main_gw\s*\n)",
            r"\1\tuse sidm_params, only: initialize_sidm_kernel\n",
            "SIDM module import",
        )
        text = replace_once(
            text,
            r"(\s*call\s+get_rh_vh_nh\(mbh,\s*rh,\s*ctl%n0,\s*ctl%v0\)\s*\n)",
            r"\1\tcall initialize_sidm_kernel(ctl%v0)\n",
            "SIDM initialization call",
        )
    md.write_text(text)
    shutil.copy2(args.module_template, sidm)

    check = md.read_text()
    if check.count("use sidm_params, only: initialize_sidm_kernel") != 1:
        raise RuntimeError("SIDM import installation is not unique")
    if check.count("call initialize_sidm_kernel(ctl%v0)") != 1:
        raise RuntimeError("SIDM initialization installation is not unique")
    print("SIDM_BORN_KERNEL_INSTALLED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
