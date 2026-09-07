#!/usr/bin/env python3
"""Prepare a private GNC build tree, install the coupling patches, and compile."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--build-dir", type=Path, default=ROOT / "build" / "GNC")
    p.add_argument("--prepare-only", action="store_true")
    p.add_argument("--fc", default=os.environ.get("MPIFC", "mpif90"))
    p.add_argument("--hdf5-libdir", default=os.environ.get("HDF5_LIBDIR"))
    p.add_argument("--hdf5-includedir", default=os.environ.get("HDF5_INCLUDEDIR"))
    args = p.parse_args()
    if not args.prepare_only:
        if not shutil.which(args.fc):
            p.error("MPI Fortran compiler not found. Set --fc or MPIFC.")
        if not args.hdf5_libdir or not args.hdf5_includedir:
            p.error("set HDF5_LIBDIR and HDF5_INCLUDEDIR for your Fortran-enabled HDF5 installation")
        if not (Path(args.hdf5_includedir) / "hdf5.mod").is_file():
            p.error("HDF5_INCLUDEDIR does not contain hdf5.mod")
    build = args.build_dir.resolve()
    if build.exists():
        p.error("build directory already exists. Choose a new directory to preserve existing work.")
    shutil.copytree(ROOT / "vendor" / "GNC", build)
    for relative in json.loads((build / "output_directories.json").read_text()):
        directory = Path(relative)
        if directory.is_absolute() or ".." in directory.parts:
            raise ValueError("invalid output skeleton directory")
        (build / directory).mkdir(parents=True, exist_ok=True)
    patch = ROOT / "fp_solver" / "integration" / "gnc_absolute_closure.patch"
    subprocess.run(["patch", "--batch", "--forward", "-p1", "-i", str(patch)], cwd=build, check=True)
    loader = ROOT / "fp_solver" / "integration" / "absolute_xj_loader.f90.inc"
    target = build / "source" / "ini_single.f90"
    with target.open("a") as stream:
        stream.write("\n" + loader.read_text())
    subprocess.run([sys.executable, str(ROOT / "fp_solver" / "install_sidm_born_kernel.py"), str(build)], check=True)
    state = {"status": "PREPARED", "compiled": False,
             "absolute_normalization_patch": True, "plunge_records_patch": True,
             "weighted_xj_loader": True, "born_kernel": True}
    if not args.prepare_only:
        opts = [f"FC={args.fc}", f"hdf5dirlib={Path(args.hdf5_libdir).resolve()}",
                f"hdf5dirinc={Path(args.hdf5_includedir).resolve()}"]
        # The upstream module dependency graph is serial. Do not use parallel make.
        subprocess.run(["make", "-j1", *opts], cwd=build / "source", check=True)
        subprocess.run(["make", "-j1", *opts], cwd=build / "main", check=True)
        state.update(status="COMPILED", compiled=True)
    (build / "build_status.json").write_text(json.dumps(state, indent=2) + "\n")
    print(json.dumps(state, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
