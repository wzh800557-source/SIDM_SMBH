# Installation and code layout

The native fluid engine is in `fluid_solver/GravothermalSIDM/SourcePy`.
`evolve.py` advances the halo and `record.py` handles its records. The copy
includes the cluster's central-point-mass gravity modification. The fluid response
adds the black-hole-dependent conductivity and the conservative inner-face energy
term through `fp_solver/measured_fluid_feedback.py`. Its location is retained for
compatibility with existing run scripts. Use `fluid_solver/run_response.py` as the
fluid entry point.

`vendor/GNC` contains the cluster's SIDM-modified GNC source, without compiled
binaries or run outputs. `tools/build_gnc.py` copies it to a separate build
directory and installs our absolute-normalization, weighted-particle, plunge-record,
and Born-scattering patches before compilation. The finite-angle collision and
steady-state calculations are in `fp_solver`. The scan drivers are in
`parameter_scan` (coarse grid) and `parameter_scan20` (20 by 20 by 3 grid).
Cluster submission scripts are in `slurm`. See `REPRODUCING.md` for production
inputs, run commands, and the limits of the completed checks.

## Python environment

Use Python 3.12 for the tested fluid environment. Run these commands from the
repository root:

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-fluid.txt
export PYTHONDONTWRITEBYTECODE=1
make audit
make test
make fluid-smoke
```

`requirements-fluid-lock.txt` pins the versions used for the local Python 3.12
checks. Install from that file instead to reproduce the tested dependency versions.

The fluid engine is found relative to the repository, so an external checkout is
not required. Set `GRAVOTHERMAL_ROOT` only to use a different copy. The smoke test
advances a pure halo and two black-hole branches, checks that shell masses stay
fixed, and checks the conduction energy budget. Its small imposed boundary current
is synthetic. It does not certify a measured FP closure or gravothermal-collapse
result.

## Generate an example halo

```sh
python fluid_solver/generate_halo.py --outdir results/example_halo --steps 20
```

This evolves a newly generated NFW halo with constant cross section and writes
radius in pc, shell-averaged density in solar masses per cubic pc, and
one-dimensional velocity dispersion in km/s. The output is an installation example,
not the production deep-collapse snapshot. Existing output directories are never
overwritten. The approved production snapshot and matched checkpoints are in
`data/production`, with reference measurements in `data/reference`.

The matched fluid-response driver takes a remapped profile, bridge metadata, and
an accepted closure file. Its full argument list is available with:

```sh
python fluid_solver/run_response.py --help
```

## Build GNC

Install an MPI Fortran compiler, GNU make, `patch`, and HDF5 with its Fortran
interface. HDF5 and MPI must be compatible with the chosen compiler.

```sh
export MPIFC=mpif90
export HDF5_LIBDIR=/path/to/hdf5/lib
export HDF5_INCLUDEDIR=/path/to/hdf5/include
python tools/build_gnc.py --build-dir build/GNC
```

The compiled programs are `build/GNC/main/ini`, `main`, `pro`, and `cfuns`.
The build helper also recreates the example output-directory skeleton. Configure
the run deck and generate the coefficient tables before launching a physical run.
The public upstream tutorial is linked in `vendor/GNC/README.md`.

To check source preparation on a machine without MPI Fortran, use:

```sh
python tools/build_gnc.py --prepare-only --build-dir build/GNC_source_check
```

This checks patch installation and writes `build_status.json` with
`compiled: false`. It does not test compilation or MPI execution. The helper
refuses an existing build directory to protect previous work.

## Source provenance and licences

`source_manifest.json` records the upstream revisions and bundled source files.
The fluid engine retains its GPL-3.0-or-later licence, copyright notice, and terms.
GNC retains its MIT licence and attribution. See `THIRD_PARTY.md` for their
locations. The top-level README is intentionally blank.
