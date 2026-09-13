# Reproducing the fixed-snapshot calculations

The repository includes the native fluid engine, the SIDM-modified GNC source,
the Python non-local collision solver, production inputs, and scan and plotting
scripts. The top-level README is intentionally blank.

The checks completed for this release are recorded in `validation/reproduction_audit.json`.
They cover the tests, a fresh production-profile preparation, re-analysis of the
archived fluid trajectories, and regeneration of the reference figures. The full
38-task collision calculation, all 1,200 scan cases, and native MPI execution have
not been rerun as part of this release audit.

## Environment and inputs

Follow `SETUP.md`, or install the recorded Python 3.12 environment:

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-fluid-lock.txt
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export REPO="$PWD"
export PYTHON_BIN="$REPO/.venv/bin/python"
make audit test fluid-smoke
python tools/verify_reproduction_inputs.py
```

`data/production/input/prof_deep.txt` is the fluid-evolved, pre-black-hole snapshot
used by the production calculation. Its columns are radius in pc, shell-averaged
density in solar masses per cubic pc, and one-dimensional velocity dispersion in
km/s. The matched post-black-hole profile, hydrostatic bridge, normalized DF, and
preparation metadata are supplied alongside it.

`data/manifest.json` records the released files and their origin. Physical tables
are unchanged. Absolute paths in JSON were replaced with filenames, and references
to the changed metadata files were updated. These changes do not alter the recorded
numbers or acceptance gates.

The archived outputs in `data/reference` allow figures and response diagnostics to
be checked without rerunning the collision Monte Carlo calculation. The large
intermediate collision operators are generated outputs and are not bundled.
The exact NFW-to-collapse restart history is also not bundled. Reproduction of the
production calculation begins at the supplied deep-collapse snapshot. The example
halo generator in `fluid_solver` is not a substitute for that input.

## Prepare the production snapshot

From the repository root, choose a new output directory:

```sh
export ROOT="$REPO/results/production"
export PKG="$REPO/fp_solver"
export PROFILE="$REPO/data/production/input/prof_deep.txt"
bash slurm/fp_solver/run_production_ej_prep.sbatch
```

The script can also be submitted with `sbatch` after selecting a suitable partition.
It adds a four-million-solar-mass central black hole through an adiabatic remap,
constructs the mass-matched inner bridge, evaluates the three reporting surfaces,
and normalizes the DF. The scattering parameters are `sigma0/m = 100 cm^2/g` and
`w = 80 km/s`. Every preparation gate must pass before proceeding.

The local replay matched the archived post-black-hole profile and bridge to maximum
relative differences of about 9e-15 and 1e-12, respectively. The DF table matched to
4e-8. These are numerical comparisons, not byte-identical-output requirements.
See `validation/preparation_comparison.json`.

To start directly from the archived matched state instead, use a different new
output directory and copy the supplied checkpoints:

```sh
export ROOT="$REPO/results/from_checkpoint"
test ! -e "$ROOT" && mkdir -p "$ROOT"
cp -R data/production/. "$ROOT/"
```

Keep the profile, bridge, and DF metadata together. Mixing a freshly regenerated
profile with an archived closure will fail the input-identity checks.

## Generate the collision operators and solve the depleted distribution

`configs/production_operators.json` contains the 38 independent-seed tasks for the
19 resolution and reporting-surface groups. Each task uses 16 radial bins and
1.6 million sampled pairs per bin. This is the production workload, not a smoke test.

```sh
for INDEX in $(seq 0 37); do
  python tools/reproduce_production.py operator --root "$ROOT" --index "$INDEX"
done
python tools/reproduce_production.py steady --root "$ROOT"
python tools/reproduce_production.py aggregate --root "$ROOT"
```

The operator tasks are independent and can be distributed across a Slurm array.
The sequential loop above defines the full task list. The steady stage combines
the two seeds within each group and solves each seed and combined operator.
The driver refuses existing stage outputs rather than overwriting earlier runs.
It avoids the historical restart dependencies in the older Slurm ladder scripts.

The aggregate stage records separate mass and capture-binding-energy convergence
gates. A mass-only closure is written only if its mass, seed, resolution,
reporting-surface, and physical-identity checks pass. A failed energy-current gate
remains failed. The current energy-moment audit is retained in the steady solver.
The thermal term used below is the specified local kinetic-energy loss per captured
mass, not the capture binding-energy current.

## Run the black-hole-aware fluid response

With a newly accepted mass closure from the previous stage:

```sh
export FLUID="$ROOT/fluid_reproduction"
python fluid_solver/run_response.py \
  --profile "$ROOT/remap/profile_bh.txt" --profile-state post-bh \
  --diagnostics "$ROOT/aggregate_reproduction/fluid_mass_closure.json" \
  --bridge-json "$ROOT/bridge/bridge.json" --outdir "$FLUID" \
  --mbh 4e6 --sigma-over-m 100 --w-units 80 \
  --branches control,sink --tau-end 0.06 --output-dtau 0.002 \
  --t-epsilon 1e-4 --wall-seconds 25000 \
  --max-unmodeled-fluid-shell-mass-fraction 0.02 \
  --max-unmodeled-gnc-mass-fraction 0.20
python fp_solver/analyze_fluid_response.py \
  --summary "$FLUID/fluid_feedback_summary.json" \
  --trajectories "$FLUID/fluid_feedback_trajectories.csv" \
  --closure-json "$ROOT/aggregate_reproduction/fluid_mass_closure.json" \
  --profile "$ROOT/remap/profile_bh.txt" \
  --bridge-json "$ROOT/bridge/bridge.json" \
  --out-json "$FLUID/fluid_response_analysis.json" \
  --out-csv "$FLUID/fluid_response_comparison.csv"
```

A successful process exit is insufficient. Inspect the analyzer's status, all
acceptance gates, stop reasons, and common evolved interval. The closure is held
fixed at the matched snapshot. This is a short fluid response, with explicit limits
on the omitted captured mass, rather than a continuously updated two-way evolution.
The density diagnostic is the mean inside the innermost resolved Lagrangian shell.
It is not an extrapolated central cusp density.

## Re-analyze the archived response and regenerate its figures

This uses the measured reference trajectories and does not rerun their evolution:

```sh
export REF="$REPO/data/reference/production"
export MASS="$REF/aggregate_v6_mass_reaudit_boundaryfix"
export ARCHIVE="$REF/fluid_mass_v6_final_20260901_0955"
export CHECK="$REPO/results/reference_check"
mkdir -p "$CHECK"
python fp_solver/analyze_fluid_response.py \
  --summary "$ARCHIVE/fluid_feedback_summary.json" \
  --trajectories "$ARCHIVE/fluid_feedback_trajectories.csv" \
  --closure-json "$MASS/fluid_mass_closure.json" \
  --profile "$REPO/data/production/remap/profile_bh.txt" \
  --bridge-json "$REPO/data/production/bridge/bridge.json" \
  --out-json "$CHECK/fluid_response_analysis.json" \
  --out-csv "$CHECK/fluid_response_comparison.csv"
python fp_solver/plot_production_results.py \
  --convergence-json "$MASS/ej_convergence.json" \
  --closure-json "$MASS/fluid_mass_closure.json" \
  --fluid-analysis-json "$CHECK/fluid_response_analysis.json" \
  --fluid-comparison-csv "$CHECK/fluid_response_comparison.csv" \
  --outdir "$CHECK/production_figures"
python figures/plot_mass_convergence.py \
  --convergence-json "$MASS/ej_convergence.json" \
  --outdir "$CHECK/manuscript_figures"
```

This archived result uses a converged mass current and a prescribed thermal sink.
It does not validate an absolute two-current energy closure. The stricter
`publication_inputs.py` loader for later audited response products has a separate
input contract. Do not rename a historical status to make that loader accept it.

## Reproduce the 20 by 20 by 3 scan

The dense scan varies black-hole mass and represented halo mass at three scattering
normalizations. It remaps each profile and measures the interface and the fluid
conductive luminosity. It does not contain 1,200 measured FP capture rates or 1,200
coupled halo evolutions. `parameter_scan` retains the older coarse workflow.

On Slurm, adjust the partition directives in `slurm/parameter_scan20` for your
cluster. The scripts inherit the activated environment and `PYTHON_BIN`.

```sh
export SCAN_ROOT="$REPO/results/scan20"
bash slurm/parameter_scan20/launch_scan20.sh
```

The launcher submits 50 batches of 24 cases, followed by fluid measurements and
aggregation. Use a new output directory. The seed defaults to the released
pre-black-hole snapshot. The released classifier labels missing critical-orbit
diagnostics as unavailable rather than treating them as a measured failure.
Archived scan tables retain their original labels and measurements.

To regenerate the manuscript heat maps from the archived measurements:

```sh
python figures/plot_manuscript_scan.py \
  --results-csv data/reference/scan20/scan20_results.csv \
  --outdir results/scan20_figures
```

## Reproduce the cross-regime calibration

The orbit-resolved current, a conductive hydrostatic spike, and Bondi inflow are
separate boundary-value problems. The calibration therefore keeps three branches
rather than blending a single FP coefficient monotonically into a Bondi rate.

```sh
make calibrate
```

This command performs two checks. First, it reads the ten vector-extracted
simulation markers from Sabarish et al. (2025) and refits their harmonic form,

```text
dotM = C + 1 / (A s + B / s).
```

The joint refit gives `C = 1.782 Msun/yr`, `A = 0.05705`, and `B = 0.27554`,
placing the conductivity turnover at `s = 2.198 cm^2/g`. Omitting one marker at a
time moves the turnover between 2.149 and 2.256 cm^2/g. These leave-one-out ranges
measure point sensitivity rather than a statistical confidence interval. The
paper prints `C = 3.67`, but that value is incompatible with its plotted fit and
its stated collisionless rate of `1.7 Msun/yr`; the code records the discrepancy
and determines all three coefficients from the vector data. The benchmark has an
LMFP excess proportional to `s` and an SMFP excess proportional to `1/s`. It
describes a conduction-supported isolated spike with different boundary
conditions from either a Bondi reservoir or the velocity-dependent production
halo.

Second, the command evaluates a nominal adiabatic Bondi benchmark and audits the
saved Yukawa profile from the ISCO to the nominal Bondi radius. For a monatomic SIDM fluid,
`gamma = 5/3`, `lambda_B = 1/4`, and `c_infinity = sqrt(gamma) sigma_1d`. At the
fiducial boundary, substituting the local density and dispersion for the
asymptotic reservoir values gives `dotM_B = 4.10e4 Msun/Myr`, about 2.92e3 times
the measured FP current. The saved hydrostatic state does not itself supply that
Bondi reservoir. The Bondi branch is also excluded by the collisionality audit,
because the saved profile has
`N_orb = 1.95e-8` at the ISCO. Orbital memory therefore survives inside the nominal
Bondi region even though `N_orb` is close to unity near its outer edge.

The output is written to `validation/cross_regime_calibration.json`, and the figure
is written to `validation/cross_regime_figures/`. The vector-extracted benchmark
and its metadata are under `data/calibration/`. The current code uses the exact
Bondi normalization in new coarse-scan summaries. Archived scan tables retain the
older dimensional comparison and are not silently rewritten.

The source figure is not redistributed. To reproduce the marker extraction from
`images/imfp.pdf` in the paper's arXiv source package, install the optional
dependency and run

```sh
python -m pip install -r requirements-calibration-extraction.txt
python tools/extract_sabarish_imfp.py /path/to/images/imfp.pdf \
  --out /tmp/sabarish_imfp.csv \
  --reference data/calibration/sabarish_2025_imfp_digitized.csv
```

The extractor identifies the ten black vector markers, applies the stored affine
axis transform, and checks the result against the released CSV.

`fp_solver/moving_interface_ledger.py` supplies conservative mass and energy
bookkeeping for a moving boundary. It requires measured specific energies,
luminosity, and work terms. The ledger does not turn the captured-mass current into
a thermal current.

To regenerate the calibration, its integrated acceptance ledger, every package
test, and the repository audit in one step, run

```sh
make final-calibration
```

The resulting `validation/final_calibration_audit.json` separates accepted
fixed-snapshot results from the energy and time-dependent coupling gates that have
not passed. A completed numerical process is never promoted across those gates.

## Native GNC calculations

The fixed-snapshot non-local operator workflow above runs in Python. The native
SIDM-modified GNC implementation and the plunge-record instrumentation are also
included under `vendor/GNC` and `fp_solver/integration`. Follow `SETUP.md` to prepare
and build them with MPI Fortran and HDF5. Compilation, coefficient-table generation,
and a physical MPI run still require validation in that environment. Historical
native-GNC run directories and generated coefficient tables are not bundled.
