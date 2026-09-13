# Final calibration of the available current branches

The released calculation now closes the captured-mass problem for the fixed
production snapshot and selects the correct transport branch without forcing an
FP coefficient to approach Bondi accretion.

The remap and hydrostatic bridge pass their mass and equilibrium checks. The
depleted finite-angle phase-space solution gives a direct captured-mass current of
14.0608 solar masses per Myr and an immediate-absorption current of 15.2254 solar
masses per Myr. The corresponding direct coefficient is
`C_M = 1.74816e-4`. The largest selected numerical or capture-rule change is
7.95 per cent.

The branch calibration uses the ten vector markers in the constant-cross-section
SIDM-spike calculation of Sabarish et al. (2025). A joint refit of
`dotM = C + 1/(A*s + B/s)` gives `C = 1.7819 solar masses per yr`,
`A = 0.057054`, and `B = 0.275543`. The conductive excess peaks at
`s = 2.1976 cm^2/g` and has an rms residual of `0.1178 solar masses per yr`.
Leave-one-out refits place the turnover between 2.1485 and 2.2556 cm2/g. The
constant-cross-section spike is an external conductivity benchmark, not a
calibration of the velocity-dependent production halo.

For the fiducial Yukawa profile, the orbital transport depth is 1.20 at the
nominal Bondi radius but only `1.95e-8` at the ISCO. The inflow path therefore
retains orbital memory. Its selected solver is the orbit-resolved FP branch. A
nominal Bondi calculation obtained by substituting the local interface state for
an asymptotic reservoir gives `4.10e4 solar masses per Myr`, 2,919 times the
measured FP current, but the saved hydrostatic profile does not satisfy the Bondi
collisionality or reservoir conditions.

The solver rule is now explicit. Use orbit-resolved FP where `N_orb < 1`. Use a
non-orbit-averaged radial kinetic solver through `N_orb` of order unity. Treat a
profile as a hydrodynamic candidate only when the complete supply-to-capture path
is continuum-like. The released operational sentinel is `N_orb >= 20*pi`, which
corresponds to `Kn <= 0.1` only when the pressure scale height equals the radius.
After this gate passes, the flow boundary condition decides between a hydrostatic
conductive spike and a maintained transonic Bondi inflow.

The moving-interface mass and energy ledger is implemented and unit tested. The
archived black-hole-aware fluid calculation also passes its numerical and
provenance checks for a prescribed local kinetic sink, reaching 0.0593 initial
relaxation times with a final inner-density ratio of 1.00108 relative to its
control. This response is not used as a physical heating or cooling claim.

The production capture-energy moment has not passed its energy-resolution gate.
An absolute two-way evolution therefore remains unauthorized. It requires a
converged capture-weighted energy exchange, a radial kinetic treatment of any
transition layer, and production use of the moving-interface work and energy
ledger.

Run `make final-calibration` to regenerate the calibration JSON, figure,
integrated acceptance ledger, tests, and repository audit.
