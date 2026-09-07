# Upstream software

## GravothermalSIDM

- Author: Kimberly Boddy.
- Upstream: https://github.com/kboddy/GravothermalSIDM
- Revision: `fe5ed30436df46c08befb0999ec86f57385eb1ac`.
- Source: `fluid_solver/GravothermalSIDM/SourcePy`.
- Licence: GPL-3.0-or-later. The complete licence and upstream terms remain in
  `fluid_solver/GravothermalSIDM/LICENSE` and `TERMS`.

The bundled `evolve.py` includes the black-hole-gravity changes used on the
cluster. Their difference from the upstream revision is preserved in
`fluid_solver/patches/black_hole_gravity.patch`. These modifications and the new
fluid-specific drivers are distributed under GPL-3.0-or-later. The top-level MIT
licence does not replace the fluid engine's licence or the terms governing its
derivatives.

## GNC

- Authors: Fupeng Zhang and Pau Amaro-Seoane.
- Upstream: https://github.com/zhangfupeng-gzhu/GNC
- Base revision: `44dffbb67fc89ab02a10e2afc9f5b4e3ed280029`.
- Source: `vendor/GNC`.
- Licence: MIT, retained in `vendor/GNC/LICENSE`.

This is the cluster's modified SIDM source tree, not an unmodified upstream
release. Its source line endings have been normalized, and machine-specific HDF5
paths have been replaced with build arguments. Obsolete patching scripts, backup
makefiles, executables, and simulation outputs are excluded. The additional
coupling patches remain in `fp_solver/integration` and are installed by
`tools/build_gnc.py` into a new build directory.

The public example input deck and empty-directory layout are retained. Their
presence does not make the example a production SIDM run. Cite the upstream
software and its associated papers when using these solvers.
