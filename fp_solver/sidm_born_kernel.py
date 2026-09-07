#!/usr/bin/env python3
"""Shared t-channel Born/Yukawa angular and conductivity moments.

The differential cross section is

    d sigma / d Omega = sigma_0/(4 pi)
        [1 + (v_rel/w)^2 sin^2(theta/2)]^-2.

``sigma_0`` is the zero-velocity total cross section in the differential kernel.
For isotropic scattering, sigma_V(0)=2 sigma_0/3.  The gravothermal K_p functions
absorb that constant into their calibrated coefficients and tend to unity as
v/w -> 0.
"""

from __future__ import annotations

import numpy as np


def tchannel_viscosity_ratio(v_over_w):
    """Return sigma_V(v_rel)/sigma_V(0) for the t-channel Born kernel.

    The exact expression is

        6 [(a+2) log(1+a) - 2a] / a^3,  a=(v_rel/w)^2.

    A series is used at small ``a`` to avoid catastrophic cancellation.
    """

    x = np.asarray(v_over_w, dtype=float)
    if np.any(x < 0.0) or np.any(~np.isfinite(x)):
        raise ValueError("v_over_w must be finite and non-negative")
    a = x * x
    out = np.empty_like(a)
    small = a < 1.0e-3
    # 1 - a + 9 a^2/10 - 4 a^3/5 + 5 a^4/7 + O(a^5)
    aa = a[small]
    out[small] = 1.0 - aa + 0.9 * aa**2 - 0.8 * aa**3 + (5.0 / 7.0) * aa**4
    bb = a[~small]
    out[~small] = 6.0 * ((bb + 2.0) * np.log1p(bb) - 2.0 * bb) / bb**3
    return float(out) if out.ndim == 0 else out


def tchannel_viscosity_over_sigma0(v_over_w):
    """Return sigma_V(v_rel)/sigma_0 for the Eq. 14 t-channel kernel."""

    return (2.0 / 3.0) * tchannel_viscosity_ratio(v_over_w)


_TCHANNEL_ORDER1 = {
    3: (8.0, 0.303941, 0.74, 0.68),
    5: (24.0, 0.826198, 0.82, 0.76),
    7: (48.0, 1.36217, 0.83, 0.79),
    9: (80.0, 1.90106, 0.84, 0.81),
}


def fluid_tchannel_kp(vdisp_over_w, p: int):
    """Return the order-one Outmezguine et al. t-channel K_p fit.

    ``vdisp_over_w`` uses the one-dimensional velocity dispersion, matching
    ``GravothermalSIDM.SourcePy.evolve``.  K_3 enters the LMFP conductivity and
    K_5 enters the SMFP conductivity.
    """

    if p not in _TCHANNEL_ORDER1:
        raise ValueError(f"unsupported K_p index {p}")
    x = np.asarray(vdisp_over_w, dtype=float)
    if np.any(x < 0.0) or np.any(~np.isfinite(x)):
        raise ValueError("vdisp_over_w must be finite and non-negative")
    par0, par1, par2, par3 = _TCHANNEL_ORDER1[p]
    s2 = x * x + 1.0e-4
    s4 = x**4 + 1.0e-8
    g = 1.5 * np.log1p((s2 * par1) ** par2) / (par0 * par2 * s4)
    out = (1.0 + 1.0 / g**par3) ** (-1.0 / par3)
    return float(out) if out.ndim == 0 else out


def gnc_effective_log_corrected(x_energy, v0_code: float, w_code: float, n: float = 4.0):
    """Exact t-channel viscosity moment mapped onto GNC's Landau prefactor.

    GNC stores x=|E|/v0^2, whereas ``w_code`` is in the absolute GNC velocity
    unit.  The physical ratio is therefore sqrt(2x)*v0_code/w_code.  With
    GNC's prefactor 2*pi*sigma0*w^4, the required multiplier is

        [(a+2) log(1+a) - 2a] / a,  a=(v/w)^2.

    This reproduces the exact viscosity cross section, including its y^4/6
    low-velocity limit and 2 log(y)-2 high-velocity limit.
    """

    x = np.asarray(x_energy, dtype=float)
    y = np.sqrt(2.0 * np.abs(x)) * float(v0_code) / float(w_code)
    a = y * y
    out = np.empty_like(a)
    small = a < 1.0e-3
    aa = a[small]
    out[small] = aa**2 / 6.0 - aa**3 / 6.0 + 3.0 * aa**4 / 20.0 - 2.0 * aa**5 / 15.0
    bb = a[~small]
    out[~small] = ((bb + 2.0) * np.log1p(bb) - 2.0 * bb) / bb
    return float(out) if out.ndim == 0 else out


def gnc_effective_log_legacy(x_energy, w_code: float, n: float = 4.0):
    """Reproduce the old, inverted GNC multiplier before its 0.1 floor."""

    x = np.asarray(x_energy, dtype=float)
    y = np.sqrt(2.0 * np.maximum(np.abs(x), 1.0e-30)) / float(w_code)
    out = 0.5 * np.log1p((1.0 / y) ** n)
    return float(out) if out.ndim == 0 else out
