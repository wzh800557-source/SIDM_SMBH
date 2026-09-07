"""Maxwellian relative-speed quadrature for diagnostic transport rates.

The background is an isotropic Gaussian with one-dimensional dispersion sigma.
This is not the general velocity distribution of a depleted nuclear DF. Existing
scan results used an RMS approximation and are not retroactively recalibrated.
"""
import math
import numpy as np
from numpy.polynomial.legendre import leggauss


def maxwellian_rate_average(test_speed, sigma, cross_section, nquad=192):
    """Return <sigma_cross(s) s> for a fixed test speed and Maxwellian background.

The caller supplies a vectorized cross section with any consistent area units.
Speed units are inherited. The tail is truncated at ten background dispersions.
"""
    v, sig = float(test_speed), float(sigma)
    if not math.isfinite(v + sig) or v < 0 or sig < 0 or nquad < 16:
        raise ValueError("Nonnegative finite speeds and nquad >= 16 are required")
    if sig == 0:
        return float(np.asarray(cross_section(np.array([v])))[0]) * v
    nodes, weights = leggauss(nquad)
    lower = max(-v / sig, -10.0)
    upper = 10.0
    z = lower + (nodes + 1) * (upper - lower) / 2
    speed = v + sig * z
    weights = weights * (upper - lower) / 2
    if v / sig < 1e-7:
        # Zero-speed Maxwell distribution; correction is second order in v/sig.
        x = speed / sig
        probability_dz = math.sqrt(2 / math.pi) * x*x * np.exp(-x*x/2)
    else:
        probability_dz = speed / (math.sqrt(2*math.pi) * v) * np.exp(-z*z/2) * (-np.expm1(-2*speed*v/(sig*sig)))
    area = np.asarray(cross_section(speed), float)
    if np.any(~np.isfinite(area)) or np.any(area < 0):
        raise ValueError("Cross section must be finite and nonnegative")
    return float(np.sum(weights * probability_dz * area * speed))
