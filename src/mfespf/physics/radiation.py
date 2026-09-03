"""Radiative primitives shared by the Stage 1 modules.

Implements EQ-PLANCK-SPECTRAL-RADIANCE and EQ-GREYBODY-EMISSION.  These are
functions rather than a module: they are constitutive physics several modules
evaluate, and giving them one implementation means the spectral and broadband
paths cannot silently disagree about what a black body does.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from mfespf.contracts.spectral import SpectralGrid
from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError

__all__ = [
    "band_integrated_exitance",
    "equilibrium_temperature",
    "greybody_exitance",
    "normalized_band_weights",
    "planck_spectral_exitance",
    "relaxation_timescale",
]


@lru_cache(maxsize=8)
def _quadrature(samples: int) -> tuple[np.ndarray, np.ndarray]:
    """Gauss-Legendre nodes and weights, computed once per sample count.

    The nodes are the roots of a Legendre polynomial and cost an eigenvalue
    solve; a column solve asks for them several times per step and they never
    change. Caching them alters no value (CMP-009).
    """
    nodes, weights = np.polynomial.legendre.leggauss(samples)
    nodes.setflags(write=False)
    weights.setflags(write=False)
    return nodes, weights


_H = value_of("planck")
_C = value_of("speed_of_light")
_KB = value_of("boltzmann")
_SIGMA = value_of("stefan_boltzmann")


def planck_spectral_exitance(wavelength_m: np.ndarray, temperature_k: float) -> np.ndarray:
    """Hemispheric spectral exitance ``M_lambda`` in W m^-2 m^-1.

    This is ``pi * B_lambda``: radiance integrated over the outward hemisphere for
    a Lambertian emitter, which is what a flux-based energy budget needs.  The
    exponent is evaluated with ``expm1`` so the Rayleigh-Jeans limit stays
    accurate instead of cancelling catastrophically.
    """
    lam = np.asarray(wavelength_m, dtype=float)
    if temperature_k <= 0.0:
        return np.zeros_like(lam)
    if np.any(lam <= 0.0):
        raise ConfigurationError("wavelengths must be positive")
    x = (_H * _C) / (lam * _KB * temperature_k)
    # Guard the exponential: above x ~ 700 the Wien tail is numerically zero at
    # double precision, and expm1 would overflow rather than underflow to it.
    out = np.zeros_like(lam)
    finite = x < 700.0
    out[finite] = (2.0 * np.pi * _H * _C**2 / lam[finite] ** 5) / np.expm1(x[finite])
    return out


def band_integrated_exitance(
    grid: SpectralGrid, temperature_k: float, *, samples: int = 64
) -> np.ndarray:
    """Planck exitance integrated over each bin of ``grid``, in W m^-2.

    Bin integration uses Gauss-Legendre quadrature in wavelength rather than a
    midpoint value times the bin width: the Planck function is strongly curved
    across a broad band, and the midpoint error would appear in the energy budget
    as a spurious residual that no physics accounts for.
    """
    nodes, weights = _quadrature(samples)
    lo = grid.edges[:-1][:, None]
    hi = grid.edges[1:][:, None]
    half = 0.5 * (hi - lo)
    mid = 0.5 * (hi + lo)
    lam = mid + half * nodes[None, :]
    values = planck_spectral_exitance(lam.reshape(-1), temperature_k).reshape(lam.shape)
    return np.sum(values * weights[None, :], axis=1) * half[:, 0]


def normalized_band_weights(grid: SpectralGrid, temperature_k: float) -> np.ndarray:
    """Fraction of a black body's total emission falling in each bin.

    The weights sum to the fraction of ``sigma T^4`` inside the grid's support,
    which is below one for a grid that does not span the whole spectrum.  The
    shortfall is returned in the weights rather than normalised away, so a caller
    scaling a measured total irradiance onto a truncated grid can see -- and
    report -- how much of the spectrum its grid omits.
    """
    return band_integrated_exitance(grid, temperature_k) / (_SIGMA * temperature_k**4)


def band_integrated_exitance_profile(
    grid: SpectralGrid, temperatures: np.ndarray, *, samples: int = 64
) -> np.ndarray:
    """``band_integrated_exitance`` for many temperatures at once.

    Returns shape ``(n_temperatures, n_bands)``.  The quadrature is identical to
    the scalar routine -- same nodes, same weights, same summation order per band
    -- so this is a vectorisation and not a different approximation.  It exists
    because a column solve evaluates the Planck integral for every layer at every
    step, and doing that one layer at a time dominates the run time (CMP-009
    permits optimisation only where physical meaning is unchanged).
    """
    T = np.atleast_1d(np.asarray(temperatures, dtype=float))
    nodes, weights = _quadrature(samples)
    lo = grid.edges[:-1][:, None]
    hi = grid.edges[1:][:, None]
    half = 0.5 * (hi - lo)
    lam = 0.5 * (hi + lo) + half * nodes[None, :]  # (n_bands, samples)

    out = np.zeros((T.size, grid.n_bins))
    positive = T > 0.0
    if not np.any(positive):
        return out

    x = (_H * _C) / (lam[None, :, :] * _KB * T[positive][:, None, None])
    radiance = np.zeros_like(x)
    finite = x < 700.0
    lam_broadcast = np.broadcast_to(lam[None, :, :], x.shape)
    radiance[finite] = (2.0 * np.pi * _H * _C**2 / lam_broadcast[finite] ** 5) / np.expm1(x[finite])
    out[positive] = np.sum(radiance * weights[None, None, :], axis=2) * half[:, 0][None, :]
    return out


def band_exitance_profile(
    grid: SpectralGrid, temperatures: np.ndarray, *, normalize: bool = True
) -> np.ndarray:
    """``band_exitance`` for many temperatures at once, shape ``(n_temps, n_bands)``."""
    raw = band_integrated_exitance_profile(grid, temperatures)
    if not normalize:
        return raw
    T = np.atleast_1d(np.asarray(temperatures, dtype=float))
    totals = np.sum(raw, axis=1)
    target = _SIGMA * np.where(T > 0.0, T, 0.0) ** 4
    scale = np.divide(target, totals, out=np.zeros_like(totals), where=totals > 0.0)
    return raw * scale[:, None]


def band_exitance(
    grid: SpectralGrid, temperature_k: float, *, normalize: bool = True
) -> np.ndarray:
    """Thermal emission per bin in W m^-2, optionally normalised to ``sigma T^4``.

    A finite grid never spans the whole Planck curve, so the raw band integrals
    sum to slightly less than ``sigma T^4``.  With ``normalize`` the shortfall is
    redistributed across the represented bands so the emitted total is exactly
    the grey-body total; without it the shortfall is simply absent.

    Normalisation is the default because an energy budget must close against the
    physical emission, not against an artifact of the grid's edges -- but it *is*
    an approximation: it places the missing far-infrared energy in bands where
    the atmosphere may have different opacity.  ``band_exitance_capture`` reports
    the fraction being redistributed so a run can state it.
    """
    raw = band_integrated_exitance(grid, temperature_k)
    if not normalize:
        return raw
    total = float(np.sum(raw))
    if total <= 0.0:
        return raw
    return raw * (greybody_exitance(temperature_k) / total)


def band_exitance_capture(grid: SpectralGrid, temperature_k: float) -> float:
    """Fraction of ``sigma T^4`` that actually falls inside the grid's support.

    One minus this is what ``band_exitance`` redistributes, and it is the
    structural uncertainty a spectrally truncated configuration carries.
    """
    if temperature_k <= 0.0:
        return 1.0
    return float(np.sum(band_integrated_exitance(grid, temperature_k))) / (
        _SIGMA * temperature_k**4
    )


def greybody_exitance(temperature_k: float, emissivity: float = 1.0) -> float:
    """EQ-GREYBODY-EMISSION: ``M = epsilon sigma T^4`` in W m^-2."""
    if not 0.0 <= emissivity <= 1.0:
        raise ConfigurationError(
            f"emissivity {emissivity} is outside the declared validity domain [0, 1]"
        )
    if temperature_k < 0.0:
        raise ConfigurationError(f"negative temperature {temperature_k} K is unphysical")
    return emissivity * _SIGMA * temperature_k**4


def equilibrium_temperature(
    irradiance_normal: float, albedo: float, emissivity: float = 1.0, geometry_factor: float = 0.25
) -> float:
    """EQ-RADIATIVE-EQUILIBRIUM-TEMPERATURE, the Stage 1 analytic target."""
    if emissivity <= 0.0:
        raise ConfigurationError("equilibrium temperature is undefined for zero emissivity")
    absorbed = (1.0 - albedo) * irradiance_normal * geometry_factor
    return float((absorbed / (emissivity * _SIGMA)) ** 0.25)


def relaxation_timescale(
    areal_heat_capacity: float, temperature_k: float, emissivity: float = 1.0
) -> float:
    """Linearised e-folding time ``tau = C / (4 epsilon sigma T^3)`` in seconds.

    The measured transient of the slab must match this to first order in the
    perturbation.  It is a prediction of the reduced equation, not an input, and
    Stage 5 exists to test whether the real system has one such timescale at all.
    """
    denominator = 4.0 * emissivity * _SIGMA * temperature_k**3
    if denominator <= 0.0:
        raise ConfigurationError(
            "relaxation timescale is undefined at zero emissivity or temperature"
        )
    return areal_heat_capacity / denominator
