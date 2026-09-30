"""The 15-micron CO2 band: Chou, Suarez, Liang and Yan (2001).

Why this exists
---------------

The reference column's only non-vapour infrared absorber was a flat caricature,
``window_absorption(kappa_opaque=1e-4)``, and it cools the mid-latitude summer
stratopause at 5.8 C/day against a line-by-line 12. The calibrated band in
``greenhouse.py`` has absorption that does not depend on pressure, and it cools
the same stratopause at 88 C/day. Real lines narrow as the pressure falls, so a
band's emissivity in thin air depends on the pressure it is at
(``docs/protocols/co2-band.md``).

Where the numbers come from
---------------------------

Chou, M.-D., Suarez, M. J., Liang, X.-Z. and Yan, M. M.-H. (2001), *A Thermal
Infrared Radiation Parameterization for Atmospheric Studies*, NASA/TM-2001-104606
vol. 15/19, a public-domain US Government work, section 4.5 and Tables 3 and 10.
They were fitted by the authors to line-by-line calculations with the 1996
HITRAN, and none is adjusted here.

- Band 3, 540-800 cm^-1: a centre (620-720) and two wings (540-620, 720-800).
- Six k-terms each, ``k_n = k_1 8^(n-1)``, with Planck-weighted fractions that
  sum to one over the band. ``k_1 / mu`` is given with Chou's diffusivity
  factor 1.66 folded in; the solver applies its own 1.66, so the vertical
  coefficient is ``(k_1 / mu) / 1.66``.
- The CO2 amount is scaled for pressure and temperature (their eq. 4.1):
  ``u (p / p_r)^m [1 + alpha (T - 250) + beta (T - 250)^2]``. The centre is
  referred to 30 hPa with ``m = 0.85``, the wings to 300 hPa with ``m = 0.5``.

How it enters the grid
----------------------

As the solar vapour's terms did: one sub-band per term, each sized to its
term's share of its sub-band's Planck flux at 250 K. The centre's fractions sum
to 0.389 and the Planck share of 620-720 cm^-1 in the band at 250 K is 0.3884,
which is the check that the fractions are what this takes them to be.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from mfespf.contracts.spectral import SpectralCoordinate, SpectralGrid
from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError

__all__ = [
    "CO2_CENTRE",
    "CO2_WINGS",
    "CO2Terms",
    "co2_band_grid",
    "co2_layer_amount",
    "co2_optical_depth",
]

_H = 6.62607015e-34
_C = 2.99792458e8
_KB = 1.380649e-23

#: Chou's diffusivity factor, folded into their k_1 / mu.
_CHOU_DIFFUSIVITY = 1.66
#: The temperature the Planck-weighted fractions and the scaling refer to.
_REFERENCE_TEMPERATURE_K = 250.0


@dataclass(frozen=True)
class CO2Terms:
    """One part of the band: Tables 3 and 10 of Chou et al. (2001)."""

    #: k_1 / mu, (cm-atm)^-1, with Chou's diffusivity 1.66 folded in.
    k1_over_mu: float
    #: Each term's fraction of the band's Planck flux at 250 K.
    fractions: tuple[float, ...]
    #: The pressure the scaling refers to, Pa, and its exponent.
    reference_pressure: float
    m: float
    #: The temperature scaling, per K and per K^2.
    alpha: float
    beta: float
    #: The wavenumber ranges its terms fill, cm^-1.
    ranges_cm: tuple[tuple[float, float], ...]


CO2_CENTRE = CO2Terms(
    k1_over_mu=2.66e-3,
    fractions=(0.0766, 0.1372, 0.1189, 0.0335, 0.0169, 0.0059),
    reference_pressure=3000.0,
    m=0.85,
    alpha=0.0042,
    beta=2.00e-5,
    ranges_cm=((620.0, 720.0),),
)
CO2_WINGS = CO2Terms(
    k1_over_mu=2.66e-5,
    fractions=(0.1395, 0.1407, 0.1549, 0.1357, 0.0182, 0.0220),
    reference_pressure=30000.0,
    m=0.50,
    alpha=0.0179,
    beta=1.02e-4,
    ranges_cm=((540.0, 620.0), (720.0, 800.0)),
)
_RATIO = 8.0


def _planck_cumulative_wavenumber(
    low: float, high: float, n: int = 20001
) -> tuple[np.ndarray, np.ndarray]:
    """Wavenumber mesh (cm^-1) and cumulative Planck flux at 250 K, normalised."""
    nu = np.linspace(low, high, n)
    x = nu * 100.0
    density = x**3 / np.expm1(_H * _C * x / (_KB * _REFERENCE_TEMPERATURE_K))
    cumulative = np.concatenate(
        ([0.0], np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(nu)))
    )
    return nu, cumulative / cumulative[-1]


def co2_band_grid(grid: SpectralGrid) -> tuple[SpectralGrid, np.ndarray, np.ndarray]:
    """Rebuild a wavelength grid's 12.5-18.5 micron region around the band.

    Returns the grid, the vertical absorption coefficient per bin in
    (cm-atm)^-1 (zero outside the band), and each bin's kind: 0 outside, 1 the
    centre, 2 the wings. Edges outside 540-800 cm^-1 are kept exactly.
    """
    if grid.coordinate is not SpectralCoordinate.WAVELENGTH:
        raise ConfigurationError("the CO2 band is built on a wavelength grid")
    edges = np.asarray(grid.edges, dtype=float)
    low_m, high_m = 1.0e-2 / 800.0, 1.0e-2 / 540.0
    if edges[0] > low_m or edges[-1] < high_m:
        raise ConfigurationError("the grid must span 12.5 to 18.5 microns to hold the CO2 band")

    pieces: list[tuple[float, float, float, int]] = []  # (lower m, upper m, k, kind)
    for kind, spec in ((1, CO2_CENTRE), (2, CO2_WINGS)):
        fractions = np.asarray(spec.fractions, dtype=float)
        shares = np.cumsum(fractions / fractions.sum())[:-1]
        ks = spec.k1_over_mu * _RATIO ** np.arange(fractions.size) / _CHOU_DIFFUSIVITY
        for low_cm, high_cm in spec.ranges_cm:
            nu, cumulative = _planck_cumulative_wavenumber(low_cm, high_cm)
            interior = np.interp(shares, cumulative, nu)
            bounds = np.concatenate(([low_cm], interior, [high_cm]))
            # Term n occupies [bounds[n], bounds[n+1]] in wavenumber; in
            # wavelength the order reverses, which the sort below handles.
            for n in range(fractions.size):
                pieces.append((1.0e-2 / bounds[n + 1], 1.0e-2 / bounds[n], float(ks[n]), kind))

    pieces.sort()
    kept_low = [float(e) for e in edges if e < low_m and not np.isclose(e, low_m)]
    kept_high = [float(e) for e in edges if e > high_m and not np.isclose(e, high_m)]
    new_edges = [*kept_low, *(p[0] for p in pieces), pieces[-1][1], *kept_high]
    k = [0.0] * len(kept_low) + [p[2] for p in pieces] + [0.0] * len(kept_high)
    kinds = [0] * len(kept_low) + [p[3] for p in pieces] + [0] * len(kept_high)
    new_grid = SpectralGrid(
        grid_id=f"{grid.grid_id}+co2_chou2001",
        coordinate=grid.coordinate,
        coordinate_units=grid.coordinate_units,
        edges=np.asarray(new_edges),
        integration_convention=grid.integration_convention,
        density_units=grid.density_units,
        resolution_identifier=f"{grid.resolution_identifier}+chou2001-co2",
    )
    k_arr, kind_arr = np.asarray(k), np.asarray(kinds, dtype=int)
    if k_arr.size != new_grid.n_bins:
        raise ConfigurationError("internal: the CO2 coefficients do not match the rebuilt grid")
    return new_grid, k_arr, kind_arr


def co2_layer_amount(co2_ppmv: float, layer_mass: np.ndarray) -> np.ndarray:
    """Each layer's CO2 in cm-atm at STP, from its mixing ratio and its air."""
    k_b = value_of("boltzmann")
    molecules_per_kg = value_of("gas_constant_dry_air") / k_b
    cm_atm = value_of("standard_surface_pressure") / (k_b * 273.15) * 0.01
    return (
        float(co2_ppmv) * 1.0e-6 * np.asarray(layer_mass, dtype=float) * molecules_per_kg / cm_atm
    )


def co2_optical_depth(
    k: np.ndarray,
    kind: np.ndarray,
    *,
    co2_ppmv: float,
    layer_mass: np.ndarray,
    pressure: np.ndarray,
    temperature: np.ndarray,
) -> np.ndarray:
    """Vertical absorption optical depth of CO2, shaped (layers, bins)."""
    amount = co2_layer_amount(co2_ppmv, layer_mass)
    p = np.asarray(pressure, dtype=float)
    dt = np.asarray(temperature, dtype=float) - _REFERENCE_TEMPERATURE_K
    tau = np.zeros((amount.size, np.asarray(k).size))
    for code, spec in ((1, CO2_CENTRE), (2, CO2_WINGS)):
        columns = np.asarray(kind) == code
        if not np.any(columns):
            continue
        scale = (p / spec.reference_pressure) ** spec.m
        scale = scale * (1.0 + spec.alpha * dt + spec.beta * dt**2)
        tau[:, columns] = (amount * scale)[:, None] * np.asarray(k)[None, columns]
    return tau
