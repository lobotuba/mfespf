"""Water vapour absorbs sunlight: the k-distributions of Chou and Suarez (1999).

Why this exists
---------------

The reference column's surface netted 183.5 W m^-2 of radiation against an
observed ~104 (Trenberth, Fasullo and Kiehl 2009) and shed the excess as a water
cycle running one and a half times too fast. Earth's atmosphere absorbs about
78 W m^-2 of sunlight and this one absorbed none: its vapour bands stopped at
4 microns. Most of the real absorption is by water vapour in the near infrared,
in bands at 0.94, 1.14, 1.38, 1.87 and 2.7 microns.

Where the numbers come from
---------------------------

Chou, M.-D. and Suarez, M. J. (1999), *A Solar Radiation Parameterization for
Atmospheric Studies*, NASA/TM-1999-104606 vol. 15, a public-domain US Government
work. Its tables were derived by its authors from line-by-line calculations and
have not been adjusted here; nothing in this module is fitted to an energy
budget, which is what leaves the budget free to test it
(``docs/protocols/solar-vapour-absorption.md``).

- Table 2: ten-term k-distributions in 0.70-1.22, 1.22-2.27 and 2.27-10 microns.
  ``h`` is each term's share of the *total* solar flux; the three bands carry
  0.32055, 0.16536 and 0.04335 of it. ``k`` is in g^-1 cm^2 at 300 hPa, 240 K.
- Table 3: a single effective ``k = 0.00075 g^-1 cm^2`` in 0.4-0.7 microns.
- Equations 3.2-3.5: the absorber amount is scaled by
  ``(p / 300 hPa)^0.8 [1 + 0.00135 (T - 240 K)]``, which folds the pressure
  broadening of the line wings and the weak temperature dependence into the
  path.

How it enters the grid
----------------------

A k-distribution sorts a band's spectrum by absorption strength. Within one of
the source's bands, a sub-band per term, sized so that it receives its term's
share of the band's sunlight and carrying that term's ``k``, is the same
transmission for anything that does not vary with wavelength inside the band --
here the surface albedo, which is flat below 4 microns, and a thermal emission
that is negligible there. So the solver is unchanged: it sees more bands, not a
new method. The bands above the split, where the thermal calibrations live, are
left exactly as they were.
"""

from __future__ import annotations

import numpy as np

from mfespf.contracts.spectral import SpectralCoordinate, SpectralGrid
from mfespf.core.errors import ConfigurationError

__all__ = [
    "SOLAR_VAPOUR_BANDS",
    "SOLAR_VAPOUR_REFERENCE_PRESSURE_PA",
    "SOLAR_VAPOUR_REFERENCE_TEMPERATURE_K",
    "scaled_vapour_path",
    "solar_vapour_bands",
    "solar_vapour_optical_depth",
]

#: g^-1 cm^2 to m^2 kg^-1.
_CM2_PER_G = 0.1

#: Table 2, the ten absorption coefficients shared by the three infrared bands,
#: g^-1 cm^2 at 300 hPa and 240 K.
_K_INFRARED = (0.0010, 0.0133, 0.0422, 0.1334, 0.4217, 1.3340, 5.6230, 31.620, 177.80, 1000.0)

#: Table 2, flux weights: fractions of the total solar flux, per band.
_H_BAND_9 = (0.20673, 0.03497, 0.03011, 0.02260, 0.01336, 0.00696, 0.00441, 0.00115, 0.00026, 0.0)
_H_BAND_10 = (
    0.08236,
    0.01157,
    0.01133,
    0.01143,
    0.01240,
    0.01258,
    0.01381,
    0.00650,
    0.00244,
    0.00094,
)
_H_BAND_11 = (
    0.01074,
    0.00360,
    0.00411,
    0.00421,
    0.00389,
    0.00326,
    0.00499,
    0.00465,
    0.00245,
    0.00145,
)

#: Table 3: the effective coefficient for 0.4-0.7 microns, g^-1 cm^2.
_K_VISIBLE = 0.00075

#: The source's bands this module represents: (lower m, upper m, ks, weights).
#: Band 11 runs to 10 microns in the source; here it stops where the grid's
#: thermal bands begin, and its k-distribution is applied over that part.
SOLAR_VAPOUR_BANDS: tuple[tuple[float, float, tuple[float, ...], tuple[float, ...]], ...] = (
    (0.40e-6, 0.70e-6, (_K_VISIBLE,), (1.0,)),
    (0.70e-6, 1.22e-6, _K_INFRARED, _H_BAND_9),
    (1.22e-6, 2.27e-6, _K_INFRARED, _H_BAND_10),
    (2.27e-6, 10.0e-6, _K_INFRARED, _H_BAND_11),
)

SOLAR_VAPOUR_REFERENCE_PRESSURE_PA = 3.0e4
SOLAR_VAPOUR_REFERENCE_TEMPERATURE_K = 240.0
_PRESSURE_EXPONENT = 0.8
_TEMPERATURE_COEFFICIENT = 0.00135

_H = 6.62607015e-34
_C = 2.99792458e8
_KB = 1.380649e-23


def _planck_cumulative(
    lower: float, upper: float, temperature: float, n: int = 20001
) -> tuple[np.ndarray, np.ndarray]:
    """Wavelength mesh and cumulative black-body exitance over [lower, upper]."""
    wl = np.geomspace(lower, upper, n)
    x = _H * _C / (wl * _KB * temperature)
    density = 1.0 / (wl**5 * np.expm1(x))
    cumulative = np.concatenate(
        ([0.0], np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(wl)))
    )
    return wl, cumulative / cumulative[-1]


def solar_vapour_bands(
    grid: SpectralGrid, *, solar_temperature: float = 5772.0
) -> tuple[SpectralGrid, np.ndarray]:
    """Rebuild a grid's solar region around the k-distributions.

    Returns the new grid and ``k``, the reference absorption coefficient per
    band in m^2 kg^-1 (zero where water vapour does not absorb sunlight here).
    The grid's edges below 0.4 microns and from the first edge at or beyond
    2.27 microns upward are kept exactly; everything between is replaced by the
    source's band edges and, inside each infrared band, one sub-band per term
    with non-zero weight. Sub-band edges are placed on the cumulative spectrum of
    a black body at ``solar_temperature`` -- the same shape the configured solar
    source uses -- so that each receives its term's share of the band's sunlight.
    """
    if grid.coordinate is not SpectralCoordinate.WAVELENGTH:
        raise ConfigurationError("solar vapour bands are built on a wavelength grid")
    edges = np.asarray(grid.edges, dtype=float)
    split = edges[edges >= 2.27e-6]
    if split.size == 0 or split[0] > 4.5e-6:
        raise ConfigurationError(
            "the grid needs an edge between 2.27 and 4.5 microns so that the thermal "
            "bands above it can be kept as they are"
        )
    top = float(split[0])

    new_edges: list[float] = [float(e) for e in edges if e < 0.40e-6]
    ks: list[float] = [0.0] * max(len(new_edges) - 1, 0)
    if new_edges and new_edges[-1] < 0.40e-6:
        new_edges.append(0.40e-6)
        ks.append(0.0)
    elif not new_edges:
        new_edges = [0.40e-6]

    for lower, upper, k_terms, weights in SOLAR_VAPOUR_BANDS:
        upper = min(upper, top)
        terms = [(k, h) for k, h in zip(k_terms, weights, strict=True) if h > 0.0]
        total = sum(h for _, h in terms)
        wl, cumulative = _planck_cumulative(lower, upper, solar_temperature)
        shares = np.cumsum([h / total for _, h in terms])[:-1]
        interior = np.interp(shares, cumulative, wl)
        for edge in (*interior, upper):
            new_edges.append(float(edge))
        ks.extend(k * _CM2_PER_G for k, _ in terms)
        if upper >= top:
            break

    kept = [float(e) for e in edges if e > top]
    new_edges.extend(kept)
    ks.extend([0.0] * len(kept))

    new_grid = SpectralGrid(
        grid_id=f"{grid.grid_id}+solar_vapour",
        coordinate=grid.coordinate,
        coordinate_units=grid.coordinate_units,
        edges=np.asarray(new_edges),
        integration_convention=grid.integration_convention,
        density_units=grid.density_units,
        resolution_identifier=f"{grid.resolution_identifier}+cs99-vapour",
    )
    k = np.asarray(ks, dtype=float)
    if k.size != new_grid.n_bins:
        raise ConfigurationError("internal: k does not match the rebuilt grid")
    return new_grid, k


def scaled_vapour_path(
    vapour_path: np.ndarray, pressure: np.ndarray, temperature: np.ndarray
) -> np.ndarray:
    """The source's scaled absorber amount per layer, kg m^-2 (its eq. 3.5)."""
    w = np.asarray(vapour_path, dtype=float)
    p = np.asarray(pressure, dtype=float)
    t = np.asarray(temperature, dtype=float)
    if np.any(w < 0.0):
        raise ConfigurationError("a vapour path cannot be negative")
    if np.any(p <= 0.0):
        raise ConfigurationError("layer pressures must be positive")
    return (
        w
        * (p / SOLAR_VAPOUR_REFERENCE_PRESSURE_PA) ** _PRESSURE_EXPONENT
        * (1.0 + _TEMPERATURE_COEFFICIENT * (t - SOLAR_VAPOUR_REFERENCE_TEMPERATURE_K))
    )


def solar_vapour_optical_depth(
    k: np.ndarray, vapour_path: np.ndarray, pressure: np.ndarray, temperature: np.ndarray
) -> np.ndarray:
    """Absorption optical depth of the vapour for sunlight, shaped (layers, bands)."""
    return scaled_vapour_path(vapour_path, pressure, temperature)[:, None] * np.asarray(k)[None, :]
