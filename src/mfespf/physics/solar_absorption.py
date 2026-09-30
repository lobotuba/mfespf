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
from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError

__all__ = [
    "OZONE_BANDS",
    "OZONE_PROFILE_HPA_PPMV",
    "SOLAR_VAPOUR_BANDS",
    "SOLAR_VAPOUR_REFERENCE_PRESSURE_PA",
    "SOLAR_VAPOUR_REFERENCE_TEMPERATURE_K",
    "ozone_bands",
    "ozone_layer_path",
    "ozone_optical_depth",
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


def _measured_cumulative(
    lower: float, upper: float, n: int = 20001
) -> tuple[np.ndarray, np.ndarray]:
    """Wavelength mesh and cumulative measured solar irradiance over [lower, upper]."""
    from mfespf.physics.solar_spectrum import hsrs_cumulative

    wl = np.geomspace(lower, upper, n)
    cumulative = hsrs_cumulative(wl) - hsrs_cumulative(np.array([lower]))[0]
    return wl, cumulative / cumulative[-1]


def solar_vapour_bands(
    grid: SpectralGrid, *, solar_temperature: float = 5772.0, measured_spectrum: bool = False
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
    With ``measured_spectrum`` they are placed on the measured spectrum instead,
    which is then the solar source's shape too.
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
        wl, cumulative = (
            _measured_cumulative(lower, upper)
            if measured_spectrum
            else _planck_cumulative(lower, upper, solar_temperature)
        )
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


# --- ozone ----------------------------------------------------------------

#: Chou and Suarez (1999) Table 1 and section 3.5: the bands in which ozone
#: absorbs sunlight, (lower m, upper m, effective coefficient in (cm-atm)^-1 at
#: STP). Band 2 is two disjoint ranges sharing one coefficient; band 8's
#: coefficient already carries ozone's weak near-infrared absorption (their eq.
#: 3.14). No pressure or temperature scaling, as in the source.
OZONE_BANDS: tuple[tuple[float, float, float], ...] = (
    (0.175e-6, 0.225e-6, 30.47),
    (0.225e-6, 0.245e-6, 187.24),
    (0.245e-6, 0.260e-6, 301.92),
    (0.260e-6, 0.280e-6, 187.24),
    (0.280e-6, 0.295e-6, 42.83),
    (0.295e-6, 0.310e-6, 7.09),
    (0.310e-6, 0.320e-6, 1.25),
    (0.320e-6, 0.400e-6, 0.0345),
    (0.400e-6, 0.700e-6, 0.0572),
)

#: Ozone mixing ratio, global area-weighted mean of SWOOSH v2.73 combined,
#: anomaly-filled, 10-degree zonal means, 1991-2020 (Davis et al. 2016,
#: doi:10.5194/essd-8-461-2016), computed for
#: ``docs/protocols/ozone-solar-heating.md``: (pressure hPa, ppmv), 316 to 1 hPa.
OZONE_PROFILE_HPA_PPMV: tuple[tuple[float, float], ...] = (
    (316.23, 0.061),
    (261.02, 0.082),
    (215.44, 0.119),
    (177.83, 0.162),
    (146.78, 0.217),
    (121.15, 0.305),
    (100.00, 0.462),
    (82.54, 0.723),
    (68.13, 1.124),
    (56.23, 1.671),
    (46.42, 2.344),
    (38.31, 3.118),
    (31.62, 3.980),
    (26.10, 4.910),
    (21.54, 5.822),
    (17.78, 6.643),
    (14.68, 7.297),
    (12.12, 7.721),
    (10.00, 7.933),
    (8.25, 7.981),
    (6.81, 7.923),
    (5.62, 7.778),
    (4.64, 7.536),
    (3.83, 7.159),
    (3.16, 6.635),
    (2.61, 6.011),
    (2.15, 5.338),
    (1.78, 4.696),
    (1.47, 4.114),
    (1.21, 3.612),
    (1.00, 3.182),
)

#: STP is 273.15 K by definition; with the standard pressure it fixes the
#: Loschmidt number that makes a cm-atm a number of molecules.
_STP_TEMPERATURE_K = 273.15


def ozone_bands(grid: SpectralGrid) -> tuple[SpectralGrid, np.ndarray]:
    """Rebuild a grid below 0.7 microns around ozone's bands.

    Returns the new grid and ozone's coefficient per band in (cm-atm)^-1, zero
    where it does not absorb. Edges below 0.175 microns are kept, the source's
    edges from 0.175 to 0.4 replace whatever the grid had there, and edges at
    0.4 and 0.7 are added if missing; everything above 0.4 microns is kept, so
    a grid already built by ``solar_vapour_bands`` keeps its sub-bands.
    """
    if grid.coordinate is not SpectralCoordinate.WAVELENGTH:
        raise ConfigurationError("ozone bands are built on a wavelength grid")
    edges = np.asarray(grid.edges, dtype=float)
    if edges[0] > 0.175e-6 or edges[-1] < 0.7e-6:
        raise ConfigurationError("the grid must span 0.175 to 0.7 microns to hold ozone's bands")
    source = {b[0] for b in OZONE_BANDS} | {b[1] for b in OZONE_BANDS}
    kept_low = {float(e) for e in edges if e < 0.175e-6}
    kept_high = {float(e) for e in edges if e > 0.4e-6 and not np.isclose(e, 0.7e-6)}
    new_edges = np.array(sorted(kept_low | source | kept_high))
    k = np.zeros(new_edges.size - 1)
    centres = 0.5 * (new_edges[:-1] + new_edges[1:])
    for lower, upper, coefficient in OZONE_BANDS:
        k[(centres > lower) & (centres < upper)] = coefficient
    new_grid = SpectralGrid(
        grid_id=f"{grid.grid_id}+ozone",
        coordinate=grid.coordinate,
        coordinate_units=grid.coordinate_units,
        edges=new_edges,
        integration_convention=grid.integration_convention,
        density_units=grid.density_units,
        resolution_identifier=f"{grid.resolution_identifier}+cs99-ozone",
    )
    return new_grid, k


def _ozone_mixing_ratio(pressure_pa: np.ndarray) -> np.ndarray:
    """SWOOSH's mixing ratio at each pressure, mol mol^-1, held below 316 hPa."""
    p = np.array([row[0] for row in OZONE_PROFILE_HPA_PPMV])[::-1] * 100.0
    v = np.array([row[1] for row in OZONE_PROFILE_HPA_PPMV])[::-1] * 1.0e-6
    return np.interp(np.log(np.asarray(pressure_pa, dtype=float)), np.log(p), v)


def ozone_layer_path(
    layer_pressure: np.ndarray, interface_pressure: np.ndarray, layer_mass: np.ndarray
) -> np.ndarray:
    """Ozone in each layer, cm-atm at STP, from the SWOOSH profile.

    Each layer carries the mixing ratio at its centre times its air. The top
    layer also carries the column between its top interface and 1 hPa, where
    SWOOSH ends; what lies above 1 hPa (about a dobson unit) is left out.
    Molecules per kilogram of dry air are ``R_d / k_B``, and a cm-atm is the
    Loschmidt number times a centimetre, so no molar mass is assumed.
    """
    p = np.asarray(layer_pressure, dtype=float)
    p_i = np.asarray(interface_pressure, dtype=float)
    m = np.asarray(layer_mass, dtype=float)
    if p_i.size != p.size + 1 or m.size != p.size:
        raise ConfigurationError("pressures and masses must describe the same layers")
    k_b = value_of("boltzmann")
    molecules_per_kg = value_of("gas_constant_dry_air") / k_b
    cm_atm = value_of("standard_surface_pressure") / (k_b * _STP_TEMPERATURE_K) * 0.01
    path = _ozone_mixing_ratio(p) * m * molecules_per_kg / cm_atm
    top = float(p_i[-1])
    if top > 100.0:
        # Gravity as the grid has it: the top layer's own dp / m.
        gravity = float((p_i[-2] - p_i[-1]) / m[-1])
        above = np.geomspace(100.0, top, 400)
        column = float(np.trapezoid(_ozone_mixing_ratio(above), above)) / gravity
        path[-1] += column * molecules_per_kg / cm_atm
    return path


def ozone_optical_depth(k: np.ndarray, path: np.ndarray) -> np.ndarray:
    """Absorption optical depth of ozone for sunlight, shaped (layers, bands)."""
    return np.asarray(path, dtype=float)[:, None] * np.asarray(k, dtype=float)[None, :]
