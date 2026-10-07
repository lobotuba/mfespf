"""Water vapour's longwave absorption: Chou, Suarez, Liang and Yan (2001).

Why this exists
---------------

The vapour was a two-wing envelope with one fitted decay, evaluated at the
centre of coarse bins (``greenhouse.vapour_mass_absorption``). On the AFGL
mid-latitude summer atmosphere that left the outgoing longwave 29.7 W m^-2
below line-by-line, almost all of it at 1100-1380 and 340-540 cm^-1, where one
wide bin applied an opaque band wing across partly transparent regions
(``docs/protocols/longwave-band-comparison.md``). The fit to all-sky OLR then
made the total right by compensation, standing in for the clouds as well.

Where the numbers come from
---------------------------

Chou et al. (2001), NASA/TM-2001-104606 vol. 19, the source of the CO2 band:

- Table 1: the bands.
- Table 4: the line k-distribution for bands 1, 2 and 4-9: k_1/mu, the ratio
  eta (k_n = k_1 eta^(n-1)) and the Planck-weighted fractions Delta-g.
- Table 10: the same for the three sub-bands of band 3.
- Table 3: the scaling of the amount, w = q dp/g (p / 500 hPa)^m h(theta).
- Table 9 and eq. 4.19-4.24: the self-broadened continuum in bands 3-7.

Nothing is fitted or adjusted. Chou's k/mu carries his diffusivity 1.66, the
solver applies its own, so each vertical coefficient is (k/mu)/1.66.

How it enters the grid
----------------------

As the CO2 and ozone terms did: one sub-band per term, sized to its fraction of
the band's Planck flux at 250 K. Where CO2 (band 3) or ozone (band 5) shares the
band, their transmittances multiply (Chou eq. 4.27). That is represented exactly
by splitting each CO2 or ozone sub-band again by the water terms, so each bin
carries one water term and one other, and holds the product of their fractions.
Everything at 4 microns and longer is rebuilt this way, from 2500 down to
10 cm^-1; everything shorter is kept edge for edge
(``docs/protocols/water-vapour-longwave.md``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from mfespf.contracts.spectral import SpectralCoordinate, SpectralGrid
from mfespf.core.errors import ConfigurationError
from mfespf.physics.co2_band import CO2_CENTRE, CO2_WINGS
from mfespf.physics.ozone_infrared import O3_IR_WEIGHTS
from mfespf.physics.water_vapour_tables import (
    WATER_TABLE_LOG10_P_HPA,
    WATER_TABLES,
    TabulatedWater,
)

__all__ = [
    "CHOU_WATER_BANDS",
    "ChouLongwave",
    "WaterBand",
    "chou_longwave_grid",
    "water_vapour_optical_depth",
]

_H = 6.62607015e-34
_C = 2.99792458e8
_KB = 1.380649e-23
_REFERENCE_TEMPERATURE_K = 250.0
_CHOU_DIFFUSIVITY = 1.66
#: Chou's reference pressure for the water scaling (Table 3), Pa.
_WATER_REFERENCE_PRESSURE = 50000.0
#: g^-1 cm^2 times kg m^-2: one kg m^-2 is 0.1 g cm^-2.
_KG_M2_TO_G_CM2 = 0.1
#: The longwave region rebuilt, cm^-1: the model's 4-micron shortwave boundary
#: down to its 1 mm long-wavelength edge.
_LONGWAVE_CM = (10.0, 2500.0)


@dataclass(frozen=True)
class WaterBand:
    """One of Chou's bands or sub-bands, for water vapour."""

    name: str
    range_cm: tuple[float, float]
    #: k_1 / mu in g^-1 cm^2, the ratio eta, and the Planck-weighted fractions.
    k1_over_mu: float
    eta: float
    delta_g: tuple[float, ...]
    #: Table 3: the scaling's exponent and temperature coefficients.
    m: float
    alpha: float
    beta: float
    #: Table 9: continuum k_c / mu in g^-1 cm^2, zero where there is none.
    continuum_over_mu: float = 0.0


_BAND3_SCALING = {"m": 1.0, "alpha": 0.0167, "beta": 8.54e-5}

#: Tables 1, 3, 4, 9 and 10 of Chou et al. (2001), as published.
CHOU_WATER_BANDS: tuple[WaterBand, ...] = (
    WaterBand(
        "1",
        (0.0, 340.0),
        2.90e1,
        6,
        (0.2747, 0.2717, 0.2752, 0.1177, 0.0352, 0.0255),
        1.0,
        0.0021,
        -1.01e-5,
    ),
    WaterBand(
        "2",
        (340.0, 540.0),
        4.17e-1,
        6,
        (0.1521, 0.3974, 0.1778, 0.1826, 0.0374, 0.0527),
        1.0,
        0.0140,
        5.57e-5,
    ),
    WaterBand(
        "3a",
        (540.0, 620.0),
        1.33e-2,
        8,
        (0.1782, 0.0593, 0.0215, 0.0068, 0.0022, 0.0),
        **_BAND3_SCALING,
        continuum_over_mu=109.6,
    ),
    WaterBand(
        "3b",
        (620.0, 720.0),
        1.33e-2,
        8,
        (0.0923, 0.1675, 0.0923, 0.0187, 0.0178, 0.0),
        **_BAND3_SCALING,
        continuum_over_mu=54.8,
    ),
    WaterBand(
        "3c",
        (720.0, 800.0),
        1.33e-2,
        8,
        (0.0, 0.1083, 0.1581, 0.0455, 0.0274, 0.0041),
        **_BAND3_SCALING,
        continuum_over_mu=27.4,
    ),
    WaterBand(
        "4",
        (800.0, 980.0),
        5.25e-4,
        6,
        (0.4654, 0.2991, 0.1343, 0.0646, 0.0226, 0.0140),
        1.0,
        0.0302,
        2.96e-4,
        15.8,
    ),
    WaterBand(
        "5",
        (980.0, 1100.0),
        5.25e-4,
        6,
        (0.5543, 0.2723, 0.1131, 0.0443, 0.0160, 0.0),
        1.0,
        0.0307,
        2.86e-4,
        9.40,
    ),
    WaterBand(
        "6",
        (1100.0, 1215.0),
        9.37e-3,
        8,
        (0.5955, 0.2693, 0.0953, 0.0335, 0.0064, 0.0),
        0.77,
        0.0195,
        1.11e-4,
        7.75,
    ),
    WaterBand(
        "7",
        (1215.0, 1380.0),
        4.72e-2,
        9,
        (0.1958, 0.3469, 0.3147, 0.1013, 0.0365, 0.0048),
        0.5,
        0.0152,
        7.61e-5,
        8.78,
    ),
    WaterBand(
        "8",
        (1380.0, 1900.0),
        1.32e0,
        6,
        (0.0740, 0.1636, 0.4174, 0.1783, 0.1101, 0.0566),
        1.0,
        0.0008,
        -3.52e-6,
    ),
    WaterBand(
        "9",
        (1900.0, 3000.0),
        5.25e-4,
        16,
        (0.1437, 0.2197, 0.3185, 0.2351, 0.0647, 0.0183),
        1.0,
        0.0096,
        1.64e-5,
    ),
)


@dataclass(frozen=True)
class ChouLongwave:
    """Per-bin coefficients on the rebuilt grid. Shortwave bins carry nothing.

    ``water_k`` is the vertical line coefficient of the bin's water term, in
    g^-1 cm^2; ``water_band`` indexes ``CHOU_WATER_BANDS`` (-1 for none) for the
    scaling; ``continuum_k`` is the vertical continuum coefficient. ``co2_k`` and
    ``co2_kind`` are in the form ``co2_band.co2_optical_depth`` takes, and
    ``ozone_term`` in the form ``ozone_infrared.ozone_infrared_optical_depth`` takes.
    ``water_term`` is the bin's term within its band, which the tabulated bands
    (1, 2 and 8) need to look up their pressure-dependent coefficient.
    """

    water_band: np.ndarray
    water_term: np.ndarray
    #: Bands whose lines come from Chou's tables rather than Table 4
    #: (water-vapour-longwave.md attempt 2); empty is attempt 1.
    tabulated: tuple[str, ...]
    water_k: np.ndarray
    continuum_k: np.ndarray
    co2_k: np.ndarray
    co2_kind: np.ndarray
    ozone_term: np.ndarray

    def carried(self, home: np.ndarray) -> ChouLongwave:
        """The same spec on a grid rebuilt elsewhere: each new bin takes its home bin's."""
        return ChouLongwave(
            water_band=self.water_band[home],
            water_term=self.water_term[home],
            tabulated=self.tabulated,
            water_k=self.water_k[home],
            continuum_k=self.continuum_k[home],
            co2_k=self.co2_k[home],
            co2_kind=self.co2_kind[home],
            ozone_term=self.ozone_term[home],
        )


def _planck_split(low_cm: float, high_cm: float, weights: np.ndarray) -> np.ndarray:
    """Boundaries splitting [low, high] cm^-1 into pieces holding these Planck shares."""
    nu = np.linspace(low_cm, high_cm, 4001)
    x = nu * 100.0
    density = x**3 / np.expm1(_H * _C * x / (_KB * _REFERENCE_TEMPERATURE_K))
    cumulative = np.concatenate(
        ([0.0], np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(nu)))
    )
    cumulative = cumulative / cumulative[-1]
    shares = np.cumsum(weights / weights.sum())[:-1]
    return np.concatenate(([low_cm], np.interp(shares, cumulative, nu), [high_cm]))


def _water_terms(band: WaterBand, tabulated: tuple[str, ...] = ()) -> tuple[np.ndarray, np.ndarray]:
    """The band's non-zero terms: vertical coefficients (g^-1 cm^2) and fractions.

    For a tabulated band the coefficient depends on pressure, so it is returned
    as NaN here and looked up per layer by term.
    """
    if band.name in tabulated:
        weights = np.asarray(WATER_TABLES[band.name].weights, dtype=float)
        return np.full(weights.size, np.nan), weights
    g = np.asarray(band.delta_g, dtype=float)
    k = band.k1_over_mu * band.eta ** np.arange(g.size) / _CHOU_DIFFUSIVITY
    keep = g > 0.0
    return k[keep], g[keep]


def chou_longwave_grid(
    grid: SpectralGrid, *, ozone: bool = False, tabulated: tuple[str, ...] = ()
) -> tuple[SpectralGrid, ChouLongwave]:
    """Rebuild a wavelength grid's longwave on Chou's bands, with CO2 (and ozone).

    Pieces are (low cm^-1, high cm^-1, band index, water k, CO2 k, CO2 kind,
    ozone term). The shortwave, below 4 microns, is kept edge for edge.
    """
    if grid.coordinate is not SpectralCoordinate.WAVELENGTH:
        raise ConfigurationError("Chou's longwave is built on a wavelength grid")
    edges = np.asarray(grid.edges, dtype=float)
    lo_cm, hi_cm = _LONGWAVE_CM
    if edges[0] > 1.0e-2 / hi_cm or edges[-1] < 1.0e-2 / lo_cm:
        raise ConfigurationError("the grid must span 4 microns to 1 mm to hold the longwave")
    pieces: list[tuple[float, float, int, float, float, int, int, int]] = []

    def add(index: int, low: float, high: float, co2: tuple[float, int], o3: int) -> None:
        k, g = _water_terms(CHOU_WATER_BANDS[index], tabulated)
        bounds = _planck_split(low, high, g)
        for n in range(k.size):
            pieces.append((bounds[n], bounds[n + 1], index, float(k[n]), co2[0], co2[1], o3, n))

    for index, band in enumerate(CHOU_WATER_BANDS):
        low, high = max(band.range_cm[0], lo_cm), min(band.range_cm[1], hi_cm)
        if band.name in ("3a", "3b", "3c"):
            spec, kind = (CO2_CENTRE, 1) if band.name == "3b" else (CO2_WINGS, 2)
            fractions = np.asarray(spec.fractions, dtype=float)
            k_co2 = spec.k1_over_mu * 8.0 ** np.arange(fractions.size) / _CHOU_DIFFUSIVITY
            bounds = _planck_split(low, high, fractions)
            for n in range(fractions.size):
                add(index, bounds[n], bounds[n + 1], (float(k_co2[n]), kind), -1)
        elif band.name == "5" and ozone:
            weights = np.asarray(O3_IR_WEIGHTS, dtype=float)
            bounds = _planck_split(low, high, weights)
            for n in range(weights.size):
                add(index, bounds[n], bounds[n + 1], (0.0, 0), n)
        else:
            add(index, low, high, (0.0, 0), -1)

    pieces.sort(key=lambda p: 1.0e-2 / p[1])  # by wavelength of the short edge
    long_edges = [1.0e-2 / p[1] for p in pieces] + [1.0e-2 / pieces[-1][0]]
    shortwave = [
        float(e) for e in edges if e < 1.0e-2 / hi_cm and not np.isclose(e, 1.0e-2 / hi_cm)
    ]
    new_edges = np.asarray([*shortwave, *long_edges])
    n_short = len(shortwave)
    new_grid = SpectralGrid(
        grid_id=f"{grid.grid_id}+chou2001-longwave",
        coordinate=grid.coordinate,
        coordinate_units=grid.coordinate_units,
        edges=new_edges,
        integration_convention=grid.integration_convention,
        density_units=grid.density_units,
        resolution_identifier=f"{grid.resolution_identifier}+chou2001-lw",
    )
    if new_grid.n_bins != n_short + len(pieces):
        raise ConfigurationError("internal: Chou's longwave pieces do not match the grid")

    def column(position: int, fill: float) -> np.ndarray:
        return np.array([fill] * n_short + [p[position] for p in pieces])

    band_index = column(2, -1).astype(int)
    continuum = np.array(
        [0.0] * n_short
        + [CHOU_WATER_BANDS[p[2]].continuum_over_mu / _CHOU_DIFFUSIVITY for p in pieces]
    )
    return new_grid, ChouLongwave(
        water_band=band_index,
        water_term=column(7, -1).astype(int),
        tabulated=tuple(tabulated),
        water_k=column(3, 0.0),
        continuum_k=continuum,
        co2_k=column(4, 0.0),
        co2_kind=column(5, 0).astype(int),
        ozone_term=column(6, -1).astype(int),
    )


def _tabulated_coefficients(table: TabulatedWater, pressure: np.ndarray) -> np.ndarray:
    """A tabulated band's coefficients at each pressure (Pa), shaped (layers, terms).

    log k is interpolated in log p; held above 0.1 hPa (the Doppler limit) and
    extrapolated log-linearly below 1000 hPa, as ozone's band 5 is.
    """
    logp = np.log10(np.asarray(pressure, dtype=float) / 100.0)
    table_p = np.asarray(WATER_TABLE_LOG10_P_HPA)
    logk = np.log(np.array(table.k_per_g_cm2, dtype=np.float64))
    clipped = np.maximum(logp, table_p[0])
    out = np.empty((logp.size, logk.shape[1]))
    for n in range(logk.shape[1]):
        inside = np.interp(clipped, table_p, logk[:, n])
        slope = (logk[-1, n] - logk[-2, n]) / (table_p[-1] - table_p[-2])
        beyond = logk[-1, n] + slope * (clipped - table_p[-1])
        out[:, n] = np.exp(np.where(clipped > table_p[-1], beyond, inside))
    return out


def water_vapour_optical_depth(
    spec: ChouLongwave,
    *,
    vapour: np.ndarray,
    layer_mass: np.ndarray,
    pressure: np.ndarray,
    temperature: np.ndarray,
) -> np.ndarray:
    """Vertical absorption optical depth of water vapour, lines and continuum.

    ``vapour`` and ``layer_mass`` in kg m^-2, ``pressure`` in Pa, ``temperature``
    in K; shaped (layers, bins).
    """
    w = np.asarray(vapour, dtype=float)
    m = np.asarray(layer_mass, dtype=float)
    p = np.asarray(pressure, dtype=float)
    theta = np.asarray(temperature, dtype=float)
    if np.any(w < 0.0):
        raise ConfigurationError("a vapour path cannot be negative")
    path = w * _KG_M2_TO_G_CM2  # g cm^-2
    q = np.divide(w, m, out=np.zeros_like(w), where=m > 0.0)
    tau = np.zeros((w.size, spec.water_band.size))
    dt = theta - _REFERENCE_TEMPERATURE_K
    for index, band in enumerate(CHOU_WATER_BANDS):
        bins = spec.water_band == index
        if not np.any(bins):
            continue
        if band.name in spec.tabulated:
            table = WATER_TABLES[band.name]
            k = _tabulated_coefficients(table, p)  # (layers, terms)
            h = 1.0 + table.alpha * dt + table.beta * dt**2
            tau[:, bins] = (path * h)[:, None] * k[:, spec.water_term[bins]]
            continue
        scaled = (
            path
            * (p / _WATER_REFERENCE_PRESSURE) ** band.m
            * (1.0 + band.alpha * dt + band.beta * dt**2)
        )
        tau[:, bins] = scaled[:, None] * spec.water_k[None, bins]
    # Self-broadened continuum (eq. 4.19-4.21): the partial pressure q p / 0.622
    # over 1013 hPa, stronger in the cold.
    self_broadening = q * p / 0.622 / 101300.0 * np.exp(1800.0 * (1.0 / theta - 1.0 / 296.0))
    tau += (path * self_broadening)[:, None] * spec.continuum_k[None, :]
    return tau
