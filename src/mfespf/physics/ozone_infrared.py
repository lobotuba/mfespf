"""Ozone's 9.6-micron band, derived from Chou et al.'s (2001) transmission tables.

Why this exists
---------------

``tropical-heating-diagnostic.md`` found heat missing from the tropical lower
stratosphere that the troposphere does not supply, and inferred ozone's absent
9.6-micron band as the likely cause. Ozone absorbs the upwelling infrared
there. Gettelman et al. (2004) put that absorption, with CO2's slight heating
at the tropopause, as the heating the CO2 and vapour cooling partly balance.

Where the numbers come from
---------------------------

Chou, Suarez, Liang and Yan (2001), NASA/TM-2001-104606 vol. 19, section 4.2
and Table 8: band 5, 980-1100 cm^-1, treated by pre-computed transmission
tables from line-by-line calculations with HITRAN 1996. The weak ozone
absorption of bands 4 and 6 is folded in (their eq. 4.17-4.18).

The tables are not printed in the report. They are carried, public domain, by
WRF's ``phys/module_ra_goddard.F`` (commit 7042598715a7; ``data/sources.yaml``)
as ``o1``, ``o2`` and ``o3``: a quadratic in temperature on a grid of log10
amount by log10 pressure.

What is derived here
--------------------

A path-transmittance table cannot be handed to a solver that takes per-layer
optical depths, so an exponential sum was fitted to it
(``docs/protocols/ozone-infrared-band.md``, attempt 2):

- six terms with weights g_n, shared by every pressure;
- coefficients k_n tabulated at the table's 21 pressures from 0.1 to 1000 hPa;
- a shared temperature scaling ``h = 1 + alpha (T - 250) + beta (T - 250)^2``.

A single power-law pressure scaling was tried first and failed (attempt 1). It
failed for the reason Chou et al. give for using tables: the band runs from
Doppler-broadened to pressure-broadened lines across these pressures.

**How well it fits.** The RMS transmittance error is 0.0015 at 250 K, with a
maximum of 0.010. At held-out pressures, with log k interpolated in log p, the
RMS is 0.0018.

**Outside the tabulated pressures:**

- above 0.1 hPa each k is held at its 0.1 hPa value: the Doppler limit, where a
  line's width stops depending on pressure;
- below 1000 hPa it is extrapolated log-linearly: the Lorentz regime.
"""

from __future__ import annotations

import numpy as np

from mfespf.contracts.spectral import SpectralCoordinate, SpectralGrid
from mfespf.core.errors import ConfigurationError

__all__ = [
    "O3_IR_ALPHA",
    "O3_IR_BAND_CM",
    "O3_IR_BETA",
    "O3_IR_K_PER_CM_ATM",
    "O3_IR_LOG10_P_HPA",
    "O3_IR_WEIGHTS",
    "ozone_infrared_coefficients",
    "ozone_infrared_grid",
    "ozone_infrared_optical_depth",
]

_H = 6.62607015e-34
_C = 2.99792458e8
_KB = 1.380649e-23
_REFERENCE_TEMPERATURE_K = 250.0

#: The band, cm^-1: Chou et al.'s band 5.
O3_IR_BAND_CM: tuple[float, float] = (980.0, 1100.0)

#: Weights of the six terms, summing to one. The first is near zero: the fit
#: found five terms enough, and the sixth is kept rather than the form changed.
O3_IR_WEIGHTS: tuple[float, ...] = (
    4.397700080985617e-06,
    0.5494647890550012,
    0.2018765862535625,
    0.1720781660770513,
    0.05804458988817751,
    0.018531471026126612,
)

#: log10 of the tabulated pressures in hPa: Chou's Table 8 grid from 0.1 hPa.
O3_IR_LOG10_P_HPA: tuple[float, ...] = (
    -1.0,
    -0.8,
    -0.6,
    -0.4,
    -0.2,
    0.0,
    0.2,
    0.4,
    0.6,
    0.8,
    1.0,
    1.2,
    1.4,
    1.6,
    1.8,
    2.0,
    2.2,
    2.4,
    2.6,
    2.8,
    3.0,
)

#: Vertical absorption coefficient of each term at each pressure, (cm-atm)^-1.
#: The solver applies the diffusivity 1.66 the fit assumed.
O3_IR_K_PER_CM_ATM: tuple[tuple[float, ...], ...] = (
    (
        0.0035433031756328935,
        0.06785609280377533,
        0.0681573629079284,
        0.06833553965477113,
        5.378410135613613,
        124.88079322246259,
    ),
    (
        0.003601097775793041,
        0.07147772337795436,
        0.07185714670961947,
        0.07204819507028101,
        5.519523169703336,
        124.26811999215231,
    ),
    (
        0.003918096612981572,
        0.07671221164781684,
        0.07724354522713268,
        0.07747963701164517,
        5.745622971638955,
        123.32129022605608,
    ),
    (
        0.003565728823598375,
        0.08374510764208305,
        0.08535312847652507,
        0.08594588400970629,
        6.104046475252462,
        121.90859761599647,
    ),
    (
        0.003282597463947263,
        0.0545535139325403,
        0.05679490073974202,
        0.3043841161088577,
        6.343452340583493,
        122.40441310115158,
    ),
    (
        0.0034367718124955946,
        0.0489053909902308,
        0.0508437614498268,
        0.46082614772047487,
        6.7129145391918605,
        123.30496896041181,
    ),
    (
        0.003963384792557368,
        0.05073372093079105,
        0.052573685489291304,
        0.6329173830043495,
        7.317772616636669,
        124.38362097151999,
    ),
    (
        0.004845208956843172,
        0.05919682009212856,
        0.060966003430751425,
        0.8411036070094717,
        8.25052524282132,
        125.51296175417349,
    ),
    (
        0.006198964973269516,
        0.07475326913971837,
        0.07644290503756886,
        1.1064971359385471,
        9.60855102046974,
        126.19520804796754,
    ),
    (
        0.007859677268168358,
        0.09809317777355733,
        0.09970946244035168,
        1.459168396035166,
        11.448721996590244,
        125.48079002323459,
    ),
    (
        0.010514667171117376,
        0.12974276447073416,
        0.13140929299315396,
        1.9443034422566643,
        13.689594664883181,
        122.29354394530776,
    ),
    (
        0.012735904212168303,
        0.16960579610923682,
        0.17218729160877766,
        2.628429227571586,
        15.945496686135295,
        116.46837102065606,
    ),
    (
        0.012337979633202401,
        0.14593245385731093,
        0.4730506281751013,
        3.4152346043986612,
        18.342385285592496,
        105.08442642018358,
    ),
    (
        0.014602399025716103,
        0.15319955119221412,
        0.7837120759556133,
        4.29579271231307,
        20.360380060570332,
        90.26833607318109,
    ),
    (
        0.019238862142113743,
        0.17275506126049883,
        1.1460117599085953,
        5.280879918114512,
        21.065454579402225,
        77.11809870065419,
    ),
    (
        0.028671324116516208,
        0.1979930684543358,
        1.5933094191844364,
        6.221162171251701,
        20.472733921027256,
        66.9285566205477,
    ),
    (
        0.04321418379989551,
        0.225524837120282,
        2.146406240434177,
        6.864370584387775,
        19.550583575568357,
        58.18224048105172,
    ),
    (
        0.05765003813374431,
        0.2521596291043291,
        2.822681203692702,
        6.945754377942939,
        19.577557841285845,
        48.89405081295847,
    ),
    (
        0.06101025722836454,
        0.2751575194147227,
        3.6698780008463308,
        6.3736591481723925,
        21.256879147174306,
        37.71566179220734,
    ),
    (
        0.043694549420828344,
        0.2938315872211517,
        5.0603961516058265,
        5.066423136297026,
        22.093917283166935,
        31.622290830855952,
    ),
    (
        0.02267770414243769,
        0.30942824086954485,
        5.3492478745690635,
        5.352790261572917,
        17.6698426864774,
        43.12268050983625,
    ),
)
O3_IR_ALPHA = 0.0046227586502563
O3_IR_BETA = 4.0957047469971285e-06


def _planck_cumulative(low: float, high: float, n: int = 20001) -> tuple[np.ndarray, np.ndarray]:
    """Wavenumber mesh (cm^-1) and cumulative Planck flux at 250 K, normalised."""
    nu = np.linspace(low, high, n)
    x = nu * 100.0
    density = x**3 / np.expm1(_H * _C * x / (_KB * _REFERENCE_TEMPERATURE_K))
    cumulative = np.concatenate(
        ([0.0], np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(nu)))
    )
    return nu, cumulative / cumulative[-1]


def ozone_infrared_grid(grid: SpectralGrid) -> tuple[SpectralGrid, np.ndarray]:
    """Rebuild a wavelength grid's 9.09-10.2 micron region around the band's terms.

    Returns the grid and each bin's term index (-1 outside the band). Each term
    occupies the share of 980-1100 cm^-1 that holds its weight of the Planck
    flux at 250 K. Edges outside the band are kept exactly.
    """
    if grid.coordinate is not SpectralCoordinate.WAVELENGTH:
        raise ConfigurationError("the ozone band is built on a wavelength grid")
    edges = np.asarray(grid.edges, dtype=float)
    low_cm, high_cm = O3_IR_BAND_CM
    low_m, high_m = 1.0e-2 / high_cm, 1.0e-2 / low_cm
    if edges[0] > low_m or edges[-1] < high_m:
        raise ConfigurationError("the grid must span 9.09 to 10.2 microns to hold the ozone band")
    weights = np.asarray(O3_IR_WEIGHTS, dtype=float)
    nu, cumulative = _planck_cumulative(low_cm, high_cm)
    shares = np.cumsum(weights)[:-1]
    bounds = np.concatenate(([low_cm], np.interp(shares, cumulative, nu), [high_cm]))
    # Term n spans [bounds[n], bounds[n+1]] in wavenumber; in wavelength the
    # order reverses, which the sort handles.
    pieces = sorted((1.0e-2 / bounds[n + 1], 1.0e-2 / bounds[n], n) for n in range(weights.size))
    kept_low = [float(e) for e in edges if e < low_m and not np.isclose(e, low_m)]
    kept_high = [float(e) for e in edges if e > high_m and not np.isclose(e, high_m)]
    new_edges = [*kept_low, *(p[0] for p in pieces), pieces[-1][1], *kept_high]
    terms = [-1] * len(kept_low) + [p[2] for p in pieces] + [-1] * len(kept_high)
    new_grid = SpectralGrid(
        grid_id=f"{grid.grid_id}+o3_9.6um_chou2001",
        coordinate=grid.coordinate,
        coordinate_units=grid.coordinate_units,
        edges=np.asarray(new_edges),
        integration_convention=grid.integration_convention,
        density_units=grid.density_units,
        resolution_identifier=f"{grid.resolution_identifier}+chou2001-o3ir",
    )
    term_index = np.asarray(terms, dtype=int)
    if term_index.size != new_grid.n_bins:
        raise ConfigurationError("internal: the ozone terms do not match the rebuilt grid")
    return new_grid, term_index


def ozone_infrared_coefficients(pressure: np.ndarray) -> np.ndarray:
    """Each term's coefficient at each pressure (Pa), shaped (layers, 6), (cm-atm)^-1."""
    logp = np.log10(np.asarray(pressure, dtype=float) / 100.0)
    table_p = np.asarray(O3_IR_LOG10_P_HPA)
    logk = np.log(np.asarray(O3_IR_K_PER_CM_ATM, dtype=float))
    # Held above 0.1 hPa (the Doppler limit); log-linear below 1000 hPa.
    clipped = np.maximum(logp, table_p[0])
    out = np.empty((logp.size, logk.shape[1]))
    for n in range(logk.shape[1]):
        inside = np.interp(clipped, table_p, logk[:, n])
        slope = (logk[-1, n] - logk[-2, n]) / (table_p[-1] - table_p[-2])
        beyond = logk[-1, n] + slope * (clipped - table_p[-1])
        out[:, n] = np.exp(np.where(clipped > table_p[-1], beyond, inside))
    return out


def ozone_infrared_optical_depth(
    term_index: np.ndarray,
    *,
    ozone_path: np.ndarray,
    pressure: np.ndarray,
    temperature: np.ndarray,
) -> np.ndarray:
    """Vertical absorption optical depth of ozone in the band, shaped (layers, bins)."""
    u = np.asarray(ozone_path, dtype=float)
    dt = np.asarray(temperature, dtype=float) - _REFERENCE_TEMPERATURE_K
    h = 1.0 + O3_IR_ALPHA * dt + O3_IR_BETA * dt**2
    k = ozone_infrared_coefficients(pressure)
    idx = np.asarray(term_index, dtype=int)
    tau = np.zeros((u.size, idx.size))
    inside = idx >= 0
    tau[:, inside] = (u * h)[:, None] * k[:, idx[inside]]
    return tau
