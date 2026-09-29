"""Water vapour absorbs sunlight, with the coefficients of Chou and Suarez (1999).

These check the implementation against the source rather than against the
model: the sub-bands must carry the source's flux weights, a uniform path must
transmit the source's own k-distribution sum, the scaling must be exactly one
at the source's reference state, and the rebuilt grid must change nothing at
all when there is no vapour to absorb.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from mfespf.physics.greenhouse import vapour_absorption_optical_depth
from mfespf.physics.radiative_transfer import solve_column_radiation
from mfespf.physics.solar import ConstantSpectrumSolarSource
from mfespf.physics.solar_absorption import (
    SOLAR_VAPOUR_BANDS,
    SOLAR_VAPOUR_REFERENCE_PRESSURE_PA,
    SOLAR_VAPOUR_REFERENCE_TEMPERATURE_K,
    scaled_vapour_path,
    solar_vapour_bands,
    solar_vapour_optical_depth,
)
from mfespf.run.column_configurations import column_grid, reference_column

pytestmark = pytest.mark.stage(3)

GRID = column_grid(16)
REBUILT, K = solar_vapour_bands(GRID)


def _incident(grid) -> np.ndarray:
    """Each band's share of the sunlight, as the configured solar source sets it."""
    return np.asarray(ConstantSpectrumSolarSource(grid=grid).spectral_shape(), dtype=float)


@pytest.mark.req("SPC-004", "EQR-002")
def test_the_sub_bands_carry_the_source_weights():
    """Each k-term's sub-band receives its term's share of its band's sunlight."""
    incident = _incident(REBUILT)
    for lower, upper, ks, weights in SOLAR_VAPOUR_BANDS[1:3]:  # 0.70-1.22, 1.22-2.27
        inside = (REBUILT.centers > lower) & (REBUILT.centers < upper)
        shares = incident[inside] / incident[inside].sum()
        expected = np.array([h for h in weights if h > 0.0])
        assert shares == pytest.approx(expected / expected.sum(), rel=2e-3)
        assert K[inside] == pytest.approx(
            [0.1 * k for k, h in zip(ks, weights, strict=True) if h > 0.0]
        )
    # The thermal bands are the original grid's, edge for edge.
    top = GRID.edges[GRID.edges >= 2.27e-6][0]
    assert np.array_equal(REBUILT.edges[REBUILT.edges >= top], GRID.edges[GRID.edges >= top])
    assert np.all(K[REBUILT.centers > top] == 0.0)


@pytest.mark.req("EQR-002", "VAL-001")
def test_a_uniform_path_transmits_the_source_k_distribution():
    """The source's equation 3.4, evaluated directly, against the sub-band sum."""
    lower, upper, ks, weights = SOLAR_VAPOUR_BANDS[1]
    inside = (REBUILT.centers > lower) & (REBUILT.centers < upper)
    incident = _incident(REBUILT)
    for path in (0.1, 1.0, 10.0, 50.0):  # kg m^-2, at the reference state
        tau = solar_vapour_optical_depth(
            K,
            np.array([path]),
            np.array([SOLAR_VAPOUR_REFERENCE_PRESSURE_PA]),
            np.array([SOLAR_VAPOUR_REFERENCE_TEMPERATURE_K]),
        )[0]
        transmitted = np.sum(incident[inside] * np.exp(-tau[inside])) / incident[inside].sum()
        w_g_cm2 = path * 0.1  # kg m^-2 to g cm^-2
        h = np.array(weights)
        source = np.sum(h * np.exp(-np.array(ks) * w_g_cm2)) / h.sum()
        assert transmitted == pytest.approx(source, rel=3e-3)


@pytest.mark.req("EQR-002")
def test_the_scaling_is_one_at_the_reference_state():
    path = np.array([3.0, 3.0, 3.0])
    at_ref = scaled_vapour_path(path, np.full(3, 3.0e4), np.full(3, 240.0))
    assert at_ref == pytest.approx(path, rel=1e-15)
    # Lower pressure narrows the line wings and weakens the absorption.
    thin = scaled_vapour_path(path, np.full(3, 1.5e4), np.full(3, 240.0))
    assert thin == pytest.approx(path * 0.5**0.8)
    # And warmer vapour absorbs slightly more, by the source's coefficient.
    warm = scaled_vapour_path(path, np.full(3, 3.0e4), np.full(3, 290.0))
    assert warm == pytest.approx(path * (1.0 + 0.00135 * 50.0))


@pytest.mark.req("SPC-005", "VAL-001")
def test_a_dry_column_on_the_rebuilt_grid_changes_nothing():
    """The null test: the regridding alone moves no energy anywhere.

    With the vapour radiatively inert, the only difference between the two
    configurations is where the solar band edges lie. The albedo is flat below
    4 microns and nothing else in the solar region absorbs, so the runs must
    agree to round-off.
    """
    base = dict(
        grid=GRID,
        areal_heat_capacity=1.0e6,
        master_dt=3600.0,
        initial_water_store=1.0e7,
        radiative_vapour=False,
    )
    plain = reference_column(**base).run(200)
    rebuilt = reference_column(solar_vapour_absorption=True, **base).run(200)
    a = float(plain.steps[-1].outputs["surface_temperature"].values)
    b = float(rebuilt.steps[-1].outputs["surface_temperature"].values)
    assert b == pytest.approx(a, abs=1e-6)
    assert float(np.sum(rebuilt.steps[-1].outputs["atmosphere_absorbed_shortwave"].values)) == (
        pytest.approx(0.0, abs=1e-9)
    )


MLS = Path(__file__).parent.parent / "data" / "afgl1986_midlatitude_summer.csv"


def _absorbed_fraction(slant: float) -> dict[str, float]:
    """Fraction of each near-infrared band's sunlight the MLS atmosphere absorbs.

    The geometry of Chou and Suarez's Table 8: sun at 60 degrees, surface albedo
    0.2, no scattering. Thermal emission is removed by differencing a lit and a
    dark solve, which is exact because nothing scatters.
    """
    rows = list(csv.DictReader(MLS.open()))
    p = np.array([float(r["p"]) for r in rows]) * 100.0
    t = np.array([float(r["t"]) for r in rows])
    q = np.array([float(r["H2O"]) for r in rows]) * 1.0e-6 * 18.015 / 28.964
    mass = (p[:-1] - p[1:]) / 9.80665
    layer_t = 0.5 * (t[:-1] + t[1:])
    layer_p = np.sqrt(p[:-1] * p[1:])
    vapour = 0.5 * (q[:-1] + q[1:]) * mass
    tau = vapour_absorption_optical_depth(REBUILT, vapour) + solar_vapour_optical_depth(
        K, vapour, layer_p, layer_t
    )
    shortwave = REBUILT.centers < 4.0e-6
    shares = _incident(REBUILT)

    def solve(total: float):
        incident = np.where(shortwave, total * shares / shares[shortwave].sum(), 0.0)
        return solve_column_radiation(
            grid=REBUILT,
            incident_spectral_irradiance=incident / REBUILT.widths,
            layer_temperature=layer_t,
            layer_mass=mass,
            absorption_coefficient=np.zeros(REBUILT.n_bins),
            surface_temperature=float(t[0]),
            surface_albedo=np.where(shortwave, 0.2, 0.0),
            absorption_optical_depth=tau,
            diffusivity_factor=slant,
        )

    lit, dark = solve(1365.0 * 0.5), solve(0.0)
    top = (lit.downward[-1] - lit.upward[-1]) - (dark.downward[-1] - dark.upward[-1])
    bottom = (lit.downward[0] - lit.upward[0]) - (dark.downward[0] - dark.upward[0])
    arriving = lit.downward[-1] - dark.downward[-1]
    out = {}
    for lower, upper, _k, _h in SOLAR_VAPOUR_BANDS[1:3]:
        inside = (REBUILT.centers > lower) & (REBUILT.centers < upper)
        out[f"{lower * 1e6:.2f}"] = float(np.sum((top - bottom)[inside]) / np.sum(arriving[inside]))
    return out


@pytest.mark.req("EQR-002", "VAL-001", "SPC-005")
def test_the_mls_atmosphere_absorbs_what_line_by_line_says():
    """Chou and Suarez's Table 8, like for like: same atmosphere, same geometry.

    Their line-by-line atmospheric absorption, as a fraction of each band's
    sunlight (their Table 3 shares of a 682.5 W m^-2 beam), is 0.224 in
    0.70-1.22 microns and 0.447 in 1.22-2.27, *including* oxygen and carbon
    dioxide, which this model does not have (their Table 7 puts them at about
    4.3 and 3.3 W m^-2). So the vapour-only model must come in below the
    line-by-line figure, and close to it. Measured at the 60-degree slant path:
    0.221 and 0.410. The lower bound, 0.9 of the line-by-line value, was set
    after that measurement as a regression guard, not derived.

    At the solver's diffusivity factor, 1.66, the fractions are 0.207 and
    0.399: the diffuse treatment of sunlight under-absorbs by a few percent for
    a 60-degree sun, which is the size of that declared approximation.
    """
    line_by_line = {"0.70": 48.9 / (0.32055 * 682.5), "1.22": 50.5 / (0.16536 * 682.5)}
    model = _absorbed_fraction(slant=2.0)
    for band, reference in line_by_line.items():
        assert model[band] <= reference, (band, model[band], reference)
        assert model[band] >= 0.9 * reference, (band, model[band], reference)
    diffuse = _absorbed_fraction(slant=1.66)
    for band in line_by_line:
        assert 0.9 * model[band] < diffuse[band] < model[band]
