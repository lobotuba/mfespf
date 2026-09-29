"""Water vapour absorbs sunlight, with the coefficients of Chou and Suarez (1999).

These check the implementation against the source rather than against the
model: the sub-bands must carry the source's flux weights, a uniform path must
transmit the source's own k-distribution sum, the scaling must be exactly one
at the source's reference state, and the rebuilt grid must change nothing at
all when there is no vapour to absorb.
"""

from __future__ import annotations

import numpy as np
import pytest

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
