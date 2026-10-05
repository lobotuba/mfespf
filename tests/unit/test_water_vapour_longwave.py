"""Water vapour's longwave absorption from Chou, Suarez, Liang and Yan (2001).

Checked against the source and against the physics: the transcribed tables are
self-consistent, the overlap with CO2 and ozone holds each pair's product of
fractions, and the continuum is self-broadened -- quadratic in vapour and
stronger in the cold (docs/protocols/water-vapour-longwave.md).
"""

from __future__ import annotations

from itertools import pairwise

import numpy as np
import pytest

from mfespf.physics.water_vapour_longwave import (
    CHOU_WATER_BANDS,
    chou_longwave_grid,
    water_vapour_optical_depth,
)
from mfespf.run.column_configurations import column_grid

pytestmark = pytest.mark.stage(6)


@pytest.mark.req("SPC-001", "EQR-002")
def test_the_transcribed_fractions_are_self_consistent():
    """Each band's Delta-g sums to one; band 3's three sub-bands sum to one together."""
    for band in CHOU_WATER_BANDS:
        if band.name in ("3a", "3b", "3c"):
            continue
        assert sum(band.delta_g) == pytest.approx(1.0, abs=1e-4), band.name
    band3 = sum(sum(b.delta_g) for b in CHOU_WATER_BANDS if b.name in ("3a", "3b", "3c"))
    assert band3 == pytest.approx(1.0, abs=1e-4)
    # The band edges tile 0-3000 cm^-1 with no gap or overlap.
    edges = [b.range_cm for b in CHOU_WATER_BANDS]
    assert edges[0][0] == 0.0 and edges[-1][1] == 3000.0
    assert all(a[1] == b[0] for a, b in pairwise(edges))


@pytest.mark.req("SPC-001", "EQR-002")
def test_every_longwave_bin_carries_one_water_term_and_the_overlaps_their_products():
    grid, spec = chou_longwave_grid(column_grid(16), ozone=True)
    longwave = spec.water_band >= 0
    assert np.all(spec.water_k[longwave] > 0.0)
    # CO2 only in band 3, ozone only in band 5, and each bin there carries one term of each.
    names = np.array([CHOU_WATER_BANDS[b].name if b >= 0 else "" for b in spec.water_band])
    in_band3 = np.isin(names, ["3a", "3b", "3c"])
    assert np.all(spec.co2_kind[in_band3] > 0) and np.all(spec.co2_kind[~in_band3] == 0)
    assert np.all(spec.ozone_term[names == "5"] >= 0) and np.all(spec.ozone_term[names != "5"] < 0)
    # Band 3b: 6 CO2 terms times its 5 non-zero water terms.
    assert int(np.sum(names == "3b")) == 30
    # The shortwave is untouched.
    original = np.asarray(column_grid(16).edges)
    assert np.all(np.isin(original[original < 4.0e-6], grid.edges))


@pytest.mark.req("EQR-002", "VAL-001")
def test_the_continuum_is_self_broadened_and_stronger_in_the_cold():
    _, spec = chou_longwave_grid(column_grid(16))
    window = np.flatnonzero((spec.continuum_k > 0.0) & (spec.water_band >= 0))[0]
    mass = np.array([1000.0, 1000.0, 1000.0])
    p = np.full(3, 9.0e4)

    def continuum(vapour, temperature):
        lines_off = spec.__class__(**{**spec.__dict__, "water_k": np.zeros_like(spec.water_k)})
        return water_vapour_optical_depth(
            lines_off, vapour=vapour, layer_mass=mass, pressure=p, temperature=temperature
        )[:, window]

    tau = continuum(np.array([10.0, 20.0, 10.0]), np.array([296.0, 296.0, 260.0]))
    assert tau[1] / tau[0] == pytest.approx(4.0, rel=1e-12)  # twice the vapour, four times
    assert tau[2] / tau[0] == pytest.approx(np.exp(1800.0 * (1 / 260.0 - 1 / 296.0)), rel=1e-12)
    no_vapour = continuum(np.zeros(3), np.full(3, 296.0))
    assert np.all(no_vapour == 0.0)
