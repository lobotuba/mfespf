"""Ozone's 9.6-micron band, derived from Chou et al.'s (2001) transmission tables.

Checked against the source: the fitted band reproduces the line-by-line
transmittances it was fitted to, its terms fill 980-1100 cm^-1 by their Planck
weights, no ozone absorbs nothing, and above 0.1 hPa the absorption per molecule
is held at the Doppler limit (docs/protocols/ozone-infrared-band.md).
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.physics.ozone_infrared import (
    O3_IR_BAND_CM,
    O3_IR_WEIGHTS,
    ozone_infrared_coefficients,
    ozone_infrared_grid,
    ozone_infrared_optical_depth,
)
from mfespf.run.column_configurations import column_grid

pytestmark = pytest.mark.stage(6)

#: Chou et al. (2001) band-5 transmittance at 250 K (table ``o1``), as carried
#: by WRF module_ra_goddard.F at commit 7042598715a7: (p hPa, u cm-atm, T).
CHOU_TABLE = (
    (10.0, 0.015849, 0.9535622),
    (10.0, 0.063096, 0.89595246),
    (10.0, 0.251189, 0.78733349),
    (10.0, 1.0, 0.61410975),
    (100.0, 0.015849, 0.92423612),
    (100.0, 0.063096, 0.80518728),
    (100.0, 0.251189, 0.62041485),
    (100.0, 1.0, 0.41349036),
    (1000.0, 0.015849, 0.9137494),
    (1000.0, 0.063096, 0.75415814),
    (1000.0, 0.251189, 0.51872611),
    (1000.0, 1.0, 0.33376986),
)


@pytest.mark.req("SPC-001", "EQR-002", "VAL-001")
def test_the_band_reproduces_chou_s_tabulated_transmittance():
    """The fit's maximum error at 250 K over the whole table was 0.0103 (Z0)."""
    g = np.asarray(O3_IR_WEIGHTS)
    assert float(np.sum(g)) == pytest.approx(1.0, abs=1e-12)
    for p_hpa, u, expected in CHOU_TABLE:
        k = ozone_infrared_coefficients(np.array([p_hpa * 100.0]))[0]
        fitted = float(np.sum(g * np.exp(-1.66 * k * u)))
        assert fitted == pytest.approx(expected, abs=0.0104), (p_hpa, u)


@pytest.mark.req("SPC-001", "EQR-002")
def test_the_terms_fill_the_band_by_their_planck_weights():
    h, c, kb = 6.62607015e-34, 2.99792458e8, 1.380649e-23
    grid, terms = ozone_infrared_grid(column_grid(16))
    inside = terms >= 0
    assert sorted(terms[inside].tolist()) == list(range(6))
    # Nothing outside 980-1100 cm^-1 is ozone's, and the grid outside is untouched.
    low_m, high_m = 1.0e-2 / O3_IR_BAND_CM[1], 1.0e-2 / O3_IR_BAND_CM[0]
    assert np.all((grid.centers[inside] > low_m) & (grid.centers[inside] < high_m))
    original = np.asarray(column_grid(16).edges)
    kept = original[(original < low_m) | (original > high_m)]
    assert np.all(np.isin(kept, grid.edges))
    # Each term's sub-band holds its weight of the band's Planck flux at 250 K.
    nu = np.linspace(*O3_IR_BAND_CM, 40001)
    b = nu**3 / np.expm1(h * c * nu * 100.0 / (kb * 250.0))
    total = np.trapezoid(b, nu)
    edges = np.asarray(grid.edges)
    for i in np.flatnonzero(inside):
        lo, hi = 1.0e-2 / edges[i + 1], 1.0e-2 / edges[i]
        m = (nu >= lo) & (nu <= hi)
        share = np.trapezoid(b[m], nu[m]) / total if m.sum() > 1 else 0.0
        assert share == pytest.approx(O3_IR_WEIGHTS[terms[i]], abs=2e-3)


@pytest.mark.req("SPC-001", "EQR-002")
def test_no_ozone_absorbs_nothing_and_thin_air_holds_the_doppler_limit():
    _, terms = ozone_infrared_grid(column_grid(16))
    p = np.array([5.0e4, 10.0, 1.0])  # 500 hPa, 0.1 hPa, 0.01 hPa
    T = np.full(3, 250.0)
    none = ozone_infrared_optical_depth(terms, ozone_path=np.zeros(3), pressure=p, temperature=T)
    assert np.all(none == 0.0)
    k = ozone_infrared_coefficients(p)
    # Above 0.1 hPa the absorption per molecule stops depending on pressure.
    np.testing.assert_array_equal(k[2], k[1])
    # Narrower lines in thinner air: the wings, which hold almost all the band's
    # weight, absorb less per molecule at 0.1 hPa than at 500 hPa, and the line
    # peaks, the strongest term, absorb more.
    assert np.all(k[1, 1:5] < k[0, 1:5])
    assert k[1, 5] > k[0, 5]
    tau = ozone_infrared_optical_depth(
        terms, ozone_path=np.array([0.01, 0.01, 0.01]), pressure=p, temperature=T
    )
    assert np.all(tau[:, terms < 0] == 0.0)
    assert np.all(tau[:, terms >= 0] >= 0.0)
