"""The 15-micron CO2 band of Chou, Suarez, Liang and Yan (2001).

Checked against the source and against the physics it exists for: the terms
carry the published coefficients, the centre's Planck-weighted fractions are
the Planck share of its wavenumbers (a consistency the table does not state
but must satisfy), and the same CO2 absorbs less per molecule in thin air --
the pressure broadening the calibrated band left out.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.physics.co2_band import (
    CO2_CENTRE,
    CO2_WINGS,
    co2_band_grid,
    co2_layer_amount,
    co2_optical_depth,
)
from mfespf.run.column_configurations import column_grid

pytestmark = pytest.mark.stage(6)


def band_bins():
    grid, k, kind = co2_band_grid(column_grid(16))
    return grid, k, kind


@pytest.mark.req("SPC-001", "EQR-002")
def test_the_band_terms_carry_the_source_coefficients():
    grid, k, kind = band_bins()
    for code, spec in ((1, CO2_CENTRE), (2, CO2_WINGS)):
        expected = spec.k1_over_mu * 8.0 ** np.arange(6) / 1.66
        present = np.unique(k[kind == code])
        np.testing.assert_allclose(np.sort(present), np.sort(expected), rtol=1e-12)
    # The centre is one sub-band of six terms; the wings are two of six each.
    assert int(np.sum(kind == 1)) == 6
    assert int(np.sum(kind == 2)) == 12
    # Nothing outside 540-800 cm^-1 is CO2's, and the grid outside is untouched.
    centres_cm = 1.0e-2 / grid.centers
    assert np.all(k[(centres_cm < 540.0) | (centres_cm > 800.0)] == 0.0)
    original = np.asarray(column_grid(16).edges)
    kept = original[(original < 1.0e-2 / 800.0) | (original > 1.0e-2 / 540.0)]
    assert np.all(np.isin(kept, grid.edges))


@pytest.mark.req("SPC-001", "VAL-001")
def test_the_centre_fractions_are_the_planck_share():
    """0.389 of the band's Planck flux at 250 K lies in 620-720 cm^-1."""
    h, c, kb = 6.62607015e-34, 2.99792458e8, 1.380649e-23
    nu = np.linspace(540.0, 800.0, 20001)
    b = (nu * 100.0) ** 3 / np.expm1(h * c * nu * 100.0 / (kb * 250.0))
    centre = (nu >= 620.0) & (nu <= 720.0)
    share = np.trapezoid(b[centre], nu[centre]) / np.trapezoid(b, nu)
    assert float(np.sum(CO2_CENTRE.fractions)) == pytest.approx(share, abs=0.002)
    total = float(np.sum(CO2_CENTRE.fractions) + np.sum(CO2_WINGS.fractions))
    assert total == pytest.approx(1.0, abs=1e-9)


@pytest.mark.req("EQR-002", "VAL-001")
def test_thin_air_absorbs_less_per_molecule():
    """The same CO2 at 1 hPa absorbs a fraction of what it does at 300 hPa."""
    _, k, kind = band_bins()
    mass = np.array([10.0, 10.0])
    tau = co2_optical_depth(
        k,
        kind,
        co2_ppmv=350.0,
        layer_mass=mass,
        pressure=np.array([3.0e4, 100.0]),
        temperature=np.array([250.0, 250.0]),
    )
    wings = kind == 2
    centre = kind == 1
    # Wings scale as (p / 300 hPa)^0.5, the centre as (p / 30 hPa)^0.85.
    ratio_wings = tau[1, wings] / tau[0, wings]
    ratio_centre = tau[1, centre] / tau[0, centre]
    np.testing.assert_allclose(ratio_wings, (100.0 / 3.0e4) ** 0.5, rtol=1e-12)
    np.testing.assert_allclose(ratio_centre, (100.0 / 3.0e4) ** 0.85, rtol=1e-12)
    # And no CO2, no absorption.
    assert np.all(
        co2_optical_depth(
            k,
            kind,
            co2_ppmv=0.0,
            layer_mass=mass,
            pressure=np.array([3.0e4, 100.0]),
            temperature=np.array([250.0, 250.0]),
        )
        == 0.0
    )
    # A standard column is about 8 km of air at STP, so 350 ppmv of CO2 is
    # about 280 cm of pure CO2 at STP: the amount is the mixing ratio times the air.
    column = co2_layer_amount(350.0, np.array([101325.0 / 9.80665]))
    assert float(column[0]) == pytest.approx(280.0, rel=0.02)
