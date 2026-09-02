"""Stage 1: limiting cases and analytic identities for the radiative primitives.

These are the charter 15.7 correctness gates for the radiation equations:
dimensional consistency, limiting behaviour, and agreement with an independently
derived analytic result.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.contracts.spectral import (
    IntegrationConvention,
    SpectralCoordinate,
    SpectralGrid,
)
from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError
from mfespf.physics.radiation import (
    band_integrated_exitance,
    equilibrium_temperature,
    greybody_exitance,
    normalized_band_weights,
    planck_spectral_exitance,
    relaxation_timescale,
)

pytestmark = pytest.mark.stage(1)

SIGMA = value_of("stefan_boltzmann")


def _wide_grid(n=400, lo=1.0e-8, hi=1.0e-2) -> SpectralGrid:
    return SpectralGrid(
        grid_id="wide",
        coordinate=SpectralCoordinate.WAVELENGTH,
        coordinate_units="m",
        edges=np.logspace(np.log10(lo), np.log10(hi), n + 1),
        integration_convention=IntegrationConvention.SPECTRAL_DENSITY,
    )


@pytest.mark.req("RAD-002", "VAL-001")
def test_planck_integrates_to_stefan_boltzmann():
    """The Planck integral must reproduce sigma T^4 without being told to."""
    for temperature in (200.0, 288.0, 5772.0):
        grid = _wide_grid()
        total = float(np.sum(band_integrated_exitance(grid, temperature)))
        assert total == pytest.approx(SIGMA * temperature**4, rel=1e-6)


@pytest.mark.req("RAD-003")
def test_planck_wien_displacement():
    """Peak wavelength times temperature is Wien's constant, 2.898e-3 m K."""
    lam = np.logspace(-8, -3, 200_001)
    for temperature in (300.0, 3000.0):
        peak = lam[int(np.argmax(planck_spectral_exitance(lam, temperature)))]
        assert peak * temperature == pytest.approx(2.897771955e-3, rel=1e-3)


@pytest.mark.req("VAL-001")
def test_planck_rayleigh_jeans_limit():
    """At long wavelength the exitance approaches 2 pi c k T / lambda^4.

    The leading correction is x/2 with x = hc/(lambda k T), so the tolerance is
    set by the wavelength: at 0.1 m and 300 K, x is 4.8e-4 and the limit holds to
    about 2.4e-4 relative. Asserting anything tighter would be testing the
    correction term, not the limit.
    """
    lam = np.array([1.0e-1])
    temperature = 300.0
    expected = (
        2.0 * np.pi * value_of("speed_of_light") * value_of("boltzmann") * temperature / lam**4
    )
    assert planck_spectral_exitance(lam, temperature)[0] == pytest.approx(expected[0], rel=1e-3)


@pytest.mark.req("VAL-001")
def test_greybody_zero_emissivity_is_null():
    assert greybody_exitance(300.0, 0.0) == 0.0
    assert greybody_exitance(0.0, 1.0) == 0.0


@pytest.mark.req("MOD-002")
def test_greybody_refuses_an_emissivity_outside_its_validity_domain():
    with pytest.raises(ConfigurationError, match="validity domain"):
        greybody_exitance(300.0, 1.5)


@pytest.mark.req("SPC-005")
def test_band_weights_sum_to_the_captured_fraction():
    """A truncated grid reports how much of the spectrum it misses."""
    full = normalized_band_weights(_wide_grid(), 5772.0)
    assert float(np.sum(full)) == pytest.approx(1.0, rel=1e-6)

    narrow = SpectralGrid(
        grid_id="visible_only",
        coordinate=SpectralCoordinate.WAVELENGTH,
        coordinate_units="m",
        edges=np.array([4.0e-7, 7.0e-7]),
    )
    captured = float(np.sum(normalized_band_weights(narrow, 5772.0)))
    assert 0.2 < captured < 0.5  # the visible band holds roughly a third


@pytest.mark.req("ORB-002", "VAL-001")
def test_equilibrium_temperature_matches_the_closed_form():
    S, albedo, emissivity = 1361.0, 0.3, 1.0
    T = equilibrium_temperature(S, albedo, emissivity)
    # Independently: absorbed = emitted at equilibrium.
    assert greybody_exitance(T, emissivity) == pytest.approx(0.25 * (1 - albedo) * S)
    assert T == pytest.approx(254.6, abs=0.5)


@pytest.mark.req("TIM-008")
def test_relaxation_timescale_is_the_linearised_slab_response():
    C = 2.1e8
    T = 254.6
    tau = relaxation_timescale(C, T)
    assert tau == pytest.approx(C / (4.0 * SIGMA * T**3))
    # A 50 m ocean mixed layer relaxes on a timescale of order a year, not a day.
    assert 1.0e7 < tau < 1.0e8
