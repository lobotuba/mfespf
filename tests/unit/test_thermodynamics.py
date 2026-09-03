"""Stage 3: correctness gates for the thermodynamic primitives."""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError
from mfespf.physics.thermodynamics import (
    EPSILON,
    KAPPA,
    dry_adiabatic_lapse_rate,
    dry_convective_adjustment,
    exner,
    moist_adiabatic_lapse_rate,
    potential_temperature,
    saturation_specific_humidity,
    saturation_vapour_pressure,
)

pytestmark = pytest.mark.stage(3)

P0 = value_of("standard_surface_pressure")


def _column(n: int = 8) -> tuple[np.ndarray, np.ndarray]:
    """Equal-mass layers, midpoint pressures from the surface upward."""
    interfaces = np.linspace(P0, 1.0e3, n + 1)
    pressure = 0.5 * (interfaces[:-1] + interfaces[1:])
    mass = -np.diff(interfaces) / value_of("earth_standard_gravity")
    return pressure, value_of("specific_heat_dry_air") * mass


@pytest.mark.req("ATM-002", "VAL-001")
def test_exner_and_potential_temperature_are_inverse():
    pressure, _ = _column()
    T = np.linspace(288.0, 220.0, pressure.size)
    theta = potential_temperature(T, pressure)
    assert np.allclose(theta * exner(pressure), T, rtol=1e-14)
    # At the reference pressure the two coincide by definition.
    assert exner(np.array([P0]))[0] == pytest.approx(1.0)


@pytest.mark.req("ATM-002")
def test_dry_adiabatic_lapse_rate_is_g_over_cp():
    expected = value_of("earth_standard_gravity") / value_of("specific_heat_dry_air")
    assert dry_adiabatic_lapse_rate() == pytest.approx(expected)
    assert dry_adiabatic_lapse_rate() == pytest.approx(9.76e-3, rel=1e-2)  # ~9.8 K/km


def _unstable_profile(pressure: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """A temperature profile built to be convectively unstable somewhere.

    Instability is constructed in potential-temperature space rather than in
    temperature space. Over a column spanning 1013 to 10 hPa the dry adiabat
    itself falls from 300 K to about 80 K, so almost any temperature profile a
    person would write down by hand is stable, and a test that perturbs
    temperature directly would exercise nothing.
    """
    theta = 300.0 + rng.normal(0.0, 15.0, size=pressure.size)
    return theta * exner(pressure)


@pytest.mark.req("ATM-002", "CON-009", "VAL-001")
def test_convective_adjustment_conserves_column_enthalpy_exactly():
    """The defining property: mixing moves energy, it does not create it."""
    pressure, C = _column(10)
    rng = np.random.default_rng(20260903)
    adjusted_any = False
    for _ in range(30):
        T = _unstable_profile(pressure, rng)
        adjusted = dry_convective_adjustment(T, C, pressure)
        adjusted_any = adjusted_any or not np.allclose(adjusted, T)
        assert float(np.sum(C * adjusted)) == pytest.approx(float(np.sum(C * T)), rel=1e-13)
    assert adjusted_any, "the profiles were all stable, so nothing was tested"


@pytest.mark.req("ATM-002", "VAL-001")
def test_convective_adjustment_leaves_a_stable_profile():
    """Potential temperature must not decrease with height afterwards."""
    pressure, C = _column(10)
    rng = np.random.default_rng(11)
    saw_instability = False
    for _ in range(30):
        T = _unstable_profile(pressure, rng)
        saw_instability = saw_instability or bool(
            np.any(np.diff(potential_temperature(T, pressure)) < 0.0)
        )
        theta = potential_temperature(dry_convective_adjustment(T, C, pressure), pressure)
        assert np.all(np.diff(theta) >= -1e-9)
    assert saw_instability, "no profile was unstable, so the criterion was never exercised"


@pytest.mark.req("ATM-002", "VAL-001")
def test_an_already_stable_profile_is_left_untouched():
    """The null case: adjustment must be inert where nothing is unstable."""
    pressure, C = _column(8)
    stable = 300.0 * exner(pressure)  # exactly a dry adiabat
    assert np.allclose(dry_convective_adjustment(stable, C, pressure), stable, rtol=1e-14)

    strongly_stable = 300.0 * exner(pressure) + np.linspace(0.0, 40.0, pressure.size)
    assert np.allclose(
        dry_convective_adjustment(strongly_stable, C, pressure), strongly_stable, rtol=1e-14
    )


@pytest.mark.req("ATM-002")
def test_a_wholly_unstable_column_mixes_to_a_single_adiabat():
    """The limiting case: potential temperature decreasing everywhere mixes to one value."""
    pressure, C = _column(6)
    theta = np.linspace(320.0, 260.0, 6)  # decreasing upward: unstable throughout
    T = theta * exner(pressure)

    adjusted = dry_convective_adjustment(T, C, pressure)
    mixed = potential_temperature(adjusted, pressure)

    assert np.allclose(mixed, mixed[0], rtol=1e-12)
    # Enthalpy conservation fixes the mixed value; it is not a free choice.
    assert mixed[0] == pytest.approx(float(np.sum(C * T) / np.sum(C * exner(pressure))), rel=1e-13)
    # The adjusted profile now cools with height, as a dry adiabat must.
    assert np.all(np.diff(adjusted) < 0.0)


@pytest.mark.req("ATM-002")
def test_adjustment_only_touches_the_unstable_group():
    """A stable upper column must not be disturbed by mixing below it."""
    pressure, C = _column(8)
    theta = 280.0 + np.linspace(0.0, 70.0, 8)  # stable everywhere
    theta[0], theta[1] = theta[1], theta[0]  # invert only the lowest pair
    T = theta * exner(pressure)

    adjusted = dry_convective_adjustment(T, C, pressure)

    assert np.allclose(adjusted[2:], T[2:], rtol=1e-14)
    assert not np.allclose(adjusted[:2], T[:2])
    # And the pair that mixed conserved its own enthalpy, not just the column's.
    assert float(np.sum(C[:2] * adjusted[:2])) == pytest.approx(
        float(np.sum(C[:2] * T[:2])), rel=1e-13
    )


@pytest.mark.req("WAT-001", "VAL-001")
def test_saturation_vapour_pressure_at_the_triple_point():
    triple = value_of("triple_point_temperature")
    assert saturation_vapour_pressure(np.array([triple]))[0] == pytest.approx(
        value_of("saturation_vapour_pressure_triple_point"), rel=1e-12
    )


@pytest.mark.req("WAT-001")
def test_saturation_vapour_pressure_matches_reference_values():
    """Against standard tables, within the constant-latent-heat approximation.

    The integrated Clausius-Clapeyron form with a fixed L is good to a percent
    or so across the terrestrial range; asserting more would be testing the
    approximation rather than the implementation.
    """
    for temperature, reference_pa in ((273.15, 611.2), (293.15, 2339.0), (303.15, 4246.0)):
        computed = float(saturation_vapour_pressure(np.array([temperature]))[0])
        assert computed == pytest.approx(reference_pa, rel=0.03)


@pytest.mark.req("WAT-001")
def test_saturation_specific_humidity_is_monotone_and_bounded():
    pressure = np.full(5, P0)
    temperatures = np.array([240.0, 260.0, 280.0, 290.0, 300.0])
    q_s = saturation_specific_humidity(temperatures, pressure)
    assert np.all(np.diff(q_s) > 0.0)
    assert np.all((q_s > 0.0) & (q_s < 1.0))
    # About 22 g/kg at 300 K and sea level.
    assert q_s[-1] == pytest.approx(0.022, rel=0.1)


@pytest.mark.req("WAT-001")
def test_specific_humidity_uses_the_full_form_not_the_thin_approximation():
    """The two differ by about 2 percent at 300 K, which is not negligible here."""
    T, p = np.array([300.0]), np.array([P0])
    full = float(saturation_specific_humidity(T, p)[0])
    thin = float(EPSILON * saturation_vapour_pressure(T)[0] / p[0])
    assert full > thin
    assert (full - thin) / full == pytest.approx(0.02, abs=0.01)


@pytest.mark.req("ATM-002", "WAT-002")
def test_moist_lapse_rate_is_below_dry_and_approaches_it_when_dry():
    """Latent release makes a saturated parcel cool more slowly as it rises."""
    warm = moist_adiabatic_lapse_rate(np.array([300.0]), np.array([P0]))[0]
    cold = moist_adiabatic_lapse_rate(np.array([220.0]), np.array([P0]))[0]
    dry = dry_adiabatic_lapse_rate()

    assert warm < cold < dry
    assert warm == pytest.approx(3.6e-3, rel=0.2)  # ~3.6 K/km in the warm tropics
    # In the cold limit the saturation humidity vanishes and the two coincide.
    assert cold == pytest.approx(dry, rel=0.02)


@pytest.mark.req("SYS-006", "MOD-002")
def test_unphysical_inputs_are_refused():
    with pytest.raises(ConfigurationError, match="positive"):
        exner(np.array([-1.0]))
    with pytest.raises(ConfigurationError, match="positive"):
        saturation_vapour_pressure(np.array([0.0]))
    with pytest.raises(ConfigurationError, match="same layers"):
        dry_convective_adjustment(np.zeros(3), np.ones(4), np.ones(3))


@pytest.mark.req("ATM-002")
def test_kappa_matches_the_declared_constants():
    assert KAPPA == pytest.approx(
        value_of("gas_constant_dry_air") / value_of("specific_heat_dry_air")
    )
    assert KAPPA == pytest.approx(0.2857, rel=1e-3)
    assert EPSILON == pytest.approx(0.622, rel=1e-3)
