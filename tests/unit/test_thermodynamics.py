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
    hydrostatic_heights,
    lapse_rate_adjustment_with_groups,
    lapse_rate_offsets,
    moist_adiabatic_lapse_rate,
    potential_temperature,
    saturation_specific_humidity,
    saturation_vapour_pressure,
    surface_coupled_adjustment,
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


# --- the surface as the bottom of the column ---------------------------------


def _surface_mixed(temperature, *, surface_temperature, surface_heat_capacity=1.0e7, n=8):
    pressure, capacity = _column(n)
    return surface_coupled_adjustment(
        np.asarray(temperature, dtype=float),
        capacity,
        pressure,
        surface_temperature=surface_temperature,
        surface_heat_capacity=surface_heat_capacity,
        surface_pressure=P0,
    )


@pytest.mark.req("ATM-002", "CON-009", "VAL-001")
def test_the_surface_mixing_conserves_the_enthalpy_of_the_surface_and_the_air():
    """The same conservation argument, now spanning a module boundary.

    Prepending the surface does not weaken the guarantee that made the
    layer-only adjustment safe: the mixed value is still defined as the one that
    leaves the total unchanged. It matters more here, not less, because this
    mixing moves energy between two modules rather than within one, and an
    exchange that did not balance would be a leak the column budget has to see.
    """
    pressure, capacity = _column()
    T = np.linspace(250.0, 220.0, pressure.size)
    Cs, Ts = 1.0e7, 320.0

    adjusted_surface, adjusted, groups, surface_mixed = _surface_mixed(
        T, surface_temperature=Ts, surface_heat_capacity=Cs
    )

    assert surface_mixed and groups
    before = Cs * Ts + float(np.sum(capacity * T))
    after = Cs * adjusted_surface + float(np.sum(capacity * adjusted))
    assert after == pytest.approx(before, rel=1e-14)
    assert adjusted_surface < Ts  # the surface paid for the warming above it


@pytest.mark.req("ATM-002", "CON-009", "CLM-003")
def test_an_inversion_between_layers_cannot_hide_a_buoyant_surface():
    """The failure this scheme exists to remove.

    The layer-only adjustment compares each layer with the one above it, so a
    profile whose potential temperature increases upward is stable by
    construction and it does nothing -- however hot the ground underneath is.
    That is the state a column with trapped radiatively active vapour settles
    into, and it is self-sustaining: the adjustment is inactive exactly where
    the water is, so nothing lifts the water, so the inversion stays.

    Launching the mixing from the surface breaks it, because the surface is
    buoyant against the whole profile rather than against its neighbour only.
    """
    pressure, capacity = _column()
    inverted = np.array([240.0, 286.0, 283.0, 280.0, 277.0, 274.0, 271.0, 268.0])

    # The layer-only scheme sees nothing to do. Compared to a tolerance rather
    # than exactly, because the adjustment returns theta * exner and that round
    # trip is not bit-exact; the claim is that it does not move the profile.
    assert dry_convective_adjustment(inverted, capacity, pressure) == pytest.approx(
        inverted, rel=1e-14
    )

    _, adjusted, groups, surface_mixed = _surface_mixed(inverted, surface_temperature=300.0)

    assert surface_mixed
    assert groups[0][0] == 0  # and it reaches the layer holding the water
    assert adjusted[0] > inverted[0] + 20.0
    theta = potential_temperature(adjusted, pressure)
    assert np.all(np.diff(theta) >= -1e-9)  # and leaves a stable profile behind


@pytest.mark.req("ATM-002", "VAL-001", "DAT-007")
def test_a_surface_colder_than_the_air_is_left_alone_and_bad_geometry_refused():
    """A stable surface must not be mixed, and the column must have one bottom."""
    pressure, capacity = _column()
    stable = np.linspace(250.0, 220.0, pressure.size)

    surface, adjusted, groups, surface_mixed = _surface_mixed(stable, surface_temperature=200.0)
    assert not surface_mixed and groups == ()
    assert surface == 200.0
    assert adjusted == pytest.approx(stable, rel=1e-14)

    with pytest.raises(ConfigurationError, match="bottom of the column"):
        surface_coupled_adjustment(
            stable,
            capacity,
            pressure,
            surface_temperature=300.0,
            surface_heat_capacity=1.0e7,
            surface_pressure=1.0e3,  # above the lowest layer, so not a surface
        )
    with pytest.raises(ConfigurationError, match="surface heat capacity"):
        surface_coupled_adjustment(
            stable,
            capacity,
            pressure,
            surface_temperature=300.0,
            surface_heat_capacity=0.0,
            surface_pressure=P0,
        )


# --- adjusting toward a target other than the dry adiabat --------------------


@pytest.mark.req("ATM-002", "VAL-001")
def test_hydrostatic_heights_match_the_isothermal_closed_form():
    """An isothermal column has an exact answer, so the integrator is checked on it.

    ``z = (R T / g) ln(p_s / p)`` is the whole of the hypsometric equation when
    the temperature does not vary, which makes it the one case where the
    layer-by-layer integration can be compared against something rather than
    against itself.
    """
    isothermal = np.full(8, 260.0)
    interfaces = np.linspace(P0, 1.0e4, 9)
    pressure = 0.5 * (interfaces[:-1] + interfaces[1:])

    heights = hydrostatic_heights(isothermal, pressure, P0)
    scale = value_of("gas_constant_dry_air") * 260.0 / value_of("earth_standard_gravity")

    assert heights == pytest.approx(scale * np.log(P0 / pressure), rel=1e-12)
    assert np.all(np.diff(heights) > 0.0)


@pytest.mark.req("ATM-002", "CON-009", "VAL-001")
def test_a_profile_already_at_the_target_lapse_rate_is_left_alone():
    """The null control for a new target: it must do nothing where nothing is wrong.

    A column built to sit exactly on the target rate has no instability with
    respect to it, however far that rate is from a dry adiabat. Getting this
    wrong would make the adjustment relax a correct profile, which is the way a
    convective scheme most easily hides a bias.
    """
    pressure, capacity = _column()
    rate = np.full(pressure.size, 6.5e-3)  # K/m, constant
    T = np.linspace(288.0, 288.0, pressure.size)
    height = hydrostatic_heights(T, pressure, P0)
    offsets = lapse_rate_offsets(rate, height)
    on_target = 288.0 - offsets

    adjusted, groups = lapse_rate_adjustment_with_groups(on_target, capacity, offsets)

    assert groups == ()
    assert adjusted == pytest.approx(on_target, rel=1e-14)


@pytest.mark.req("ATM-002", "CON-009", "CLM-003")
def test_a_steeper_profile_relaxes_to_the_target_and_keeps_its_enthalpy():
    """The adjustment's one guarantee has to survive the change of target.

    Mixing to a prescribed lapse rate conserves the group's enthalpy for the
    same reason mixing to a dry adiabat does -- the mixed value is defined as
    the one that leaves ``sum(C T)`` unchanged -- so the new target must not
    weaken it. It is checked here against a profile steeper than the target,
    which is the case the scheme exists to act on.
    """
    pressure, capacity = _column()
    rate = np.full(pressure.size, 6.5e-3)
    steep = 300.0 - 9.0e-3 * np.linspace(0.0, 1.2e4, pressure.size)  # steeper than target
    height = hydrostatic_heights(steep, pressure, P0)
    offsets = lapse_rate_offsets(rate, height)

    adjusted, groups = lapse_rate_adjustment_with_groups(steep, capacity, offsets)

    assert groups  # there was something to do
    before = float(np.sum(capacity * steep))
    after = float(np.sum(capacity * adjusted))
    assert after == pytest.approx(before, rel=1e-14)

    # Inside a mixed group the profile now follows the target exactly.
    low, high = groups[0]
    span = slice(low, high + 1)
    assert np.allclose((adjusted + offsets)[span], (adjusted + offsets)[low], rtol=1e-12)

    # And it is less steep than it was, in K per metre through the real column.
    # Per-layer temperature steps will not do: the layers are equally spaced in
    # pressure and so grow thicker with height, which makes the step size grow
    # even where the lapse rate is constant.
    def rate_over(profile: np.ndarray) -> np.ndarray:
        return -np.diff(profile[span]) / np.diff(height[span])

    assert np.mean(rate_over(steep)) > 6.5e-3
    assert rate_over(adjusted) == pytest.approx(6.5e-3, rel=1e-9)


@pytest.mark.req("ATM-002", "DAT-007")
def test_the_target_offsets_must_cover_the_surface_and_every_layer():
    """A mismatched target is a silent bias, so it is refused rather than broadcast."""
    pressure, capacity = _column()
    T = np.linspace(290.0, 240.0, pressure.size)

    with pytest.raises(ConfigurationError, match="one more value than there are layers"):
        surface_coupled_adjustment(
            T,
            capacity,
            pressure,
            surface_temperature=300.0,
            surface_heat_capacity=1.0e7,
            surface_pressure=P0,
            target_offsets=np.zeros(T.size),  # missing the surface element
        )
    with pytest.raises(ConfigurationError, match="same layers"):
        lapse_rate_adjustment_with_groups(T, capacity, np.zeros(T.size - 1))
