"""Water vapour as a spectral absorber, and the feedback that follows.

Requirements EQR-002, EQR-004, RAD-001, SPC-001, VAL-001, VAL-002, INF-003.

The claim worth testing is not that a moist atmosphere is warmer -- any added
absorber does that -- but that the *feedback* is a consequence rather than an
input. Nothing in the constitutive relation knows about temperature. The water
module supplies a vapour path, Clausius-Clapeyron makes that path grow when the
column warms, and the amplification is what the radiative transfer computes from
those two facts. So the tests below hold the vapour fixed and let it follow, and
measure the difference.

Nothing here runs a coupled model or reads `data/raw`.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.contracts.module import StepContext
from mfespf.core.errors import ConfigurationError
from mfespf.core.quantities import (
    Field,
    Provenance,
    QuantityKind,
    SpatialSupport,
    TimeSemantics,
)
from mfespf.physics.atmosphere import ATMOSPHERE_TEMPERATURE_STATE, hydrostatic_layer_masses
from mfespf.physics.convection import ConvectiveAdjustmentModule
from mfespf.physics.greenhouse import (
    CO2_WING_DECAY_PER_METRE,
    CO2_WING_DECAY_WITH_VAPOUR_PER_METRE,
    VAPOUR_BAND_ABSORPTION,
    co2_band_absorption,
    co2_wing_decay_for,
    greenhouse_grid,
    vapour_absorption_optical_depth,
    vapour_mass_absorption,
)
from mfespf.physics.radiative_transfer import solve_column_radiation
from mfespf.physics.surface import SURFACE_TEMPERATURE_STATE
from mfespf.physics.thermodynamics import saturation_specific_humidity
from mfespf.physics.water import VAPOUR_STATE, PrescribedHumidityModule

pytestmark = pytest.mark.stage(6)

GRID = greenhouse_grid()
N_LAYERS = 10
INTERFACES, MASS = hydrostatic_layer_masses(
    n_layers=N_LAYERS, surface_pressure=101325.0, top_pressure=1.0e3, gravity=9.80665
)
PRESSURE = 0.5 * (INTERFACES[:-1] + INTERFACES[1:])
PROFILE = np.linspace(288.0, 220.0, N_LAYERS)
SHORTWAVE = GRID.centers < 4.0e-6
RELATIVE_HUMIDITY = 0.7

_WEIGHTS = np.where(SHORTWAVE, GRID.widths, 0.0)
INCIDENT = np.where(
    SHORTWAVE, 340.25 * _WEIGHTS / np.sum(_WEIGHTS) / np.maximum(GRID.widths, 1e-300), 0.0
)


def _vapour(temperature: np.ndarray) -> np.ndarray:
    """The vapour path a column at constant relative humidity carries."""
    return RELATIVE_HUMIDITY * saturation_specific_humidity(temperature, PRESSURE) * MASS


def _band(wavelength_m: float) -> int:
    return int(np.argmin(np.abs(GRID.centers - wavelength_m)))


def _transmittance(vapour_path: np.ndarray, wavelength_m: float) -> float:
    tau = vapour_absorption_optical_depth(GRID, vapour_path).sum(axis=0)
    return float(np.exp(-tau[_band(wavelength_m)]))


def _solve(shift: float, *, follow: bool, ppm: float = 284.7):
    temperature = PROFILE + shift
    vapour = _vapour(temperature) if follow else _vapour(PROFILE)
    return solve_column_radiation(
        grid=GRID,
        incident_spectral_irradiance=INCIDENT,
        layer_temperature=temperature,
        layer_mass=MASS,
        absorption_coefficient=co2_band_absorption(
            GRID, ppm, wing_decay_per_metre=co2_wing_decay_for(active_vapour=True), background=0.0
        ),
        surface_temperature=288.0 + shift,
        surface_albedo=np.where(SHORTWAVE, 0.3, 0.0),
        absorption_optical_depth=vapour_absorption_optical_depth(GRID, vapour),
    )


def _olr_with(vapour_path: np.ndarray, shift: float, ppm: float = 284.7) -> float:
    """Outgoing longwave for a vapour distribution supplied directly."""
    solution = solve_column_radiation(
        grid=GRID,
        incident_spectral_irradiance=INCIDENT,
        layer_temperature=PROFILE + shift,
        layer_mass=MASS,
        absorption_coefficient=co2_band_absorption(
            GRID, ppm, wing_decay_per_metre=co2_wing_decay_for(active_vapour=True), background=0.0
        ),
        surface_temperature=288.0 + shift,
        surface_albedo=np.where(SHORTWAVE, 0.3, 0.0),
        absorption_optical_depth=vapour_absorption_optical_depth(GRID, vapour_path),
    )
    return float(np.sum(solution.upward[-1][~SHORTWAVE]))


def _olr(shift: float, *, follow: bool, ppm: float = 284.7) -> float:
    return float(np.sum(_solve(shift, follow=follow, ppm=ppm).upward[-1][~SHORTWAVE]))


# --- the constitutive relation ----------------------------------------------


@pytest.mark.req("VAL-001", "SPC-001")
def test_a_dry_column_absorbs_nothing_and_the_window_is_open():
    """The null case, and it is the one that says the effect is the vapour's.

    A column with no vapour has no vapour opacity anywhere -- not a small
    residue, exactly none -- so anything the moist column does differently is
    attributable to the water and not to a floor someone left in.
    """
    tau = vapour_absorption_optical_depth(GRID, np.zeros(N_LAYERS))

    assert tau.shape == (N_LAYERS, GRID.n_bins)
    assert np.all(tau == 0.0)
    assert _transmittance(np.zeros(N_LAYERS), 10.0e-6) == 1.0


@pytest.mark.req("SPC-001", "RAD-001")
def test_the_bands_are_saturated_and_confined_to_the_infrared():
    """Two sign statements, not bookkeeping.

    Vapour absorbing sunlight would warm the atmosphere directly rather than
    trap the surface's emission, which is a different mechanism with a different
    signature; this configuration does not represent it and so must not produce
    a trace of it. And the response has to come from the window and the wings,
    which is only true if the band cores are already saturated.
    """
    kappa = vapour_mass_absorption(GRID)
    assert np.all(kappa[SHORTWAVE] == 0.0)

    column = _vapour(PROFILE).sum()
    assert _transmittance(_vapour(PROFILE), 6.3e-6) < 0.02  # the bending band
    assert _transmittance(_vapour(PROFILE), 30.0e-6) < 0.02  # the rotational band
    # Saturated means saturated: an optical depth of order a thousand, so the
    # exact band-centre coefficient stops mattering.
    assert VAPOUR_BAND_ABSORPTION * column > 1.0e3


@pytest.mark.req("EQR-002", "CLM-003")
def test_the_window_closes_as_the_column_moistens():
    """Where the feedback actually comes from.

    The bands are saturated, so they cannot trap more. What changes with vapour
    amount is the window: it is the only part of the spectrum through which the
    surface still reaches space, and the continuum closes it as the atmosphere
    moistens. That is the mechanism, and it is measured here rather than assumed.
    """
    dry, moist = _vapour(PROFILE), _vapour(PROFILE + 4.0)
    assert moist.sum() > dry.sum()

    open_window = _transmittance(dry, 10.0e-6)
    closed_window = _transmittance(moist, 10.0e-6)
    assert 0.5 < open_window < 0.8  # a window, not a hole and not a wall
    assert closed_window < open_window

    # Clausius-Clapeyron, which nothing here fitted: about 7 percent more vapour
    # per kelvin.
    per_kelvin = (moist.sum() / dry.sum()) ** 0.25 - 1.0
    assert 0.06 < per_kelvin < 0.09


@pytest.mark.req("EQR-002")
def test_optical_depth_is_linear_in_the_vapour_the_water_module_supplies():
    single = vapour_absorption_optical_depth(GRID, _vapour(PROFILE))
    double = vapour_absorption_optical_depth(GRID, 2.0 * _vapour(PROFILE))
    assert np.allclose(double, 2.0 * single, rtol=1e-12)


@pytest.mark.req("SYS-006", "MOD-002")
def test_unphysical_vapour_parameters_are_refused():
    with pytest.raises(ConfigurationError, match="cannot be negative"):
        vapour_absorption_optical_depth(GRID, -np.ones(N_LAYERS))
    with pytest.raises(ConfigurationError, match="non-negative"):
        vapour_mass_absorption(GRID, band_absorption=-1.0)
    with pytest.raises(ConfigurationError, match="non-negative"):
        vapour_mass_absorption(GRID, window_continuum=-1.0)
    with pytest.raises(ConfigurationError, match="wing decay"):
        vapour_mass_absorption(GRID, wing_decay_per_metre=0.0)


# --- what it was calibrated against -----------------------------------------


@pytest.mark.req("EQR-004", "INF-003", "VAL-002")
def test_the_outgoing_longwave_matches_the_flux_the_column_absorbs():
    """The calibration target is a radiation observation, not a temperature.

    CERES measures an all-sky outgoing longwave flux near 240 W m^-2. Fitting to
    that leaves the surface-temperature record untouched and so still available
    as independent validation, which is the whole of INF-003.

    **All-sky and not clear-sky, because of the albedo these columns carry.**
    An albedo of 0.3 is the all-sky planetary value, so the column absorbs
    340.25 * 0.7 = 238 W m^-2 and in equilibrium must emit exactly that. It was
    fitted to the clear-sky 265 first, which gave the column the cloud's cooling
    without the cloud's warming: it settled 7.5 K below its own calibration
    point, and the feedback measured there was half the right size, because
    Clausius-Clapeyron is exponential and a cold column is a dry one.
    """
    absorbed = 340.25 * (1.0 - 0.3)
    assert _olr(0.0, follow=True) == pytest.approx(absorbed, abs=2.0)

    # The greenhouse effect follows from it rather than being a second target:
    # sigma T^4 at 288 K is 390 W m^-2, and the observed all-sky value is ~155.
    assert 390.1 - _olr(0.0, follow=True) == pytest.approx(155.0, abs=6.0)


@pytest.mark.req("EQR-004", "SYS-006")
def test_the_co2_band_has_a_calibration_for_each_background_it_runs_against():
    """A band's fitted width depends on what else absorbs around it.

    With vapour active the 15 micron band is partly masked, so the same
    published forcing per doubling needs a narrower band than it does against
    the fixed background. Two calibrations coexist and neither is correct for
    the other configuration, so the selection is a function rather than a
    default somebody has to remember.
    """
    assert co2_wing_decay_for(active_vapour=False) == CO2_WING_DECAY_PER_METRE
    assert co2_wing_decay_for(active_vapour=True) == CO2_WING_DECAY_WITH_VAPOUR_PER_METRE
    assert CO2_WING_DECAY_WITH_VAPOUR_PER_METRE < CO2_WING_DECAY_PER_METRE

    def net(ppm: float) -> float:
        solution = _solve(0.0, follow=True, ppm=ppm)
        return float(np.sum(solution.downward[-1] - solution.upward[-1]))

    assert net(569.4) - net(284.7) == pytest.approx(3.70, abs=0.05)


# --- the feedback ------------------------------------------------------------


@pytest.mark.req("CLM-003", "VAL-002", "EQR-002")
def test_letting_the_vapour_follow_the_temperature_amplifies_the_response():
    """The result the whole mechanism exists to produce, and its null control.

    Nothing in the absorption relation knows about temperature. Holding the
    vapour path fixed while the column warms gives the response without the
    feedback; letting it follow constant relative humidity gives the response
    with it. The difference is the water-vapour feedback, and it is a
    consequence of two independent pieces of physics meeting -- Clausius-
    Clapeyron and a window that closes -- rather than a term anyone added.
    """
    held = _olr(0.5, follow=False) - _olr(-0.5, follow=False)
    free = _olr(0.5, follow=True) - _olr(-0.5, follow=True)

    assert held > 0.0  # a warmer column radiates more: the Planck response
    assert free > 0.0  # still stable -- the feedback does not run away
    assert free < held  # and it is amplifying

    feedback = held - free
    assert feedback == pytest.approx(2.29, abs=0.25)
    # Published clear-sky estimates are near 1.8 W m^-2 K^-1. This is larger
    # because the configuration has no lapse-rate feedback, which in the real
    # atmosphere cancels a substantial part of it -- declared, not tuned away.
    assert 1.5 < feedback < 3.0


@pytest.mark.req("CLM-003", "VAL-005", "WAT-002")
def test_vapour_confined_to_the_lowest_layer_has_no_feedback_at_all():
    """Why the coupled column does not get the feedback the relation can produce.

    The water module has no vertical transport -- it says so itself -- so
    evaporation enters the lowest layer and nothing carries it higher. In
    equilibrium the entire column's vapour sits in that layer at saturation.

    That is fatal, and not because the amount is small. A greenhouse absorber
    warms the surface only by emitting from somewhere colder than the surface.
    Vapour held a few kelvin below the surface re-emits at almost the
    temperature it absorbed at, so it traps nothing however much of it there is.
    The same column mass gives half the greenhouse effect and an amplification
    of exactly one.

    Pinned here rather than left in a run, because it is the difference between
    a feedback that is weak and a feedback that is structurally impossible, and
    only the second is worth changing the water module for.
    """
    distributed = _vapour(PROFILE)
    confined = np.zeros(N_LAYERS)
    confined[0] = distributed.sum()
    assert confined.sum() == pytest.approx(distributed.sum())  # the same water

    def amplification(build) -> float:
        held = _olr_with(build(PROFILE), 0.5) - _olr_with(build(PROFILE), -0.5)
        free = _olr_with(build(PROFILE + 0.5), 0.5) - _olr_with(build(PROFILE - 0.5), -0.5)
        return held / free

    spread = 390.1 - _olr_with(distributed, 0.0)
    stacked = 390.1 - _olr_with(confined, 0.0)
    assert spread > 1.9 * stacked  # roughly halved by being stacked at the bottom

    def confine(temperature: np.ndarray) -> np.ndarray:
        out = np.zeros(N_LAYERS)
        out[0] = _vapour(temperature).sum()
        return out

    assert amplification(_vapour) > 2.0  # the feedback the relation can produce
    assert amplification(confine) == pytest.approx(1.0, abs=1e-3)  # and none of it survives


# --- the prescribed-humidity control ----------------------------------------


def _humidity_inputs(temperature: np.ndarray) -> dict[str, Field]:
    """The three fields the control module reads, built directly."""
    common = {
        "quantity_kind": QuantityKind.DIAGNOSTIC,
        "time_semantics": TimeSemantics.INSTANTANEOUS,
        "timestamp": 0.0,
        "spatial_support": SpatialSupport.LAYER_AVERAGE,
        "dimensions": ("layer",),
        "provenance": Provenance(producer="test", derivation="test fixture"),
    }
    return {
        ATMOSPHERE_TEMPERATURE_STATE: Field(
            field_id=ATMOSPHERE_TEMPERATURE_STATE,
            physical_name="layer temperature",
            symbol="T",
            values=np.asarray(temperature, dtype=float),
            units="K",
            **common,  # type: ignore[arg-type]
        ),
        "atmosphere_layer_mass": Field(
            field_id="atmosphere_layer_mass",
            physical_name="layer mass",
            symbol="m",
            values=MASS,
            units="kg m^-2",
            **common,  # type: ignore[arg-type]
        ),
        "atmosphere_layer_pressure": Field(
            field_id="atmosphere_layer_pressure",
            physical_name="layer pressure",
            symbol="p",
            values=PRESSURE,
            units="Pa",
            **common,  # type: ignore[arg-type]
        ),
    }


def _stepped(module: PrescribedHumidityModule, temperature: np.ndarray, step: int) -> np.ndarray:
    module.begin_step(StepContext(time=float(step) * 3600.0, dt=3600.0, step_index=step))
    module.receive_inputs(_humidity_inputs(temperature))
    module.compute_tendencies()
    return np.asarray(module.emit_outputs()[VAPOUR_STATE].values, dtype=float)


@pytest.mark.req("EQR-002", "WAT-002", "CLM-003")
def test_the_prescribed_profile_follows_clausius_clapeyron():
    """What is prescribed is the distribution; what is not is the amount.

    Relative humidity is imposed uniformly, so the vapour a layer carries is a
    fixed fraction of what it *can* carry -- and how much it can carry comes
    from Clausius-Clapeyron. The seven percent per kelvin that results is
    therefore not imposed, which is the whole reason this counts as a control
    for a feedback rather than a way of writing one in.
    """
    module = PrescribedHumidityModule(n_layers=N_LAYERS, relative_humidity=0.7)
    module.validate_configuration()
    module.initialize({}, {"master_start_time": 0.0})

    cool = _stepped(module, PROFILE, 0)
    warm = _stepped(module, PROFILE + 4.0, 1)

    assert warm.sum() > cool.sum()
    per_kelvin = (warm.sum() / cool.sum()) ** 0.25 - 1.0
    assert 0.06 < per_kelvin < 0.09

    # And the relative humidity really is what was asked for, in every layer.
    published = np.asarray(module.emit_outputs()["atmosphere_relative_humidity"].values)
    assert np.allclose(published, 0.7)

    with pytest.raises(ConfigurationError, match="fraction"):
        PrescribedHumidityModule(relative_humidity=1.4).validate_configuration()


@pytest.mark.req("VAL-001", "CLM-003", "HYP-007")
def test_holding_the_vapour_removes_the_feedback_and_nothing_else():
    """The null control for the control.

    Freezing the vapour has to change exactly one thing: whether the humidity
    is allowed to answer the temperature. The frozen module sees the same warmer
    column and publishes the same vapour it published before, so a difference
    between the two runs is the feedback and cannot be anything else.
    """
    frozen = PrescribedHumidityModule(
        n_layers=N_LAYERS, relative_humidity=0.7, respond_to_temperature=False
    )
    frozen.initialize({}, {"master_start_time": 0.0})
    first = _stepped(frozen, PROFILE, 0)
    later = _stepped(frozen, PROFILE + 4.0, 1)
    assert np.array_equal(first, later)  # the column warmed; the vapour did not

    free = PrescribedHumidityModule(n_layers=N_LAYERS, relative_humidity=0.7)
    free.initialize({}, {"master_start_time": 0.0})
    assert np.array_equal(_stepped(free, PROFILE, 0), first)  # identical at the start
    assert _stepped(free, PROFILE + 4.0, 1).sum() > later.sum()  # and only then differ

    assert frozen.parameters()["respond_to_temperature"] is False
    assert frozen.contract.conserved_quantities == ()  # it moves no water and says so


@pytest.mark.req("MOD-006", "AUD-002")
def test_the_prescribed_profile_is_republished_every_step():
    """The trap a PRESCRIBED module falls into, and this one must not.

    ``emit_state`` runs at initialisation. A module that published its profile
    only there would freeze the humidity at whatever the first step implied and
    supply no feedback at all, while looking exactly like a module that did --
    which is how the cloud module once went wrong.
    """
    module = PrescribedHumidityModule(n_layers=N_LAYERS, relative_humidity=0.7)
    module.initialize({}, {"master_start_time": 0.0})

    warmer = [_stepped(module, PROFILE + shift, i) for i, shift in enumerate((0.0, 2.0, 4.0))]
    assert warmer[0].sum() < warmer[1].sum() < warmer[2].sum()

    # And it owns nothing: the vapour it publishes is not a state it advances.
    assert module.contract.owned_state_ids == frozenset()


# --- convective moisture transport ------------------------------------------


def _convection_inputs(
    temperature: np.ndarray, vapour: np.ndarray, surface: float | None = None
) -> dict[str, Field]:
    common = {
        "quantity_kind": QuantityKind.DIAGNOSTIC,
        "time_semantics": TimeSemantics.INSTANTANEOUS,
        "timestamp": 0.0,
        "spatial_support": SpatialSupport.LAYER_AVERAGE,
        "dimensions": ("layer",),
        "provenance": Provenance(producer="test", derivation="test fixture"),
    }
    scalar = {**common, "spatial_support": SpatialSupport.GLOBAL_SCALAR, "dimensions": ()}

    def layer(field_id, values, units):
        return Field(
            field_id=field_id,
            physical_name=field_id,
            symbol="x",
            values=np.asarray(values, dtype=float),
            units=units,
            **common,  # type: ignore[arg-type]
        )

    return {
        ATMOSPHERE_TEMPERATURE_STATE: layer(ATMOSPHERE_TEMPERATURE_STATE, temperature, "K"),
        "atmosphere_layer_heat_capacity": layer(
            "atmosphere_layer_heat_capacity", MASS * 1004.0, "J m^-2 K^-1"
        ),
        "atmosphere_layer_pressure": layer("atmosphere_layer_pressure", PRESSURE, "Pa"),
        "atmosphere_layer_mass": layer("atmosphere_layer_mass", MASS, "kg m^-2"),
        "layer_radiative_heating": layer("layer_radiative_heating", np.zeros(N_LAYERS), "W m^-2"),
        VAPOUR_STATE: layer(VAPOUR_STATE, vapour, "kg m^-2"),
        SURFACE_TEMPERATURE_STATE: Field(
            field_id=SURFACE_TEMPERATURE_STATE,
            physical_name="surface temperature",
            symbol="Ts",
            values=np.array(float(temperature[0]) if surface is None else surface),
            units="K",
            **scalar,  # type: ignore[arg-type]
        ),
        "surface_heat_capacity": Field(
            field_id="surface_heat_capacity",
            physical_name="surface heat capacity",
            symbol="Cs",
            values=np.array(1.0e7),
            units="J m^-2 K^-1",
            **scalar,  # type: ignore[arg-type]
        ),
    }


def _transport(
    temperature: np.ndarray,
    vapour: np.ndarray,
    *,
    surface_based: bool = False,
    surface: float | None = None,
):
    module = ConvectiveAdjustmentModule(
        transport_vapour=True,
        sensible_heat_coefficient=0.0,
        surface_based=surface_based,
        surface_pressure=float(INTERFACES[0]),
    )
    module.initialize({}, {"master_start_time": 0.0})
    module.begin_step(StepContext(time=0.0, dt=3600.0, step_index=0))
    module.receive_inputs(_convection_inputs(temperature, vapour, surface))
    tendencies = module.compute_tendencies()
    moved = [t for t in tendencies if t.state_id == VAPOUR_STATE]
    return module, (moved[0].as_array() * 3600.0 if moved else np.zeros(N_LAYERS))


@pytest.mark.req("CON-001", "WAT-002", "EQR-002")
def test_convection_carries_water_over_the_layers_it_mixes_and_conserves_it():
    """The same closure applied to a second conserved scalar.

    An overturning circulation homogenises whatever it carries, so the vapour is
    mixed to a uniform specific humidity over exactly the layers the temperature
    adjustment mixed -- no new group boundaries, and therefore no parameter the
    temperature adjustment did not already have.
    """
    unstable = np.array(PROFILE)
    # Twelve kelvin across layers about 880 m thick: steeper than the 9.8 K/km
    # dry adiabat, so the bottom four layers overturn.
    unstable[:4] = [300.0, 288.0, 276.0, 264.0]
    vapour = np.zeros(N_LAYERS)
    vapour[0] = 20.0

    module, moved = _transport(unstable, vapour)

    assert module._mixing_groups[0][0] == 0  # and the group reaches the layer holding it
    assert float(np.sum(moved)) == pytest.approx(0.0, abs=1e-12)  # water is redistributed
    assert moved[0] < 0.0  # out of the layer that held it
    assert np.any(moved[1:] > 0.0)  # and into the ones above

    # Uniform specific humidity over the mixed group, which is what "mixed" means.
    low, high = module._mixing_groups[0]
    final = (vapour + moved)[low : high + 1]
    specific = final / MASS[low : high + 1]
    assert np.allclose(specific, specific[0])


@pytest.mark.req("CLM-003", "VAL-005", "CON-001")
def test_the_transport_moves_nothing_when_the_adjustment_is_inactive():
    """The negative result this mechanism ran into, pinned so it cannot be lost.

    Keying moisture transport to the convective closure reuses a scheme that has
    no free parameters, which is why it was worth trying first. It also inherits
    the closure's silence: where the profile is stable the adjustment does
    nothing, and so does the transport.

    That is not a corner case in this column, it is the equilibrium. Vapour with
    nothing to carry it upward accumulates in the lowest layer; radiatively
    active, it makes that layer an efficient radiator; it cools 45 K below the
    layer above; and the inversion that results decouples it from the neutral
    region higher up. The adjustment is then inactive exactly where the water
    is. The missing transport sustains the stratification that blocks the
    transport, and adding the transport alone does not break the loop.
    """
    stable = 290.0 - 4.0 * np.arange(N_LAYERS)  # gentler than the adiabat: stable
    vapour = np.zeros(N_LAYERS)
    vapour[0] = 20.0

    module, moved = _transport(stable, vapour)

    assert module._mixing_groups == ()
    assert np.all(moved == 0.0)
    assert module._vapour_moved == 0.0


@pytest.mark.req("CON-001", "WAT-002", "CLM-003")
def test_a_buoyant_surface_lifts_the_water_the_layer_scheme_could_not_reach():
    """The mechanism that answers the negative result above.

    Same column, same trapped water, same inversion -- and the only change is
    that the mixing is launched from the surface rather than from a comparison
    between neighbouring layers. The surface is buoyant against the whole
    profile, so it cannot be hidden by a step between two layers, and the group
    it opens reaches the layer holding the water.

    Both schemes are checked on the identical input, because the claim is about
    the difference between them and not about either one alone.
    """
    inverted = np.array(PROFILE)
    inverted[:2] = [240.0, 286.0]  # trapped vapour cools layer 0 below layer 1
    vapour = np.zeros(N_LAYERS)
    vapour[0] = 20.0

    _, stuck = _transport(inverted, vapour, surface=300.0)
    assert np.all(stuck == 0.0)  # the layer-only scheme, for the same column

    module, moved = _transport(inverted, vapour, surface_based=True, surface=300.0)

    assert module._surface_mixed
    assert float(np.sum(moved)) == pytest.approx(0.0, abs=1e-12)  # still a redistribution
    assert moved[0] < -1.0  # and now the water actually leaves layer 0
    assert np.sum(moved[1:]) > 1.0
    assert module._surface_convective_flux > 0.0  # paid for by the surface, upward


@pytest.mark.req("VAL-002", "CLM-003", "WAT-002")
def test_complete_mixing_oversaturates_the_top_of_the_group_it_moistens():
    """Why the feedback this scheme produces is a property of the layer grid.

    Mixing to a uniform specific humidity is what an overturning circulation
    does to a conserved scalar, but a real plume entrains dry air as it rises
    and detrains before it saturates, which is why the observed mid-troposphere
    sits near thirty to fifty percent rather than at saturation. Mixing
    completely instead drives the top of the group past saturation, and the
    water module's condensation then clips it back.

    So the water that survives aloft is set by where the group's top happens to
    fall, which is a property of the discretisation rather than of the physics.
    Measured in the coupled column, the amplification runs 1.10x, 1.22x, 1.88x,
    1.28x and 1.30x at 10, 20, 30, 40 and 60 layers, with the spin-up converged
    to three decimals at every one. The two finest grids agree, so a limit near
    1.3x is likely; the 45 percent excursion at 30 layers is not explained, and
    until it is the feedback is recorded rather than quoted.
    """
    inverted = np.array(PROFILE)
    inverted[:2] = [240.0, 286.0]
    vapour = np.zeros(N_LAYERS)
    vapour[0] = 20.0

    module, moved = _transport(inverted, vapour, surface_based=True, surface=300.0)
    low, high = module._mixing_groups[0]

    specific = (vapour + moved)[low : high + 1] / MASS[low : high + 1]
    saturation = saturation_specific_humidity(inverted[low : high + 1], PRESSURE[low : high + 1])

    assert np.allclose(specific, specific[0])  # mixed, by construction
    assert specific[-1] > saturation[-1]  # and past saturation at the top


@pytest.mark.req("CON-001", "CON-007", "WAT-001")
def test_over_drawing_a_layer_is_refused_rather_than_clipped():
    """Clipping a negative vapour to zero creates water, and bought nothing.

    Two modules can each compute a vapour tendency from the same starting state
    without seeing the other's, and together draw more from a layer than it
    holds. That happened on a column driven to a realistic evaporation rate:
    convective transport moved vapour out of a layer while the saturation
    adjustment condensed more out of it, and a single clip of 0.62 kg m^-2
    surfaced 1475 steps later as a 5.4 percent water residual with nothing to
    say where it came from.

    The reconciliation catches the invented mass either way, so the clip only
    ever obscured the cause. The state is now refused at the step that creates
    it, naming the layer, the over-draw and the labels responsible.
    """
    from mfespf.core.errors import ConservationError
    from mfespf.physics.equations import EQUATION_IDS
    from mfespf.physics.water import WaterModule
    from mfespf.registries.states import TendencyContribution

    module = WaterModule(n_layers=3)
    module.vapour = np.array([2.0, 1.0, 0.0])
    dt = 100.0

    def draw(rate: np.ndarray, label: str) -> TendencyContribution:
        return TendencyContribution(
            state_id=VAPOUR_STATE,
            provider_module_id="test",
            label=label,
            value=rate,
            units="kg m^-2",
            equation_id=EQUATION_IDS["condensation"],
        )

    # Comfortably inside what the layers hold: committed unchanged.
    module.advance_state(dt, [draw(np.array([-0.005, -0.005, 0.0]), "condensation")])
    assert module.vapour == pytest.approx([1.5, 0.5, 0.0])

    # Two contributions that are each affordable and jointly are not.
    module.vapour = np.array([2.0, 1.0, 0.0])
    with pytest.raises(ConservationError) as raised:
        module.advance_state(
            dt,
            [
                draw(np.array([0.0, -0.006, 0.006]), "convective_vapour_transport"),
                draw(np.array([0.0, -0.006, 0.0]), "evaporation_minus_condensation"),
            ],
        )

    message = str(raised.value)
    assert "layer 1" in message
    assert "convective_vapour_transport" in message
    assert "evaporation_minus_condensation" in message
    # The vapour is left untouched rather than half-committed.
    assert module.vapour == pytest.approx([2.0, 1.0, 0.0])
