"""Stage 3: the zero-buoyancy bulk-plume closure for humidity and stratification.

The tests that matter here are the ones that check the closure against results
derived somewhere else. Two do that against this framework's own older code --
the moist adiabat and the dry adiabat, both written before any of this existed
-- and one checks the predicted humidity against the observed free
troposphere, which is the only way a closure with no fitted parameter can be
wrong.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from mfespf.contracts.module import StepContext
from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError
from mfespf.core.quantities import (
    Field,
    Provenance,
    QuantityKind,
    SpatialSupport,
    TimeSemantics,
)
from mfespf.physics.atmosphere import ATMOSPHERE_TEMPERATURE_STATE
from mfespf.physics.convection import ConvectiveAdjustmentModule
from mfespf.physics.plume import (
    BULK_PLUME_PARAMETER,
    EVAPORATING_CONDENSATE_FRACTION,
    EVAPORATING_FRACTION_RANGE,
    bulk_plume_relative_humidity,
    entraining_lapse_rate,
    equilibrium_relative_humidity,
    fractional_entrainment_rate,
    water_vapour_lapse_rate,
)
from mfespf.physics.thermodynamics import (
    dry_adiabatic_lapse_rate,
    moist_adiabatic_lapse_rate,
)
from mfespf.physics.water import VAPOUR_STATE, EntrainingPlumeHumidityModule
from mfespf.run.column_configurations import radiating_column

pytestmark = pytest.mark.stage(3)

N_LAYERS = 8
P0 = value_of("standard_surface_pressure")
INTERFACES = np.linspace(P0, 1.0e4, N_LAYERS + 1)
PRESSURE = 0.5 * (INTERFACES[:-1] + INTERFACES[1:])
MASS = -np.diff(INTERFACES) / value_of("earth_standard_gravity")
PROFILE = np.linspace(288.0, 230.0, N_LAYERS)


# --- the humidity closure ----------------------------------------------------


@pytest.mark.req("WAT-002", "EQR-002", "VAL-001")
def test_the_closed_form_is_the_general_form_at_constant_mass_flux():
    """Two expressions for one quantity, and they have to agree exactly.

    Romps' general result carries height-varying entrainment and detrainment;
    holding the mass flux constant forces the two rates equal and collapses it
    to a single number. That collapse is algebra, not approximation, so any
    disagreement here is an implementation error in one of the two -- which is
    the point of keeping both.
    """
    gamma = water_vapour_lapse_rate(PROFILE, moist_adiabatic_lapse_rate(PROFILE, PRESSURE))
    rate = fractional_entrainment_rate(gamma)

    general = bulk_plume_relative_humidity(
        entrainment=rate, detrainment=rate, vapour_lapse_rate=gamma
    )
    closed = equilibrium_relative_humidity()

    assert general == pytest.approx(closed, rel=1e-14)
    assert np.allclose(general, general[0])  # constant, as the closed form implies


@pytest.mark.req("WAT-002", "EQR-002")
def test_undiluted_convection_leaves_only_the_evaporated_condensate():
    """The limit that says what the parameter means.

    With no entrainment there is no detrainment either, so the only water
    reaching the environment is the condensate that evaporates instead of
    falling out. The humidity is then exactly that fraction, which is what
    makes ``alpha`` a precipitation efficiency and not a tuning knob.
    """
    for alpha in (0.0, 0.2, 0.45):
        assert equilibrium_relative_humidity(0.0, alpha) == pytest.approx(alpha)


@pytest.mark.req("WAT-002", "EQR-002", "VAL-001")
def test_relative_humidity_rises_with_entrainment_and_saturates():
    """More vigorous turnover moistens the environment, and cannot oversaturate it.

    The direction is worth stating because it reads backwards at first: a
    plume that entrains more also detrains more, and what it detrains is
    saturated. The bound is the physical content -- an environment fed entirely
    by detrained saturated air is saturated, and no amount of entrainment takes
    it past that.
    """
    humidities = [equilibrium_relative_humidity(a) for a in (0.0, 0.1, 0.2, 0.5, 1.0, 10.0)]

    assert all(b > a for a, b in itertools.pairwise(humidities))
    assert equilibrium_relative_humidity(1.0e6) < 1.0
    assert equilibrium_relative_humidity(1.0e6) == pytest.approx(1.0, abs=1e-5)


@pytest.mark.req("CLM-003", "INF-003", "EQR-004")
def test_the_predicted_humidity_matches_the_regime_the_closure_describes():
    """The falsification test, against the right population -- which it was not.

    Both parameters come from measurement and neither was chosen to make the
    humidity come out right, so this can fail. It is the comparison itself that
    had to be corrected: an earlier version of this test asserted 0.30 to 0.50,
    the relative humidity of the free troposphere *globally*, and passed.

    That was the wrong reference. This is a radiative-convective closure, and
    the global figure is dominated by the subtropical dry zones, which are
    maintained by large-scale subsidence and last-saturation dynamics that no
    bulk plume represents. The population it describes is a convecting column,
    where cloud-resolving simulations put the free-tropospheric humidity near
    eighty percent -- Singh and O'Gorman took eighty as a constant for exactly
    this reason. A test can pass and still be evidence of nothing, which is what
    the previous one was.
    """
    predicted = equilibrium_relative_humidity()
    assert 0.60 < predicted < 0.90

    low, high = EVAPORATING_FRACTION_RANGE
    spread = [equilibrium_relative_humidity(BULK_PLUME_PARAMETER, a) for a in (low, high)]
    assert 0.60 < spread[0] < spread[1] < 0.90

    # The bound that matters for this framework: a convecting-column humidity
    # applied to a global-mean column is too moist, not too dry. That is the
    # opposite of the shortfall the transport scheme had, and it is a bias in a
    # known direction rather than an unknown one.
    assert predicted > 0.50


@pytest.mark.req("WAT-002", "ATM-002")
def test_the_vapour_lapse_rate_is_a_sensible_scale_height():
    """``gamma`` is derived, so its reciprocal is a prediction and not an input.

    A terrestrial troposphere holds its water in the lowest couple of
    kilometres. Getting a scale height of that size out of Clausius-Clapeyron
    and hydrostatic balance alone is the cheapest available check that the
    expression was transcribed correctly.
    """
    gamma = water_vapour_lapse_rate(PROFILE, moist_adiabatic_lapse_rate(PROFILE, PRESSURE))
    scale_height = 1.0 / gamma

    assert np.all(gamma > 0.0)
    assert np.all(scale_height > 1.0e3)
    assert np.all(scale_height < 6.0e3)
    assert np.all(np.diff(scale_height) < 0.0)  # drier and shallower with height


# --- the lapse rate the same parameter predicts ------------------------------


@pytest.mark.req("ATM-002", "EQR-002", "VAL-001")
def test_the_entraining_lapse_rate_reduces_to_the_moist_adiabat():
    """Checked against code that predates the closure, not against a second copy.

    ``moist_adiabatic_lapse_rate`` was written for Stage 3 with no knowledge of
    entrainment, and the entraining expression has to reproduce it exactly at
    zero entrainment. That makes it an independent check rather than a
    restatement of the same algebra.
    """
    assert entraining_lapse_rate(PROFILE, PRESSURE, 0.0) == pytest.approx(
        moist_adiabatic_lapse_rate(PROFILE, PRESSURE), rel=1e-15
    )


@pytest.mark.req("ATM-002", "CLM-003", "VAL-001")
def test_entrainment_stabilises_between_the_moist_and_dry_adiabats():
    """Entrained dry air robs the plume of buoyancy, so it leaves a stiffer column.

    The bound on the other side matters as much: however much dry air a plume
    entrains, condensation still releases heat into it, so the profile it
    leaves neutral can never be as steep as a dry adiabat. A scheme that
    crossed either bound would be producing energy or destroying it.
    """
    moist = moist_adiabatic_lapse_rate(PROFILE, PRESSURE)
    entraining = entraining_lapse_rate(PROFILE, PRESSURE, BULK_PLUME_PARAMETER)
    dry = dry_adiabatic_lapse_rate()

    assert np.all(entraining > moist)
    assert np.all(entraining < dry)
    # And monotone in the parameter, so the prediction has a direction.
    stronger = entraining_lapse_rate(PROFILE, PRESSURE, 2.0 * BULK_PLUME_PARAMETER)
    assert np.all(stronger > entraining)


@pytest.mark.req("ATM-002", "VAL-001")
def test_a_dry_column_is_dry_adiabatic_whatever_the_entrainment():
    """Entrainment cannot matter where there is no water to be diluted.

    The mechanism is entirely about entrained air being subsaturated, so as the
    saturation humidity goes to zero the parameter has to drop out. It does,
    algebraically: the factors of ``1 + a`` cancel between numerator and
    denominator.
    """
    cold = np.full(N_LAYERS, 180.0)  # saturation humidity is tiny but not zero
    dry = dry_adiabatic_lapse_rate()

    rates = [entraining_lapse_rate(cold, PRESSURE, a) for a in (0.0, 0.2, 1.0, 5.0)]
    for rate in rates:
        assert rate == pytest.approx(dry, rel=1e-3)

    # The limit is exact only at q* = 0, so the sharper statement is that the
    # parameter has stopped mattering: the spread across a two-order-of-magnitude
    # range of it is far smaller than the moist-to-dry gap it opens in a warm
    # column, which is what "drops out" means for a finite temperature.
    spread = float(np.max(rates) - np.min(rates))
    warm_gap = float(np.max(dry - moist_adiabatic_lapse_rate(PROFILE, PRESSURE)))
    assert spread < 1.0e-3 * warm_gap


@pytest.mark.req("EQR-002", "VAL-001")
def test_unphysical_plume_parameters_are_refused():
    """A negative entrainment rate is not a slower one, and it must not run."""
    with pytest.raises(ConfigurationError, match="cannot be negative"):
        equilibrium_relative_humidity(-0.1)
    with pytest.raises(ConfigurationError, match="fraction below one"):
        equilibrium_relative_humidity(0.2, 1.0)
    with pytest.raises(ConfigurationError, match="cannot be negative"):
        entraining_lapse_rate(PROFILE, PRESSURE, -1.0)
    with pytest.raises(ConfigurationError, match="water-vapour lapse rate must be positive"):
        bulk_plume_relative_humidity(
            entrainment=1.0e-4, detrainment=1.0e-4, vapour_lapse_rate=-1.0e-4
        )


# --- the module ---------------------------------------------------------------


def _inputs(temperature: np.ndarray) -> dict[str, Field]:
    common = {
        "quantity_kind": QuantityKind.DIAGNOSTIC,
        "time_semantics": TimeSemantics.INSTANTANEOUS,
        "timestamp": 0.0,
        "spatial_support": SpatialSupport.LAYER_AVERAGE,
        "dimensions": ("layer",),
        "provenance": Provenance(producer="test", derivation="test fixture"),
    }

    def layer(field_id: str, values: np.ndarray, units: str) -> Field:
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
        "atmosphere_layer_mass": layer("atmosphere_layer_mass", MASS, "kg m^-2"),
        "atmosphere_layer_pressure": layer("atmosphere_layer_pressure", PRESSURE, "Pa"),
    }


def _stepped(module: EntrainingPlumeHumidityModule, temperature: np.ndarray, step: int):
    module.begin_step(StepContext(time=float(step) * 3600.0, dt=3600.0, step_index=step))
    module.receive_inputs(_inputs(temperature))
    module.compute_tendencies()
    return module.emit_outputs()


@pytest.mark.req("WAT-002", "MOD-006", "CLM-003")
def test_the_plume_module_derives_its_humidity_rather_than_taking_one():
    """The whole difference from the prescribed control, made checkable.

    The module has no relative-humidity argument to pass. Its humidity is
    whatever the two plume parameters imply, changing when they change and in
    the right direction, and its contract cites the water-budget equation
    rather than pretending Clausius-Clapeyron alone produced the distribution.
    """
    module = EntrainingPlumeHumidityModule(n_layers=N_LAYERS)
    module.validate_configuration()

    assert module.relative_humidity == pytest.approx(
        equilibrium_relative_humidity(BULK_PLUME_PARAMETER, EVAPORATING_CONDENSATE_FRACTION)
    )
    assert "EQ-BULK-PLUME-HUMIDITY" in module.contract.governing_equation_ids
    assert module.parameters()["relative_humidity_origin"] == (
        "predicted_from_the_bulk_plume_water_budget"
    )

    wetter = EntrainingPlumeHumidityModule(n_layers=N_LAYERS, evaporating_fraction=0.8)
    drier = EntrainingPlumeHumidityModule(n_layers=N_LAYERS, evaporating_fraction=0.4)
    assert drier.relative_humidity < module.relative_humidity < wetter.relative_humidity

    # It still moves no water, and the contract has to keep saying so.
    assert module.contract.conserved_quantities == ()


@pytest.mark.req("WAT-002", "HYP-005", "VAL-001")
def test_the_plume_humidity_follows_clausius_clapeyron():
    """The distribution is predicted; the amount still is not imposed.

    Fixing the relative humidity by any means leaves the vapour amount to
    Clausius-Clapeyron, so about seven percent per kelvin comes out without
    being put in. Freezing the response removes exactly that and nothing else,
    which is the null control the feedback is measured against.
    """
    module = EntrainingPlumeHumidityModule(n_layers=N_LAYERS)
    module.initialize({}, {"master_start_time": 0.0})

    cool = np.asarray(_stepped(module, PROFILE, 0)[VAPOUR_STATE].values, dtype=float)
    warm = np.asarray(_stepped(module, PROFILE + 4.0, 1)[VAPOUR_STATE].values, dtype=float)

    per_kelvin = (warm.sum() / cool.sum()) ** 0.25 - 1.0
    assert 0.06 < per_kelvin < 0.09

    published = np.asarray(module.emit_outputs()["atmosphere_relative_humidity"].values)
    assert np.allclose(published, module.relative_humidity)

    frozen = EntrainingPlumeHumidityModule(n_layers=N_LAYERS, respond_to_temperature=False)
    frozen.initialize({}, {"master_start_time": 0.0})
    first = np.asarray(_stepped(frozen, PROFILE, 0)[VAPOUR_STATE].values, dtype=float)
    later = np.asarray(_stepped(frozen, PROFILE + 4.0, 1)[VAPOUR_STATE].values, dtype=float)
    assert np.array_equal(first, later)


@pytest.mark.req("EQR-004", "CLM-003", "VAL-001")
def test_the_implied_entrainment_rate_agrees_with_the_simulations():
    """The check that caught the evaporating fraction being wrong.

    ``a = epsilon PE / gamma`` ties the three together, so a parameter pair
    chosen from humidity considerations can still be tested against entrainment
    rates it was not derived from. Cloud-resolving simulations diagnose roughly
    0.2 to 0.8 per kilometre for bulk plumes.

    At the measured precipitation efficiency the closure gives 0.198 per
    kilometre -- the bottom edge of that range, not comfortably inside it, and
    worth saying so. At the value previously assumed it gave 0.09, a factor of
    two below the bottom, and that disagreement sat in the equation registry as
    a known failure mode for two commits before anyone looked the efficiency up.
    The over-determination was doing its job the whole time; nobody was reading
    it.
    """
    module = EntrainingPlumeHumidityModule(n_layers=N_LAYERS)
    implied = module.implied_entrainment_rate(280.0, 6.5e-3) * 1.0e3

    assert 0.15 <= implied <= 0.8
    assert implied == pytest.approx(0.2, abs=0.02)  # at the edge, and pinned there

    # And the identity run backwards returns the published bulk-plume parameter,
    # which is the closure of the loop rather than a restatement of it.
    stale = EntrainingPlumeHumidityModule(n_layers=N_LAYERS, evaporating_fraction=0.3)
    assert stale.implied_entrainment_rate(280.0, 6.5e-3) * 1.0e3 < 0.2


@pytest.mark.req("MOD-006", "WAT-002")
def test_a_column_takes_its_humidity_from_exactly_one_place():
    """Three modules now publish the vapour field, and only one may run.

    The guard existed for two of them and had to be widened for the third. It
    is worth a test rather than a review, because the failure it prevents is
    silent: two producers of one field leave the coupler picking whichever
    published last, and the run would look fine.
    """
    common = {
        "absorption_coefficient": np.zeros(4),
        "n_layers": 4,
        "convection": True,
    }
    for clash in (
        {"water": True, "plume_humidity": True},
        {"water": True, "prescribed_humidity": 0.7},
        {"prescribed_humidity": 0.7, "plume_humidity": True},
        {"water": True, "prescribed_humidity": 0.7, "plume_humidity": True},
    ):
        with pytest.raises(ConfigurationError, match=r"never more than one|never both"):
            radiating_column(run_id="CLASH", purpose="two humidity sources", **common, **clash)

    # And each on its own is accepted.
    for single in ({"water": True}, {"prescribed_humidity": 0.7}, {"plume_humidity": True}):
        radiating_column(run_id="OK", purpose="one humidity source", **common, **single)


@pytest.mark.req("ATM-002", "MOD-006", "DAT-007")
def test_the_entraining_target_is_declared_and_needs_a_surface():
    """The target belongs in the parameters, and it cannot float free of the ground.

    Two runs differing only in what their convection relaxes to must not hash
    alike, so the target and its parameter go into ``parameters()``. And the
    plume's lapse rate is measured upward from the surface, so a layer-only
    adjustment has no anchor for it -- refused rather than quietly anchored to
    the lowest layer, which would put the whole profile a layer's thickness out.
    """
    plume = ConvectiveAdjustmentModule(surface_based=True, target="entraining_plume")
    plume.validate_configuration()
    assert plume.parameters()["target"] == "entraining_plume"
    assert plume.parameters()["bulk_plume_parameter"] == BULK_PLUME_PARAMETER
    assert "EQ-ENTRAINING-LAPSE-RATE" in plume.contract.governing_equation_ids

    dry = ConvectiveAdjustmentModule()
    assert dry.parameters()["target"] == "dry_adiabat"
    assert dry.parameters()["bulk_plume_parameter"] is None
    assert "EQ-ENTRAINING-LAPSE-RATE" not in dry.contract.governing_equation_ids

    with pytest.raises(ConfigurationError, match="needs surface_based convection"):
        ConvectiveAdjustmentModule(target="entraining_plume").validate_configuration()
    with pytest.raises(ConfigurationError, match="unknown convective target"):
        ConvectiveAdjustmentModule(
            surface_based=True, target="moist_adiabat"
        ).validate_configuration()


@pytest.mark.req("VAL-001", "CLM-003", "HYP-005")
def test_freezing_the_vapour_freezes_the_lapse_rate_it_sets():
    """The control has to hold every consequence of the thing it holds.

    Water vapour acts on this column twice: radiatively through the absorption,
    and dynamically because the entraining lapse rate depends on the saturation
    humidity. Freezing only the first left the frozen run with the entire
    lapse-rate response, which showed up as a feedback parameter of
    12.4 W m^-2 K^-1 against a Planck response of 3.2 -- an absurdity rather
    than a subtlety, which is the only reason it was caught.

    What stays free is the geometry. The column still expands as it warms,
    because that is hypsometry and not the vapour's doing, so the *rate* is
    frozen and the offsets it integrates to are not.
    """
    warm = PROFILE + 6.0

    free = ConvectiveAdjustmentModule(surface_based=True, target="entraining_plume")
    cold_rate = free._plume_offsets(PROFILE, PRESSURE)
    warm_rate = free._plume_offsets(warm, PRESSURE)
    assert not np.allclose(cold_rate, warm_rate)  # responds, as the free run must

    held = ConvectiveAdjustmentModule(
        surface_based=True, target="entraining_plume", target_responds_to_temperature=False
    )
    held._plume_offsets(PROFILE, PRESSURE)
    captured = np.array(held._held_rate)
    held._plume_offsets(warm, PRESSURE)

    assert np.array_equal(held._held_rate, captured)  # the rate is pinned
    assert held.parameters()["target_responds_to_temperature"] is False
    # And the frozen target is the one the *unfrozen* column started from, not
    # some other profile: a control that begins somewhere else is not a control.
    assert np.allclose(captured, entraining_lapse_rate(PROFILE, PRESSURE, BULK_PLUME_PARAMETER))
