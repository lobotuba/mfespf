"""Condensing what convection lifts, in the step it lifts it.

docs/protocols/convective-condensation-order.md, test O1. In the default order
the water module condenses first and convection then mixes the vapour to a
uniform specific humidity over each group, so the top of a group ends the step
far above its own saturation -- 63166 layer-steps in the tropical column, up to
300 times saturation (tropical-supersaturation-cause.md). With
``condense_after_convection`` convection runs first and publishes its mixing,
and the water module's saturation adjustment acts on the mixed profile, as it
already does on the upwelling's.
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
from mfespf.physics.surface import SURFACE_TEMPERATURE_STATE
from mfespf.physics.thermodynamics import (
    deposition_latent_heat_excess,
    saturation_humidity_within_domain,
)
from mfespf.physics.water import VAPOUR_STATE, WaterModule
from mfespf.run.column_configurations import reference_column

pytestmark = pytest.mark.stage(6)

N = 30
DT = 3600.0
INTERFACES, MASS = hydrostatic_layer_masses(
    n_layers=N, surface_pressure=101325.0, top_pressure=1.0e3, gravity=9.80665
)
PRESSURE = 0.5 * (INTERFACES[:-1] + INTERFACES[1:])
C = MASS * 1004.0


def _field(field_id: str, values, units: str, *, layer: bool = True, **extra) -> Field:
    return Field(
        field_id=field_id,
        physical_name=field_id,
        symbol="x",
        values=np.asarray(values, dtype=float),
        units=units,
        quantity_kind=extra.pop("quantity_kind", QuantityKind.DIAGNOSTIC),
        time_semantics=extra.pop("time_semantics", TimeSemantics.INSTANTANEOUS),
        timestamp=0.0,
        spatial_support=SpatialSupport.LAYER_AVERAGE if layer else SpatialSupport.GLOBAL_SCALAR,
        dimensions=("layer",) if layer else (),
        provenance=Provenance(producer="test", derivation="test fixture"),
        **extra,
    )


def _unstable_moist_column() -> tuple[np.ndarray, np.ndarray]:
    """Four layers steeper than the dry adiabat, all the water in the lowest one."""
    temperature = np.linspace(288.0, 220.0, N)
    temperature[:4] = [300.0, 288.0, 276.0, 264.0]
    vapour = np.zeros(N)
    vapour[0] = 20.0
    return temperature, vapour


def _convect(temperature: np.ndarray, vapour: np.ndarray) -> ConvectiveAdjustmentModule:
    module = ConvectiveAdjustmentModule(
        transport_vapour=True,
        condensation_follows=True,
        sensible_heat_coefficient=0.0,
        surface_pressure=float(INTERFACES[0]),
    )
    module.initialize({}, {"master_start_time": 0.0})
    module.begin_step(StepContext(time=0.0, dt=DT, step_index=0))
    module.receive_inputs(
        {
            ATMOSPHERE_TEMPERATURE_STATE: _field(ATMOSPHERE_TEMPERATURE_STATE, temperature, "K"),
            "atmosphere_layer_heat_capacity": _field(
                "atmosphere_layer_heat_capacity", C, "J m^-2 K^-1"
            ),
            "atmosphere_layer_pressure": _field("atmosphere_layer_pressure", PRESSURE, "Pa"),
            "atmosphere_layer_mass": _field("atmosphere_layer_mass", MASS, "kg m^-2"),
            "layer_radiative_heating": _field("layer_radiative_heating", np.zeros(N), "W m^-2"),
            VAPOUR_STATE: _field(VAPOUR_STATE, vapour, "kg m^-2"),
            SURFACE_TEMPERATURE_STATE: _field(
                SURFACE_TEMPERATURE_STATE, float(temperature[0]), "K", layer=False
            ),
            "surface_heat_capacity": _field(
                "surface_heat_capacity", 1.0e7, "J m^-2 K^-1", layer=False
            ),
        }
    )
    module.compute_tendencies()
    return module


@pytest.mark.req("CON-001", "WAT-002", "SYS-006")
def test_what_convection_lifts_above_saturation_condenses_in_the_same_step():
    """O1: the lifted layer ends at its saturation, and the excess rains out."""
    temperature, vapour = _unstable_moist_column()
    convection = _convect(temperature, vapour)
    published = convection.emit_outputs()["convective_vapour_transport_rate"]
    moved = np.asarray(published.values, dtype=float) * DT
    assert float(np.sum(moved)) == pytest.approx(0.0, abs=1e-12)
    low, high = convection._mixing_groups[0]
    assert low == 0 and high >= 1

    water = WaterModule(n_layers=N, evaporation_coefficient=0.0, receives_convective_transport=True)
    water.initialize({VAPOUR_STATE: vapour}, {"master_start_time": 0.0})
    water.begin_step(StepContext(time=0.0, dt=DT, step_index=0))
    water.receive_inputs(
        {
            ATMOSPHERE_TEMPERATURE_STATE: _field(ATMOSPHERE_TEMPERATURE_STATE, temperature, "K"),
            "atmosphere_layer_mass": _field("atmosphere_layer_mass", MASS, "kg m^-2"),
            "atmosphere_layer_pressure": _field("atmosphere_layer_pressure", PRESSURE, "Pa"),
            "atmosphere_layer_heat_capacity": _field(
                "atmosphere_layer_heat_capacity", C, "J m^-2 K^-1"
            ),
            SURFACE_TEMPERATURE_STATE: _field(
                SURFACE_TEMPERATURE_STATE, float(temperature[0]), "K", layer=False
            ),
            "convective_vapour_transport_rate": published,
        }
    )
    contributions = water.compute_tendencies()

    saturation, _ = saturation_humidity_within_domain(temperature, PRESSURE)
    capacity = saturation * MASS
    mixed = vapour + moved
    expected = np.maximum(0.0, mixed - capacity)
    # The mixing lifted water above saturation somewhere above the layer that
    # held it -- the case the default order leaves for the next step.
    assert np.any(expected[1:] > 0.0)
    np.testing.assert_allclose(water._condensation, expected, rtol=1e-12, atol=1e-15)

    # The committed vapour: convection's own contribution plus the water
    # module's, which excludes convection's so it is not counted twice.
    own = sum(c.as_array() for c in contributions if c.state_id == VAPOUR_STATE)
    committed = vapour + moved + np.asarray(own) * DT
    condensing = expected > 0.0
    np.testing.assert_allclose(committed[condensing], capacity[condensing], rtol=1e-12)
    assert np.all(committed[~condensing] <= capacity[~condensing] * (1 + 1e-12) + 1e-300)

    # Water closes: what left the vapour is what precipitates.
    assert float(np.sum(vapour) - np.sum(committed)) == pytest.approx(
        water._precipitation * DT, rel=1e-12
    )
    # Its latent heat is released where it condensed, and nowhere else.
    latent = (
        water.latent_heat * expected / DT
        + deposition_latent_heat_excess(temperature) * expected / DT
    )
    np.testing.assert_allclose(water._latent_heating, latent, rtol=1e-12, atol=1e-15)
    assert np.all(water._latent_heating[~condensing] == 0.0)


@pytest.mark.req("CON-001", "WAT-002")
def test_condensing_after_convection_needs_convection_to_carry_vapour():
    with pytest.raises(ConfigurationError):
        ConvectiveAdjustmentModule(transport_vapour=False, condensation_follows=True)
    with pytest.raises(ConfigurationError):
        reference_column(condense_after_convection=True, vapour_transport=False)


@pytest.mark.req("CON-001", "WAT-002", "SYS-006")
def test_the_coupled_column_runs_convection_first_and_closes_its_budgets():
    """The option reverses the order and the conservation gates still close."""
    temperature, vapour = _unstable_moist_column()
    kernel = reference_column(
        condense_after_convection=True,
        n_layers=N,
        initial_air_temperature=temperature,
        initial_surface_temperature=300.0,
    )
    order = [s.module_id for s in kernel.schedule.ordered()]
    assert order.index("convection") < order.index("water")
    kernel.initial_state[VAPOUR_STATE] = vapour
    run = kernel.run(24)
    assert str(run.manifest.status) == "complete", run.manifest.invalidation_reasons
    assert run.worst_relative_residual() < 1e-11


@pytest.mark.req("CON-001", "SYS-006")
def test_the_tropical_column_keeps_the_order_the_column_chose():
    """The exchange is spliced in by splitting order, not by append order."""
    from mfespf.run.tropical_column import tropical_column

    on = [s.module_id for s in tropical_column(condense_after_convection=True).schedule.ordered()]
    off = [s.module_id for s in tropical_column().schedule.ordered()]
    assert on.index("tropical_exchange") < on.index("convection") < on.index("water")
    assert off.index("tropical_exchange") < off.index("water") < off.index("convection")
