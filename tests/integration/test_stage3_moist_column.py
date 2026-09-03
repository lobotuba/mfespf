"""Stage 3 of the validation ladder: thermodynamics, convection and water.

Charter section 19 stage 3; requirements ATM-002, WAT-001, WAT-002, CON-001,
CON-007, MOD-011, STA-003, STA-006.

Two things are new here and both are architectural as much as physical. Three
modules now contribute labeled tendencies to one state, and the framework moves
a conserved quantity that is not energy -- water -- with both budgets required to
close at once across a phase change.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.coupling.transfers import ConservedQuantity
from mfespf.physics.atmosphere import ATMOSPHERE_TEMPERATURE_STATE
from mfespf.physics.surface import SURFACE_TEMPERATURE_STATE, SURFACE_WATER_STATE
from mfespf.physics.thermodynamics import potential_temperature
from mfespf.physics.water import VAPOUR_STATE
from mfespf.run.column_configurations import (
    column_grid,
    moist_radiative_convective,
    radiative_convective,
    window_atmosphere,
)

pytestmark = pytest.mark.stage(3)

GRID = column_grid(16)
C_SMALL = 1.0e6
BASE = dict(grid=GRID, n_layers=8, areal_heat_capacity=C_SMALL, master_dt=3600.0)


def _layer_pressure(kernel) -> np.ndarray:
    atmosphere = next(s.module for s in kernel.schedule.ordered() if s.module_id == "atmosphere")
    return atmosphere.layer_pressure


# --- convection ------------------------------------------------------------


@pytest.mark.req("ATM-002", "CON-009", "KER-012")
def test_radiative_convective_equilibrium_is_stable_and_closes():
    """The Stage 3 configuration must reach a dry-stable profile with a closed budget."""
    kernel = radiative_convective(
        1.0e-4, initial_surface_temperature=270.0, initial_air_temperature=250.0, **BASE
    )
    result = kernel.run(6000)

    assert result.may_support_a_claim
    temperatures = result.steps[-1].outputs[ATMOSPHERE_TEMPERATURE_STATE].values
    theta = potential_temperature(temperatures, _layer_pressure(kernel))
    assert np.all(np.diff(theta) >= -1e-6)
    assert result.worst_relative_residual() < 1e-11


@pytest.mark.req("ATM-002", "CON-009")
def test_convection_adds_no_energy_to_the_column():
    """Mixing redistributes; only the sensible flux crosses a boundary.

    The convective contribution to every layer, weighted by heat capacity, must
    sum to zero. If it did not, the scheme would be a heat source disguised as a
    transport, and the column budget would close only because the error was
    hidden inside a term nobody separated out.
    """
    kernel = radiative_convective(
        1.0e-4, initial_surface_temperature=290.0, initial_air_temperature=230.0, **BASE
    )
    result = kernel.run(40)
    atmosphere = next(s.module for s in kernel.schedule.ordered() if s.module_id == "atmosphere")
    capacity = atmosphere.layer_heat_capacity

    saw_mixing = False
    for step in result.steps:
        budget = step.state_budgets[ATMOSPHERE_TEMPERATURE_STATE]
        mixing = [c for c in budget.contributions if c.label == "dry_convective_adjustment"]
        assert len(mixing) == 1
        energy = float(np.sum(capacity * mixing[0].as_array()))
        assert energy == pytest.approx(0.0, abs=1e-6)
        saw_mixing = saw_mixing or bool(np.any(np.abs(mixing[0].as_array()) > 1e-12))
    assert saw_mixing, "the adjustment never activated, so nothing was tested"


@pytest.mark.req("ATM-002", "MOD-011", "STA-003", "STA-006")
def test_three_modules_contribute_labeled_tendencies_to_one_state():
    """A1.3: several processes act, one owner commits, and the labels survive."""
    kernel = moist_radiative_convective(
        1.0e-4,
        initial_surface_temperature=295.0,
        initial_air_temperature=280.0,
        evaporation_coefficient=5.0e-3,
        **BASE,
    )
    result = kernel.run(20)

    for step in result.steps:
        budget = step.state_budgets[ATMOSPHERE_TEMPERATURE_STATE]
        providers = {c.provider_module_id for c in budget.contributions}
        labels = {c.label for c in budget.contributions}
        assert providers == {"radiation", "convection", "water"}
        assert labels == {
            "radiative_flux_divergence",
            "surface_sensible_heat",
            "dry_convective_adjustment",
            "latent_heat_release",
        }
        assert np.allclose(budget.residual(), 0.0, atol=1e-10)


@pytest.mark.req("ATM-002", "CLM-003")
def test_convection_changes_the_profile_it_is_given():
    """Convection must cool the surface and warm the air beneath it.

    The declared sign, before the run: carrying heat upward out of the surface
    lowers the surface temperature and raises the lowest layer's. The comparison
    is against the same optics with convection absent, so the difference is
    attributable to the process and not to a different atmosphere.
    """
    settings = dict(initial_surface_temperature=270.0, initial_air_temperature=250.0, **BASE)
    radiative = window_atmosphere(kappa_opaque=1.0e-4, **settings).run(12000)
    convective = radiative_convective(1.0e-4, **settings).run(12000)

    radiative_surface = radiative.steps[-1].outputs["surface_temperature"].scalar()
    convective_surface = convective.steps[-1].outputs["surface_temperature"].scalar()
    radiative_air = radiative.steps[-1].outputs[ATMOSPHERE_TEMPERATURE_STATE].values[0]
    convective_air = convective.steps[-1].outputs[ATMOSPHERE_TEMPERATURE_STATE].values[0]

    assert convective_surface < radiative_surface - 1.0
    assert convective_air > radiative_air + 1.0
    # And the surface-to-air jump the flux feeds on is much reduced.
    assert (convective_surface - convective_air) < (radiative_surface - radiative_air)


@pytest.mark.req("ATM-002", "VAL-001")
def test_sensible_heat_flux_vanishes_when_isothermal():
    """The null case for the bulk formula: no gradient, no flux."""
    kernel = radiative_convective(
        1.0e-4, initial_surface_temperature=260.0, initial_air_temperature=260.0, **BASE
    )
    coupler, _ = kernel.prepare(n_steps=2)
    coupler.step()
    assert float(coupler.bus["surface_sensible_heat_flux"].values) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.req("CON-002", "CON-004", "ATM-002")
def test_sensible_heat_is_transferred_not_created():
    """What the surface loses to the flux is exactly what the column gains."""
    kernel = radiative_convective(
        1.0e-4, initial_surface_temperature=290.0, initial_air_temperature=250.0, **BASE
    )
    result = kernel.run(20)

    for report in result.conservation_reports:
        matched = [
            r
            for r in report.coupling_residuals
            if r.sender == "surface" and r.receiver == "atmosphere"
        ]
        assert matched
        for residual in matched:
            assert residual.relative_mismatch == pytest.approx(0.0, abs=1e-12)
        assert report.passed, report.describe()


# --- water -----------------------------------------------------------------


@pytest.mark.req("WAT-001", "WAT-002", "CON-001", "CON-007")
def test_moist_column_budgets_close_for_water_and_energy():
    """Two conserved quantities, both closing, across a phase change.

    This is CON-007 in full: the matter transferred by evaporation and
    precipitation and the latent energy that travels with it must reconcile at
    the same time, between the same modules, over the same interval.
    """
    kernel = moist_radiative_convective(
        1.0e-4,
        initial_surface_temperature=300.0,
        initial_air_temperature=285.0,
        evaporation_coefficient=1.0e-2,
        **BASE,
    )
    result = kernel.run(300)

    assert result.may_support_a_claim
    for report in result.conservation_reports:
        assert report.passed, report.describe()
        quantities = {b.quantity for b in report.budgets}
        assert ConservedQuantity.WATER in quantities
        assert ConservedQuantity.ENERGY in quantities
        # Water is a closed loop between the surface and the column.
        assert report.system_residual(ConservedQuantity.WATER) == pytest.approx(0.0, abs=1e-9)


@pytest.mark.req("WAT-001", "CON-001")
def test_total_water_is_conserved_between_the_surface_and_the_column():
    """WAT-001 at the system level: the loop neither gains nor loses water."""
    kernel = moist_radiative_convective(
        1.0e-4,
        initial_surface_temperature=300.0,
        initial_air_temperature=285.0,
        evaporation_coefficient=1.0e-2,
        **BASE,
    )
    result = kernel.run(300)

    first, last = result.steps[0], result.steps[-1]
    initial_total = float(np.sum(first.state_budgets[VAPOUR_STATE].initial_value)) + float(
        first.state_budgets[SURFACE_WATER_STATE].initial_value
    )
    final_total = float(np.sum(last.outputs[VAPOUR_STATE].values)) + float(
        last.outputs[SURFACE_WATER_STATE].values
    )
    assert final_total == pytest.approx(initial_total, rel=1e-12)
    # And some water actually moved, or the conservation is vacuous.
    assert float(last.outputs["atmosphere_precipitable_water"].values) > 0.1


@pytest.mark.req("WAT-002", "CON-007")
def test_evaporation_conserves_water_and_carries_its_latent_energy():
    """Every kilogram leaving the surface takes L_v joules with it."""
    kernel = moist_radiative_convective(
        1.0e-4,
        initial_surface_temperature=300.0,
        initial_air_temperature=285.0,
        evaporation_coefficient=1.0e-2,
        **BASE,
    )
    result = kernel.run(20)
    water = next(s.module for s in kernel.schedule.ordered() if s.module_id == "water")

    for report in result.conservation_reports:
        water_transfer = next(
            r
            for r in report.coupling_residuals
            if r.sender == "surface"
            and r.receiver == "water"
            and r.quantity is ConservedQuantity.WATER
        )
        energy_transfer = next(
            r
            for r in report.coupling_residuals
            if r.sender == "surface"
            and r.receiver == "water"
            and r.quantity is ConservedQuantity.ENERGY
        )
        assert water_transfer.imported > 0.0
        assert energy_transfer.imported == pytest.approx(
            water.latent_heat * water_transfer.imported, rel=1e-12
        )


@pytest.mark.req("WAT-001", "VAL-001")
def test_evaporation_stops_when_the_air_is_saturated():
    """The limiting case of the bulk formula: no humidity deficit, no flux."""
    from mfespf.contracts.module import StepContext

    kernel = moist_radiative_convective(
        1.0e-4,
        initial_surface_temperature=290.0,
        initial_air_temperature=290.0,
        evaporation_coefficient=1.0e-2,
        **BASE,
    )
    coupler, _ = kernel.prepare(n_steps=4)
    water = next(s.module for s in kernel.schedule.ordered() if s.module_id == "water")

    coupler.step()
    assert float(coupler.bus["surface_evaporation_flux"].values) > 0.0

    # Saturate the lowest layer by hand and re-evaluate: the deficit is gone.
    from mfespf.physics.thermodynamics import saturation_specific_humidity

    atmosphere = next(s.module for s in kernel.schedule.ordered() if s.module_id == "atmosphere")
    saturated = saturation_specific_humidity(
        np.array([290.0]), np.array([atmosphere.layer_pressure[0]])
    )[0]
    water.vapour[0] = saturated * atmosphere.layer_mass[0] * 1.5

    water.begin_step(StepContext(time=0.0, dt=3600.0, step_index=0))
    water.receive_inputs(dict(coupler.bus))
    water.compute_tendencies()
    assert water._evaporation == pytest.approx(0.0, abs=1e-15)


@pytest.mark.req("WAT-001", "WAT-002", "CON-007")
def test_condensation_conserves_water_and_energy_together():
    """Vapour above saturation must leave as precipitation and release exactly L_v per kg."""
    kernel = moist_radiative_convective(
        1.0e-4,
        initial_surface_temperature=305.0,
        initial_air_temperature=250.0,
        evaporation_coefficient=5.0e-2,
        **BASE,
    )
    result = kernel.run(200)
    water = next(s.module for s in kernel.schedule.ordered() if s.module_id == "water")

    condensing = [s for s in result.steps if float(s.outputs["precipitation_flux"].values) > 0.0]
    assert condensing, "no condensation occurred, so the phase change was never tested"

    for step in condensing:
        precipitation = float(step.outputs["precipitation_flux"].values)
        latent = float(np.sum(step.outputs["latent_heating_rate"].values))
        assert latent == pytest.approx(water.latent_heat * precipitation, rel=1e-12)

    for report in result.conservation_reports:
        assert report.passed, report.describe()


@pytest.mark.req("VAL-004", "MOD-014")
def test_a_dry_column_with_no_evaporation_is_inert():
    """The water module's null control: switch the source off and nothing moves."""
    kernel = moist_radiative_convective(
        1.0e-4,
        initial_surface_temperature=290.0,
        initial_air_temperature=270.0,
        evaporation_coefficient=0.0,
        **BASE,
    )
    result = kernel.run(50)

    for step in result.steps:
        assert float(step.outputs["atmosphere_precipitable_water"].values) == 0.0
        assert float(step.outputs["surface_evaporation_flux"].values) == 0.0
        assert float(step.outputs["precipitation_flux"].values) == 0.0
    assert float(result.steps[-1].outputs[SURFACE_WATER_STATE].values) == pytest.approx(
        1.0e4, rel=1e-14
    )
    for report in result.conservation_reports:
        assert report.passed, report.describe()


@pytest.mark.req("STA-006")
def test_the_surface_budget_separates_radiative_from_turbulent_cooling():
    """Charter 9.9 wants the surface energy terms apart, not netted."""
    kernel = moist_radiative_convective(
        1.0e-4,
        initial_surface_temperature=300.0,
        initial_air_temperature=280.0,
        evaporation_coefficient=1.0e-2,
        **BASE,
    )
    result = kernel.run(20)

    for step in result.steps:
        budget = step.state_budgets[SURFACE_TEMPERATURE_STATE]
        labels = {c.label for c in budget.contributions}
        assert labels == {
            "net_radiative_flux",
            "evaporative_cooling",
            "surface_sensible_heat",
        }
        assert np.allclose(budget.residual(), 0.0, atol=1e-10)


@pytest.mark.req("MOD-008")
def test_moist_column_checkpoint_round_trips_both_stores():
    kernel = moist_radiative_convective(
        1.0e-4, evaporation_coefficient=1.0e-2, initial_surface_temperature=295.0, **BASE
    )
    coupler, _ = kernel.prepare(n_steps=20)
    coupler.run(10)

    water = next(s.module for s in kernel.schedule.ordered() if s.module_id == "water")
    surface = next(s.module for s in kernel.schedule.ordered() if s.module_id == "surface")
    saved_vapour = water.checkpoint()
    saved_surface = surface.checkpoint()

    coupler.run(5)
    assert not np.array_equal(water.vapour, saved_vapour[VAPOUR_STATE])

    water.restore(saved_vapour)
    surface.restore(saved_surface)
    assert np.array_equal(water.vapour, saved_vapour[VAPOUR_STATE])
    assert surface.water_store == saved_surface[SURFACE_WATER_STATE]
