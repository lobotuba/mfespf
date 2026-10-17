"""A surface held at a prescribed temperature (coupled-ascent-subsidence.md).

The energy that holds it crosses the system boundary as a declared source in the
surface's own budget, so the run still closes while its temperature never moves.
"""

from __future__ import annotations

import pytest

from mfespf.core.errors import ConfigurationError
from mfespf.physics.column_surface import FIXED_SURFACE_TRANSFER, ColumnSurfaceModule
from mfespf.physics.surface import SURFACE_TEMPERATURE_STATE
from mfespf.run.column_configurations import column_grid, moist_radiative_convective

pytestmark = pytest.mark.stage(3)

GRID = column_grid(16)
BASE = dict(grid=GRID, n_layers=8, areal_heat_capacity=1.0e6, master_dt=3600.0)


@pytest.mark.req("CON-001", "CON-009")
def test_a_held_surface_stays_put_and_declares_the_energy_that_holds_it():
    """Every step the surface ends exactly at its prescribed temperature; one
    labelled contribution cancels all the others; the power it represents is
    published and enters the surface's energy budget as the declared boundary
    source, and the whole system still reconciles."""
    kernel = moist_radiative_convective(
        1.0e-4,
        prescribed_surface_temperature=295.0,
        initial_air_temperature=280.0,
        evaporation_coefficient=1.0e-2,
        **BASE,
    )
    result = kernel.run(120)
    assert result.may_support_a_claim
    for report in result.conservation_reports:
        assert report.passed, report.describe()
    held_any = False
    for step, report in zip(result.steps, result.conservation_reports, strict=True):
        budget = step.state_budgets[SURFACE_TEMPERATURE_STATE]
        assert float(budget.final_value) == 295.0
        hold = [c for c in budget.contributions if c.label == FIXED_SURFACE_TRANSFER]
        others = [c for c in budget.contributions if c.label != FIXED_SURFACE_TRANSFER]
        assert len(hold) == 1
        assert float(hold[0].value) == -sum(float(c.value) for c in others)
        power = float(step.outputs[FIXED_SURFACE_TRANSFER].values)
        assert power == pytest.approx(float(hold[0].value) * 1.0e6, rel=1e-12, abs=1e-12)
        surface = next(
            b for b in report.budgets if b.module_id == "surface" and b.quantity.value == "energy"
        )
        assert surface.source_labels[FIXED_SURFACE_TRANSFER] == pytest.approx(
            power * 3600.0, rel=1e-12, abs=1e-9
        )
        held_any = held_any or abs(power) > 1.0
    assert held_any, "the surface was never pushed away from its held temperature"


@pytest.mark.req("CON-009")
def test_an_unheld_surface_is_unchanged_and_a_held_one_starts_where_it_is_held():
    """Without the option nothing new is published or declared; with it, a start
    elsewhere than the held temperature is refused."""
    kernel = moist_radiative_convective(
        1.0e-4,
        initial_surface_temperature=295.0,
        initial_air_temperature=280.0,
        evaporation_coefficient=1.0e-2,
        **BASE,
    )
    result = kernel.run(3)
    assert FIXED_SURFACE_TRANSFER not in result.steps[-1].outputs
    surface = ColumnSurfaceModule(
        grid=GRID, initial_temperature=290.0, prescribed_temperature=290.0
    )
    with pytest.raises(ConfigurationError, match="start it where it is held"):
        surface.initialize({SURFACE_TEMPERATURE_STATE: 280.0}, {})
