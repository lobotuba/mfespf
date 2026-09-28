"""Stage 1 of the validation ladder: the airless radiating body.

Charter section 19 stage 1, requirements KER-012 and KER-013.

Three things are checked, and they are independent of one another: the coupled
system reaches the analytically known equilibrium, it approaches it on the
timescale the reduced equation predicts, and its energy budget closes to
round-off at every coupling interval.  A model can pass any one of these while
failing the others -- an equilibrium reached with a leaking budget is the
classic case -- so all three are gates, not one.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from mfespf.core.constants import value_of
from mfespf.coupling.transfers import ConservedQuantity
from mfespf.physics.radiation import equilibrium_temperature, relaxation_timescale
from mfespf.physics.surface import SURFACE_TEMPERATURE_STATE
from mfespf.run.configurations import airless_radiating_body, null_control

pytestmark = pytest.mark.stage(1)

TSI = value_of("total_solar_irradiance_reference")
SMALL_HEAT_CAPACITY = 1.0e6  # a thin slab, so equilibrium is reached quickly


@pytest.mark.req("KER-012", "KER-013", "VAL-001")
def test_slab_reaches_analytic_equilibrium(report_bias):
    """The coupled run must land on the closed-form temperature, not near it.

    Run for roughly 27 e-folding times from well below equilibrium, so what is
    being tested is the fixed point the coupled system settles on rather than how
    far it happened to get. Forward Euler's fixed point is exact -- the scheme
    stops moving precisely where absorbed equals emitted -- so a discrepancy here
    would be a real error in the coupling, not truncation.
    """
    albedo, emissivity = 0.3, 1.0
    kernel = airless_radiating_body(
        albedo=albedo,
        emissivity=emissivity,
        areal_heat_capacity=SMALL_HEAT_CAPACITY,
        initial_temperature=100.0,
        master_dt=600.0,
        solar_step_multiple=1,
    )
    result = kernel.run(12000)

    expected = equilibrium_temperature(TSI, albedo, emissivity)
    final = result.steps[-1].outputs["surface_temperature"].scalar()

    assert result.may_support_a_claim
    assert final == pytest.approx(expected, rel=1e-9)

    # Reported for the same reason as everywhere else, and it is the useful end
    # of the scale: this column has no atmosphere, so its 33 K deficit against
    # the observed surface temperature is the greenhouse effect, measured.
    bias = report_bias(result, label="Stage 1 airless body (no atmosphere)")
    assert bias.equilibrated
    assert bias.error_k == pytest.approx(expected - 288.0, abs=1e-6)


@pytest.mark.req("KER-012", "TIM-007", "TIM-008")
def test_slab_relaxation_timescale_matches_the_linearised_prediction():
    """A small perturbation must decay with tau = C / (4 eps sigma T_eq^3).

    Started 1 K below equilibrium so the linearisation is valid; the measured
    e-folding is compared with the analytic value the reduced equation predicts.
    """
    albedo, emissivity, C = 0.3, 1.0, SMALL_HEAT_CAPACITY
    equilibrium = equilibrium_temperature(TSI, albedo, emissivity)
    tau = relaxation_timescale(C, equilibrium, emissivity)

    kernel = airless_radiating_body(
        albedo=albedo,
        emissivity=emissivity,
        areal_heat_capacity=C,
        initial_temperature=equilibrium - 1.0,
        master_dt=tau / 400.0,
        solar_step_multiple=1,
    )
    result = kernel.run(400)

    times = np.array([s.time for s in result.steps])
    temperatures = np.array([s.outputs["surface_temperature"].scalar() for s in result.steps])
    departure = equilibrium - temperatures
    usable = departure > departure[0] * 1e-3
    slope = np.polyfit(times[usable], np.log(departure[usable]), 1)[0]
    measured_tau = -1.0 / slope

    assert measured_tau == pytest.approx(tau, rel=0.02)


@pytest.mark.req("CON-009", "CON-004", "KER-013", "ACC-004")
def test_stage1_energy_budget_closes_at_every_coupling_interval():
    kernel = airless_radiating_body(
        areal_heat_capacity=SMALL_HEAT_CAPACITY, master_dt=600.0, solar_step_multiple=6
    )
    result = kernel.run(120)

    reports = result.conservation_reports
    assert reports, "the run produced no reconciled coupling interval"
    for report in reports:
        assert report.passed, report.describe()
        assert report.worst_relative_residual() < 1e-12
        # Every module's budget is published, not just the failing ones.
        assert {b.module_id for b in report.budgets} == {"orbit", "surface", "space"}


@pytest.mark.req("CON-001", "CON-009")
def test_system_energy_is_conserved_across_the_boundary():
    """Storage gained by the surface and space must equal what the Sun supplied."""
    kernel = airless_radiating_body(
        areal_heat_capacity=SMALL_HEAT_CAPACITY, master_dt=600.0, solar_step_multiple=1
    )
    result = kernel.run(50)

    supplied = 0.0
    stored = 0.0
    for report in result.conservation_reports:
        for budget in report.budgets_for(ConservedQuantity.ENERGY):
            if budget.module_id == "orbit":
                supplied += budget.net_physical
            else:
                stored += budget.stored_change

    assert stored == pytest.approx(supplied, rel=1e-12)


@pytest.mark.req("VAL-004", "MOD-014", "KER-012")
def test_zero_insolation_null_control():
    """With the star switched off the surface may only cool.

    The disabled-coupling control MOD-014 requires: if this run warms, the
    coupling is injecting energy that no physics accounts for.
    """
    kernel = null_control(
        areal_heat_capacity=SMALL_HEAT_CAPACITY, initial_temperature=250.0, master_dt=600.0
    )
    result = kernel.run(120)

    temperatures = [s.outputs["surface_temperature"].scalar() for s in result.steps]
    assert temperatures[-1] < temperatures[0]
    assert all(later <= earlier for earlier, later in itertools.pairwise(temperatures))
    for report in result.conservation_reports:
        assert report.passed, report.describe()


@pytest.mark.req("STA-006", "STA-001")
def test_owner_commits_state_from_labeled_contributions_only():
    kernel = airless_radiating_body(
        areal_heat_capacity=SMALL_HEAT_CAPACITY, master_dt=600.0, solar_step_multiple=1
    )
    result = kernel.run(10)

    for step in result.steps:
        budget = step.state_budgets[SURFACE_TEMPERATURE_STATE]
        assert float(budget.residual()) == pytest.approx(0.0, abs=1e-12)
        assert [c.label for c in budget.contributions] == ["net_radiative_flux"]


@pytest.mark.req("MOD-008")
def test_checkpoint_and_restore_round_trip_exactly():
    kernel = airless_radiating_body(
        areal_heat_capacity=SMALL_HEAT_CAPACITY, master_dt=600.0, solar_step_multiple=1
    )
    coupler, _ = kernel.prepare(n_steps=20)
    coupler.run(10)

    surface = next(s.module for s in kernel.schedule.ordered() if s.module_id == "surface")
    saved = surface.checkpoint()
    coupler.run(5)
    assert surface.temperature != saved[SURFACE_TEMPERATURE_STATE]

    surface.restore(saved)
    assert surface.temperature == saved[SURFACE_TEMPERATURE_STATE]
