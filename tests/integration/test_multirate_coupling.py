"""KER-011 and TIM-002 to TIM-004: multirate coupling.

Two modules advancing at different rates must still produce a budget that
closes, and the interpolation the schedule uses to bridge them must be declared
and reported rather than silently absorbed.
"""

from __future__ import annotations

import pytest

from mfespf.run.configurations import airless_radiating_body

pytestmark = pytest.mark.stage(1)

C_SMALL = 1.0e6


@pytest.mark.req("KER-011", "TIM-002", "TIM-003")
def test_modules_advance_at_different_rates_within_one_schedule():
    kernel = airless_radiating_body(
        areal_heat_capacity=C_SMALL, master_dt=3600.0, solar_step_multiple=24
    )
    schedule = kernel.schedule

    assert schedule.module_dt("surface") == 3600.0
    assert schedule.module_dt("solar") == 86400.0
    assert schedule.coupling_period == 24
    assert "synchronization every 24" in schedule.describe()


@pytest.mark.req("KER-011", "CON-004")
def test_the_budget_closes_across_a_multirate_coupling_interval():
    """The slow module writes one record where the fast one writes twenty-four.

    Reconciling per master step would compare a daily export against the first
    hourly import and report a mismatch that is an artifact of the schedule.
    Reconciling at the synchronization point compares like with like.
    """
    kernel = airless_radiating_body(
        areal_heat_capacity=C_SMALL, master_dt=3600.0, solar_step_multiple=24
    )
    result = kernel.run(96)

    reports = result.conservation_reports
    assert len(reports) == 4  # 96 steps / 24-step coupling period
    for report in reports:
        assert report.passed, report.describe()


@pytest.mark.req("TIM-003")
def test_only_scheduled_steps_run_the_slow_modules():
    kernel = airless_radiating_body(
        areal_heat_capacity=C_SMALL, master_dt=3600.0, solar_step_multiple=6
    )
    result = kernel.run(12)

    ran_solar = [i for i, s in enumerate(result.steps) if "solar" in s.ran]
    held_solar = [i for i, s in enumerate(result.steps) if "solar" in s.held]
    assert ran_solar == [0, 6]
    assert held_solar == [1, 2, 3, 4, 5, 7, 8, 9, 10, 11]
    # The surface never stalls: it advances on every master step.
    assert all("surface" in s.ran for s in result.steps)


@pytest.mark.req("TIM-004")
def test_the_zero_order_hold_is_reported_with_its_staleness():
    kernel = airless_radiating_body(
        areal_heat_capacity=C_SMALL, master_dt=3600.0, solar_step_multiple=6
    )
    result = kernel.run(12)

    holds = [w for w in result.warnings if "zero-order hold" in w]
    assert holds
    assert "up to 18000 s" in holds[0]  # five stale master steps of 3600 s
    assert "TIM-004" in holds[0]


@pytest.mark.req("CON-009")
def test_an_unreconciled_trailing_interval_is_reported_not_ignored():
    """Stopping mid-interval leaves transfers unaudited, and the run says so."""
    kernel = airless_radiating_body(
        areal_heat_capacity=C_SMALL, master_dt=3600.0, solar_step_multiple=24
    )
    result = kernel.run(30)  # 24 + 6: one closed interval, six unaudited steps

    assert len(result.conservation_reports) == 1
    assert any("unreconciled coupling interval" in w for w in result.warnings)


@pytest.mark.req("TIM-004", "KER-011")
def test_a_faster_exchange_reduces_the_interpolation_error():
    """The hold is a real approximation, and refining the exchange shrinks it.

    Under an eccentric orbit the incident flux varies through the year, so a
    daily exchange and an hourly one give measurably different answers. The
    difference is the coupling error TIM-004 requires be characterised rather
    than assumed negligible.
    """

    def final_temperature(step_multiple: int) -> float:
        kernel = airless_radiating_body(
            areal_heat_capacity=C_SMALL,
            master_dt=3600.0,
            solar_step_multiple=step_multiple,
            eccentricity=0.6,
            initial_temperature=250.0,
        )
        result = kernel.run(24 * 60)
        return result.steps[-1].outputs["surface_temperature"].scalar()

    hourly = final_temperature(1)
    six_hourly = final_temperature(6)
    daily = final_temperature(24)

    assert abs(six_hourly - hourly) < abs(daily - hourly)
