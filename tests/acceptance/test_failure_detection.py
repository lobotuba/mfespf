"""ACC-006 and ACC-007: the framework must detect its own failures.

These are the acceptance criteria that distinguish a framework from a model.  A
model that only reports success when it succeeds is untested; what has to be
demonstrated is that a deliberately broken mechanism is caught and correctly
classified, and that a numerically invalid run is labeled as such rather than
quietly reported.
"""

from __future__ import annotations

from typing import Any

import pytest

from mfespf.contracts.module import ModuleStatus
from mfespf.core.errors import ConservationError, OwnershipError
from mfespf.coupling.transfers import ConservedQuantity, Transfer
from mfespf.evaluation.convergence import time_step_convergence
from mfespf.physics.surface import SURFACE_TEMPERATURE_STATE, SlabSurfaceModule
from mfespf.registries.states import TendencyContribution
from mfespf.run.configurations import airless_radiating_body
from mfespf.run.manifest import RunStatus

pytestmark = pytest.mark.stage(0)

C_SMALL = 1.0e6


class LeakySurface(SlabSurfaceModule):
    """A surface that heats itself without accounting for where the energy came from.

    This is the archetypal undetected error: the model warms, the temperature
    looks plausible, and nothing in the physics says anything is wrong. The
    conservation audit is what catches it.
    """

    def compute_tendencies(self) -> list[TendencyContribution]:
        contributions = super().compute_tendencies()
        contributions.append(
            TendencyContribution(
                state_id=SURFACE_TEMPERATURE_STATE,
                provider_module_id=self.module_id,
                label="undeclared_heating",
                value=1.0e-5,
                units="K s^-1",
            )
        )
        return contributions


class UnownedWriter(SlabSurfaceModule):
    """A module that tries to contribute to a state it is not authorised for."""

    def compute_tendencies(self) -> list[TendencyContribution]:
        return [
            TendencyContribution(
                state_id="space_cumulative_energy_received",
                provider_module_id=self.module_id,
                label="reaching_into_another_module",
                value=1.0,
                units="W m^-2",
            )
        ]


class DoubleCountingSurface(SlabSurfaceModule):
    """A surface that reports its emitted energy twice over the same interval."""

    def report_transfers(self) -> tuple[list[Transfer], list[Transfer]]:
        exports, imports = super().report_transfers()
        start, end = self._interval
        exports.append(
            Transfer(
                transfer_id=f"emitted_longwave_again@{start:g}",
                sender_module=self.module_id,
                receiver_module=self.space_module_id,
                quantity=ConservedQuantity.ENERGY,
                amount=self._emitted * (end - start),
                units="J m^-2",
                start_time=start,
                end_time=end,
                channel="emitted_longwave",
            )
        )
        return exports, imports


def _kernel_with(surface: SlabSurfaceModule, **kwargs: Any):
    kernel = airless_radiating_body(areal_heat_capacity=C_SMALL, master_dt=600.0, **kwargs)
    for schedule in kernel.schedule.schedules:
        if schedule.module_id == "surface":
            schedule.module = surface
    return kernel


@pytest.mark.req("ACC-006", "CON-009", "CON-010")
def test_an_unaccounted_energy_source_is_detected_and_not_reported_as_support():
    """ACC-006: an intentional mechanism failure, correctly detected."""
    kernel = _kernel_with(
        LeakySurface(areal_heat_capacity=C_SMALL, initial_temperature=250.0),
        solar_step_multiple=1,
    )
    result = kernel.run(10, assert_conservation=False)

    assert not result.may_support_a_claim
    assert result.manifest.status in (
        RunStatus.FAILED,
        RunStatus.INVALID_FOR_SCIENTIFIC_INFERENCE,
    )
    report = result.conservation_reports[0]
    assert not report.passed
    assert any("surface" in f and "residual" in f for f in report.failures())


@pytest.mark.req("ACC-006", "SYS-006")
def test_the_leaky_run_stops_when_conservation_is_enforced():
    kernel = _kernel_with(
        LeakySurface(areal_heat_capacity=C_SMALL, initial_temperature=250.0),
        solar_step_multiple=1,
    )
    coupler, _ = kernel.prepare(n_steps=10)
    with pytest.raises(ConservationError, match="reconciliation failed"):
        coupler.run(5)


@pytest.mark.req("MOD-010", "STA-003", "ACC-002")
def test_writing_another_modules_state_is_refused_at_run_time():
    kernel = _kernel_with(UnownedWriter(areal_heat_capacity=C_SMALL), solar_step_multiple=1)
    coupler, _ = kernel.prepare(n_steps=2)
    with pytest.raises(OwnershipError, match="not an authorized tendency provider"):
        coupler.step()


@pytest.mark.req("CON-006", "ACC-006")
def test_a_double_counted_transfer_is_detected():
    kernel = _kernel_with(
        DoubleCountingSurface(areal_heat_capacity=C_SMALL, initial_temperature=250.0),
        solar_step_multiple=1,
    )
    result = kernel.run(4, assert_conservation=False)
    report = result.conservation_reports[0]
    assert report.double_counts
    assert not report.passed
    assert not result.may_support_a_claim


@pytest.mark.req("ACC-007", "SYS-008", "AUD-007")
def test_a_numerically_unstable_run_is_labeled_invalid_for_inference():
    """ACC-007: forward Euler beyond its stability limit produces nonsense.

    The slab's declared stability limit is dt < 2C / (4 eps sigma T^3). A day-long
    step on a slab whose relaxation time is minutes overshoots so far that the
    temperature goes negative on the first step. The framework refuses the
    impossible state (charter 15.2) rather than integrating on from it, and the
    run is not admissible as evidence.
    """
    kernel = airless_radiating_body(
        areal_heat_capacity=1.0e4,  # very thin slab: tau of order minutes
        master_dt=86400.0,  # a day per step, far beyond the limit
        solar_step_multiple=1,
        initial_temperature=250.0,
    )
    result = kernel.run(20, assert_conservation=False)

    assert result.manifest.status == RunStatus.FAILED
    assert not result.may_support_a_claim
    assert any("unphysical" in w for w in result.manifest.warnings)

    # The steps taken before the refusal are kept rather than discarded. They
    # are not evidence -- the status above is what a consumer checks, and it
    # says so -- but a failed run whose trajectory has been thrown away cannot
    # be diagnosed, and charter 3.5 asks for negative results to be preserved.
    assert result.steps
    assert len(result.steps) < 20
    assert result.steps[0].outputs["surface_temperature"].scalar() > 0.0


def _final_temperature_at(dt: float, *, end_time: float, heat_capacity: float) -> float:
    kernel = airless_radiating_body(
        areal_heat_capacity=heat_capacity,
        master_dt=dt,
        solar_step_multiple=1,
        initial_temperature=250.0,
    )
    result = kernel.run(round(end_time / dt), assert_conservation=False)
    if not result.steps:
        return float("nan")
    return result.steps[-1].outputs["surface_temperature"].scalar()


@pytest.mark.req("ACC-007", "VAL-001", "SYS-008")
def test_a_marginally_unstable_run_fails_the_time_step_convergence_gate():
    """The harder failure: finite, plausible numbers that mean nothing.

    Just past the forward-Euler stability limit the temperature oscillates with
    growing amplitude while every individual value stays physically reasonable.
    Nothing about a single trajectory reveals this. Refining the time step does:
    the answer does not settle, so the result is a property of the solver rather
    than of the physics (charter 15.2).
    """
    end_time = 8.0e6
    report = time_step_convergence(
        lambda dt: _final_temperature_at(dt, end_time=end_time, heat_capacity=1.0e6),
        base_dt=7.0e5,  # the slab's stability limit is about 5.3e5 s
        refinements=3,
        tolerance=1e-4,
        scale=250.0,
    )
    assert not report.converged
    assert (
        "not converging" in report.diagnosis or "above the declared tolerance" in report.diagnosis
    )


@pytest.mark.req("VAL-001", "UNC-003")
def test_a_stable_configuration_passes_the_convergence_gate_at_first_order():
    """The same gate, applied to a configuration inside its validity domain.

    Forward Euler is first order, and the study must measure that. A gate that
    only ever fails is not evidence of anything.
    """
    # Sampled mid-transient, at roughly one e-folding time. At the fixed point
    # forward Euler is exact, so an equilibrium sample would measure round-off
    # rather than the order of the scheme.
    end_time = 2.7e5
    report = time_step_convergence(
        lambda dt: _final_temperature_at(dt, end_time=end_time, heat_capacity=1.0e6),
        base_dt=1.35e4,
        refinements=4,
        tolerance=1e-3,
        expected_order=1.0,
        scale=250.0,
    )
    assert report.converged, report.describe()
    assert report.observed_order == pytest.approx(1.0, abs=0.25)
    assert report.richardson_estimate is not None


@pytest.mark.req("SYS-004")
def test_stub_output_is_marked_and_cannot_support_a_claim():
    from mfespf.physics.solar import ConstantSpectrumSolarSource, broadband_grid

    stub = ConstantSpectrumSolarSource(grid=broadband_grid(), status=ModuleStatus.STUB)
    assert not stub.status.may_support_a_claim
    stub.initialize({}, {})
    from mfespf.contracts.module import StepContext

    stub.begin_step(StepContext(time=0.0, dt=1.0, step_index=0))
    outputs = stub.emit_outputs()
    assert all(field.stub for field in outputs.values())
    assert "STUB - NOT PHYSICAL" in outputs["solar_total_irradiance_1au"].describe()
