"""Stage 5: the forcing protocols TIM-007 asks for.

Every test here is about the same property, approached four ways: the energy a
forcing delivers is a declaration, not an artifact of the time step. A forcing
sampled at the start of each step would deliver a different amount of energy at
every resolution, and a convergence study run against it would be measuring the
forcing rather than the solver.
"""

from __future__ import annotations

import math

import pytest

from mfespf.core.errors import ConfigurationError
from mfespf.experiments.forcing import (
    ImpulseForcing,
    NoForcing,
    PeriodicForcing,
    PulseForcing,
    StepForcing,
)

pytestmark = pytest.mark.stage(5)

DAY = 86400.0


def _delivered(forcing, span: float, dt: float) -> float:
    """Total energy the forcing delivers over ``span`` in steps of ``dt``."""
    n = round(span / dt)
    return sum(forcing.energy_delivered(i * dt, (i + 1) * dt) for i in range(n))


# --- the null control ------------------------------------------------------


@pytest.mark.req("TIM-007", "VAL-001")
def test_no_forcing_delivers_exactly_zero():
    """The null case, run through the same path as every other forcing."""
    forcing = NoForcing()
    assert forcing.regime == "none"
    assert forcing.mean_amplitude(0.0, 1.0e6) == 0.0
    assert _delivered(forcing, 100 * DAY, DAY) == 0.0


# --- step ------------------------------------------------------------------


@pytest.mark.req("TIM-007")
def test_step_is_off_before_its_onset_and_on_after():
    forcing = StepForcing(amplitude=4.0, onset=10 * DAY)
    assert forcing.mean_amplitude(0.0, 10 * DAY) == 0.0
    assert forcing.mean_amplitude(10 * DAY, 11 * DAY) == pytest.approx(4.0)
    # A step whose onset falls inside a step is applied for the part of the
    # interval it covers, not for all of it or none of it.
    assert forcing.mean_amplitude(9.5 * DAY, 10.5 * DAY) == pytest.approx(2.0)


@pytest.mark.req("TIM-007", "CMP-009")
def test_step_energy_is_independent_of_the_time_step():
    forcing = StepForcing(amplitude=4.0, onset=3.5 * DAY)
    reference = _delivered(forcing, 100 * DAY, DAY)
    for dt in (0.25 * DAY, 2 * DAY, 5 * DAY):
        assert _delivered(forcing, 100 * DAY, dt) == pytest.approx(reference, rel=1e-12)
    assert reference == pytest.approx(4.0 * 96.5 * DAY, rel=1e-12)


# --- pulse -----------------------------------------------------------------


@pytest.mark.req("TIM-007")
def test_pulse_delivers_its_declared_energy_whatever_the_step():
    """Both edges land inside a step at two of these resolutions."""
    forcing = PulseForcing(amplitude=8.0, onset=2.5 * DAY, duration=7.5 * DAY)
    declared = 8.0 * 7.5 * DAY
    assert forcing.parameters()["energy_J_m2"] == pytest.approx(declared)
    for dt in (0.5 * DAY, DAY, 2 * DAY, 3 * DAY):
        assert _delivered(forcing, 60 * DAY, dt) == pytest.approx(declared, rel=1e-12)


@pytest.mark.req("TIM-007", "SYS-006")
def test_a_pulse_must_declare_a_duration():
    """A zero-width pulse is an impulse, and the two are not the same object."""
    with pytest.raises(ConfigurationError, match="ImpulseForcing"):
        PulseForcing(amplitude=1.0, duration=0.0)
    with pytest.raises(ConfigurationError, match="onset"):
        PulseForcing(amplitude=1.0, duration=DAY, onset=-1.0)


# --- impulse ---------------------------------------------------------------


@pytest.mark.req("TIM-007", "CMP-009")
def test_impulse_energy_is_exact_at_every_step_and_its_flux_is_not():
    """An impulse fixes the energy; the flux is whatever the step makes it.

    This is the honest form of a delta function in discrete time. Fixing the
    flux instead would make the delivered energy proportional to the step, so an
    impulse response would change with resolution for reasons that have nothing
    to do with the physics.
    """
    energy = 5.0e7
    forcing = ImpulseForcing(energy=energy, onset=0.0)
    for dt in (0.5 * DAY, DAY, 4 * DAY):
        assert _delivered(forcing, 40 * DAY, dt) == pytest.approx(energy, rel=1e-12)
        # The flux does depend on the step, which is the part that is not
        # resolved and is declared as such rather than hidden.
        assert forcing.mean_amplitude(0.0, dt) == pytest.approx(energy / dt, rel=1e-12)
    assert "not resolved below the step" in forcing.parameters()["shape"]


@pytest.mark.req("TIM-007")
def test_impulse_has_no_instantaneous_value():
    """Asking for one would invent a flux from a step the caller did not give."""
    assert ImpulseForcing(energy=1.0e7).mean_amplitude(0.0, 0.0) == 0.0


# --- periodic --------------------------------------------------------------


@pytest.mark.req("TIM-007")
def test_periodic_forcing_integrates_to_zero_over_whole_periods():
    """A periodic experiment must not warm the system by accident.

    True at any step size because the mean is the exact integral rather than a
    midpoint sample: at a tenth of a period the two differ by about 0.4 percent
    of the amplitude, which over many cycles is a drift, not a rounding error.
    """
    period = 30 * DAY
    forcing = PeriodicForcing(amplitude=6.0, period=period)
    for dt in (period / 4.0, period / 7.0, period / 40.0):
        assert _delivered(forcing, 12 * period, dt) == pytest.approx(0.0, abs=1e-6)


@pytest.mark.req("TIM-007")
def test_periodic_mean_over_a_step_is_the_exact_integral():
    period = 30 * DAY
    forcing = PeriodicForcing(amplitude=6.0, period=period)
    omega = 2.0 * math.pi / period
    start, end = 3 * DAY, 8 * DAY
    expected = (6.0 / omega) * (math.cos(omega * start) - math.cos(omega * end)) / (end - start)
    assert forcing.mean_amplitude(start, end) == pytest.approx(expected, rel=1e-14)
    # And the midpoint sample it deliberately is not.
    assert forcing.mean_amplitude(start, end) != pytest.approx(
        6.0 * math.sin(omega * 0.5 * (start + end)), rel=1e-6
    )


@pytest.mark.req("TIM-007", "SYS-006")
def test_unphysical_forcing_declarations_are_refused():
    with pytest.raises(ConfigurationError, match="positive period"):
        PeriodicForcing(amplitude=1.0, period=0.0)
    with pytest.raises(ConfigurationError, match="onset"):
        StepForcing(amplitude=1.0, onset=-DAY)
    with pytest.raises(ConfigurationError, match="ends before it starts"):
        StepForcing(amplitude=1.0).mean_amplitude(10.0, 5.0)


@pytest.mark.req("TIM-007", "AUD-002")
def test_every_forcing_declares_itself_for_the_manifest():
    """The forcing enters the configuration hash, so it must be serialisable."""
    for forcing in (
        NoForcing(),
        StepForcing(amplitude=4.0),
        PulseForcing(amplitude=4.0, duration=DAY),
        ImpulseForcing(energy=1.0e7),
        PeriodicForcing(amplitude=4.0, period=DAY),
    ):
        parameters = forcing.parameters()
        assert parameters["regime"] == forcing.regime
        assert forcing.describe()
