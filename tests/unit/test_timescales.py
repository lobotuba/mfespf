"""Stage 5: measuring timescales, and refusing to measure what is not there.

Charter 11.1, 11.3 and requirement TIM-008.

The analysis is tested against series whose answers are known exactly -- one
exponential, two exponentials, one sinusoid -- rather than against model output,
so a failure here is a defect in the measurement and not in the physics. The
coupled runs in the integration suite then use these on output whose answer is
not known in advance.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError, InvalidForScientificInferenceError
from mfespf.evaluation.timescales import (
    TimescaleMeasurement,
    TimescaleReport,
    efolding_time,
    fit_exponential_modes,
    lag_correlation,
    periodic_response,
)

pytestmark = pytest.mark.stage(5)


def _measurement(tau: float, mechanism: str = "m", **kwargs) -> TimescaleMeasurement:
    return TimescaleMeasurement(
        observable="surface_temperature",
        subsystem="slab",
        mechanism=mechanism,
        forcing_regime="step",
        timescale_s=tau,
        method="test",
        **kwargs,
    )


# --- measurement without a model -------------------------------------------


@pytest.mark.req("TIM-008", "VAL-002")
def test_efolding_time_recovers_a_single_exponential_exactly():
    """The definition and the answer coincide for the one case they must."""
    tau = 3.7e6
    t = np.linspace(0.0, 20.0 * tau, 20001)
    y = 5.0 * (1.0 - np.exp(-t / tau))
    assert efolding_time(t, y) == pytest.approx(tau, rel=1e-5)


@pytest.mark.req("TIM-008")
def test_efolding_time_assumes_nothing_about_the_shape():
    """A two-mode response still has a well-defined 1 - 1/e crossing.

    It is not "the" timescale of that response and the report must not present it
    as one, but it is a measured property of the trajectory, and unlike a fit it
    cannot be wrong about a model it never assumed.
    """
    t = np.linspace(0.0, 4.0e8, 4001)
    y = 1.0 - 0.5 * np.exp(-t / 1.0e7) - 0.5 * np.exp(-t / 1.0e8)
    crossing = efolding_time(t, y, final=1.0)
    assert 1.0e7 < crossing < 1.0e8


@pytest.mark.req("TIM-008", "SYS-006")
def test_a_response_that_does_not_move_has_no_timescale():
    t = np.linspace(0.0, 100.0, 101)
    with pytest.raises(ConfigurationError, match="does not change"):
        efolding_time(t, np.full_like(t, 3.0))
    with pytest.raises(ConfigurationError, match="never reaches"):
        efolding_time(t, 1.0 - np.exp(-t / 1.0e6), final=1.0)


# --- measurement under an assumed model ------------------------------------


@pytest.mark.req("TIM-008", "VAL-002")
def test_a_one_mode_fit_recovers_a_known_single_exponential():
    tau, asymptote, amplitude = 4.2e6, 288.0, -3.5
    t = np.linspace(0.0, 8.0 * tau, 4001)
    y = asymptote + amplitude * np.exp(-t / tau)
    fit = fit_exponential_modes(t, y, n_modes=1)
    assert fit.timescales[0] == pytest.approx(tau, rel=2e-3)
    assert fit.asymptote == pytest.approx(asymptote, rel=1e-6)
    assert fit.amplitudes[0] == pytest.approx(amplitude, rel=5e-3)
    assert fit.relative_residual < 1e-4


@pytest.mark.req("TIM-008", "VAL-002")
def test_a_two_mode_fit_recovers_two_well_separated_timescales():
    fast, slow = 2.0e6, 6.0e7
    t = np.linspace(0.0, 4.0 * slow, 6001)
    y = 290.0 - 2.0 * np.exp(-t / fast) - 1.0 * np.exp(-t / slow)
    fit = fit_exponential_modes(t, y, n_modes=2)
    assert fit.timescales[0] == pytest.approx(fast, rel=5e-3)
    assert fit.timescales[1] == pytest.approx(slow, rel=5e-3)
    assert fit.amplitudes[0] == pytest.approx(-2.0, rel=1e-2)
    assert fit.amplitudes[1] == pytest.approx(-1.0, rel=1e-2)
    # Not zero: the search quantises each timescale to about a part in a
    # thousand, and that shows up here rather than being hidden.
    assert fit.relative_residual < 5e-4


@pytest.mark.req("TIM-008", "HYP-006")
def test_the_residual_is_what_says_whether_a_second_mode_is_there():
    """Granting a fit a second mode always lowers its residual; the question is
    by how much.

    On a single exponential a two-mode fit buys almost nothing and returns two
    near-degenerate timescales. On a genuine two-mode series it buys orders of
    magnitude. Reading a second timescale off a fit without checking this is how
    a model with one storage term comes to be described as having two.
    """
    t = np.linspace(0.0, 3.0e7, 3001)
    single = 280.0 - 4.0 * np.exp(-t / 4.0e6)
    double = 280.0 - 2.0 * np.exp(-t / 8.0e5) - 2.0 * np.exp(-t / 1.5e7)

    one_single = fit_exponential_modes(t, single, n_modes=1).relative_residual
    two_single = fit_exponential_modes(t, single, n_modes=2)
    one_double = fit_exponential_modes(t, double, n_modes=1).relative_residual
    two_double = fit_exponential_modes(t, double, n_modes=2).relative_residual

    assert one_single < 1e-4  # one mode already describes it
    assert one_double > 100.0 * one_single  # one mode does not describe this
    assert two_double < one_double / 100.0  # two modes do
    # The spurious pair is degenerate, which is the other signature.
    assert two_single.timescales[1] / two_single.timescales[0] < 3.0


@pytest.mark.req("TIM-008", "TIM-006")
def test_a_timescale_longer_than_the_record_is_marked_as_extrapolated():
    """A fit will return a number the data cannot support, and say so.

    The number is not useless -- it is the best estimate under the assumed model
    -- but it was inferred from the first part of a curve the run never
    finished, and reporting it in the language of a measurement is the confusion
    TIM-006 forbids between operating modes.
    """
    fast, slow = 1.0e6, 5.0e8
    t = np.linspace(0.0, 4.0e7, 2001)  # a twelfth of the slow mode
    y = 280.0 - 1.0 * np.exp(-t / fast) - 1.0 * np.exp(-t / slow)
    fit = fit_exponential_modes(t, y, n_modes=2)
    assert fit.constrained[0] is True
    assert fit.constrained[1] is False
    assert "extrapolated" in fit.describe()
    assert "EXTRAPOLATED" in _measurement(fit.timescales[1], constrained=False).describe()


@pytest.mark.req("TIM-008", "SYS-006")
def test_a_fit_refuses_input_it_cannot_fit():
    t = np.linspace(0.0, 1.0e6, 50)
    with pytest.raises(ConfigurationError, match="constant"):
        fit_exponential_modes(t, np.full_like(t, 2.0), n_modes=1)
    with pytest.raises(ConfigurationError, match="at least"):
        fit_exponential_modes(t[:4], np.arange(4.0), n_modes=2)
    with pytest.raises(ConfigurationError, match="at least one mode"):
        fit_exponential_modes(t, np.exp(-t / 1e5), n_modes=0)


@pytest.mark.req("TIM-008", "CMP-009")
def test_the_fit_is_deterministic():
    """No random start, so the same series gives the same answer every time."""
    t = np.linspace(0.0, 1.0e8, 2001)
    y = 285.0 - 3.0 * np.exp(-t / 4.0e6) - 1.5 * np.exp(-t / 5.0e7)
    first = fit_exponential_modes(t, y, n_modes=2)
    second = fit_exponential_modes(t, y, n_modes=2)
    assert first.timescales == second.timescales
    assert first.amplitudes == second.amplitudes


# --- frequency response ----------------------------------------------------


@pytest.mark.req("TIM-008", "VAL-002")
def test_periodic_response_recovers_a_known_gain_and_phase():
    period = 1.0e7
    omega = 2.0 * math.pi / period
    t = np.linspace(0.0, 10.0 * period, 4001)
    forcing = 4.0 * np.sin(omega * t)
    lag = 0.9
    response = 280.0 + 0.25 * 4.0 * np.sin(omega * t - lag)
    measured = periodic_response(t, forcing, response, period=period)
    assert measured.gain == pytest.approx(0.25, rel=1e-6)
    assert measured.phase_lag_rad == pytest.approx(lag, rel=1e-6)
    assert measured.phase_lag_s == pytest.approx(lag * period / (2.0 * math.pi), rel=1e-5)


@pytest.mark.req("TIM-008", "SYS-006")
def test_periodic_response_refuses_a_record_shorter_than_a_period():
    period = 1.0e7
    t = np.linspace(0.0, 0.8 * period, 200)
    signal = np.sin(2.0 * math.pi * t / period)
    with pytest.raises(ConfigurationError, match="less than one"):
        periodic_response(t, signal, signal, period=period, discard_fraction=0.0)


@pytest.mark.req("TIM-008")
def test_lag_correlation_finds_a_known_delay():
    t = np.linspace(0.0, 1000.0, 1001)
    driver = np.sin(2.0 * math.pi * t / 200.0)
    response = np.sin(2.0 * math.pi * (t - 25.0) / 200.0)
    lag, correlation = lag_correlation(t, driver, response, max_lag_samples=100)
    assert lag == pytest.approx(25.0, abs=1.0)
    assert correlation > 0.999


# --- the prohibition -------------------------------------------------------


@pytest.mark.req("TIM-008")
def test_measurements_carry_all_four_labels():
    """TIM-008 is about attribution as much as about numbers."""
    m = _measurement(1.0e7)
    assert m.observable and m.subsystem and m.mechanism and m.forcing_regime
    assert m.timescale_days == pytest.approx(1.0e7 / 86400.0)
    assert m.mechanism in m.describe()
    assert m.forcing_regime in m.describe()


@pytest.mark.req("TIM-008", "INF-006")
def test_a_single_climate_index_is_refused_when_the_timescales_disagree():
    """The one number the charter's timescale philosophy exists to prevent.

    Refusing rather than averaging is the point: an average of a two-week
    adjustment and a two-century one describes no mechanism in the system, and a
    reader given it has no way to know that.
    """
    report = TimescaleReport()
    report.add(_measurement(4.0e6, "radiative_relaxation"))
    report.add(_measurement(3.0e9, "deep_ocean_uptake"))
    assert report.spread == pytest.approx(750.0)
    with pytest.raises(InvalidForScientificInferenceError, match="single climate index"):
        report.single_index()
    # And the refusal names what to report instead.
    with pytest.raises(InvalidForScientificInferenceError, match="deep_ocean_uptake"):
        report.single_index()


@pytest.mark.req("TIM-008")
def test_one_index_is_allowed_when_the_measurements_actually_agree():
    """The rule is about honesty, not about never quoting a number."""
    report = TimescaleReport()
    report.add(_measurement(4.0e6, "a"))
    report.add(_measurement(4.4e6, "b"))
    assert report.single_index() == pytest.approx(4.2e6)


@pytest.mark.req("TIM-008")
def test_measurements_are_selectable_by_every_label():
    report = TimescaleReport()
    report.add(_measurement(1.0e6, "fast"))
    report.add(_measurement(1.0e9, "slow"))
    assert len(report.select(mechanism="fast")) == 1
    assert len(report.select(forcing_regime="step")) == 2
    assert report.select(forcing_regime="periodic") == []
