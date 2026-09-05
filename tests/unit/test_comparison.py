"""Putting a model series and an observed one on comparable footing.

Requirements DAT-004, DAT-008, INF-007, TIM-004, VAL-003, CLM-002, VAL-007.

Two operations stand between a model's output and an observational record.
Getting the first wrong produces an error of order 250 K, which is obvious.
Getting it subtly wrong -- differencing against the wrong window -- produces a
constant offset that looks like a physical bias, which is not.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError, ProvenanceError
from mfespf.data.records import (
    ProvenanceClass,
    SamplingConvention,
    TimeSeries,
    UncertaintyKind,
)
from mfespf.evaluation.comparison import (
    align_to_record,
    anomaly_against_baseline,
    model_series,
)
from mfespf.registries.claims import ClaimStatus
from mfespf.run.modes import OperatingMode

pytestmark = pytest.mark.stage(6)

YEAR = 365.25 * 86400.0


def _record(n: int = 200, offset: float = 0.0, units: str = "K") -> TimeSeries:
    times = np.arange(n, dtype=float) * (YEAR / 12.0)
    values = offset + 0.3 * np.sin(2.0 * np.pi * times / (11.0 * YEAR))
    return TimeSeries(
        variable="surface_temperature_anomaly",
        times=times,
        values=values,
        units=units,
        sampling=SamplingConvention.INTERVAL_MEAN,
        provenance_class=ProvenanceClass.OBSERVATIONAL,
        uncertainty=np.full(n, 0.1),
        uncertainty_kind=UncertaintyKind.HALF_WIDTH_95,
    )


# --- differencing to a baseline --------------------------------------------


@pytest.mark.req("INF-007", "DAT-008")
def test_the_baseline_subtracted_is_the_series_own_mean_and_is_recorded():
    absolute = _record(offset=288.0)
    start, end = 5.0 * YEAR, 10.0 * YEAR
    anomaly, step = anomaly_against_baseline(absolute, start=start, end=end, label="test")

    inside = np.logical_and(absolute.times >= start, absolute.times <= end)
    expected = float(np.mean(absolute.values[inside]))
    assert float(np.mean(anomaly.values[inside])) == pytest.approx(0.0, abs=1e-12)
    assert anomaly.values[0] == pytest.approx(absolute.values[0] - expected)

    assert step.operation == "anomaly:surface_temperature_anomaly"
    assert f"subtracted={expected:.10g}" in step.parameters
    assert "baseline_window=" in step.parameters
    assert "not propagated" in step.uncertainty_effect


@pytest.mark.req("VAL-003")
def test_differencing_against_the_wrong_window_shows_up_as_a_bias():
    """The failure mode worth having a test for.

    Using a convenient window instead of the record's own shifts the model by a
    constant. It is not a large error and it looks exactly like a physical
    offset, which is why the window is a required argument rather than a
    default.
    """
    absolute = _record(offset=288.0)
    right, _ = anomaly_against_baseline(absolute, start=0.0, end=5.0 * YEAR)
    wrong, _ = anomaly_against_baseline(absolute, start=8.0 * YEAR, end=13.0 * YEAR)
    difference = wrong.values - right.values
    assert np.allclose(difference, difference[0])  # a pure constant offset
    assert abs(difference[0]) > 0.05


@pytest.mark.req("SYS-006")
def test_an_unusable_baseline_window_is_refused():
    absolute = _record(offset=288.0)
    with pytest.raises(ConfigurationError, match="is empty"):
        anomaly_against_baseline(absolute, start=YEAR, end=YEAR)
    with pytest.raises(ProvenanceError, match="beyond its coverage"):
        anomaly_against_baseline(absolute, start=0.0, end=500.0 * YEAR)
    with pytest.raises(ConfigurationError, match="not a baseline"):
        anomaly_against_baseline(absolute, start=0.0, end=YEAR / 24.0)


# --- alignment --------------------------------------------------------------


@pytest.mark.req("TIM-004", "VAL-003")
def test_the_model_is_moved_onto_the_record_and_never_the_reverse():
    """Resampling the record onto the model's clock would let the model choose
    how much of the observations survive."""
    observed = _record(n=240)
    modelled = TimeSeries(
        variable="surface_temperature_anomaly",
        times=observed.times[::2] + 1.0e5,
        values=np.zeros(observed.times[::2].size),
        units="K",
        sampling=SamplingConvention.INSTANTANEOUS,
        provenance_class=ProvenanceClass.MODEL_OUTPUT,
    )
    comparison = align_to_record(modelled, observed)

    # Every paired time is one of the record's own, not one of the model's.
    assert set(comparison.times.tolist()) <= set(observed.times.tolist())
    assert comparison.n > 100
    assert any(t.operation.startswith("align:") for t in comparison.transformations)
    assert any("direction=model_onto_record" in t.parameters for t in comparison.transformations)


@pytest.mark.req("DAT-004", "SYS-006")
def test_series_in_different_units_or_with_no_overlap_are_refused():
    observed = _record()
    wrong_units = _record(units="W m^-2")
    with pytest.raises(ProvenanceError, match="convert at a declared boundary"):
        align_to_record(wrong_units, observed)

    elsewhere = observed.with_values(observed.times + 1000.0 * YEAR, observed.values)
    with pytest.raises(ProvenanceError, match="do not overlap"):
        align_to_record(elsewhere, observed)


@pytest.mark.req("VAL-003")
def test_a_comparison_reports_what_the_model_did_and_did_not_explain():
    observed = _record()
    perfect = observed.with_values(observed.times, observed.values)
    flat = observed.with_values(observed.times, np.zeros_like(observed.values))

    good = align_to_record(perfect, observed)
    assert good.rms_error == pytest.approx(0.0, abs=1e-9)
    assert good.explained_fraction() == pytest.approx(1.0, abs=1e-9)
    assert good.within_observational_uncertainty() == pytest.approx(1.0)

    # A model that predicts nothing explains nothing. Not a negative number and
    # not a correlation: zero.
    none = align_to_record(flat, observed)
    assert none.explained_fraction() == pytest.approx(0.0, abs=1e-9)
    assert none.modelled_range == pytest.approx(0.0)
    assert "variance explained" in none.describe()


@pytest.mark.req("VAL-003")
def test_a_model_worse_than_a_flat_line_reports_a_negative_explained_fraction():
    """Worth being able to state plainly rather than reporting a correlation."""
    observed = _record()
    inverted = observed.with_values(observed.times, -observed.values)
    comparison = align_to_record(inverted, observed)
    assert comparison.explained_fraction() < -1.0


# --- model output as a series ------------------------------------------------


@pytest.mark.req("TIM-001")
def test_a_model_series_is_timestamped_at_the_end_of_its_step():
    """A committed state in a step's outputs is the value at the step's end.

    Using the start time would put the whole series one step early, which at a
    monthly step against a two-year response is a twentieth of a phase.
    """

    class _Step:
        def __init__(self, time, dt, value):
            self.time, self.dt = time, dt
            self.outputs = {"surface_temperature": type("F", (), {"values": np.array(value)})()}

    steps = [_Step(i * 100.0, 100.0, 280.0 + i) for i in range(5)]
    series = model_series(steps, "surface_temperature", units="K")
    assert series.times[0] == pytest.approx(100.0)
    assert series.times[-1] == pytest.approx(500.0)
    assert series.provenance_class is ProvenanceClass.MODEL_OUTPUT
    assert not series.provenance_class.constrains_the_world

    with pytest.raises(ConfigurationError, match="no steps"):
        model_series([], "surface_temperature", units="K")


# --- what a result is worth --------------------------------------------------


@pytest.mark.req("CLM-002", "VAL-007", "TIM-006")
def test_a_failure_is_decisive_and_a_reconstruction_pass_is_only_qualified():
    """The asymmetry this whole comparison path exists to encode.

    A mechanism handed the boundary conditions of the period it is judged on and
    still unable to account for it has been refuted. The same mechanism
    succeeding under those conditions has not been supported -- it was told part
    of the answer.
    """
    for mode in (OperatingMode.PRESCRIBED_HISTORY, OperatingMode.HISTORICAL_RECONSTRUCTION):
        assert mode.status_for_result(False) is ClaimStatus.REFUTED
        assert mode.status_for_result(True) is ClaimStatus.QUALIFIED

    assert (
        OperatingMode.OUT_OF_SAMPLE_PREDICTION.status_for_result(True)
        is ClaimStatus.SUPPORTED_WITHIN_SCOPE
    )
    assert OperatingMode.HINDCAST.status_for_result(False) is ClaimStatus.REFUTED
    # A mode that compared against nothing settles nothing, either way.
    for mode in (OperatingMode.EQUILIBRIUM, OperatingMode.CONTROLLED_EXPERIMENT):
        assert mode.status_for_result(True) is ClaimStatus.EMPIRICALLY_UNRESOLVED
        assert mode.status_for_result(False) is ClaimStatus.EMPIRICALLY_UNRESOLVED
