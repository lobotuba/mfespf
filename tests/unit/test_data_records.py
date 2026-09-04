"""Stage 6: declared records and what resampling does to them.

Charter 11.4 and requirements DAT-007, DAT-008, INF-002, INF-007, TIM-004.

Charter 11.4 asks a resampling to declare aggregation, interpolation,
uncertainty and lost frequency content. The first three are bookkeeping and are
checked here as such. The fourth gets most of the file, because it is the one
that decides whether a later comparison means anything: a record averaged onto a
coarse clock has had variability removed, and a model that reproduces the
smoothed version has reproduced less than it appears to.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError, ProvenanceError
from mfespf.data.records import (
    ProvenanceClass,
    SamplingConvention,
    TimeSeries,
    resample,
)

pytestmark = pytest.mark.stage(6)

DAY = 86400.0


def _series(
    n: int = 400,
    interval: float = DAY,
    sampling: SamplingConvention = SamplingConvention.INTERVAL_MEAN,
    fast_amplitude: float = 0.3,
) -> TimeSeries:
    times = np.arange(n, dtype=float) * interval
    values = (
        1361.0
        + 0.8 * np.sin(2.0 * np.pi * times / (200.0 * DAY))
        + fast_amplitude * np.sin(2.0 * np.pi * times / (4.0 * DAY))
    )
    return TimeSeries(
        variable="total_solar_irradiance",
        times=times,
        values=values,
        units="W m^-2",
        sampling=sampling,
        dataset_id="TEST-TSI",
        dataset_version="1.0.0",
        provenance_class=ProvenanceClass.SYNTHETIC,
    )


# --- what a record is ------------------------------------------------------


@pytest.mark.req("DAT-003", "TIM-001")
def test_a_record_reports_what_it_can_resolve():
    """The shortest period a sampling can represent is twice its interval.

    Reported because it is the number that says which of a record's features a
    model could possibly be asked to reproduce; a claim about anything faster is
    about the interpolation, not the data.
    """
    series = _series(n=100, interval=DAY)
    assert series.median_interval == pytest.approx(DAY)
    assert series.resolved_period == pytest.approx(2.0 * DAY)
    assert series.is_uniform
    assert series.duration == pytest.approx(99.0 * DAY)
    assert series.coverage == (0.0, 99.0 * DAY)


@pytest.mark.req("DAT-007", "SYS-006")
def test_a_record_refuses_to_be_read_outside_its_coverage():
    """Extrapolating a record would invent a boundary condition."""
    series = _series(n=10)
    series.require_coverage(0.0, 9.0 * DAY)
    with pytest.raises(ProvenanceError, match="beyond its coverage"):
        series.require_coverage(0.0, 20.0 * DAY)
    with pytest.raises(ProvenanceError, match="beyond its coverage"):
        series.value_at(20.0 * DAY)


@pytest.mark.req("DAT-007", "SYS-006")
def test_a_gap_is_a_shorter_series_not_a_sentinel():
    """A missing value encoded as a number is a number something will multiply."""
    times = np.arange(5.0) * DAY
    with pytest.raises(ConfigurationError, match="non-finite"):
        TimeSeries(
            variable="x",
            times=times,
            values=np.array([1.0, 2.0, np.nan, 4.0, 5.0]),
            units="W m^-2",
            sampling=SamplingConvention.INTERVAL_MEAN,
        )


@pytest.mark.req("SYS-006")
def test_malformed_records_are_refused():
    times = np.arange(5.0) * DAY
    with pytest.raises(ConfigurationError, match="same length"):
        TimeSeries(
            variable="x",
            times=times,
            values=np.zeros(4),
            units="W m^-2",
            sampling=SamplingConvention.INTERVAL_MEAN,
        )
    with pytest.raises(ConfigurationError, match="non-increasing"):
        TimeSeries(
            variable="x",
            times=np.array([0.0, 2.0, 1.0]),
            values=np.zeros(3),
            units="W m^-2",
            sampling=SamplingConvention.INTERVAL_MEAN,
        )
    with pytest.raises(Exception, match="not registered"):
        TimeSeries(
            variable="x",
            times=times,
            values=np.zeros(5),
            units="furlongs per fortnight",
            sampling=SamplingConvention.INTERVAL_MEAN,
        )


# --- resampling ------------------------------------------------------------


@pytest.mark.req("TIM-004", "INF-007")
def test_coarsening_reports_the_variability_it_removed():
    """The charter's fourth demand, and the reason the other three matter.

    A four-day component cannot survive a fifty-day clock. The report says so in
    two ways a reader can act on: which band of periods is gone, and how much of
    the record's variance went with it.
    """
    series = _series(n=400, interval=DAY, fast_amplitude=0.6)
    targets = np.arange(25.0, 375.0, 50.0) * DAY
    _, report = resample(series, targets)

    assert report.is_coarsening
    assert report.target_resolved_period == pytest.approx(100.0 * DAY)
    assert report.variance_retained < 0.95
    assert "not represented after resampling" in report.lost_frequency_note
    assert f"{report.source_resolved_period:g}" in report.lost_frequency_note


@pytest.mark.req("TIM-004")
def test_refining_reports_that_nothing_was_lost():
    """The null case: a finer clock removes no frequency content."""
    series = _series(n=60, interval=10.0 * DAY, fast_amplitude=0.0)
    targets = np.arange(5.0, 585.0, 5.0) * DAY
    _, report = resample(series, targets)
    assert not report.is_coarsening
    assert "no frequency content lost" in report.lost_frequency_note
    assert report.variance_retained > 0.98


@pytest.mark.req("DAT-008", "INF-007")
def test_the_transformation_record_replays_to_the_same_numbers():
    """DAT-008: preprocessing reproducible from a versioned transformation record.

    The record names the operation, its parameters and the software; running it
    again from those must give the identical array, or the record documents
    something other than what happened.
    """
    series = _series(n=200)
    targets = np.arange(10.0, 190.0, 20.0) * DAY

    first, report = resample(series, targets)
    step = report.as_transformation_step()
    assert step.operation == "resample:total_solar_irradiance"
    assert "aggregation=interval_mean" in step.parameters
    assert "interpolation=piecewise_constant" in step.parameters
    assert step.software.endswith("resample")
    assert "UNRESOLVED" in step.uncertainty_effect

    second, _ = resample(series, targets, method=report.aggregation)
    assert np.array_equal(first.values, second.values)


@pytest.mark.req("DAT-003", "SYS-006")
def test_a_method_the_sampling_convention_forbids_is_refused():
    """Averaging instantaneous samples and interpolating means are both wrong.

    Both produce a plausible array and leave no trace of the mistake in it,
    which is exactly why the convention has to be carried and checked rather
    than remembered.
    """
    means = _series(sampling=SamplingConvention.INTERVAL_MEAN)
    instants = _series(sampling=SamplingConvention.INSTANTANEOUS)
    targets = np.arange(10.0, 390.0, 20.0) * DAY

    with pytest.raises(ConfigurationError, match="not valid for a interval_mean"):
        resample(means, targets, method="linear")
    with pytest.raises(ConfigurationError, match="not valid for a instantaneous"):
        resample(instants, targets, method="interval_mean")
    # And the defaults are the ones each convention permits.
    assert resample(means, targets)[1].aggregation == "interval_mean"
    assert resample(instants, targets)[1].interpolation == "linear"


@pytest.mark.req("TIM-004", "CON-009")
def test_an_interval_mean_resampling_preserves_the_record_mean():
    """Averaging is conservative in the mean, which is what makes it an average.

    Not exact at the edges, where the outermost target cells extend half a step
    past the record and hold the end values; the interior is exact and that is
    what the tolerance below reflects.
    """
    series = _series(n=400, interval=DAY, fast_amplitude=0.0)
    targets = np.arange(20.0, 380.0, 40.0) * DAY
    resampled, _ = resample(series, targets)
    interior = np.logical_and(series.times >= targets[0], series.times <= targets[-1])
    assert float(np.mean(resampled.values)) == pytest.approx(
        float(np.mean(series.values[interior])), rel=2e-4
    )


@pytest.mark.req("DAT-007")
def test_resampling_past_the_record_coverage_is_refused():
    series = _series(n=50)
    with pytest.raises(ProvenanceError, match="beyond its coverage"):
        resample(series, np.arange(0.0, 100.0, 5.0) * DAY)


# --- provenance class ------------------------------------------------------


@pytest.mark.req("INF-002")
def test_provenance_class_says_what_a_record_is_not_what_it_was_used_for():
    """A role and a provenance class answer different questions.

    A synthetic record can be correctly declared as withheld evaluation data.
    The role is right; what it cannot tell a reader is that agreement with the
    record says nothing about the world.
    """
    assert ProvenanceClass.OBSERVATIONAL.constrains_the_world
    assert ProvenanceClass.PROXY.constrains_the_world
    assert ProvenanceClass.REANALYSIS.constrains_the_world
    for generated in (
        ProvenanceClass.MODEL_OUTPUT,
        ProvenanceClass.SYNTHETIC,
        ProvenanceClass.IDEALISED,
    ):
        assert generated.is_generated
        assert not generated.constrains_the_world
    # Undeclared is neither, which is why it is refused in an evidence role.
    assert not ProvenanceClass.UNDECLARED.constrains_the_world
    assert not ProvenanceClass.UNDECLARED.is_generated
