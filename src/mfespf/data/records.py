"""Historical input records, and what happens to them on the way in.

Charter 11.4 and 10.4; requirements DAT-008, INF-002, INF-007, TIM-004.

Charter 11.4 is one sentence with four demands in it: *any resampling shall
declare aggregation, interpolation, uncertainty, and lost frequency content*.
The first three are bookkeeping. The fourth is the one that decides whether a
reconstruction means anything, because a record averaged onto a coarse model
clock has had variance removed, and a model that then reproduces the smoothed
record has reproduced less than it appears to.

So ``resample`` returns the resampled series *and* a report, and the report says
how much variance survived and what the shortest resolvable period became. The
report converts to a ``TransformationStep`` for the dataset registry, which is
what DAT-008 means by a versioned transformation record: the operation, its
parameters, the software that ran it, and its effect on uncertainty, in a form
that can be replayed to reproduce the same numbers.

**Provenance class is separate from role.** A dataset's *role* (INF-002) says
what it was used for; its *provenance class* says what kind of thing it is. A
synthetic record used as withheld evidence is a legitimate test of the
evaluation pathway and is not evidence about the world, and only the second of
those facts is invisible from the role alone.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from enum import StrEnum
from typing import Any

import numpy as np

from mfespf.core.errors import ConfigurationError, ProvenanceError
from mfespf.core.units import parse_unit
from mfespf.registries.datasets import ProvenanceClass, TransformationStep

__all__ = [
    "ProvenanceClass",
    "ResamplingReport",
    "SamplingConvention",
    "TimeSeries",
    "UncertaintyKind",
    "cell_edges",
    "resample",
]

SOFTWARE = "mfespf.data.records.resample"


class UncertaintyKind(StrEnum):
    """What a record's uncertainty column means (UNC-002).

    Added because real records disagree. HadCRUT5 publishes a 95 percent
    interval, Berkeley Earth a 95 percent interval it calls an uncertainty,
    NOAA's heat content a standard error, and NRLTSI2 a standard uncertainty.
    Combining them as though they were the same quantity would understate some
    by a factor of two, and nothing in the numbers says which is which -- so the
    meaning travels with the array rather than living in whoever loaded it.
    """

    UNDECLARED = "undeclared"
    STANDARD_UNCERTAINTY = "standard_uncertainty"
    STANDARD_ERROR = "standard_error"
    HALF_WIDTH_95 = "half_width_95"

    @property
    def is_one_sigma(self) -> bool:
        """Whether the value may be used directly in a quadrature combination.

        A half-width of a 95 percent interval may not: it is roughly two sigma
        for a normal distribution, and treating it as one would halve the
        combined uncertainty of anything it entered (UNC-006).
        """
        return self in (
            UncertaintyKind.STANDARD_UNCERTAINTY,
            UncertaintyKind.STANDARD_ERROR,
        )


class SamplingConvention(StrEnum):
    """What one sample of a series means (DAT-003, charter 11.4).

    The difference decides how the series may be resampled: an interval mean can
    be averaged again, an instantaneous sample cannot without an interpolation
    assumption, and an interval total must be summed rather than averaged. A
    record that does not say which it is cannot be resampled correctly, which is
    why there is no default.
    """

    INSTANTANEOUS = "instantaneous"
    INTERVAL_MEAN = "interval_mean"
    INTERVAL_TOTAL = "interval_total"


@dataclass(frozen=True)
class TimeSeries:
    """One variable of a declared record, on its own clock.

    Times are seconds on the master coordinate (TIM-001), strictly increasing.
    Missing values are absent, not encoded: a gap is a shorter series, because
    DAT-007 forbids a sentinel that a consumer might arithmetic on by mistake.
    """

    variable: str
    times: np.ndarray
    values: np.ndarray
    units: str
    sampling: SamplingConvention
    dataset_id: str = ""
    dataset_version: str = ""
    provenance_class: ProvenanceClass = ProvenanceClass.UNDECLARED
    uncertainty: np.ndarray | None = None
    uncertainty_kind: UncertaintyKind = UncertaintyKind.UNDECLARED
    note: str = ""

    def __post_init__(self) -> None:
        times = np.asarray(self.times, dtype=float)
        values = np.asarray(self.values, dtype=float)
        object.__setattr__(self, "times", times)
        object.__setattr__(self, "values", values)
        parse_unit(self.units)
        if times.ndim != 1 or values.shape != times.shape:
            raise ConfigurationError(
                f"series {self.variable!r} has {times.shape} times and {values.shape} values; "
                f"they must be one-dimensional and the same length"
            )
        if times.size < 2:
            raise ConfigurationError(f"series {self.variable!r} needs at least two samples")
        if np.any(np.diff(times) <= 0.0):
            raise ConfigurationError(
                f"series {self.variable!r} has non-increasing times; a record whose order is "
                f"not monotone cannot be resampled without an undeclared sort"
            )
        if not np.all(np.isfinite(values)):
            raise ConfigurationError(
                f"series {self.variable!r} contains non-finite values; a gap is represented by "
                f"a shorter series, never by a sentinel (DAT-007)"
            )
        if self.uncertainty is not None:
            uncertainty = np.asarray(self.uncertainty, dtype=float)
            object.__setattr__(self, "uncertainty", uncertainty)
            if uncertainty.shape != values.shape:
                raise ConfigurationError(
                    f"series {self.variable!r} has {uncertainty.shape} uncertainties for "
                    f"{values.shape} values"
                )
            if self.uncertainty_kind is UncertaintyKind.UNDECLARED:
                raise ConfigurationError(
                    f"series {self.variable!r} carries uncertainties but does not say what "
                    f"they are; a standard error and a 95 percent half-width differ by a "
                    f"factor of two and nothing in the numbers distinguishes them (UNC-002)"
                )

    # --- coverage ---------------------------------------------------------
    @property
    def coverage(self) -> tuple[float, float]:
        return (float(self.times[0]), float(self.times[-1]))

    @property
    def duration(self) -> float:
        return float(self.times[-1] - self.times[0])

    @property
    def median_interval(self) -> float:
        return float(np.median(np.diff(self.times)))

    @property
    def is_uniform(self) -> bool:
        intervals = np.diff(self.times)
        return bool(np.allclose(intervals, intervals[0], rtol=1e-9, atol=0.0))

    @property
    def resolved_period(self) -> float:
        """The shortest period the sampling can represent: twice the interval.

        Reported rather than assumed, because it is the number that says which
        of a record's features a model could possibly be asked to reproduce.
        """
        return 2.0 * self.median_interval

    def covers(self, start: float, end: float) -> bool:
        low, high = self.coverage
        return low <= start and end <= high

    def require_coverage(self, start: float, end: float) -> None:
        if not self.covers(start, end):
            low, high = self.coverage
            raise ProvenanceError(
                f"series {self.variable!r} covers [{low:g}, {high:g}] s but was asked for "
                f"[{start:g}, {end:g}] s; extrapolating a record beyond its coverage would "
                f"invent a boundary condition (DAT-007)"
            )

    # --- reading ----------------------------------------------------------
    def value_at(self, time: float) -> float:
        """The sample in force at ``time``, held from the last one at or before it.

        A zero-order hold, not an interpolation: interpolating would smooth a
        step the record may intend. It is deliberately *not* the cell
        reconstruction ``resample`` uses, which puts the boundaries halfway
        between samples because that is what an interval mean asserts. The two
        agree exactly at the sample times, which is where a replay reads them,
        and differ by at most half an interval elsewhere.
        """
        self.require_coverage(time, time)
        index = int(np.searchsorted(self.times, time, side="right") - 1)
        return float(self.values[max(index, 0)])

    def with_values(self, times: np.ndarray, values: np.ndarray, **overrides: Any) -> TimeSeries:
        fields: dict[str, Any] = {
            "variable": self.variable,
            "times": times,
            "values": values,
            "units": self.units,
            "sampling": self.sampling,
            "dataset_id": self.dataset_id,
            "dataset_version": self.dataset_version,
            "provenance_class": self.provenance_class,
            "uncertainty_kind": self.uncertainty_kind,
            "note": self.note,
        }
        fields.update(overrides)
        return TimeSeries(**fields)

    def describe(self) -> str:
        low, high = self.coverage
        return (
            f"{self.variable} [{self.units}] {self.sampling.value}, {self.times.size} samples "
            f"over [{low:g}, {high:g}] s, resolving periods down to {self.resolved_period:g} s "
            f"({self.provenance_class.value})"
        )


@dataclass(frozen=True)
class ResamplingReport:
    """What a resampling did, and what it cost (charter 11.4, INF-007).

    ``variance_retained`` is the honest headline. A record resampled onto a
    coarser clock has had variance removed, and a model compared against the
    result is being asked to reproduce a smoother signal than the instrument
    recorded. TIM-004 requires that be stated rather than assumed harmless.
    """

    variable: str
    dataset_id: str
    dataset_version: str
    source_coverage: tuple[float, float]
    source_samples: int
    target_samples: int
    aggregation: str
    interpolation: str
    source_resolved_period: float
    target_resolved_period: float
    input_variance: float
    output_variance: float
    uncertainty_effect: str
    software: str = SOFTWARE

    @property
    def variance_retained(self) -> float:
        if self.input_variance == 0.0:
            return 1.0
        return self.output_variance / self.input_variance

    @property
    def is_coarsening(self) -> bool:
        return self.target_resolved_period > self.source_resolved_period

    @property
    def lost_frequency_note(self) -> str:
        if not self.is_coarsening:
            return (
                "no frequency content lost: the target clock resolves at least as short a "
                "period as the record does"
            )
        return (
            f"variability at periods between {self.source_resolved_period:g} s and "
            f"{self.target_resolved_period:g} s is not represented after resampling; "
            f"{100.0 * (1.0 - self.variance_retained):.2f} percent of the record's variance "
            f"was removed, and a model reproducing the result has reproduced a smoother "
            f"signal than the record holds"
        )

    def as_transformation_step(self) -> TransformationStep:
        """The DAT-008 record: enough to replay this exact operation."""
        return TransformationStep(
            operation=f"resample:{self.variable}",
            parameters=(
                f"aggregation={self.aggregation}; interpolation={self.interpolation}; "
                f"source_coverage=[{self.source_coverage[0]:g},{self.source_coverage[1]:g}]; "
                f"source_samples={self.source_samples}; target_samples={self.target_samples}; "
                f"target_resolved_period={self.target_resolved_period:g}"
            ),
            software=self.software,
            uncertainty_effect=f"{self.uncertainty_effect}; {self.lost_frequency_note}",
        )

    def describe(self) -> str:
        return (
            f"{self.variable}: {self.source_samples} -> {self.target_samples} samples by "
            f"{self.aggregation}/{self.interpolation}; variance retained "
            f"{self.variance_retained:.4f}; {self.lost_frequency_note}"
        )


def resample(
    series: TimeSeries,
    target_times: np.ndarray,
    *,
    method: str | None = None,
    target_interval: float | None = None,
) -> tuple[TimeSeries, ResamplingReport]:
    """Put a record on a new clock, and report what that cost.

    ``method`` defaults to the one the sampling convention permits: an interval
    mean is averaged, an interval total is summed, and an instantaneous sample is
    linearly interpolated. Passing a method the convention does not support is
    refused rather than silently accepted, because averaging instantaneous
    samples and interpolating interval means are both wrong in ways that leave no
    trace in the output.
    """
    targets = np.asarray(target_times, dtype=float)
    if targets.ndim != 1 or targets.size < 2:
        raise ConfigurationError("target_times must be a one-dimensional array of two or more")
    if np.any(np.diff(targets) <= 0.0):
        raise ConfigurationError("target_times must be strictly increasing")

    default = {
        SamplingConvention.INTERVAL_MEAN: "interval_mean",
        SamplingConvention.INTERVAL_TOTAL: "interval_total",
        SamplingConvention.INSTANTANEOUS: "linear",
    }[series.sampling]
    chosen = method or default
    permitted = {
        SamplingConvention.INTERVAL_MEAN: {"interval_mean", "zero_order_hold"},
        SamplingConvention.INTERVAL_TOTAL: {"interval_total", "zero_order_hold"},
        SamplingConvention.INSTANTANEOUS: {"linear", "zero_order_hold"},
    }[series.sampling]
    if chosen not in permitted:
        raise ConfigurationError(
            f"method {chosen!r} is not valid for a {series.sampling.value} series; "
            f"permitted methods are {sorted(permitted)}. Averaging instantaneous samples and "
            f"interpolating interval means are both wrong in ways the output does not show"
        )

    interval = (
        target_interval if target_interval is not None else float(np.median(np.diff(targets)))
    )
    series.require_coverage(float(targets[0]), float(targets[-1]))

    if chosen == "linear":
        values = np.interp(targets, series.times, series.values)
        aggregation, interpolation = "none", "linear"
    elif chosen == "zero_order_hold":
        indices = np.clip(np.searchsorted(series.times, targets, side="right") - 1, 0, None)
        values = series.values[indices]
        aggregation, interpolation = "none", "zero_order_hold"
    else:
        # Overlap-weighted integration of the source cells over each target cell.
        # Both sets of cells are centred on their sample times, which is what an
        # interval mean asserts, so resampling a series onto its own clock
        # returns the series exactly.
        edges = np.empty(targets.size + 1)
        edges[1:-1] = 0.5 * (targets[:-1] + targets[1:])
        edges[0] = targets[0] - 0.5 * interval
        edges[-1] = targets[-1] + 0.5 * interval
        source = cell_edges(series.times)
        values = np.array(
            [_piecewise_mean(series.values, source, lo, hi) for lo, hi in itertools.pairwise(edges)]
        )
        if chosen == "interval_total":
            values = values * interval
        aggregation = chosen
        interpolation = "piecewise_constant"

    resampled = series.with_values(
        targets,
        values,
        note=(f"{series.note}; resampled onto a {interval:g} s clock" if series.note else ""),
    )
    report = ResamplingReport(
        variable=series.variable,
        dataset_id=series.dataset_id,
        dataset_version=series.dataset_version,
        source_coverage=series.coverage,
        source_samples=int(series.times.size),
        target_samples=int(targets.size),
        aggregation=aggregation,
        interpolation=interpolation,
        source_resolved_period=series.resolved_period,
        target_resolved_period=2.0 * interval,
        input_variance=float(np.var(series.values)),
        output_variance=float(np.var(values)),
        uncertainty_effect=(
            "sample uncertainties are not propagated through the resampling and remain "
            "UNRESOLVED; averaging correlated samples would not reduce them the way "
            "independent ones would (UNC-006)"
        ),
    )
    return resampled, report


def cell_edges(times: np.ndarray) -> np.ndarray:
    """The interval boundaries an interval-mean series' samples assert.

    A sample of an interval mean sits at the *middle* of the interval it
    describes, so the boundaries are halfway between consecutive samples, with
    half-cells at each end.

    Reconstructing the series with breakpoints at the sample times instead --
    which is what a zero-order hold does -- shifts every cell by half an
    interval. Resampling a series onto its own clock then fails to return it,
    and the discrepancy comes out as a variance loss that looks like the
    resampling removed real variability. It had been doing exactly that, and
    reporting two percent of it as though it were physical.
    """
    if times.size < 2:
        raise ConfigurationError("cell edges need at least two samples")
    edges = np.empty(times.size + 1)
    edges[1:-1] = 0.5 * (times[:-1] + times[1:])
    edges[0] = times[0] - 0.5 * (times[1] - times[0])
    edges[-1] = times[-1] + 0.5 * (times[-1] - times[-2])
    return edges


def _piecewise_mean(values: np.ndarray, source_edges: np.ndarray, low: float, high: float) -> float:
    """Mean over ``[low, high]`` of the step function those cells describe.

    The overlap of the target cell with each source cell, weighted by that
    source cell's value. Where the target cells are the source cells this
    returns the source values unchanged, which is what makes resampling onto
    the same clock the identity.
    """
    if high <= low:
        index = np.clip(np.searchsorted(source_edges, low, side="right") - 1, 0, values.size - 1)
        return float(values[int(index)])
    left = np.maximum(source_edges[:-1], low)
    right = np.minimum(source_edges[1:], high)
    overlap = np.clip(right - left, 0.0, None)
    covered = float(np.sum(overlap))
    if covered <= 0.0:
        index = np.clip(np.searchsorted(source_edges, low, side="right") - 1, 0, values.size - 1)
        return float(values[int(index)])
    # Divided by the width actually covered rather than by the target width: a
    # target cell reaching past the outermost source edge is averaged over the
    # part the record covers instead of being diluted towards zero.
    return float(np.sum(values * overlap) / covered)
