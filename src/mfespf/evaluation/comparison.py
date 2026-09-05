"""Putting a model series and an observed series on comparable footing.

Requirements DAT-008, INF-007, TIM-004, UNC-002, VAL-003.

Two operations stand between a model's output and an observational record, and
both are transformations that change what a comparison means. Neither is
performed silently here.

**Differencing to a baseline.** An observational anomaly is a temperature minus
the mean of that record over a stated window. A model producing absolute
temperature must be differenced against *its own* mean over *the same* window
before the two are comparable. Skipping it compares an absolute temperature with
a difference and produces an error of order 250 K; doing it against a different
window produces a smaller error that looks like a physical result.

**Alignment.** The two series are on different clocks and cover different
periods. The overlap is computed, the model is resampled onto the record's times
through the same reporting path as any other resampling, and the record is never
resampled onto the model -- a comparison should move the model to the data, not
the data to the model.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from mfespf.core.errors import ConfigurationError, ProvenanceError
from mfespf.data.records import (
    ProvenanceClass,
    SamplingConvention,
    TimeSeries,
    UncertaintyKind,
    resample,
)
from mfespf.registries.datasets import TransformationStep

__all__ = ["Comparison", "align_to_record", "anomaly_against_baseline", "model_series"]


def model_series(
    steps: list,
    field_id: str,
    *,
    units: str,
    variable: str | None = None,
    sampling: SamplingConvention = SamplingConvention.INSTANTANEOUS,
) -> TimeSeries:
    """A model output as a series, timestamped at the end of the step it belongs to.

    ``StepReport.time`` is when the step *started*; a committed state in its
    outputs is the value at the end. Using the start time would put the whole
    series one step early -- negligible when the step is a hundredth of the
    response time, and a twentieth of a phase when it is a month against a
    two-year timescale, which is where this is used.
    """
    if not steps:
        raise ConfigurationError("the run produced no steps to build a series from")
    times = np.array([s.time + s.dt for s in steps])
    values = np.array([float(s.outputs[field_id].values) for s in steps])
    return TimeSeries(
        variable=variable or field_id,
        times=times,
        values=values,
        units=units,
        sampling=sampling,
        provenance_class=ProvenanceClass.MODEL_OUTPUT,
        note=f"{field_id} from a model run, at end-of-step instants",
    )


def anomaly_against_baseline(
    series: TimeSeries, *, start: float, end: float, label: str = ""
) -> tuple[TimeSeries, TransformationStep]:
    """Subtract the series' own mean over ``[start, end]``, and say so.

    The window must be the observational record's baseline window, not a
    convenient one. Differencing against a different window shifts the whole
    model series by a constant, which shows up as a bias a reader may mistake
    for a physical offset.
    """
    if end <= start:
        raise ConfigurationError(f"baseline window [{start}, {end}] is empty")
    series.require_coverage(start, end)
    inside = np.logical_and(series.times >= start, series.times <= end)
    if int(np.count_nonzero(inside)) < 2:
        raise ConfigurationError(
            f"the baseline window [{start:g}, {end:g}] contains "
            f"{int(np.count_nonzero(inside))} samples of {series.variable!r}; a mean taken "
            f"from one sample is that sample, not a baseline"
        )

    baseline = float(np.mean(series.values[inside]))
    anomaly = series.with_values(
        series.times,
        series.values - baseline,
        note=(
            f"anomaly relative to this series' own mean over [{start:g}, {end:g}] s "
            f"({baseline:.6g} {series.units}){f'; baseline: {label}' if label else ''}"
        ),
    )
    step = TransformationStep(
        operation=f"anomaly:{series.variable}",
        parameters=(
            f"baseline_window=[{start:g},{end:g}]; samples_in_window="
            f"{int(np.count_nonzero(inside))}; subtracted={baseline:.10g} {series.units}"
            + (f"; baseline_label={label}" if label else "")
        ),
        software="mfespf.evaluation.comparison.anomaly_against_baseline",
        uncertainty_effect=(
            "the baseline mean carries its own uncertainty, which is not propagated; a "
            "comparison of two anomalies is therefore insensitive to any constant offset "
            "between them and cannot test one (UNC-002)"
        ),
    )
    return anomaly, step


@dataclass(frozen=True)
class Comparison:
    """A model series and an observed one, on the same times, with the record of how."""

    variable: str
    times: np.ndarray
    modelled: np.ndarray
    observed: np.ndarray
    units: str
    observed_uncertainty: np.ndarray | None
    observed_uncertainty_kind: UncertaintyKind
    transformations: tuple[TransformationStep, ...]

    @property
    def n(self) -> int:
        return int(self.times.size)

    @property
    def residual(self) -> np.ndarray:
        return self.modelled - self.observed

    @property
    def rms_error(self) -> float:
        return float(np.sqrt(np.mean(self.residual**2)))

    @property
    def bias(self) -> float:
        return float(np.mean(self.residual))

    @property
    def observed_range(self) -> float:
        return float(np.ptp(self.observed))

    @property
    def modelled_range(self) -> float:
        return float(np.ptp(self.modelled))

    def explained_fraction(self) -> float:
        """Fraction of the observed variance the model accounts for.

        Negative when the model is worse than a flat line through the
        observations, which is a perfectly possible outcome and one worth being
        able to state plainly rather than reporting a correlation instead.
        """
        variance = float(np.var(self.observed))
        if variance == 0.0:
            raise ConfigurationError("the observed series does not vary; nothing to explain")
        return 1.0 - float(np.var(self.residual)) / variance

    def within_observational_uncertainty(self) -> float:
        """Fraction of points where the residual is inside the record's own error.

        Uses the record's uncertainty at face value, whatever kind it is, and the
        caller has to know which: a 95 percent half-width is a laxer test than a
        standard uncertainty by about a factor of two.
        """
        if self.observed_uncertainty is None:
            raise ConfigurationError(
                "the observed record carries no uncertainty, so there is nothing to compare "
                "the residual against"
            )
        return float(np.mean(np.abs(self.residual) <= self.observed_uncertainty))

    def describe(self) -> str:
        lines = [
            f"{self.variable}: {self.n} paired samples in {self.units}",
            f"  observed range  {self.observed_range:.4g}",
            f"  modelled range  {self.modelled_range:.4g}",
            f"  bias            {self.bias:+.4g}",
            f"  rms error       {self.rms_error:.4g}",
            f"  variance explained {self.explained_fraction():+.4f}",
        ]
        if self.observed_uncertainty is not None:
            lines.append(
                f"  within the record's own {self.observed_uncertainty_kind.value}: "
                f"{100.0 * self.within_observational_uncertainty():.1f} percent of points"
            )
        return "\n".join(lines)


def align_to_record(
    modelled: TimeSeries, observed: TimeSeries, *, method: str | None = None
) -> Comparison:
    """Resample the model onto the record's times over their overlap.

    The model moves to the data and never the other way round. Resampling an
    observational record onto a model's clock would make the comparison depend
    on a choice belonging to the model, and would smooth the record by an amount
    that varies with the model's time step.
    """
    if modelled.units != observed.units:
        raise ProvenanceError(
            f"cannot compare {modelled.variable!r} in {modelled.units!r} with "
            f"{observed.variable!r} in {observed.units!r}; convert at a declared boundary "
            f"(DAT-004)"
        )
    low = max(modelled.coverage[0], observed.coverage[0])
    high = min(modelled.coverage[1], observed.coverage[1])
    if high <= low:
        raise ProvenanceError(
            f"{modelled.variable!r} covers [{modelled.coverage[0]:g}, "
            f"{modelled.coverage[1]:g}] and {observed.variable!r} covers "
            f"[{observed.coverage[0]:g}, {observed.coverage[1]:g}]; they do not overlap"
        )

    inside = np.logical_and(observed.times >= low, observed.times <= high)
    targets = observed.times[inside]
    if targets.size < 2:
        raise ProvenanceError("the overlap contains fewer than two observations")

    resampled, report = resample(
        modelled, targets, method=method, target_interval=observed.median_interval
    )
    return Comparison(
        variable=observed.variable,
        times=targets,
        modelled=resampled.values,
        observed=observed.values[inside],
        units=observed.units,
        observed_uncertainty=(
            observed.uncertainty[inside] if observed.uncertainty is not None else None
        ),
        observed_uncertainty_kind=observed.uncertainty_kind,
        transformations=(
            report.as_transformation_step(),
            TransformationStep(
                operation=f"align:{observed.variable}",
                parameters=(
                    f"overlap=[{low:g},{high:g}]; paired={targets.size}; "
                    f"observations_outside_overlap={int(observed.times.size - targets.size)}; "
                    f"direction=model_onto_record"
                ),
                software="mfespf.evaluation.comparison.align_to_record",
                uncertainty_effect=(
                    "the record's uncertainties are carried through unchanged; the model "
                    "carries none, so the residual's uncertainty is the record's alone"
                ),
            ),
        ),
    )
