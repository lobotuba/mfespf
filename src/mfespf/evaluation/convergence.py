"""Numerical convergence gates.

Charter 15.2 and 15.7; requirements VAL-001, VAL-002, VAL-003, UNC-003, ACC-007.

A result that changes when the time step is halved is a property of the solver,
not of the physics.  The charter treats that as a numerical failure, and this
module is the gate: run the same configuration at successively finer steps to
the same end time, and require the answer to settle.

The gate reports an observed order of convergence as well as a verdict.  That
number is the useful part -- a first-order scheme that measures first order is
working, while one that measures zeroth order is not converging at all, and the
difference is invisible if only the final value is compared.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from dataclasses import field as dc_field

from mfespf.core.errors import ConfigurationError

__all__ = ["ConvergenceReport", "time_step_convergence"]


@dataclass
class ConvergenceReport:
    """The record of a time-step refinement study."""

    time_steps: list[float] = dc_field(default_factory=list)
    values: list[float] = dc_field(default_factory=list)
    differences: list[float] = dc_field(default_factory=list)
    observed_order: float | None = None
    expected_order: float | None = None
    tolerance: float = 1e-3
    converged: bool = False
    diagnosis: str = ""

    @property
    def richardson_estimate(self) -> float | None:
        """Extrapolated value at zero step size, when an order was measurable.

        Richardson extrapolation on the two finest solutions.  It is reported as
        an estimate of the exact solution, and is the natural quantity against
        which to state the discretisation component of the uncertainty budget
        (UNC-003) rather than a value to use in place of the run.
        """
        if self.observed_order is None or len(self.values) < 2:
            return None
        ratio = self.time_steps[-2] / self.time_steps[-1]
        denominator = ratio**self.observed_order - 1.0
        if abs(denominator) < 1e-12:
            return None
        return self.values[-1] + (self.values[-1] - self.values[-2]) / denominator

    def describe(self) -> str:
        lines = ["time-step convergence study"]
        for dt, value in zip(self.time_steps, self.values, strict=True):
            lines.append(f"  dt = {dt:12.6g}  value = {value:.10g}")
        for dt, diff in zip(self.time_steps[1:], self.differences, strict=True):
            lines.append(f"  |value(dt) - value(2dt)| at dt = {dt:.6g}: {diff:.6g}")
        if self.observed_order is not None:
            lines.append(f"  observed order: {self.observed_order:.3f}")
            if self.expected_order is not None:
                lines.append(f"  expected order: {self.expected_order:.3f}")
            estimate = self.richardson_estimate
            if estimate is not None:
                lines.append(f"  Richardson estimate at dt -> 0: {estimate:.10g}")
        lines.append("  PASS" if self.converged else f"  FAIL: {self.diagnosis}")
        return "\n".join(lines)


def time_step_convergence(
    run_to_end_time: Callable[[float], float],
    *,
    base_dt: float,
    refinements: int = 3,
    tolerance: float = 1e-3,
    expected_order: float | None = None,
    scale: float | None = None,
) -> ConvergenceReport:
    """Refine the time step and measure whether the answer settles.

    ``run_to_end_time`` takes a time step and returns one scalar observable, run
    to the *same physical end time* each call -- halving the step means doubling
    the number of steps, and it is the caller's job to do that, because only the
    caller knows what the configuration's end time means.

    The verdict is about ``base_dt``: the question a run needs answered is not
    whether some finer step would converge, but whether *its own* step gives an
    answer the refinement agrees with.  So convergence requires the base value to
    sit within ``tolerance`` of the finest, and the successive differences to
    shrink.  A marginally unstable configuration fails the first test -- its
    coarse answer oscillates around the right one while every individual number
    looks plausible -- and a non-converging scheme fails the second.

    ``tolerance`` is relative to ``scale`` when given, else to the finest value.
    Differences at the round-off floor are treated as converged rather than
    ranked, since their ratios are noise.
    """
    if refinements < 2:
        raise ConfigurationError("a convergence study needs at least two refinements")

    report = ConvergenceReport(tolerance=tolerance, expected_order=expected_order)
    for level in range(refinements):
        dt = base_dt / (2**level)
        report.time_steps.append(dt)
        report.values.append(float(run_to_end_time(dt)))

    report.differences = [
        abs(later - earlier)
        for earlier, later in zip(report.values, report.values[1:], strict=False)
    ]

    reference = scale if scale is not None else abs(report.values[-1])
    reference = reference if reference > 0 else 1.0

    if any(not math.isfinite(v) for v in report.values):
        report.diagnosis = (
            "at least one refinement produced a non-finite result; the configuration left "
            "its numerical validity domain (charter 15.2)"
        )
        return report

    relative = [d / reference for d in report.differences]
    # Below this, successive values differ only by floating-point round-off and
    # their ratios carry no information about the order of the scheme.
    round_off_floor = 1e-11
    significant = [d for d in relative if d > round_off_floor]

    if not significant:
        report.observed_order = float("inf")
        report.converged = True
        return report

    if len(report.differences) >= 2 and report.differences[-1] > 0:
        ratio = report.differences[-2] / report.differences[-1]
        if ratio > 0:
            report.observed_order = math.log2(ratio)

    base_error = abs(report.values[0] - report.values[-1]) / reference
    shrinking = all(
        later <= earlier * (1.0 + 1e-9)
        for earlier, later in zip(report.differences, report.differences[1:], strict=False)
    )

    if not shrinking:
        report.diagnosis = (
            "successive differences do not shrink under refinement; the solution is not "
            "converging and the result is a property of the time step (charter 15.2)"
        )
    elif base_error > tolerance:
        report.diagnosis = (
            f"the answer at dt = {report.time_steps[0]:.6g} differs from the refined limit by "
            f"{base_error:.3e} relative, above the declared tolerance {tolerance:.3e}; this "
            f"time step is outside the configuration's numerical validity domain"
        )
    elif (
        expected_order is not None
        and report.observed_order is not None
        and math.isfinite(report.observed_order)
        and report.observed_order < expected_order * 0.5
    ):
        report.diagnosis = (
            f"observed order {report.observed_order:.2f} is far below the declared "
            f"{expected_order:.2f}; the scheme is not behaving as documented"
        )
    else:
        report.converged = True

    return report
