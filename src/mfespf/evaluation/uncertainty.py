"""Parameter sweeps, sensitivity derivatives and declared intervals.

Charter 7.6 and section 15; requirements UNC-005, UNC-006, UNC-007.

UNC-005 asks the framework to support bounds, sensitivity derivatives, parameter
sweeps, ensembles, Monte Carlo propagation and model-form comparison. There were
interfaces and no implementations, and the evidence that this mattered was
accidental: establishing what the bulk-plume closure's assumed parameter was
worth needed a three-point sweep at two resolutions, run by hand in scratch
scripts, and the answer -- a factor of 2.7 in the feedback and eleven kelvin of
surface temperature -- was the most decision-relevant number produced that day.
Four more sweeps have been run by hand since, for the evaporating fraction, the
timestep, the exchange coefficient and gravity.

**The six methods are one mechanism.** Every one of them runs a configuration
several times with something changed and collects an observable; they differ only
in *how the members are chosen*. So there is one ``Ensemble``, and six
constructors that build its members:

============================ ==================================================
``bounds``                   two members, at a declared range's ends
``sweep``                    members along one parameter
``ensemble``                 members at joint parameter assignments
``monte_carlo``              members sampled from declared distributions
``model_forms``              members that are different configurations
sensitivity derivative       a finite difference over a numeric ``sweep``
============================ ==================================================

Treating them as one thing is not a simplification; it is what they are. Writing
six mechanisms would have made the difference between an ensemble and a Monte
Carlo sample a matter of which function was called rather than a matter of how
the members were chosen, and the second is the only statement that carries any
information about what the spread means.

What the spread means
---------------------

Which is UNC-007, and the reason this module refuses more than it computes. A
three-point sweep across an *assumed* range gives a range. It does not give a 95
percent interval, and calling it one would be the false precision this project
has warned about since the evaporating condensate fraction was quoted at its
midpoint. So a ``PredictiveInterval`` carries its meaning and its nominal
coverage, and **a method that cannot support a coverage must state ``None`` and
say why** rather than picking a number that reads like a probability.

The distinction is load-bearing in both directions. A Monte Carlo sample from
declared distributions *can* state a coverage, and refusing it one would be as
wrong as granting one to a sweep.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from dataclasses import field as dc_field
from enum import StrEnum
from typing import Any

import numpy as np

from mfespf.core.errors import ConfigurationError
from mfespf.data.records import UncertaintyKind

__all__ = [
    "Ensemble",
    "EnsembleResult",
    "IntervalMeaning",
    "Member",
    "PredictiveInterval",
    "UncertaintyMethod",
    "bounds",
    "ensemble",
    "model_forms",
    "monte_carlo",
    "root_sum_square",
    "run_ensemble",
    "sweep",
]


class UncertaintyMethod(StrEnum):
    """UNC-005's six, named so a result says which one produced it."""

    BOUNDS = "bounds"
    SENSITIVITY_DERIVATIVE = "sensitivity_derivative"
    PARAMETER_SWEEP = "parameter_sweep"
    ENSEMBLE = "ensemble"
    MONTE_CARLO = "monte_carlo"
    MODEL_FORM_COMPARISON = "model_form_comparison"

    @property
    def supports_a_coverage(self) -> bool:
        """Whether a spread from this method can carry a nominal coverage.

        Only where the members were drawn from a declared distribution. A sweep
        across a range someone asserted tells you what the model does across that
        range and nothing about how likely any part of it is, however many points
        it has. An ensemble of joint assignments is in the same position unless
        the assignments were sampled, which is what makes it a Monte Carlo.
        """
        return self is UncertaintyMethod.MONTE_CARLO


class IntervalMeaning(StrEnum):
    """What the two numbers of an interval are (UNC-007).

    Parallel to ``UncertaintyKind``, which says what a *record's* uncertainty
    column means, and deliberately a separate vocabulary: a record's column is a
    measurement's uncertainty, and this is a prediction's. ``FULL_RANGE`` has no
    counterpart there because no dataset publishes one, and it is the honest
    answer for most of what this module produces.
    """

    FULL_RANGE = "full_range"
    STANDARD_DEVIATION = "standard_deviation"
    PERCENTILE_INTERVAL = "percentile_interval"


@dataclass(frozen=True, slots=True)
class Member:
    """One run of an ensemble: a label and what was changed to produce it."""

    label: str
    overrides: Mapping[str, Any] = dc_field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.label.strip():
            raise ConfigurationError("an ensemble member needs a label")


@dataclass(frozen=True, slots=True)
class Ensemble:
    """A set of runs to be made, and the declarations that make them readable.

    ``basis`` is where the values came from, and it is required for the same
    reason INF-003 requires a calibration to name its target: a range asserted
    from taste and a range read out of the literature produce identical numbers
    and mean different things. The evaporating condensate fraction was asserted
    over 0.1 to 0.5 when published precipitation efficiencies put it near two
    thirds, so every sweep across it was measuring the wrong interval -- and
    nothing in the sweep's output said where 0.1 to 0.5 had come from.
    """

    method: UncertaintyMethod
    members: tuple[Member, ...]
    parameters: tuple[str, ...]
    basis: str
    distributional_assumptions: str = "unresolved"
    correlations: str = "unresolved"

    def __post_init__(self) -> None:
        if len(self.members) < 2:
            raise ConfigurationError(
                f"an ensemble needs at least two members to be a spread; got {len(self.members)}"
            )
        labels = [m.label for m in self.members]
        if len(set(labels)) != len(labels):
            raise ConfigurationError(f"ensemble member labels must be distinct: {labels}")
        if not self.basis.strip():
            raise ConfigurationError(
                "an ensemble must declare where its values came from; a range asserted "
                "and a range measured produce identical numbers and mean different things"
            )
        if self.method.supports_a_coverage and self.distributional_assumptions == "unresolved":
            raise ConfigurationError(
                f"{self.method.value} draws from distributions, so it must declare them "
                "(UNC-002); an undeclared distribution makes the coverage meaningless"
            )


@dataclass(frozen=True, slots=True)
class PredictiveInterval:
    """An interval that states what it means and what it covers (UNC-007).

    ``nominal_coverage`` is ``None`` when the method cannot support one, and then
    ``coverage_note`` must say why. That is not a formality: the alternative is a
    number between 0 and 1 sitting beside a range, which every reader will take
    for a probability.
    """

    low: float
    high: float
    meaning: IntervalMeaning
    nominal_coverage: float | None
    method: UncertaintyMethod
    coverage_note: str = ""

    def __post_init__(self) -> None:
        if self.high < self.low:
            raise ConfigurationError(f"interval [{self.low}, {self.high}] is inverted")
        if self.nominal_coverage is None:
            if not self.coverage_note.strip():
                raise ConfigurationError(
                    "an interval with no nominal coverage must say why it has none "
                    "(UNC-007); otherwise the absence reads as an oversight"
                )
        else:
            if not 0.0 < self.nominal_coverage <= 1.0:
                raise ConfigurationError(
                    f"nominal coverage {self.nominal_coverage} is not a probability"
                )
            if not self.method.supports_a_coverage:
                raise ConfigurationError(
                    f"{self.method.value} cannot support a nominal coverage: its members "
                    "were not drawn from a declared distribution, so the spread says what "
                    "the model does across the range and nothing about how likely any "
                    "part of it is"
                )

    @property
    def width(self) -> float:
        return self.high - self.low

    def describe(self) -> str:
        coverage = (
            f"{self.nominal_coverage:.0%} nominal coverage"
            if self.nominal_coverage is not None
            else f"no nominal coverage ({self.coverage_note})"
        )
        return f"[{self.low:.6g}, {self.high:.6g}] {self.meaning.value}, {coverage}"


@dataclass(frozen=True, slots=True)
class EnsembleResult:
    """One observable, over every member of an ensemble."""

    design: Ensemble
    observable: str
    values: tuple[float, ...]

    @property
    def labels(self) -> tuple[str, ...]:
        return tuple(m.label for m in self.design.members)

    @property
    def span(self) -> float:
        return max(self.values) - min(self.values)

    @property
    def mean(self) -> float:
        return float(np.mean(self.values))

    def interval(self, coverage: float | None = None) -> PredictiveInterval:
        """The spread, as an interval that declares what it is.

        With no ``coverage`` this is the full range of what was run, which is the
        honest reading of a sweep and is what the constructor will accept from
        one. Asking a non-sampling method for a coverage is refused rather than
        approximated.
        """
        if coverage is None:
            return PredictiveInterval(
                low=min(self.values),
                high=max(self.values),
                meaning=IntervalMeaning.FULL_RANGE,
                nominal_coverage=None,
                method=self.design.method,
                coverage_note=(
                    f"a {self.design.method.value} spans the values it was given; its "
                    "members were not drawn from a distribution, so no fraction of "
                    "anything lies inside this range"
                ),
            )
        quantiles = np.quantile(self.values, [0.5 * (1.0 - coverage), 0.5 * (1.0 + coverage)])
        return PredictiveInterval(
            low=float(quantiles[0]),
            high=float(quantiles[1]),
            meaning=IntervalMeaning.PERCENTILE_INTERVAL,
            nominal_coverage=coverage,
            method=self.design.method,
        )

    def sensitivity(self) -> float:
        """d(observable)/d(parameter), by central difference over the sweep.

        Only defined for a sweep along one numeric parameter: a derivative with
        respect to "which configuration was used" is not a derivative. Computed
        over the whole span rather than between adjacent points, because a
        three-point sweep across a wide asserted range is the usual case here and
        the endpoint slope is what it supports.
        """
        if self.design.method is not UncertaintyMethod.PARAMETER_SWEEP:
            raise ConfigurationError(
                f"a sensitivity derivative needs a parameter sweep, not a "
                f"{self.design.method.value}"
            )
        if len(self.design.parameters) != 1:
            raise ConfigurationError(
                f"a sensitivity derivative needs one parameter, not {self.design.parameters}"
            )
        name = self.design.parameters[0]
        x = np.array([float(m.overrides[name]) for m in self.design.members], dtype=float)
        y = np.array(self.values, dtype=float)
        order = np.argsort(x)
        x, y = x[order], y[order]
        if x[-1] == x[0]:
            raise ConfigurationError("the swept parameter did not vary")
        return float((y[-1] - y[0]) / (x[-1] - x[0]))

    def describe(self) -> str:
        lines = [
            f"{self.observable} over a {self.design.method.value} "
            f"of {len(self.values)} members ({self.design.basis})"
        ]
        for label, value in zip(self.labels, self.values, strict=True):
            lines.append(f"  {label:32s} {value:.6g}")
        lines.append(f"  span {self.span:.6g}")
        return "\n".join(lines)


# --- how the members get chosen --------------------------------------------


def bounds(parameter: str, low: float, high: float, *, basis: str) -> Ensemble:
    """Two members, at the ends of a declared range."""
    if not high > low:
        raise ConfigurationError(f"bounds [{low}, {high}] must be increasing")
    return Ensemble(
        method=UncertaintyMethod.BOUNDS,
        members=(
            Member(f"{parameter}={low:g}", {parameter: low}),
            Member(f"{parameter}={high:g}", {parameter: high}),
        ),
        parameters=(parameter,),
        basis=basis,
    )


def sweep(parameter: str, values: Sequence[float], *, basis: str) -> Ensemble:
    """Members along one parameter, which is what a sensitivity is taken over."""
    return Ensemble(
        method=UncertaintyMethod.PARAMETER_SWEEP,
        members=tuple(Member(f"{parameter}={v:g}", {parameter: v}) for v in values),
        parameters=(parameter,),
        basis=basis,
    )


def ensemble(assignments: Sequence[Mapping[str, Any]], *, basis: str) -> Ensemble:
    """Members at joint parameter assignments, chosen rather than sampled."""
    names = tuple(sorted({k for a in assignments for k in a}))
    return Ensemble(
        method=UncertaintyMethod.ENSEMBLE,
        members=tuple(
            Member(", ".join(f"{k}={a[k]!r}" for k in sorted(a)), dict(a)) for a in assignments
        ),
        parameters=names,
        basis=basis,
    )


def monte_carlo(
    distributions: Mapping[str, tuple[float, float]],
    *,
    n_members: int,
    seed: int,
    basis: str,
    distributional_assumptions: str,
    correlations: str = "assumed independent",
) -> Ensemble:
    """Members sampled from declared normal distributions, as ``(mean, sigma)``.

    The seed is required rather than defaulted. An unseeded sample makes two runs
    of the same declared ensemble produce different numbers with identical
    manifests, which is the reproducibility failure the whole audit trail exists
    to prevent.
    """
    if n_members < 2:
        raise ConfigurationError("a Monte Carlo sample needs at least two members")
    rng = np.random.default_rng(seed)
    names = tuple(sorted(distributions))
    draws = {
        name: rng.normal(distributions[name][0], distributions[name][1], n_members)
        for name in names
    }
    for name, (_mean, sigma) in distributions.items():
        if sigma <= 0.0:
            raise ConfigurationError(f"{name} has a non-positive standard deviation {sigma}")
    members = tuple(
        Member(f"sample {i:04d}", {name: float(draws[name][i]) for name in names})
        for i in range(n_members)
    )
    return Ensemble(
        method=UncertaintyMethod.MONTE_CARLO,
        members=members,
        parameters=names,
        basis=f"{basis} (seed {seed})",
        distributional_assumptions=distributional_assumptions,
        correlations=correlations,
    )


def model_forms(variants: Mapping[str, Mapping[str, Any]], *, basis: str) -> Ensemble:
    """Members that are different configurations rather than different numbers.

    Structural uncertainty, which UNC-003 requires to stay separately traceable
    from parameter uncertainty -- hence its own method rather than an ensemble
    with a discrete parameter.
    """
    return Ensemble(
        method=UncertaintyMethod.MODEL_FORM_COMPARISON,
        members=tuple(Member(label, dict(overrides)) for label, overrides in variants.items()),
        parameters=tuple(sorted({k for o in variants.values() for k in o})),
        basis=basis,
    )


# --- running one ------------------------------------------------------------


def run_ensemble(
    design: Ensemble,
    *,
    build: Callable[..., Any],
    observe: Callable[[Any], float],
    n_steps: int,
) -> EnsembleResult:
    """Run every member and collect one observable.

    ``build`` takes the member's overrides as keyword arguments and returns
    something with a ``run(n_steps)``; ``observe`` reduces that run to a number.
    Deliberately thin -- what this module is for is the declarations around the
    loop, not the loop.
    """
    values = []
    for member in design.members:
        result = build(**member.overrides).run(n_steps)
        values.append(float(observe(result)))
    return EnsembleResult(design=design, observable=observe.__name__, values=tuple(values))


# --- UNC-006 ----------------------------------------------------------------


def root_sum_square(
    contributions: Mapping[str, tuple[float, UncertaintyKind]],
    *,
    independent: bool,
    justification: str,
) -> float:
    """Combine uncertainties in quadrature, refusing when its assumptions are not met.

    UNC-006: root-sum-square is only valid when the contributions are independent
    and each is a standard uncertainty. Both conditions are checked rather than
    assumed, because both fail silently and in the same direction -- understating
    the combined figure -- and neither leaves a trace in the number.

    A 95 percent half-width is roughly two sigma, so passing one in unconverted
    halves the combined uncertainty of anything it enters. ``UncertaintyKind``
    exists because real records disagree about which they publish, and this is
    the place that disagreement would do its damage.
    """
    if not independent:
        raise ConfigurationError(
            "root-sum-square assumes independent contributions; correlated ones need "
            "their covariance, and combining them in quadrature understates the total "
            "(UNC-006). Declare the covariance or combine them another way."
        )
    if not justification.strip():
        raise ConfigurationError(
            "independence must be justified, not merely asserted (UNC-006, UNC-002)"
        )
    if not contributions:
        raise ConfigurationError("nothing to combine")

    wrong_kind = {
        name: kind.value for name, (_value, kind) in contributions.items() if not kind.is_one_sigma
    }
    if wrong_kind:
        raise ConfigurationError(
            f"these contributions are not standard uncertainties and cannot enter a "
            f"quadrature sum as they stand: {wrong_kind}. A 95 percent half-width is "
            "about two sigma, so including one unconverted understates the combined "
            "figure by close to a factor of two (UNC-006)."
        )
    negative = {name: v for name, (v, _k) in contributions.items() if v < 0.0}
    if negative:
        raise ConfigurationError(f"an uncertainty cannot be negative: {negative}")
    return math.sqrt(sum(value * value for value, _kind in contributions.values()))
