"""Sweeps, sensitivities and intervals that say what they mean.

Charter 7.6; requirements UNC-005, UNC-006, UNC-007.

The sensitivity test is the one that matters. A sweep that merely runs is a loop;
what makes it a measurement is that the derivative it produces can be checked
against a closed form, and the airless body supplies one -- its equilibrium
temperature is analytic in the irradiance, so ``dT/dS = T / 4S`` exactly, and a
finite difference over a sweep has to land on it.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError
from mfespf.data.records import UncertaintyKind
from mfespf.evaluation.uncertainty import (
    IntervalMeaning,
    Member,
    PredictiveInterval,
    UncertaintyMethod,
    bounds,
    ensemble,
    model_forms,
    monte_carlo,
    root_sum_square,
    run_ensemble,
    sweep,
)
from mfespf.physics.radiation import equilibrium_temperature
from mfespf.run.configurations import airless_radiating_body

pytestmark = pytest.mark.stage(1)


def surface_temperature(result) -> float:
    """The observable, named: ``run_ensemble`` takes the name from the function."""
    return float(np.asarray(result.steps[-1].outputs["surface_temperature"].values).sum())


# --- UNC-005: the sweep, and the derivative that checks it -----------------


@pytest.mark.req("UNC-005", "VAL-001")
def test_a_sweep_recovers_a_sensitivity_that_has_a_closed_form():
    """The check that makes the machinery a measurement rather than a loop.

    An airless body in radiative balance sits at ``T = (S(1-a)/4 e sigma)^(1/4)``,
    so ``dT/dS = T / 4S`` exactly. A secant over a sweep must land on that, and it
    must *miss* by a predictable amount -- which is the stronger half of the
    check, because landing close could be a coincidence and missing by the right
    term cannot.

    The direction is worth getting right, and getting it wrong here is what
    produced this test. ``T`` is concave in ``S``, so the first guess is that the
    secant falls below the tangent. It does not: the mean of ``f'`` over the
    interval is what a secant measures, ``f'`` is *convex* (``f'' < 0`` but
    ``f''' > 0`` for a fourth root), and by Jensen the mean of a convex function
    exceeds its midpoint value. The excess is ``h^2 f'''(S_0) / 6``, and that is
    what is asserted.
    """
    s_0 = 1361.0
    width = 40.0
    design = sweep(
        "total_irradiance",
        [s_0 - width, s_0, s_0 + width],
        basis="a +/-3 percent excursion about the reference irradiance, chosen to "
        "be wide enough to resolve the difference and narrow enough that the "
        "fourth-root curvature stays small",
    )
    assert design.method is UncertaintyMethod.PARAMETER_SWEEP

    result = run_ensemble(
        design,
        build=lambda **kw: airless_radiating_body(
            initial_temperature=equilibrium_temperature(kw["total_irradiance"], 0.3, 1.0), **kw
        ),
        observe=surface_temperature,
        n_steps=48,
    )

    assert result.observable == "surface_temperature"
    assert len(result.values) == 3
    t_mid = equilibrium_temperature(s_0, 0.3, 1.0)
    analytic = t_mid / (4.0 * s_0)
    assert result.sensitivity() == pytest.approx(analytic, rel=2e-3)

    # f = c S^(1/4), so f' = c S^(-3/4) / 4 = analytic and
    # f''' = (21/64) c S^(-11/4) = (21/64) (4 analytic) / S^2.
    third_derivative = (21.0 / 64.0) * (4.0 * analytic) / s_0**2
    assert result.sensitivity() - analytic == pytest.approx(
        width**2 * third_derivative / 6.0, rel=0.02
    )


@pytest.mark.req("UNC-005")
def test_the_six_methods_differ_only_in_how_the_members_are_chosen():
    """Which is the design, and is worth pinning so it is not quietly undone.

    Each constructor produces the same kind of object and declares a different
    method. That is what keeps "this was sampled from a distribution" a property
    of the ensemble rather than of which function happened to be called.
    """
    got = {
        bounds("alpha", 0.1, 0.5, basis="asserted").method: 2,
        sweep("alpha", [0.1, 0.3, 0.5], basis="asserted").method: 3,
        ensemble([{"a": 1, "b": 2}, {"a": 3, "b": 4}], basis="chosen pairs").method: 2,
        monte_carlo(
            {"alpha": (0.3, 0.05)},
            n_members=16,
            seed=7,
            basis="published spread",
            distributional_assumptions="normal, truncated nowhere",
        ).method: 16,
        model_forms(
            {"dry adiabat": {"convective_target": "dry"}, "plume": {"convective_target": "plume"}},
            basis="the two closures implemented",
        ).method: 2,
    }
    assert set(got) == set(UncertaintyMethod) - {UncertaintyMethod.SENSITIVITY_DERIVATIVE}
    assert got[UncertaintyMethod.MONTE_CARLO] == 16

    # A Monte Carlo is reproducible from its declaration, seed included.
    first = monte_carlo(
        {"alpha": (0.3, 0.05)},
        n_members=8,
        seed=11,
        basis="published spread",
        distributional_assumptions="normal",
    )
    again = monte_carlo(
        {"alpha": (0.3, 0.05)},
        n_members=8,
        seed=11,
        basis="published spread",
        distributional_assumptions="normal",
    )
    assert [m.overrides for m in first.members] == [m.overrides for m in again.members]


@pytest.mark.req("UNC-005", "UNC-002", "INF-003")
def test_an_ensemble_must_say_where_its_values_came_from():
    """A range asserted and a range measured produce identical numbers.

    The evaporating condensate fraction was swept over 0.1 to 0.5 on the
    reasoning that a stated range beats a fitted number, and published
    precipitation efficiencies put it near two thirds -- so the sweep excluded
    the answer, and nothing in its output said where the bounds had come from.
    """
    with pytest.raises(ConfigurationError, match="where its values came from"):
        sweep("alpha", [0.1, 0.5], basis="   ")
    with pytest.raises(ConfigurationError, match="at least two members"):
        sweep("alpha", [0.3], basis="a single point is not a spread")
    with pytest.raises(ConfigurationError, match="must be distinct"):
        ensemble([{"a": 1}, {"a": 1}], basis="duplicated by accident")

    # A sampling method has to declare the distribution it sampled, or its
    # coverage would be a number with nothing behind it.
    with pytest.raises(ConfigurationError, match="must declare them"):
        monte_carlo(
            {"alpha": (0.3, 0.05)},
            n_members=4,
            seed=1,
            basis="published spread",
            distributional_assumptions="unresolved",
        )


@pytest.mark.req("UNC-005", "UNC-003")
def test_a_derivative_is_refused_where_it_would_not_be_one():
    """A derivative with respect to "which configuration" is not a derivative.

    Model-form uncertainty stays separately traceable from parameter uncertainty
    (UNC-003), and this is where the two would otherwise get mixed.
    """
    forms = model_forms({"a": {"x": 1}, "b": {"x": 2}}, basis="two closures")
    result = run_ensemble(forms, build=lambda **_: _Stub(1.0), observe=_value, n_steps=1)
    with pytest.raises(ConfigurationError, match="needs a parameter sweep"):
        result.sensitivity()

    joint = ensemble([{"a": 1.0, "b": 2.0}, {"a": 2.0, "b": 1.0}], basis="pairs")
    with pytest.raises(ConfigurationError, match="needs a parameter sweep"):
        run_ensemble(joint, build=lambda **_: _Stub(1.0), observe=_value, n_steps=1).sensitivity()


# --- UNC-007: what the interval means ---------------------------------------


@pytest.mark.req("UNC-007")
def test_a_sweep_gives_a_range_and_refuses_to_call_it_a_coverage():
    """The refusal this module exists for.

    A three-point sweep across an asserted range tells you what the model does
    across that range. It says nothing about how likely any part of it is, so a
    number between 0 and 1 beside it would be read as a probability it has not
    earned -- the same false precision as quoting an assumed parameter's midpoint.
    """
    design = sweep("alpha", [0.1, 0.3, 0.5], basis="asserted, pending a literature value")
    result = run_ensemble(
        design, build=lambda **kw: _Stub(10.0 * kw["alpha"]), observe=_value, n_steps=1
    )

    interval = result.interval()
    assert interval.meaning is IntervalMeaning.FULL_RANGE
    assert interval.nominal_coverage is None
    assert "not drawn from a distribution" in interval.coverage_note
    assert (interval.low, interval.high) == pytest.approx((1.0, 5.0))
    assert "no nominal coverage" in interval.describe()

    with pytest.raises(ConfigurationError, match="cannot support a nominal coverage"):
        result.interval(0.95)


@pytest.mark.req("UNC-007")
def test_a_sampled_ensemble_may_state_a_coverage():
    """And refusing it one would be as wrong as granting one to a sweep.

    The distinction is about where the members came from, not about how many
    there are: these were drawn from a declared distribution, so a fraction of
    them really does lie inside the interval.
    """
    design = monte_carlo(
        {"alpha": (0.3, 0.05)},
        n_members=400,
        seed=3,
        basis="a published spread",
        distributional_assumptions="normal about the published central value",
    )
    result = run_ensemble(design, build=lambda **kw: _Stub(kw["alpha"]), observe=_value, n_steps=1)

    interval = result.interval(0.95)
    assert interval.meaning is IntervalMeaning.PERCENTILE_INTERVAL
    assert interval.nominal_coverage == 0.95
    assert "95% nominal coverage" in interval.describe()
    # Two sigma of a normal, to the accuracy 400 samples support.
    assert interval.width == pytest.approx(2 * 1.96 * 0.05, rel=0.15)
    assert result.interval().nominal_coverage is None  # still available, still honest


@pytest.mark.req("UNC-007")
def test_an_interval_without_a_coverage_must_say_why():
    """Otherwise the absence reads as an oversight rather than a statement."""
    with pytest.raises(ConfigurationError, match="must say why"):
        PredictiveInterval(
            low=1.0,
            high=2.0,
            meaning=IntervalMeaning.FULL_RANGE,
            nominal_coverage=None,
            method=UncertaintyMethod.PARAMETER_SWEEP,
        )
    with pytest.raises(ConfigurationError, match="not a probability"):
        PredictiveInterval(
            low=1.0,
            high=2.0,
            meaning=IntervalMeaning.PERCENTILE_INTERVAL,
            nominal_coverage=95.0,
            method=UncertaintyMethod.MONTE_CARLO,
        )
    with pytest.raises(ConfigurationError, match="inverted"):
        PredictiveInterval(
            low=2.0,
            high=1.0,
            meaning=IntervalMeaning.FULL_RANGE,
            nominal_coverage=None,
            method=UncertaintyMethod.BOUNDS,
            coverage_note="a range",
        )


# --- UNC-006: quadrature, and when it is not allowed ------------------------


@pytest.mark.req("UNC-006")
def test_quadrature_refuses_contributions_that_are_not_standard_uncertainties():
    """The factor of two that leaves no trace in the number.

    A 95 percent half-width is about two sigma. Passing one into a quadrature sum
    unconverted halves the combined figure, and nothing about the result looks
    wrong -- which is exactly why ``UncertaintyKind`` travels with the value.
    """
    good = {
        "radiative": (0.4, UncertaintyKind.STANDARD_UNCERTAINTY),
        "sampling": (0.3, UncertaintyKind.STANDARD_ERROR),
    }
    combined = root_sum_square(good, independent=True, justification="different instruments")
    assert combined == pytest.approx(0.5)  # 3-4-5

    bad = {**good, "record": (0.2, UncertaintyKind.HALF_WIDTH_95)}
    with pytest.raises(ConfigurationError, match="close to a factor of two"):
        root_sum_square(bad, independent=True, justification="different instruments")

    undeclared = {**good, "other": (0.2, UncertaintyKind.UNDECLARED)}
    with pytest.raises(ConfigurationError, match="not standard uncertainties"):
        root_sum_square(undeclared, independent=True, justification="different instruments")


@pytest.mark.req("UNC-006", "UNC-002")
def test_quadrature_refuses_correlated_contributions_and_bare_assertions():
    """Independence declared is the requirement; independence asserted is not."""
    contributions = {"a": (0.3, UncertaintyKind.STANDARD_UNCERTAINTY)}
    with pytest.raises(ConfigurationError, match="assumes independent"):
        root_sum_square(contributions, independent=False, justification="they are correlated")
    with pytest.raises(ConfigurationError, match="justified, not merely asserted"):
        root_sum_square(contributions, independent=True, justification="")
    with pytest.raises(ConfigurationError, match="cannot be negative"):
        root_sum_square(
            {"a": (-0.3, UncertaintyKind.STANDARD_UNCERTAINTY)},
            independent=True,
            justification="unrelated measurements",
        )


# --- stubs ------------------------------------------------------------------


class _Stub:
    """Stands in for a kernel where the point is the declaration, not the physics."""

    def __init__(self, value: float) -> None:
        self.value = value

    def run(self, n_steps: int) -> _Stub:
        return self


def _value(result: _Stub) -> float:
    return result.value


def test_a_member_needs_a_label():
    with pytest.raises(ConfigurationError, match="needs a label"):
        Member("  ")
