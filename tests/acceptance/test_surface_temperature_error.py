"""How far this column is from the observed surface temperature, in kelvin.

The number the project is judged on. It had no test: it was recomputed by hand in
a scratch script whenever anyone wanted it, so it moved between commits without
anyone noticing and the figure quoted in the documentation was whatever the last
person to run a script had seen.

The configurations are run in the order they gain physics, because the shape of
the sequence is the result and no single number is. Each one is an equilibrium --
the report refuses to state an error otherwise -- and each is a time-mean over
the tail of the run, because the column limit-cycles and an endpoint is a point
on a cycle.
"""

from __future__ import annotations

import pytest

from mfespf.evaluation.bias import OBSERVED_TEMPERATURE_SPREAD_K
from mfespf.run.column_configurations import column_grid, moist_radiative_convective

pytestmark = pytest.mark.stage(3)

GRID = column_grid(16)
BASE = dict(
    grid=GRID,
    n_layers=8,
    areal_heat_capacity=1.0e6,
    master_dt=3600.0,
    initial_water_store=1.0e7,
)

#: Long enough that every configuration below reaches an imbalance well inside
#: the equilibrium threshold, which was measured rather than guessed: the full
#: column carries +5.94 W m^-2 at 3000 steps, +1.44 at 6000 and +0.096 at 12000.
STEPS = 16_000


@pytest.fixture
def equilibrate(report_bias):
    """Run one configuration to equilibrium and report its error in kelvin."""

    def run(label: str, *, steps: int = STEPS, **kwargs):
        result = moist_radiative_convective(**kwargs).run(steps)
        return report_bias(result, label=label)

    return run


@pytest.mark.req("KER-012", "VAL-005", "CLM-003")
def test_the_water_is_the_whole_warm_bias(equilibrate):
    """The headline error, and the sequence that says where it comes from.

    Three configurations, differing by one mechanism each:

    * vapour present but radiatively inert -- a water cycle with no greenhouse
      contribution from its own water;
    * vapour radiatively active -- the absorption switched on, transport still off;
    * vapour transported -- convection carries moisture upward.

    The dry column is *cold*. Every step toward a complete water cycle warms it,
    which is the thing about this bias that is easy to get backwards: the vapour
    is a greenhouse gas and its absorption beats its own evaporative cooling. So
    this is not a missing cooling term to be found, and a change that makes the
    surface temperature worse on the way to a more complete water cycle is not a
    regression.
    """
    inert = equilibrate("vapour present, radiatively inert", **BASE)
    radiative = equilibrate("vapour radiatively active", radiative_vapour=True, **BASE)
    transported = equilibrate(
        "vapour radiatively active and transported",
        radiative_vapour=True,
        vapour_transport=True,
        **BASE,
    )

    for report in (inert, radiative, transported):
        assert report.equilibrated, report.describe()

    # Each mechanism warms the column, and the ordering is the claim.
    assert inert.error_k < radiative.error_k < transported.error_k

    # The dry-side configurations are too cold and the complete one too warm, so
    # the bias changes sign along the sequence. A configuration that happened to
    # land near zero in the middle of that would be passing through, not right.
    assert inert.error_k < 0.0 < transported.error_k

    # The full column, which is the number to quote. Held to a wide band on
    # purpose: this asserts that the bias has not moved by more than a kelvin
    # without anyone saying so, not that +6.7 K is correct.
    assert transported.error_k == pytest.approx(6.7, abs=1.0)
    assert transported.distinguishable, "a bias inside the observational spread"

    # And the process quantities, because a temperature on its own has been
    # wrong here before while looking right. The excess evaporation is not an
    # independent error: Clausius-Clapeyron over a surface 6.7 K too warm is a
    # factor of about 1.6, and 4.0 against an observed 2.7 is 1.48.
    assert transported.evaporation_mm_day is not None
    assert transported.bowen_ratio is not None
    assert 2.5 < transported.evaporation_mm_day < 5.5
    assert 0.2 < transported.bowen_ratio < 0.5
    assert transported.precipitable_water_kg_m2 is not None
    assert 20.0 < transported.precipitable_water_kg_m2 < 35.0


@pytest.mark.req("VAL-005")
def test_an_unfinished_run_reports_no_error_rather_than_a_small_one(equilibrate):
    """The guard that makes every other number here mean something.

    A bias read from a run still carrying an imbalance is not an uncertain
    result, it is a different quantity: the dry radiative-convective column at a
    large heat capacity reads -16.3 K at 20000 steps while carrying +43.8 W m^-2,
    and it is unfinished rather than cold. Deep-ocean configurations have looked
    3 K cooler for the same reason.
    """
    early = equilibrate(
        "the same column, stopped early",
        radiative_vapour=True,
        vapour_transport=True,
        steps=2000,
        **BASE,
    )
    assert not early.equilibrated
    assert early.error_k is None
    assert "no error to report" in early.describe()
    assert abs(early.imbalance_w_m2) > OBSERVED_TEMPERATURE_SPREAD_K
