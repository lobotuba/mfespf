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
    areal_heat_capacity=1.0e6,
    master_dt=3600.0,
    initial_water_store=1.0e7,
)

#: Long enough that every configuration below reaches an imbalance well inside
#: the equilibrium threshold at the default 30 layers, which was measured rather
#: than guessed: the full column carries +0.016 W m^-2 at 24000 steps and +0.008
#: at 32000, with the time-mean moving 0.03 K between them.
STEPS = 24_000


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
    inert = equilibrate("30 layers: vapour present, radiatively inert", **BASE)
    radiative = equilibrate("30 layers: vapour radiatively active", radiative_vapour=True, **BASE)
    transported = equilibrate(
        "30 layers: vapour radiatively active and transported",
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

    # The full column at the default 30 layers, which the vapour band is now
    # calibrated at. Held to a wide band on purpose: this asserts that the bias
    # has not moved by more than a kelvin without anyone saying so, not that
    # +9.8 K is correct.
    assert transported.error_k == pytest.approx(9.8, abs=1.0)
    assert transported.distinguishable, "a bias inside the observational spread"

    # And the process quantities, because a temperature on its own has been
    # wrong here before while looking right. The excess evaporation is not an
    # independent error: Clausius-Clapeyron over a surface 9.8 K too warm is a
    # factor of about 1.9, and 4.5 against an observed 2.7 is 1.66. The Bowen
    # ratio is *below* the observed 0.26 at this resolution (0.18): it passes
    # through the observed value near 12 layers and keeps falling, so it is
    # bounded here, not claimed as an agreement.
    assert transported.evaporation_mm_day is not None
    assert transported.bowen_ratio is not None
    assert 2.5 < transported.evaporation_mm_day < 5.5
    assert 0.1 < transported.bowen_ratio < 0.3
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
        "30 layers: the same column, stopped early",
        radiative_vapour=True,
        vapour_transport=True,
        steps=2000,
        **BASE,
    )
    assert not early.equilibrated
    assert early.error_k is None
    assert "no error to report" in early.describe()
    assert abs(early.imbalance_w_m2) > OBSERVED_TEMPERATURE_SPREAD_K


#: Enough for 40 layers to settle with margin: at 36000 steps the imbalance is
#: 0.013 W m^-2 and the time-mean moves 0.01 K between a tail of a fifth and of
#: two fifths. At 24000 it carries 0.086, too close to the threshold to rely on.
RESOLVED_STEPS = 36_000


@pytest.mark.req("VAL-002", "KER-012", "VAL-005")
def test_the_headline_error_is_quoted_against_its_resolution_limit(equilibrate):
    """What the default resolution reports, and how far it is from the limit.

    The headline was once quoted at 8 layers with no layer study behind it. With
    the vapour band calibrated at the default 30 layers
    (``CALIBRATION_N_LAYERS``), the full column reads:

    ======  =========  =====
    layers  error (K)  Bowen
    ======  =========  =====
    8       +5.58      0.363
    20      +9.70      0.211
    30      +9.81      0.177
    40      +10.62     0.158
    60      +10.63     0.142
    ======  =========  =====

    Converged at +10.6 K: 40 and 60 layers agree to 0.02 K. The default reads
    0.8 K low, and not because it is creeping toward the limit -- 20 and 30
    agree, 40 and 60 agree, and there is a step between them nobody has
    explained yet. This test holds both facts: the limit, and the size of the
    step, so that it cannot grow or vanish without the test noticing.
    """
    coarse = equilibrate(
        "full water cycle, 8 layers (NOT resolution-converged)",
        radiative_vapour=True,
        vapour_transport=True,
        **{**BASE, "n_layers": 8},
    )
    default = equilibrate(
        "full water cycle, 30 layers (the default)",
        radiative_vapour=True,
        vapour_transport=True,
        **BASE,
    )
    limit = equilibrate(
        "full water cycle, 40 layers (resolution limit)",
        radiative_vapour=True,
        vapour_transport=True,
        steps=RESOLVED_STEPS,
        **{**BASE, "n_layers": 40},
    )
    for report in (coarse, default, limit):
        assert report.error_k is not None, report.describe()
    assert coarse.error_k is not None and default.error_k is not None
    assert limit.error_k is not None

    # The number to quote.
    assert limit.error_k == pytest.approx(10.6, abs=0.5)
    # The default is below it by the unexplained step, and by no more than a
    # kelvin: enough to state, not enough to change any conclusion drawn here.
    assert 0.3 < limit.error_k - default.error_k < 1.2
    # And 8 layers is not a small correction away from either.
    assert default.error_k - coarse.error_k > 3.0
