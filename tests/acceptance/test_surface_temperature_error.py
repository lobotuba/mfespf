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

import numpy as np
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


#: Enough for 30 and 40 layers to settle with margin: at 28000 steps both carry
#: under 0.02 W m^-2 and the time-means move under 0.03 K with the tail length.
RESOLVED_STEPS = 28_000


def _mean_convective_top_hpa(result, tail_fraction: float = 0.2) -> float:
    """Time-mean pressure of the highest layer the adjustment mixed, in hPa."""
    tops = []
    for step in result.steps[-int(len(result.steps) * tail_fraction) :]:
        mixed = np.asarray(step.outputs["convectively_mixed_layer_temperature"].values)
        pressure = np.asarray(step.outputs["atmosphere_layer_pressure"].values)
        levels = np.flatnonzero(~np.isnan(mixed))
        if levels.size:
            tops.append(pressure[levels.max()])
    return float(np.mean(tops)) / 100.0


@pytest.mark.req("VAL-002", "KER-012", "VAL-005")
def test_the_headline_error_is_quoted_with_its_resolution_scatter(report_bias):
    """The error does not converge with layer count; it scatters, and why.

    With the vapour band calibrated at the default 30 layers
    (``CALIBRATION_N_LAYERS``), the full column reads:

    ======  =========  =====================
    layers  error (K)  top of mixing (hPa)
    ======  =========  =====================
    8       +5.58
    20      +9.70
    25      +9.89      462
    30      +9.81      446
    32      +10.53     400
    36      +10.23     413
    40      +10.64     391
    45      +10.16     415
    60      +10.63
    ======  =========  =====================

    From 25 layers up there is no trend -- the error scatters from +9.8 to
    +10.6 K around +10.3 -- and it follows the convective top almost exactly:
    the deeper the mixing reaches, the warmer the surface, about 0.015 K per hPa,
    at an outgoing longwave of 240.0 W m^-2 in every case. The adjustment mixes
    whole layers or none, so its top snaps to a layer edge, and the edges fall
    somewhere different at each resolution. The scatter is the all-or-nothing
    closure, seen through the grid.

    It was first read, from 30, 40 and 60 layers alone, as a step between two
    plateaus. More points showed a scatter, and the diagnostic showed its cause.
    """
    base = dict(radiative_vapour=True, vapour_transport=True, **BASE)
    coarse = report_bias(
        moist_radiative_convective(**{**base, "n_layers": 8}).run(STEPS),
        label="full water cycle, 8 layers (NOT resolution-converged)",
    )
    runs = {
        n: moist_radiative_convective(**{**base, "n_layers": n}).run(RESOLVED_STEPS)
        for n in (30, 40)
    }
    reports = {
        n: report_bias(result, label=f"full water cycle, {n} layers") for n, result in runs.items()
    }
    errors = {n: r.error_k for n, r in reports.items()}
    assert coarse.error_k is not None
    assert all(e is not None for e in errors.values()), [r.describe() for r in reports.values()]
    e30, e40 = errors[30], errors[40]
    assert e30 is not None and e40 is not None

    # Both inside the scatter band the layer study measured, about +10.3 K.
    for error in (e30, e40):
        assert error == pytest.approx(10.3, abs=0.6)
    # The scatter is bounded: resolution moves the answer by under a kelvin.
    assert abs(e40 - e30) < 1.2
    # And its cause: whichever resolution mixes deeper is the warmer one.
    tops = {n: _mean_convective_top_hpa(result) for n, result in runs.items()}
    deeper = min(tops, key=lambda n: tops[n])
    warmer = max(errors, key=lambda n: errors[n] or 0.0)
    assert deeper == warmer, (tops, errors)
    # 8 layers is not a small correction away from any of it.
    assert min(e30, e40) - coarse.error_k > 3.0
