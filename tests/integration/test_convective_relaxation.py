"""The relaxed convective adjustment: its limits, its fraction, its budgets.

The relaxation applies ``1 - exp(-dt / tau)`` of the instantaneous adjustment
each step. That is chosen to keep every property the instantaneous closure had
-- exact per-step conservation above all -- and these tests check that it does,
against independent statements rather than against itself: the fraction against
its closed form, the limit against the unrelaxed scheme run separately, and the
budgets against the conservation auditor.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from mfespf.physics.water import VAPOUR_STATE
from mfespf.run.column_configurations import column_grid, moist_radiative_convective

pytestmark = pytest.mark.stage(3)

COLUMN = dict(
    grid=column_grid(16),
    n_layers=12,
    areal_heat_capacity=1.0e6,
    master_dt=3600.0,
    initial_water_store=1.0e7,
    radiative_vapour=True,
    vapour_transport=True,
)


#: A tendency below this is round-off, not mixing: 1e-9 K/s is 3.6e-6 K per step.
SIGNIFICANT = 1.0e-9


def _tendencies(result) -> list[np.ndarray]:
    return [
        np.asarray(step.outputs["convective_temperature_tendency"].values, dtype=float)
        for step in result.steps
    ]


@pytest.mark.req("EQR-002", "VAL-001")
def test_a_zero_timescale_is_the_instantaneous_adjustment_exactly():
    """The limit, checked against the unrelaxed scheme run on its own."""
    instantaneous = moist_radiative_convective(**COLUMN).run(300)
    vanishing = moist_radiative_convective(convective_relaxation_timescale=1.0e-9, **COLUMN).run(
        300
    )
    a = instantaneous.steps[-1].outputs["surface_temperature"].values
    b = vanishing.steps[-1].outputs["surface_temperature"].values
    assert float(b) == pytest.approx(float(a), abs=1e-9)
    assert any(np.any(t != 0.0) for t in _tendencies(instantaneous)), "nothing was adjusted"


@pytest.mark.req("EQR-002", "VAL-001")
def test_the_relaxed_adjustment_applies_the_exact_exponential_fraction():
    """At the first step anything mixes, the tendency is the closed-form share.

    Until then the two runs are the same run -- no mixing, identical states --
    so the first step where the instantaneous scheme adjusts is a step where
    both start from the same provisional profile, and their tendencies must
    stand in the ratio ``1 - exp(-dt / tau)`` layer by layer.
    """
    tau = 7200.0
    probe = moist_radiative_convective(**COLUMN).run(400)
    # A real adjustment, not the 1e-17 K/s round-off a group already at
    # neutrality leaves behind.
    first = next(k for k, t in enumerate(_tendencies(probe)) if np.max(np.abs(t)) > SIGNIFICANT)

    full = _tendencies(moist_radiative_convective(**COLUMN).run(first + 1))[first]
    relaxed = _tendencies(
        moist_radiative_convective(convective_relaxation_timescale=tau, **COLUMN).run(first + 1)
    )[first]

    expected = 1.0 - math.exp(-3600.0 / tau)
    mixed = np.abs(full) > SIGNIFICANT
    assert np.any(mixed)
    assert np.allclose(relaxed[mixed] / full[mixed], expected, rtol=1e-9)
    assert np.all(np.abs(relaxed[~mixed]) <= SIGNIFICANT)


@pytest.mark.req("CON-001", "CON-004", "EQR-002")
def test_the_relaxed_adjustment_conserves_energy_and_water_every_step():
    """Scaling each group's adjustment by one number keeps its sum at zero.

    So the budgets close per step to the same tolerance as the instantaneous
    scheme's, and the vapour stays non-negative because the relaxed profile is
    a convex combination of two non-negative ones.
    """
    result = moist_radiative_convective(
        convective_relaxation_timescale=7200.0, **{**COLUMN, "n_layers": 30}
    ).run(1500)
    assert result.may_support_a_claim
    assert result.worst_relative_residual() < 1e-10
    assert all(
        np.all(np.asarray(step.outputs[VAPOUR_STATE].values) >= 0.0) for step in result.steps
    )
    assert any(np.any(t != 0.0) for t in _tendencies(result))
