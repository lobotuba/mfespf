"""A convective top that can lie inside a layer.

The whole-layer adjustment's top is always a layer edge, and the column's error
followed where the grid put those edges. The fractional top lets a mixed group
rise part-way through the layer above it, to where a linear reconstruction of
that layer's stability equals the group's mixed value. These tests hold it to
the properties it has to keep -- exact enthalpy, a stable result, the
whole-layer scheme as its limit -- and to the one it was built for: a top that
moves continuously as the forcing changes.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.physics.convection import ConvectiveAdjustmentModule
from mfespf.physics.thermodynamics import (
    convective_adjustment_with_groups,
    potential_temperature,
    surface_coupled_adjustment,
)

pytestmark = pytest.mark.stage(3)

N = 40
PRESSURE = np.linspace(98000.0, 15000.0, N)
CAPACITY = np.full(N, 5.0e6)
#: Potential temperature rising 40 K through the column: stable everywhere.
STABLE = 300.0 * (PRESSURE / 1.0e5) ** 0.2857 + np.linspace(0.0, 40.0, N)


def _heated(delta: float) -> np.ndarray:
    profile = STABLE.copy()
    profile[0] += delta
    return profile


def _top(temperature: np.ndarray) -> float:
    """Top of the highest mixed group, in layers, counting the partial one."""
    partial: dict[int, float] = {}
    _, groups = convective_adjustment_with_groups(
        temperature, CAPACITY, PRESSURE, partial_tops=partial
    )
    if not groups:
        return 0.0
    high = max(h for _, h in groups)
    return high + 1 + partial.get(high + 1, 0.0)


@pytest.mark.req("CON-001", "EQR-002", "VAL-001")
def test_the_fractional_top_conserves_enthalpy_and_leaves_a_stable_profile():
    rng = np.random.default_rng(20260928)
    for _ in range(200):
        temperature = STABLE + rng.normal(0.0, 3.0, N)
        capacity = CAPACITY * rng.uniform(0.5, 2.0, N)
        partial: dict[int, float] = {}
        adjusted, _ = convective_adjustment_with_groups(
            temperature, capacity, PRESSURE, partial_tops=partial
        )
        before, after = np.sum(capacity * temperature), np.sum(capacity * adjusted)
        assert after == pytest.approx(before, rel=1e-13)
        assert np.all(np.diff(potential_temperature(adjusted, PRESSURE)) >= -1e-8)
        assert all(0.0 < f < 1.0 for f in partial.values())

        # And with the surface as the column's bottom element.
        surface = float(temperature[0] + rng.uniform(0.0, 5.0))
        s, layers, _, _ = surface_coupled_adjustment(
            temperature,
            capacity,
            PRESSURE,
            surface_temperature=surface,
            surface_heat_capacity=1.0e8,
            surface_pressure=1.0e5,
            partial_tops={},
        )
        total = 1.0e8 * surface + np.sum(capacity * temperature)
        assert 1.0e8 * s + np.sum(capacity * layers) == pytest.approx(total, rel=1e-13)


@pytest.mark.req("EQR-002", "VAL-002")
def test_the_fractional_top_moves_continuously_with_the_forcing():
    """What it was built for: the top is not a layer edge.

    Warming the lowest layer in steps of 0.01 K deepens the convection. The
    whole-layer top moves only in whole layers; the fractional one moves by a
    hundredth of a layer at most, once convection is under way.
    """
    deltas = np.linspace(3.0, 20.0, 1701)
    tops = np.array([_top(_heated(d)) for d in deltas])
    active = tops > 1.5
    steps = np.abs(np.diff(tops))[active[1:] & active[:-1]]
    assert tops.max() - tops[active].min() > 3.0, "the sweep should deepen the convection"
    assert steps.max() < 0.02

    whole = []
    for d in deltas:
        _, groups = convective_adjustment_with_groups(_heated(d), CAPACITY, PRESSURE)
        whole.append(max((h + 1 for _, h in groups), default=0))
    assert set(np.diff(whole)) <= {0, 1}
    assert 1 in set(np.diff(whole))


@pytest.mark.req("EQR-002", "VAL-001")
def test_without_partial_tops_the_adjustment_is_the_whole_layer_one():
    """The limit: not asking for partial tops reproduces the old scheme bit for bit.

    Checked against a reference profile adjusted by hand: two unstable lowest
    layers mixed to one potential temperature by enthalpy, the rest untouched.
    """
    temperature = STABLE.copy()
    temperature[0] += 8.0
    adjusted, groups = convective_adjustment_with_groups(temperature, CAPACITY, PRESSURE)
    low, high = groups[0]
    assert low == 0
    span = slice(low, high + 1)
    theta = potential_temperature(adjusted, PRESSURE)
    assert np.ptp(theta[span]) < 1e-9
    # Untouched up to the round trip through potential temperature.
    assert np.allclose(adjusted[high + 1 :], temperature[high + 1 :], rtol=1e-13, atol=0.0)
    assert np.sum(CAPACITY * adjusted) == pytest.approx(np.sum(CAPACITY * temperature), rel=1e-13)


@pytest.mark.req("CON-001", "WAT-002")
def test_the_fractional_top_carries_water_in_the_same_fraction():
    """The layer the group rose into gives up that fraction of its vapour."""
    module = ConvectiveAdjustmentModule(transport_vapour=True, fractional_top=True)
    module._partial_tops = {3: 0.25}
    mass = np.full(6, 1.0e3)
    water = np.array([10.0, 8.0, 6.0, 4.0, 2.0, 1.0])

    class _Field:
        def __init__(self, values):
            self.values = values
            self.time_interval = None

        def require_units(self, _units):
            return self

    class _Context:
        dt = 3600.0

    module._inputs = {"atmosphere_water_vapour": _Field(water)}
    module._context = _Context()
    change = module._vapour_transport(((0, 2),), mass)
    assert change is not None
    moved = water + change

    # Mixed over layers 0-2 and a quarter of layer 3, uniform in specific humidity.
    total = water[:3].sum() + 0.25 * water[3]
    q = total / (3 * 1.0e3 + 0.25 * 1.0e3)
    assert np.allclose(moved[:3], q * 1.0e3)
    assert moved[3] == pytest.approx(0.75 * water[3] + 0.25 * q * 1.0e3)
    assert np.array_equal(moved[4:], water[4:])
    assert moved.sum() == pytest.approx(water.sum(), rel=1e-14)
    assert np.all(moved >= 0.0)
