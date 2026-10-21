"""The descending boundary layer's diagnostics (descending-boundary-layer-ventilation.md)."""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.run.column_configurations import column_grid, reference_column
from mfespf.run.two_column import two_column

pytestmark = pytest.mark.stage(3)


@pytest.mark.req("CON-009")
def test_the_surface_mixed_layer_top_is_a_layer_of_the_surface_group():
    """A surface-based column publishes the top of the group its surface mixed
    with: a layer index inside the column, and the lowest layer whenever the
    surface took no part. Diagnostic only."""
    kernel = reference_column(
        grid=column_grid(16),
        n_layers=8,
        master_dt=3600.0,
        areal_heat_capacity=1.0e6,
        initial_surface_temperature=300.0,
        initial_air_temperature=270.0,
    )
    result = kernel.run(48)
    tops = [float(s.outputs["surface_mixed_layer_top"].values) for s in result.steps]
    assert all(t == int(t) and 0 <= t < 8 for t in tops)
    assert max(tops) > 0, "a warm surface under cold air must mix upward at least once"


@pytest.mark.req("CON-009")
def test_the_published_mixing_groups_replay_the_convective_vapour_transport():
    """With the groups and partial tops published (shallow-return-flow.md),
    mixing the provisional vapour offline -- the start-of-step vapour plus every
    transport's rate this step -- reproduces convection's vapour transport
    exactly. Output only: the shared numbers are identical with it on."""
    settings = dict(
        grid=column_grid(16),
        areal_heat_capacity=1.0e6,
        tropical_latitude=20.0,
        ocean_heat_contrast=40.0,
        condense_after_convection=True,
        end_of_step_saturation=True,
    )
    plain = two_column(**settings).run(30)
    run = two_column(diagnose_mixing_groups=True, **settings).run(30)
    for a, b in zip(plain.steps, run.steps, strict=True):
        for name, field in a.outputs.items():
            assert np.array_equal(
                np.asarray(field.values), np.asarray(b.outputs[name].values), equal_nan=True
            ), name

    mixed_steps = 0
    for previous, step in zip(run.steps[:-1], run.steps[1:], strict=True):
        out, dt = step.outputs, 3600.0
        for region in ("ascent", "descent"):
            w = np.asarray(
                previous.outputs[f"{region}:atmosphere_water_vapour"].values, dtype=float
            )
            m = np.asarray(out[f"{region}:atmosphere_layer_mass"].values, dtype=float)
            dt = float(out[f"{region}:convective_mixing_group"].time_interval)
            advection = np.asarray(out[f"{region}:vapour_advection_rate"].values, dtype=float) + (
                np.asarray(out[f"{region}:circulation_vapour_advection_rate"].values, dtype=float)
            )
            provisional = w + advection * dt
            group = np.asarray(out[f"{region}:convective_mixing_group"].values)
            partial = np.asarray(out[f"{region}:convective_partial_top_fraction"].values)
            mixed = np.array(provisional)
            for g in np.unique(group[group >= 0]):
                span = np.flatnonzero(group == g)
                low, high = int(span[0]), int(span[-1])
                top = high + 1
                f = float(partial[top]) if top < w.size else 0.0
                total_water = float(np.sum(provisional[low : high + 1])) + (
                    f * float(provisional[top]) if f else 0.0
                )
                total_mass = float(np.sum(m[low : high + 1])) + (f * float(m[top]) if f else 0.0)
                mixed[low : high + 1] = total_water * m[low : high + 1] / total_mass
                if f:
                    mixed[top] = (1.0 - f) * provisional[top] + f * total_water * m[
                        top
                    ] / total_mass
                mixed_steps += 1
            label = step.state_budgets[f"{region}:atmosphere_water_vapour"]
            transport = sum(
                np.asarray(c.as_array(), dtype=float)
                for c in label.contributions
                if c.label == "convective_vapour_transport"
            )
            assert np.allclose(
                (mixed - provisional) / dt, transport, rtol=0.0, atol=1e-15 * np.abs(w).max() / dt
            ), (region, step.step_index)
    assert mixed_steps > 0, "the columns must convect at least once"
