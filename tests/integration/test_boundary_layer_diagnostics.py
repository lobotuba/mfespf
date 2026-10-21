"""The descending boundary layer's diagnostics (descending-boundary-layer-ventilation.md)."""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError
from mfespf.physics.convection import ConvectiveAdjustmentModule
from mfespf.physics.thermodynamics import (
    KAPPA,
    exner,
    lifting_condensation_pressure,
    saturation_specific_humidity,
)
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


@pytest.mark.req("CON-009", "ATM-002")
def test_subcloud_vapour_mixing_stops_at_the_surface_parcels_condensation_level():
    """With vapour_mixing_top="lifting_condensation_level" (rising-mixing-variants.md)
    the surface group's vapour is mixed only up to the layer holding the
    surface parcel's lifting condensation level, while the heat is still mixed
    over the whole group. The published vapour groups replay the transport, the
    transport moves no water into or out of the column, and the cut happens."""
    settings = dict(
        grid=column_grid(16),
        areal_heat_capacity=1.0e6,
        tropical_latitude=20.0,
        ocean_heat_contrast=40.0,
        condense_after_convection=True,
        end_of_step_saturation=True,
        diagnose_mixing_groups=True,
    )
    steps = 60
    run = two_column(vapour_mixing_top="lifting_condensation_level", **settings).run(steps)
    assert len(run.steps) == steps, "the run must complete"
    cuts = 0
    for previous, step in zip(run.steps[:-1], run.steps[1:], strict=True):
        out = step.outputs
        for region in ("ascent", "descent"):
            w = np.asarray(
                previous.outputs[f"{region}:atmosphere_water_vapour"].values, dtype=float
            )
            m = np.asarray(out[f"{region}:atmosphere_layer_mass"].values, dtype=float)
            p = np.asarray(out[f"{region}:atmosphere_layer_pressure"].values, dtype=float)
            dt = float(out[f"{region}:convective_mixing_group"].time_interval)
            provisional = w + dt * (
                np.asarray(out[f"{region}:vapour_advection_rate"].values, dtype=float)
                + np.asarray(out[f"{region}:circulation_vapour_advection_rate"].values, dtype=float)
            )
            group = np.asarray(out[f"{region}:convective_mixing_group"].values)
            partial = np.asarray(out[f"{region}:convective_partial_top_fraction"].values)
            adjusted = np.asarray(
                out[f"{region}:atmosphere_potential_temperature"].values, dtype=float
            ) * exner(p)
            heat_top = int(np.asarray(out[f"{region}:surface_mixed_layer_top"].values))
            vapour_top = int(np.flatnonzero(group == group[0])[-1]) if group[0] >= 0 else 0
            # The subcloud layer as a fixed point (erratum 1): the mean parcel of
            # the layers kept must be unsaturated at their top, or the top moves
            # down to the layer holding its condensation level.
            T0, p0, top = float(adjusted[0]), float(p[0]), heat_top
            while top >= 1:
                q = float(np.sum(provisional[: top + 1]) / np.sum(m[: top + 1]))
                T_top = T0 * (p[top] / p0) ** KAPPA
                qs = float(saturation_specific_humidity(np.array([T_top]), np.array([p[top]]))[0])
                if q < qs:
                    break
                lower = int(np.flatnonzero(p >= lifting_condensation_pressure(T0, q, p0))[-1])
                if lower >= top:
                    break
                top = lower
            assert vapour_top == max(top, 0), (region, step.step_index, vapour_top, top)
            cuts += vapour_top < heat_top
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
            label = step.state_budgets[f"{region}:atmosphere_water_vapour"]
            transport = sum(
                np.asarray(c.as_array(), dtype=float)
                for c in label.contributions
                if c.label == "convective_vapour_transport"
            )
            scale = np.abs(w).sum() / dt
            assert np.allclose((mixed - provisional) / dt, transport, rtol=0.0, atol=1e-15 * scale)
            assert abs(float(np.sum(transport))) <= 1e-12 * scale
    assert cuts > 0, "the condensation level must cut the surface group at least once"


@pytest.mark.req("CON-009")
def test_subcloud_vapour_mixing_refuses_a_column_without_surface_based_convection():
    """The surface parcel's level only bounds a group the surface takes part in."""
    with pytest.raises(ConfigurationError):
        ConvectiveAdjustmentModule(
            transport_vapour=True,
            surface_based=False,
            vapour_mixing_top="lifting_condensation_level",
        )
    with pytest.raises(ConfigurationError):
        ConvectiveAdjustmentModule(vapour_mixing_top="cloud_top")
