"""The descending boundary layer's diagnostics (descending-boundary-layer-ventilation.md)."""

from __future__ import annotations

import pytest

from mfespf.run.column_configurations import column_grid, reference_column

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
