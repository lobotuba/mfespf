"""Several columns in one run: the framework holds regions without cross-talk.

A second column is what the circulation physics needs, and the framework could
not hold one: every field and state name was a single global constant.
``regionalize`` files a kernel under a region and the coupler routes its fields
and tendencies there, with the modules unchanged.

The tests are the ones that can catch a region leaking into another. Two
unconnected copies of a column must each reproduce it run alone, to the bit;
two *different* columns must each reproduce their own; and ownership must hold
per region, so one column's module cannot write the other's state.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import OwnershipError
from mfespf.registries.states import TendencyContribution
from mfespf.run.column_configurations import column_grid, reference_column
from mfespf.run.regions import combine_regions, regionalize

pytestmark = pytest.mark.stage(3)

BASE = dict(grid=column_grid(16), areal_heat_capacity=1.0e6, initial_water_store=1.0e7)
STEPS = 200
FIELDS = ("surface_temperature", "outgoing_longwave_radiation", "atmosphere_precipitable_water")


def _value(result, name: str) -> float:
    return float(np.sum(result.steps[-1].outputs[name].values))


@pytest.mark.req("STA-001", "VER-006", "MOD-010")
def test_two_unconnected_copies_each_reproduce_the_column_alone():
    alone = reference_column(**BASE).run(STEPS)
    kernel = combine_regions(
        [regionalize(reference_column(**BASE), "a"), regionalize(reference_column(**BASE), "b")],
        run_id="REGIONS-NULL",
        purpose="Two unconnected copies of one column: the regions must not interact.",
    )
    both = kernel.run(STEPS)

    assert both.may_support_a_claim, both.manifest.invalidation_reasons
    for name in FIELDS:
        assert _value(both, f"a:{name}") == _value(alone, name)
        assert _value(both, f"b:{name}") == _value(alone, name)
    # Each region negotiated its own interfaces, with its own providers.
    assert kernel.coupler is not None
    pairs = {(n.consumer_module_id, n.provider_module_id) for n in kernel.coupler.negotiations}
    assert all(c.split("/")[0] == p.split("/")[0] for c, p in pairs)
    assert len(kernel.coupler.negotiations) % 2 == 0


@pytest.mark.req("STA-001", "MOD-010")
def test_two_different_columns_each_reproduce_their_own():
    """A leak would show here: region b is warmer-coupled and must not pull on a."""
    warm = dict(BASE, sensible_heat_coefficient=14.0)
    alone_a = reference_column(**BASE).run(STEPS)
    alone_b = reference_column(**warm).run(STEPS)
    both = combine_regions(
        [regionalize(reference_column(**BASE), "a"), regionalize(reference_column(**warm), "b")],
        run_id="REGIONS-DIFFERENT",
        purpose="Two different unconnected columns: each must match its own run alone.",
    ).run(STEPS)
    for name in FIELDS:
        assert _value(both, f"a:{name}") == _value(alone_a, name)
        assert _value(both, f"b:{name}") == _value(alone_b, name)
    assert _value(alone_a, "surface_temperature") != _value(alone_b, "surface_temperature")


@pytest.mark.req("MOD-010", "STA-004")
def test_a_module_in_one_region_cannot_write_the_other_regions_state():
    kernel = combine_regions(
        [regionalize(reference_column(**BASE), "a"), regionalize(reference_column(**BASE), "b")],
        run_id="REGIONS-OWNERSHIP",
        purpose="Ownership holds per region.",
    )
    trespass = TendencyContribution(
        state_id="b:atmosphere_layer_temperature",
        provider_module_id="a/convection",
        label="cross_region_write",
        value=np.zeros(30),
        units="K s^-1",
    )
    with pytest.raises(OwnershipError, match="not an authorized tendency provider"):
        kernel.states.assert_may_contribute(trespass)
    permitted = TendencyContribution(
        state_id="b:atmosphere_layer_temperature",
        provider_module_id="b/convection",
        label="own_region",
        value=np.zeros(30),
        units="K s^-1",
    )
    kernel.states.assert_may_contribute(permitted)
