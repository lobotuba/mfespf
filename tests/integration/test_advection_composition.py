"""Two producers of one field, and the transports a column composes with
(advection-composition.md).

The tropical two-column configuration had the circulation and each column's
tropical exchange publishing the same advection fields, and the coupler kept
whichever came later: the water and convection modules composed with the
exchange's transport and never with the circulation's, in the step it was
made. Now the coupler refuses two producers of one field, each transport
publishes under its own name, and the consumers compose with every source
they declare.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError
from mfespf.run.column_configurations import column_grid
from mfespf.run.tropical_column import tropical_column
from mfespf.run.two_column import two_column

pytestmark = pytest.mark.stage(3)


def _modules(kernel):
    return {s.module_id: s.module for s in kernel.schedule.ordered()}


def _tropical_pair():
    # An ocean heat contrast makes the two columns differ, so the circulation
    # has something to carry; identical columns would give it nothing.
    return two_column(
        grid=column_grid(16),
        areal_heat_capacity=1.0e6,
        tropical_latitude=20.0,
        ocean_heat_contrast=40.0,
    )


@pytest.mark.req("STA-001", "CON-004")
def test_two_producers_of_one_field_are_refused():
    """The collision as it was: the circulation publishing the exchange's
    names into the columns. Refused, naming the field and both producers,
    rather than keeping whichever published last."""
    kernel = _tropical_pair()
    circulation = _modules(kernel)["circulation"]
    circulation.contract = replace(
        circulation.contract,
        produced_outputs=(*circulation.contract.produced_outputs, "ascent:vapour_advection_rate"),
    )
    with pytest.raises(ConfigurationError, match="ascent:vapour_advection_rate") as refused:
        kernel.prepare(n_steps=2)
    assert "circulation" in str(refused.value)
    assert "ascent/tropical_exchange" in str(refused.value)


@pytest.mark.req("CON-004", "MOD-011")
def test_the_tropical_columns_compose_with_both_transports():
    """Each column's water and convection modules receive the exchange's rate
    and the circulation's, and compose with their sum."""
    kernel = _tropical_pair()
    coupler, _ = kernel.prepare(n_steps=30)
    coupler.run(24)
    modules = _modules(kernel)
    circulation = modules["circulation"]
    for region in ("ascent", "descent"):
        exchange = modules[f"{region}/tropical_exchange"]
        for kind, own_name in (
            ("vapour", "circulation_vapour_advection_rate"),
            ("heating", "circulation_advective_heating_rate"),
        ):
            own = circulation._rates[f"{region}:{own_name}"]
            theirs = exchange._rates[
                "vapour_advection_rate" if kind == "vapour" else "advective_heating_rate"
            ]
            assert np.abs(own).max() > 0.0
            for consumer in ("water", "convection"):
                composed = modules[f"{region}/{consumer}"].advection_rate(kind)
                assert composed is not None
                assert np.array_equal(composed, theirs + own), (region, consumer, kind)


@pytest.mark.req("CON-004")
@pytest.mark.parametrize(
    "build",
    [
        lambda: two_column(grid=column_grid(16), areal_heat_capacity=1.0e6),
        lambda: tropical_column(grid=column_grid(16), areal_heat_capacity=1.0e6),
    ],
    ids=["two_column", "tropical_column"],
)
def test_a_column_with_one_transport_composes_with_exactly_that_rate(build):
    """One source: the composed rate is that source's array, not a sum that
    could round differently, so the configurations that had a single producer
    cannot move."""
    kernel = build()
    coupler, _ = kernel.prepare(n_steps=4)
    coupler.run(3)
    modules = _modules(kernel)
    for module_id, module in modules.items():
        if not module_id.endswith(("water", "convection")):
            continue
        for kind in ("vapour", "heating"):
            sources = module.advection_fields(kind)
            assert len(sources) == 1, (module_id, sources)
            composed = module.advection_rate(kind)
            if composed is None:
                continue
            assert composed is module._inputs[sources[0]].values or np.array_equal(
                composed, np.asarray(module._inputs[sources[0]].values, dtype=float)
            )
