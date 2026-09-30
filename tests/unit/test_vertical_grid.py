"""The vertical grid with a log-pressure cap aloft.

Checked against its own definition rather than against a run: without the cap
it is today's grid exactly; with it no layer is deeper than the cap (the merged
top layer by less than half a step more), everything below the first capped
interface is untouched, and the column's mass is still exactly
``(p_s - p_top) / g``.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError
from mfespf.physics.atmosphere import (
    AtmosphericColumnModule,
    hydrostatic_layer_masses,
    layer_interface_pressures,
)

pytestmark = pytest.mark.stage(2)

PS, PT, G = 101325.0, 1.0e3, 9.80665


def grid(alpha: float | None) -> np.ndarray:
    return layer_interface_pressures(
        n_layers=30, surface_pressure=PS, top_pressure=PT, max_log_step=alpha
    )


@pytest.mark.req("ATM-003")
def test_without_a_cap_the_grid_is_todays():
    assert np.array_equal(grid(None), np.linspace(PS, PT, 31))


@pytest.mark.req("ATM-003")
def test_no_layer_is_deeper_than_the_cap():
    for alpha in (0.4, 0.2, 0.1):
        interfaces = grid(alpha)
        depth = np.log(interfaces[:-1] / interfaces[1:])
        assert interfaces[0] == PS and interfaces[-1] == PT
        assert np.all(np.diff(interfaces) < 0.0)
        # Every layer but the top respects the cap exactly; the top one may
        # carry a merged remainder of under half a step.
        assert np.all(depth[:-1] <= alpha + 1e-12)
        assert depth[-1] <= 1.5 * alpha + 1e-12


@pytest.mark.req("ATM-003")
def test_below_the_cap_nothing_moves():
    today = grid(None)
    uniform = (PS - PT) / 30
    for alpha in (0.4, 0.2, 0.1):
        capped = grid(alpha)
        # Wherever today's step is thinner than the cap, the interfaces agree.
        n = int(np.argmax(today[:-1] * (1.0 - np.exp(-alpha)) < uniform))
        # Built by repeated steps rather than linspace, so equal to rounding.
        np.testing.assert_allclose(capped[: n + 1], today[: n + 1], rtol=1e-12, atol=0)


@pytest.mark.req("ATM-003", "CON-001")
def test_the_capped_column_has_exactly_its_mass():
    for alpha in (None, 0.4, 0.2, 0.1):
        interfaces, masses = hydrostatic_layer_masses(
            n_layers=30, surface_pressure=PS, top_pressure=PT, gravity=G, max_log_step=alpha
        )
        assert masses.size == interfaces.size - 1
        assert float(masses.sum()) == pytest.approx((PS - PT) / G, rel=1e-14)
    atmosphere = AtmosphericColumnModule(n_layers=30, top_pressure=PT, max_log_step=0.2)
    assert atmosphere.n_layers == 39 and atmosphere.reference_layers == 30
    assert atmosphere.layer_mass.size == 39


@pytest.mark.req("ATM-003")
def test_a_cap_that_cannot_be_built_is_refused():
    with pytest.raises(ConfigurationError, match="positive"):
        grid(0.0)
    with pytest.raises(ConfigurationError, match="above zero"):
        layer_interface_pressures(
            n_layers=30, surface_pressure=PS, top_pressure=0.0, max_log_step=0.2
        )


@pytest.mark.req("ATM-003")
def test_the_standard_atmosphere_reproduces_its_own_layer_bases():
    """US Standard Atmosphere 1976: each tabulated base is hit exactly."""
    from mfespf.physics.atmosphere import standard_atmosphere_temperature

    bases = {
        101325.0: 288.15,
        22632.06: 216.65,
        5474.889: 216.65,
        868.0187: 228.65,
        110.9063: 270.65,
    }
    computed = standard_atmosphere_temperature(np.array(list(bases)))
    np.testing.assert_allclose(computed, list(bases.values()), atol=1e-6)
    # Isothermal from 11 to 20 km, warming above it: the standard's stratosphere.
    assert standard_atmosphere_temperature(np.array([10000.0]))[0] == pytest.approx(216.65)
    profile = standard_atmosphere_temperature(np.geomspace(5000.0, 120.0, 20))
    assert np.all(np.diff(profile) > 0.0)
