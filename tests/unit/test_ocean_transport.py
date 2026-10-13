"""The ocean's heat transport between a rising and a sinking region.

ocean-heat-transport-circulation.md: the ocean moves heat and makes none; each
surface warms at its delivered flux over its heat capacity; with no contrast the
two columns stay identical, so the transport is the only thing that breaks
their symmetry.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError
from mfespf.physics.ocean_transport import (
    OBSERVED_OCEAN_HEAT_CONTRAST_W_M2,
    OceanHeatTransport,
    ocean_heat_delivered,
)
from mfespf.run.column_configurations import column_grid
from mfespf.run.two_column import two_column

pytestmark = pytest.mark.stage(3)

SIGMA = 0.558
CAPACITY = 1.0e6


def _surfaces(run):
    out = run.steps[-1].outputs
    return tuple(float(out[f"{r}:surface_temperature"].values) for r in ("ascent", "descent"))


@pytest.mark.req("CON-001", "CON-004", "EQR-002")
def test_the_ocean_moves_heat_and_makes_none():
    """sigma Q_r + (1 - sigma) Q_s = 0 for any contrast, the transfer carries exactly
    what leaves one surface, and a coupled run closes its budgets."""
    rng = np.random.default_rng(11)
    for _ in range(100):
        sigma, contrast = rng.uniform(0.05, 0.95), rng.uniform(-30.0, 30.0)
        q_r, q_s = ocean_heat_delivered(contrast, sigma)
        assert sigma * q_r + (1.0 - sigma) * q_s == pytest.approx(0.0, abs=1e-13)
        assert q_r - q_s == pytest.approx(contrast, rel=1e-12)
    with pytest.raises(ConfigurationError, match="fraction"):
        ocean_heat_delivered(10.0, 1.0)

    contrast = OBSERVED_OCEAN_HEAT_CONTRAST_W_M2["oaflux"]
    run = two_column(
        grid=column_grid(16),
        areal_heat_capacity=CAPACITY,
        ascent_fraction=SIGMA,
        ocean_heat_contrast=contrast,
    ).run(20)
    assert run.may_support_a_claim, run.manifest.invalidation_reasons
    assert run.worst_relative_residual() < 1e-10
    out = run.steps[-1].outputs
    q_r = float(out["ascent:ocean_heat_convergence"].values)
    q_s = float(out["descent:ocean_heat_convergence"].values)
    assert q_r == pytest.approx((1.0 - SIGMA) * contrast, rel=1e-12)
    assert SIGMA * q_r + (1.0 - SIGMA) * q_s == pytest.approx(0.0, abs=1e-12)


@pytest.mark.req("CON-004", "EQR-002")
def test_each_surface_warms_at_its_flux_over_its_heat_capacity():
    """One step from identical states: the surfaces differ from the run without
    transport by exactly Q dt / C, everything else being the same at the start."""
    contrast, dt = 13.39, 3600.0
    options = dict(grid=column_grid(16), areal_heat_capacity=CAPACITY, ascent_fraction=SIGMA)
    with_ocean = _surfaces(two_column(ocean_heat_contrast=contrast, **options).run(1))
    without = _surfaces(two_column(**options).run(1))
    q_r, q_s = ocean_heat_delivered(contrast, SIGMA)
    assert with_ocean[0] - without[0] == pytest.approx(q_r * dt / CAPACITY, rel=1e-9)
    assert with_ocean[1] - without[1] == pytest.approx(q_s * dt / CAPACITY, rel=1e-9)


@pytest.mark.req("EQR-002", "VAL-001")
def test_no_contrast_leaves_the_two_columns_identical():
    """With dQ = 0 the module is present and the columns stay bit-for-bit equal: no
    contrast, no circulation, as spectral-wtg-circulation.md Test A found."""
    run = two_column(
        grid=column_grid(16),
        areal_heat_capacity=CAPACITY,
        ascent_fraction=SIGMA,
        ocean_heat_contrast=0.0,
    ).run(30)
    out = run.steps[-1].outputs
    up, down = _surfaces(run)
    assert up == down
    assert np.array_equal(
        np.asarray(out["ascent:atmosphere_water_vapour"].values),
        np.asarray(out["descent:atmosphere_water_vapour"].values),
    )
    assert float(out["circulation_subsidence_omega_500"].values) == 0.0


@pytest.mark.req("SPC-001", "EQR-002")
def test_the_tropical_two_column_builds_and_closes():
    """Each region can be the tropical column of a band, with the transport between
    them; the run closes through the lateral boundaries and the crossing."""
    run = two_column(
        grid=column_grid(16),
        areal_heat_capacity=CAPACITY,
        ascent_fraction=SIGMA,
        tropical_latitude=20.0,
        ocean_heat_contrast=13.39,
        relaxation_distance=7.34e5,
    ).run(10)
    assert run.may_support_a_claim, run.manifest.invalidation_reasons
    assert run.worst_relative_residual() < 1e-10
    module = OceanHeatTransport(ascent_fraction=SIGMA, heat_contrast=13.39)
    assert module.parameters()["into_ascent_W_m2"] == pytest.approx(0.442 * 13.39, rel=1e-12)
