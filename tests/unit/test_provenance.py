"""Water provenance: molecule-once milestones survive any amount of recycling.

fall-path-closure-candidate.md, F6. Checked against the definition itself: a
parcel taken through repeated condensation, evaporation and precipitation
crosses each milestone exactly once, whatever its gross phase-change totals.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.physics.provenance import (
    Milestones,
    class_count,
    condense,
    evaporate,
    form_precipitation,
    proportional,
    surface_arrival,
)

pytestmark = pytest.mark.stage(3)


@pytest.mark.req("CON-001", "VAL-001")
def test_a_recycled_parcel_crosses_each_milestone_once():
    """condense, evaporate, condense, evaporate, condense, form rain, evaporate,
    condense, form rain, reach the surface: C_unique = C_precip,unique =
    P_surface,unique = m, so FE = SE = PE = 1, while the gross condensation is
    4 m and the gross precipitation formation 2 m."""
    n, layer, mass = 4, 2, 3.7
    milestones = Milestones.zeros(n)
    vapour = np.zeros(class_count(n))
    vapour[0] = mass
    gross_condensed = gross_formed = 0.0
    for _ in range(2):
        vapour, cloud = condense(vapour, mass, layer, milestones, large_scale=False)
        gross_condensed += mass
        cloud, vapour = evaporate(cloud, mass)
    vapour, cloud = condense(vapour, mass, layer, milestones, large_scale=False)
    gross_condensed += mass
    cloud, rain = form_precipitation(cloud, mass, layer, milestones)
    gross_formed += mass
    rain, vapour = evaporate(rain, mass)
    vapour, cloud = condense(vapour, mass, layer + 1, milestones, large_scale=False)
    gross_condensed += mass
    cloud, rain = form_precipitation(cloud, mass, layer + 1, milestones)
    gross_formed += mass
    surface_arrival(rain, milestones)

    assert gross_condensed == pytest.approx(4 * mass) and gross_formed == pytest.approx(2 * mass)
    first_condensed = float(np.sum(milestones.first_condensed))
    first_precip = float(np.sum(milestones.first_precip))
    first_surface = float(np.sum(milestones.first_surface))
    assert first_condensed == pytest.approx(mass, rel=1e-15)
    assert first_precip == pytest.approx(mass, rel=1e-15)
    assert first_surface == pytest.approx(mass, rel=1e-15)
    # All of it is counted at the layer it first condensed in, and formed its
    # first precipitation there.
    assert milestones.first_condensed[layer] == pytest.approx(mass, rel=1e-15)
    assert milestones.first_precip[layer] == pytest.approx(mass, rel=1e-15)
    assert milestones.first_precip_location[layer] == pytest.approx(mass, rel=1e-15)
    assert milestones.first_surface[layer] == pytest.approx(mass, rel=1e-15)
    assert not np.any(vapour) and not np.any(cloud)


@pytest.mark.req("CON-001")
def test_taking_from_a_pool_is_proportional_and_empties_it_exactly():
    pool = np.array([1.0, 0.0, 2.0, 0.0, 0.0, 1.0, 0.0])
    taken = proportional(pool, 2.0)
    np.testing.assert_allclose(taken, pool / 2.0, rtol=1e-15)
    whole = proportional(pool, 10.0)
    assert np.array_equal(pool - whole, np.zeros_like(pool))
    assert not np.any(proportional(pool, 0.0))


@pytest.mark.req("VAL-001")
def test_only_precipitation_reaches_the_surface():
    """A cohort ends at the surface, and only water that formed precipitation
    gets there: vapour or cloud arriving would be a bookkeeping error."""
    milestones = Milestones.zeros(2)
    with pytest.raises(ValueError, match="formed precipitation"):
        surface_arrival(np.array([1.0, 0.0, 0.0, 0.0, 0.0]), milestones)
