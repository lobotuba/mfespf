"""Per-region settings on the two-column kernel (descending-boundary-layer-ventilation.md, V1)."""

from __future__ import annotations

import pytest

from mfespf.run.column_configurations import column_grid
from mfespf.run.two_column import two_column

pytestmark = pytest.mark.stage(3)


def _water(kernel, region):
    return next(s.module for s in kernel.schedule.ordered() if s.module_id == f"{region}/water")


@pytest.mark.req("CON-009")
def test_regional_settings_change_one_region_only_and_refuse_unknown_regions():
    """An override names its region and reaches that region's modules alone."""
    common = dict(grid=column_grid(16), n_layers=8, master_dt=3600.0, areal_heat_capacity=1.0e6)
    kernel = two_column(
        regional_settings={"descent": {"eddy_diffusivity": 0.0}}, eddy_diffusivity=0.39, **common
    )
    assert _water(kernel, "descent").eddy_diffusivity == 0.0
    assert _water(kernel, "ascent").eddy_diffusivity == 0.39
    with pytest.raises(ValueError, match="unknown regions"):
        two_column(regional_settings={"sideways": {}}, **common)
