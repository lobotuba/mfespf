"""The observed cloud field: ISCCP's clouds, at ISCCP's heights.

Checked against the source and against the clear column, not against the
CERES quantities the field is later judged by: the types must partition the
sky exactly as ISCCP reports them, each cloud must sit where its top pressure
says, the infrared optical thickness must follow ISCCP's phase relations, and a
field with no cloud in it must be the clear column to round-off.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError
from mfespf.physics.cloud_field import (
    ISCCP_CLOUD_TYPES,
    CloudType,
    ObservedCloudFieldModule,
    cloud_type_layers,
)
from mfespf.run.column_configurations import column_grid, reference_column

pytestmark = pytest.mark.stage(4)


@pytest.mark.req("CLD-002", "EQR-002")
def test_the_cloud_types_partition_the_sky():
    """ISCCP's D2 amounts sum to its observed total cloud cover, 0.675 to 0.687."""
    total = sum(t.amount for t in ISCCP_CLOUD_TYPES)
    assert total == pytest.approx(0.687, abs=1e-9)
    assert 0.65 < total < 0.70
    too_much = ObservedCloudFieldModule(
        cloud_types=(CloudType("a", 0.7, 5.0e4, 3.0, 0.0), CloudType("b", 0.4, 8.0e4, 3.0, 0.0))
    )
    with pytest.raises(ConfigurationError, match="partition the sky"):
        too_much.validate_configuration()


@pytest.mark.req("CLD-002")
def test_each_cloud_sits_in_the_layer_containing_its_top():
    interfaces = np.array([101325.0, 90000.0, 70000.0, 50000.0, 30000.0, 20000.0, 1000.0])
    tops = np.array([26700.0, 32600.0, 55700.0, 82600.0])
    # 267 hPa lies in 300-200, 326 in 500-300, 557 in 700-500, 826 in 900-700.
    assert cloud_type_layers(tops, interfaces).tolist() == [4, 3, 2, 1]
    with pytest.raises(ConfigurationError, match="outside the column"):
        cloud_type_layers(np.array([500.0]), interfaces)


@pytest.mark.req("CLD-002", "EQR-002")
def test_the_infrared_optical_thickness_follows_the_phase():
    """Rossow and Schiffer's section 3: tau_ir = tau / 2.13 for ice, / 2.56 for liquid."""
    ice = CloudType("ice", 0.1, 3.0e4, 4.26, 1.0)
    liquid = CloudType("liquid", 0.1, 8.0e4, 5.12, 0.0)
    mixed = CloudType("mixed", 0.1, 5.5e4, 1.0, 0.5)
    assert ice.infrared_optical_thickness == pytest.approx(2.0)
    assert liquid.infrared_optical_thickness == pytest.approx(2.0)
    assert mixed.infrared_optical_thickness == pytest.approx(0.5 / 2.13 + 0.5 / 2.56)


@pytest.mark.req("CLD-001", "VAL-001")
def test_an_empty_field_is_the_clear_column():
    """Zero amounts: the mixed sky is the clear sky, and the effects are zero.

    Run on the same surface reflectance as the cloud field uses, so the only
    difference between the two columns is the machinery of the sub-columns.
    """
    empty = tuple(
        CloudType(t.name, 0.0, t.top_pressure_pa, t.visible_optical_thickness, t.ice_fraction)
        for t in ISCCP_CLOUD_TYPES
    )
    base = dict(grid=column_grid(16), areal_heat_capacity=1.0e6, master_dt=3600.0)
    with_field = reference_column(observed_cloud_field=True, **base)
    for slot in with_field.schedule.ordered():
        if isinstance(slot.module, ObservedCloudFieldModule):
            slot.module.cloud_types = empty
    cloudy = with_field.run(50)
    outputs = cloudy.steps[-1].outputs
    assert float(outputs["cloud_shortwave_radiative_effect"].values) == 0.0
    assert float(outputs["cloud_longwave_radiative_effect"].values) == 0.0
    assert float(outputs["reflected_shortwave_radiation"].values) == pytest.approx(
        float(outputs["clear_sky_reflected_shortwave_radiation"].values), rel=1e-14
    )
