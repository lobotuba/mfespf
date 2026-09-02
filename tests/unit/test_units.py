"""Stage 0: units, dimensions and the coercion prohibition."""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError, DimensionalityError, UnitError
from mfespf.core.quantities import (
    Field,
    MissingValuePolicy,
    Provenance,
    QuantityKind,
    SpatialSupport,
    TimeSemantics,
)
from mfespf.core.units import (
    ENERGY,
    IRRADIANCE,
    POWER,
    Dimension,
    convert,
    parse_unit,
    require_compatible,
)

pytestmark = pytest.mark.stage(0)


@pytest.mark.req("DAT-003", "KER-002")
def test_registered_units_carry_correct_dimensions():
    assert parse_unit("W m^-2").dimension == IRRADIANCE
    assert parse_unit("J").dimension == ENERGY
    assert parse_unit("W").dimension == POWER
    assert parse_unit("J").dimension / parse_unit("s").dimension == POWER


@pytest.mark.req("DAT-005")
def test_incompatible_dimensions_fail_rather_than_coerce():
    with pytest.raises(DimensionalityError, match="incompatible dimensions"):
        require_compatible(parse_unit("W m^-2"), parse_unit("K"), context="test")


@pytest.mark.req("DAT-005", "KER-002")
def test_unknown_unit_is_an_error_not_an_assumption():
    with pytest.raises(UnitError, match="not registered"):
        parse_unit("furlongs per fortnight")


@pytest.mark.req("DAT-004")
def test_conversion_returns_a_logged_record():
    value, record = convert(1.0, parse_unit("hour"), parse_unit("s"), context="boundary")
    assert value == pytest.approx(3600.0)
    assert record.from_unit == "hour" and record.to_unit == "s"
    assert "boundary" in record.describe()


@pytest.mark.req("DAT-004")
def test_affine_temperature_conversion_uses_the_declared_offset():
    kelvin, _ = convert(0.0, parse_unit("degC"), parse_unit("K"), context="ingest")
    assert kelvin == pytest.approx(273.15)
    assert parse_unit("degC").is_affine
    assert not parse_unit("K").is_affine


@pytest.mark.req("DAT-003")
def test_dimension_algebra_is_exact_with_rational_exponents():
    half = Dimension.of(length=1) ** 0.5
    assert half * half == Dimension.of(length=1)
    assert (Dimension.of(length=2) / Dimension.of(length=2)).is_dimensionless


@pytest.mark.req("DAT-001", "DAT-002")
def test_field_requires_complete_metadata(scalar_field):
    field = scalar_field("toa_flux", 340.25, "W m^-2", quantity_kind=QuantityKind.FLUX)
    assert field.unit.dimension == IRRADIANCE
    assert field.quantity_kind is QuantityKind.FLUX
    assert field.provenance.producer == "test"


@pytest.mark.req("DAT-006")
def test_interval_semantics_require_an_interval():
    with pytest.raises(ConfigurationError, match="DAT-006"):
        Field(
            field_id="x",
            physical_name="x",
            values=np.array(1.0),
            units="W m^-2",
            quantity_kind=QuantityKind.FLUX,
            time_semantics=TimeSemantics.INTERVAL_AVERAGE,
            spatial_support=SpatialSupport.GLOBAL_SCALAR,
            provenance=Provenance(producer="test"),
        )


@pytest.mark.req("DAT-007")
def test_missing_data_is_not_silently_zero_filled():
    with pytest.raises(ConfigurationError, match="DAT-007"):
        Field(
            field_id="x",
            physical_name="x",
            values=np.array([1.0, np.nan]),
            units="K",
            quantity_kind=QuantityKind.INTENSIVE_STATE,
            time_semantics=TimeSemantics.INSTANTANEOUS,
            spatial_support=SpatialSupport.GLOBAL_SCALAR,
            provenance=Provenance(producer="test"),
        )
    # Declaring the policy makes the same data admissible: the requirement is
    # that the behaviour be explicit, not that gaps be forbidden.
    field = Field(
        field_id="x",
        physical_name="x",
        values=np.array([1.0, np.nan]),
        units="K",
        quantity_kind=QuantityKind.INTENSIVE_STATE,
        time_semantics=TimeSemantics.INSTANTANEOUS,
        spatial_support=SpatialSupport.GLOBAL_SCALAR,
        provenance=Provenance(producer="test"),
        missing_value_policy=MissingValuePolicy.PROPAGATE_NAN,
    )
    assert np.isnan(field.values[1])


@pytest.mark.req("DAT-005")
def test_require_units_rejects_a_compatible_but_differently_scaled_unit(scalar_field):
    field = scalar_field("s", 1.0, "hPa")
    with pytest.raises(DimensionalityError, match="convert at a declared boundary"):
        field.require_units("Pa")
