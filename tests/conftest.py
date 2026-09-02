"""Shared fixtures.

Tests bind themselves to requirements with ``@pytest.mark.req("ID", ...)``.
``tools/traceability.py`` reads those markers to build the acceptance matrix, so
a marker citing a requirement that does not exist is a build error rather than a
harmless typo.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.contracts.spectral import (
    IntegrationConvention,
    SpectralCoordinate,
    SpectralGrid,
)
from mfespf.core.quantities import (
    Field,
    Provenance,
    QuantityKind,
    SpatialSupport,
    TimeSemantics,
)
from mfespf.physics.equations import build_equation_registry
from mfespf.registries.states import StateRegistry
from mfespf.run.configurations import build_state_registry


@pytest.fixture
def equations():
    return build_equation_registry()


@pytest.fixture
def states() -> StateRegistry:
    return build_state_registry()


@pytest.fixture
def linear_grid() -> SpectralGrid:
    """A simple uniform grid, easy to reason about in regridding tests."""
    return SpectralGrid(
        grid_id="uniform_4",
        coordinate=SpectralCoordinate.WAVELENGTH,
        coordinate_units="m",
        edges=np.array([1.0e-6, 2.0e-6, 3.0e-6, 4.0e-6, 5.0e-6]),
        integration_convention=IntegrationConvention.SPECTRAL_DENSITY,
        resolution_identifier="uniform-4",
    )


@pytest.fixture
def scalar_field():
    def _make(field_id: str, value: float, units: str, **kwargs) -> Field:
        defaults = {
            "physical_name": field_id.replace("_", " "),
            "quantity_kind": QuantityKind.DIAGNOSTIC,
            "time_semantics": TimeSemantics.INSTANTANEOUS,
            "spatial_support": SpatialSupport.GLOBAL_SCALAR,
            "provenance": Provenance(producer="test"),
        }
        defaults.update(kwargs)
        return Field(field_id=field_id, values=np.array(value), units=units, **defaults)

    return _make
