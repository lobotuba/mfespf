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
from mfespf.evaluation.bias import BiasReport, surface_temperature_bias
from mfespf.physics.equations import build_equation_registry
from mfespf.registries.states import StateRegistry
from mfespf.run.configurations import build_state_registry

#: Collected across the whole session so the summary can print them together.
_BIAS_REPORTS: list[BiasReport] = []


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


# --- the error in kelvin, reported from every test that measures one --------


@pytest.fixture(scope="session")
def _bias_log() -> list[BiasReport]:
    return _BIAS_REPORTS


@pytest.fixture
def report_bias(request):
    """Reduce a completed run to its error in kelvin and report it at the end.

    Any test that equilibrates a column should call this, whatever else it is
    checking. The error against the observed surface temperature is the number
    the project is judged on, and it was previously recomputed by hand in scratch
    scripts each time anyone wanted it -- so it moved without anyone noticing
    between one commit and the next.

    The returned report refuses to state an error for a run that has not settled.
    That is not caution: the dry column reads -16.3 K at 20000 steps while
    carrying +43.8 W m^-2, and it is unfinished rather than cold.
    """

    def _report(result, *, label: str | None = None, tail_fraction: float = 0.2) -> BiasReport:
        report = surface_temperature_bias(
            result, label=label or request.node.name, tail_fraction=tail_fraction
        )
        _BIAS_REPORTS.append(report)
        return report

    return _report


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    """Print every error in kelvin the run measured, together, at the end.

    Together is the point. One configuration's bias is a number; the set of them
    is the shape of the problem -- which configurations are equilibria, which are
    unfinished, and which look good only because two large errors cancel.
    """
    if not _BIAS_REPORTS:
        return
    terminalreporter.write_sep("=", "surface temperature error (K)")
    settled = [r for r in _BIAS_REPORTS if r.equilibrated]
    for report in _BIAS_REPORTS:
        terminalreporter.write_line(report.describe())
    if settled:
        worst = max(settled, key=lambda r: abs(r.error_k or 0.0))
        terminalreporter.write_line("")
        terminalreporter.write_line(
            f"{len(settled)} of {len(_BIAS_REPORTS)} runs equilibrated; "
            f"largest error {worst.error_k:+.3f} K ({worst.label})"
        )
