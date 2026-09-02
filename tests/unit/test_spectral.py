"""Stage 0 and 2: spectral grids, conservative regridding, and negotiation."""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.contracts.spectral import (
    IntegrationConvention,
    SpectralCoordinate,
    SpectralGrid,
    negotiate_spectral_grid,
    regrid_spectral_density,
)
from mfespf.core.errors import ConfigurationError, DimensionalityError

pytestmark = pytest.mark.stage(0)


def _grid(edges, grid_id="g", convention=IntegrationConvention.SPECTRAL_DENSITY) -> SpectralGrid:
    return SpectralGrid(
        grid_id=grid_id,
        coordinate=SpectralCoordinate.WAVELENGTH,
        coordinate_units="m",
        edges=np.asarray(edges, dtype=float),
        integration_convention=convention,
        resolution_identifier=grid_id,
    )


@pytest.mark.req("SPC-004")
def test_grid_declares_edges_widths_and_ordering():
    grid = _grid([1.0, 2.0, 4.0])
    assert grid.n_bins == 2
    assert grid.widths == pytest.approx([1.0, 2.0])
    assert grid.support == (1.0, 4.0)


@pytest.mark.req("SPC-004")
def test_unordered_edges_are_rejected_rather_than_sorted():
    with pytest.raises(ConfigurationError, match="strictly increasing"):
        _grid([2.0, 1.0, 3.0])


@pytest.mark.req("SPC-004")
def test_integration_convention_decides_whether_width_is_applied():
    density = _grid([0.0, 2.0], convention=IntegrationConvention.SPECTRAL_DENSITY)
    integrated = _grid([0.0, 2.0], convention=IntegrationConvention.BIN_INTEGRATED)
    assert density.integrate(np.array([3.0])) == pytest.approx(6.0)
    assert integrated.integrate(np.array([3.0])) == pytest.approx(3.0)


@pytest.mark.req("SPC-005", "SPC-006", "KER-009")
def test_regridding_conserves_the_integrated_quantity():
    source = _grid([1.0, 2.0, 3.0, 4.0], "source")
    target = _grid([1.0, 2.5, 4.0], "target")
    values = np.array([10.0, 20.0, 30.0])

    result = regrid_spectral_density(values, source, target)

    assert result.source_integral == pytest.approx(60.0)
    assert result.target_integral == pytest.approx(60.0)
    assert result.relative_residual < 1e-15
    assert "residual" in result.describe()


@pytest.mark.req("SPC-005")
def test_regridding_to_a_single_bin_preserves_the_total():
    source = _grid(np.linspace(1.0, 5.0, 9), "fine")
    target = _grid([1.0, 5.0], "broadband")
    values = np.linspace(1.0, 8.0, 8)

    result = regrid_spectral_density(values, source, target)
    assert result.target_integral == pytest.approx(source.integrate(values))


@pytest.mark.req("SPC-005", "SPC-006")
def test_source_content_outside_the_target_support_is_reported_not_discarded():
    source = _grid([1.0, 2.0, 3.0], "source")
    target = _grid([1.0, 2.0], "narrow")
    values = np.array([10.0, 40.0])

    result = regrid_spectral_density(values, source, target)
    assert result.target_integral == pytest.approx(10.0)
    assert result.unmatched_source_integral == pytest.approx(40.0)
    assert "unmatched" in result.describe()


@pytest.mark.req("SPC-003")
def test_regridding_across_incompatible_coordinates_is_refused():
    wavelength = _grid([1.0, 2.0], "lam")
    wavenumber = SpectralGrid(
        grid_id="nu",
        coordinate=SpectralCoordinate.WAVENUMBER,
        coordinate_units="cm^-1",
        edges=np.array([1.0, 2.0]),
    )
    with pytest.raises(DimensionalityError, match="cannot regrid"):
        regrid_spectral_density(np.array([1.0]), wavelength, wavenumber)


@pytest.mark.req("SPC-007", "SPC-008")
def test_negotiation_records_an_unmet_resolution_request():
    broadband = _grid([1.0e-7, 1.0e-4], "broadband")
    selected, unmet = negotiate_spectral_grid(provider_grids=[broadband], requested_min_bins=1000)
    assert selected is broadband
    assert any("1000 bins" in u for u in unmet)


@pytest.mark.req("SPC-007", "CMP-004")
def test_negotiation_prefers_the_cheapest_adequate_grid():
    coarse = _grid(np.linspace(1.0e-7, 1.0e-4, 5), "coarse")
    fine = _grid(np.linspace(1.0e-7, 1.0e-4, 101), "fine")
    selected, unmet = negotiate_spectral_grid(provider_grids=[fine, coarse], requested_min_bins=4)
    assert selected is coarse
    assert unmet == ()
