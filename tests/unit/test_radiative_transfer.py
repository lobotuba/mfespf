"""Stage 2: correctness gates for the column radiative transfer solver.

These test the solver in isolation, before any coupling.  The most valuable one
is the analytic greenhouse identity: the closed form for a single grey layer is
derived independently of this code, so agreement is evidence about the solver
rather than about its own arithmetic.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from mfespf.contracts.spectral import (
    IntegrationConvention,
    SpectralCoordinate,
    SpectralGrid,
)
from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError
from mfespf.physics.radiation import band_exitance, greybody_exitance
from mfespf.physics.radiative_transfer import (
    grey_absorption,
    solve_column_radiation,
    window_absorption,
)
from mfespf.run.column_configurations import column_grid

pytestmark = pytest.mark.stage(2)

SIGMA = value_of("stefan_boltzmann")
TSI = value_of("total_solar_irradiance_reference")
SOLAR_BAND_EDGE_M = 4.0e-6


def _two_band_grid() -> SpectralGrid:
    """One band for sunlight, one for terrestrial radiation, split at 4 microns."""
    return SpectralGrid(
        grid_id="solar_terrestrial_2band",
        coordinate=SpectralCoordinate.WAVELENGTH,
        coordinate_units="m",
        edges=np.array([1.0e-7, SOLAR_BAND_EDGE_M, 1.0e-3]),
        integration_convention=IntegrationConvention.SPECTRAL_DENSITY,
        resolution_identifier="solar-terrestrial-2band",
    )


def _incident_in_solar_band(grid: SpectralGrid, total: float) -> np.ndarray:
    """Place all incident energy in band 0 as a spectral irradiance."""
    values = np.zeros(grid.n_bins)
    values[0] = total / grid.widths[0]
    return values


@pytest.mark.req("RAD-002", "CON-009")
def test_radiation_flux_divergence_closes_identically():
    """What enters the column and does not leave is what heated the layers."""
    grid = column_grid(16)
    solution = solve_column_radiation(
        grid=grid,
        incident_spectral_irradiance=_incident_in_solar_band(grid, 340.25),
        layer_temperature=np.array([260.0, 240.0, 220.0, 210.0]),
        layer_mass=np.full(4, 2500.0),
        absorption_coefficient=grey_absorption(grid, 1.0e-4),
        surface_temperature=288.0,
        surface_albedo=np.full(grid.n_bins, 0.2),
    )
    assert solution.energy_residual() == pytest.approx(0.0, abs=1e-9)


@pytest.mark.req("RAD-002", "CON-006")
def test_incident_equals_absorbed_plus_reflected_plus_transmitted():
    """Charter 9.3's conservation role, checked with emission switched off.

    At zero temperature nothing emits, so the only energy in the column is what
    came in at the top; it must reappear exactly as layer absorption plus surface
    absorption plus what leaves the top.
    """
    grid = column_grid(16)
    incident_total = 340.25
    solution = solve_column_radiation(
        grid=grid,
        incident_spectral_irradiance=_incident_in_solar_band(grid, incident_total),
        layer_temperature=np.zeros(4),
        layer_mass=np.full(4, 2500.0),
        absorption_coefficient=grey_absorption(grid, 1.0e-4),
        surface_temperature=0.0,
        surface_albedo=np.full(grid.n_bins, 0.25),
    )
    absorbed_in_layers = float(np.sum(solution.layer_heating))
    escaping = solution.interface_up(solution.toa_index)
    assert absorbed_in_layers + solution.surface_absorbed + escaping == pytest.approx(
        incident_total, rel=1e-12
    )


@pytest.mark.req("RAD-002", "VAL-001", "KER-012")
def test_transparent_column_transmits_the_incident_field_unchanged():
    """The null case: zero opacity must leave the incident field untouched."""
    grid = column_grid(16)
    incident = _incident_in_solar_band(grid, 340.25)
    solution = solve_column_radiation(
        grid=grid,
        incident_spectral_irradiance=incident,
        layer_temperature=np.full(6, 250.0),
        layer_mass=np.full(6, 1700.0),
        absorption_coefficient=grey_absorption(grid, 0.0),
        surface_temperature=288.0,
        surface_albedo=np.zeros(grid.n_bins),
    )
    assert np.allclose(solution.layer_heating, 0.0, atol=1e-12)
    assert solution.interface_down(0) == pytest.approx(340.25, rel=1e-12)
    # A black surface at 288 K radiates sigma T^4 straight to space.
    assert solution.interface_up(solution.toa_index) == pytest.approx(
        greybody_exitance(288.0), rel=1e-12
    )


@pytest.mark.req("RAD-002", "VAL-001", "EQR-002")
def test_single_grey_layer_matches_the_analytic_greenhouse():
    """Evaluate the solver at the closed-form temperatures; residuals must vanish.

    For one layer transparent to sunlight and grey with emissivity ``eps`` in the
    thermal infrared, the steady state is ``T_s^4 = 2 T_e^4 / (2 - eps)`` and
    ``T_a^4 = T_s^4 / 2``. That identity is derived independently of this code, so
    checking it by *evaluation* rather than by integration isolates the solver
    from every other part of the framework -- no coupler, no time stepping, no
    convergence question.

    The residual is bounded by one real physical effect rather than by tolerance
    shopping: the surface's own emission below 4 microns is not zero, and the
    analytic model assumes it is. The test measures that leak and requires the
    residual to sit inside it.
    """
    grid = _two_band_grid()
    absorbed_solar = 238.175  # (1 - 0.3) * 1361 / 4
    T_e = (absorbed_solar / SIGMA) ** 0.25

    for eps in (0.25, 0.5, 0.75, 0.9):
        T_s = T_e * (2.0 / (2.0 - eps)) ** 0.25
        T_a = T_s / 2.0**0.25

        layer_mass = np.array([1.0e4])
        diffusivity = 1.66
        kappa_lw = -np.log1p(-eps) / (diffusivity * layer_mass[0])
        absorption = np.array([0.0, kappa_lw])

        solution = solve_column_radiation(
            grid=grid,
            incident_spectral_irradiance=_incident_in_solar_band(grid, absorbed_solar),
            layer_temperature=np.array([T_a]),
            layer_mass=layer_mass,
            absorption_coefficient=absorption,
            surface_temperature=T_s,
            surface_albedo=np.zeros(grid.n_bins),
            diffusivity_factor=diffusivity,
        )

        # The one assumption of the analytic model this configuration violates.
        solar_band_leak = float(band_exitance(grid, T_s)[0])
        bound = max(2.0 * solar_band_leak, 1e-9)

        surface_net = solution.interface_down(0) - solution.interface_up(0)
        assert surface_net == pytest.approx(0.0, abs=bound), f"eps={eps}"
        assert float(solution.layer_heating[0]) == pytest.approx(0.0, abs=bound), f"eps={eps}"
        # The gate is only meaningful if the term it allows for is small. It is
        # measured, not assumed: the surface emits 0.07 percent of its output
        # below 4 microns at eps = 0.25 and 0.34 percent at eps = 0.9, the growth
        # being the Wien tail moving into the solar band as the surface warms.
        # So the identity is verified to within half a percent, and the bound is
        # a physical quantity rather than a fitted tolerance.
        assert solar_band_leak / absorbed_solar < 5.0e-3


@pytest.mark.req("RAD-002", "VAL-001")
def test_zero_emissivity_layer_recovers_the_airless_balance():
    """eps = 0 must reduce the greenhouse identity to T_s = T_e exactly."""
    grid = _two_band_grid()
    absorbed_solar = 238.175
    T_e = (absorbed_solar / SIGMA) ** 0.25
    solution = solve_column_radiation(
        grid=grid,
        incident_spectral_irradiance=_incident_in_solar_band(grid, absorbed_solar),
        layer_temperature=np.array([200.0]),
        layer_mass=np.array([1.0e4]),
        absorption_coefficient=np.zeros(grid.n_bins),
        surface_temperature=T_e,
        surface_albedo=np.zeros(grid.n_bins),
    )
    assert solution.interface_down(0) - solution.interface_up(0) == pytest.approx(0.0, abs=1e-9)
    assert float(solution.layer_heating[0]) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.req("SPC-009")
def test_equal_opacity_distributed_differently_gives_a_different_answer():
    """SPC-009: equal integrated opacity, different spectrum, different response.

    A grey coefficient and a window coefficient are chosen to carry the same
    band-summed opacity. The outgoing longwave differs, which is exactly the
    behaviour SPC-009 forbids assuming away.
    """
    grid = column_grid(32)
    windowed = window_absorption(grid, kappa_opaque=1.0e-3, kappa_window=0.0)
    grey = grey_absorption(grid, float(np.mean(windowed)))
    assert np.sum(grey) == pytest.approx(np.sum(windowed), rel=1e-12)

    common = dict(
        grid=grid,
        incident_spectral_irradiance=_incident_in_solar_band(grid, 340.25),
        layer_temperature=np.array([270.0, 250.0, 230.0, 215.0]),
        layer_mass=np.full(4, 2500.0),
        surface_temperature=288.0,
        surface_albedo=np.zeros(grid.n_bins),
    )
    a = solve_column_radiation(absorption_coefficient=grey, **common)  # type: ignore[arg-type]
    b = solve_column_radiation(absorption_coefficient=windowed, **common)  # type: ignore[arg-type]

    _, olr_a = a.partition(a.upward[a.toa_index])
    _, olr_b = b.partition(b.upward[b.toa_index])
    assert abs(olr_a - olr_b) > 1.0  # W m^-2, far beyond numerical noise


@pytest.mark.req("SPC-007", "VAL-003", "KER-009")
def test_a_grey_atmosphere_is_exactly_band_count_invariant():
    """With a grey coefficient the answer must not depend on the band count at all.

    Every band then shares one emissivity, and the thermal emission is normalised
    to sigma T^4 whatever the grid, so the transfer is linear in band-summed
    quantities. Agreement to the last digit -- not merely convergence -- is the
    strongest available statement that the spectral machinery adds no bias of its
    own, and it is what makes the *spectrally structured* case below interpretable
    as physics rather than discretisation.
    """

    def toa_upward(n_bands: int) -> float:
        grid = column_grid(n_bands)
        solution = solve_column_radiation(
            grid=grid,
            incident_spectral_irradiance=_incident_in_solar_band(grid, 340.25),
            layer_temperature=np.array([270.0, 250.0, 230.0, 215.0]),
            layer_mass=np.full(4, 2500.0),
            absorption_coefficient=grey_absorption(grid, 1.0e-4),
            surface_temperature=288.0,
            surface_albedo=np.zeros(grid.n_bins),
        )
        return solution.interface_up(solution.toa_index)

    reference = toa_upward(8)
    for n_bands in (16, 32, 64, 128, 256):
        assert toa_upward(n_bands) == pytest.approx(reference, rel=1e-12)


@pytest.mark.req("SPC-007", "VAL-003")
def test_spectral_convergence_of_the_column_solution():
    """A spectrally structured atmosphere must converge as the bands refine.

    VAL-003 asks spectral modules for a convergence test. Here the coefficient is
    discontinuous at the window edges, so a coarse grid genuinely cannot resolve
    it and the outgoing longwave is badly wrong at 16 bands. What must be true is
    that refining fixes it: the successive differences shrink and the finest pair
    agrees to a fraction of a watt.
    """

    def olr(n_bands: int) -> float:
        grid = column_grid(n_bands)
        solution = solve_column_radiation(
            grid=grid,
            incident_spectral_irradiance=_incident_in_solar_band(grid, 340.25),
            layer_temperature=np.array([270.0, 250.0, 230.0, 215.0]),
            layer_mass=np.full(4, 2500.0),
            absorption_coefficient=window_absorption(grid, kappa_opaque=1.0e-3),
            surface_temperature=288.0,
            surface_albedo=np.zeros(grid.n_bins),
        )
        return solution.partition(solution.upward[solution.toa_index])[1]

    values = [olr(n) for n in (64, 128, 256, 512)]
    differences = [abs(b - a) for a, b in itertools.pairwise(values)]
    assert differences[-1] < differences[0]
    assert differences[-1] < 0.1  # W m^-2 between the two finest grids

    # And the coarse grid is genuinely inadequate here, which is the reason the
    # configuration's declared spectral resolution belongs in its manifest.
    assert abs(olr(16) - values[-1]) > 10.0


@pytest.mark.req("RAD-002", "CLD-001")
def test_shortwave_and_longwave_are_separately_available():
    """RAD-002 and CLD-001: the components exist before any net is formed."""
    grid = column_grid(32)
    solution = solve_column_radiation(
        grid=grid,
        incident_spectral_irradiance=_incident_in_solar_band(grid, 340.25),
        layer_temperature=np.full(4, 250.0),
        layer_mass=np.full(4, 2500.0),
        absorption_coefficient=window_absorption(grid, kappa_opaque=1.0e-3),
        surface_temperature=288.0,
        surface_albedo=np.full(grid.n_bins, 0.3),
    )
    toa = solution.upward[solution.toa_index]
    reflected, outgoing_longwave = solution.partition(toa)
    assert reflected > 0.0 and outgoing_longwave > 0.0
    assert reflected + outgoing_longwave == pytest.approx(
        solution.interface_up(solution.toa_index), rel=1e-12
    )


@pytest.mark.req("MOD-002", "SYS-006")
def test_mismatched_array_lengths_are_refused():
    grid = column_grid(8)
    with pytest.raises(ConfigurationError, match="absorption_coefficient has"):
        solve_column_radiation(
            grid=grid,
            incident_spectral_irradiance=np.zeros(grid.n_bins),
            layer_temperature=np.full(3, 250.0),
            layer_mass=np.full(3, 2500.0),
            absorption_coefficient=np.zeros(grid.n_bins + 1),
            surface_temperature=280.0,
            surface_albedo=np.zeros(grid.n_bins),
        )


@pytest.mark.req("EQR-003", "UNC-003")
def test_diffusivity_is_a_declared_parameterization():
    """The diffusivity factor must be classified as parameterized, not measured."""
    from mfespf.core.constants import ParameterClass, get

    constant = get("radiative_diffusivity_factor")
    assert constant.parameter_class is ParameterClass.PARAMETERIZED
    assert constant.standard_uncertainty and constant.standard_uncertainty > 0.0

    # And it changes the answer, so it is not a free stylistic choice.
    grid = column_grid(16)
    common = dict(
        grid=grid,
        incident_spectral_irradiance=_incident_in_solar_band(grid, 340.25),
        layer_temperature=np.full(4, 250.0),
        layer_mass=np.full(4, 2500.0),
        absorption_coefficient=grey_absorption(grid, 1.0e-4),
        surface_temperature=288.0,
        surface_albedo=np.zeros(grid.n_bins),
    )
    a = solve_column_radiation(diffusivity_factor=1.66, **common)  # type: ignore[arg-type]
    b = solve_column_radiation(diffusivity_factor=2.00, **common)  # type: ignore[arg-type]
    assert a.interface_up(a.toa_index) != pytest.approx(b.interface_up(b.toa_index), rel=1e-6)
