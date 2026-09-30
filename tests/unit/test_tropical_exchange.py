"""The tropical column's exchange with the rest of the planet.

Checked against physics rather than against the code: mass is conserved layer
by layer, rising air in an isothermal layer cools at the dry adiabatic rate
g / c_p per metre it rises, the band's sunlight is geometry that reduces to one
quarter for the whole globe, no exchange changes nothing, and a coupled run's
budgets close through the lateral boundary.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.constants import value_of
from mfespf.physics.orbit import band_insolation_factor
from mfespf.physics.tropical_exchange import (
    export_heating,
    upward_advection_rate,
    upwelling_mass_flux,
)

pytestmark = pytest.mark.stage(4)

G = value_of("earth_standard_gravity")
R = value_of("gas_constant_dry_air")
CP = value_of("specific_heat_dry_air")
F0 = 2.737e-5


def log_grid(n: int = 400, top: float = 100.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    interfaces = np.geomspace(101325.0, top, n + 1)
    centres = np.sqrt(interfaces[:-1] * interfaces[1:])
    mass = (interfaces[:-1] - interfaces[1:]) / G
    return interfaces, centres, mass


@pytest.mark.req("CON-001", "EQR-002")
def test_the_upwelling_conserves_mass_layer_by_layer():
    interfaces, _, _ = log_grid()
    flux = upwelling_mass_flux(interfaces, flux_at_reference=F0)
    assert flux[0] == 0.0 and flux[-1] == 0.0
    assert np.all(flux >= 0.0)
    assert np.all(flux[(interfaces <= 1.0e4) & (interfaces >= 7.0e3)] == F0)
    assert np.all(flux[interfaces > 1.0e4] == 0.0)
    # Whatever enters laterally leaves laterally: the net lateral mass is zero.
    lateral = flux[1:] - flux[:-1]
    assert float(np.sum(lateral)) == pytest.approx(0.0, abs=1e-18)
    # A uniform quantity is not changed by being moved.
    assert np.all(upward_advection_rate(flux, np.full(400, 3.0)) == 0.0)


@pytest.mark.req("EQR-002", "VAL-001")
def test_rising_air_cools_at_the_dry_adiabatic_rate_in_an_isothermal_layer():
    """dT/dt = -w g / c_p with w = F / rho, where the flux is uniform."""
    temperature = 200.0
    interfaces, centres, mass = log_grid(4000)
    height = R * temperature / G * np.log(101325.0 / centres)
    s = CP * temperature + G * height
    flux = upwelling_mass_flux(interfaces, flux_at_reference=F0)
    warming = upward_advection_rate(flux, s) / (CP * mass)
    # Well inside the uniform-flux band, away from its edges.
    inside = (centres < 9.5e3) & (centres > 7.5e3)
    rho = centres[inside] / (R * temperature)
    expected = -(F0 / rho) * G / CP
    np.testing.assert_allclose(warming[inside], expected, rtol=2e-3)


@pytest.mark.req("SOL-001", "EQR-002")
def test_the_band_insolation_is_geometry():
    assert band_insolation_factor(90.0) == pytest.approx(0.25, abs=1e-5)
    # Nearer the equator a band gets more; these are the protocol's numbers.
    assert band_insolation_factor(25.0) == pytest.approx(0.29681, abs=2e-5)
    assert (
        band_insolation_factor(20.0) > band_insolation_factor(25.0) > band_insolation_factor(30.0)
    )


@pytest.mark.req("EQR-002", "VAL-001")
def test_no_exchange_leaves_the_column_unchanged():
    interfaces, centres, mass = log_grid(50)
    flux = upwelling_mass_flux(interfaces, flux_at_reference=0.0)
    assert np.all(flux == 0.0)
    assert np.all(upward_advection_rate(flux, np.linspace(0.0, 1.0, 50)) == 0.0)
    assert np.all(export_heating(mass, centres, 0.0) == 0.0)
    # And a real export removes exactly what it says, from below 100 hPa only.
    heating = export_heating(mass, centres, 38.464)
    assert float(np.sum(heating)) == pytest.approx(-38.464, rel=1e-12)
    assert np.all(heating[centres < 1.0e4] == 0.0)


@pytest.mark.req("CON-001", "CON-004")
def test_the_exchange_budget_closes_through_the_lateral_boundary():
    from mfespf.run.column_configurations import column_grid
    from mfespf.run.tropical_column import tropical_column

    run = tropical_column(
        grid=column_grid(16),
        areal_heat_capacity=1.0e6,
        measured_solar_spectrum=True,
        ozone=True,
        eddy_diffusivity=0.39,
        max_log_step=0.2,
        top_pressure=100.0,
    ).run(200)
    assert run.may_support_a_claim, run.manifest.invalidation_reasons
    assert run.worst_relative_residual() < 1e-10
    exported = float(np.asarray(run.steps[-1].outputs["lateral_energy_export"].values))
    assert exported > 0.0
