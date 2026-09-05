"""Stage 4: two-stream scattering and cloud optical properties.

The two-stream limits come first. Scattering was introduced underneath work that
had already been validated without it, so the first thing that must be shown is
that it changes nothing where it should not.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError
from mfespf.physics.clouds import (
    CloudProfile,
    cloud_optical_properties,
    single_layer_cloud,
)
from mfespf.physics.radiative_transfer import (
    grey_absorption,
    solve_column_radiation,
    two_stream_layer,
)
from mfespf.run.column_configurations import column_grid

pytestmark = pytest.mark.stage(4)

GRID = column_grid(16)
DIFFUSIVITY = 1.66


def _incident(grid, total: float = 340.25) -> np.ndarray:
    values = np.zeros(grid.n_bins)
    values[0] = total / grid.widths[0]
    return values


# --- two-stream limits -----------------------------------------------------


@pytest.mark.req("RAD-002", "VAL-001")
def test_two_stream_reduces_to_pure_absorption():
    """With no scattering the layer coefficients must be the old ones exactly."""
    tau = np.array([0.0, 1e-6, 0.01, 0.5, 2.0, 50.0, 500.0])
    zero = np.zeros_like(tau)
    reflectance, transmittance, emissivity = two_stream_layer(tau, zero, zero, DIFFUSIVITY)

    assert np.array_equal(reflectance, zero)
    assert np.allclose(transmittance, np.exp(-DIFFUSIVITY * tau), rtol=1e-14, atol=0.0)
    assert np.allclose(emissivity, 1.0 - transmittance, rtol=0.0, atol=0.0)


@pytest.mark.req("RAD-002", "VAL-001")
def test_conservative_scattering_neither_absorbs_nor_emits():
    """A purely scattering layer must satisfy r + t = 1 and emit nothing."""
    tau = np.array([0.01, 0.5, 2.0, 50.0, 500.0])
    ones = np.ones_like(tau)
    reflectance, transmittance, emissivity = two_stream_layer(
        tau, ones, np.zeros_like(tau), DIFFUSIVITY
    )

    assert np.max(np.abs(reflectance + transmittance - 1.0)) <= 2.3e-16
    assert np.max(emissivity) <= 2.3e-16
    assert np.all(np.diff(reflectance) > 0.0)  # thicker scatters more


@pytest.mark.req("RAD-002")
def test_delta_scaling_makes_forward_scattering_less_reflective():
    """Forward-peaked scattering must reflect less than isotropic scattering.

    Without delta scaling a two-stream treats a droplet's strong forward lobe as
    though half of it went backwards, and a cloud comes out far too bright.
    """
    tau = np.array([10.0])
    isotropic = two_stream_layer(tau, np.array([1.0]), np.array([0.0]), DIFFUSIVITY)[0]
    forward = two_stream_layer(tau, np.array([1.0]), np.array([0.85]), DIFFUSIVITY)[0]
    assert forward[0] < isotropic[0]
    assert 0.3 < forward[0] < 0.8  # a plausible albedo for an optical depth of ten


@pytest.mark.req("RAD-002", "CON-009")
def test_scattering_solver_reproduces_the_absorption_only_results():
    """The adding method must agree with the sweeps it replaced.

    Checked against an independently computed absorption-only solution rather
    than against a stored number, so it tests the algebra and not a snapshot.
    """
    n_layers = 5
    temperature = np.array([280.0, 265.0, 250.0, 235.0, 220.0])
    mass = np.full(n_layers, 2000.0)
    kappa = grey_absorption(GRID, 2.0e-4)
    albedo = np.full(GRID.n_bins, 0.15)

    solution = solve_column_radiation(
        grid=GRID,
        incident_spectral_irradiance=_incident(GRID),
        layer_temperature=temperature,
        layer_mass=mass,
        absorption_coefficient=kappa,
        surface_temperature=290.0,
        surface_albedo=albedo,
    )

    # Independent reference: the explicit downward-then-upward sweeps, which are
    # exact when nothing scatters.
    from mfespf.physics.radiation import band_exitance, band_exitance_profile

    optical_depth = DIFFUSIVITY * kappa[None, :] * mass[:, None]
    emissivity = -np.expm1(-optical_depth)
    transmittance = 1.0 - emissivity
    emission = band_exitance_profile(GRID, temperature)

    downward = np.zeros((n_layers + 1, GRID.n_bins))
    downward[n_layers] = _incident(GRID) * GRID.widths
    for i in range(n_layers - 1, -1, -1):
        downward[i] = downward[i + 1] * transmittance[i] + emissivity[i] * emission[i]
    upward = np.zeros_like(downward)
    upward[0] = (1.0 - albedo) * band_exitance(GRID, 290.0) + albedo * downward[0]
    for i in range(n_layers):
        upward[i + 1] = upward[i] * transmittance[i] + emissivity[i] * emission[i]

    assert np.allclose(solution.downward, downward, rtol=1e-12)
    assert np.allclose(solution.upward, upward, rtol=1e-12)
    assert solution.energy_residual() == pytest.approx(0.0, abs=1e-9)


@pytest.mark.req("RAD-002", "CON-009")
def test_flux_divergence_still_closes_with_scattering_present():
    """Adding scattering must not disturb the identity the budget relies on."""
    profile = single_layer_cloud(n_layers=6, layer=3, liquid_water_path=0.08, fraction=1.0)
    optics = cloud_optical_properties(GRID, profile)
    solution = solve_column_radiation(
        grid=GRID,
        incident_spectral_irradiance=_incident(GRID),
        layer_temperature=np.array([285.0, 275.0, 262.0, 250.0, 238.0, 225.0]),
        layer_mass=np.full(6, 1700.0),
        absorption_coefficient=grey_absorption(GRID, 1.0e-4),
        surface_temperature=292.0,
        surface_albedo=np.full(GRID.n_bins, 0.1),
        absorption_optical_depth=optics.absorption,
        scattering_optical_depth=optics.scattering,
        asymmetry=optics.asymmetry,
    )
    assert solution.energy_residual() == pytest.approx(0.0, abs=1e-9)


# --- cloud optics ----------------------------------------------------------


@pytest.mark.req("CLD-001", "CLD-002")
def test_cloud_state_is_a_separate_variable_set():
    """CLD-001: liquid, ice, fraction and vertical support, all addressable."""
    profile = CloudProfile(
        liquid_water_path=np.array([0.0, 0.05, 0.0]),
        ice_water_path=np.array([0.0, 0.0, 0.01]),
        fraction=0.5,
    )
    assert profile.n_layers == 3
    assert profile.total_water_path == pytest.approx(0.06)
    assert not profile.is_clear
    # The cloud lives in identifiable layers, not in a column-mean number.
    assert profile.liquid_water_path[1] > 0.0 and profile.liquid_water_path[0] == 0.0


@pytest.mark.req("CLD-002", "EQR-002")
def test_cloud_optical_depth_scales_with_water_path():
    """Optical depth must be linear in water path at fixed effective radius."""
    depths = []
    for water in (0.01, 0.02, 0.04):
        profile = single_layer_cloud(n_layers=3, layer=1, liquid_water_path=water, fraction=1.0)
        optics = cloud_optical_properties(GRID, profile)
        depths.append(float(np.max(optics.scattering[1])))
    assert depths[1] == pytest.approx(2.0 * depths[0], rel=1e-12)
    assert depths[2] == pytest.approx(4.0 * depths[0], rel=1e-12)

    # And inversely proportional to the effective radius, which is why that
    # parameter sets the cloud albedo as much as the water path does.
    small = cloud_optical_properties(
        GRID,
        single_layer_cloud(
            n_layers=3,
            layer=1,
            liquid_water_path=0.02,
            fraction=1.0,
            droplet_effective_radius=5.0e-6,
        ),
    )
    assert float(np.max(small.scattering[1])) == pytest.approx(2.0 * depths[1], rel=1e-12)


@pytest.mark.req("CLD-002", "SPC-010")
def test_cloud_scatters_in_the_shortwave_and_absorbs_in_the_longwave():
    """The two effects come from different physics and must stay distinguishable."""
    profile = single_layer_cloud(n_layers=3, layer=1, liquid_water_path=0.05, fraction=1.0)
    optics = cloud_optical_properties(GRID, profile)
    shortwave = GRID.centers < 4.0e-6

    assert np.all(optics.scattering[1][shortwave] > 0.0)
    assert np.all(optics.scattering[1][~shortwave] == 0.0)
    assert np.all(optics.absorption[1][~shortwave] > 0.0)
    assert np.all(optics.absorption[1][shortwave] == 0.0)
    # Cloud-free layers stay cloud-free.
    assert np.all(optics.scattering[0] == 0.0) and np.all(optics.absorption[0] == 0.0)


@pytest.mark.req("CLD-002", "VAL-001")
def test_a_clear_profile_has_no_optical_effect():
    """The null case: zero fraction or zero water gives identically nothing."""
    for profile in (
        single_layer_cloud(n_layers=4, layer=2, liquid_water_path=0.05, fraction=0.0),
        single_layer_cloud(n_layers=4, layer=2, liquid_water_path=0.0, fraction=1.0),
    ):
        optics = cloud_optical_properties(GRID, profile)
        assert optics.is_empty
        assert np.array_equal(optics.scattering, np.zeros_like(optics.scattering))
        assert np.array_equal(optics.absorption, np.zeros_like(optics.absorption))


@pytest.mark.req("CLD-002")
def test_ice_is_less_extinguishing_per_kilogram_than_liquid():
    """Larger crystals give less shortwave extinction for the same mass."""
    liquid = cloud_optical_properties(
        GRID,
        single_layer_cloud(n_layers=3, layer=1, liquid_water_path=0.02, fraction=1.0),
    )
    ice = cloud_optical_properties(
        GRID,
        single_layer_cloud(n_layers=3, layer=1, liquid_water_path=0.02, fraction=1.0, ice=True),
    )
    assert float(np.max(ice.scattering[1])) < float(np.max(liquid.scattering[1]))
    # But their longwave absorption is the same: it depends on mass, not radius.
    assert float(np.max(ice.absorption[1])) == pytest.approx(
        float(np.max(liquid.absorption[1])), rel=1e-12
    )


@pytest.mark.req("SYS-006", "MOD-002")
def test_unphysical_cloud_states_are_refused():
    with pytest.raises(ConfigurationError, match="fraction"):
        CloudProfile(liquid_water_path=np.zeros(3), fraction=1.5)
    with pytest.raises(ConfigurationError, match="negative"):
        CloudProfile(liquid_water_path=np.array([-0.1, 0.0]))
    with pytest.raises(ConfigurationError, match="same shape"):
        CloudProfile(liquid_water_path=np.zeros(3), ice_water_path=np.zeros(4))
    with pytest.raises(ConfigurationError, match="outside a column"):
        single_layer_cloud(n_layers=3, layer=7, liquid_water_path=0.01, fraction=1.0)
