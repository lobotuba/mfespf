"""Stratospheric volcanic aerosol as a scattering layer.

Requirements RAD-001, SPC-001, EQR-002, MOD-001, MOD-002, VAL-001, VAL-002.

The forcing is checkable against a published number, which the carbon dioxide
band's was too -- but here the record *is* the optical property, so far less is
constructed and the check is correspondingly sharper. What the tests below fix
is the sign, the confinement of the effect to where it belongs, the combination
with any other scatterer in the same layer, and the magnitude against Pinatubo.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError
from mfespf.physics.atmosphere import hydrostatic_layer_masses
from mfespf.physics.greenhouse import co2_band_absorption, greenhouse_grid
from mfespf.physics.radiative_transfer import combine_scatterers, solve_column_radiation
from mfespf.physics.volcanic import (
    SULFATE_ASYMMETRY,
    PrescribedVolcanicModule,
    stratospheric_aerosol_optics,
)

pytestmark = pytest.mark.stage(6)

GRID = greenhouse_grid(n_background=24, n_band=32, band_lower_m=6.5e-6, band_upper_m=5.0e-5)
N_LAYERS = 8
_, MASS = hydrostatic_layer_masses(
    n_layers=N_LAYERS, surface_pressure=101325.0, top_pressure=1.0e3, gravity=9.80665
)
PROFILE = np.linspace(310.0, 240.0, N_LAYERS)
SHORTWAVE = GRID.centers < 4.0e-6


def _toa_net(optical_depth: float) -> float:
    weights = np.where(SHORTWAVE, GRID.widths, 0.0)
    incident = np.where(
        SHORTWAVE, 340.25 * weights / np.sum(weights) / np.maximum(GRID.widths, 1e-300), 0.0
    )
    absorption, scattering, asymmetry = stratospheric_aerosol_optics(GRID, optical_depth, N_LAYERS)
    solution = solve_column_radiation(
        grid=GRID,
        incident_spectral_irradiance=incident,
        layer_temperature=PROFILE,
        layer_mass=MASS,
        absorption_coefficient=co2_band_absorption(GRID, 350.0, background=3.0e-5),
        surface_temperature=315.0,
        surface_albedo=np.where(SHORTWAVE, 0.3, 0.0),
        absorption_optical_depth=absorption,
        scattering_optical_depth=scattering,
        asymmetry=asymmetry,
    )
    return float(np.sum(solution.downward[-1] - solution.upward[-1]))


# --- the optical properties -------------------------------------------------


@pytest.mark.req("VAL-001")
def test_zero_optical_depth_scatters_nothing():
    absorption, scattering, asymmetry = stratospheric_aerosol_optics(GRID, 0.0, N_LAYERS)
    assert np.all(scattering == 0.0)
    assert np.all(absorption == 0.0)
    assert np.all(asymmetry == 0.0)  # undefined where nothing scatters, so not misleading
    assert _toa_net(0.0) == _toa_net(0.0)


@pytest.mark.req("SPC-001", "RAD-001")
def test_the_aerosol_is_confined_to_the_shortwave_and_the_stratosphere():
    """Two confinements, and both are sign statements rather than bookkeeping.

    Scattering in the longwave would reflect the surface's own emission back
    down and *warm* it, reversing the effect. Aerosol in the lowest layers
    would be tropospheric aerosol, which is a different mechanism with a
    different lifetime.
    """
    absorption, scattering, _g = stratospheric_aerosol_optics(
        GRID, 0.15, N_LAYERS, stratospheric_layers=2
    )
    assert np.all(scattering[:, ~SHORTWAVE] == 0.0)
    assert np.all(scattering[:-2, :] == 0.0)
    assert np.max(scattering[-2:, SHORTWAVE]) > 0.0
    # Pure scattering: sulfate absorbs almost nothing in the shortwave, and
    # representing it as partly absorbing would put energy into the layer that
    # in reality leaves the planet.
    assert np.all(absorption == 0.0)


@pytest.mark.req("EQR-002")
def test_extinction_does_not_keep_rising_into_the_ultraviolet():
    """Post-eruption droplets are large compared with visible light.

    Their extinction efficiency has reached the geometric-optics limit by 550
    nm and stops rising as the wavelength falls. An uncapped power law put two
    or three times the reference optical depth into the blue, where there is a
    great deal of solar energy, and the modelled cooling came out about twice
    what the record implies.
    """
    _a, scattering, _g = stratospheric_aerosol_optics(GRID, 0.2, N_LAYERS)
    top = scattering[-1]
    assert float(np.max(top)) == pytest.approx(0.2, rel=1e-9)
    # Longer wavelengths still fall away.
    longer = np.logical_and(SHORTWAVE, GRID.centers > 1.5e-6)
    assert np.all(top[longer] < 0.2)


@pytest.mark.req("EQR-002")
def test_the_optical_depth_is_linear_in_the_record():
    _a, single, _g = stratospheric_aerosol_optics(GRID, 0.1, N_LAYERS)
    _a2, double, _g2 = stratospheric_aerosol_optics(GRID, 0.2, N_LAYERS)
    assert np.allclose(double, 2.0 * single, rtol=1e-12)


@pytest.mark.req("SYS-006", "MOD-002")
def test_unphysical_aerosol_parameters_are_refused():
    with pytest.raises(ConfigurationError, match="cannot be negative"):
        stratospheric_aerosol_optics(GRID, -0.1, N_LAYERS)
    with pytest.raises(ConfigurationError, match="layers"):
        stratospheric_aerosol_optics(GRID, 0.1, N_LAYERS, stratospheric_layers=99)
    with pytest.raises(ConfigurationError, match="asymmetry"):
        stratospheric_aerosol_optics(GRID, 0.1, N_LAYERS, asymmetry=1.0)


# --- combining two scatterers ----------------------------------------------


@pytest.mark.req("RAD-001", "CLD-002")
def test_two_scatterers_in_one_layer_add_and_their_asymmetry_is_weighted():
    """Optical depths add; asymmetry does not.

    Asymmetry is a property of the scattered light, so the combination is
    weighted by how much each contributor scatters. A thin aerosol beneath a
    thick cloud must not drag the layer's asymmetry towards its own value.
    """
    shape = (N_LAYERS, GRID.n_bins)
    cloud = (np.zeros(shape), np.full(shape, 4.0), np.full(shape, 0.85))
    aerosol = (np.zeros(shape), np.full(shape, 1.0), np.full(shape, 0.6))

    absorption, scattering, asymmetry = combine_scatterers(cloud, aerosol)
    assert np.all(scattering == 5.0)
    assert np.all(absorption == 0.0)
    expected = (4.0 * 0.85 + 1.0 * 0.6) / 5.0
    assert np.allclose(asymmetry, expected)
    # Nearer the cloud's value than the aerosol's, because the cloud scatters more.
    assert abs(expected - 0.85) < abs(expected - 0.6)


@pytest.mark.req("RAD-001")
def test_combining_where_nothing_scatters_gives_zero_rather_than_a_division():
    shape = (2, 3)
    empty = (np.zeros(shape), np.zeros(shape), np.full(shape, 0.85))
    _a, scattering, asymmetry = combine_scatterers(empty, empty)
    assert np.all(scattering == 0.0)
    assert np.all(asymmetry == 0.0)
    with pytest.raises(ConfigurationError, match="at least one"):
        combine_scatterers()


# --- the forcing ------------------------------------------------------------


@pytest.mark.req("VAL-002", "CLM-003", "RAD-001")
def test_a_volcanic_layer_reflects_and_cools():
    """Sign first, then magnitude against a published number.

    An optical depth of 0.15 is roughly what Pinatubo produced, and the
    accepted forcing from it is near -3 W m^-2. Nothing here was calibrated to
    that: the record supplies the optical depth and the two-stream solver does
    the rest.
    """
    base = _toa_net(0.0)
    forcings = [_toa_net(tau) - base for tau in (0.05, 0.10, 0.15)]

    assert all(f < 0.0 for f in forcings)  # it cools
    assert forcings[0] > forcings[1] > forcings[2]  # more aerosol, more cooling
    assert forcings[2] == pytest.approx(-2.9, abs=0.5)
    # Nearly linear at these optical depths, which is what a thin scattering
    # layer does; the saturation a thick one would show is far away.
    assert forcings[2] / forcings[0] == pytest.approx(3.0, rel=0.05)


@pytest.mark.req("EQR-002", "CLM-003")
def test_a_cloud_like_asymmetry_would_halve_the_reflection():
    """Why the asymmetry is a declared parameter and not an afterthought.

    The fraction scattered back to space goes as roughly ``(1-g)/2``, so taking
    a cloud droplet's 0.85 instead of sulfate's 0.6 would lose most of the
    effect for the same optical depth.
    """
    weights = np.where(SHORTWAVE, GRID.widths, 0.0)
    incident = np.where(
        SHORTWAVE, 340.25 * weights / np.sum(weights) / np.maximum(GRID.widths, 1e-300), 0.0
    )

    def net(asymmetry: float) -> float:
        absorption, scattering, g = stratospheric_aerosol_optics(
            GRID, 0.15, N_LAYERS, asymmetry=asymmetry
        )
        solution = solve_column_radiation(
            grid=GRID,
            incident_spectral_irradiance=incident,
            layer_temperature=PROFILE,
            layer_mass=MASS,
            absorption_coefficient=co2_band_absorption(GRID, 350.0, background=3.0e-5),
            surface_temperature=315.0,
            surface_albedo=np.where(SHORTWAVE, 0.3, 0.0),
            absorption_optical_depth=absorption,
            scattering_optical_depth=scattering,
            asymmetry=g,
        )
        return float(np.sum(solution.downward[-1] - solution.upward[-1]))

    sulfate = net(SULFATE_ASYMMETRY) - _toa_net(0.0)
    cloudlike = net(0.85) - _toa_net(0.0)
    assert abs(cloudlike) < 0.6 * abs(sulfate)


# --- the module -------------------------------------------------------------


@pytest.mark.req("MOD-001", "MOD-002")
def test_the_module_declares_the_longwave_effect_it_omits():
    module = PrescribedVolcanicModule(optical_depth=0.1)
    module.validate_configuration()
    limitations = " ".join(module.contract.known_limitations)

    assert "longwave" in limitations
    assert "stratospheric temperature" in limitations
    assert "hemispheric" in limitations
    assert module.contract.owned_state_ids == frozenset()
    assert module.contract.produced_outputs == ("stratospheric_aerosol_optical_depth",)
    with pytest.raises(ConfigurationError, match="cannot be negative"):
        PrescribedVolcanicModule(optical_depth=-0.1).validate_configuration()


@pytest.mark.req("MOD-006", "AUD-002")
def test_a_time_varying_optical_depth_is_republished_every_step():
    from mfespf.contracts.module import StepContext

    module = PrescribedVolcanicModule(optical_depth=lambda t: 0.01 + t / 1.0e11)
    module.initialize({}, {"master_start_time": 0.0})
    assert float(module.emit_state()["stratospheric_aerosol_optical_depth"].values) == (
        pytest.approx(0.01)
    )

    module.begin_step(StepContext(time=1.0e11, dt=1.0e6, step_index=1))
    module.compute_tendencies()
    published = module.emit_outputs()["stratospheric_aerosol_optical_depth"]
    assert float(published.values) == pytest.approx(1.01)
    assert module.parameters()["longwave_effect"] == "omitted"
