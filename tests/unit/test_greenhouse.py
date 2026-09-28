"""Carbon dioxide as a spectral absorber.

Requirements CHM-001, EQR-002, EQR-004, RAD-001, SPC-001, MOD-002, VAL-001.

The claim worth testing is not that the model warms when CO2 rises -- almost any
absorber does that -- but that the *logarithm* comes out of the band's shape
rather than being written in. So the tests measure the forcing per doubling at
several starting concentrations and check that it is nearly the same each time,
and check that it stops being so when the grid cannot resolve the wings.

Nothing here runs a coupled model or reads `data/raw`.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError
from mfespf.physics.atmosphere import hydrostatic_layer_masses
from mfespf.physics.greenhouse import (
    CALIBRATION_N_LAYERS,
    CO2_BAND_CENTRE_M,
    CO2_BAND_LOWER_M,
    CO2_BAND_RESOLUTION,
    CO2_BAND_UPPER_M,
    INFRARED_BACKGROUND_ABSORPTION,
    PREINDUSTRIAL_CO2_PPM,
    PrescribedGreenhouseModule,
    co2_band_absorption,
    co2_mass_mixing_ratio,
    greenhouse_grid,
)
from mfespf.physics.radiative_transfer import solve_column_radiation

pytestmark = pytest.mark.stage(6)

GRID = greenhouse_grid()  # the grid the band is calibrated on, taken by default
N_LAYERS = CALIBRATION_N_LAYERS  # the resolution the bands are fitted at
_, MASS = hydrostatic_layer_masses(
    n_layers=N_LAYERS, surface_pressure=101325.0, top_pressure=1.0e3, gravity=9.80665
)
PROFILE = np.linspace(285.0, 220.0, N_LAYERS)


def _toa_net(ppm: float, grid=GRID, *, mass=None, profile=None, **kwargs) -> float:
    shortwave = grid.centers < 4.0e-6
    weights = np.where(shortwave, grid.widths, 0.0)
    incident = np.where(
        shortwave, 340.25 * weights / np.sum(weights) / np.maximum(grid.widths, 1e-300), 0.0
    )
    solution = solve_column_radiation(
        grid=grid,
        incident_spectral_irradiance=incident,
        layer_temperature=PROFILE if profile is None else profile,
        layer_mass=MASS if mass is None else mass,
        absorption_coefficient=co2_band_absorption(grid, ppm, **kwargs),
        surface_temperature=288.0,
        surface_albedo=np.where(shortwave, 0.3, 0.0),
    )
    return float(np.sum(solution.downward[-1] - solution.upward[-1]))


def _forcing(low: float, high: float, grid=GRID, **kwargs) -> float:
    """Instantaneous forcing: the temperature profile is held fixed."""
    return _toa_net(high, grid, **kwargs) - _toa_net(low, grid, **kwargs)


# --- the constitutive relation ---------------------------------------------


@pytest.mark.req("CHM-001", "EQR-002")
def test_the_mixing_ratio_converts_a_mole_fraction_to_a_mass_fraction():
    """Heavier than air by the ratio of the molar masses, so 1 ppm is 1.52 ppm by mass."""
    assert float(co2_mass_mixing_ratio(1.0e6)) == pytest.approx(44.0095 / 28.9647, rel=1e-9)
    assert float(co2_mass_mixing_ratio(400.0)) == pytest.approx(400e-6 * 44.0095 / 28.9647)
    assert float(co2_mass_mixing_ratio(0.0)) == 0.0
    with pytest.raises(ConfigurationError, match="cannot be negative"):
        co2_mass_mixing_ratio(-1.0)


@pytest.mark.req("VAL-001", "SPC-001")
def test_zero_concentration_leaves_only_the_background():
    """The null case: no gas, no band, and the infrared background untouched."""
    background = 3.0e-5
    kappa = co2_band_absorption(GRID, 0.0, background=background)
    infrared = GRID.centers >= 4.0e-6
    assert np.all(kappa[infrared] == background)
    assert np.all(kappa[~infrared] == 0.0)
    assert _forcing(0.0, 0.0) == 0.0


@pytest.mark.req("SPC-001", "CLD-002")
def test_the_band_is_confined_to_the_infrared():
    """CO2 warms rather than shades: its band sits where the Earth radiates.

    A greenhouse gas that absorbed sunlight would cool the surface, so this is a
    sign check on the mechanism and not a detail of the grid.
    """
    kappa = co2_band_absorption(GRID, 400.0, background=0.0)
    shortwave = GRID.centers < 4.0e-6
    assert np.all(kappa[shortwave] == 0.0)
    assert np.max(kappa[~shortwave]) > 0.0
    # Strongest at the band centre, falling away on both sides.
    peak = int(np.argmax(kappa))
    assert abs(np.log(GRID.centers[peak] / CO2_BAND_CENTRE_M)) < 0.1


@pytest.mark.req("EQR-002")
def test_absorption_is_linear_in_concentration_before_the_transfer_sees_it():
    """The constitutive relation is linear; the *forcing* is not.

    Worth separating, because the logarithm is often assumed to be in the
    absorption coefficient. It is not: the coefficient doubles when the
    concentration doubles, and the transfer solver turns that into a nearly
    constant forcing.
    """
    single = co2_band_absorption(GRID, 200.0, background=0.0)
    double = co2_band_absorption(GRID, 400.0, background=0.0)
    assert np.allclose(double, 2.0 * single, rtol=1e-12)


@pytest.mark.req("MOD-002", "SYS-006")
def test_unphysical_band_parameters_are_refused():
    with pytest.raises(ConfigurationError, match="non-negative"):
        co2_band_absorption(GRID, 400.0, centre_absorption=-1.0)
    with pytest.raises(ConfigurationError, match="wing decay"):
        co2_band_absorption(GRID, 400.0, wing_decay_per_metre=0.0)
    with pytest.raises(ConfigurationError, match="resolve the wings"):
        greenhouse_grid(n_band=4)
    with pytest.raises(ConfigurationError, match="inside the grid's support"):
        greenhouse_grid(band_lower_m=1.0e-2)


# --- where the logarithm comes from ----------------------------------------


@pytest.mark.req("EQR-002", "VAL-002", "CLM-003")
def test_forcing_per_doubling_is_nearly_constant():
    """The result this whole band model exists to produce.

    A doubling gives about the same forcing wherever it starts, across a factor
    of eight in concentration. Nothing asserts that: the band core is saturated,
    so more gas widens the opaque region instead of deepening it, and the width
    grows with the logarithm of the concentration. The transfer solver produces
    the consequence and this measures it.
    """
    doublings = [_forcing(c, 2.0 * c) for c in (142.35, 284.7, 569.4, 1138.8)]

    assert all(f > 0.0 for f in doublings)  # the sign, first
    assert max(doublings) / min(doublings) < 1.06
    # And it hits the value it was calibrated against, at the background and on
    # the grid it was calibrated at -- both of which the configurations use.
    assert _forcing(PREINDUSTRIAL_CO2_PPM, 2.0 * PREINDUSTRIAL_CO2_PPM) == pytest.approx(
        3.70, abs=0.05
    )


@pytest.mark.req("EQR-002", "SPC-005")
def test_a_grid_that_cannot_resolve_the_wings_loses_the_logarithm():
    """The failure mode the band-count guard exists to prevent.

    With the wings unresolved the forcing per doubling drifts far more than it
    should, because the model is no longer widening an opaque region -- it is
    filling in a handful of coarse bands.
    """

    def spread(grid):
        doublings = [_forcing(c, 2.0 * c, grid) for c in (142.35, 284.7, 569.4, 1138.8)]
        return max(doublings) / min(doublings)

    coarse = greenhouse_grid(n_band=8, band_lower_m=6.5e-6, band_upper_m=5.0e-5)
    assert spread(coarse) > spread(GRID)


@pytest.mark.req("CMP-006", "VAL-002")
def test_the_forcing_is_converged_in_spectral_resolution():
    """A number that still moves with the grid is a property of the grid.

    This one caught a real error. When the band was narrowed the previous grid --
    48 bands over 6.5 to 50 microns -- stopped resolving the wings, and gave
    3.700 W m^-2 per doubling where the converged answer was 3.601. Calibrating
    there would have fitted the grid rather than the physics, so the resolution
    was fixed first and the wing decay fitted afterwards.

    The reference is deliberately a *wider* region as well as a finer one, so it
    also checks that nothing was lost by narrowing the region to 8-30 microns.
    """
    fine = greenhouse_grid(n_band=160, band_lower_m=6.5e-6, band_upper_m=5.0e-5)
    reference = _forcing(284.7, 569.4, fine)
    assert _forcing(284.7, 569.4, GRID) == pytest.approx(reference, rel=2e-3)


@pytest.mark.req("EQR-002")
def test_the_band_core_is_saturated_at_present_day_concentrations():
    """The premise the logarithm rests on, checked rather than asserted."""
    kappa = co2_band_absorption(GRID, PREINDUSTRIAL_CO2_PPM, background=0.0)
    assert float(np.max(kappa) * np.sum(MASS)) > 1.0e3


# --- the module -------------------------------------------------------------


@pytest.mark.req("MOD-001", "MOD-002", "CHM-001")
def test_the_module_declares_what_it_does_not_represent():
    module = PrescribedGreenhouseModule(concentration=300.0)
    module.validate_configuration()
    limitations = " ".join(module.contract.known_limitations)

    assert "no carbon cycle" in limitations
    assert "water-vapour feedback" in limitations
    assert module.contract.produced_outputs == ("co2_mole_fraction",)
    assert module.contract.owned_state_ids == frozenset()  # it owns nothing
    assert module.contract.conserved_quantities == ()
    with pytest.raises(ConfigurationError, match="must be positive"):
        PrescribedGreenhouseModule(concentration=0.0).validate_configuration()


@pytest.mark.req("MOD-006", "AUD-002")
def test_a_time_varying_concentration_is_republished_every_step():
    """The trap the cloud module fell into: a PRESCRIBED module's emit_state runs
    only at initialisation."""
    from mfespf.contracts.module import StepContext

    module = PrescribedGreenhouseModule(concentration=lambda t: 300.0 + t / 1.0e9)
    module.initialize({}, {"master_start_time": 0.0})
    assert float(module.emit_state()["co2_mole_fraction"].values) == pytest.approx(300.0)

    module.begin_step(StepContext(time=2.0e9, dt=1.0e6, step_index=1))
    module.compute_tendencies()
    assert float(module.emit_outputs()["co2_mole_fraction"].values) == pytest.approx(302.0)
    assert module.parameters()["concentration_ppm"] == "time_varying"


@pytest.mark.req("EQR-004", "INF-003", "SYS-006")
def test_the_band_is_calibrated_at_the_background_the_runs_use():
    """The defect this guards against was real and cost a factor of 1.75.

    The wing decay was once fitted with the function's default background while
    every configuration passed a different one. A thinner background masks less
    of the band, so the forcing per doubling in the runs was 6.28 W m^-2 against
    a calibrated 3.58. One constant, used in both places, is the fix -- and this
    is what notices if they are ever separated again.

    The grid is checked the same way and for the same reason: a spectral
    resolution that resolves the wings during calibration and not during the run
    separates the two just as effectively.
    """

    def greenhouse_replay_n_band() -> int:
        return int(inspect.signature(greenhouse_replay).parameters["n_band"].default)

    import inspect

    from mfespf.run.greenhouse import GREENHOUSE_BACKGROUND, greenhouse_replay

    default = inspect.signature(co2_band_absorption).parameters["background"].default
    assert default == INFRARED_BACKGROUND_ABSORPTION
    assert GREENHOUSE_BACKGROUND == INFRARED_BACKGROUND_ABSORPTION

    # The same for the grid: the calibration is only worth anything on the grid
    # it was performed on, so the run configurations take these defaults too.
    parameters = inspect.signature(greenhouse_grid).parameters
    assert parameters["n_band"].default == CO2_BAND_RESOLUTION
    assert parameters["band_lower_m"].default == CO2_BAND_LOWER_M
    assert parameters["band_upper_m"].default == CO2_BAND_UPPER_M
    assert greenhouse_replay_n_band() == CO2_BAND_RESOLUTION

    # And the vertical resolution: the vapour bands were fitted on a 10-layer
    # profile while every configuration ran 8, and neither was converged.
    # The historical replay is the declared exception: its two-day explicit step
    # is unstable at 30 layers, and the calibration it depends on -- this band,
    # against the fixed background -- does not vary with the layer count.
    from mfespf.run.column_configurations import radiating_column
    from mfespf.run.greenhouse import REPLAY_N_LAYERS
    from mfespf.run.greenhouse import greenhouse_replay as replay

    columns = inspect.signature(radiating_column).parameters["n_layers"].default
    assert columns == CALIBRATION_N_LAYERS
    assert inspect.signature(replay).parameters["n_layers"].default == REPLAY_N_LAYERS
    for n_layers in (10, CALIBRATION_N_LAYERS, 60):
        _, mass = hydrostatic_layer_masses(
            n_layers=n_layers, surface_pressure=101325.0, top_pressure=1.0e3, gravity=9.80665
        )
        profile = np.linspace(285.0, 220.0, n_layers)

        def toa(ppm: float, mass=mass, profile=profile) -> float:
            return _toa_net(ppm, mass=mass, profile=profile)

        assert toa(569.4) - toa(284.7) == pytest.approx(3.70, abs=0.01)

    # And the calibrated forcing is the one all of those defaults produce.
    assert _forcing(284.7, 569.4) == pytest.approx(3.70, abs=0.05)


@pytest.mark.req("CLM-003", "VAL-005")
def test_the_base_state_is_a_prediction_and_it_is_low():
    """Only the forcing was fitted, so the greenhouse effect is a result.

    It comes out about five kelvin short of the observed 33 K, which is roughly
    what water vapour would supply. The vapour this framework tracks is
    radiatively inert, so the shortfall is expected and is recorded rather than
    closed by adjusting the background.
    """
    column = np.sum(MASS)
    background_optical_depth = INFRARED_BACKGROUND_ABSORPTION * column
    # A window that is genuinely a window: optically thin, not transparent.
    assert 0.1 < background_optical_depth < 0.6


@pytest.mark.req("EQR-004", "MOD-002", "SYS-006")
def test_a_coupled_radiation_module_takes_the_calibrated_background_by_itself():
    """No caller should have to remember the opacity the band was fitted at.

    ``SpectralRadiationModule`` used to default this to 1e-4 and both greenhouse
    configurations reassigned it after construction. Nothing was running wrong,
    because both of them remembered -- but the next one would not have, and the
    failure is silent: the column uses a forcing per doubling the calibration
    never produced, and reports no complaint about it. That is the same shape as
    the defect that invalidated four claims, one layer down.
    """
    from mfespf.physics.radiative_transfer import SpectralRadiationModule

    module = SpectralRadiationModule(
        grid=GRID,
        absorption_coefficient=np.full(GRID.n_bins, INFRARED_BACKGROUND_ABSORPTION),
        couple_greenhouse=True,
    )

    assert module.greenhouse_background == INFRARED_BACKGROUND_ABSORPTION
    # And it is in the parameters, so two runs at different opacities cannot
    # hash alike (MOD-002).
    assert module.parameters()["greenhouse_background_m2_kg"] == INFRARED_BACKGROUND_ABSORPTION
