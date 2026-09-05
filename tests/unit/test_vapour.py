"""Water vapour as a spectral absorber, and the feedback that follows.

Requirements EQR-002, EQR-004, RAD-001, SPC-001, VAL-001, VAL-002, INF-003.

The claim worth testing is not that a moist atmosphere is warmer -- any added
absorber does that -- but that the *feedback* is a consequence rather than an
input. Nothing in the constitutive relation knows about temperature. The water
module supplies a vapour path, Clausius-Clapeyron makes that path grow when the
column warms, and the amplification is what the radiative transfer computes from
those two facts. So the tests below hold the vapour fixed and let it follow, and
measure the difference.

Nothing here runs a coupled model or reads `data/raw`.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError
from mfespf.physics.atmosphere import hydrostatic_layer_masses
from mfespf.physics.greenhouse import (
    CO2_WING_DECAY_PER_METRE,
    CO2_WING_DECAY_WITH_VAPOUR_PER_METRE,
    VAPOUR_BAND_ABSORPTION,
    co2_band_absorption,
    co2_wing_decay_for,
    greenhouse_grid,
    vapour_absorption_optical_depth,
    vapour_mass_absorption,
)
from mfespf.physics.radiative_transfer import solve_column_radiation
from mfespf.physics.thermodynamics import saturation_specific_humidity

pytestmark = pytest.mark.stage(6)

GRID = greenhouse_grid()
N_LAYERS = 10
INTERFACES, MASS = hydrostatic_layer_masses(
    n_layers=N_LAYERS, surface_pressure=101325.0, top_pressure=1.0e3, gravity=9.80665
)
PRESSURE = 0.5 * (INTERFACES[:-1] + INTERFACES[1:])
PROFILE = np.linspace(288.0, 220.0, N_LAYERS)
SHORTWAVE = GRID.centers < 4.0e-6
RELATIVE_HUMIDITY = 0.7

_WEIGHTS = np.where(SHORTWAVE, GRID.widths, 0.0)
INCIDENT = np.where(
    SHORTWAVE, 340.25 * _WEIGHTS / np.sum(_WEIGHTS) / np.maximum(GRID.widths, 1e-300), 0.0
)


def _vapour(temperature: np.ndarray) -> np.ndarray:
    """The vapour path a column at constant relative humidity carries."""
    return RELATIVE_HUMIDITY * saturation_specific_humidity(temperature, PRESSURE) * MASS


def _band(wavelength_m: float) -> int:
    return int(np.argmin(np.abs(GRID.centers - wavelength_m)))


def _transmittance(vapour_path: np.ndarray, wavelength_m: float) -> float:
    tau = vapour_absorption_optical_depth(GRID, vapour_path).sum(axis=0)
    return float(np.exp(-tau[_band(wavelength_m)]))


def _solve(shift: float, *, follow: bool, ppm: float = 284.7):
    temperature = PROFILE + shift
    vapour = _vapour(temperature) if follow else _vapour(PROFILE)
    return solve_column_radiation(
        grid=GRID,
        incident_spectral_irradiance=INCIDENT,
        layer_temperature=temperature,
        layer_mass=MASS,
        absorption_coefficient=co2_band_absorption(
            GRID, ppm, wing_decay_per_metre=co2_wing_decay_for(active_vapour=True), background=0.0
        ),
        surface_temperature=288.0 + shift,
        surface_albedo=np.where(SHORTWAVE, 0.3, 0.0),
        absorption_optical_depth=vapour_absorption_optical_depth(GRID, vapour),
    )


def _olr_with(vapour_path: np.ndarray, shift: float, ppm: float = 284.7) -> float:
    """Outgoing longwave for a vapour distribution supplied directly."""
    solution = solve_column_radiation(
        grid=GRID,
        incident_spectral_irradiance=INCIDENT,
        layer_temperature=PROFILE + shift,
        layer_mass=MASS,
        absorption_coefficient=co2_band_absorption(
            GRID, ppm, wing_decay_per_metre=co2_wing_decay_for(active_vapour=True), background=0.0
        ),
        surface_temperature=288.0 + shift,
        surface_albedo=np.where(SHORTWAVE, 0.3, 0.0),
        absorption_optical_depth=vapour_absorption_optical_depth(GRID, vapour_path),
    )
    return float(np.sum(solution.upward[-1][~SHORTWAVE]))


def _olr(shift: float, *, follow: bool, ppm: float = 284.7) -> float:
    return float(np.sum(_solve(shift, follow=follow, ppm=ppm).upward[-1][~SHORTWAVE]))


# --- the constitutive relation ----------------------------------------------


@pytest.mark.req("VAL-001", "SPC-001")
def test_a_dry_column_absorbs_nothing_and_the_window_is_open():
    """The null case, and it is the one that says the effect is the vapour's.

    A column with no vapour has no vapour opacity anywhere -- not a small
    residue, exactly none -- so anything the moist column does differently is
    attributable to the water and not to a floor someone left in.
    """
    tau = vapour_absorption_optical_depth(GRID, np.zeros(N_LAYERS))

    assert tau.shape == (N_LAYERS, GRID.n_bins)
    assert np.all(tau == 0.0)
    assert _transmittance(np.zeros(N_LAYERS), 10.0e-6) == 1.0


@pytest.mark.req("SPC-001", "RAD-001")
def test_the_bands_are_saturated_and_confined_to_the_infrared():
    """Two sign statements, not bookkeeping.

    Vapour absorbing sunlight would warm the atmosphere directly rather than
    trap the surface's emission, which is a different mechanism with a different
    signature; this configuration does not represent it and so must not produce
    a trace of it. And the response has to come from the window and the wings,
    which is only true if the band cores are already saturated.
    """
    kappa = vapour_mass_absorption(GRID)
    assert np.all(kappa[SHORTWAVE] == 0.0)

    column = _vapour(PROFILE).sum()
    assert _transmittance(_vapour(PROFILE), 6.3e-6) < 0.02  # the bending band
    assert _transmittance(_vapour(PROFILE), 30.0e-6) < 0.02  # the rotational band
    # Saturated means saturated: an optical depth of order a thousand, so the
    # exact band-centre coefficient stops mattering.
    assert VAPOUR_BAND_ABSORPTION * column > 1.0e3


@pytest.mark.req("EQR-002", "CLM-003")
def test_the_window_closes_as_the_column_moistens():
    """Where the feedback actually comes from.

    The bands are saturated, so they cannot trap more. What changes with vapour
    amount is the window: it is the only part of the spectrum through which the
    surface still reaches space, and the continuum closes it as the atmosphere
    moistens. That is the mechanism, and it is measured here rather than assumed.
    """
    dry, moist = _vapour(PROFILE), _vapour(PROFILE + 4.0)
    assert moist.sum() > dry.sum()

    open_window = _transmittance(dry, 10.0e-6)
    closed_window = _transmittance(moist, 10.0e-6)
    assert 0.5 < open_window < 0.8  # a window, not a hole and not a wall
    assert closed_window < open_window

    # Clausius-Clapeyron, which nothing here fitted: about 7 percent more vapour
    # per kelvin.
    per_kelvin = (moist.sum() / dry.sum()) ** 0.25 - 1.0
    assert 0.06 < per_kelvin < 0.09


@pytest.mark.req("EQR-002")
def test_optical_depth_is_linear_in_the_vapour_the_water_module_supplies():
    single = vapour_absorption_optical_depth(GRID, _vapour(PROFILE))
    double = vapour_absorption_optical_depth(GRID, 2.0 * _vapour(PROFILE))
    assert np.allclose(double, 2.0 * single, rtol=1e-12)


@pytest.mark.req("SYS-006", "MOD-002")
def test_unphysical_vapour_parameters_are_refused():
    with pytest.raises(ConfigurationError, match="cannot be negative"):
        vapour_absorption_optical_depth(GRID, -np.ones(N_LAYERS))
    with pytest.raises(ConfigurationError, match="non-negative"):
        vapour_mass_absorption(GRID, band_absorption=-1.0)
    with pytest.raises(ConfigurationError, match="non-negative"):
        vapour_mass_absorption(GRID, window_continuum=-1.0)
    with pytest.raises(ConfigurationError, match="wing decay"):
        vapour_mass_absorption(GRID, wing_decay_per_metre=0.0)


# --- what it was calibrated against -----------------------------------------


@pytest.mark.req("EQR-004", "INF-003", "VAL-002")
def test_the_clear_sky_outgoing_longwave_matches_the_calibration_target():
    """The calibration target is a radiation observation, not a temperature.

    CERES measures a clear-sky outgoing longwave flux near 265 W m^-2. Fitting
    to that leaves the surface-temperature record untouched and therefore still
    available as independent validation, which is the whole of INF-003.
    """
    assert _olr(0.0, follow=True) == pytest.approx(265.0, abs=2.0)

    # The clear-sky greenhouse effect follows from it rather than being a second
    # target: sigma T^4 at 288 K is 390 W m^-2, and the observed value is ~125.
    assert 390.1 - _olr(0.0, follow=True) == pytest.approx(125.0, abs=5.0)


@pytest.mark.req("EQR-004", "SYS-006")
def test_the_co2_band_has_a_calibration_for_each_background_it_runs_against():
    """A band's fitted width depends on what else absorbs around it.

    With vapour active the 15 micron band is partly masked, so the same
    published forcing per doubling needs a narrower band than it does against
    the fixed background. Two calibrations coexist and neither is correct for
    the other configuration, so the selection is a function rather than a
    default somebody has to remember.
    """
    assert co2_wing_decay_for(active_vapour=False) == CO2_WING_DECAY_PER_METRE
    assert co2_wing_decay_for(active_vapour=True) == CO2_WING_DECAY_WITH_VAPOUR_PER_METRE
    assert CO2_WING_DECAY_WITH_VAPOUR_PER_METRE < CO2_WING_DECAY_PER_METRE

    def net(ppm: float) -> float:
        solution = _solve(0.0, follow=True, ppm=ppm)
        return float(np.sum(solution.downward[-1] - solution.upward[-1]))

    assert net(569.4) - net(284.7) == pytest.approx(3.70, abs=0.05)


# --- the feedback ------------------------------------------------------------


@pytest.mark.req("CLM-003", "VAL-002", "EQR-002")
def test_letting_the_vapour_follow_the_temperature_amplifies_the_response():
    """The result the whole mechanism exists to produce, and its null control.

    Nothing in the absorption relation knows about temperature. Holding the
    vapour path fixed while the column warms gives the response without the
    feedback; letting it follow constant relative humidity gives the response
    with it. The difference is the water-vapour feedback, and it is a
    consequence of two independent pieces of physics meeting -- Clausius-
    Clapeyron and a window that closes -- rather than a term anyone added.
    """
    held = _olr(0.5, follow=False) - _olr(-0.5, follow=False)
    free = _olr(0.5, follow=True) - _olr(-0.5, follow=True)

    assert held > 0.0  # a warmer column radiates more: the Planck response
    assert free > 0.0  # still stable -- the feedback does not run away
    assert free < held  # and it is amplifying

    feedback = held - free
    assert feedback == pytest.approx(2.3, abs=0.4)
    # Published clear-sky estimates are near 1.8 W m^-2 K^-1. This is larger
    # because the configuration has no lapse-rate feedback, which in the real
    # atmosphere cancels a substantial part of it -- declared, not tuned away.
    assert 1.5 < feedback < 3.0


@pytest.mark.req("CLM-003", "VAL-005", "WAT-002")
def test_vapour_confined_to_the_lowest_layer_has_no_feedback_at_all():
    """Why the coupled column does not get the feedback the relation can produce.

    The water module has no vertical transport -- it says so itself -- so
    evaporation enters the lowest layer and nothing carries it higher. In
    equilibrium the entire column's vapour sits in that layer at saturation.

    That is fatal, and not because the amount is small. A greenhouse absorber
    warms the surface only by emitting from somewhere colder than the surface.
    Vapour held a few kelvin below the surface re-emits at almost the
    temperature it absorbed at, so it traps nothing however much of it there is.
    The same column mass gives half the greenhouse effect and an amplification
    of exactly one.

    Pinned here rather than left in a run, because it is the difference between
    a feedback that is weak and a feedback that is structurally impossible, and
    only the second is worth changing the water module for.
    """
    distributed = _vapour(PROFILE)
    confined = np.zeros(N_LAYERS)
    confined[0] = distributed.sum()
    assert confined.sum() == pytest.approx(distributed.sum())  # the same water

    def amplification(build) -> float:
        held = _olr_with(build(PROFILE), 0.5) - _olr_with(build(PROFILE), -0.5)
        free = _olr_with(build(PROFILE + 0.5), 0.5) - _olr_with(build(PROFILE - 0.5), -0.5)
        return held / free

    spread = 390.1 - _olr_with(distributed, 0.0)
    stacked = 390.1 - _olr_with(confined, 0.0)
    assert spread > 1.9 * stacked  # roughly halved by being stacked at the bottom

    def confine(temperature: np.ndarray) -> np.ndarray:
        out = np.zeros(N_LAYERS)
        out[0] = _vapour(temperature).sum()
        return out

    assert amplification(_vapour) > 2.0  # the feedback the relation can produce
    assert amplification(confine) == pytest.approx(1.0, abs=1e-3)  # and none of it survives
