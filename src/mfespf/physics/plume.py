"""The zero-buoyancy bulk-plume closure for humidity and stratification.

Charter section 19 stage 3 and section 20; requirements WAT-002, EQR-002,
EQR-004, INF-003, CLM-003.

This exists because of a measured failure. Convective moisture transport keyed
to the dry adjustment mixes completely, which drives the top of the mixed group
past saturation, and condensation then clips it back; the water that survives
aloft is decided by which layer the group's top happens to fall in, so the
vapour feedback came out a property of the layer grid. Above the group the
column was bone dry, because evaporation is the only source of water in this
framework and convection its only transport.

The literature's answer to that is entrainment. A real plume dilutes as it
rises and detrains air that is not saturated, which is why the observed
mid-troposphere sits near a third to a half of saturation rather than at it.
The cheapest form of the idea is a bulk plume in the limit of zero cloud
buoyancy, and it has a property worth the whole implementation: **one
dimensionless parameter predicts both the humidity and the lapse rate**. The
relation is therefore over-determined -- calibrate against one and the other is
a prediction that can fail, which is the difference between a closure and a
knob.

The equations
-------------

Zero buoyancy means the plume is neutrally buoyant with respect to its
environment at every level, so plume and environment share a temperature and
the plume's moist static energy declines only through entrainment of
subsaturated air::

    dh*/dz = -epsilon (1 - RH) L q*                       (Singh & O'Gorman)

Closing the environment's water budget -- subsidence drying against detrainment
of saturated air and evaporating condensate -- gives the relative humidity in
terms of the plume's own rates::

    RH = (delta + alpha gamma - alpha epsilon)
         / (delta + gamma - alpha epsilon)                (Romps 2014)

with ``gamma = -d ln q* / dz`` the water-vapour lapse rate, which follows from
Clausius-Clapeyron and hydrostatic balance alone::

    gamma = L Gamma / (R_v T^2) - g / (R T)

``alpha`` is the fraction of condensate that evaporates rather than
precipitating, so the precipitation efficiency is ``1 - alpha``.

Holding the mass flux constant with height forces ``epsilon = delta`` and
collapses both to one dimensionless number ``a = epsilon (1 - alpha) / gamma``::

    RH = (a + alpha) / (1 + a)                            (Romps 2016)
    Gamma = (g/c_p) [1 + a + q* L/(R T)]
                   / [1 + a + q* L^2/(c_p R_v T^2)]

and that lapse rate reduces to the moist adiabat exactly as ``a -> 0``, which
is checked against this framework's own ``moist_adiabatic_lapse_rate`` rather
than against a second copy of the formula.

What is fitted here
-------------------

Nothing. ``a = 0.2`` is the published value for Earth's tropics and ``alpha``
follows from measured precipitation efficiencies, so the humidity this produces
is a prediction and not a knob.

It was not always so, and the way that went wrong is worth keeping. ``alpha``
was first *assumed* -- 0.3, asserted over 0.1 to 0.5, on the reasoning that a
declared range was more honest than a fitted number. The reasoning was fine and
the range was wrong: it placed the precipitation efficiency at 0.5 to 0.9 where
the literature places it at 0.2 to 0.4, so it excluded the measured value
instead of bracketing it. A declared range is only as honest as its bounds, and
asserting one is not a substitute for looking the number up.

What makes the corrected value trustworthy is that it is over-determined too.
``a = epsilon PE / gamma`` ties the bulk-plume parameter to the entrainment rate
and the precipitation efficiency, so the measured ``alpha`` can be checked
against entrainment rates it was not derived from. It passes and the assumed one
failed, which is how the error was found.

Sources
-------

Singh, M. S. and O'Gorman, P. A. (2013), "Influence of entrainment on the
thermal stratification in simulations of radiative-convective equilibrium",
Geophysical Research Letters 40, 4398-4403, doi:10.1002/grl.50796.

Romps, D. M. (2014), "An analytical model for tropical relative humidity",
Journal of Climate 27(19), 7432-7449, doi:10.1175/JCLI-D-14-00255.1.

Romps, D. M. (2016), "Clausius-Clapeyron scaling of CAPE from analytical
solutions to RCE", Journal of the Atmospheric Sciences 73(9), 3719-3737,
doi:10.1175/JAS-D-15-0327.1.

Seeley, J. T. and Wordsworth, R. D. (2023), "Moist convection is most vigorous
at intermediate atmospheric humidity", arXiv:2301.03669 -- the equations above
were read from this reproduction of them, the primary sources being paywalled.
"""

from __future__ import annotations

import numpy as np

from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError
from mfespf.physics.thermodynamics import saturation_specific_humidity

__all__ = [
    "BULK_PLUME_PARAMETER",
    "EVAPORATING_CONDENSATE_FRACTION",
    "EVAPORATING_FRACTION_RANGE",
    "bulk_plume_relative_humidity",
    "convective_mass_flux",
    "detrainment_rate",
    "entraining_lapse_rate",
    "equilibrium_relative_humidity",
    "fractional_entrainment_rate",
    "water_vapour_lapse_rate",
]

_CP_DRY = value_of("specific_heat_dry_air")
_R_DRY = value_of("gas_constant_dry_air")
_R_VAPOUR = value_of("gas_constant_water_vapour")
_L_VAPORIZATION = value_of("water_latent_heat_vaporization_273")
_GRAVITY = value_of("earth_standard_gravity")

#: The dry adiabatic lapse rate, and the stability floor below which no steady
#: subsidence balances a given radiative cooling.
_DRY_LAPSE_RATE = _GRAVITY / _CP_DRY
_MINIMUM_STABILITY = 1.0e-4

#: Half-depth, in metres, over which the mass flux's logarithmic gradient is
#: taken. A physical scale rather than a numerical one: detrainment belongs to a
#: cloud layer, and Romps puts the moistening lengthscale near two kilometres.
#: Fixing it in metres is what makes the closure grid-convergent.
_DETRAINMENT_LENGTHSCALE = 2000.0

#: The dimensionless bulk-plume parameter ``a = epsilon (1 - alpha) / gamma``.
#: Romps (2016) for Earth's tropics; Seeley and Wordsworth (2023) state that
#: this value "corresponds to typical values of precipitation efficiency and
#: entrainment rate as diagnosed in cloud-resolving simulations", and that a
#: generally falls in O(0.1) to O(1). Not fitted here, and not to be fitted
#: against a temperature record under any circumstances (INF-003).
BULK_PLUME_PARAMETER = 0.2

#: Fraction of condensate that evaporates instead of precipitating, so the
#: precipitation efficiency is ``1 - alpha``. Taken from measurement rather than
#: assumed: Langhans, Yeo and Romps (2015) find PE of 33 and 38 percent in
#: large-eddy simulations of cumulus congestus, using the definition this
#: closure needs -- the fraction of condensate that reaches the surface -- and
#: report Braham (1952) at 19 percent and other numerical studies near 20.
#: A precipitation efficiency near a third gives ``alpha`` near two thirds.
#:
#: This replaced an assumed 0.3 asserted over 0.1 to 0.5. That range was not
#: conservative, it was inverted: it put the precipitation efficiency at 0.5 to
#: 0.9 where the literature puts it at 0.2 to 0.4, and so excluded the measured
#: value entirely. Two checks now agree that this one is right, which the old
#: one failed -- see ``fractional_entrainment_rate``.
EVAPORATING_CONDENSATE_FRACTION = 0.67

#: The range the measurements span, for sensitivity studies. Wide, because
#: precipitation efficiency genuinely is: definitions differ, and estimates
#: differ within a definition.
EVAPORATING_FRACTION_RANGE = (0.62, 0.81)


def water_vapour_lapse_rate(temperature: np.ndarray, lapse_rate: np.ndarray | float) -> np.ndarray:
    """``gamma = L Gamma / (R_v T^2) - g / (R T)`` in m^-1.

    The fractional rate at which saturation specific humidity falls off with
    height, and the only place Clausius-Clapeyron enters the plume closure. It
    is a derived quantity, not a parameter: given the temperature and the lapse
    rate the column already has, gamma follows from the saturation vapour
    pressure and hydrostatic balance with nothing left to choose.

    Its reciprocal is the water-vapour scale height, a couple of kilometres for
    a terrestrial troposphere, which is the sanity check worth remembering.
    """
    T = np.asarray(temperature, dtype=float)
    if np.any(T <= 0.0):
        raise ConfigurationError("temperatures must be positive")
    gamma_env = np.asarray(lapse_rate, dtype=float)
    return (_L_VAPORIZATION * gamma_env) / (_R_VAPOUR * T**2) - _GRAVITY / (_R_DRY * T)


def bulk_plume_relative_humidity(
    *,
    entrainment: np.ndarray | float,
    detrainment: np.ndarray | float,
    vapour_lapse_rate: np.ndarray | float,
    evaporating_fraction: float = EVAPORATING_CONDENSATE_FRACTION,
) -> np.ndarray:
    """``RH = (delta + alpha gamma - alpha epsilon) / (delta + gamma - alpha epsilon)``.

    Romps (2014), from the environment's water budget: subsidence dries the
    environment, detrainment of saturated plume air and evaporating condensate
    moisten it, and the balance fixes the relative humidity. Nothing about
    temperature appears -- the temperature dependence enters only through
    ``gamma``, which is why the humidity this predicts is not a restatement of
    the warming it will go on to produce.
    """
    if not 0.0 <= evaporating_fraction < 1.0:
        raise ConfigurationError(
            f"the evaporating condensate fraction {evaporating_fraction} must be a "
            f"fraction below one; at one, no condensate precipitates at all"
        )
    eps = np.asarray(entrainment, dtype=float)
    delta = np.asarray(detrainment, dtype=float)
    gamma = np.asarray(vapour_lapse_rate, dtype=float)
    if np.any(eps < 0.0) or np.any(delta < 0.0):
        raise ConfigurationError("entrainment and detrainment rates cannot be negative")
    if np.any(gamma <= 0.0):
        raise ConfigurationError(
            "the water-vapour lapse rate must be positive; a column whose saturation "
            "humidity increases with height is outside this closure's domain"
        )
    denominator = delta + gamma - evaporating_fraction * eps
    if np.any(denominator <= 0.0):
        raise ConfigurationError(
            "the bulk-plume water budget has no positive solution for these rates"
        )
    numerator = delta + evaporating_fraction * gamma - evaporating_fraction * eps
    return np.asarray(numerator / denominator, dtype=float)


def equilibrium_relative_humidity(
    bulk_plume_parameter: float = BULK_PLUME_PARAMETER,
    evaporating_fraction: float = EVAPORATING_CONDENSATE_FRACTION,
) -> float:
    """``RH = (a + alpha) / (1 + a)`` -- the constant-mass-flux case.

    Romps (2016) holds the mass flux constant with height, which forces
    ``epsilon = delta`` and makes the relative humidity a single number. The
    two limits are worth stating because they are the whole physical content:
    undiluted convection (``a = 0``) leaves the environment at the humidity
    that evaporating condensate alone can supply, and infinitely dilute
    convection saturates it, because everything the plume carries is detrained
    into the environment before it can precipitate.
    """
    if bulk_plume_parameter < 0.0:
        raise ConfigurationError(
            f"the bulk-plume parameter {bulk_plume_parameter} cannot be negative; it is "
            f"proportional to the entrainment rate"
        )
    if not 0.0 <= evaporating_fraction < 1.0:
        raise ConfigurationError(
            f"the evaporating condensate fraction {evaporating_fraction} must be a "
            f"fraction below one"
        )
    return (bulk_plume_parameter + evaporating_fraction) / (1.0 + bulk_plume_parameter)


def fractional_entrainment_rate(
    vapour_lapse_rate: np.ndarray | float,
    bulk_plume_parameter: float = BULK_PLUME_PARAMETER,
    evaporating_fraction: float = EVAPORATING_CONDENSATE_FRACTION,
) -> np.ndarray:
    """``epsilon = a gamma / (1 - alpha)`` in m^-1.

    The dimensional entrainment rate the dimensionless parameter stands for.
    Reported rather than used internally, because it is the number the
    cloud-resolving literature quotes and therefore the one a reader can check
    this closure against: a few tenths per kilometre for deep convection.

    It is also the check that caught the evaporating fraction being wrong. The
    three quantities are not independent -- ``a = epsilon PE / gamma`` -- so a
    parameter pair can be tested against entrainment rates even when it was
    chosen from humidity considerations. At the measured precipitation
    efficiency this gives 0.20 to 0.34 per kilometre, inside the 0.2 to 0.8 that
    cloud-resolving simulations diagnose; at the value previously assumed it
    gave 0.09, outside it. Running the identity the other way, an entrainment
    rate of 0.2 per kilometre with a precipitation efficiency of 0.33 returns
    ``a = 0.204``, against the 0.2 Romps quotes, which closes the loop.
    """
    if not 0.0 <= evaporating_fraction < 1.0:
        raise ConfigurationError(
            f"the evaporating condensate fraction {evaporating_fraction} must be a "
            f"fraction below one"
        )
    gamma = np.asarray(vapour_lapse_rate, dtype=float)
    return bulk_plume_parameter * gamma / (1.0 - evaporating_fraction)


def entraining_lapse_rate(
    temperature: np.ndarray,
    pressure: np.ndarray,
    bulk_plume_parameter: float = BULK_PLUME_PARAMETER,
) -> np.ndarray:
    """The lapse rate of an atmosphere neutral to an entraining plume, in K m^-1.

    ``Gamma = (g/c_p) [1 + a + q* L/(R T)] / [1 + a + q* L^2/(c_p R_v T^2)]``

    This is the second thing the one parameter predicts, and the reason the
    closure can be falsified rather than tuned: fix ``a`` against the humidity
    and the stratification is no longer free. Entrainment of subsaturated air
    robs the plume of buoyancy, so the environment it leaves neutral is more
    stable than a moist adiabat and less stable than a dry one, and the
    difference is what shows up as convective available potential energy.

    At ``a = 0`` it reduces exactly to ``moist_adiabatic_lapse_rate``, which is
    the limiting case the tests check against this framework's existing
    implementation rather than against a second copy of the same algebra.
    """
    if bulk_plume_parameter < 0.0:
        raise ConfigurationError(
            f"the bulk-plume parameter {bulk_plume_parameter} cannot be negative"
        )
    T = np.asarray(temperature, dtype=float)
    q_s = saturation_specific_humidity(T, pressure)
    a = float(bulk_plume_parameter)
    numerator = 1.0 + a + (_L_VAPORIZATION * q_s) / (_R_DRY * T)
    denominator = 1.0 + a + (_L_VAPORIZATION**2 * q_s) / (_CP_DRY * _R_VAPOUR * T**2)
    return (_GRAVITY / _CP_DRY) * numerator / denominator


def convective_mass_flux(
    radiative_cooling: np.ndarray,
    temperature: np.ndarray,
    pressure: np.ndarray,
    lapse_rate: np.ndarray,
) -> np.ndarray:
    """The convective mass flux implied by clear-sky radiative cooling, kg m^-2 s^-1.

    In radiative-convective equilibrium the environment is not free to choose
    how fast it subsides. It sinks at whatever rate makes adiabatic compression
    balance its radiative cooling::

        -w (Gamma_d - Gamma) + Q = 0

    and mass continuity then requires the convective updraught to carry the same
    mass back up, so ``M = rho |Q| / (Gamma_d - Gamma)``.

    This is the constraint that turns the bulk-plume closure from a
    one-parameter relation into a profile. Romps' derivation notes that the mass
    flux -- and therefore ``d_z M = e - d`` -- is "largely constrained by
    clear-sky radiative cooling"; this is that constraint, evaluated rather than
    assumed. Nothing here is fitted: the cooling comes from the radiative
    transfer the column already solves.

    ``radiative_cooling`` is a positive-definite cooling rate in K s^-1. Where
    the column is being heated rather than cooled, or where the stratification
    approaches the dry adiabat, no steady subsidence balances it and the mass
    flux is zero -- which is the top of the convecting layer.
    """
    cooling = np.asarray(radiative_cooling, dtype=float)
    T = np.asarray(temperature, dtype=float)
    p = np.asarray(pressure, dtype=float)
    stability = _DRY_LAPSE_RATE - np.asarray(lapse_rate, dtype=float)

    density = p / (_R_DRY * T)
    flux = np.zeros_like(cooling)
    active = (cooling > 0.0) & (stability > _MINIMUM_STABILITY)
    flux[active] = density[active] * cooling[active] / stability[active]
    return flux


def detrainment_rate(
    mass_flux: np.ndarray,
    entrainment: np.ndarray,
    height: np.ndarray,
) -> np.ndarray:
    """``delta = epsilon - d_z ln M`` in m^-1, from ``d_z M = (epsilon - delta) M``.

    The bulk-plume mass budget gives the *difference* between the rates, not
    each of them, so one more relation is needed and the entrainment rate
    supplies it. What comes out is the reason the observed humidity profile is
    C-shaped: where the mass flux is falling off with height -- the upper
    troposphere, as convection runs out of buoyancy -- detrainment exceeds
    entrainment, the environment is moistened by saturated air, and the relative
    humidity climbs back towards one.

    Detrainment cannot be negative: a plume that is not detraining is simply not
    detraining. Where the mass budget would ask for less than zero it is clipped,
    and that clipping is the closure admitting the profile is entraining faster
    than the mass flux supports.
    """
    M = np.asarray(mass_flux, dtype=float)
    eps = np.broadcast_to(np.asarray(entrainment, dtype=float), M.shape)
    z = np.asarray(height, dtype=float)

    log_gradient = np.zeros_like(M)
    active = M > 0.0
    if np.count_nonzero(active) >= 2:
        # Differentiated in the log, so the result is the fractional rate the
        # bulk-plume equations are written in rather than an absolute gradient
        # -- and over a fixed physical depth rather than layer by layer.
        #
        # That distinction decides whether the closure converges. Detrainment is
        # a property of a cloud layer, whose depth is set by the physics at
        # around a kilometre or two (Romps puts the moistening lengthscale 1/δ
        # near 2 km), not by whatever thickness the grid happens to use. Taken
        # layer by layer, the derivative differentiates the noise in a mass flux
        # that is itself derived from the radiative cooling the humidity
        # controls, so a kink feeds back on itself: refining the grid then makes
        # the answer worse rather than better, with more and more layers driven
        # onto the zero-detrainment floor. Measured, the surface temperature ran
        # away from 299.9 K at ten layers to 305.9 K at sixty.
        za, la = z[active], np.log(M[active])
        slopes = np.zeros_like(la)
        for i, centre in enumerate(za):
            window = np.abs(za - centre) <= _DETRAINMENT_LENGTHSCALE
            if np.count_nonzero(window) >= 2:
                slopes[i] = np.polyfit(za[window], la[window], 1)[0]
        log_gradient[active] = slopes

    return np.maximum(eps - log_gradient, 0.0)
