"""Dry and moist thermodynamic primitives.

Charter section 19 stage 3; requirements ATM-002, WAT-001, WAT-002, VAL-001.

Functions rather than a module, for the same reason as the radiative
primitives: several modules evaluate them, and one implementation means the
convective and the moist paths cannot disagree about what a dry adiabat is.

The important routine is ``dry_convective_adjustment``.  It conserves the
column's enthalpy exactly, by construction rather than to a tolerance: the mixed
potential temperature of each unstable group is *defined* as the value that
leaves ``sum(C_i T_i)`` unchanged.  That matters because convection is the first
process in this framework that moves energy between layers without any transfer
crossing a module boundary, and its net contribution to the column budget must
therefore be exactly zero.
"""

from __future__ import annotations

import numpy as np

from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError

__all__ = [
    "convective_adjustment_with_groups",
    "dry_adiabatic_lapse_rate",
    "dry_convective_adjustment",
    "exner",
    "moist_adiabatic_lapse_rate",
    "potential_temperature",
    "saturation_specific_humidity",
    "saturation_vapour_pressure",
    "surface_coupled_adjustment",
]

_CP_DRY = value_of("specific_heat_dry_air")
_R_DRY = value_of("gas_constant_dry_air")
_R_VAPOUR = value_of("gas_constant_water_vapour")
_L_VAPORIZATION = value_of("water_latent_heat_vaporization_273")
_T_TRIPLE = value_of("triple_point_temperature")
_E_TRIPLE = value_of("saturation_vapour_pressure_triple_point")
_GRAVITY = value_of("earth_standard_gravity")

#: Ratio of the gas constants, R_d / R_v. Appears wherever vapour pressure is
#: converted to a mass mixing ratio.
EPSILON = _R_DRY / _R_VAPOUR

#: R_d / c_p, the exponent in the Exner function.
KAPPA = _R_DRY / _CP_DRY


def exner(pressure: np.ndarray, reference_pressure: float | None = None) -> np.ndarray:
    """``(p / p_ref)^kappa`` -- the factor relating temperature to potential temperature."""
    p_ref = (
        reference_pressure
        if reference_pressure is not None
        else value_of("standard_surface_pressure")
    )
    if p_ref <= 0.0:
        raise ConfigurationError("the reference pressure must be positive")
    p = np.asarray(pressure, dtype=float)
    if np.any(p <= 0.0):
        raise ConfigurationError("pressures must be positive to form the Exner function")
    return (p / p_ref) ** KAPPA


def potential_temperature(
    temperature: np.ndarray, pressure: np.ndarray, reference_pressure: float | None = None
) -> np.ndarray:
    """``theta = T / exner(p)``.

    Conserved under dry adiabatic displacement, which is what makes it the right
    variable for a stability criterion: an atmosphere is stable against dry
    convection exactly when ``theta`` does not decrease with height.
    """
    return np.asarray(temperature, dtype=float) / exner(pressure, reference_pressure)


def dry_adiabatic_lapse_rate() -> float:
    """``g / c_p`` in K m^-1, about 9.8 K per kilometre."""
    return _GRAVITY / _CP_DRY


def _merged(groups: list[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    """Collapse overlapping or touching mixed regions into maximal ones."""
    if not groups:
        return ()
    ordered = sorted(groups)
    out = [ordered[0]]
    for low, high in ordered[1:]:
        last_low, last_high = out[-1]
        if low <= last_high + 1:
            out[-1] = (last_low, max(last_high, high))
        else:
            out.append((low, high))
    return tuple(out)


def dry_convective_adjustment(
    temperature: np.ndarray,
    heat_capacity: np.ndarray,
    layer_pressure: np.ndarray,
    *,
    reference_pressure: float | None = None,
    tolerance: float = 1e-12,
    max_passes: int = 1000,
) -> np.ndarray:
    """Mix convectively unstable layers to a dry adiabat, conserving enthalpy.

    Layer index 0 is the lowest layer.  Where potential temperature decreases
    with height the column is unstable; those layers are mixed to a single
    potential temperature chosen so that ``sum(C_i T_i)`` is unchanged, and the
    group is grown until it is stable against both neighbours.

    This is the Manabe-Wetherald convective adjustment. It is a *closure*, not a
    derivation: real convection has finite speed, entrains, and does not
    instantaneously homogenise. What it buys is a profile that is stable against
    dry displacement, reached without creating or destroying energy.

    Returns a new temperature array; the input is not modified.
    """
    return convective_adjustment_with_groups(
        temperature,
        heat_capacity,
        layer_pressure,
        reference_pressure=reference_pressure,
        tolerance=tolerance,
        max_passes=max_passes,
    )[0]


def convective_adjustment_with_groups(
    temperature: np.ndarray,
    heat_capacity: np.ndarray,
    layer_pressure: np.ndarray,
    *,
    reference_pressure: float | None = None,
    tolerance: float = 1e-12,
    max_passes: int = 1000,
) -> tuple[np.ndarray, tuple[tuple[int, int], ...]]:
    """The same adjustment, also reporting which layers it mixed together.

    The groups are what convection *did*: each ``(low, high)`` pair is an
    inclusive run of layers the closure homogenised. They are worth returning
    because heat is not the only thing convection carries. A moisture transport
    that mixes over exactly these layers reuses the closure rather than
    inventing a second one, and so introduces no parameter the temperature
    adjustment did not already have -- which is the difference between deriving
    a humidity profile and fitting one.

    Overlapping groups from successive passes are merged, so what comes back is
    the final set of mixed regions rather than a history of how it got there.
    """
    T = np.array(temperature, dtype=float)
    C = np.asarray(heat_capacity, dtype=float)
    n = T.size
    if C.size != n or np.asarray(layer_pressure).size != n:
        raise ConfigurationError(
            "temperature, heat capacity and pressure must describe the same layers"
        )
    if np.any(C <= 0.0):
        raise ConfigurationError("layer heat capacities must be positive")
    if n < 2:
        return T, ()

    pi = exner(layer_pressure, reference_pressure)
    theta = T / pi

    groups: list[tuple[int, int]] = []
    for _ in range(max_passes):
        unstable = np.nonzero(theta[1:] < theta[:-1] - tolerance)[0]
        if unstable.size == 0:
            return theta * pi, _merged(groups)
        low = int(unstable[0])
        high = low + 1
        while True:
            # The mixed value is *defined* by enthalpy conservation, so the
            # adjustment cannot create or destroy energy however the group grows.
            weights = C[low : high + 1] * pi[low : high + 1]
            mixed = float(np.sum(C[low : high + 1] * T[low : high + 1]) / np.sum(weights))
            if low > 0 and theta[low - 1] > mixed + tolerance:
                low -= 1
                continue
            if high < n - 1 and theta[high + 1] < mixed - tolerance:
                high += 1
                continue
            break
        theta[low : high + 1] = mixed
        T[low : high + 1] = mixed * pi[low : high + 1]
        groups.append((low, high))

    raise ConfigurationError(
        f"convective adjustment did not stabilise in {max_passes} passes; the profile or "
        f"the tolerance is outside the scheme's validity domain"
    )


def surface_coupled_adjustment(
    temperature: np.ndarray,
    heat_capacity: np.ndarray,
    layer_pressure: np.ndarray,
    *,
    surface_temperature: float,
    surface_heat_capacity: float,
    surface_pressure: float,
    reference_pressure: float | None = None,
    tolerance: float = 1e-12,
    max_passes: int = 1000,
) -> tuple[float, np.ndarray, tuple[tuple[int, int], ...], bool]:
    """The same adjustment with the surface as the column's bottom element.

    ``convective_adjustment_with_groups`` mixes air with air, so it can only act
    on an instability *between layers*. But the instability that drives real
    convection is not between layers, it is between the ground and the air on
    top of it: the surface is heated, the air in contact with it becomes
    buoyant, and the plume that rises carries heat and moisture through whatever
    stratification happens to lie above. A scheme that cannot see the surface
    cannot see the cause.

    So the surface is prepended as element zero of the profile, with its own
    heat capacity and its own Exner factor at surface pressure, and the existing
    algorithm runs unchanged. Nothing is added: the mixed value is still defined
    by enthalpy conservation, now over the surface and the air together, and no
    parameter enters that the layer-only adjustment did not already have. This
    is what makes it worth trying before an eddy diffusivity -- the diffusivity
    would have to be fitted, and this does not.

    It is the parameter-free counterpart of the bulk sensible heat flux, and the
    two differ in exactly the way that matters here. The flux moves energy at a
    rate set by a fitted coefficient, so a layer that radiates hard enough can
    stay far colder than the surface for ever. The adjustment asserts that where
    the ground is buoyant relative to the air, mixing has already removed the
    instability, at whatever rate that took.

    Returns the adjusted surface temperature, the adjusted layer temperatures,
    the mixed groups in *layer* indices, and whether the surface took part in
    any of them.
    """
    T = np.asarray(temperature, dtype=float)
    C = np.asarray(heat_capacity, dtype=float)
    p = np.asarray(layer_pressure, dtype=float)
    if surface_heat_capacity <= 0.0:
        raise ConfigurationError("the surface heat capacity must be positive")
    if surface_pressure <= 0.0:
        raise ConfigurationError("the surface pressure must be positive")
    if T.size and surface_pressure < float(np.max(p)):
        raise ConfigurationError(
            f"the surface pressure {surface_pressure} lies above the lowest layer "
            f"at {float(np.max(p))}; the surface must be the bottom of the column"
        )

    extended, groups = convective_adjustment_with_groups(
        np.concatenate(([float(surface_temperature)], T)),
        np.concatenate(([float(surface_heat_capacity)], C)),
        np.concatenate(([float(surface_pressure)], p)),
        reference_pressure=reference_pressure,
        tolerance=tolerance,
        max_passes=max_passes,
    )

    # Back to layer indices. A group that reaches the surface starts at layer 0
    # by construction, since the surface is element zero and groups are runs.
    surface_mixed = any(low == 0 for low, _ in groups)
    layer_groups = tuple((max(low - 1, 0), high - 1) for low, high in groups if high >= 1)
    return float(extended[0]), extended[1:], layer_groups, surface_mixed


def saturation_vapour_pressure(temperature: np.ndarray) -> np.ndarray:
    """Clausius-Clapeyron saturation vapour pressure over liquid water, in Pa.

    Integrated with a constant latent heat, which is the standard reduced form:

        e_s(T) = e_triple * exp( (L / R_v) * (1/T_triple - 1/T) )

    Holding L constant is an approximation worth about a percent between 250 and
    300 K and degrading outside it; the equation registry carries that as the
    declared validity domain rather than leaving it to the caller to remember.
    """
    T = np.asarray(temperature, dtype=float)
    if np.any(T <= 0.0):
        raise ConfigurationError("temperatures must be positive")
    return _E_TRIPLE * np.exp((_L_VAPORIZATION / _R_VAPOUR) * (1.0 / _T_TRIPLE - 1.0 / T))


def saturation_specific_humidity(temperature: np.ndarray, pressure: np.ndarray) -> np.ndarray:
    """Saturation specific humidity in kg of vapour per kg of moist air.

    ``q_s = eps e_s / (p - (1 - eps) e_s)``. The full form rather than the
    ``eps e_s / p`` approximation: near the surface at 300 K the two differ by
    about 2 percent, which is larger than several of the effects this framework
    is meant to be able to resolve.
    """
    e_s = saturation_vapour_pressure(temperature)
    p = np.asarray(pressure, dtype=float)
    if np.any(p <= 0.0):
        raise ConfigurationError("pressures must be positive")
    denominator = p - (1.0 - EPSILON) * e_s
    # Above the boiling point the denominator vanishes: the air cannot be
    # saturated at that pressure at all. Clamp rather than divide by zero, and
    # let the validity domain of the equation record carry the caveat.
    denominator = np.maximum(denominator, 1e-6 * p)
    return EPSILON * e_s / denominator


def moist_adiabatic_lapse_rate(temperature: np.ndarray, pressure: np.ndarray) -> np.ndarray:
    """Saturated adiabatic lapse rate in K m^-1.

    Always smaller than the dry rate because condensation releases latent heat
    into a rising parcel; it approaches the dry rate as the saturation humidity
    goes to zero, which is the limiting case the tests check.
    """
    T = np.asarray(temperature, dtype=float)
    q_s = saturation_specific_humidity(T, pressure)
    numerator = 1.0 + (_L_VAPORIZATION * q_s) / (_R_DRY * T)
    denominator = 1.0 + (_L_VAPORIZATION**2 * q_s) / (_CP_DRY * _R_VAPOUR * T**2)
    return (_GRAVITY / _CP_DRY) * numerator / denominator
