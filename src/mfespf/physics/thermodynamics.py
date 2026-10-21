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
    "SATURATION_DOMAIN_FLOOR_K",
    "convective_adjustment_with_groups",
    "deposition_latent_heat_excess",
    "dry_adiabatic_lapse_rate",
    "dry_convective_adjustment",
    "evaporation_coefficient_for",
    "exner",
    "hydrostatic_heights",
    "lapse_rate_adjustment_with_groups",
    "lapse_rate_offsets",
    "lifting_condensation_pressure",
    "moist_adiabatic_lapse_rate",
    "potential_temperature",
    "saturation_humidity_within_domain",
    "saturation_specific_humidity",
    "saturation_vapour_pressure",
    "surface_coupled_adjustment",
]

#: Coldest temperature at which `EQ-CLAUSIUS-CLAPEYRON-SATURATION` is declared
#: valid, and therefore the coldest at which a humidity closure in this
#: framework will evaluate it. It lives here so the relation, its declared
#: domain and every caller share one number.
#:
#: The floor is measured rather than asserted. Against Murphy and Koop (2005)
#: the constant-heat-capacity form is within 2 percent to 200 K, 6.3 percent at
#: 172 K and 13.6 percent by 150 K, and 170 K is where the record puts the line.
SATURATION_DOMAIN_FLOOR_K = 170.0

_CP_DRY = value_of("specific_heat_dry_air")
_R_DRY = value_of("gas_constant_dry_air")
_R_VAPOUR = value_of("gas_constant_water_vapour")
_L_VAPORIZATION = value_of("water_latent_heat_vaporization_273")
_L_FUSION = value_of("water_latent_heat_fusion")
#: Latent heat of sublimation: a molecule leaving ice for vapour must be both
#: melted and evaporated, so the two latent heats add.
_L_SUBLIMATION = _L_VAPORIZATION + _L_FUSION
_CP_VAPOUR = value_of("specific_heat_water_vapour")
_CP_LIQUID = value_of("specific_heat_liquid_water")
_CP_ICE = value_of("specific_heat_ice")
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
    partial_tops: dict[int, float] | None = None,
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
    T = np.asarray(temperature, dtype=float)
    if np.asarray(layer_pressure).size != T.size:
        raise ConfigurationError(
            "temperature, heat capacity and pressure must describe the same layers"
        )
    pi = exner(layer_pressure, reference_pressure) if T.size else np.ones(0)
    return _adjust_to_neutral(
        T,
        heat_capacity,
        scale=pi,
        offset=np.zeros(T.size),
        tolerance=tolerance,
        max_passes=max_passes,
        partial_tops=partial_tops,
    )


def _adjust_to_neutral(
    temperature: np.ndarray,
    heat_capacity: np.ndarray,
    *,
    scale: np.ndarray,
    offset: np.ndarray,
    tolerance: float = 1e-12,
    max_passes: int = 1000,
    partial_tops: dict[int, float] | None = None,
) -> tuple[np.ndarray, tuple[tuple[int, int], ...]]:
    """Mix to uniform ``s = (T - offset) / scale``, conserving ``sum(C T)`` exactly.

    One algorithm for two targets, because they differ only in what "neutral"
    means. The dry adiabat is neutral when ``T / pi`` is uniform, so
    ``scale = pi`` and ``offset = 0``. Any prescribed lapse rate is neutral when
    ``T`` plus the cumulative drop that rate implies is uniform, so ``scale = 1``
    and ``offset`` is minus that drop. Writing it once means a second target
    cannot quietly acquire different conservation behaviour from the first.

    The mixed value stays *defined* by enthalpy conservation rather than
    computed and then checked:

        s = (sum(C T) - sum(C offset)) / sum(C scale)

    so the adjustment cannot create or destroy energy however the group grows,
    for either target.

    **A top that is not a layer edge.** When ``partial_tops`` is given, a group
    that stops growing penetrates the next layer by a mass fraction ``f`` rather
    than stopping at its lower edge. Inside that layer the stability variable is
    reconstructed as linear in mass with the gradient it has to the layer above,
    keeping the layer's mean; the mixed part is what lies below the level where
    that reconstruction equals the group's mixed value, i.e. the level of
    neutral buoyancy. ``f`` reaches 1 exactly where the whole layer would have
    joined, so the top moves continuously in pressure instead of snapping to an
    edge -- which is what made the column's error follow the grid. The fraction
    of each partially mixed layer is written into ``partial_tops``, keyed by
    layer index. Enthalpy stays exact: the reconstruction integrates to the
    layer mean, and the unmixed remainder keeps its reconstructed mean, which
    lies above the group's value, so no instability is created.
    """
    T = np.array(temperature, dtype=float)
    C = np.asarray(heat_capacity, dtype=float)
    n = T.size
    if C.size != n:
        raise ConfigurationError(
            "temperature, heat capacity and pressure must describe the same layers"
        )
    if np.any(C <= 0.0):
        raise ConfigurationError("layer heat capacities must be positive")
    if np.any(scale <= 0.0):
        raise ConfigurationError("the stability variable's scale must be positive")
    if n < 2:
        return T, ()

    stability = (T - offset) / scale
    fractional = partial_tops is not None
    partial: dict[int, float] = {}

    groups: list[tuple[int, int]] = []
    for _ in range(max_passes):
        unstable = np.nonzero(stability[1:] < stability[:-1] - tolerance)[0]
        if unstable.size == 0:
            final = _merged(groups)
            if partial_tops is not None:
                partial_tops.clear()
                partial_tops.update(
                    {
                        k: f
                        for k, f in partial.items()
                        if not any(low <= k <= high for low, high in final)
                    }
                )
            return stability * scale + offset, final
        low = int(unstable[0])
        # Start from the buoyant layer alone and let it rise: in whole-layer
        # mode the unstable neighbour joins on the first pass exactly as
        # before, and in fractional mode it joins only as far as the parcel
        # is buoyant in it, so the mixing depth grows from zero at onset.
        high = low
        f = 0.0
        while True:
            span = slice(low, high + 1)
            numerator = float(np.sum(C[span] * T[span]) - np.sum(C[span] * offset[span]))
            denominator = float(np.sum(C[span] * scale[span]))
            mixed = numerator / denominator
            f = 0.0
            if low > 0 and stability[low - 1] > mixed + tolerance:
                low -= 1
                continue
            if high < n - 1:
                k = high + 1
                gradient = stability[k + 1] - stability[k] if k + 1 < n else 0.0
                if fractional and gradient > 0.0:
                    f = _penetration(
                        numerator, denominator, C[k] * scale[k], stability[k], gradient
                    )
                    if f >= 1.0:
                        high += 1
                        continue
                elif stability[k] < mixed - tolerance:
                    high += 1
                    continue
            break
        if f > 0.0:
            k = high + 1
            gradient = stability[k + 1] - stability[k]
            bottom = stability[k] - 0.5 * gradient
            weight = C[k] * scale[k]
            part = bottom + 0.5 * gradient * f
            mixed = (numerator + f * weight * part) / (denominator + f * weight)
            # The unmixed remainder keeps its reconstructed mean, which lies
            # above the level of neutral buoyancy and so above ``mixed``.
            remainder = bottom + 0.5 * gradient * (1.0 + f)
            stability[k] = f * mixed + (1.0 - f) * remainder
            T[k] = stability[k] * scale[k] + offset[k]
            partial[k] = f
        stability[low : high + 1] = mixed
        T[low : high + 1] = mixed * scale[low : high + 1] + offset[low : high + 1]
        groups.append((low, high))

    raise ConfigurationError(
        f"convective adjustment did not stabilise in {max_passes} passes; the profile or "
        f"the tolerance is outside the scheme's validity domain"
    )


def _penetration(
    numerator: float, denominator: float, weight: float, mean: float, gradient: float
) -> float:
    """Mass fraction of the next layer a mixed group rises through, in [0, 1].

    The group's mixed value with a fraction ``f`` of the layer included is
    ``(numerator + f w m(f)) / (denominator + f w)``, where ``m(f)`` is the
    mean of the reconstructed profile over the lowest ``f`` of the layer. The
    group stops where that equals the reconstructed profile itself, at
    ``bottom + gradient f``. Monotone in ``f`` for a positive gradient, so
    bisection finds it; 0 when the group is not buoyant at the layer's bottom
    edge, 1 when it is buoyant at its top.
    """
    bottom = mean - 0.5 * gradient

    def excess(f: float) -> float:
        mixed = (numerator + f * weight * (bottom + 0.5 * gradient * f)) / (
            denominator + f * weight
        )
        return mixed - (bottom + gradient * f)

    if excess(0.0) <= 0.0:
        return 0.0
    if excess(1.0) >= 0.0:
        return 1.0
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if excess(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def hydrostatic_heights(
    temperature: np.ndarray, pressure: np.ndarray, surface_pressure: float
) -> np.ndarray:
    """Geometric height of each layer midpoint above the surface, in metres.

    The hypsometric equation, integrated upward from the surface through the
    column's own temperatures. It is a diagnostic and not a state: this
    framework's columns are defined in pressure, and height appears only where a
    relation is written in it -- a lapse rate in kelvin per metre being the case
    that forced the issue.
    """
    T = np.asarray(temperature, dtype=float)
    p = np.asarray(pressure, dtype=float)
    if T.size != p.size:
        raise ConfigurationError("temperature and pressure must describe the same layers")
    if surface_pressure <= 0.0 or np.any(p <= 0.0):
        raise ConfigurationError("pressures must be positive to integrate hydrostatically")
    if T.size == 0:
        return np.zeros(0)

    scale = _R_DRY / _GRAVITY
    # From the surface to the lowest midpoint, then midpoint to midpoint with the
    # mean temperature of each pair, which is the standard layer-by-layer form.
    steps = np.empty(T.size)
    steps[0] = scale * T[0] * np.log(surface_pressure / p[0])
    if T.size > 1:
        mean_T = 0.5 * (T[:-1] + T[1:])
        steps[1:] = scale * mean_T * np.log(p[:-1] / p[1:])
    return np.cumsum(steps)


def lapse_rate_offsets(lapse_rate: np.ndarray, height: np.ndarray) -> np.ndarray:
    """Cumulative temperature drop a target lapse rate implies, in K.

    ``offset[i]`` is how much colder than the bottom of the column the target
    profile is at layer ``i``. Subtracting it turns "the profile follows this
    lapse rate" into "this quantity is uniform", which is what lets one mixing
    routine serve any target.
    """
    G = np.asarray(lapse_rate, dtype=float)
    z = np.asarray(height, dtype=float)
    if G.size != z.size:
        raise ConfigurationError("the lapse rate and the height must describe the same layers")
    if G.size == 0:
        return np.zeros(0)
    drops = np.zeros(G.size)
    if G.size > 1:
        drops[1:] = 0.5 * (G[:-1] + G[1:]) * np.diff(z)
    return np.cumsum(drops)


def lapse_rate_adjustment_with_groups(
    temperature: np.ndarray,
    heat_capacity: np.ndarray,
    target_offsets: np.ndarray,
    *,
    tolerance: float = 1e-12,
    max_passes: int = 1000,
    partial_tops: dict[int, float] | None = None,
) -> tuple[np.ndarray, tuple[tuple[int, int], ...]]:
    """Adjust toward a prescribed lapse rate instead of toward the dry adiabat.

    Same closure, same conservation, different target. Where the profile is
    steeper than the target it is mixed until it is not, and the group's
    enthalpy is unchanged by construction. The dry adjustment is the special
    case whose target is ``g / c_p``, and it keeps its own routine because the
    Exner form is exact in pressure coordinates and this one is not.
    """
    T = np.asarray(temperature, dtype=float)
    offsets = np.asarray(target_offsets, dtype=float)
    if offsets.size != T.size:
        raise ConfigurationError(
            "the target offsets and the temperature must describe the same layers"
        )
    return _adjust_to_neutral(
        T,
        heat_capacity,
        scale=np.ones(T.size),
        offset=-offsets,
        tolerance=tolerance,
        max_passes=max_passes,
        partial_tops=partial_tops,
    )


def surface_coupled_adjustment(
    temperature: np.ndarray,
    heat_capacity: np.ndarray,
    layer_pressure: np.ndarray,
    *,
    surface_temperature: float,
    surface_heat_capacity: float,
    surface_pressure: float,
    target_offsets: np.ndarray | None = None,
    reference_pressure: float | None = None,
    tolerance: float = 1e-12,
    max_passes: int = 1000,
    partial_tops: dict[int, float] | None = None,
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

    all_T = np.concatenate(([float(surface_temperature)], T))
    all_C = np.concatenate(([float(surface_heat_capacity)], C))
    extended_partials: dict[int, float] | None = {} if partial_tops is not None else None
    if target_offsets is None:
        extended, groups = convective_adjustment_with_groups(
            all_T,
            all_C,
            np.concatenate(([float(surface_pressure)], p)),
            reference_pressure=reference_pressure,
            tolerance=tolerance,
            max_passes=max_passes,
            partial_tops=extended_partials,
        )
    else:
        # A prescribed lapse rate instead of the dry adiabat. The offsets are
        # measured from the surface, so element zero is zero by construction and
        # the surface enters on the same footing it does above.
        offsets = np.asarray(target_offsets, dtype=float)
        if offsets.size != all_T.size:
            raise ConfigurationError(
                "the target offsets must cover the surface and every layer, so one "
                "more value than there are layers"
            )
        extended, groups = lapse_rate_adjustment_with_groups(
            all_T,
            all_C,
            offsets,
            tolerance=tolerance,
            max_passes=max_passes,
            partial_tops=extended_partials,
        )

    # Back to layer indices. A group that reaches the surface starts at layer 0
    # by construction, since the surface is element zero and groups are runs.
    surface_mixed = any(low == 0 for low, _ in groups)
    layer_groups = tuple((max(low - 1, 0), high - 1) for low, high in groups if high >= 1)
    if partial_tops is not None and extended_partials is not None:
        # A partial layer always lies above a group, so never at the surface.
        partial_tops.clear()
        partial_tops.update({k - 1: f for k, f in extended_partials.items() if k >= 1})
    return float(extended[0]), extended[1:], layer_groups, surface_mixed


def saturation_vapour_pressure(temperature: np.ndarray) -> np.ndarray:
    """Clausius-Clapeyron saturation vapour pressure, in Pa, over the stable phase.

    Two corrections to the reduced form this replaced, both of which matter in a
    column whose top is at 172 K.

    **The latent heat is not constant.** Kirchhoff's relation makes it vary with
    temperature at a rate set by the heat capacities of the two phases::

        dL/dT = c_pv - c_condensed        so     L(T) = L_0 + dC (T - T_t)

    Integrating Clausius-Clapeyron ``d ln(e_s)/dT = L(T) / (R_v T^2)`` with that
    linear form is exact rather than numerical, and gives::

        e_s(T) = e_t (T/T_t)^(dC/R_v)
                 exp{ [L_0 - dC T_t] / R_v * (1/T_t - 1/T) }

    Setting ``dC = 0`` recovers the constant-latent-heat expression exactly,
    which is the limiting case the tests use.

    **Below the triple point the stable phase is ice.** Vapour is then in
    equilibrium with ice rather than with supercooled liquid, so the relevant
    latent heat is that of sublimation, ``L_v + L_f``, and the relevant heat
    capacity is the ice one. Ice holds its water more tightly: at 253 K the
    saturation pressure over ice is about 82 percent of that over supercooled
    liquid, which is the whole basis of the Wegener-Bergeron-Findeisen process.
    Using the liquid branch down to 172 K, as this did before, overstates the
    saturation humidity of the upper troposphere and the stratosphere.

    Both branches pass through the triple point by construction, so the function
    is continuous there and needs no blending.
    """
    T = np.asarray(temperature, dtype=float)
    if np.any(T <= 0.0):
        raise ConfigurationError("temperatures must be positive")

    over_ice = T < _T_TRIPLE
    latent = np.where(over_ice, _L_SUBLIMATION, _L_VAPORIZATION)
    heat_capacity_difference = np.where(over_ice, _CP_VAPOUR - _CP_ICE, _CP_VAPOUR - _CP_LIQUID)

    exponent = heat_capacity_difference / _R_VAPOUR
    offset = (latent - heat_capacity_difference * _T_TRIPLE) / _R_VAPOUR
    return _E_TRIPLE * (T / _T_TRIPLE) ** exponent * np.exp(offset * (1.0 / _T_TRIPLE - 1.0 / T))


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


def lifting_condensation_pressure(
    temperature: float, specific_humidity: float, pressure: float
) -> float:
    """Pressure in Pa at which a parcel lifted dry-adiabatically saturates.

    The parcel keeps its specific humidity and follows ``T (p / p_0)^kappa``
    until ``q`` reaches ``saturation_specific_humidity`` -- the same relation
    the water module condenses by, so the level is where this model would
    condense the parcel, not where another saturation formula would
    (EQ-LIFTING-CONDENSATION-LEVEL). A parcel already saturated is at its own
    level. Solved by bisection, because the saturation humidity falls
    monotonically along a dry adiabat.

    The lift stops where the parcel would leave the saturation relation's
    declared domain, and a parcel still unsaturated there is refused rather
    than given a level the relation cannot vouch for.
    """
    T0, q0, p0 = float(temperature), float(specific_humidity), float(pressure)
    if T0 <= 0.0 or p0 <= 0.0 or q0 < 0.0:
        raise ConfigurationError("a lifted parcel needs a positive temperature and pressure")

    def saturated(p: float) -> bool:
        T = T0 * (p / p0) ** KAPPA
        return bool(q0 >= saturation_specific_humidity(np.array([T]), np.array([p]))[0])

    if saturated(p0):
        return p0
    floor = p0 * (SATURATION_DOMAIN_FLOOR_K / T0) ** (1.0 / KAPPA)
    if not saturated(floor):
        raise ConfigurationError(
            "the parcel does not saturate before leaving the saturation relation's domain"
        )
    low, high = floor, p0  # saturated at low, not at high
    for _ in range(60):
        middle = 0.5 * (low + high)
        if saturated(middle):
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def evaporation_coefficient_for(sensible_heat_coefficient: float) -> float:
    """The bulk evaporation coefficient implied by the sensible-heat one.

    These were two independent free parameters and they are not independent.
    The same turbulent eddies carry heat and moisture across the surface layer,
    so their exchange coefficients are equal to within the accuracy either is
    known to -- the Reynolds analogy, and the reason bulk formulae are written
    with a single ``C`` in most texts::

        H = rho c_p C |U| (T_s - T_1)     =  k_H (T_s - T_1)
        E = rho     C |U| (q_s - q_1)     =  k_E (q_s - q_1)

    so ``k_E = k_H / c_p`` exactly, with no new number. The configured pair was
    ``k_H = 10 W m^-2 K^-1`` against ``k_E = 2.0e-4 kg m^-2 s^-1``, which is
    **fifty times** apart: the column transported heat fifty times more
    efficiently than moisture using eddies it had declared to be the same ones.
    That is the whole of the Bowen-ratio error, and tying the two together
    turns the ratio from a symptom into a prediction.

    Over-determined, which is the point: one coefficient now fixes both surface
    fluxes, so it can be caught being wrong. Measured on the moist column at
    ``k_H = 10``, the Bowen ratio comes out 0.30 against an observed 0.26,
    without being fitted to it.
    """
    if sensible_heat_coefficient < 0.0:
        raise ConfigurationError("the sensible heat coefficient cannot be negative")
    return sensible_heat_coefficient / _CP_DRY


def deposition_latent_heat_excess(temperature: np.ndarray) -> np.ndarray:
    """The extra latent heat released where the condensate forms as ice, in J kg^-1.

    Below the triple point the stable condensed phase is ice, so vapour deposits
    rather than condenses and releases the latent heat of sublimation, ``L_v +
    L_f``. `saturation_vapour_pressure` has saturated over the stable phase
    since the ice branch was added; the energy bookkeeping did not follow, and
    the gap is the whole of ``L_f`` -- 11.8 percent of ``L_v``, discontinuously,
    the moment the condensation is deposition.

    Returned as the *excess* over ``L_v`` rather than as the total, because that
    is what has to be accounted separately. The reference latent heat stays
    constant, so the module's energy store does not change when a layer's
    temperature changes with no water moving -- which a temperature-dependent
    store would do, and which no transfer would balance.

    The excess is real energy and it comes from somewhere. Water leaves the
    surface as liquid, absorbing ``L_v``; it deposits aloft as ice, releasing
    ``L_v + L_f``; and the ice that falls back melts at the surface, absorbing
    ``L_f`` again. Over the cycle it is the surface that pays the fusion energy
    and the atmosphere that receives it, which is what freezing aloft and
    melting at the ground physically do.
    """
    T = np.asarray(temperature, dtype=float)
    return np.where(T < _T_TRIPLE, _L_FUSION, 0.0)


def saturation_humidity_within_domain(
    temperature: np.ndarray, pressure: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Saturation humidity where the relation holds, and exactly zero below it.

    Returns the humidity and the mask of layers it was actually evaluated at, so
    a module can publish the temperatures the relation was applied to and let the
    domain check judge it there rather than over the whole profile.

    This exists because refining the grid walks the column out of the relation's
    own declared domain. A radiative-convective top layer gets colder as it gets
    thinner, and at 20 layers the greenhouse column reaches 159 K against a floor
    of 170. VAL-002 asks for the refinement and `EQ-CLAUSIUS-CLAPEYRON-
    SATURATION` refuses the result, and both are right: the relation really is
    13.6 percent wrong by 150 K.

    Declining to place vapour below the floor is a statement about the closure
    rather than a repair of the relation, and the price of it is measured: the
    water in layers below 170 K is one part in ten million of the column, from
    12 to 40 layers, because the saturation mixing ratio at the floor is already
    of order 1e-9. The real upper atmosphere is not in equilibrium with ice
    anyway -- its few parts per million of water are set by transport and by
    methane oxidation, neither of which this column contains -- so a saturation
    closure has nothing to say up there and now says nothing.
    """
    T = np.asarray(temperature, dtype=float)
    evaluated = T >= SATURATION_DOMAIN_FLOOR_K
    humidity = np.zeros_like(T)
    if np.any(evaluated):
        humidity[evaluated] = saturation_specific_humidity(
            T[evaluated], np.asarray(pressure, dtype=float)[evaluated]
        )
    return humidity, evaluated


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
