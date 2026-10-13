"""Mass-flux convection: an entraining plume and the environment's compensating subsidence.

Why this exists
---------------

The column's free troposphere sat at RH 0.89 to 0.95 because convection only
mixed water upward and nothing in the column descended
(``docs/protocols/free-troposphere-moisture-budget.md``): the moist convection
was a lapse-rate adjustment with uniform vapour mixing. Real convection rises in
narrow plumes, and the rest of the column sinks to compensate, carrying drier air
down from above. Only the plume's detrainment and its evaporating condensate
re-moisten it (Yanai et al. 1973; Arakawa & Schubert 1974).

The physics (``docs/protocols/mass-flux-moisture-transport.md``)
------------------------------------------------------------------

One entraining plume per step, launched from the layer of largest moist static
energy below 850 hPa, which is also its cloud base (Zhang & McFarlane 1995). It
takes all its mass from that layer, rises at constant mass flux entraining and
detraining at the same fractional rate (Tiedtke 1989; Romps 2010 for the rate),
and detrains everything in the layer where its integrated buoyancy work runs out.
Within each layer it mixes its inflow with entrained air, conserving moist static
energy and total water exactly, then condenses to saturation. The condensate leaves
it, and a measured fraction evaporates back into the environment.

The environment sinks at exactly the plume's mass flux at every interface, so
subsidence has no amplitude of its own. Dry static energy and vapour move by one
flux-form operator, upwind. Entrainment draws the start-of-step environment, so
the plume's budget identities hold exactly under subcycling.

The amplitude is a dilute-CAPE closure. The cloud-base mass flux removes the
plume's own CAPE over a declared time, at the rate the same operator removes it
per unit mass flux. Nothing in it reads a humidity, a water amount or a
precipitation rate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from mfespf.core.errors import ConfigurationError
from mfespf.physics.provenance import (
    Milestones,
    condense,
    evaporate,
    form_precipitation,
    surface_arrival,
)
from mfespf.physics.thermodynamics import saturation_specific_humidity

__all__ = [
    "CAPE_TRIGGER_J_KG",
    "CLOUD_WATER_CONVERSION_LAND_PER_M",
    "CLOUD_WATER_CONVERSION_PER_M",
    "CONDENSATE_PATHWAYS",
    "CONVECTIVE_ADJUSTMENT_TIME_S",
    "DETRAINED_STATES",
    "MAX_SUBSTEPS",
    "PLUME_ENTRAINMENT_PER_M",
    "PRECIPITATION_EVAPORATIONS",
    "PRECIPITATION_EVAPORATION_EFFICIENCY",
    "SOURCE_LAYER_FLOOR_PA",
    "Convection",
    "MassFluxStep",
    "Plume",
    "apply_mass_flux",
    "cloud_base_mass_flux",
    "convect",
    "deposit_evaporation",
    "entraining_plume",
    "evaporation_capacity",
    "fall_and_evaporate",
    "provenance_substep",
    "unit_rates",
]

#: The plume's fractional entrainment rate, m^-1: Romps (2010, J. Atmos. Sci.
#: 67, 1908), measured directly in LES of deep convection, "~1 km^-1 throughout
#: the depth of the troposphere". Declared range 0.5-2 km^-1: Romps finds bulk
#: methods underestimate it by about a factor of 2.
PLUME_ENTRAINMENT_PER_M = 1.0e-3

#: Deep convection's trigger, J kg^-1: CAM5 ``zm_conv.F90`` ``capelmt``.
CAPE_TRIGGER_J_KG = 70.0

#: The time the closure removes dilute CAPE over, s: CAM5 ``zm_conv.F90``
#: ``tau``. A closure hypothesis, swept over 1800-14400 s in the protocol.
CONVECTIVE_ADJUSTMENT_TIME_S = 3600.0

#: The source layer is sought among layers centred below this pressure, Pa.
SOURCE_LAYER_FLOOR_PA = 85000.0

#: The state detrained air enters the environment in
#: (``docs/protocols/detrainment-thermodynamic-state.md``):
#:
#: - ``"plume"``: the plume's own s_u and q_u (D0, the original scheme);
#: - ``"mixed_equilibrium"``: as D0, then each receiving layer brought to
#:   saturation at constant moist static energy after every substep, the
#:   condensate detrained explicitly (D1);
#: - ``"environment_saturation"``: vapour capped at the receiving layer's
#:   saturation, the excess detrained as condensate with its latent heat left
#:   in s (D2, a diagnostic limit).
#:
#: In D1 and D2 the detrained condensate joins the plume's condensate in that
#: layer and follows the same evaporating-fraction pathway.
DETRAINED_STATES = ("plume", "mixed_equilibrium", "environment_saturation")

#: Where the evaporating fraction of convective condensate evaporates
#: (``docs/protocols/convective-reevaporation-placement.md``):
#:
#: - ``"cascade"``: where it formed, up to the layer's deficit before the
#:   evaporation's own cooling, the rest passed down (E0, the original);
#: - ``"wet_bulb_cascade"``: the same, limited by the wet-bulb capacity (E1);
#: - ``"falling_rate_law"``: precipitation falls interface by interface and
#:   evaporates under CAM5 Zhang-McFarlane's rate law (E2);
#: - ``"below_cloud_base"``: E0's rule, only below cloud base (E3, a
#:   diagnostic limit).
#:
#: In every one, the evaporating fraction is the most that may evaporate.
PRECIPITATION_EVAPORATIONS = (
    "cascade",
    "wet_bulb_cascade",
    "falling_rate_law",
    "below_cloud_base",
    "none",
    "rate_law_only",
)
# ``"none"`` evaporates no precipitation, and ``"rate_law_only"`` is
# ``"falling_rate_law"`` with no evaporating-fraction cap
# (``docs/protocols/convective-condensate-ledger.md``).

#: What becomes of the plume's condensate (convective-condensate-ledger.md):
#:
#: - ``"immediate"``: it leaves the plume in the layer it forms (the original);
#: - ``"carried"``: the plume carries it as cloud water, converts it to
#:   precipitation at ``c0`` per metre and detrains the rest, as CAM5
#:   Zhang-McFarlane's ``cldprp``.
CONDENSATE_PATHWAYS = ("immediate", "carried")

#: c0, m^-1: cloud water converted to precipitation per metre of ascent,
#: ``rprd = c0 M ql`` with ``ql = ql1 / (1 + dz c0)``, CAM5 Zhang-McFarlane
#: ``cldprp`` at the values its WRF port sets (``c0_ocn = 0.0450``,
#: ``c0_lnd = 0.0059``). The ocean value is nominal for a column over ocean;
#: the land value is the other end of the declared range.
CLOUD_WATER_CONVERSION_PER_M = 0.045
CLOUD_WATER_CONVERSION_LAND_PER_M = 0.0059

#: ke in e = ke (1 - q/q_s) sqrt(F), (kg m^-2 s^-1)^-1/2 s^-1: the evaporation
#: of convective precipitation in CAM5's Zhang-McFarlane scheme
#: (``zm_conv.F90``, ``zm_conv_evap``, "a Sundqvist type algorithm"), at the
#: value the WRF port of that scheme sets; CAM reads it from ``zmconv_ke``.
PRECIPITATION_EVAPORATION_EFFICIENCY = 1.0e-6


@dataclass(frozen=True)
class Plume:
    """One entraining plume, per unit cloud-base mass flux.

    Layers run from the surface (0) up; interfaces from the surface (0) to the
    top (n). ``flux`` is the plume's mass flux at each interface, and
    ``entrainment`` and ``detrainment`` are each layer's exchange, all per unit
    cloud-base flux. ``s_u`` and ``q_u`` are the dry static energy and vapour
    of the air the plume leaves each layer with, NaN where it does not reach.
    ``condensate`` is the water it condenses in each layer, per unit cloud-base
    flux. Its ``buoyancy`` includes the loading of that condensate before it
    falls out.

    Carrying cloud water (``condensate_pathway="carried"``), ``condensate`` is
    the net condensation in each layer, which can be negative where entrained
    air evaporates carried water. ``liquid`` is the cloud water the plume leaves
    each layer with, per unit mass; ``conversion`` the precipitation it forms,
    ``formation`` the condensate it forms and ``cloud_evaporation`` the carried
    water it evaporates, each per unit cloud-base flux. ``None`` when the
    condensate leaves at once.
    """

    source: int
    lfc: int
    top: int
    flux: np.ndarray
    entrainment: np.ndarray
    detrainment: np.ndarray
    s_u: np.ndarray
    q_u: np.ndarray
    temperature: np.ndarray
    condensate: np.ndarray
    buoyancy: np.ndarray
    cape: float
    liquid: np.ndarray | None = field(default=None)
    top_pressure_continuous: float = field(default=float("nan"))
    conversion: np.ndarray | None = field(default=None)
    formation: np.ndarray | None = field(default=None)
    cloud_evaporation: np.ndarray | None = field(default=None)

    @property
    def active(self) -> bool:
        return self.top >= 0


def _layer_thickness(
    temperature: np.ndarray, interface_pressure: np.ndarray, gas_constant: float, gravity: float
) -> np.ndarray:
    return (
        gas_constant
        * temperature
        / gravity
        * np.log(interface_pressure[:-1] / interface_pressure[1:])
    )


def _saturate(
    h: float, q_total: float, z: float, p: float, cp: float, latent: float, gravity: float
) -> tuple[float, float]:
    """The temperature and vapour of air with moist static energy ``h`` and total
    water ``q_total`` at height ``z`` and pressure ``p``, after condensing to
    saturation if it is above it. ``h`` is conserved: the condensate takes no
    energy with it."""
    temperature = (h - gravity * z - latent * q_total) / cp
    saturation = float(saturation_specific_humidity(np.array([temperature]), np.array([p]))[0])
    if q_total <= saturation:
        return temperature, q_total
    # Newton on cp T + L q_s(T) = h - g z, bracketed between the unsaturated
    # temperature (too cold: condensing warms it) and the temperature with all
    # the water condensed.
    target = h - gravity * z
    low, high = temperature, (target) / cp
    t = temperature
    for _ in range(60):
        qs = float(saturation_specific_humidity(np.array([t]), np.array([p]))[0])
        f = cp * t + latent * qs - target
        if abs(f) < 1.0e-9 * abs(target):
            break
        dq = (
            float(saturation_specific_humidity(np.array([t + 1.0e-3]), np.array([p]))[0]) - qs
        ) / 1.0e-3
        step = t - f / (cp + latent * dq)
        if f > 0.0:
            high = t
        else:
            low = t
        t = step if low < step < high else 0.5 * (low + high)
    qs = float(saturation_specific_humidity(np.array([t]), np.array([p]))[0])
    return t, min(qs, q_total)


def entraining_plume(
    *,
    temperature: np.ndarray,
    vapour: np.ndarray,
    layer_mass: np.ndarray,
    heat_capacity: np.ndarray,
    pressure: np.ndarray,
    interface_pressure: np.ndarray,
    entrainment_rate: float,
    latent_heat: float,
    gravity: float,
    gas_constant: float,
    source_floor: float = SOURCE_LAYER_FLOOR_PA,
    condensate_pathway: str = "immediate",
    conversion_rate: float = CLOUD_WATER_CONVERSION_PER_M,
) -> Plume:
    """The plume of mass-flux-moisture-transport.md and its second erratum.

    Launched from the layer of largest moist static energy below
    ``source_floor``, which is its cloud base: it takes its whole mass from that
    layer. It rises at constant mass flux with δ = ε, and detrains everything in
    the layer where its integrated buoyancy work from the level of free
    convection first returns to zero or below after its level of neutral
    buoyancy. CAPE_d is the positive area between the level of free convection
    and the level of neutral buoyancy.

    ``condensate_pathway="carried"`` carries the condensate as cloud water,
    brings each layer's mixture of moist static energy and total water to
    equilibrium, converts ``ql / (1 + c0 dz)`` of it implicitly and loads the
    buoyancy with what it carries on.
    """
    if condensate_pathway not in CONDENSATE_PATHWAYS:
        raise ConfigurationError(
            f"condensate_pathway {condensate_pathway!r} is not one of {CONDENSATE_PATHWAYS}"
        )
    carried = condensate_pathway == "carried"
    latent = latent_heat
    T = np.asarray(temperature, dtype=float)
    m = np.asarray(layer_mass, dtype=float)
    q = np.asarray(vapour, dtype=float) / m
    cp = np.asarray(heat_capacity, dtype=float) / m
    p = np.asarray(pressure, dtype=float)
    pi = np.asarray(interface_pressure, dtype=float)
    n = T.size
    dz = _layer_thickness(T, pi, gas_constant, gravity)
    z = np.concatenate(([0.0], np.cumsum(dz)))[:-1] + 0.5 * dz
    s_e = cp * T + gravity * z
    h_e = s_e + latent * q

    nan = np.full(n, np.nan)
    empty = Plume(
        source=-1,
        lfc=-1,
        top=-1,
        flux=np.zeros(n + 1),
        entrainment=np.zeros(n),
        detrainment=np.zeros(n),
        s_u=nan.copy(),
        q_u=nan.copy(),
        temperature=nan.copy(),
        condensate=np.zeros(n),
        buoyancy=nan.copy(),
        cape=0.0,
    )
    candidates = np.flatnonzero(p > source_floor)
    if candidates.size == 0:
        return empty
    source = int(candidates[np.argmax(h_e[candidates])])

    # The plume's properties layer by layer, per unit cloud-base flux.
    s_u, q_u, t_u = nan.copy(), nan.copy(), nan.copy()
    condensate = np.zeros(n)
    buoyancy = nan.copy()
    entrained = np.zeros(n)
    liquid = nan.copy()
    conversion = np.zeros(n)
    formation = np.zeros(n)
    cloud_evaporation = np.zeros(n)

    def leave(k: int, h_mix: float, q_mix: float, mass: float, carried_in: float = 0.0) -> None:
        t, qv = _saturate(h_mix, q_mix, z[k], p[k], cp[k], latent, gravity)
        t_u[k], q_u[k] = t, qv
        s_u[k] = h_mix - latent * qv
        condensed = q_mix - qv
        if not carried:
            condensate[k] = mass * condensed
            buoyancy[k] = gravity * ((t - T[k]) / T[k] + 0.61 * (qv - q[k]) - condensed)
            return
        # ``condensed`` is the equilibrium cloud water of the mixture, whose
        # total water includes what was carried in, diluted.
        net = condensed - carried_in
        formation[k] = mass * max(net, 0.0)
        cloud_evaporation[k] = mass * max(-net, 0.0)
        kept = condensed / (1.0 + conversion_rate * dz[k])
        conversion[k] = mass * (condensed - kept)
        liquid[k] = kept
        condensate[k] = mass * net
        buoyancy[k] = gravity * ((t - T[k]) / T[k] + 0.61 * (qv - q[k]) - kept)

    entrained[source] = 1.0
    leave(source, h_e[source], q[source], 1.0)
    h_u, q_t = h_e[source], q_u[source]
    if carried:
        q_t = q_u[source] + liquid[source]
    lfc = lnb = top = -1
    work = 0.0
    before = 0.0
    cape = 0.0
    for k in range(source + 1, n):
        e = entrainment_rate * dz[k]
        h_mix = (h_u + e * h_e[k]) / (1.0 + e)
        q_mix = (q_t + e * q[k]) / (1.0 + e)
        entrained[k] = e
        if carried:
            leave(k, h_mix, q_mix, 1.0 + e, liquid[k - 1] / (1.0 + e))
            h_u, q_t = h_mix, q_u[k] + liquid[k]
        else:
            leave(k, h_mix, q_mix, 1.0 + e)
            h_u, q_t = h_mix, q_u[k]
        b = buoyancy[k]
        if lfc < 0:
            if b > 0.0:
                lfc = k
            else:
                continue
        before = work
        work += b * dz[k]
        if lnb < 0:
            if b > 0.0:
                cape += b * dz[k]
                continue
            lnb = k
        if work <= 0.0:
            top = k
            break
    if lfc < 0:
        return empty
    # Where the buoyancy work crosses zero between the last layer it was
    # positive in and the terminal layer, interpolated in pressure between
    # their centres (convective-temporal-convergence.md, T0). A diagnostic:
    # the detrainment is still deposited in the terminal layer.
    if top >= 0 and before > 0.0:
        crossing = before / (before - work)
        top_continuous = float(p[top - 1] + crossing * (p[top] - p[top - 1]))
    else:
        top_continuous = float(p[top]) if top >= 0 else float(p[n - 1])
    if top < 0:
        top = n - 1
    # Above the termination layer the plume does not exist.
    for arr in (s_u, q_u, t_u, buoyancy, liquid):
        arr[top + 1 :] = np.nan
    for arr in (condensate, entrained, conversion, formation, cloud_evaporation):
        arr[top + 1 :] = 0.0

    flux = np.zeros(n + 1)
    flux[source + 1 : top + 1] = 1.0
    detrained = np.zeros(n)
    detrained[source + 1 : top] = entrained[source + 1 : top]
    detrained[top] = 1.0 + entrained[top]
    return Plume(
        source=source,
        lfc=lfc,
        top=top,
        flux=flux,
        entrainment=entrained,
        detrainment=detrained,
        s_u=s_u,
        q_u=q_u,
        temperature=t_u,
        condensate=condensate,
        buoyancy=buoyancy,
        cape=float(cape),
        top_pressure_continuous=top_continuous,
        **(
            {
                "liquid": liquid,
                "conversion": conversion,
                "formation": formation,
                "cloud_evaporation": cloud_evaporation,
            }
            if carried
            else {}
        ),
    )


def deposit_evaporation(amounts: np.ndarray, deficits: np.ndarray) -> tuple[np.ndarray, float]:
    """Evaporate ``amounts`` (kg m^-2 per layer) where they formed, up to each
    layer's saturation deficit; what exceeds it falls to the next layer below and
    evaporates there under the same limit. Returns the evaporation per layer and
    what reaches the surface."""
    evaporated = np.zeros_like(np.asarray(amounts, dtype=float))
    falling = 0.0
    for k in range(evaporated.size - 1, -1, -1):
        falling += float(amounts[k])
        take = min(falling, max(float(deficits[k]), 0.0))
        evaporated[k] = take
        falling -= take
    return evaporated, falling


def evaporation_capacity(
    temperature: np.ndarray,
    specific_humidity: np.ndarray,
    pressure: np.ndarray,
    cp: np.ndarray,
    latent_heat: float,
) -> np.ndarray:
    """The most water each layer can take by evaporation, per unit mass: the
    largest d with q + d <= q_s(T - L d / c_p, p), the wet-bulb limit. Bisection
    to 2^-60 of the bracket [0, min(max(q_s(T) - q, 0), c_p T / (2 L))]: the
    excess is negative at 0 and not negative at the top.

    The second bound keeps every trial temperature at T/2 or above. Near the
    model top q_s(T) can exceed the layer's whole capacity to cool (e_s
    approaches p), and a bracket of q_s(T) - q alone drove the trial
    temperature below zero (convective-reevaporation-placement.md, erratum).
    At T/2, q_s is negligible against q + c_p T / (2 L); that is checked, and
    refused if it fails, rather than assumed."""
    T = np.asarray(temperature, dtype=float)
    q = np.asarray(specific_humidity, dtype=float)
    p = np.asarray(pressure, dtype=float)
    low = np.zeros_like(T)
    high = np.minimum(
        np.maximum(saturation_specific_humidity(T, p) - q, 0.0), cp * T / (2.0 * latent_heat)
    )
    unbracketed = (high > 0.0) & (
        q + high - saturation_specific_humidity(T - latent_heat * high / cp, p) < 0.0
    )
    if np.any(unbracketed):
        raise ConfigurationError(
            "the wet-bulb capacity is not bracketed in layers "
            f"{np.flatnonzero(unbracketed).tolist()}"
        )
    for _ in range(60):
        mid = 0.5 * (low + high)
        over = q + mid - saturation_specific_humidity(T - latent_heat * mid / cp, p) > 0.0
        high = np.where(over, mid, high)
        low = np.where(over, low, mid)
    return low


def fall_and_evaporate(
    formed: np.ndarray,
    *,
    potential: float,
    capacity: np.ndarray,
    specific_humidity: np.ndarray,
    saturation: np.ndarray,
    layer_mass: np.ndarray,
    dt: float,
    efficiency: float = PRECIPITATION_EVAPORATION_EFFICIENCY,
) -> tuple[np.ndarray, float]:
    """E2: condensate falls interface by interface over ``dt``. In each layer the
    flux from above, F, evaporates at ke (1 - q/q_s) sqrt(F) per unit mass,
    limited by the layer's capacity (kg m^-2), by what falls into it, and by
    what the evaporating fraction ``potential`` still allows of the condensate
    formed above. Rain formed in a layer is not evaporated there. Returns the
    evaporation per layer and what reaches the surface."""
    evaporated = np.zeros_like(np.asarray(formed, dtype=float))
    falling = 0.0
    allowed = 0.0
    for k in range(evaporated.size - 1, -1, -1):
        if falling > 0.0:
            dryness = max(0.0, 1.0 - float(specific_humidity[k]) / float(saturation[k]))
            rate = efficiency * dryness * np.sqrt(falling / dt) * float(layer_mass[k]) * dt
            take = min(rate, max(float(capacity[k]), 0.0), falling, max(allowed, 0.0))
            evaporated[k] = take
            falling -= take
            allowed -= take
        falling += float(formed[k])
        allowed += potential * float(formed[k])
    return evaporated, falling


def _operator(
    plume: Plume,
    s: np.ndarray,
    q: np.ndarray,
    s_det: np.ndarray | None = None,
    q_det: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """m dphi/dt per unit cloud-base flux for phi = s and q: subsidence upwind, and
    the exchange with a plume launched into this same environment, so the plume's
    budget identities hold exactly (the protocol's third erratum). The detrained
    air carries ``s_det`` and ``q_det``, the plume's own unless given."""
    M = plume.flux
    E, D = plume.entrainment, plume.detrainment
    s_u = np.nan_to_num(plume.s_u) if s_det is None else s_det
    q_u = np.nan_to_num(plume.q_u) if q_det is None else q_det
    rate_s = np.zeros_like(s)
    rate_q = np.zeros_like(q)
    rate_s[:-1] += M[1:-1] * s[1:]
    rate_s -= M[:-1] * s
    rate_q[:-1] += M[1:-1] * q[1:]
    rate_q -= M[:-1] * q
    rate_s += -E * s + D * s_u
    rate_q += -E * q + D * q_u
    return rate_s, rate_q


def unit_rates(plume: Plume, s: np.ndarray, q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The environment's m ds/dt and m dq/dt at unit cloud-base flux, now."""
    return _operator(plume, s, q)


def _advective_split(
    plume: Plume, phi: np.ndarray, phi_det: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """The operator in advective form, per unit cloud-base flux: subsidence
    M_(k+1) (phi_(k+1) - phi_k) and net detrainment D_k (phi_det,k - phi_k).
    Their sum is the flux form because M_(k+1) - M_k = E_k - D_k; they are
    computed independently, so a test can check that."""
    M, D = plume.flux, plume.detrainment
    subsidence = np.zeros_like(phi)
    subsidence[:-1] = M[1:-1] * (phi[1:] - phi[:-1])
    return subsidence, D * (phi_det - phi)


@dataclass(frozen=True)
class MassFluxStep:
    """One step's convective change, as the protocol's erratum divides it.

    ``energy`` is each layer's change in m s (J m^-2) less the latent heat of
    its net condensation: a redistribution, summing to zero. ``transport`` is
    the vapour moved (kg m^-2), also summing to zero. ``condensation`` is the
    net condensation c = P - e in each layer (kg m^-2), which the water module
    turns into precipitation and latent heat; its column sum is the surface
    precipitation. ``formed`` includes the detrained condensate.

    The rest decompose it (detrainment-thermodynamic-state.md), all kg m^-2 or
    J m^-2 over the step: the operator's subsidence and net-detrainment parts
    for vapour and dry static energy; the vapour and condensate detrained,
    gross; and the supersaturation detrainment leaves in a receiving layer at
    constant moist static energy, removed by convection only in D1.
    """

    energy: np.ndarray
    transport: np.ndarray
    condensation: np.ndarray
    evaporation: np.ndarray
    formed: np.ndarray
    surface_precipitation: float
    substeps: int
    subsidence_vapour: np.ndarray
    detrainment_vapour: np.ndarray
    subsidence_energy: np.ndarray
    detrainment_energy: np.ndarray
    detrained_vapour: np.ndarray
    detrained_condensate: np.ndarray
    post_detrainment_excess: np.ndarray
    precipitation_entering: np.ndarray
    evaporation_limit_time: np.ndarray
    formation_gross: np.ndarray
    cloud_evaporation: np.ndarray
    rain_evaporation: np.ndarray
    precipitation_formed: np.ndarray
    detrained_cloud: np.ndarray
    detrained_mass: np.ndarray
    plume_mass: np.ndarray


#: The decomposition fields, summed over substeps. ``precipitation_entering``
#: is the convective condensate falling into each layer from above (kg m^-2);
#: ``evaporation_limit_time`` is each layer's evaporation limit times the
#: substep's length (kg m^-2 s), so that its sum over a step divided by the step
#: is the step-mean limit.
_PARTS = (
    "subsidence_vapour",
    "detrainment_vapour",
    "subsidence_energy",
    "detrainment_energy",
    "detrained_vapour",
    "detrained_condensate",
    "post_detrainment_excess",
    "precipitation_entering",
    "evaporation_limit_time",
    "formation_gross",
    "cloud_evaporation",
    "rain_evaporation",
    "precipitation_formed",
    "detrained_cloud",
    "detrained_mass",
    "plume_mass",
)
# The ledger (convective-condensate-ledger.md), kg m^-2 over the step:
# ``formation_gross`` is C_form; ``cloud_evaporation`` E_cloud (in the plume and
# on detrainment, or, when condensate leaves at once, the local-first share of
# the cascade); ``rain_evaporation`` E_precip; ``precipitation_formed`` C_precip
# where it forms; ``detrained_cloud`` the cloud water detrained.


def apply_mass_flux(
    plume: Plume,
    *,
    cloud_base_flux: float,
    dt: float,
    temperature: np.ndarray,
    vapour: np.ndarray,
    layer_mass: np.ndarray,
    heat_capacity: np.ndarray,
    pressure: np.ndarray,
    interface_pressure: np.ndarray,
    latent_heat: float,
    gravity: float,
    gas_constant: float,
    evaporating_fraction: float,
    detrained_state: str = "plume",
    precipitation_evaporation: str = "cascade",
) -> MassFluxStep:
    """One exact, conservative step of the plume's transport at a fixed mass flux,
    then the measured fraction of its condensate evaporated where it formed and
    below. Refused if any layer would lose more than half its mass: the
    integrator (``convect``) chooses the step so that it never does.

    ``detrained_state`` is one of ``DETRAINED_STATES``, and
    ``precipitation_evaporation`` one of ``PRECIPITATION_EVAPORATIONS``.
    """
    if detrained_state not in DETRAINED_STATES:
        raise ConfigurationError(
            f"detrained_state {detrained_state!r} is not one of {DETRAINED_STATES}"
        )
    if precipitation_evaporation not in PRECIPITATION_EVAPORATIONS:
        raise ConfigurationError(
            f"precipitation_evaporation {precipitation_evaporation!r} is not one of "
            f"{PRECIPITATION_EVAPORATIONS}"
        )
    carried = plume.liquid is not None
    if carried and (
        detrained_state != "plume" or precipitation_evaporation not in ("none", "rate_law_only")
    ):
        raise ConfigurationError(
            "carried cloud water detrains with the plume's own state and evaporates its "
            "precipitation by the rate law or not at all"
        )
    n = np.asarray(layer_mass).size
    zero = np.zeros(n)
    if not plume.active or cloud_base_flux <= 0.0 or dt <= 0.0:
        return MassFluxStep(
            zero,
            zero.copy(),
            zero.copy(),
            zero.copy(),
            zero.copy(),
            0.0,
            0,
            *(zero.copy() for _ in _PARTS),
        )
    T = np.asarray(temperature, dtype=float)
    m = np.asarray(layer_mass, dtype=float)
    C = np.asarray(heat_capacity, dtype=float)
    cp = C / m
    pi = np.asarray(interface_pressure, dtype=float)
    dz = _layer_thickness(T, pi, gas_constant, gravity)
    z = np.concatenate(([0.0], np.cumsum(dz)))[:-1] + 0.5 * dz
    s0 = cp * T + gravity * z
    q0 = np.asarray(vapour, dtype=float) / m
    outflow = cloud_base_flux * (plume.flux[:-1] + plume.entrainment) * dt
    if np.any(outflow > 0.5 * m * (1.0 + 1.0e-12)):
        raise ConfigurationError(
            "a mass-flux step would take more than half of a layer's mass; the "
            "integrator must shorten it"
        )
    p = np.asarray(pressure, dtype=float)
    scale = cloud_base_flux * dt
    D = plume.detrainment
    detraining = (D > 0.0) & ~np.isnan(plume.q_u)
    s_det = np.nan_to_num(plume.s_u)
    q_det = np.nan_to_num(plume.q_u)
    held = np.zeros(n)
    if detrained_state == "environment_saturation":
        # D2: no more vapour than the receiving layer holds at its own
        # temperature; the rest detrains as condensate, its latent heat in s so
        # that D h_u still enters the layer exactly.
        capped = np.where(detraining, np.minimum(q_det, saturation_specific_humidity(T, p)), q_det)
        excess = q_det - capped
        s_det = s_det + latent_heat * excess
        q_det = capped
        held = scale * D * excess
    rate_s, rate_q = _operator(plume, s0, q0, s_det, q_det)
    s = s0 + scale * rate_s / m
    q = q0 + scale * rate_q / m
    if np.any(q < 0.0):
        raise ConfigurationError("mass-flux transport drew a layer's vapour below zero")
    subsidence_s, detrainment_s = _advective_split(plume, s0, s_det)
    subsidence_q, detrainment_q = _advective_split(plume, q0, q_det)
    # What detrainment leaves above saturation in a receiving layer, brought to
    # saturation at constant h = s + L q by the plume's own solve. Measured in
    # every mode; removed here only in D1, where it is detrained condensate.
    post = np.zeros(n)
    for k in (int(j) for j in np.flatnonzero(detraining)):
        h = s[k] + latent_heat * q[k]
        _, held_vapour = _saturate(h, q[k], z[k], p[k], cp[k], latent_heat, gravity)
        post[k] = (q[k] - held_vapour) * m[k]
        if detrained_state == "mixed_equilibrium" and held_vapour < q[k]:
            s[k] = h - latent_heat * held_vapour
            q[k] = held_vapour
    if detrained_state == "mixed_equilibrium":
        held = post.copy()
    formed = scale * plume.condensate + held
    t_after = (s - gravity * z) / cp
    saturation = saturation_specific_humidity(t_after, p)
    deficits = saturation * m - q * m
    detrained_cloud = np.zeros(n)
    if carried:
        # Cloud water detrains into the receiving layer and evaporates there up
        # to its wet-bulb capacity; what it cannot take becomes precipitation in
        # that layer, as there is no prognostic cloud water to hold it.
        assert plume.liquid is not None and plume.conversion is not None
        assert plume.formation is not None
        assert plume.cloud_evaporation is not None
        limit = evaporation_capacity(t_after, q, p, cp, latent_heat) * m
        detrained_cloud = scale * D * np.nan_to_num(plume.liquid)
        taken = np.minimum(detrained_cloud, limit)
        precipitation_formed = scale * plume.conversion + (detrained_cloud - taken)
        if precipitation_evaporation == "rate_law_only":
            moistened = q + taken / m
            cooled = t_after - latent_heat * taken / (cp * m)
            rain, _ = fall_and_evaporate(
                precipitation_formed,
                potential=1.0,
                capacity=limit - taken,
                specific_humidity=moistened,
                saturation=saturation_specific_humidity(cooled, p),
                layer_mass=m,
                dt=dt,
            )
        else:
            rain = np.zeros(n)
        evaporated = taken + rain
        formation_gross = scale * plume.formation
        cloud_evaporation = scale * plume.cloud_evaporation + taken
        rain_evaporation = rain
    elif precipitation_evaporation == "none":
        limit = np.zeros(n)
        evaporated = np.zeros(n)
    elif precipitation_evaporation == "rate_law_only":
        limit = evaporation_capacity(t_after, q, p, cp, latent_heat) * m
        evaporated, _ = fall_and_evaporate(
            formed,
            potential=1.0,
            capacity=limit,
            specific_humidity=q,
            saturation=saturation,
            layer_mass=m,
            dt=dt,
        )
    elif precipitation_evaporation == "cascade":
        limit = deficits
        evaporated, _ = deposit_evaporation(evaporating_fraction * formed, limit)
    elif precipitation_evaporation == "below_cloud_base":
        condensing = np.flatnonzero(plume.condensate > 0.0)
        base = int(condensing[0]) if condensing.size else n
        limit = np.where(np.arange(n) < base, deficits, 0.0)
        evaporated, _ = deposit_evaporation(evaporating_fraction * formed, limit)
    else:
        limit = evaporation_capacity(t_after, q, p, cp, latent_heat) * m
        if precipitation_evaporation == "wet_bulb_cascade":
            evaporated, _ = deposit_evaporation(evaporating_fraction * formed, limit)
        else:
            evaporated, _ = fall_and_evaporate(
                formed,
                potential=evaporating_fraction,
                capacity=limit,
                specific_humidity=q,
                saturation=saturation,
                layer_mass=m,
                dt=dt,
            )
    if not carried:
        formation_gross = formed
        if precipitation_evaporation in ("cascade", "wet_bulb_cascade", "below_cloud_base"):
            # Local first: a layer's own share evaporates before what fell into it.
            cloud_evaporation = np.minimum(evaporating_fraction * formed, np.maximum(limit, 0.0))
        else:
            cloud_evaporation = np.zeros(n)
        rain_evaporation = evaporated - cloud_evaporation
        precipitation_formed = formed - cloud_evaporation
    falling = precipitation_formed - rain_evaporation
    entering = np.concatenate((np.cumsum(falling[::-1])[::-1][1:], [0.0]))
    condensation = formed - evaporated
    vapour_change = (q - q0) * m + evaporated
    energy_change = (s - s0) * m - latent_heat * evaporated
    return MassFluxStep(
        energy=energy_change - latent_heat * condensation,
        transport=vapour_change + condensation,
        condensation=condensation,
        evaporation=evaporated,
        formed=formed,
        surface_precipitation=float(np.sum(formed) - np.sum(evaporated)),
        substeps=1,
        subsidence_vapour=scale * subsidence_q,
        detrainment_vapour=scale * detrainment_q,
        subsidence_energy=scale * subsidence_s,
        detrainment_energy=scale * detrainment_s,
        detrained_vapour=scale * D * q_det,
        detrained_condensate=held,
        post_detrainment_excess=post,
        precipitation_entering=entering,
        evaporation_limit_time=np.maximum(limit, 0.0) * dt,
        formation_gross=formation_gross,
        cloud_evaporation=cloud_evaporation,
        rain_evaporation=rain_evaporation,
        precipitation_formed=precipitation_formed,
        detrained_cloud=detrained_cloud,
        detrained_mass=scale * D,
        plume_mass=scale * plume.flux[:-1],
    )


def cloud_base_mass_flux(
    plume: Plume,
    *,
    temperature: np.ndarray,
    vapour: np.ndarray,
    layer_mass: np.ndarray,
    heat_capacity: np.ndarray,
    pressure: np.ndarray,
    interface_pressure: np.ndarray,
    entrainment_rate: float,
    latent_heat: float,
    gravity: float,
    gas_constant: float,
    evaporating_fraction: float,
    adjustment_time: float = CONVECTIVE_ADJUSTMENT_TIME_S,
    trigger: float = CAPE_TRIGGER_J_KG,
    residual: float = 0.0,
    detrained_state: str = "plume",
    precipitation_evaporation: str = "cascade",
    condensate_pathway: str = "immediate",
    conversion_rate: float = CLOUD_WATER_CONVERSION_PER_M,
) -> tuple[float, float]:
    """M_b = max[0, (CAPE_d - C_residual) / (τ_c Γ_M)], and Γ_M.

    Γ_M is the rate unit cloud-base flux removes CAPE_d (erratum 1): the state is
    advanced along the operator's tendencies at M_b = 1 by M_b Δt = 1 kg m^-2,
    including the latent heat of what condenses and evaporates, and CAPE_d is
    re-evaluated on the plume launched into the result. Zero below the trigger,
    and zero if convection would not reduce CAPE_d.
    """
    if not plume.active or plume.cape <= trigger:
        return 0.0, 0.0
    displacement = 1.0  # kg m^-2 of cloud-base mass, M_b dt at M_b = 1 for 1 s
    step = apply_mass_flux(
        plume,
        cloud_base_flux=1.0,
        dt=displacement,
        temperature=temperature,
        vapour=vapour,
        layer_mass=layer_mass,
        heat_capacity=heat_capacity,
        pressure=pressure,
        interface_pressure=interface_pressure,
        latent_heat=latent_heat,
        gravity=gravity,
        gas_constant=gas_constant,
        evaporating_fraction=evaporating_fraction,
        detrained_state=detrained_state,
        precipitation_evaporation=precipitation_evaporation,
    )
    C = np.asarray(heat_capacity, dtype=float)
    # The full change the environment sees: the redistribution and the latent
    # heat the water module will add for the net condensation.
    t_new = (
        np.asarray(temperature, dtype=float) + (step.energy + latent_heat * step.condensation) / C
    )
    w_new = np.asarray(vapour, dtype=float) + step.transport - step.condensation
    displaced = entraining_plume(
        temperature=t_new,
        vapour=w_new,
        layer_mass=layer_mass,
        heat_capacity=heat_capacity,
        pressure=pressure,
        interface_pressure=interface_pressure,
        entrainment_rate=entrainment_rate,
        latent_heat=latent_heat,
        gravity=gravity,
        gas_constant=gas_constant,
        condensate_pathway=condensate_pathway,
        conversion_rate=conversion_rate,
    )
    rate = -(displaced.cape - plume.cape) / displacement
    if rate <= 0.0:
        return 0.0, rate
    return max(0.0, (plume.cape - residual) / (adjustment_time * rate)), rate


def provenance_substep(
    plume: Plume,
    step: MassFluxStep,
    *,
    scale: float,
    classes: np.ndarray,
    layer_mass: np.ndarray,
) -> tuple[np.ndarray, Milestones]:
    """Carry the environment's water provenance classes through one substep
    (fall-path-closure-candidate.md).

    ``classes`` is the environment's vapour by class, (2n + 1, n) kg m^-2, at
    the start of the substep; ``step`` the substep's bulk result, at cloud-base
    mass ``scale`` (M_b dt). The plume carries its vapour and cloud water by
    class, mixing in entrained air and detraining in proportion; it condenses,
    evaporates and converts by the same amounts as the bulk plume. The
    environment's classes then take the same flux-form operator as its vapour,
    the detrained cloud water's evaporation and the falling precipitation's.
    What reaches the surface ends its cohort. Returns the classes after the
    substep and the milestones it crossed.
    """
    if plume.liquid is None or plume.conversion is None:
        raise ConfigurationError("water provenance is carried with the cloud-water pathway only")
    K, n = classes.shape
    m = np.asarray(layer_mass, dtype=float)
    q = classes / m
    unit = Milestones.zeros(n)
    scaled = Milestones.zeros(n)
    detrained_vapour = np.zeros((K, n))
    detrained_cloud = np.zeros((K, n))
    rain_formed = np.zeros((K, n))
    vapour = np.zeros(K)
    cloud = np.zeros(K)
    E, D, M = plume.entrainment, plume.detrainment, plume.flux
    if plume.active:
        for k in range(plume.source, plume.top + 1):
            vapour = vapour + E[k] * q[:, k]
            net = float(plume.condensate[k])
            if net > 0.0:
                vapour, formed = condense(vapour, net, k, unit, large_scale=False)
                cloud = cloud + formed
            elif net < 0.0:
                cloud, freed = evaporate(cloud, -net)
                vapour = vapour + freed
            if plume.conversion[k] > 0.0:
                cloud, rain = form_precipitation(cloud, float(plume.conversion[k]), k, unit)
                rain_formed[:, k] += scale * rain
            mass = 1.0 if k == plume.source else 1.0 + float(E[k])
            share = min(float(D[k]) / mass, 1.0)
            detrained_vapour[:, k] = vapour * share
            detrained_cloud[:, k] = cloud * share
            vapour = vapour - detrained_vapour[:, k]
            cloud = cloud - detrained_cloud[:, k]
    for name in ("first_condensed", "first_precip", "first_precip_location"):
        setattr(scaled, name, scale * getattr(unit, name))

    rate = np.zeros((K, n))
    rate[:, :-1] += M[1:-1] * q[:, 1:]
    rate -= M[:-1] * q
    rate += -E * q + detrained_vapour
    after = classes + scale * rate

    # The detrained cloud water evaporates into the receiving layer as far as
    # the bulk step found it could, and the rest becomes precipitation there.
    plume_evaporation = scale * np.asarray(plume.cloud_evaporation, dtype=float)
    for k in range(n):
        arriving = scale * detrained_cloud[:, k]
        if not np.any(arriving):
            continue
        taken = float(step.cloud_evaporation[k] - plume_evaporation[k])
        left, freed = evaporate(arriving, taken)
        after[:, k] += freed
        if np.sum(left) > 0.0:
            _, rain = form_precipitation(left, float(np.sum(left)), k, scaled)
            rain_formed[:, k] += rain

    # The fall: what comes from above evaporates as the bulk fall did, and rain
    # formed in a layer joins the flux to the layer below.
    falling = np.zeros(K)
    for k in range(n - 1, -1, -1):
        if step.rain_evaporation[k] > 0.0:
            falling, freed = evaporate(falling, float(step.rain_evaporation[k]))
            after[:, k] += freed
        falling = falling + rain_formed[:, k]
    surface_arrival(falling, scaled)
    return after, scaled


#: A step needing more substeps than this is refused, not limited (third erratum).
MAX_SUBSTEPS = 500

#: What follows a window whose closure finds no mass flux
#: (convective-window-integration.md): "stop" ends convection for the rest of
#: the host step, "evaluate" ends only that window (I1).
WINDOW_CONTINUATIONS = ("stop", "evaluate")

#: A convective window's status (convective-window-integration.md, I0).
WINDOW_ACTIVE = 0.0
WINDOW_STOPPED_WITHIN = 1.0
WINDOW_BELOW_TRIGGER = 2.0
WINDOW_NO_PLUME = 3.0
WINDOW_NO_REMOVAL = 4.0
WINDOW_SKIPPED = 5.0


@dataclass(frozen=True)
class Convection:
    """What one step of mass-flux convection did, and how it got there."""

    step: MassFluxStep
    plume: Plume
    cape: float
    cape_removal: float
    first_cloud_base_flux: float
    mean_cloud_base_flux: float
    provenance: np.ndarray | None = None
    milestones: Milestones | None = None
    cape_after: float = float("nan")
    cape_time_mean: float = float("nan")
    cape_window_start: np.ndarray | None = None
    active_time_fraction: float = float("nan")
    window_evaluated: np.ndarray | None = None
    window_status: np.ndarray | None = None
    window_mass: np.ndarray | None = None

    @property
    def cape_window_end(self) -> np.ndarray | None:
        """CAPE at each window's end: the next window's start, the last the
        CAPE convection leaves (convective-window-integration.md)."""
        if self.cape_window_start is None:
            return None
        return np.append(self.cape_window_start[1:], self.cape_after)

    @property
    def cape_pre_convection(self) -> float:
        """The mean over the step's windows of CAPE at each window's start
        (convective-substep-convergence.md)."""
        if self.cape_window_start is None:
            return float("nan")
        return float(np.mean(self.cape_window_start))

    @property
    def cape_post_convection(self) -> float:
        """The mean over the step's windows of CAPE at each window's end:
        nothing else acts within convection, so a window ends at the next
        one's start, and the last at the CAPE convection leaves."""
        if self.cape_window_start is None:
            return float("nan")
        return float(np.mean(np.append(self.cape_window_start[1:], self.cape_after)))


def convect(
    *,
    dt: float,
    temperature: np.ndarray,
    vapour: np.ndarray,
    layer_mass: np.ndarray,
    heat_capacity: np.ndarray,
    pressure: np.ndarray,
    interface_pressure: np.ndarray,
    entrainment_rate: float,
    latent_heat: float,
    gravity: float,
    gas_constant: float,
    evaporating_fraction: float,
    adjustment_time: float = CONVECTIVE_ADJUSTMENT_TIME_S,
    trigger: float = CAPE_TRIGGER_J_KG,
    detrained_state: str = "plume",
    precipitation_evaporation: str = "cascade",
    condensate_pathway: str = "immediate",
    conversion_rate: float = CLOUD_WATER_CONVERSION_PER_M,
    provenance: np.ndarray | None = None,
    substep: float | None = None,
    continuation: str = "stop",
) -> Convection:
    """Integrate the closure over ``dt`` as the rate law it is (third erratum).

    Each substep launches the plume into the current state, closes its mass flux
    on its own CAPE, and takes one exact step of the operator, as long as no
    layer loses more than half its mass. It stops when CAPE falls to the
    trigger or convection stops removing it.
    """
    T = np.array(temperature, dtype=float)
    w = np.array(vapour, dtype=float)
    m = np.asarray(layer_mass, dtype=float)
    C = np.asarray(heat_capacity, dtype=float)
    n = m.size
    keys = ("energy", "transport", "condensation", "evaporation", "formed", *_PARTS)
    totals = {k: np.zeros(n) for k in keys}
    fixed: dict[str, Any] = dict(
        layer_mass=m,
        heat_capacity=C,
        pressure=pressure,
        interface_pressure=interface_pressure,
        latent_heat=latent_heat,
        gravity=gravity,
        gas_constant=gas_constant,
    )
    elapsed, substeps, moved = 0.0, 0, 0.0
    # The step's time-mean CAPE (convective-cape-sampling.md, C0): each
    # substep's starting CAPE times its length, a diagnostic never fed back.
    cape_time = 0.0
    # Each window's starting CAPE and the time convection integrated
    # (convective-substep-convergence.md): diagnostics, never fed back.
    windows = 1 if substep is None else round(dt / substep)
    window_start = np.full(windows, np.nan)
    window, fresh = 0, True
    active_time = 0.0
    # Per window (convective-window-integration.md, I0): whether the closure
    # was evaluated at its start, what it found, and the mass it moved.
    window_evaluated = np.zeros(windows)
    window_status = np.full(windows, WINDOW_SKIPPED)
    window_mass = np.zeros(windows)
    if continuation not in WINDOW_CONTINUATIONS:
        raise ConfigurationError(
            f"unknown convective window continuation {continuation!r}; "
            f"expected one of {WINDOW_CONTINUATIONS}"
        )
    idle_from: float | None = None
    idle_cape = 0.0
    first: tuple[Plume, float, float] | None = None
    # A common convective step (convective-temporal-convergence.md, T1): the
    # host step is integrated in windows of ``substep``, every substep ending at
    # the next window boundary at the latest, so the closure is re-evaluated at
    # least that often whatever the host step.
    if substep is not None:
        if not substep > 0.0 or abs(dt / substep - round(dt / substep)) > 1.0e-9:
            raise ConfigurationError(
                f"the convective substep {substep:g} s must divide the step {dt:g} s exactly"
            )
    window_end = dt if substep is None else min(dt, substep)
    classes = None if provenance is None else np.array(provenance, dtype=float)
    milestones = None if provenance is None else Milestones.zeros(n)
    if classes is not None and (
        condensate_pathway != "carried"
        or detrained_state != "plume"
        or precipitation_evaporation not in ("none", "rate_law_only")
    ):
        raise ConfigurationError(
            "water provenance follows the carried cloud-water pathway with the rate law or "
            "no precipitation evaporation"
        )
    while elapsed < dt * (1.0 - 1.0e-12):
        plume = entraining_plume(
            temperature=T,
            vapour=w,
            entrainment_rate=entrainment_rate,
            condensate_pathway=condensate_pathway,
            conversion_rate=conversion_rate,
            **fixed,
        )
        flux, removal = cloud_base_mass_flux(
            plume,
            temperature=T,
            vapour=w,
            entrainment_rate=entrainment_rate,
            evaporating_fraction=evaporating_fraction,
            adjustment_time=adjustment_time,
            trigger=trigger,
            detrained_state=detrained_state,
            precipitation_evaporation=precipitation_evaporation,
            condensate_pathway=condensate_pathway,
            conversion_rate=conversion_rate,
            **fixed,
        )
        if first is None:
            first = (plume, removal, flux)
        if fresh and window < windows:
            window_start[window] = plume.cape
            fresh = False
            window_evaluated[window] = 1.0
            if flux > 0.0:
                window_status[window] = WINDOW_ACTIVE
            elif not plume.active:
                window_status[window] = WINDOW_NO_PLUME
            elif plume.cape <= trigger:
                window_status[window] = WINDOW_BELOW_TRIGGER
            else:
                window_status[window] = WINDOW_NO_REMOVAL
        elif flux <= 0.0 and window < windows:
            window_status[window] = WINDOW_STOPPED_WITHIN
        if flux <= 0.0:
            # "stop" ends convection for the rest of the host step; "evaluate"
            # ends only this window, and the next is closed on the state as it
            # stands (convective-window-integration.md, I1). Nothing else acts
            # within this call, so the idle time holds this CAPE.
            if continuation == "evaluate" and substep is not None and window_end < dt:
                if idle_from is None:
                    idle_from, idle_cape = elapsed, plume.cape
                elapsed = window_end
                window_end = min(dt, window_end + substep)
                window, fresh = window + 1, True
                continue
            break
        if idle_from is not None:
            cape_time += (elapsed - idle_from) * idle_cape
            idle_from = None
        outflow = flux * (plume.flux[:-1] + plume.entrainment)
        active = outflow > 0.0
        limit = float(np.min(0.5 * m[active] / outflow[active])) if np.any(active) else dt
        h = min(dt - elapsed, limit)
        if substep is not None:
            h = min(h, window_end - elapsed)
        cape_time += plume.cape * h
        active_time += h
        if window < windows:
            window_mass[window] += flux * h
        step = apply_mass_flux(
            plume,
            cloud_base_flux=flux,
            dt=h,
            temperature=T,
            vapour=w,
            evaporating_fraction=evaporating_fraction,
            detrained_state=detrained_state,
            precipitation_evaporation=precipitation_evaporation,
            **fixed,
        )
        if classes is not None and milestones is not None:
            classes, crossed = provenance_substep(
                plume, step, scale=flux * h, classes=classes, layer_mass=m
            )
            milestones.add(crossed)
        T = T + (step.energy + latent_heat * step.condensation) / C
        w = w + step.transport - step.condensation
        for key in keys:
            totals[key] += getattr(step, key)
        moved += flux * h
        elapsed += h
        substeps += 1
        if substep is not None and elapsed >= window_end * (1.0 - 1.0e-12):
            elapsed = window_end
            window_end = min(dt, window_end + substep)
            window, fresh = window + 1, True
        if substeps > MAX_SUBSTEPS:
            raise ConfigurationError(
                f"mass-flux convection needed more than {MAX_SUBSTEPS} substeps in one step "
                f"(cloud-base flux {flux:.3g} kg m^-2 s^-1, CAPE {plume.cape:.3g} J kg^-1, "
                f"removal {removal:.3g} J m^2 kg^-2); refused, not limited"
            )
    if first is None:
        raise ConfigurationError("mass-flux convection was asked to integrate no time")
    plume0, removal0, flux0 = first
    cape_after = entraining_plume(
        temperature=T,
        vapour=w,
        entrainment_rate=entrainment_rate,
        condensate_pathway=condensate_pathway,
        conversion_rate=conversion_rate,
        **fixed,
    ).cape
    result = MassFluxStep(
        energy=totals["energy"],
        transport=totals["transport"],
        condensation=totals["condensation"],
        evaporation=totals["evaporation"],
        formed=totals["formed"],
        surface_precipitation=float(np.sum(totals["formed"]) - np.sum(totals["evaporation"])),
        substeps=substeps,
        **{k: totals[k] for k in _PARTS},
    )
    return Convection(
        step=result,
        plume=plume0,
        cape=plume0.cape,
        cape_removal=removal0,
        first_cloud_base_flux=flux0,
        mean_cloud_base_flux=moved / dt,
        provenance=classes,
        milestones=milestones,
        cape_after=cape_after,
        # The idle remainder of the step holds the CAPE convection left.
        cape_time_mean=(
            cape_time + max(dt - (elapsed if idle_from is None else idle_from), 0.0) * cape_after
        )
        / dt,
        # Windows after convection stopped start, and end, where it left CAPE.
        cape_window_start=np.where(np.isnan(window_start), cape_after, window_start),
        active_time_fraction=active_time / dt,
        window_evaluated=window_evaluated,
        window_status=window_status,
        window_mass=window_mass,
    )
