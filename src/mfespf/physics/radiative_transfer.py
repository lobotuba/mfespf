"""Spectral radiative transfer through a layered column.

Charter 9.3 and section 19 stage 2; requirements RAD-001, RAD-002, RAD-003,
SPC-001, SPC-002, SPC-007, SPC-009.

One solve handles the whole spectrum.  There is no shortwave/longwave split in
the physics -- there is radiation in band ``b``, with two sources: what arrives
at the top from the Sun, and what the layers and surface emit thermally.  A
"greenhouse" configuration is then simply an absorption coefficient that is
small where the solar spectrum lives and large where the terrestrial spectrum
does, which is what the real one is. The shortwave and longwave *diagnostics*
RAD-002 and CLD-001 ask for are recovered by partitioning bands at a declared
cutoff, after the solve rather than before it.

Scattering is not represented. Absorption, emission and surface reflection are.
The angular integral over the hemisphere is replaced by a single slant path
through the declared diffusivity factor, which is a parameterization and is
registered as one.

The fluxes are computed on interfaces and the layer heating is their divergence,
so the column's energy budget closes **by construction**: whatever leaves one
interface and does not arrive at the next is exactly what heated the layer
between them. That is a deliberate structural choice -- it means a conservation
failure in a run points at the coupling or at another module, never at
bookkeeping inside the solver.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

from mfespf.contracts.cost import CostDeclaration, ResourceClass
from mfespf.contracts.module import (
    FidelityLevel,
    Module,
    ModuleContract,
    ModuleStatus,
    UncertaintyReport,
)
from mfespf.contracts.spectral import SpectralGrid
from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError
from mfespf.core.quantities import (
    Field,
    Provenance,
    QuantityKind,
    SignConvention,
    SpatialSupport,
    TimeSemantics,
    Uncertainty,
)
from mfespf.core.versioning import (
    InterfaceContract,
    RequiredInterface,
    Version,
    VersionRange,
)
from mfespf.physics.atmosphere import ATMOSPHERE_TEMPERATURE_STATE
from mfespf.physics.equations import EQUATION_IDS
from mfespf.physics.radiation import (
    band_exitance,
    band_exitance_capture,
    band_exitance_profile,
)
from mfespf.registries.states import TendencyContribution

if TYPE_CHECKING:
    from mfespf.physics.rrtmg_longwave import RRTMGLongwave
    from mfespf.physics.water_vapour_longwave import ChouLongwave

__all__ = [
    "RADIATION_INTERFACE",
    "ColumnRadiationSolution",
    "SpectralRadiationModule",
    "grey_absorption",
    "solve_column_radiation",
    "two_stream_layer",
    "window_absorption",
]

V1 = Version(1, 0, 0)
_DEFAULT_SW_LW_CUTOFF_M = 4.0e-6

RADIATION_INTERFACE = InterfaceContract(
    interface_id="spectral_radiation",
    version=V1,
    physical_semantics={
        "flux_units": "W m^-2",
        "spectral_flux_units": "W m^-2 m^-1",
        "interface_ordering": "index_0_is_the_surface",
        "layer_ordering": "index_0_is_the_lowest_layer",
        "heating_units": "W m^-2 per layer",
        "sign_convention": "separate_upward_and_downward_components",
        "scattering": "not_represented",
    },
    optional_capabilities=frozenset({"spectral_flux_profiles", "shortwave_longwave_partition"}),
)


def grey_absorption(grid: SpectralGrid, kappa: float) -> np.ndarray:
    """A mass absorption coefficient identical in every band, in m^2 kg^-1.

    The null case is ``kappa = 0``: a transparent atmosphere, which must
    reproduce the airless result exactly.
    """
    if kappa < 0.0:
        raise ConfigurationError("a mass absorption coefficient cannot be negative")
    return np.full(grid.n_bins, float(kappa))


def window_absorption(
    grid: SpectralGrid,
    *,
    kappa_opaque: float,
    kappa_window: float = 0.0,
    window: tuple[float, float] = (8.0e-6, 12.0e-6),
    transparent_below: float = _DEFAULT_SW_LW_CUTOFF_M,
) -> np.ndarray:
    """Absorption that is transparent to sunlight and opaque outside a window.

    A deliberately crude caricature of the real spectral structure: nothing
    absorbs below ``transparent_below`` (so sunlight reaches the surface), the
    8-12 micron window is comparatively open, and the rest of the terrestrial
    spectrum is opaque. It exists so that SPC-009 can be tested with something
    other than a grey coefficient -- two configurations carrying equal integrated
    power but different spectra must be allowed to respond differently.
    """
    if kappa_opaque < 0.0 or kappa_window < 0.0:
        raise ConfigurationError("absorption coefficients cannot be negative")
    if not window[0] < window[1]:
        raise ConfigurationError(f"window {window} must be an increasing interval")

    lower, upper = grid.edges[:-1], grid.edges[1:]
    widths = grid.widths

    def overlap_fraction(low: float, high: float) -> np.ndarray:
        span = np.clip(np.minimum(upper, high) - np.maximum(lower, low), 0.0, None)
        return span / widths

    # Each band gets the mass-weighted mean coefficient over its own wavelength
    # range rather than the value at its centre. Sampling at the centre makes the
    # coefficient jump as band edges move across the window's edges, which shows
    # up as a sawtooth in any spectral convergence study -- an artifact of how the
    # coefficient was built, not of the transfer solver being tested.
    transparent = overlap_fraction(grid.edges[0], transparent_below)
    in_window = overlap_fraction(*window)
    opaque = np.clip(1.0 - transparent - in_window, 0.0, None)
    return opaque * float(kappa_opaque) + in_window * float(kappa_window)


def two_stream_layer(
    optical_depth: np.ndarray,
    single_scattering_albedo: np.ndarray,
    asymmetry: np.ndarray,
    diffusivity: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reflectance, transmittance and emissivity of a homogeneous layer.

    The standard two-stream solution for a slab that absorbs, scatters and emits.
    Its two limits are what make it safe to introduce underneath work that was
    validated without it:

    * with no scattering it reduces to ``t = exp(-D tau)``, ``r = 0`` and
      ``e = 1 - t`` -- algebraically identical to the absorption-only sweeps this
      replaces, so every Stage 2 and Stage 3 result is unchanged;
    * with conservative scattering it gives ``r + t = 1`` and ``e = 0``, so a
      purely scattering layer neither absorbs nor emits.

    Scattering is needed for clouds specifically. A cloud's shortwave effect is
    reflection, and representing it as absorption would warm the column where the
    real one cools it -- a wrong response sign, which the charter's failure
    taxonomy (15.4) treats as a refutation rather than an inaccuracy.

    Emissivity follows from ``e = 1 - r - t`` rather than being supplied
    independently: that is Kirchhoff's law plus energy conservation, and deriving
    it means the layer cannot emit energy it did not absorb.
    """
    tau = np.asarray(optical_depth, dtype=float)
    # Not clipped below one: a conservatively scattering layer must come out
    # with r + t exactly one and e exactly zero, and a guard of 1 - 1e-12 leaves
    # a few picowatts of absorption that no physics put there.
    omega = np.clip(np.asarray(single_scattering_albedo, dtype=float), 0.0, 1.0)
    g = np.clip(np.asarray(asymmetry, dtype=float), 0.0, 1.0 - 1e-12)

    # Delta scaling removes the forward-scattering peak, which a two-stream
    # cannot represent, and rescales the remainder so the transmitted energy is
    # right. Without it a cloud with g ~ 0.85 would reflect far too much.
    forward = g**2
    tau = tau * (1.0 - omega * forward)
    omega = omega * (1.0 - forward) / np.maximum(1.0 - omega * forward, 1e-300)
    g = g / (1.0 + g)

    backscatter = 0.5 * (1.0 - g)
    gamma1 = diffusivity * (1.0 - omega * (1.0 - backscatter))
    gamma2 = diffusivity * omega * backscatter
    k = np.sqrt(np.maximum(gamma1**2 - gamma2**2, 0.0))

    x = k * tau
    reflectance = np.zeros_like(tau)
    transmittance = np.ones_like(tau)

    # Two regimes, for accuracy rather than taste: the series form is the k -> 0
    # limit where the exponential form is 0/0, and the scaled exponential form
    # avoids the overflow that cosh and sinh reach around tau ~ 400.
    small = x < 1e-6
    if np.any(small):
        denominator = 1.0 + gamma1[small] * tau[small]
        reflectance[small] = gamma2[small] * tau[small] / denominator
        transmittance[small] = 1.0 / denominator
    large = ~small
    if np.any(large):
        decay = np.exp(-2.0 * x[large])
        denominator = (k[large] + gamma1[large]) + (k[large] - gamma1[large]) * decay
        reflectance[large] = gamma2[large] * (1.0 - decay) / denominator
        transmittance[large] = 2.0 * k[large] * np.exp(-x[large]) / denominator

    emissivity = np.clip(1.0 - reflectance - transmittance, 0.0, 1.0)
    return reflectance, transmittance, emissivity


@dataclass(frozen=True)
class PlanckFractions:
    """Correlated-k emission (longwave-source-function.md (p)).

    A g-point is the set of wavenumbers across its band whose absorption falls
    in one range, not a wavelength interval, so it emits its Planck fraction of
    the band's Planck flux (Lacis & Oinas 1991; Mlawer et al. 1997). ``band``
    names each bin's band, -1 for a bin that emits its own Planck flux;
    ``layer`` is each layer's fraction per bin and ``surface`` the surface's.
    """

    band: np.ndarray  # (n_bins,)
    layer: np.ndarray  # (n_layers, n_bins)
    surface: np.ndarray  # (n_bins,)

    def weighted(self, own: np.ndarray, fractions: np.ndarray) -> np.ndarray:
        """Each g-point bin's fraction of its band's total; other bins unchanged.

        ``own`` is the Planck flux per bin, shaped (rows, bins), and
        ``fractions`` the fractions for the same rows.
        """
        own = np.atleast_2d(own)
        grouped = self.band >= 0
        if not np.any(grouped):
            return own
        member = self.band[None, :] == np.arange(int(self.band.max()) + 1)[:, None]
        totals = own @ member.T.astype(float)
        out = own.copy()
        out[:, grouped] = np.atleast_2d(fractions)[:, grouped] * totals[:, self.band[grouped]]
        return out


def linear_in_tau_weight(x: np.ndarray) -> np.ndarray:
    """RRTMG's linear-in-tau weight f(x) (Clough et al. 1992; ``rrtmg_lw_init``).

    A layer of emitting depth ``x`` emits ``B_lay + f(x) (B_edge - B_lay)`` toward
    that edge: f -> 0 when transparent, so it emits at its mean temperature, and
    f -> 1 when opaque, so it emits at the edge's.
    """
    x = np.asarray(x, dtype=float)
    out = x / 6.0
    thick = x > 0.06
    xt = x[thick]
    out[thick] = 1.0 - 2.0 * (1.0 / xt - np.exp(-xt) / -np.expm1(-xt))
    return out


def interface_temperatures(
    layer_temperature: np.ndarray,
    layer_pressure: np.ndarray,
    interface_pressure: np.ndarray,
    surface_temperature: float,
) -> np.ndarray:
    """Interface temperatures by CAM's RRTMG convention (``rrtmg_state.F90``).

    Interior interfaces lie linearly in ln p between the adjacent layer
    temperatures; the top interface takes the top layer's temperature; the
    bottom is CAM's (LW_up / sigma)^(1/4), the surface temperature for a black
    surface. Indexed from the surface (0) to the model top.
    """
    t = np.asarray(layer_temperature, dtype=float)
    lnp = np.log(np.asarray(layer_pressure, dtype=float))
    lnpi = np.log(np.asarray(interface_pressure, dtype=float))
    if lnpi.size != t.size + 1 or lnp.size != t.size:
        raise ConfigurationError("interface temperatures need n layers and n + 1 interfaces")
    out = np.empty(t.size + 1)
    out[0] = float(surface_temperature)
    out[-1] = t[-1]
    dy = (lnpi[1:-1] - lnp[:-1]) / (lnp[1:] - lnp[:-1])
    out[1:-1] = t[:-1] + dy * (t[1:] - t[:-1])
    return out


@dataclass
class ColumnRadiationSolution:
    """Interface fluxes, layer heating and the diagnostics RAD-002 requires.

    Interfaces are indexed from the surface (0) to the model top (``n_layers``);
    layer ``i`` lies between interfaces ``i`` and ``i+1``.
    """

    upward: np.ndarray  # (n_interfaces, n_bands)
    downward: np.ndarray  # (n_interfaces, n_bands)
    layer_heating: np.ndarray  # (n_layers,) W m^-2
    surface_absorbed: float
    surface_emitted: float
    grid: SpectralGrid
    shortwave_mask: np.ndarray
    # The optical properties the solve was actually performed with, carried out
    # rather than recomputed. Three registered equations declare ranges on the
    # optical depth and one more on the single-scattering albedo, and none of
    # them could be checked while these stayed local to the solver. Recomputing
    # them outside would check a reconstruction rather than the thing used.
    # None on a mixed sky, where there is no single set: see ``_mix_columns``.
    optical_depth: np.ndarray | None  # (n_layers, n_bands), extinction
    single_scattering_albedo: np.ndarray | None  # (n_layers, n_bands)
    layer_emissivity: np.ndarray | None  # (n_layers, n_bands)

    @property
    def n_layers(self) -> int:
        return int(self.layer_heating.size)

    def interface_up(self, index: int) -> float:
        return float(np.sum(self.upward[index]))

    def interface_down(self, index: int) -> float:
        return float(np.sum(self.downward[index]))

    def partition(self, values: np.ndarray) -> tuple[float, float]:
        """Split a per-band flux into (shortwave, longwave) at the declared cutoff."""
        return (
            float(np.sum(values[self.shortwave_mask])),
            float(np.sum(values[~self.shortwave_mask])),
        )

    @property
    def toa_index(self) -> int:
        return self.n_layers

    def energy_residual(self) -> float:
        """Net radiative energy the column gains minus what the layers absorbed.

        Zero to round-off by construction; it is checked anyway, because "by
        construction" is a claim about the code and this is the measurement of it.
        """
        net_in = (
            self.interface_down(self.toa_index)
            - self.interface_up(self.toa_index)
            + self.interface_up(0)
            - self.interface_down(0)
        )
        return float(np.sum(self.layer_heating)) - net_in


def _mix_many(
    clear: ColumnRadiationSolution, cloudy: list[tuple[float, ColumnRadiationSolution]]
) -> ColumnRadiationSolution:
    """The independent-column mean of a clear column and disjoint cloudy ones.

    ``_mix_columns`` generalised to several cloud types that partition the sky.
    Linear in every flux for the same reason, so the mix closes its budget.
    """
    total = sum(amount for amount, _ in cloudy)
    if any(a < 0.0 for a, _ in cloudy) or total > 1.0 + 1e-12:
        raise ConfigurationError(f"cloud amounts {[a for a, _ in cloudy]} must partition the sky")
    parts = [(1.0 - total, clear), *cloudy]

    def mix(name: str) -> Any:
        return sum(w * getattr(s, name) for w, s in parts)

    return ColumnRadiationSolution(
        upward=mix("upward"),
        downward=mix("downward"),
        layer_heating=mix("layer_heating"),
        surface_absorbed=float(mix("surface_absorbed")),
        surface_emitted=float(mix("surface_emitted")),
        grid=clear.grid,
        shortwave_mask=clear.shortwave_mask,
        optical_depth=None,
        single_scattering_albedo=None,
        layer_emissivity=None,
    )


def combine_scatterers(
    *contributions: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Add two or more scatterers occupying the same layers.

    Each contribution is ``(absorption, scattering, asymmetry)`` on the same
    ``(n_layers, n_bands)`` grid. Optical depths add; the asymmetry does not.
    It is a property of the scattered light, so the combination is weighted by
    how much each contributor scatters -- a thin, strongly forward-scattering
    aerosol beneath a thick cloud must not drag the layer's asymmetry towards
    its own value.

    Where nothing scatters the weighted mean is undefined, and zero is returned
    rather than a division by zero: with no scattering the asymmetry multiplies
    nothing and any value would do, so the one that cannot mislead is best.
    """
    if not contributions:
        raise ConfigurationError("combine_scatterers needs at least one contribution")
    absorption = np.sum([np.asarray(a, dtype=float) for a, _s, _g in contributions], axis=0)
    scattering = np.sum([np.asarray(s, dtype=float) for _a, s, _g in contributions], axis=0)
    weighted = np.sum(
        [np.asarray(s, dtype=float) * np.asarray(g, dtype=float) for _a, s, g in contributions],
        axis=0,
    )
    asymmetry = np.divide(weighted, scattering, out=np.zeros_like(weighted), where=scattering > 0.0)
    return absorption, scattering, asymmetry


def solve_column_radiation(
    *,
    grid: SpectralGrid,
    incident_spectral_irradiance: np.ndarray,
    layer_temperature: np.ndarray,
    layer_mass: np.ndarray,
    absorption_coefficient: np.ndarray,
    surface_temperature: float,
    surface_albedo: np.ndarray,
    absorption_optical_depth: np.ndarray | None = None,
    scattering_optical_depth: np.ndarray | None = None,
    asymmetry: np.ndarray | None = None,
    diffusivity_factor: float | None = None,
    shortwave_cutoff_m: float = _DEFAULT_SW_LW_CUTOFF_M,
    planck_fractions: PlanckFractions | None = None,
    interface_temperature: np.ndarray | None = None,
) -> ColumnRadiationSolution:
    """Solve absorption, scattering, emission and surface reflection in the column.

    Solved by the adding method: sweeping upward, accumulate the reflectance and
    the emitted source of everything below each interface, then sweep down once
    with the incident field. That handles multiple reflection between layers
    exactly, which a pair of independent sweeps cannot once anything scatters.

    The optical-property arrays are optional and shaped ``(n_layers, n_bands)``.
    Omitting them leaves the column non-scattering, in which case the result is
    the absorption-only solution to floating-point round-off.

    They describe *a layer*, not a cloud. Cloud droplets and volcanic aerosol
    scatter in the same layers and their contributions add, so the caller
    combines them with ``combine_scatterers`` before calling here -- the solver
    has no business knowing which scatterer is which.

    Two source-function options (longwave-source-function.md), both off by
    default: ``planck_fractions`` makes g-point bins emit their fraction of the
    band's Planck flux, and ``interface_temperature`` (n_layers + 1, from the
    surface) lets each layer's source vary linearly in optical depth toward the
    edge it emits from.
    """
    n_layers = int(np.asarray(layer_mass).size)
    n_bands = grid.n_bins
    if np.asarray(layer_temperature).size != n_layers:
        raise ConfigurationError(
            f"layer_temperature has {np.asarray(layer_temperature).size} entries for "
            f"{n_layers} layers"
        )
    if np.asarray(absorption_coefficient).size != n_bands:
        raise ConfigurationError(
            f"absorption_coefficient has {np.asarray(absorption_coefficient).size} entries "
            f"for {n_bands} bands"
        )
    if np.asarray(surface_albedo).size != n_bands:
        raise ConfigurationError(
            f"surface_albedo has {np.asarray(surface_albedo).size} entries for {n_bands} bands"
        )

    D = (
        diffusivity_factor
        if diffusivity_factor is not None
        else value_of("radiative_diffusivity_factor")
    )
    kappa = np.asarray(absorption_coefficient, dtype=float)
    mass = np.asarray(layer_mass, dtype=float)
    temperature = np.asarray(layer_temperature, dtype=float)
    albedo = np.clip(np.asarray(surface_albedo, dtype=float), 0.0, 1.0)

    shape = (n_layers, n_bands)

    def _optional(values: np.ndarray | None) -> np.ndarray:
        if values is None:
            return np.zeros(shape)
        array = np.asarray(values, dtype=float)
        if array.shape != shape:
            raise ConfigurationError(
                f"cloud optical property has shape {array.shape}, expected {shape}"
            )
        if np.any(array < 0.0):
            # Every absorber and scatterer arrives here, gases included, so the
            # refusal names the quantity rather than guessing its source.
            raise ConfigurationError("absorption and scattering optical depths cannot be negative")
        return array

    cloud_absorption = _optional(absorption_optical_depth)
    cloud_scattering = _optional(scattering_optical_depth)
    asymmetry = (
        np.zeros(shape)
        if asymmetry is None
        else np.broadcast_to(np.asarray(asymmetry, dtype=float), shape)
    )

    gas_absorption = kappa[None, :] * mass[:, None]
    total_optical_depth = gas_absorption + cloud_absorption + cloud_scattering
    single_scattering_albedo = np.divide(
        cloud_scattering,
        total_optical_depth,
        out=np.zeros(shape),
        where=total_optical_depth > 0.0,
    )

    reflectance, transmittance, emissivity = two_stream_layer(
        total_optical_depth, single_scattering_albedo, asymmetry, D
    )
    layer_emission = band_exitance_profile(grid, temperature)
    surface_emission = band_exitance(grid, float(surface_temperature))
    if planck_fractions is not None:
        fractions = planck_fractions
        if (
            fractions.band.shape != (n_bands,)
            or fractions.layer.shape != shape
            or fractions.surface.shape != (n_bands,)
        ):
            raise ConfigurationError("Planck fractions must match the layers and the grid")
        grouped = fractions.band >= 0
        closure = np.bincount(
            fractions.band[grouped], weights=fractions.surface[grouped], minlength=1
        )
        present = np.bincount(fractions.band[grouped], minlength=1) > 0
        if np.any(np.abs(closure[present] - 1.0) > 1.0e-6):
            # The surface module emits its own greybody total; the radiation's
            # surface emission carries the same total only if the fractions close.
            raise ConfigurationError("each band's surface Planck fractions must sum to one")
        layer_emission = fractions.weighted(layer_emission, fractions.layer)
        surface_emission = fractions.weighted(surface_emission, fractions.surface)[0]

    if interface_temperature is None:
        source_up = emissivity * layer_emission
        source_down = source_up
    else:
        edges = np.asarray(interface_temperature, dtype=float)
        if edges.shape != (n_layers + 1,):
            raise ConfigurationError(
                f"interface_temperature has {edges.size} entries for {n_layers + 1} interfaces"
            )
        edge_emission = band_exitance_profile(grid, edges)
        top, bottom = edge_emission[1:], edge_emission[:-1]
        if planck_fractions is not None:
            # Each layer weights both its edges by its own fractions, as rtrnmr does.
            top = planck_fractions.weighted(top, planck_fractions.layer)
            bottom = planck_fractions.weighted(bottom, planck_fractions.layer)
        # Emission comes from absorption, so the emitting depth is the
        # absorption optical depth along the diffuse path, not the extinction.
        weight = linear_in_tau_weight(D * (gas_absorption + cloud_absorption))
        source_up = emissivity * (layer_emission + weight * (top - layer_emission))
        source_down = emissivity * (layer_emission + weight * (bottom - layer_emission))

    # Upward pass: reflectance and emitted source of everything below interface i.
    below_reflectance = np.zeros((n_layers + 1, n_bands))
    below_source = np.zeros((n_layers + 1, n_bands))
    surface_emissivity = 1.0 - albedo
    below_reflectance[0] = albedo
    below_source[0] = surface_emissivity * surface_emission

    for i in range(n_layers):
        coupling = 1.0 - reflectance[i] * below_reflectance[i]
        below_reflectance[i + 1] = (
            reflectance[i] + (transmittance[i] ** 2 * below_reflectance[i]) / coupling
        )
        # What leaves the layer's bottom is reflected back up from below; what
        # leaves its top goes straight up.
        below_source[i + 1] = (
            transmittance[i]
            * below_reflectance[i]
            * (reflectance[i] * below_source[i] + source_down[i])
            / coupling
            + transmittance[i] * below_source[i]
            + source_up[i]
        )

    # Downward pass with the incident field, then the upward field follows.
    downward = np.zeros((n_layers + 1, n_bands))
    downward[n_layers] = np.asarray(incident_spectral_irradiance, dtype=float) * grid.widths
    for i in range(n_layers - 1, -1, -1):
        coupling = 1.0 - reflectance[i] * below_reflectance[i]
        downward[i] = (
            transmittance[i] * downward[i + 1] + reflectance[i] * below_source[i] + source_down[i]
        ) / coupling

    upward = below_reflectance * downward + below_source

    # Heating is still the flux divergence across the layer, so the column budget
    # closes identically whatever the layer coefficients turned out to be.
    layer_heating = np.array(
        [
            float(np.sum((upward[i] - upward[i + 1]) + (downward[i + 1] - downward[i])))
            for i in range(n_layers)
        ]
    )

    shortwave_mask = grid.centers < shortwave_cutoff_m
    return ColumnRadiationSolution(
        upward=upward,
        downward=downward,
        layer_heating=layer_heating,
        surface_absorbed=float(np.sum(surface_emissivity * downward[0])),
        surface_emitted=float(np.sum(surface_emissivity * surface_emission)),
        grid=grid,
        shortwave_mask=shortwave_mask,
        optical_depth=total_optical_depth,
        single_scattering_albedo=single_scattering_albedo,
        layer_emissivity=emissivity,
    )


def _extent_over(solutions: list[ColumnRadiationSolution], attribute: str) -> np.ndarray:
    """The smallest and largest value an optical property took across the skies.

    **Why the extent and not the field.** Three equations declare ranges on the
    optical depth and one on the single-scattering albedo, and a declared range
    is a statement about the extent a relation was evaluated over -- which is
    two numbers, whatever the resolution behind them. Publishing the whole
    (sky, layer, band) array instead was measured at 80 MB across the 25,916
    steps of the historical replay, enough to exhaust memory on a run that was
    already near the limit, in order to carry information that is consumed as a
    minimum and a maximum. The full arrays remain on each
    ``ColumnRadiationSolution`` for the current step, where anyone who wants to
    know *which* layer is thick can still ask.
    """
    values = [np.asarray(getattr(s, attribute), dtype=float) for s in solutions]
    stacked = np.concatenate([v.reshape(-1) for v in values])
    return np.array([float(stacked.min()), float(stacked.max())])


def _mix_columns(
    clear: ColumnRadiationSolution, overcast: ColumnRadiationSolution, fraction: float
) -> ColumnRadiationSolution:
    """Combine a clear and an overcast column under maximum overlap.

    The independent-column approximation: solve the two skies separately and
    weight them. Every quantity mixes linearly because the radiative transfer is
    linear in the source terms, so the mixed fluxes remain a consistent solution
    of the same energy budget -- which is what keeps the flux divergence closing
    identically after the mix.
    """
    if not 0.0 <= fraction <= 1.0:
        raise ConfigurationError(f"cloud fraction {fraction} must lie in [0, 1]")
    weight = fraction
    return ColumnRadiationSolution(
        upward=(1.0 - weight) * clear.upward + weight * overcast.upward,
        downward=(1.0 - weight) * clear.downward + weight * overcast.downward,
        layer_heating=(1.0 - weight) * clear.layer_heating + weight * overcast.layer_heating,
        surface_absorbed=(1.0 - weight) * clear.surface_absorbed
        + weight * overcast.surface_absorbed,
        surface_emitted=(1.0 - weight) * clear.surface_emitted + weight * overcast.surface_emitted,
        grid=clear.grid,
        shortwave_mask=clear.shortwave_mask,
        # The fluxes mix; the optical properties do not. A fraction-weighted
        # optical depth is not the optical depth of any column that was solved,
        # and the transfer is linear in the *source* terms, not in tau. Carrying
        # the weighted mean would publish a number that looks checkable and
        # describes nothing, so the mixed sky declines to claim one and the
        # module reports the two skies it actually solved instead.
        optical_depth=None,
        single_scattering_albedo=None,
        layer_emissivity=None,
    )


class SpectralRadiationModule(Module):
    """Solves the column radiation field and hands out labeled heating tendencies.

    The module owns no state.  It reads temperatures and optical properties from
    their owners, solves, publishes fluxes, and submits a per-layer tendency the
    atmosphere commits (MOD-011). It is ACTIVE rather than DIAGNOSTIC because it
    participates in the evolution -- a diagnostic module by charter 6.4 does not.
    """

    def __init__(
        self,
        *,
        module_id: str = "radiation",
        atmosphere_module_id: str = "atmosphere",
        grid: SpectralGrid,
        absorption_coefficient: np.ndarray,
        diffusivity_factor: float | None = None,
        shortwave_cutoff_m: float = _DEFAULT_SW_LW_CUTOFF_M,
        couple_clouds: bool = False,
        cloud_module_id: str = "cloud",
        couple_greenhouse: bool = False,
        greenhouse_module_id: str = "greenhouse",
        couple_volcanic: bool = False,
        volcanic_module_id: str = "volcanic",
        couple_water: bool = False,
        water_module_id: str = "water",
        stratospheric_layers: int = 1,
        greenhouse_background: float | None = None,
        solar_vapour_k: np.ndarray | None = None,
        ozone_k: np.ndarray | None = None,
        co2_band: tuple[np.ndarray, np.ndarray, float] | None = None,
        ozone_infrared: np.ndarray | None = None,
        water_longwave: ChouLongwave | None = None,
        rrtmg: RRTMGLongwave | None = None,
        k_distribution_source: bool = False,
        linear_in_tau_source: bool = False,
        couple_cloud_field: bool = False,
        cloud_field_module_id: str = "cloud_field",
        status: ModuleStatus = ModuleStatus.ACTIVE,
    ) -> None:
        n_bands = grid.n_bins
        solar_vapour = solar_vapour_k is not None
        ozone = ozone_k is not None
        co2 = co2_band is not None
        ozone_ir = ozone_infrared is not None
        chou_vapour = water_longwave is not None and couple_water
        use_rrtmg = rrtmg is not None
        contract = ModuleContract(
            module_id=module_id,
            module_name=f"spectral column radiative transfer ({n_bands} bands)",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.REDUCED,
            physical_basis=(
                "Plane-parallel two-stream radiative transfer with absorption, "
                "scattering, thermal emission under Kirchhoff's law, and surface "
                "reflection, solved by the adding method."
            ),
            governing_equation_ids=(
                EQUATION_IDS["schwarzschild"],
                EQUATION_IDS["beer_lambert"],
                EQUATION_IDS["planck"],
                EQUATION_IDS["diffusivity"],
                EQUATION_IDS["two_stream"],
                *((EQUATION_IDS["solar_vapour_absorption"],) if solar_vapour else ()),
                *((EQUATION_IDS["ozone_solar_absorption"],) if ozone else ()),
                *((EQUATION_IDS["co2_15_micron_band"],) if co2 else ()),
                *((EQUATION_IDS["ozone_9_6_micron_band"],) if ozone_ir else ()),
                *((EQUATION_IDS["water_vapour_longwave"],) if chou_vapour else ()),
                *((EQUATION_IDS["rrtmg_longwave"],) if use_rrtmg else ()),
                *((EQUATION_IDS["correlated_k_emission"],) if k_distribution_source else ()),
                *((EQUATION_IDS["linear_in_tau_source"],) if linear_in_tau_source else ()),
                *(
                    (EQUATION_IDS["observed_cloud_field"], EQUATION_IDS["cloud_radiative_effect"])
                    if couple_cloud_field
                    else ()
                ),
                *(
                    (EQUATION_IDS["cloud_optics"], EQUATION_IDS["cloud_radiative_effect"])
                    if couple_clouds
                    else ()
                ),
            ),
            equation_variables=(
                # RRTMG's tables are tabulated at T_ref(p) and +/-15, +/-30 K; the
                # departure is what its domain bounds (rrtmg-domain-and-refusals.md).
                *(
                    (
                        (
                            EQUATION_IDS["rrtmg_longwave"],
                            (("dT_ref", "rrtmg_reference_temperature_departure", "K"),),
                        ),
                    )
                    if use_rrtmg
                    else ()
                ),
                # The cloud field's optical thickness is the cloud's, published
                # by the cloud field. Inferred, it matched this module's total
                # layer optical depth, saturated vapour bands and all.
                *(
                    (
                        (
                            EQUATION_IDS["observed_cloud_field"],
                            (("tau_c", "cloud_type_visible_optical_thickness", "1"),),
                        ),
                    )
                    if couple_cloud_field
                    else ()
                ),
                (
                    EQUATION_IDS["schwarzschild"],
                    (
                        # The source function of the medium is the layer
                        # temperature. The surface enters this solve as the
                        # lower boundary condition, not as the transfer
                        # equation's T, so binding it here would judge the
                        # wrong quantity against this record's range.
                        ("T", ATMOSPHERE_TEMPERATURE_STATE, "K"),
                    ),
                ),
                *(
                    (
                        (
                            EQUATION_IDS["cloud_optics"],
                            (
                                ("W", "cloud_liquid_water_path", "kg m^-2"),
                                ("W", "cloud_ice_water_path", "kg m^-2"),
                                ("r_eff", "cloud_effective_radius", "m"),
                            ),
                        ),
                    )
                    if couple_clouds
                    else ()
                ),
                (
                    EQUATION_IDS["planck"],
                    (
                        # Planck, unlike Schwarzschild, genuinely is evaluated
                        # at both: once per layer for the emission source and
                        # once at the surface for the boundary. Both are bound
                        # rather than one chosen, because choosing would leave
                        # the other evaluated and unjudged.
                        ("T", ATMOSPHERE_TEMPERATURE_STATE, "K"),
                        ("T", "surface_temperature", "K"),
                        ("lambda", "spectral_edges", "m"),
                    ),
                ),
            ),
            assumptions=(
                "plane-parallel, horizontally homogeneous column",
                "gases absorb and emit but do not scatter; only clouds scatter",
                "local thermodynamic equilibrium, so emissivity equals absorptivity",
                "each layer emits at its single isothermal temperature",
                "the hemispheric angular integral is replaced by one diffusivity factor",
                *(
                    (
                        "maximum overlap: all-sky is a linear mix of a clear and an "
                        "overcast column weighted by the cloud fraction",
                    )
                    if couple_clouds
                    else ()
                ),
            ),
            validity_domain=(
                "Optical depths up to about 10 per layer. The lower end used to be "
                "declared as 0.1, which measurement contradicts: the diffusivity "
                "approximation is most accurate in the thin limit, not least, and its "
                "worst absolute error falls near unit optical depth in the middle of "
                "the range. What survives of the thin-limit concern is a fixed 17 "
                "percent underestimate of a vanishing absorptance, reported as the "
                "diffusivity_approximation uncertainty rather than as a bound. "
                "Strongly inhomogeneous layers are outside this domain."
            ),
            required_inputs=(
                "toa_spectral_irradiance",
                ATMOSPHERE_TEMPERATURE_STATE,
                "atmosphere_layer_mass",
                "atmosphere_layer_heat_capacity",
                "surface_temperature",
                "surface_spectral_albedo",
                # All four, not just the fraction. This module derives the cloud
                # optical properties itself from the water paths and the radii,
                # so it consumes them as surely as it consumes the fraction.
                # Declaring only the fraction left the other three undeclared
                # inputs: the coupler could not check the wiring, and a missing
                # one surfaced as a run-time error from inside the solve rather
                # than as an unsatisfied contract before the run started.
                *(
                    (
                        "cloud_fraction",
                        "cloud_liquid_water_path",
                        "cloud_ice_water_path",
                        "cloud_effective_radius",
                    )
                    if couple_clouds
                    else ()
                ),
                *(("co2_mole_fraction",) if couple_greenhouse else ()),
                *(("stratospheric_aerosol_optical_depth",) if couple_volcanic else ()),
                *(("atmosphere_water_vapour",) if couple_water else ()),
                # The solar k-distributions scale the path by pressure, so the
                # pressure is consumed, and declared, only when they are on.
                *(
                    ("atmosphere_layer_pressure",)
                    if (solar_vapour and couple_water)
                    or ozone
                    or co2
                    or ozone_ir
                    or chou_vapour
                    or use_rrtmg
                    or linear_in_tau_source
                    else ()
                ),
                # Ozone is prescribed on pressure, and the column above the
                # model top is placed in the top layer from its top interface.
                # The linear source's interface temperatures lie in ln p between the layers.
                *(
                    ("atmosphere_interface_pressure",)
                    if ozone or ozone_ir or use_rrtmg or linear_in_tau_source
                    else ()
                ),
                *(
                    (
                        "cloud_type_amount",
                        "cloud_type_top_pressure",
                        "cloud_type_visible_optical_thickness",
                        "cloud_type_infrared_optical_thickness",
                        "atmosphere_interface_pressure",
                    )
                    if couple_cloud_field
                    else ()
                ),
            ),
            produced_outputs=(
                *(
                    (
                        "rrtmg_reference_temperature_departure",
                        "rrtmg_negative_optical_depths_floored",
                        "rrtmg_most_negative_optical_depth",
                    )
                    if use_rrtmg
                    else ()
                ),
                "layer_radiative_heating",
                "layer_optical_depth_extent",
                "layer_single_scattering_albedo_extent",
                "layer_emissivity_extent",
                "toa_downward_radiation",
                "toa_upward_radiation",
                "surface_downward_radiation",
                "surface_upward_radiation",
                "surface_absorbed_radiation",
                "outgoing_longwave_radiation",
                "reflected_shortwave_radiation",
                "atmosphere_absorbed_shortwave",
                "surface_absorbed_shortwave",
                "atmospheric_radiative_forcing_diagnostics",
                *(
                    (
                        "clear_sky_outgoing_longwave_radiation",
                        "clear_sky_reflected_shortwave_radiation",
                        "cloud_shortwave_radiative_effect",
                        "cloud_longwave_radiative_effect",
                        "cloud_net_radiative_effect",
                    )
                    if couple_clouds or couple_cloud_field
                    else ()
                ),
            ),
            conserved_quantities=("energy",),
            supported_spectral_representations=(grid.resolution_identifier,),
            provided_interfaces=(RADIATION_INTERFACE,),
            required_interfaces=(
                RequiredInterface(
                    interface_id="atmospheric_state",
                    version_range=VersionRange.parse(">=1.0.0,<2.0.0"),
                    expects={
                        "layer_ordering": "index_0_is_the_lowest_layer",
                        "vertical_coordinate": "pressure",
                        "temperature_units": "K",
                        "layer_mass_units": "kg m^-2",
                    },
                ),
                RequiredInterface(
                    interface_id="surface",
                    version_range=VersionRange.parse(">=1.0.0,<2.0.0"),
                    expects={
                        "temperature_units": "K",
                        "flux_sign_convention": "positive_upward_away_from_surface",
                    },
                ),
            ),
            uncertainty_outputs=(
                "spectral_discretization",
                "diffusivity_approximation",
                "absorption_coefficient_parameter",
                "omitted_scattering_structural",
            ),
            correctness_tests=(
                "test_transparent_column_recovers_the_airless_result",
                "test_single_grey_layer_matches_the_analytic_greenhouse",
                "test_radiation_flux_divergence_closes_identically",
                "test_incident_equals_absorbed_plus_reflected_plus_transmitted",
                "test_spectral_convergence_of_the_column_solution",
                "test_two_stream_reduces_to_pure_absorption",
                *(
                    (
                        "test_no_cloud_control_matches_the_clear_column",
                        "test_cloud_shortwave_effect_cools_and_longwave_effect_warms",
                    )
                    if couple_clouds
                    else ()
                ),
            ),
            known_limitations=(
                "Gas scattering is absent, so Rayleigh scattering and its contribution "
                "to planetary albedo are not represented at all.",
                "One diffusivity factor for both solar beam and thermal diffuse paths.",
                *(
                    (
                        "All-sky fluxes are a linear mix of two independent columns, so "
                        "horizontal photon transport between the clear and cloudy parts "
                        "is absent. The cost is two full solves per step.",
                    )
                    if couple_clouds
                    else ()
                ),
            ),
            fidelity_justification=(
                "A grey atmosphere has no window and no band wings, so it cannot "
                "produce a logarithmic forcing in carbon dioxide or a water-vapour "
                "continuum. Measured: a grid too coarse to resolve the wings gives "
                "linear growth or complete saturation instead, which is why "
                "greenhouse_grid refuses one."
            ),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(n_layers * n_bands)",
                memory_complexity="O(n_layers * n_bands)",
                scaling_variables=("n_bands", "n_layers"),
                work_units_per_step=8.0,
                benchmark_configuration=f"{n_bands} bands, global scalar column",
                peak_memory_bytes=16 * 8 * n_bands + 8192,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
                lower_fidelity_alternatives=("grey_broadband_radiation",),
            ),
        )
        super().__init__(contract, status=status)
        self.atmosphere_module_id = atmosphere_module_id
        self.grid = grid
        self.absorption_coefficient = np.asarray(absorption_coefficient, dtype=float)
        self.diffusivity_factor = (
            diffusivity_factor
            if diffusivity_factor is not None
            else value_of("radiative_diffusivity_factor")
        )
        self.shortwave_cutoff_m = shortwave_cutoff_m
        self.couple_clouds = couple_clouds
        self.couple_volcanic = couple_volcanic
        self.couple_water = couple_water
        self.couple_cloud_field = couple_cloud_field
        self.cloud_field_module_id = cloud_field_module_id
        # The sub-columns of an observed cloud field, as (amount, solution),
        # solved every step alongside the clear column.
        self.subcolumns: list[tuple[float, ColumnRadiationSolution]] = []
        # Reference absorption coefficients for sunlight, one per band, from
        # Chou and Suarez (1999); None leaves the atmosphere transparent to it.
        self.solar_vapour_k = (
            np.asarray(solar_vapour_k, dtype=float) if solar_vapour_k is not None else None
        )
        if self.solar_vapour_k is not None and self.solar_vapour_k.size != n_bands:
            raise ConfigurationError(
                f"solar vapour coefficients cover {self.solar_vapour_k.size} bands and the "
                f"grid has {n_bands}; build both with solar_vapour_bands"
            )
        # Ozone's effective coefficients for sunlight, (cm-atm)^-1 per band
        # (Chou and Suarez 1999); None leaves ozone out. Its amount is the
        # SWOOSH profile placed on this column's layers on the first step.
        self.ozone_k = np.asarray(ozone_k, dtype=float) if ozone_k is not None else None
        if self.ozone_k is not None and self.ozone_k.size != n_bands:
            raise ConfigurationError(
                f"ozone coefficients cover {self.ozone_k.size} bands and the grid has "
                f"{n_bands}; build both with ozone_bands"
            )
        self.ozone_path: np.ndarray | None = None
        # Chou et al.'s (2001) 15-micron CO2 band: (vertical k per bin in
        # (cm-atm)^-1, each bin's kind, CO2 ppmv). None leaves CO2 to whatever
        # absorption_coefficient carries.
        self.co2_band = co2_band
        if co2_band is not None and np.asarray(co2_band[0]).size != n_bands:
            raise ConfigurationError(
                f"CO2 band coefficients cover {np.asarray(co2_band[0]).size} bands and the "
                f"grid has {n_bands}; build both with co2_band_grid"
            )
        # Ozone's 9.6-micron band (ozone-infrared-band.md): each bin's term
        # index, -1 outside the band. Its amount is the same SWOOSH placement
        # the solar ozone uses.
        self.ozone_infrared = (
            np.asarray(ozone_infrared, dtype=int) if ozone_infrared is not None else None
        )
        if self.ozone_infrared is not None and self.ozone_infrared.size != n_bands:
            raise ConfigurationError(
                f"ozone band terms cover {self.ozone_infrared.size} bands and the grid has "
                f"{n_bands}; build both with ozone_infrared_grid"
            )
        # Chou et al.'s (2001) water vapour k-distributions on this grid
        # (global-column-chou-longwave.md); None keeps the fitted envelope.
        self.water_longwave = water_longwave
        if water_longwave is not None and water_longwave.water_band.size != n_bands:
            raise ConfigurationError(
                f"water vapour terms cover {water_longwave.water_band.size} bands and the "
                f"grid has {n_bands}; build both with chou_longwave_grid"
            )
        # RRTMG's gas optics (rrtmg-longwave.md): every longwave gas at once.
        self.rrtmg = rrtmg
        if rrtmg is not None:
            if rrtmg.gpoint.size != n_bands:
                raise ConfigurationError(
                    f"RRTMG g-points cover {rrtmg.gpoint.size} bands and the grid has "
                    f"{n_bands}; build both with rrtmg_longwave_grid"
                )
            if not couple_water:
                raise ConfigurationError("RRTMG's longwave needs the radiative vapour")
        # The longwave source function (longwave-source-function.md): (p) each
        # g-point emits its Planck fraction of its band, (s) each layer's source
        # varies linearly in optical depth toward the edge it emits from.
        if k_distribution_source and rrtmg is None:
            raise ConfigurationError(
                "correlated-k emission needs g-points; it applies to RRTMG's longwave"
            )
        self.k_distribution_source = k_distribution_source
        self.linear_in_tau_source = linear_in_tau_source
        self._rrtmg_floor: tuple[int, float] = (0, 0.0)
        self._rrtmg_departure: np.ndarray | None = None
        self.water_module_id = water_module_id
        self.vapour_optical_depth: np.ndarray | None = None
        self.volcanic_optics: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None
        self.volcanic_module_id = volcanic_module_id
        self.stratospheric_layers = stratospheric_layers
        self.volcanic_optical_depth = 0.0
        self.couple_greenhouse = couple_greenhouse
        self.greenhouse_module_id = greenhouse_module_id
        if greenhouse_background is None:
            # The calibrated value, not a literal repeated here. The CO2 band's
            # wing decay is fitted at this opacity, so a column that quietly ran
            # at another one would use a forcing per doubling the calibration
            # never produced -- which is exactly what happened once already.
            from mfespf.physics.greenhouse import INFRARED_BACKGROUND_ABSORPTION

            greenhouse_background = INFRARED_BACKGROUND_ABSORPTION
        self.greenhouse_background = greenhouse_background
        self.co2_mole_fraction = 0.0
        self.cloud_module_id = cloud_module_id
        self.solution: ColumnRadiationSolution | None = None
        self.clear_sky: ColumnRadiationSolution | None = None
        # Kept so the published optical properties cover every sky the
        # diffusivity approximation was applied to, not just the clear one.
        self.overcast: ColumnRadiationSolution | None = None
        self._cloud_fraction = 0.0
        self._layer_heat_capacity: np.ndarray | None = None
        # A property of the grid, not of the run: computed once and reported
        # every step as the structural uncertainty of the spectral truncation.
        self._capture_at_255k = band_exitance_capture(grid, 255.0)

    def validate_configuration(self) -> None:
        super().validate_configuration()
        if self.absorption_coefficient.size != self.grid.n_bins:
            raise ConfigurationError(
                f"absorption coefficient has {self.absorption_coefficient.size} entries for "
                f"{self.grid.n_bins} bands"
            )
        if np.any(self.absorption_coefficient < 0.0):
            raise ConfigurationError("absorption coefficients cannot be negative")
        if self.diffusivity_factor <= 0.0:
            raise ConfigurationError("the diffusivity factor must be positive")

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        self.mark_initialized()

    def compute_tendencies(self) -> list[TendencyContribution]:
        incident = self.input_field("toa_spectral_irradiance").require_units("W m^-2 m^-1")
        temperature = self.input_field(ATMOSPHERE_TEMPERATURE_STATE).require_units("K")
        mass = self.input_field("atmosphere_layer_mass").require_units("kg m^-2")
        surface_temperature = self.input_field("surface_temperature").require_units("K").scalar()
        albedo = self.input_field("surface_spectral_albedo").require_units("1")

        incident_grid = incident.metadata.get("spectral_grid")
        if not isinstance(incident_grid, SpectralGrid):
            raise ConfigurationError(
                "toa_spectral_irradiance arrived without its spectral grid (SPC-004)"
            )
        if incident_grid.n_bins != self.grid.n_bins:
            raise ConfigurationError(
                f"incident radiation is on {incident_grid.n_bins} bands but this module was "
                f"configured for {self.grid.n_bins}; regrid at a declared boundary (SPC-005)"
            )

        # A greenhouse gas enters as a constitutive relation on the absorption
        # coefficient, not as a forcing added to a budget. The transfer solver is
        # unchanged and unaware: it is handed a different kappa and computes the
        # consequence, which is what makes the greenhouse effect a result of this
        # configuration rather than an input to it.
        kappa = self.absorption_coefficient
        if self.couple_greenhouse:
            from mfespf.physics.greenhouse import co2_band_absorption

            self.co2_mole_fraction = (
                self.input_field("co2_mole_fraction").require_units("ppm").scalar()
            )
            from mfespf.physics.greenhouse import co2_wing_decay_for

            # With vapour active the infrared background *is* the vapour, layer
            # by layer, so the flat background is switched off rather than added
            # to it -- and the band's fitted width is the one calibrated against
            # that opacity, which is not the one calibrated against the flat
            # background.
            kappa = co2_band_absorption(
                self.grid,
                self.co2_mole_fraction,
                wing_decay_per_metre=co2_wing_decay_for(active_vapour=self.couple_water),
                background=0.0 if self.couple_water else self.greenhouse_background,
            )

        common = dict(
            grid=self.grid,
            incident_spectral_irradiance=incident.values,
            layer_temperature=temperature.values,
            layer_mass=mass.values,
            absorption_coefficient=kappa,
            surface_temperature=surface_temperature,
            surface_albedo=albedo.values,
            diffusivity_factor=self.diffusivity_factor,
            shortwave_cutoff_m=self.shortwave_cutoff_m,
        )

        if self.linear_in_tau_source:
            common["interface_temperature"] = interface_temperatures(
                temperature.values,
                self.input_field("atmosphere_layer_pressure").require_units("Pa").values,
                self.input_field("atmosphere_interface_pressure").require_units("Pa").values,
                surface_temperature,
            )

        # Water vapour and volcanic aerosol are both in the clear sky as well as
        # under the cloud, so they go into the base optical properties rather
        # than into the cloudy solve. A clear-sky diagnostic that excluded them
        # would attribute their effect to the cloud.
        contributions: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []

        if self.couple_water:
            vapour = self.input_field("atmosphere_water_vapour").require_units("kg m^-2").values
            if self.rrtmg is not None:
                # RRTMG carries the vapour's longwave with every other gas, below;
                # here only the solar part, which is disjoint by band.
                self.vapour_optical_depth = np.zeros((np.asarray(vapour).size, self.grid.n_bins))
            elif self.water_longwave is not None:
                from mfespf.physics.water_vapour_longwave import water_vapour_optical_depth

                # Lines and continuum, scaled for each layer's pressure and
                # temperature every step; zero in every band below 4 microns.
                self.vapour_optical_depth = water_vapour_optical_depth(
                    self.water_longwave,
                    vapour=vapour,
                    layer_mass=mass.values,
                    pressure=self.input_field("atmosphere_layer_pressure")
                    .require_units("Pa")
                    .values,
                    temperature=temperature.values,
                )
            else:
                from mfespf.physics.greenhouse import vapour_absorption_optical_depth

                self.vapour_optical_depth = vapour_absorption_optical_depth(self.grid, vapour)
            if self.solar_vapour_k is not None:
                from mfespf.physics.solar_absorption import solar_vapour_optical_depth

                pressure = self.input_field("atmosphere_layer_pressure").require_units("Pa").values
                # Added, not substituted: the thermal bands and the solar ones
                # are disjoint, so each band gets exactly one of the two.
                self.vapour_optical_depth = self.vapour_optical_depth + (
                    solar_vapour_optical_depth(
                        self.solar_vapour_k, vapour, pressure, temperature.values
                    )
                )
            # Pure absorption: vapour does not scatter at these wavelengths, and
            # representing it as partly scattering would reflect the surface's
            # own emission back down and warm it, which is the wrong sign.
            zeros = np.zeros_like(self.vapour_optical_depth)
            contributions.append((self.vapour_optical_depth, zeros, zeros))

        if self.rrtmg is not None:
            from mfespf.physics.rrtmg_longwave import rrtmg_optics, rrtmg_reference_departure
            from mfespf.physics.solar_absorption import ozone_layer_path

            pressure = self.input_field("atmosphere_layer_pressure").require_units("Pa").values
            if self.ozone_path is None:
                self.ozone_path = ozone_layer_path(
                    pressure,
                    self.input_field("atmosphere_interface_pressure").require_units("Pa").values,
                    mass.values,
                )
            rrtmg_optics_now = rrtmg_optics(
                self.rrtmg,
                vapour=self.input_field("atmosphere_water_vapour").require_units("kg m^-2").values,
                layer_mass=mass.values,
                pressure=pressure,
                temperature=temperature.values,
                ozone_path=self.ozone_path,
            )
            gas_tau, fractions = rrtmg_optics_now.tau, rrtmg_optics_now.fractions
            self._rrtmg_floor = (rrtmg_optics_now.floored, rrtmg_optics_now.most_negative)
            self._rrtmg_departure = rrtmg_reference_departure(pressure, temperature.values)
            if self.k_distribution_source:
                from mfespf.physics.rrtmg_longwave import rrtmg_planck_fractions

                common["planck_fractions"] = rrtmg_planck_fractions(self.rrtmg, fractions)
            contributions.append((gas_tau, np.zeros_like(gas_tau), np.zeros_like(gas_tau)))

        if self.co2_band is not None:
            from mfespf.physics.co2_band import co2_optical_depth

            # Scaled for each layer's pressure and temperature every step: the
            # band's emissivity in thin air is what it exists to get right.
            k_co2, kind_co2, ppmv = self.co2_band
            co2_tau = co2_optical_depth(
                k_co2,
                kind_co2,
                co2_ppmv=ppmv,
                layer_mass=mass.values,
                pressure=self.input_field("atmosphere_layer_pressure").require_units("Pa").values,
                temperature=temperature.values,
            )
            contributions.append((co2_tau, np.zeros_like(co2_tau), np.zeros_like(co2_tau)))

        if self.ozone_infrared is not None:
            from mfespf.physics.ozone_infrared import ozone_infrared_optical_depth
            from mfespf.physics.solar_absorption import ozone_layer_path

            if self.ozone_path is None:
                self.ozone_path = ozone_layer_path(
                    self.input_field("atmosphere_layer_pressure").require_units("Pa").values,
                    self.input_field("atmosphere_interface_pressure").require_units("Pa").values,
                    mass.values,
                )
            # Tabulated in pressure and scaled in temperature, so evaluated
            # per layer every step.
            o3_tau = ozone_infrared_optical_depth(
                self.ozone_infrared,
                ozone_path=self.ozone_path,
                pressure=self.input_field("atmosphere_layer_pressure").require_units("Pa").values,
                temperature=temperature.values,
            )
            contributions.append((o3_tau, np.zeros_like(o3_tau), np.zeros_like(o3_tau)))

        if self.ozone_k is not None:
            from mfespf.physics.solar_absorption import ozone_layer_path, ozone_optical_depth

            # In the clear sky and under cloud alike, like the vapour. The
            # layers do not move, so the amount is placed once.
            if self.ozone_path is None:
                self.ozone_path = ozone_layer_path(
                    self.input_field("atmosphere_layer_pressure").require_units("Pa").values,
                    self.input_field("atmosphere_interface_pressure").require_units("Pa").values,
                    mass.values,
                )
            ozone_tau = ozone_optical_depth(self.ozone_k, self.ozone_path)
            contributions.append((ozone_tau, np.zeros_like(ozone_tau), np.zeros_like(ozone_tau)))

        if self.couple_volcanic:
            from mfespf.physics.volcanic import stratospheric_aerosol_optics

            self.volcanic_optical_depth = (
                self.input_field("stratospheric_aerosol_optical_depth").require_units("1").scalar()
            )
            self.volcanic_optics = stratospheric_aerosol_optics(
                self.grid,
                self.volcanic_optical_depth,
                int(np.asarray(mass.values).size),
                stratospheric_layers=self.stratospheric_layers,
            )
            contributions.append(self.volcanic_optics)

        # Kept, not just written into ``common``: the cloudy solve needs these
        # too. Everything in the clear sky is also present under the cloud.
        base_optics: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None
        if contributions:
            base_optics = combine_scatterers(*contributions)
            absorption, scattering, asymmetry = base_optics
            common["absorption_optical_depth"] = absorption
            common["scattering_optical_depth"] = scattering
            common["asymmetry"] = asymmetry

        # The clear-sky column is always solved. It is the reference the cloud
        # radiative effect is defined against, and when no cloud is active it is
        # also the answer -- so a no-cloud run costs exactly what it did before.
        self.clear_sky = solve_column_radiation(**common)  # type: ignore[arg-type]

        cloud = self._cloud_optics() if self.couple_clouds else None
        if self.couple_cloud_field:
            # An observed field: several disjoint cloud types, each its own
            # sub-column, weighted by the amount seen from above.
            self.subcolumns = self._solve_cloud_field(common, base_optics)
            self._cloud_fraction = sum(a for a, _ in self.subcolumns)
            self.overcast = None
            self.solution = _mix_many(self.clear_sky, self.subcolumns)
        elif cloud is None or cloud[0] <= 0.0:
            self._cloud_fraction = 0.0
            self.overcast = None
            self.solution = self.clear_sky
        else:
            fraction, optics = cloud
            self._cloud_fraction = fraction
            # Cloud droplets sit in the same layers as the water vapour and the
            # volcanic aerosol, so all three optical depths add and the layer's
            # asymmetry is the scattering-weighted mean. The cloudy column must
            # be the clear column *plus* the cloud: passing the cloud's
            # properties alone overwrites the base optics rather than adding to
            # them, which deletes the vapour from under the cloud and lets the
            # cloudy fraction radiate to space through an atmosphere that is not
            # there. That reverses the sign of the longwave cloud effect.
            absorption, scattering, asymmetry = combine_scatterers(
                (optics.absorption, optics.scattering, optics.asymmetry),
                *((base_optics,) if base_optics is not None else ()),
            )
            overcast = solve_column_radiation(
                **{
                    **common,
                    "absorption_optical_depth": absorption,
                    "scattering_optical_depth": scattering,
                    "asymmetry": asymmetry,
                }  # type: ignore[arg-type]
            )
            self.overcast = overcast
            self.solution = _mix_columns(self.clear_sky, overcast, fraction)

        heat_capacity = self.input_field("atmosphere_layer_heat_capacity").require_units(
            "J m^-2 K^-1"
        )
        self._layer_heat_capacity = heat_capacity.values
        return [
            TendencyContribution(
                state_id=ATMOSPHERE_TEMPERATURE_STATE,
                provider_module_id=self.module_id,
                label="radiative_flux_divergence",
                value=self.solution.layer_heating / heat_capacity.values,
                units="K s^-1",
                equation_id=EQUATION_IDS["schwarzschild"],
            )
        ]

    def emit_outputs(self) -> dict[str, Field]:
        ctx = self.context
        solution = self._require_solution()
        provenance = Provenance(producer=self.module_id, derivation=EQUATION_IDS["schwarzschild"])
        flux = {
            "quantity_kind": QuantityKind.FLUX,
            "time_semantics": TimeSemantics.INTERVAL_AVERAGE,
            "time_interval": ctx.dt,
            "timestamp": ctx.time,
            "spatial_support": SpatialSupport.GLOBAL_SCALAR,
            "provenance": provenance,
            "units": "W m^-2",
        }
        toa = solution.toa_index
        toa_up_sw, toa_up_lw = solution.partition(solution.upward[toa])
        toa_down_sw, _ = solution.partition(solution.downward[toa])
        surface_down_sw, _ = solution.partition(solution.downward[0])
        surface_up_sw, _ = solution.partition(solution.upward[0])
        surface_net_sw = surface_down_sw - surface_up_sw

        # Every sky this step was solved for, stacked. Two when a cloud is
        # active, because the diffusivity approximation and the two-stream
        # coefficients were applied to both of them and a diagnostic that
        # covered only the clear one would leave the more demanding column --
        # the optically thick one -- unreported.
        skies = [s for s in (self.clear_sky, self.overcast) if s is not None]
        skies.extend(solution for _, solution in self.subcolumns)
        optics = {
            "quantity_kind": QuantityKind.DIAGNOSTIC,
            "time_semantics": TimeSemantics.INSTANTANEOUS,
            "timestamp": ctx.time,
            "spatial_support": SpatialSupport.LAYER_AVERAGE,
            "dimensions": ("bound",),
            "sign_convention": SignConvention.POSITIVE_DEFINITE,
            "provenance": provenance,
        }

        outputs = {
            "layer_optical_depth_extent": Field(
                field_id="layer_optical_depth_extent",
                physical_name="smallest and largest layer extinction optical depth solved",
                symbol="tau",
                values=_extent_over(skies, "optical_depth"),
                units="1",
                **optics,  # type: ignore[arg-type]
            ),
            "layer_single_scattering_albedo_extent": Field(
                field_id="layer_single_scattering_albedo_extent",
                physical_name="smallest and largest layer single-scattering albedo solved",
                symbol="omega",
                values=_extent_over(skies, "single_scattering_albedo"),
                units="1",
                **optics,  # type: ignore[arg-type]
            ),
            "layer_emissivity_extent": Field(
                field_id="layer_emissivity_extent",
                physical_name="smallest and largest layer emissivity solved",
                symbol="e_i",
                values=_extent_over(skies, "layer_emissivity"),
                units="1",
                **optics,  # type: ignore[arg-type]
            ),
            "layer_radiative_heating": Field(
                field_id="layer_radiative_heating",
                physical_name="net radiative heating of each layer",
                symbol="H_i",
                values=np.array(solution.layer_heating),
                units="W m^-2",
                quantity_kind=QuantityKind.SOURCE_RATE,
                time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                time_interval=ctx.dt,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.LAYER_AVERAGE,
                dimensions=("layer",),
                coordinate_system="pressure_layers",
                sign_convention=SignConvention.POSITIVE_INTO_RECEIVER,
                provenance=provenance,
                uncertainty=Uncertainty(
                    method="none",
                    unresolved=(
                        "spectral_discretization",
                        "diffusivity_approximation",
                        "omitted_scattering",
                    ),
                ),
            ),
            "toa_downward_radiation": Field(
                field_id="toa_downward_radiation",
                physical_name="downward radiation at the top of the column",
                values=np.array(solution.interface_down(toa)),
                sign_convention=SignConvention.POSITIVE_DOWNWARD,
                **flux,  # type: ignore[arg-type]
            ),
            "toa_upward_radiation": Field(
                field_id="toa_upward_radiation",
                physical_name="upward radiation at the top of the column",
                values=np.array(solution.interface_up(toa)),
                sign_convention=SignConvention.POSITIVE_UPWARD,
                **flux,  # type: ignore[arg-type]
            ),
            "surface_downward_radiation": Field(
                field_id="surface_downward_radiation",
                physical_name="downward radiation arriving at the surface",
                values=np.array(solution.interface_down(0)),
                sign_convention=SignConvention.POSITIVE_DOWNWARD,
                **flux,  # type: ignore[arg-type]
            ),
            "surface_upward_radiation": Field(
                field_id="surface_upward_radiation",
                physical_name="upward radiation leaving the surface",
                values=np.array(solution.interface_up(0)),
                sign_convention=SignConvention.POSITIVE_UPWARD,
                **flux,  # type: ignore[arg-type]
            ),
            "surface_absorbed_radiation": Field(
                field_id="surface_absorbed_radiation",
                physical_name="radiation absorbed by the surface",
                values=np.array(solution.surface_absorbed),
                sign_convention=SignConvention.POSITIVE_INTO_RECEIVER,
                **flux,  # type: ignore[arg-type]
            ),
            "surface_emitted_radiation": Field(
                field_id="surface_emitted_radiation",
                physical_name="thermal radiation emitted by the surface",
                values=np.array(solution.surface_emitted),
                sign_convention=SignConvention.POSITIVE_UPWARD,
                **flux,  # type: ignore[arg-type]
            ),
            # RAD-002 and CLD-001: the shortwave and longwave components are
            # separately available before any net radiative effect is formed.
            "outgoing_longwave_radiation": Field(
                field_id="outgoing_longwave_radiation",
                physical_name="outgoing longwave radiation at the top of the column",
                values=np.array(toa_up_lw),
                sign_convention=SignConvention.POSITIVE_UPWARD,
                **flux,  # type: ignore[arg-type]
            ),
            "atmosphere_absorbed_shortwave": Field(
                field_id="atmosphere_absorbed_shortwave",
                physical_name="sunlight absorbed within the atmosphere",
                values=np.array(toa_down_sw - toa_up_sw - surface_net_sw),
                **flux,  # type: ignore[arg-type]
            ),
            "surface_absorbed_shortwave": Field(
                field_id="surface_absorbed_shortwave",
                physical_name="sunlight absorbed at the surface",
                values=np.array(surface_net_sw),
                **flux,  # type: ignore[arg-type]
            ),
            "reflected_shortwave_radiation": Field(
                field_id="reflected_shortwave_radiation",
                physical_name="reflected shortwave radiation at the top of the column",
                values=np.array(toa_up_sw),
                sign_convention=SignConvention.POSITIVE_UPWARD,
                **flux,  # type: ignore[arg-type]
            ),
        }

        if (self.couple_clouds or self.couple_cloud_field) and self.clear_sky is not None:
            clear_sw, clear_lw = self.clear_sky.partition(
                self.clear_sky.upward[self.clear_sky.toa_index]
            )
            effect_sw, effect_lw = self.cloud_radiative_effect()
            outputs.update(
                {
                    "clear_sky_reflected_shortwave_radiation": Field(
                        field_id="clear_sky_reflected_shortwave_radiation",
                        physical_name="reflected shortwave with the cloud removed",
                        values=np.array(clear_sw),
                        sign_convention=SignConvention.POSITIVE_UPWARD,
                        **flux,  # type: ignore[arg-type]
                    ),
                    "clear_sky_outgoing_longwave_radiation": Field(
                        field_id="clear_sky_outgoing_longwave_radiation",
                        physical_name="outgoing longwave with the cloud removed",
                        values=np.array(clear_lw),
                        sign_convention=SignConvention.POSITIVE_UPWARD,
                        **flux,  # type: ignore[arg-type]
                    ),
                    "cloud_shortwave_radiative_effect": Field(
                        field_id="cloud_shortwave_radiative_effect",
                        physical_name="shortwave cloud radiative effect at the top",
                        values=np.array(effect_sw),
                        sign_convention=SignConvention.POSITIVE_DOWNWARD,
                        **flux,  # type: ignore[arg-type]
                    ),
                    "cloud_longwave_radiative_effect": Field(
                        field_id="cloud_longwave_radiative_effect",
                        physical_name="longwave cloud radiative effect at the top",
                        values=np.array(effect_lw),
                        sign_convention=SignConvention.POSITIVE_DOWNWARD,
                        **flux,  # type: ignore[arg-type]
                    ),
                    # Formed only after both components are published above.
                    "cloud_net_radiative_effect": Field(
                        field_id="cloud_net_radiative_effect",
                        physical_name="net cloud radiative effect at the top",
                        values=np.array(effect_sw + effect_lw),
                        sign_convention=SignConvention.POSITIVE_DOWNWARD,
                        **flux,  # type: ignore[arg-type]
                    ),
                }
            )
        if self.rrtmg is not None and self._rrtmg_departure is not None:
            # What RRTMG's tables were asked for, and what its zero-floor did:
            # the floor is the published scheme's own, and it is counted, not hidden.
            diagnostic = {
                "quantity_kind": QuantityKind.DIAGNOSTIC,
                "time_semantics": TimeSemantics.INSTANTANEOUS,
                "timestamp": ctx.time,
                "provenance": Provenance(
                    producer=self.module_id, derivation=EQUATION_IDS["rrtmg_longwave"]
                ),
            }
            outputs.update(
                {
                    "rrtmg_reference_temperature_departure": Field(
                        field_id="rrtmg_reference_temperature_departure",
                        physical_name="layer temperature minus RRTMG's reference at its pressure",
                        symbol="dT_ref",
                        values=np.array(self._rrtmg_departure),
                        units="K",
                        spatial_support=SpatialSupport.LAYER_AVERAGE,
                        dimensions=("layer",),
                        **diagnostic,  # type: ignore[arg-type]
                    ),
                    "rrtmg_negative_optical_depths_floored": Field(
                        field_id="rrtmg_negative_optical_depths_floored",
                        physical_name="layer-bin optical depths RRTMG's zero-floor raised to zero",
                        values=np.array(float(self._rrtmg_floor[0])),
                        units="1",
                        spatial_support=SpatialSupport.GLOBAL_SCALAR,
                        **diagnostic,  # type: ignore[arg-type]
                    ),
                    "rrtmg_most_negative_optical_depth": Field(
                        field_id="rrtmg_most_negative_optical_depth",
                        physical_name="lowest RRTMG optical depth before its zero-floor",
                        values=np.array(self._rrtmg_floor[1]),
                        units="1",
                        spatial_support=SpatialSupport.GLOBAL_SCALAR,
                        **diagnostic,  # type: ignore[arg-type]
                    ),
                }
            )
        return outputs

    def report_uncertainty(self) -> UncertaintyReport:
        solution = self.solution
        quantified: dict[str, float] = {
            "flux_divergence_closure_residual_W_m2": (
                abs(solution.energy_residual()) if solution else 0.0
            ),
            "diffusivity_factor": self.diffusivity_factor,
        }
        quantified["spectral_capture_fraction_at_255K"] = self._capture_at_255k
        return UncertaintyReport(
            module_id=self.module_id,
            quantified=quantified,
            method=(
                "closure residual measured directly; the physical contributions are "
                "declared but not propagated"
            ),
            distributional_assumptions="unresolved",
            correlations="unresolved",
            unresolved=(
                "spectral_discretization_error",
                "diffusivity_approximation_error",
                "absorption_coefficient_parameter_uncertainty",
                "structural_uncertainty_of_omitting_scattering",
            ),
        )

    def _co2_wing_decay(self) -> float:
        """The CO2 band width calibrated for the background this run carries."""
        from mfespf.physics.greenhouse import co2_wing_decay_for

        return co2_wing_decay_for(active_vapour=self.couple_water)

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "spectral_grid_id": self.grid.grid_id,
            "spectral_resolution_identifier": self.grid.resolution_identifier,
            "n_bands": self.grid.n_bins,
            "spectral_edges": self.grid.edges.tolist(),
            # Recorded band by band: a band-mean would make two spectrally
            # different atmospheres hash alike.
            "absorption_coefficient_m2_kg": self.absorption_coefficient.tolist(),
            "diffusivity_factor": self.diffusivity_factor,
            "shortwave_longwave_cutoff_m": self.shortwave_cutoff_m,
            "gas_scattering": "not_represented",
            "k_distribution_source": self.k_distribution_source,
            "linear_in_tau_source": self.linear_in_tau_source,
            "cloud_scattering": "two_stream_with_delta_scaling"
            if self.couple_clouds
            else "no_clouds_active",
            "couple_clouds": self.couple_clouds,
            "couple_greenhouse": self.couple_greenhouse,
            "couple_volcanic": self.couple_volcanic,
            "couple_water": self.couple_water,
            "solar_vapour_absorption": (
                {
                    "source": "Chou and Suarez 1999, NASA/TM-1999-104606 vol. 15",
                    "k_m2_kg": self.solar_vapour_k.tolist(),
                }
                if self.solar_vapour_k is not None
                else "not_represented"
            ),
            **(
                {
                    "co2_15_micron_band": {
                        "source": "Chou, Suarez, Liang and Yan 2001, Tables 3 and 10",
                        "co2_ppmv": float(self.co2_band[2]),
                    }
                }
                if self.co2_band is not None
                else {}
            ),
            **(
                {
                    "rrtmg_longwave": {
                        "source": "AER-RC/RRTMG_LW 958edb1a (rrtmg-longwave.md)",
                        "co2_ppmv": self.rrtmg.co2_ppmv,
                        "ch4_ppb": self.rrtmg.ch4_ppb,
                        "n2o_ppb": self.rrtmg.n2o_ppb,
                        "o2_fraction": self.rrtmg.o2_fraction,
                    }
                }
                if self.rrtmg is not None
                else {}
            ),
            **(
                {
                    "water_vapour_longwave": {
                        "source": (
                            "Chou, Suarez, Liang and Yan 2001, Tables 1, 3, 4, 9 and 10 "
                            "(global-column-chou-longwave.md)"
                        ),
                        "tabulated_bands": list(self.water_longwave.tabulated),
                    }
                }
                if self.water_longwave is not None and self.couple_water
                else {}
            ),
            **(
                {
                    "ozone_9_6_micron_band": {
                        "source": (
                            "Chou, Suarez, Liang and Yan 2001 band-5 transmission tables "
                            "(WRF module_ra_goddard.F, commit 7042598715a7), fitted as a "
                            "pressure-tabulated k-distribution (ozone-infrared-band.md)"
                        ),
                    }
                }
                if self.ozone_infrared is not None
                else {}
            ),
            # Only when on, so a column without ozone hashes as before.
            **(
                {
                    "ozone_solar_absorption": {
                        "source": "Chou and Suarez 1999 Table 3; SWOOSH v2.73 1991-2020 profile",
                        "k_per_cm_atm": self.ozone_k.tolist(),
                    }
                }
                if self.ozone_k is not None
                else {}
            ),
            **({"stratospheric_layers": self.stratospheric_layers} if self.couple_volcanic else {}),
            **(
                {
                    "greenhouse_background_m2_kg": (
                        "water_vapour" if self.couple_water else self.greenhouse_background
                    ),
                    "co2_wing_decay_per_metre": self._co2_wing_decay(),
                }
                if self.couple_greenhouse
                else {}
            ),
        }

    def _cloud_optics(self) -> tuple[float, Any] | None:
        """Read the cloud state, if a cloud module published one this step."""
        fraction_field = self._inputs.get("cloud_fraction")
        if fraction_field is None:
            return None
        fraction = float(fraction_field.require_units("1").scalar())
        if fraction <= 0.0:
            return None
        liquid = self._inputs.get("cloud_liquid_water_path")
        ice = self._inputs.get("cloud_ice_water_path")
        radii = self._inputs.get("cloud_effective_radius")
        if liquid is None or ice is None or radii is None:
            raise ConfigurationError(
                "a cloud fraction was published without the water paths and effective "
                "radii the optical properties are derived from (CLD-002)"
            )
        from mfespf.physics.clouds import CloudProfile, cloud_optical_properties

        profile = CloudProfile(
            liquid_water_path=liquid.require_units("kg m^-2").values,
            ice_water_path=ice.require_units("kg m^-2").values,
            fraction=fraction,
            droplet_effective_radius=float(radii.require_units("m").values[0]),
            crystal_effective_radius=float(radii.require_units("m").values[1]),
        )
        return fraction, cloud_optical_properties(
            self.grid, profile, shortwave_cutoff_m=self.shortwave_cutoff_m
        )

    def _solve_cloud_field(
        self,
        common: dict[str, Any],
        base_optics: tuple[np.ndarray, np.ndarray, np.ndarray] | None,
    ) -> list[tuple[float, ColumnRadiationSolution]]:
        """Solve one sub-column per observed cloud type (EQ-ISCCP-CLOUD-FIELD).

        Each cloud is added to the clear column's own optics in the layer that
        contains its top: scattering for sunlight, absorption for thermal
        radiation. Added, not substituted, for the reason the prescribed cloud
        adds: the vapour and aerosol under a cloud are still there.
        """
        from mfespf.physics.cloud_field import CLOUD_ASYMMETRY, cloud_type_layers

        amounts = self.input_field("cloud_type_amount").require_units("1").values
        tops = self.input_field("cloud_type_top_pressure").require_units("Pa").values
        tau = self.input_field("cloud_type_visible_optical_thickness").require_units("1").values
        tau_ir = self.input_field("cloud_type_infrared_optical_thickness").require_units("1").values
        interfaces = self.input_field("atmosphere_interface_pressure").require_units("Pa").values
        layers = cloud_type_layers(tops, interfaces)
        n_layers = int(np.asarray(common["layer_mass"]).size)
        shortwave = self.grid.centers < self.shortwave_cutoff_m
        solved: list[tuple[float, ColumnRadiationSolution]] = []
        for amount, layer, visible, infrared in zip(amounts, layers, tau, tau_ir, strict=True):
            if amount <= 0.0:
                continue
            absorption = np.zeros((n_layers, self.grid.n_bins))
            scattering = np.zeros_like(absorption)
            asymmetry = np.zeros_like(absorption)
            scattering[layer, shortwave] = visible
            asymmetry[layer, shortwave] = CLOUD_ASYMMETRY
            absorption[layer, ~shortwave] = infrared
            combined = combine_scatterers(
                (absorption, scattering, asymmetry),
                *((base_optics,) if base_optics is not None else ()),
            )
            solution = solve_column_radiation(
                **{
                    **common,
                    "absorption_optical_depth": combined[0],
                    "scattering_optical_depth": combined[1],
                    "asymmetry": combined[2],
                }
            )
            solved.append((float(amount), solution))
        return solved

    def cloud_radiative_effect(self) -> tuple[float, float]:
        """Shortwave and longwave cloud radiative effect at the top, in W m^-2.

        Returned as a pair and never pre-summed. The two have opposite signs and
        comparable magnitudes, so the net is a small difference of large numbers;
        SPC-010 and CLD-001 both require the components to be available first.
        """
        if self.clear_sky is None or self.solution is None:
            return (0.0, 0.0)
        toa = self.solution.toa_index
        all_sky_sw, all_sky_lw = self.solution.partition(self.solution.upward[toa])
        clear_sw, clear_lw = self.clear_sky.partition(self.clear_sky.upward[toa])
        # Net downward is incident minus outgoing, and the incident is identical
        # in both columns, so the effect is the difference in what leaves.
        return (clear_sw - all_sky_sw, clear_lw - all_sky_lw)

    def _require_solution(self) -> ColumnRadiationSolution:
        if self.solution is None:
            raise ConfigurationError(
                f"module {self.module_id!r} was asked for outputs before it solved"
            )
        return self.solution
