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
from typing import Any

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
from mfespf.core.versioning import InterfaceContract, Version, VersionRange
from mfespf.physics.atmosphere import ATMOSPHERE_TEMPERATURE_STATE
from mfespf.physics.equations import EQUATION_IDS
from mfespf.physics.radiation import (
    band_exitance,
    band_exitance_capture,
    band_exitance_profile,
)
from mfespf.registries.states import TendencyContribution

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


def solve_column_radiation(
    *,
    grid: SpectralGrid,
    incident_spectral_irradiance: np.ndarray,
    layer_temperature: np.ndarray,
    layer_mass: np.ndarray,
    absorption_coefficient: np.ndarray,
    surface_temperature: float,
    surface_albedo: np.ndarray,
    cloud_absorption_optical_depth: np.ndarray | None = None,
    cloud_scattering_optical_depth: np.ndarray | None = None,
    cloud_asymmetry: np.ndarray | None = None,
    diffusivity_factor: float | None = None,
    shortwave_cutoff_m: float = _DEFAULT_SW_LW_CUTOFF_M,
) -> ColumnRadiationSolution:
    """Solve absorption, scattering, emission and surface reflection in the column.

    Solved by the adding method: sweeping upward, accumulate the reflectance and
    the emitted source of everything below each interface, then sweep down once
    with the incident field. That handles multiple reflection between layers
    exactly, which a pair of independent sweeps cannot once anything scatters.

    ``cloud_*`` arrays are optional and shaped ``(n_layers, n_bands)``. Omitting
    them leaves the column non-scattering, in which case the result is the
    absorption-only solution to floating-point round-off.
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
            raise ConfigurationError("cloud optical depths cannot be negative")
        return array

    cloud_absorption = _optional(cloud_absorption_optical_depth)
    cloud_scattering = _optional(cloud_scattering_optical_depth)
    asymmetry = (
        np.zeros(shape)
        if cloud_asymmetry is None
        else np.broadcast_to(np.asarray(cloud_asymmetry, dtype=float), shape)
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

    # Upward pass: reflectance and emitted source of everything below interface i.
    below_reflectance = np.zeros((n_layers + 1, n_bands))
    below_source = np.zeros((n_layers + 1, n_bands))
    surface_emissivity = 1.0 - albedo
    below_reflectance[0] = albedo
    below_source[0] = surface_emissivity * band_exitance(grid, float(surface_temperature))

    for i in range(n_layers):
        coupling = 1.0 - reflectance[i] * below_reflectance[i]
        source = emissivity[i] * layer_emission[i]
        below_reflectance[i + 1] = (
            reflectance[i] + (transmittance[i] ** 2 * below_reflectance[i]) / coupling
        )
        below_source[i + 1] = (
            transmittance[i]
            * below_reflectance[i]
            * (reflectance[i] * below_source[i] + source)
            / coupling
            + transmittance[i] * below_source[i]
            + source
        )

    # Downward pass with the incident field, then the upward field follows.
    downward = np.zeros((n_layers + 1, n_bands))
    downward[n_layers] = np.asarray(incident_spectral_irradiance, dtype=float) * grid.widths
    for i in range(n_layers - 1, -1, -1):
        coupling = 1.0 - reflectance[i] * below_reflectance[i]
        source = emissivity[i] * layer_emission[i]
        downward[i] = (
            transmittance[i] * downward[i + 1] + reflectance[i] * below_source[i] + source
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
        surface_emitted=float(
            np.sum(surface_emissivity * band_exitance(grid, float(surface_temperature)))
        ),
        grid=grid,
        shortwave_mask=shortwave_mask,
    )


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
        greenhouse_background: float = 1.0e-4,
        status: ModuleStatus = ModuleStatus.ACTIVE,
    ) -> None:
        n_bands = grid.n_bins
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
                *(
                    (EQUATION_IDS["cloud_optics"], EQUATION_IDS["cloud_radiative_effect"])
                    if couple_clouds
                    else ()
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
                "Optical depths where the diffusivity approximation holds to a few "
                "percent, roughly 0.1 to 10 per layer. Scattering atmospheres, cloudy "
                "columns and strongly inhomogeneous layers are outside it."
            ),
            required_inputs=(
                "toa_spectral_irradiance",
                ATMOSPHERE_TEMPERATURE_STATE,
                "atmosphere_layer_mass",
                "atmosphere_layer_heat_capacity",
                "surface_temperature",
                "surface_spectral_albedo",
                *(("cloud_fraction",) if couple_clouds else ()),
                *(("co2_mole_fraction",) if couple_greenhouse else ()),
            ),
            produced_outputs=(
                "layer_radiative_heating",
                "toa_downward_radiation",
                "toa_upward_radiation",
                "surface_downward_radiation",
                "surface_upward_radiation",
                "surface_absorbed_radiation",
                "outgoing_longwave_radiation",
                "reflected_shortwave_radiation",
                "atmospheric_radiative_forcing_diagnostics",
                *(
                    (
                        "clear_sky_outgoing_longwave_radiation",
                        "clear_sky_reflected_shortwave_radiation",
                        "cloud_shortwave_radiative_effect",
                        "cloud_longwave_radiative_effect",
                        "cloud_net_radiative_effect",
                    )
                    if couple_clouds
                    else ()
                ),
            ),
            conserved_quantities=("energy",),
            supported_spectral_representations=(grid.resolution_identifier,),
            provided_interfaces=(RADIATION_INTERFACE,),
            required_interfaces=(
                ("atmospheric_state", VersionRange.parse(">=1.0.0,<2.0.0")),
                ("surface", VersionRange.parse(">=1.0.0,<2.0.0")),
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
        self.couple_greenhouse = couple_greenhouse
        self.greenhouse_module_id = greenhouse_module_id
        self.greenhouse_background = greenhouse_background
        self.co2_mole_fraction = 0.0
        self.cloud_module_id = cloud_module_id
        self.solution: ColumnRadiationSolution | None = None
        self.clear_sky: ColumnRadiationSolution | None = None
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
            kappa = co2_band_absorption(
                self.grid, self.co2_mole_fraction, background=self.greenhouse_background
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

        # The clear-sky column is always solved. It is the reference the cloud
        # radiative effect is defined against, and when no cloud is active it is
        # also the answer -- so a no-cloud run costs exactly what it did before.
        self.clear_sky = solve_column_radiation(**common)  # type: ignore[arg-type]

        cloud = self._cloud_optics() if self.couple_clouds else None
        if cloud is None or cloud[0] <= 0.0:
            self._cloud_fraction = 0.0
            self.solution = self.clear_sky
        else:
            fraction, optics = cloud
            self._cloud_fraction = fraction
            overcast = solve_column_radiation(
                cloud_absorption_optical_depth=optics.absorption,
                cloud_scattering_optical_depth=optics.scattering,
                cloud_asymmetry=optics.asymmetry,
                **common,  # type: ignore[arg-type]
            )
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

        outputs = {
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
            "reflected_shortwave_radiation": Field(
                field_id="reflected_shortwave_radiation",
                physical_name="reflected shortwave radiation at the top of the column",
                values=np.array(toa_up_sw),
                sign_convention=SignConvention.POSITIVE_UPWARD,
                **flux,  # type: ignore[arg-type]
            ),
        }

        if self.couple_clouds and self.clear_sky is not None:
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

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "spectral_grid_id": self.grid.grid_id,
            "spectral_resolution_identifier": self.grid.resolution_identifier,
            "n_bands": self.grid.n_bins,
            # Recorded band by band: a band-mean would make two spectrally
            # different atmospheres hash alike.
            "absorption_coefficient_m2_kg": self.absorption_coefficient.tolist(),
            "diffusivity_factor": self.diffusivity_factor,
            "shortwave_longwave_cutoff_m": self.shortwave_cutoff_m,
            "gas_scattering": "not_represented",
            "cloud_scattering": "two_stream_with_delta_scaling"
            if self.couple_clouds
            else "no_clouds_active",
            "couple_clouds": self.couple_clouds,
            "couple_greenhouse": self.couple_greenhouse,
            **(
                {"greenhouse_background_m2_kg": self.greenhouse_background}
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
