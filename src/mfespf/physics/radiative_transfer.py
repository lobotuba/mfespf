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
    diffusivity_factor: float | None = None,
    shortwave_cutoff_m: float = _DEFAULT_SW_LW_CUTOFF_M,
) -> ColumnRadiationSolution:
    """Solve absorption, emission and surface reflection through the column.

    Two sweeps per band, downward then upward, because without scattering the
    downward field does not depend on the upward one. Surface reflection is
    computed between them, which is why the order matters and is fixed here
    rather than left to a caller.

    ``incident_spectral_irradiance`` is the spectral irradiance arriving at the
    top of the column in W m^-2 m^-1 on ``grid``; ``surface_albedo`` is the
    per-band reflectance, and the surface's emissivity is its complement by
    Kirchhoff's law.
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

    # Emissivity of each layer in each band. Kirchhoff's law makes absorptivity
    # and emissivity the same number, which is what keeps the column from
    # manufacturing energy at equilibrium.
    optical_depth = D * kappa[None, :] * mass[:, None]
    layer_emissivity = -np.expm1(-optical_depth)  # 1 - exp(-tau), accurate when tau is small
    transmittance = 1.0 - layer_emissivity

    layer_emission = band_exitance_profile(grid, temperature)

    n_interfaces = n_layers + 1
    downward = np.zeros((n_interfaces, n_bands))
    upward = np.zeros((n_interfaces, n_bands))

    downward[n_layers] = np.asarray(incident_spectral_irradiance, dtype=float) * grid.widths
    for i in range(n_layers - 1, -1, -1):
        downward[i] = downward[i + 1] * transmittance[i] + layer_emissivity[i] * layer_emission[i]

    surface_emission = band_exitance(grid, float(surface_temperature))
    surface_emissivity = 1.0 - albedo
    upward[0] = surface_emissivity * surface_emission + albedo * downward[0]
    for i in range(n_layers):
        upward[i + 1] = upward[i] * transmittance[i] + layer_emissivity[i] * layer_emission[i]

    # Heating is the flux divergence across the layer, so the column budget
    # closes identically rather than approximately.
    layer_heating = np.array(
        [
            float(
                np.sum(
                    (upward[i] - upward[i + 1]) + (downward[i + 1] - downward[i]),
                )
            )
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
                "Plane-parallel radiative transfer with absorption, thermal emission "
                "under Kirchhoff's law, and surface reflection; no scattering."
            ),
            governing_equation_ids=(
                EQUATION_IDS["schwarzschild"],
                EQUATION_IDS["beer_lambert"],
                EQUATION_IDS["planck"],
                EQUATION_IDS["diffusivity"],
            ),
            assumptions=(
                "plane-parallel, horizontally homogeneous column",
                "no scattering in any band",
                "local thermodynamic equilibrium, so emissivity equals absorptivity",
                "each layer emits at its single isothermal temperature",
                "the hemispheric angular integral is replaced by one diffusivity factor",
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
            ),
            known_limitations=(
                "No scattering, so a shortwave-scattering atmosphere cannot be "
                "represented at all -- not merely represented badly.",
                "One diffusivity factor for both solar beam and thermal diffuse paths.",
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
        self.solution: ColumnRadiationSolution | None = None
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

        self.solution = solve_column_radiation(
            grid=self.grid,
            incident_spectral_irradiance=incident.values,
            layer_temperature=temperature.values,
            layer_mass=mass.values,
            absorption_coefficient=self.absorption_coefficient,
            surface_temperature=surface_temperature,
            surface_albedo=albedo.values,
            diffusivity_factor=self.diffusivity_factor,
            shortwave_cutoff_m=self.shortwave_cutoff_m,
        )

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
            "scattering": "not_represented",
        }

    def _require_solution(self) -> ColumnRadiationSolution:
        if self.solution is None:
            raise ConfigurationError(
                f"module {self.module_id!r} was asked for outputs before it solved"
            )
        return self.solution
