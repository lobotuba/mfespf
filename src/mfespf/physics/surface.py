"""Slab surface module: the Stage 1 airless radiating body.

Charter 9.9 and section 19 stage 1; requirements LND-001, RAD-002, STA-001,
STA-006, CON-007-adjacent energy closure, VAL-001.

The surface owns one prognostic state, its temperature.  It absorbs incident
shortwave according to a spectral albedo, reflects the remainder, and emits as a
grey body.  Absorbed and emitted energy are accounted through ordinary transfers,
so the Stage 1 energy budget closes to round-off rather than to a tolerance
chosen to make it close.

Spectral absorption is integrated over whatever grid the incident field arrives
on.  With a grey albedo the result is identical on a one-band and an N-band grid,
which is the KER-009 interchangeability test; with a wavelength-dependent albedo
the two differ, and that difference is a physical result rather than a numerical
artifact.
"""

from __future__ import annotations

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
from mfespf.coupling.conservation import ModuleBudget
from mfespf.coupling.transfers import ConservedQuantity, Transfer
from mfespf.physics.equations import EQUATION_IDS
from mfespf.physics.radiation import greybody_exitance
from mfespf.registries.states import (
    StateBudget,
    StateEntry,
    StateKind,
    TendencyContribution,
    UpdateAuthority,
)

__all__ = [
    "SURFACE_INTERFACE",
    "SURFACE_TEMPERATURE_STATE",
    "SURFACE_WATER_STATE",
    "SlabSurfaceModule",
    "surface_state_entry",
    "surface_water_state_entry",
]

V1 = Version(1, 0, 0)

SURFACE_TEMPERATURE_STATE = "surface_temperature"

SURFACE_INTERFACE = InterfaceContract(
    interface_id="surface",
    version=V1,
    physical_semantics={
        "temperature_units": "K",
        "flux_units": "W m^-2",
        "flux_sign_convention": "positive_upward_away_from_surface",
        "spatial_support": "global_scalar",
        "state_owner": "surface",
    },
)


def surface_state_entry(
    *,
    owner: str = "surface",
    readers: frozenset[str] = frozenset(),
    tendency_providers: frozenset[str] = frozenset(),
) -> StateEntry:
    """The registry row for the surface temperature (STA-002, STA-003)."""
    return StateEntry(
        record_id=SURFACE_TEMPERATURE_STATE,
        version=V1,
        physical_name="slab surface temperature",
        symbol="T_s",
        owner_module_id=owner,
        units="K",
        state_kind=StateKind.PROGNOSTIC,
        conservation_participation=("energy",),
        authorized_readers=readers,
        authorized_tendency_providers=frozenset({owner}) | tendency_providers,
        update_authority=UpdateAuthority.OWNER_ONLY,
        initialization_source="manifest",
        checkpoint_required=True,
        uncertainty_representation="initial_condition_and_parameter_uncertainty",
    )


SURFACE_WATER_STATE = "surface_water_store"


def surface_water_state_entry(*, owner: str = "surface") -> StateEntry:
    """Registry row for the surface water reservoir (STA-002, WAT-001).

    A single bulk store: this fidelity has no soil, runoff or groundwater, so the
    charter's hydrology carve-outs (9.6) stay reserved and inactive rather than
    being approximated by this one number.
    """
    return StateEntry(
        record_id=SURFACE_WATER_STATE,
        version=V1,
        physical_name="liquid water stored at the surface per unit area",
        symbol="W_s",
        owner_module_id=owner,
        units="kg m^-2",
        state_kind=StateKind.PROGNOSTIC,
        phase_or_species="liquid_water",
        conservation_participation=("water",),
        authorized_tendency_providers=frozenset({owner}),
        update_authority=UpdateAuthority.OWNER_ONLY,
        initialization_source="manifest",
        checkpoint_required=True,
        uncertainty_representation="initial_condition_uncertainty",
    )


class SlabSurfaceModule(Module):
    """A single isothermal layer of declared areal heat capacity."""

    def __init__(
        self,
        *,
        module_id: str = "surface",
        space_module_id: str = "space",
        source_module_id: str = "orbit",
        areal_heat_capacity: float = 2.1e8,
        emissivity: float = 1.0,
        albedo: float | np.ndarray = 0.3,
        initial_temperature: float = 250.0,
        status: ModuleStatus = ModuleStatus.ACTIVE,
    ) -> None:
        contract = ModuleContract(
            module_id=module_id,
            module_name="slab surface with spectral albedo and grey emission",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.REDUCED,
            physical_basis=(
                "First law applied to an isothermal layer exchanging energy only by "
                "radiation, with spectrally resolved absorption and grey thermal emission."
            ),
            governing_equation_ids=(
                EQUATION_IDS["slab"],
                EQUATION_IDS["spectral_absorption"],
                EQUATION_IDS["greybody"],
                EQUATION_IDS["forward_euler"],
            ),
            assumptions=(
                "the layer is isothermal and well mixed",
                "no atmosphere, conduction, convection or latent exchange",
                "emissivity independent of wavelength",
                "heat capacity constant over the simulated temperature range",
            ),
            validity_domain=(
                "Airless or transparent-atmosphere configurations, temperatures within the "
                "grey-body equation's declared range, and time steps below the forward-Euler "
                "stability limit for the configured heat capacity."
            ),
            owned_state_ids=frozenset({SURFACE_TEMPERATURE_STATE}),
            required_inputs=("toa_incident_irradiance", "toa_spectral_irradiance"),
            produced_outputs=(
                "surface_temperature",
                "surface_absorbed_shortwave",
                "surface_reflected_shortwave",
                "surface_emitted_longwave",
                "surface_net_radiative_flux",
            ),
            imported_transfers=("toa_incident_energy",),
            exported_transfers=("reflected_shortwave_energy", "emitted_longwave_energy"),
            conserved_quantities=("energy",),
            internal_sources_and_sinks=(),
            provided_interfaces=(SURFACE_INTERFACE,),
            required_interfaces=(("orbital_geometry", VersionRange.parse(">=1.0.0,<2.0.0")),),
            uncertainty_outputs=(
                "initial_condition",
                "areal_heat_capacity_parameter",
                "albedo_parameter",
                "time_discretization",
            ),
            correctness_tests=(
                "test_slab_reaches_analytic_equilibrium",
                "test_slab_relaxation_timescale",
                "test_slab_energy_closure",
                "test_grey_albedo_matches_broadband",
                "test_zero_insolation_null_control",
            ),
            known_limitations=(
                "One temperature and one timescale by construction; it cannot exhibit the "
                "multiple response timescales the charter's timescale philosophy expects a "
                "real system to have.",
            ),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(n_bands)",
                memory_complexity="O(n_bands)",
                scaling_variables=("n_bands",),
                work_units_per_step=3.0,
                benchmark_configuration="global scalar, 16 bands",
                peak_memory_bytes=8192,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
                lower_fidelity_alternatives=("zero_heat_capacity_instant_equilibrium",),
                checkpoint_cost_bytes=64,
            ),
        )
        super().__init__(contract, status=status)
        self.space_module_id = space_module_id
        self.source_module_id = source_module_id
        self.areal_heat_capacity = areal_heat_capacity
        self.emissivity = emissivity
        self.albedo = albedo
        self.initial_temperature = initial_temperature

        self.temperature = initial_temperature
        self._incident = 0.0
        self._absorbed = 0.0
        self._reflected = 0.0
        self._emitted = 0.0
        self._interval = (0.0, 0.0)
        self._stored_initial = 0.0
        self._budget: StateBudget | None = None

    # --- configuration ----------------------------------------------------
    def validate_configuration(self) -> None:
        super().validate_configuration()
        if self.areal_heat_capacity <= 0.0:
            raise ConfigurationError(
                f"areal heat capacity {self.areal_heat_capacity} must be positive"
            )
        if not 0.0 <= self.emissivity <= 1.0:
            raise ConfigurationError(f"emissivity {self.emissivity} outside [0, 1]")
        albedo = np.asarray(self.albedo, dtype=float)
        if np.any(albedo < 0.0) or np.any(albedo > 1.0):
            raise ConfigurationError("albedo must lie in [0, 1] at every wavelength")
        if self.initial_temperature < 0.0:
            raise ConfigurationError("initial temperature must be non-negative")

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        self.temperature = float(
            initial_state.get(SURFACE_TEMPERATURE_STATE, self.initial_temperature)
        )
        self.mark_initialized()

    # --- physics ----------------------------------------------------------
    def stored_energy(self, temperature: float | None = None) -> float:
        """Areal energy content ``C * T`` in J m^-2, the budget's stored quantity."""
        T = self.temperature if temperature is None else temperature
        return self.areal_heat_capacity * T

    def _spectral_albedo(self, grid: SpectralGrid) -> np.ndarray:
        albedo = np.asarray(self.albedo, dtype=float)
        if albedo.ndim == 0:
            return np.full(grid.n_bins, float(albedo))
        if albedo.shape != (grid.n_bins,):
            raise ConfigurationError(
                f"spectral albedo has {albedo.shape} entries but the incident grid "
                f"{grid.grid_id!r} has {grid.n_bins} bins; regrid the albedo at a declared "
                f"boundary rather than broadcasting it here"
            )
        return albedo

    def compute_tendencies(self) -> list[TendencyContribution]:
        ctx = self.context
        incident_field = self.input_field("toa_incident_irradiance").require_units("W m^-2")
        spectral_field = self.input_field("toa_spectral_irradiance").require_units("W m^-2 m^-1")
        grid = spectral_field.metadata.get("spectral_grid")
        if not isinstance(grid, SpectralGrid):
            raise ConfigurationError(
                "toa_spectral_irradiance arrived without its spectral grid (SPC-004)"
            )

        albedo = self._spectral_albedo(grid)
        self._incident = incident_field.scalar()
        self._absorbed = float(np.sum((1.0 - albedo) * spectral_field.values * grid.widths))
        # Reflection is defined as the complement of absorption on the *same*
        # grid, so absorbed + reflected equals the spectrally integrated incident
        # energy by construction. Any mismatch against the broadband incident
        # field is a real disagreement between the two representations and is
        # reported below rather than absorbed into the budget.
        spectral_incident = float(np.sum(spectral_field.values * grid.widths))
        self._reflected = spectral_incident - self._absorbed
        self._emitted = greybody_exitance(self.temperature, self.emissivity)
        self._interval = (ctx.time, ctx.time + ctx.dt)
        self._stored_initial = self.stored_energy()

        net = self._absorbed - self._emitted
        return [
            TendencyContribution(
                state_id=SURFACE_TEMPERATURE_STATE,
                provider_module_id=self.module_id,
                label="net_radiative_flux",
                value=net / self.areal_heat_capacity,
                units="K s^-1",
                equation_id=EQUATION_IDS["slab"],
            )
        ]

    def advance_state(self, dt: float, contributions: list[TendencyContribution]) -> None:
        """Forward Euler, with the labeled contributions retained (STA-006)."""
        initial = np.asarray(self.temperature, dtype=float)
        total = sum(float(np.asarray(c.value)) for c in contributions)
        self.temperature = float(initial + total * dt)
        self._budget = StateBudget(
            state_id=SURFACE_TEMPERATURE_STATE,
            initial_value=initial,
            final_value=np.asarray(self.temperature, dtype=float),
            dt=dt,
            contributions=list(contributions),
        )

    def state_budgets(self) -> dict[str, Any]:
        return {SURFACE_TEMPERATURE_STATE: self._budget} if self._budget else {}

    # --- outputs ----------------------------------------------------------
    def emit_outputs(self) -> dict[str, Field]:
        ctx = self.context
        provenance = Provenance(producer=self.module_id, derivation=EQUATION_IDS["slab"])
        flux_common = {
            "quantity_kind": QuantityKind.FLUX,
            "time_semantics": TimeSemantics.INTERVAL_AVERAGE,
            "time_interval": ctx.dt,
            "timestamp": ctx.time,
            "spatial_support": SpatialSupport.GLOBAL_SCALAR,
            "sign_convention": SignConvention.POSITIVE_UPWARD,
            "provenance": provenance,
            "units": "W m^-2",
        }
        return {
            "surface_temperature": Field(
                field_id="surface_temperature",
                physical_name="slab surface temperature",
                symbol="T_s",
                values=np.array(self.temperature),
                units="K",
                quantity_kind=QuantityKind.INTENSIVE_STATE,
                time_semantics=TimeSemantics.START_OF_STEP,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                provenance=provenance,
                uncertainty=Uncertainty(
                    method="none",
                    unresolved=(
                        "initial_condition",
                        "areal_heat_capacity_parameter",
                        "time_discretization",
                    ),
                ),
            ),
            "surface_absorbed_shortwave": Field(
                field_id="surface_absorbed_shortwave",
                physical_name="absorbed shortwave flux density",
                symbol="A",
                values=np.array(self._absorbed),
                **flux_common,  # type: ignore[arg-type]
            ),
            "surface_reflected_shortwave": Field(
                field_id="surface_reflected_shortwave",
                physical_name="reflected shortwave flux density",
                symbol="R",
                values=np.array(self._reflected),
                **flux_common,  # type: ignore[arg-type]
            ),
            "surface_emitted_longwave": Field(
                field_id="surface_emitted_longwave",
                physical_name="emitted longwave flux density",
                symbol="M",
                values=np.array(self._emitted),
                **flux_common,  # type: ignore[arg-type]
            ),
            "surface_net_radiative_flux": Field(
                field_id="surface_net_radiative_flux",
                physical_name="net radiative flux into the surface",
                symbol="F_net",
                values=np.array(self._absorbed - self._emitted),
                quantity_kind=QuantityKind.DIAGNOSTIC,
                time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                time_interval=ctx.dt,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                sign_convention=SignConvention.POSITIVE_DOWNWARD,
                provenance=provenance,
                units="W m^-2",
            ),
        }

    def report_transfers(self) -> tuple[list[Transfer], list[Transfer]]:
        start, end = self._interval
        dt = end - start
        exports = [
            Transfer(
                transfer_id=f"reflected_shortwave@{start:g}",
                sender_module=self.module_id,
                receiver_module=self.space_module_id,
                quantity=ConservedQuantity.ENERGY,
                amount=self._reflected * dt,
                units="J m^-2",
                start_time=start,
                end_time=end,
                channel="reflected_shortwave",
                equation_id=EQUATION_IDS["spectral_absorption"],
            ),
            Transfer(
                transfer_id=f"emitted_longwave@{start:g}",
                sender_module=self.module_id,
                receiver_module=self.space_module_id,
                quantity=ConservedQuantity.ENERGY,
                amount=self._emitted * dt,
                units="J m^-2",
                start_time=start,
                end_time=end,
                channel="emitted_longwave",
                equation_id=EQUATION_IDS["greybody"],
            ),
        ]
        imports = [
            Transfer(
                transfer_id=f"toa_incident@{start:g}",
                sender_module=self.source_module_id,
                receiver_module=self.module_id,
                quantity=ConservedQuantity.ENERGY,
                amount=self._incident * dt,
                units="J m^-2",
                start_time=start,
                end_time=end,
                spatial_mapping="global_mean_over_planetary_surface",
                channel="radiation_down_toa",
            )
        ]
        return (exports, imports)

    def report_conservation(self) -> list[ModuleBudget]:
        return [
            ModuleBudget(
                module_id=self.module_id,
                quantity=ConservedQuantity.ENERGY,
                units="J m^-2",
                stored_initial=self._stored_initial,
                stored_final=self.stored_energy(),
            )
        ]

    def report_uncertainty(self) -> UncertaintyReport:
        return UncertaintyReport(
            module_id=self.module_id,
            quantified={},
            method="none propagated at this fidelity",
            distributional_assumptions="unresolved",
            correlations="unresolved",
            unresolved=(
                "initial_condition_uncertainty",
                "areal_heat_capacity_parameter_uncertainty",
                "albedo_parameter_uncertainty",
                "forward_euler_time_discretization_error",
                "structural_uncertainty_of_the_single_slab_representation",
            ),
        )

    def parameters(self) -> dict[str, Any]:
        albedo = np.asarray(self.albedo, dtype=float)
        return {
            **dict(self.configuration),
            "areal_heat_capacity_J_m2_K": self.areal_heat_capacity,
            "emissivity": self.emissivity,
            # Spectral albedo is recorded band by band: a band-mean would make
            # two physically different surfaces hash alike.
            "albedo": albedo.reshape(-1).tolist(),
            "initial_temperature_K": self.initial_temperature,
            "time_integration": "forward_euler",
        }

    def checkpoint(self) -> dict[str, Any]:
        return {SURFACE_TEMPERATURE_STATE: self.temperature}

    def restore(self, checkpoint: dict[str, Any]) -> None:
        self.temperature = float(checkpoint[SURFACE_TEMPERATURE_STATE])
