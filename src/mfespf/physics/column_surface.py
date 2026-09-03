"""Slab surface beneath a radiating atmospheric column.

Charter 9.9 and section 19 stage 2; requirements LND-001, RAD-002, STA-001,
STA-006, CON-009.

The same slab thermodynamics as the airless surface, with one difference that
matters: it no longer computes its own radiation.  It publishes its optical
properties -- spectral albedo, and hence by Kirchhoff its emissivity -- and reads
back the fluxes the radiation module computed for the whole column.  That split
follows the charter's division of labour: the surface owns its thermal and
material state (9.9), the radiation module owns transport (9.3), and neither
duplicates the other's numbers.

Its transfers are with the atmosphere, not with space: energy leaving this
surface enters the column above it, and what reaches space is the column's
business. Under a transparent atmosphere those are numerically the same thing,
which is exactly what the Stage 1 regression test checks.
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
from mfespf.core.versioning import VersionRange
from mfespf.coupling.conservation import ModuleBudget
from mfespf.coupling.transfers import ConservedQuantity, Transfer
from mfespf.physics.equations import EQUATION_IDS
from mfespf.physics.surface import (
    SURFACE_INTERFACE,
    SURFACE_TEMPERATURE_STATE,
    V1,
)
from mfespf.registries.states import StateBudget, TendencyContribution

__all__ = ["ColumnSurfaceModule"]


class ColumnSurfaceModule(Module):
    """A slab surface whose radiation is solved by the column radiation module.

    A second implementation of the ``surface`` interface alongside the airless
    slab (KER-010): same contract, same owned state, different physical
    configuration and different required inputs.
    """

    def __init__(
        self,
        *,
        module_id: str = "surface",
        atmosphere_module_id: str = "atmosphere",
        grid: SpectralGrid,
        areal_heat_capacity: float = 2.1e8,
        albedo: float | np.ndarray = 0.3,
        initial_temperature: float = 250.0,
        status: ModuleStatus = ModuleStatus.ACTIVE,
    ) -> None:
        contract = ModuleContract(
            module_id=module_id,
            module_name="slab surface coupled to a radiating column",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.REDUCED,
            physical_basis=(
                "First law applied to an isothermal layer exchanging energy with the "
                "atmosphere above it by radiation only. Spectral emissivity follows from "
                "the albedo by Kirchhoff's law."
            ),
            governing_equation_ids=(
                EQUATION_IDS["slab"],
                EQUATION_IDS["planck"],
                EQUATION_IDS["forward_euler"],
            ),
            assumptions=(
                "the layer is isothermal and well mixed",
                "emissivity equals one minus albedo in every band (Kirchhoff)",
                "no conduction to depth, no turbulent or latent exchange",
                "heat capacity constant over the simulated temperature range",
            ),
            validity_domain=(
                "Radiative-only surface energy balance beneath a non-scattering column. "
                "Sensible and latent heat fluxes are absent, so the surface-to-air "
                "temperature difference this produces is not comparable with observations."
            ),
            owned_state_ids=frozenset({SURFACE_TEMPERATURE_STATE}),
            required_inputs=(
                "surface_downward_radiation",
                "surface_upward_radiation",
                "surface_absorbed_radiation",
            ),
            produced_outputs=(
                SURFACE_TEMPERATURE_STATE,
                "surface_spectral_albedo",
                "surface_net_radiative_flux",
            ),
            imported_transfers=("radiation_down_surface",),
            exported_transfers=("radiation_up_surface",),
            conserved_quantities=("energy",),
            provided_interfaces=(SURFACE_INTERFACE,),
            required_interfaces=(("spectral_radiation", VersionRange.parse(">=1.0.0,<2.0.0")),),
            uncertainty_outputs=(
                "initial_condition",
                "areal_heat_capacity_parameter",
                "albedo_parameter",
                "time_discretization",
                "missing_turbulent_fluxes_structural",
            ),
            correctness_tests=(
                "test_transparent_column_recovers_the_airless_result",
                "test_single_grey_layer_matches_the_analytic_greenhouse",
                "test_column_energy_budget_closes",
                "test_zero_insolation_column_only_cools",
            ),
            known_limitations=(
                "No sensible or latent heat flux, so this surface runs warmer relative to "
                "the air than a real one and the difference is not a model result.",
            ),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(n_bands)",
                memory_complexity="O(n_bands)",
                scaling_variables=("n_bands",),
                work_units_per_step=2.0,
                benchmark_configuration=f"{grid.n_bins} bands, global scalar",
                peak_memory_bytes=8 * grid.n_bins + 4096,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
                lower_fidelity_alternatives=("airless_slab_surface",),
                checkpoint_cost_bytes=8,
            ),
        )
        super().__init__(contract, status=status)
        self.atmosphere_module_id = atmosphere_module_id
        self.grid = grid
        self.areal_heat_capacity = areal_heat_capacity
        self.albedo = albedo
        self.initial_temperature = initial_temperature

        self.temperature = initial_temperature
        self._downward = 0.0
        self._upward = 0.0
        self._absorbed = 0.0
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
        albedo = np.asarray(self.albedo, dtype=float)
        if np.any(albedo < 0.0) or np.any(albedo > 1.0):
            raise ConfigurationError("albedo must lie in [0, 1] at every wavelength")
        if albedo.ndim == 1 and albedo.size != self.grid.n_bins:
            raise ConfigurationError(
                f"spectral albedo has {albedo.size} entries but the configured grid has "
                f"{self.grid.n_bins} bands"
            )
        if self.initial_temperature <= 0.0:
            raise ConfigurationError("initial temperature must be positive")

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        self.temperature = float(
            initial_state.get(SURFACE_TEMPERATURE_STATE, self.initial_temperature)
        )
        self.mark_initialized()

    # --- physics ----------------------------------------------------------
    @property
    def spectral_albedo(self) -> np.ndarray:
        albedo = np.asarray(self.albedo, dtype=float)
        if albedo.ndim == 0:
            return np.full(self.grid.n_bins, float(albedo))
        return albedo

    def stored_energy(self, temperature: float | None = None) -> float:
        T = self.temperature if temperature is None else temperature
        return self.areal_heat_capacity * T

    def compute_tendencies(self) -> list[TendencyContribution]:
        ctx = self.context
        self._downward = (
            self.input_field("surface_downward_radiation").require_units("W m^-2").scalar()
        )
        self._upward = self.input_field("surface_upward_radiation").require_units("W m^-2").scalar()
        self._absorbed = (
            self.input_field("surface_absorbed_radiation").require_units("W m^-2").scalar()
        )
        self._interval = (ctx.time, ctx.time + ctx.dt)
        self._stored_initial = self.stored_energy()

        # Net into the surface is what arrives minus what leaves. Reflection is
        # already inside the upward flux, so subtracting the whole upward stream
        # from the whole downward one counts it exactly once.
        net = self._downward - self._upward
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
    def emit_state(self) -> dict[str, Field]:
        """Temperature and optical properties, published before radiation solves.

        The radiation module runs ahead of this one in the splitting order, so it
        reads both from here rather than from ``emit_outputs``.
        """
        provenance = Provenance(producer=self.module_id, derivation=EQUATION_IDS["slab"])
        timestamp = self._context.time if self._context else 0.0
        return {
            SURFACE_TEMPERATURE_STATE: Field(
                field_id=SURFACE_TEMPERATURE_STATE,
                physical_name="slab surface temperature",
                symbol="T_s",
                values=np.array(self.temperature),
                units="K",
                quantity_kind=QuantityKind.INTENSIVE_STATE,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                timestamp=timestamp,
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
            "surface_spectral_albedo": Field(
                field_id="surface_spectral_albedo",
                physical_name="spectral albedo of the surface",
                symbol="alpha",
                values=self.spectral_albedo,
                units="1",
                quantity_kind=QuantityKind.BOUNDARY_VALUE,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                timestamp=timestamp,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                dimensions=("band",),
                spectral_coordinate=self.grid.grid_id,
                metadata={"spectral_grid": self.grid},
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                provenance=Provenance(
                    producer=self.module_id,
                    notes="emissivity is one minus this by Kirchhoff's law",
                ),
            ),
        }

    def emit_outputs(self) -> dict[str, Field]:
        ctx = self.context
        return {
            "surface_net_radiative_flux": Field(
                field_id="surface_net_radiative_flux",
                physical_name="net radiative flux into the surface",
                symbol="F_net",
                values=np.array(self._downward - self._upward),
                units="W m^-2",
                quantity_kind=QuantityKind.DIAGNOSTIC,
                time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                time_interval=ctx.dt,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                sign_convention=SignConvention.POSITIVE_DOWNWARD,
                provenance=Provenance(producer=self.module_id, derivation=EQUATION_IDS["slab"]),
            )
        }

    def report_transfers(self) -> tuple[list[Transfer], list[Transfer]]:
        start, end = self._interval
        dt = end - start
        common = {
            "quantity": ConservedQuantity.ENERGY,
            "units": "J m^-2",
            "start_time": start,
            "end_time": end,
            "spatial_mapping": "global_mean_over_planetary_surface",
        }
        exports = [
            Transfer(
                transfer_id=f"radiation_up_surface@{start:g}",
                sender_module=self.module_id,
                receiver_module=self.atmosphere_module_id,
                amount=self._upward * dt,
                channel="radiation_up_surface",
                equation_id=EQUATION_IDS["planck"],
                **common,  # type: ignore[arg-type]
            )
        ]
        imports = [
            Transfer(
                transfer_id=f"radiation_down_surface@{start:g}",
                sender_module=self.atmosphere_module_id,
                receiver_module=self.module_id,
                amount=self._downward * dt,
                channel="radiation_down_surface",
                **common,  # type: ignore[arg-type]
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
            method="none propagated at this fidelity",
            distributional_assumptions="unresolved",
            correlations="unresolved",
            unresolved=(
                "initial_condition_uncertainty",
                "areal_heat_capacity_parameter_uncertainty",
                "albedo_parameter_uncertainty",
                "forward_euler_time_discretization_error",
                "structural_uncertainty_of_omitting_turbulent_and_latent_fluxes",
            ),
        )

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "areal_heat_capacity_J_m2_K": self.areal_heat_capacity,
            "albedo": self.spectral_albedo.tolist(),
            "initial_temperature_K": self.initial_temperature,
            "emissivity_rule": "one_minus_albedo_kirchhoff",
            "time_integration": "forward_euler",
        }

    def checkpoint(self) -> dict[str, Any]:
        return {SURFACE_TEMPERATURE_STATE: self.temperature}

    def restore(self, checkpoint: dict[str, Any]) -> None:
        self.temperature = float(checkpoint[SURFACE_TEMPERATURE_STATE])
