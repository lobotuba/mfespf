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
from mfespf.core.versioning import VersionRange
from mfespf.coupling.conservation import ModuleBudget
from mfespf.coupling.transfers import ConservedQuantity, Transfer
from mfespf.coupling.validity import INTEGRATED_OVER
from mfespf.physics.equations import EQUATION_IDS
from mfespf.physics.surface import (
    SURFACE_INTERFACE,
    SURFACE_TEMPERATURE_STATE,
    SURFACE_WATER_STATE,
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
        couple_water: bool = False,
        water_module_id: str = "water",
        initial_water_store: float = 1.0e4,
        latent_heat: float | None = None,
        couple_deep_ocean: bool = False,
        ocean_module_id: str = "ocean",
        status: ModuleStatus = ModuleStatus.ACTIVE,
    ) -> None:
        # The contract genuinely differs between a dry and a moist configuration:
        # a moist surface owns a second state, needs three more inputs and moves
        # a second conserved quantity. MOD-001 requires the contract to describe
        # the module as configured, so it is built rather than fixed.
        moist_inputs = (
            (
                "surface_evaporation_flux",
                "precipitation_flux",
                "surface_latent_heat_flux",
                "deposition_fusion_flux",
            )
            if couple_water
            else ()
        )
        contract = ModuleContract(
            module_id=module_id,
            module_name=(
                "slab surface coupled to a radiating column"
                + (" with evaporation and precipitation" if couple_water else "")
                + (" over a deep ocean reservoir" if couple_deep_ocean else "")
            ),
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
            equation_variables=(
                (
                    EQUATION_IDS["planck"],
                    (
                        # This module cites Planck for the provenance of the
                        # longwave flux it exchanges, but it never evaluates the
                        # function at a wavelength: the radiation module does
                        # that band by band and hands over a broadband flux.
                        # Binding lambda to this module's grid edges would state
                        # a wavelength support for a quantity that has none left.
                        ("lambda", INTEGRATED_OVER, "m"),
                        ("T", SURFACE_TEMPERATURE_STATE, "K"),
                    ),
                ),
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
            owned_state_ids=frozenset(
                {SURFACE_TEMPERATURE_STATE} | ({SURFACE_WATER_STATE} if couple_water else set())
            ),
            required_inputs=(
                "surface_downward_radiation",
                "surface_upward_radiation",
                "surface_absorbed_radiation",
                *moist_inputs,
                *(("deep_ocean_heat_flux",) if couple_deep_ocean else ()),
            ),
            produced_outputs=(
                SURFACE_TEMPERATURE_STATE,
                "surface_spectral_albedo",
                "surface_heat_capacity",
                "surface_net_radiative_flux",
                *((SURFACE_WATER_STATE, "surface_energy_balance") if couple_water else ()),
            ),
            imported_transfers=(
                "radiation_down_surface",
                *(("precipitation_water",) if couple_water else ()),
            ),
            exported_transfers=(
                "radiation_up_surface",
                *(
                    (
                        "evaporation_water",
                        "evaporation_latent_energy",
                        # The fusion energy the surface pays when the ice that
                        # formed aloft lands and melts. Declared on both sides,
                        # because a mechanism that moves a conserved quantity
                        # between two modules has to be.
                        "deposition_fusion_energy",
                    )
                    if couple_water
                    else ()
                ),
                *(("deep_ocean_heat",) if couple_deep_ocean else ()),
            ),
            conserved_quantities=("energy", "water") if couple_water else ("energy",),
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
                *(
                    (
                        "test_moist_column_budgets_close_for_water_and_energy",
                        "test_evaporation_conserves_water_and_carries_its_latent_energy",
                    )
                    if couple_water
                    else ()
                ),
            ),
            known_limitations=(
                (
                    "No sensible heat flux, so the surface-to-air temperature difference "
                    "is not comparable with observations."
                    if couple_water
                    else "No sensible or latent heat flux, so this surface runs warmer "
                    "relative to the air than a real one and the difference is not a "
                    "model result."
                ),
                *(
                    (
                        "The water store does not limit evaporation, so the surface "
                        "cannot dry out however much it loses.",
                    )
                    if couple_water
                    else ()
                ),
                *(
                    (
                        "The mixed layer is one isothermal reservoir; the slower "
                        "response comes from the deep store below it, not from any "
                        "structure within this one.",
                    )
                    if couple_deep_ocean
                    else ()
                ),
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

        self.couple_deep_ocean = couple_deep_ocean
        self.ocean_module_id = ocean_module_id
        self._deep_ocean_flux = 0.0

        self.couple_water = couple_water
        self.water_module_id = water_module_id
        self.initial_water_store = initial_water_store
        self.latent_heat = (
            latent_heat
            if latent_heat is not None
            else value_of("water_latent_heat_vaporization_273")
        )

        self.temperature = initial_temperature
        self.water_store = initial_water_store
        self._evaporation = 0.0
        self._precipitation = 0.0
        self._latent_flux = 0.0
        self._fusion_flux = 0.0
        self._water_initial = initial_water_store
        self._water_budget: StateBudget | None = None
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
        self.water_store = float(initial_state.get(SURFACE_WATER_STATE, self.initial_water_store))
        if self.water_store < 0.0:
            raise ConfigurationError("the surface water store cannot start negative")
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
        contributions = [
            TendencyContribution(
                state_id=SURFACE_TEMPERATURE_STATE,
                provider_module_id=self.module_id,
                label="net_radiative_flux",
                value=net / self.areal_heat_capacity,
                units="K s^-1",
                equation_id=EQUATION_IDS["slab"],
            )
        ]
        if self.couple_deep_ocean:
            # A labelled contribution of its own, next to the radiative one: the
            # surface energy balance is the thing being examined here, and a sum
            # would hide which of the two terms moved (STA-006, charter 9.9).
            self._deep_ocean_flux = (
                self.input_field("deep_ocean_heat_flux").require_units("W m^-2").scalar()
            )
            contributions.append(
                TendencyContribution(
                    state_id=SURFACE_TEMPERATURE_STATE,
                    provider_module_id=self.module_id,
                    label="deep_ocean_exchange",
                    value=-self._deep_ocean_flux / self.areal_heat_capacity,
                    units="K s^-1",
                    equation_id=EQUATION_IDS["deep_ocean_exchange"],
                )
            )

        if not self.couple_water:
            return contributions

        self._evaporation = (
            self.input_field("surface_evaporation_flux").require_units("kg m^-2 s^-1").scalar()
        )
        self._precipitation = (
            self.input_field("precipitation_flux").require_units("kg m^-2 s^-1").scalar()
        )
        self._latent_flux = (
            self.input_field("surface_latent_heat_flux").require_units("W m^-2").scalar()
        )
        # The ice that formed aloft melts when it lands, and the surface pays
        # for it. This is the other end of the fusion energy the deposition
        # released into the atmosphere; without it that energy would have come
        # from nowhere.
        self._fusion_flux = (
            self.input_field("deposition_fusion_flux").require_units("W m^-2").scalar()
        )
        self._water_initial = self.water_store

        # Evaporative cooling is a separate labeled contribution rather than
        # being folded into the net radiative term, so the surface energy budget
        # shows the radiative and turbulent parts apart (charter 9.9).
        contributions.append(
            TendencyContribution(
                state_id=SURFACE_TEMPERATURE_STATE,
                provider_module_id=self.module_id,
                label="evaporative_cooling",
                value=-self._latent_flux / self.areal_heat_capacity,
                units="K s^-1",
                equation_id=EQUATION_IDS["latent_energy"],
            )
        )
        if self._fusion_flux:
            contributions.append(
                TendencyContribution(
                    state_id=SURFACE_TEMPERATURE_STATE,
                    provider_module_id=self.module_id,
                    label="melting_of_frozen_precipitation",
                    value=-self._fusion_flux / self.areal_heat_capacity,
                    units="K s^-1",
                    equation_id=EQUATION_IDS["latent_energy"],
                )
            )
        contributions.append(
            TendencyContribution(
                state_id=SURFACE_WATER_STATE,
                provider_module_id=self.module_id,
                label="precipitation_minus_evaporation",
                value=self._precipitation - self._evaporation,
                units="kg m^-2 s^-1",
                equation_id=EQUATION_IDS["surface_evaporation"],
            )
        )
        return contributions

    def advance_state(self, dt: float, contributions: list[TendencyContribution]) -> None:
        thermal = [c for c in contributions if c.state_id == SURFACE_TEMPERATURE_STATE]
        initial = np.asarray(self.temperature, dtype=float)
        self.temperature = float(initial + sum(float(np.asarray(c.value)) for c in thermal) * dt)
        if self.temperature <= 0.0:
            raise ConfigurationError(
                f"surface temperature reached {self.temperature:.4g} K, which is "
                f"unphysical; the run is stopped rather than continued from an "
                f"impossible state (charter 15.2)"
            )
        self._budget = StateBudget(
            state_id=SURFACE_TEMPERATURE_STATE,
            initial_value=initial,
            final_value=np.asarray(self.temperature, dtype=float),
            dt=dt,
            contributions=thermal,
        )

        hydrological = [c for c in contributions if c.state_id == SURFACE_WATER_STATE]
        if not hydrological:
            return
        water_initial = np.asarray(self.water_store, dtype=float)
        self.water_store = float(
            water_initial + sum(float(np.asarray(c.value)) for c in hydrological) * dt
        )
        self._water_budget = StateBudget(
            state_id=SURFACE_WATER_STATE,
            initial_value=water_initial,
            final_value=np.asarray(self.water_store, dtype=float),
            dt=dt,
            contributions=hydrological,
        )

    def state_budgets(self) -> dict[str, Any]:
        budgets: dict[str, Any] = {}
        if self._budget:
            budgets[SURFACE_TEMPERATURE_STATE] = self._budget
        if self._water_budget:
            budgets[SURFACE_WATER_STATE] = self._water_budget
        return budgets

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
            **(
                {
                    SURFACE_WATER_STATE: Field(
                        field_id=SURFACE_WATER_STATE,
                        physical_name="liquid water stored at the surface",
                        symbol="W_s",
                        values=np.array(self.water_store),
                        units="kg m^-2",
                        quantity_kind=QuantityKind.EXTENSIVE_STATE,
                        time_semantics=TimeSemantics.INSTANTANEOUS,
                        timestamp=timestamp,
                        spatial_support=SpatialSupport.GLOBAL_SCALAR,
                        phase_or_species="liquid_water",
                        sign_convention=SignConvention.POSITIVE_DEFINITE,
                        provenance=Provenance(
                            producer=self.module_id,
                            derivation=EQUATION_IDS["surface_evaporation"],
                        ),
                    )
                }
                if self.couple_water
                else {}
            ),
            "surface_heat_capacity": Field(
                field_id="surface_heat_capacity",
                physical_name="areal heat capacity of the surface slab",
                symbol="C_s",
                values=np.array(self.areal_heat_capacity),
                units="J m^-2 K^-1",
                quantity_kind=QuantityKind.EXTENSIVE_STATE,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                timestamp=timestamp,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                provenance=provenance,
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
            ),
            **(
                {
                    "surface_energy_balance": Field(
                        field_id="surface_energy_balance",
                        physical_name="net energy into the surface after the latent loss",
                        values=np.array(self._downward - self._upward - self._latent_flux),
                        units="W m^-2",
                        quantity_kind=QuantityKind.DIAGNOSTIC,
                        time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                        time_interval=ctx.dt,
                        timestamp=ctx.time,
                        spatial_support=SpatialSupport.GLOBAL_SCALAR,
                        sign_convention=SignConvention.POSITIVE_DOWNWARD,
                        provenance=Provenance(
                            producer=self.module_id, derivation=EQUATION_IDS["slab"]
                        ),
                    ),
                }
                if self.couple_water
                else {}
            ),
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
        if self.couple_deep_ocean:
            exports.append(
                Transfer(
                    transfer_id=f"deep_ocean_heat@{start:g}",
                    sender_module=self.module_id,
                    receiver_module=self.ocean_module_id,
                    amount=self._deep_ocean_flux * dt,
                    channel="deep_ocean_heat",
                    equation_id=EQUATION_IDS["deep_ocean_exchange"],
                    **common,  # type: ignore[arg-type]
                )
            )
        sensible = self._inputs.get("surface_sensible_heat_flux")
        if sensible is not None:
            exports.append(
                Transfer(
                    transfer_id=f"sensible_heat@{start:g}",
                    sender_module=self.module_id,
                    receiver_module=self.atmosphere_module_id,
                    amount=sensible.require_units("W m^-2").scalar() * dt,
                    channel="sensible_heat",
                    equation_id=EQUATION_IDS["sensible_heat"],
                    **common,  # type: ignore[arg-type]
                )
            )
        # Present when convection mixes the surface into the column rather than
        # exchanging with it through a bulk coefficient. Both sides declare the
        # transfer, so an exchange that only one module believed in shows up as
        # an unmatched record instead of quietly balancing.
        convective = self._inputs.get("surface_convective_heat_flux")
        if convective is not None:
            exports.append(
                Transfer(
                    transfer_id=f"surface_convective_heat@{start:g}",
                    sender_module=self.module_id,
                    receiver_module=self.atmosphere_module_id,
                    amount=convective.require_units("W m^-2").scalar() * dt,
                    channel="surface_convective_heat",
                    equation_id=EQUATION_IDS["convective_adjustment"],
                    **common,  # type: ignore[arg-type]
                )
            )

        if not self.couple_water:
            return (exports, imports)

        water_common = {
            "start_time": start,
            "end_time": end,
            "spatial_mapping": "global_mean_over_planetary_surface",
        }
        exports.extend(
            [
                Transfer(
                    transfer_id=f"evaporation_water@{start:g}",
                    sender_module=self.module_id,
                    receiver_module=self.water_module_id,
                    quantity=ConservedQuantity.WATER,
                    amount=self._evaporation * dt,
                    units="kg m^-2",
                    channel="evaporation_water",
                    species="water_vapour",
                    equation_id=EQUATION_IDS["surface_evaporation"],
                    **water_common,  # type: ignore[arg-type]
                ),
                Transfer(
                    transfer_id=f"evaporation_latent_energy@{start:g}",
                    sender_module=self.module_id,
                    receiver_module=self.water_module_id,
                    quantity=ConservedQuantity.ENERGY,
                    amount=self._latent_flux * dt,
                    units="J m^-2",
                    channel="evaporation_latent_energy",
                    equation_id=EQUATION_IDS["latent_energy"],
                    **water_common,  # type: ignore[arg-type]
                ),
                Transfer(
                    transfer_id=f"deposition_fusion_energy@{start:g}",
                    sender_module=self.module_id,
                    receiver_module=self.water_module_id,
                    quantity=ConservedQuantity.ENERGY,
                    amount=self._fusion_flux * dt,
                    units="J m^-2",
                    channel="deposition_fusion_energy",
                    equation_id=EQUATION_IDS["latent_energy"],
                    **water_common,  # type: ignore[arg-type]
                ),
            ]
        )
        imports.append(
            Transfer(
                transfer_id=f"precipitation_water@{start:g}",
                sender_module=self.water_module_id,
                receiver_module=self.module_id,
                quantity=ConservedQuantity.WATER,
                amount=self._precipitation * dt,
                units="kg m^-2",
                channel="precipitation_water",
                species="liquid_water",
                equation_id=EQUATION_IDS["condensation"],
                **water_common,  # type: ignore[arg-type]
            )
        )
        return (exports, imports)

    def report_conservation(self) -> list[ModuleBudget]:
        budgets = [
            ModuleBudget(
                module_id=self.module_id,
                quantity=ConservedQuantity.ENERGY,
                units="J m^-2",
                stored_initial=self._stored_initial,
                stored_final=self.stored_energy(),
            )
        ]
        if self.couple_water:
            budgets.append(
                ModuleBudget(
                    module_id=self.module_id,
                    quantity=ConservedQuantity.WATER,
                    units="kg m^-2",
                    stored_initial=self._water_initial,
                    stored_final=self.water_store,
                )
            )
        return budgets

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
            "couple_water": self.couple_water,
            "couple_deep_ocean": self.couple_deep_ocean,
            **(
                {
                    "initial_water_store_kg_m2": self.initial_water_store,
                    "latent_heat_J_kg": self.latent_heat,
                }
                if self.couple_water
                else {}
            ),
        }

    def checkpoint(self) -> dict[str, Any]:
        return {
            SURFACE_TEMPERATURE_STATE: self.temperature,
            SURFACE_WATER_STATE: self.water_store,
        }

    def restore(self, checkpoint: dict[str, Any]) -> None:
        self.temperature = float(checkpoint[SURFACE_TEMPERATURE_STATE])
        self.water_store = float(checkpoint.get(SURFACE_WATER_STATE, self.water_store))
