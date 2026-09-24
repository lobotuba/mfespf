"""Atmospheric water vapour, evaporation, condensation and precipitation.

Charter 9.6 and section 19 stage 3; requirements WAT-001, WAT-002, CON-001,
CON-007, CON-009.

This is the first module to move a conserved quantity that is not energy, and
the first to move two at once. Both budgets must close simultaneously, which is
what CON-007 asks for when it says a phase change shall reconcile transferred
matter *and* latent energy.

The accounting choice that makes that possible: **the water module owns the
latent energy**, not the atmosphere. Its energy store is ``L_v`` times the
vapour it holds. Evaporation is then an ordinary transfer of two quantities from
the surface -- water mass, and the energy that was spent lifting it into the
vapour phase -- and condensation is a transfer of that energy onward to the
atmosphere as sensible heat, while the water mass returns to the surface as
precipitation carrying none.

Without that store the budget would not close between the two events: energy
would leave the surface at evaporation and appear in the atmosphere at
condensation with nothing holding it in between.
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
from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError
from mfespf.core.quantities import (
    Field,
    MissingValuePolicy,
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
from mfespf.physics.atmosphere import ATMOSPHERE_TEMPERATURE_STATE
from mfespf.physics.equations import EQUATION_IDS
from mfespf.physics.plume import (
    BULK_PLUME_PARAMETER,
    EVAPORATING_CONDENSATE_FRACTION,
    EVAPORATING_FRACTION_RANGE,
    bulk_plume_relative_humidity,
    convective_mass_flux,
    detrainment_rate,
    equilibrium_relative_humidity,
    fractional_entrainment_rate,
    water_vapour_lapse_rate,
)
from mfespf.physics.thermodynamics import saturation_specific_humidity
from mfespf.registries.states import (
    StateBudget,
    StateEntry,
    StateKind,
    TendencyContribution,
    UpdateAuthority,
)

__all__ = [
    "VAPOUR_STATE",
    "WATER_INTERFACE",
    "EntrainingPlumeHumidityModule",
    "PrescribedHumidityModule",
    "WaterModule",
    "vapour_state_entry",
]

_GAS_CONSTANT_DRY_AIR = value_of("gas_constant_dry_air")

V1 = Version(1, 0, 0)

VAPOUR_STATE = "atmosphere_water_vapour"

WATER_INTERFACE = InterfaceContract(
    interface_id="atmospheric_water",
    version=V1,
    physical_semantics={
        "vapour_units": "kg m^-2",
        "vapour_basis": "mass_per_unit_area_in_each_layer",
        "layer_ordering": "index_0_is_the_lowest_layer",
        "latent_energy_owner": "water",
        "precipitation_energy_content": "neglected",
        "state_owner": "water",
    },
)


def vapour_state_entry(
    *,
    owner: str = "water",
    n_layers: int,
    readers: frozenset[str] = frozenset(),
    tendency_providers: frozenset[str] = frozenset(),
) -> StateEntry:
    """Registry row for the layer vapour masses (STA-002, STA-003)."""
    return StateEntry(
        record_id=VAPOUR_STATE,
        version=V1,
        physical_name="water vapour mass per unit area in each layer",
        symbol="w",
        owner_module_id=owner,
        units="kg m^-2",
        state_kind=StateKind.PROGNOSTIC,
        dimensions=("layer",),
        phase_or_species="water_vapour",
        conservation_participation=("water", "energy"),
        authorized_readers=readers,
        authorized_tendency_providers=frozenset({owner}) | tendency_providers,
        update_authority=UpdateAuthority.OWNER_ONLY,
        initialization_source="manifest",
        checkpoint_required=True,
        uncertainty_representation="initial_condition_and_parameter_uncertainty",
        validity_flags=(f"n_layers={n_layers}",),
    )


class WaterModule(Module):
    """Owns atmospheric vapour and the latent energy it carries."""

    def __init__(
        self,
        *,
        module_id: str = "water",
        atmosphere_module_id: str = "atmosphere",
        surface_module_id: str = "surface",
        n_layers: int = 8,
        evaporation_coefficient: float = 2.0e-4,
        initial_relative_humidity: float = 0.0,
        latent_heat: float | None = None,
        status: ModuleStatus = ModuleStatus.ACTIVE,
    ) -> None:
        contract = ModuleContract(
            module_id=module_id,
            module_name="atmospheric water vapour with saturation adjustment",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.REDUCED,
            physical_basis=(
                "Conservation of water mass between a surface reservoir and atmospheric "
                "vapour, with bulk evaporation, saturation-adjustment condensation, and "
                "immediate precipitation."
            ),
            governing_equation_ids=(
                EQUATION_IDS["clausius_clapeyron"],
                EQUATION_IDS["specific_humidity"],
                EQUATION_IDS["surface_evaporation"],
                EQUATION_IDS["condensation"],
                EQUATION_IDS["latent_energy"],
            ),
            equation_variables=(
                (
                    EQUATION_IDS["clausius_clapeyron"],
                    (
                        # Saturation is evaluated twice per step and at two
                        # different temperatures: once per layer for the
                        # adjustment, and once at the surface for the
                        # evaporative deficit. Both are bound, because binding
                        # one would leave the other evaluated and unjudged.
                        ("T", ATMOSPHERE_TEMPERATURE_STATE, "K"),
                        ("T", "surface_temperature", "K"),
                    ),
                ),
                (
                    EQUATION_IDS["latent_energy"],
                    (
                        # Not the layer temperature. This relation sums over
                        # w_i, so it is evaluated only where there is vapour,
                        # and binding the whole temperature field flagged a
                        # column whose every out-of-range layer carried a weight
                        # of exactly zero. The masked field says where the sum
                        # actually has terms.
                        ("T", "vapour_bearing_layer_temperature", "K"),
                    ),
                ),
            ),
            assumptions=(
                "no vertical or horizontal transport of vapour between layers",
                "the surface is saturated with unlimited water available",
                "condensate precipitates immediately, so no cloud water is stored",
                "precipitation carries no sensible heat",
                "constant latent heat of vaporisation",
            ),
            validity_domain=(
                "Columns where vapour transport and cloud storage may be neglected. With "
                "no transport, the vapour profile is set entirely by local evaporation "
                "and condensation, so no claim about humidity structure aloft is in scope."
            ),
            owned_state_ids=frozenset({VAPOUR_STATE}),
            required_inputs=(
                ATMOSPHERE_TEMPERATURE_STATE,
                "atmosphere_layer_mass",
                "atmosphere_layer_pressure",
                "atmosphere_layer_heat_capacity",
                "surface_temperature",
            ),
            produced_outputs=(
                VAPOUR_STATE,
                "atmosphere_specific_humidity",
                "atmosphere_relative_humidity",
                "surface_evaporation_flux",
                "precipitation_flux",
                "surface_latent_heat_flux",
                "latent_heating_rate",
                "atmosphere_precipitable_water",
                "vapour_bearing_layer_temperature",
            ),
            imported_transfers=("evaporation_water", "evaporation_latent_energy"),
            exported_transfers=("precipitation_water", "latent_heat_release"),
            conserved_quantities=("water", "energy"),
            internal_sources_and_sinks=(),
            provided_interfaces=(WATER_INTERFACE,),
            required_interfaces=(
                ("atmospheric_state", VersionRange.parse(">=1.0.0,<2.0.0")),
                ("surface", VersionRange.parse(">=1.0.0,<2.0.0")),
            ),
            uncertainty_outputs=(
                "evaporation_coefficient_parameter",
                "missing_vapour_transport_structural",
                "missing_cloud_storage_structural",
                "latent_heat_constant_approximation",
            ),
            correctness_tests=(
                "test_evaporation_conserves_water_and_carries_its_latent_energy",
                "test_evaporation_stops_when_the_air_is_saturated",
                "test_condensation_conserves_water_and_energy_together",
                "test_moist_column_budgets_close_for_water_and_energy",
                "test_a_dry_column_with_no_evaporation_is_inert",
            ),
            known_limitations=(
                "No vapour transport, so humidity aloft comes only from local "
                "processes. Evaporation enters the lowest layer and nothing carries it "
                "higher, so in equilibrium the whole column's vapour sits in that layer "
                "at saturation and every layer above is exactly dry.",
                "No cloud water, and therefore no cloud radiative effect whatsoever.",
                "Whether the vapour is radiatively active is the radiation module's "
                "configuration, not this module's. When it is not, this column has no "
                "water-vapour feedback at all. When it is, it still has almost none, and "
                "the reason is the transport limitation above rather than the radiation: "
                "vapour confined to the lowest layer emits within a few kelvin of the "
                "surface temperature, and an absorber that re-emits at the temperature "
                "it absorbed at traps nothing. The same 23 kg m^-2 gives a greenhouse "
                "effect of 63 W m^-2 confined to the lowest layer against 125 spread "
                "through the troposphere, and an amplification of exactly 1.00 against "
                "2.23. The feedback is not weakened here; it is absent.",
            ),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(n_layers)",
                memory_complexity="O(n_layers)",
                scaling_variables=("n_layers",),
                work_units_per_step=1.5,
                benchmark_configuration=f"{n_layers} layers, global scalar column",
                peak_memory_bytes=8 * n_layers + 4096,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
                lower_fidelity_alternatives=("dry_column_no_water",),
                checkpoint_cost_bytes=8 * n_layers,
            ),
        )
        super().__init__(contract, status=status)
        self.atmosphere_module_id = atmosphere_module_id
        self.surface_module_id = surface_module_id
        self.n_layers = n_layers
        self.evaporation_coefficient = evaporation_coefficient
        self.initial_relative_humidity = initial_relative_humidity
        self.latent_heat = (
            latent_heat
            if latent_heat is not None
            else value_of("water_latent_heat_vaporization_273")
        )

        self.vapour = np.zeros(n_layers)
        self._evaporation = 0.0
        self._condensation = np.zeros(n_layers)
        self._precipitation = 0.0
        self._latent_heating = np.zeros(n_layers)
        self._specific_humidity = np.zeros(n_layers)
        self._relative_humidity = np.zeros(n_layers)
        self._interval = (0.0, 0.0)
        self._stored_initial = 0.0
        self._budget: StateBudget | None = None
        self._layer_mass = np.ones(n_layers)
        self._layer_temperature = np.zeros(n_layers)

    # --- configuration ----------------------------------------------------
    def validate_configuration(self) -> None:
        super().validate_configuration()
        if self.evaporation_coefficient < 0.0:
            raise ConfigurationError(
                f"evaporation coefficient {self.evaporation_coefficient} cannot be negative"
            )
        if not 0.0 <= self.initial_relative_humidity <= 1.0:
            raise ConfigurationError("initial relative humidity must lie in [0, 1]")
        if self.latent_heat <= 0.0:
            raise ConfigurationError("the latent heat of vaporisation must be positive")

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        supplied = initial_state.get(VAPOUR_STATE)
        if supplied is None:
            # Left at zero until the first step, when the layer masses and
            # temperatures needed to convert a relative humidity into a mass are
            # available. Initialising from a guess would put a number in the
            # manifest that the run never used.
            self.vapour = np.zeros(self.n_layers)
        else:
            values = np.asarray(supplied, dtype=float)
            self.vapour = (
                np.full(self.n_layers, float(values)) if values.ndim == 0 else np.array(values)
            )
            if self.vapour.size != self.n_layers:
                raise ConfigurationError(
                    f"initial vapour has {self.vapour.size} layers, expected {self.n_layers}"
                )
        if np.any(self.vapour < 0.0):
            raise ConfigurationError("vapour mass cannot be negative")
        self.mark_initialized()

    # --- physics ----------------------------------------------------------
    @property
    def precipitable_water(self) -> float:
        """Total column vapour in kg m^-2, the quantity WAT-001 conserves."""
        return float(np.sum(self.vapour))

    def latent_energy(self, vapour: np.ndarray | None = None) -> float:
        """``L_v`` times the column vapour, in J m^-2."""
        w = self.vapour if vapour is None else vapour
        return self.latent_heat * float(np.sum(w))

    def compute_tendencies(self) -> list[TendencyContribution]:
        ctx = self.context
        temperature = self.input_field(ATMOSPHERE_TEMPERATURE_STATE).require_units("K").values
        layer_mass = self.input_field("atmosphere_layer_mass").require_units("kg m^-2").values
        pressure = self.input_field("atmosphere_layer_pressure").require_units("Pa").values
        heat_capacity = (
            self.input_field("atmosphere_layer_heat_capacity").require_units("J m^-2 K^-1").values
        )
        surface_temperature = self.input_field("surface_temperature").require_units("K").scalar()

        self._layer_mass = np.asarray(layer_mass, dtype=float)
        self._layer_temperature = np.asarray(temperature, dtype=float)
        self._interval = (ctx.time, ctx.time + ctx.dt)
        self._stored_initial = self.latent_energy()

        saturation = saturation_specific_humidity(temperature, pressure)
        specific_humidity = self.vapour / self._layer_mass

        # Evaporation into the lowest layer only: with no vertical transport there
        # is no mechanism to carry vapour higher, and pretending otherwise would
        # be transport by another name.
        surface_saturation = float(
            saturation_specific_humidity(np.array([surface_temperature]), np.array([pressure[0]]))[
                0
            ]
        )
        deficit = max(0.0, surface_saturation - float(specific_humidity[0]))
        self._evaporation = self.evaporation_coefficient * deficit

        provisional = np.array(self.vapour)
        provisional[0] += self._evaporation * ctx.dt

        # Saturation adjustment: everything above saturation condenses and leaves.
        capacity = saturation * self._layer_mass
        self._condensation = np.maximum(0.0, provisional - capacity)
        self._precipitation = float(np.sum(self._condensation)) / ctx.dt
        self._latent_heating = self.latent_heat * self._condensation / ctx.dt

        final = provisional - self._condensation
        self._specific_humidity = final / self._layer_mass
        self._relative_humidity = np.divide(
            self._specific_humidity,
            saturation,
            out=np.zeros_like(saturation),
            where=saturation > 0.0,
        )

        return [
            TendencyContribution(
                state_id=VAPOUR_STATE,
                provider_module_id=self.module_id,
                label="evaporation_minus_condensation",
                value=(final - self.vapour) / ctx.dt,
                units="kg m^-2",
                equation_id=EQUATION_IDS["condensation"],
            ),
            # Latent heat released into the layer where it condensed. The column
            # owns its temperature; this is a labeled contribution to it, exactly
            # as the radiative and convective ones are.
            TendencyContribution(
                state_id=ATMOSPHERE_TEMPERATURE_STATE,
                provider_module_id=self.module_id,
                label="latent_heat_release",
                value=self._latent_heating / heat_capacity,
                units="K s^-1",
                equation_id=EQUATION_IDS["condensation"],
            ),
        ]

    def advance_state(self, dt: float, contributions: list[TendencyContribution]) -> None:
        initial = np.array(self.vapour, dtype=float)
        total = np.zeros(self.n_layers)
        for contribution in contributions:
            total = total + np.broadcast_to(contribution.as_array(), (self.n_layers,))
        self.vapour = np.maximum(initial + total * dt, 0.0)
        self._budget = StateBudget(
            state_id=VAPOUR_STATE,
            initial_value=initial,
            final_value=np.array(self.vapour, dtype=float),
            dt=dt,
            contributions=list(contributions),
        )

    def state_budgets(self) -> dict[str, Any]:
        return {VAPOUR_STATE: self._budget} if self._budget else {}

    # --- outputs ----------------------------------------------------------
    def emit_state(self) -> dict[str, Field]:
        provenance = Provenance(producer=self.module_id, derivation=EQUATION_IDS["latent_energy"])
        timestamp = self._context.time if self._context else 0.0
        return {
            "vapour_bearing_layer_temperature": Field(
                field_id="vapour_bearing_layer_temperature",
                physical_name="layer temperature where the column holds vapour",
                symbol="T_w",
                # NaN where there is no vapour, which is a statement rather than
                # a gap: the latent energy sums over w_i, so a layer holding
                # none contributes no term, and its temperature is not one the
                # relation was evaluated at. The masked profile is published
                # rather than a vapour-weighted mean because a declared range is
                # about the extent -- a mean would hide a hot layer and a cold
                # one by averaging them into the middle.
                values=np.where(
                    np.asarray(self.vapour, dtype=float) > 0.0,
                    np.asarray(self._layer_temperature, dtype=float),
                    np.nan,
                ),
                units="K",
                quantity_kind=QuantityKind.DIAGNOSTIC,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                timestamp=timestamp,
                spatial_support=SpatialSupport.LAYER_AVERAGE,
                dimensions=("layer",),
                coordinate_system="pressure_layers",
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                # The mask is the content of this field, so it declares the
                # policy DAT-007 already provides for it rather than smuggling
                # the NaNs past a default that exists to catch them.
                missing_value_policy=MissingValuePolicy.MASKED,
                provenance=provenance,
            ),
            VAPOUR_STATE: Field(
                field_id=VAPOUR_STATE,
                physical_name="water vapour mass per unit area in each layer",
                symbol="w",
                values=np.array(self.vapour),
                units="kg m^-2",
                quantity_kind=QuantityKind.EXTENSIVE_STATE,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                timestamp=timestamp,
                spatial_support=SpatialSupport.LAYER_AVERAGE,
                dimensions=("layer",),
                coordinate_system="pressure_layers",
                phase_or_species="water_vapour",
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                provenance=provenance,
                uncertainty=Uncertainty(
                    method="none",
                    unresolved=("initial_condition", "evaporation_coefficient_parameter"),
                ),
            ),
        }

    def emit_outputs(self) -> dict[str, Field]:
        ctx = self.context
        provenance = Provenance(producer=self.module_id, derivation=EQUATION_IDS["condensation"])
        scalar_flux = {
            "quantity_kind": QuantityKind.FLUX,
            "time_semantics": TimeSemantics.INTERVAL_AVERAGE,
            "time_interval": ctx.dt,
            "timestamp": ctx.time,
            "spatial_support": SpatialSupport.GLOBAL_SCALAR,
            "provenance": provenance,
            "sign_convention": SignConvention.POSITIVE_DEFINITE,
        }
        profile = {
            "quantity_kind": QuantityKind.DIAGNOSTIC,
            "time_semantics": TimeSemantics.END_OF_STEP,
            "timestamp": ctx.time,
            "spatial_support": SpatialSupport.LAYER_AVERAGE,
            "dimensions": ("layer",),
            "coordinate_system": "pressure_layers",
            "provenance": provenance,
        }
        return {
            "surface_evaporation_flux": Field(
                field_id="surface_evaporation_flux",
                physical_name="evaporative water mass flux from the surface",
                symbol="E",
                values=np.array(self._evaporation),
                units="kg m^-2 s^-1",
                **scalar_flux,  # type: ignore[arg-type]
            ),
            "precipitation_flux": Field(
                field_id="precipitation_flux",
                physical_name="precipitation water mass flux to the surface",
                symbol="P",
                values=np.array(self._precipitation),
                units="kg m^-2 s^-1",
                **scalar_flux,  # type: ignore[arg-type]
            ),
            "surface_latent_heat_flux": Field(
                field_id="surface_latent_heat_flux",
                physical_name="latent heat flux leaving the surface with evaporation",
                symbol="LE",
                values=np.array(self.latent_heat * self._evaporation),
                units="W m^-2",
                **scalar_flux,  # type: ignore[arg-type]
            ),
            "atmosphere_precipitable_water": Field(
                field_id="atmosphere_precipitable_water",
                physical_name="total column water vapour",
                symbol="W",
                values=np.array(self.precipitable_water),
                units="kg m^-2",
                quantity_kind=QuantityKind.EXTENSIVE_STATE,
                time_semantics=TimeSemantics.END_OF_STEP,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.DOMAIN_INTEGRAL,
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                provenance=provenance,
            ),
            "latent_heating_rate": Field(
                field_id="latent_heating_rate",
                physical_name="latent heat released in each layer by condensation",
                symbol="LC",
                values=np.array(self._latent_heating),
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
            ),
            "atmosphere_specific_humidity": Field(
                field_id="atmosphere_specific_humidity",
                physical_name="specific humidity of each layer",
                symbol="q",
                values=np.array(self._specific_humidity),
                units="kg kg^-1",
                **profile,  # type: ignore[arg-type]
            ),
            "atmosphere_relative_humidity": Field(
                field_id="atmosphere_relative_humidity",
                physical_name="relative humidity of each layer",
                symbol="RH",
                values=np.array(self._relative_humidity),
                units="1",
                **profile,  # type: ignore[arg-type]
            ),
        }

    def report_transfers(self) -> tuple[list[Transfer], list[Transfer]]:
        """Two quantities across two boundaries, recorded from this side.

        Water enters from the surface and leaves to it; energy enters from the
        surface with the evaporated vapour and leaves to the atmosphere when that
        vapour condenses. Precipitation carries water but no energy, which is the
        declared approximation of neglecting its sensible heat.
        """
        start, end = self._interval
        dt = end - start
        common = {
            "start_time": start,
            "end_time": end,
            "spatial_mapping": "global_mean_over_planetary_surface",
        }
        condensed = float(np.sum(self._condensation))

        imports = [
            Transfer(
                transfer_id=f"evaporation_water@{start:g}",
                sender_module=self.surface_module_id,
                receiver_module=self.module_id,
                quantity=ConservedQuantity.WATER,
                amount=self._evaporation * dt,
                units="kg m^-2",
                channel="evaporation_water",
                species="water_vapour",
                equation_id=EQUATION_IDS["surface_evaporation"],
                **common,  # type: ignore[arg-type]
            ),
            Transfer(
                transfer_id=f"evaporation_latent_energy@{start:g}",
                sender_module=self.surface_module_id,
                receiver_module=self.module_id,
                quantity=ConservedQuantity.ENERGY,
                amount=self.latent_heat * self._evaporation * dt,
                units="J m^-2",
                channel="evaporation_latent_energy",
                equation_id=EQUATION_IDS["latent_energy"],
                **common,  # type: ignore[arg-type]
            ),
        ]
        exports = [
            Transfer(
                transfer_id=f"precipitation_water@{start:g}",
                sender_module=self.module_id,
                receiver_module=self.surface_module_id,
                quantity=ConservedQuantity.WATER,
                amount=condensed,
                units="kg m^-2",
                channel="precipitation_water",
                species="liquid_water",
                equation_id=EQUATION_IDS["condensation"],
                **common,  # type: ignore[arg-type]
            ),
            Transfer(
                transfer_id=f"latent_heat_release@{start:g}",
                sender_module=self.module_id,
                receiver_module=self.atmosphere_module_id,
                quantity=ConservedQuantity.ENERGY,
                amount=self.latent_heat * condensed,
                units="J m^-2",
                channel="latent_heat_release",
                equation_id=EQUATION_IDS["condensation"],
                **common,  # type: ignore[arg-type]
            ),
        ]
        return (exports, imports)

    def report_conservation(self) -> list[ModuleBudget]:
        return [
            ModuleBudget(
                module_id=self.module_id,
                quantity=ConservedQuantity.WATER,
                units="kg m^-2",
                stored_initial=self._stored_initial / self.latent_heat,
                stored_final=self.precipitable_water,
            ),
            ModuleBudget(
                module_id=self.module_id,
                quantity=ConservedQuantity.ENERGY,
                units="J m^-2",
                stored_initial=self._stored_initial,
                stored_final=self.latent_energy(),
            ),
        ]

    def report_uncertainty(self) -> UncertaintyReport:
        return UncertaintyReport(
            module_id=self.module_id,
            quantified={
                "precipitable_water_kg_m2": self.precipitable_water,
                "evaporation_coefficient_kg_m2_s": self.evaporation_coefficient,
            },
            method="none propagated at this fidelity",
            distributional_assumptions="unresolved",
            correlations="unresolved",
            unresolved=(
                "evaporation_coefficient_parameter_uncertainty",
                "structural_uncertainty_of_omitting_vapour_transport",
                "structural_uncertainty_of_omitting_cloud_storage",
                "constant_latent_heat_approximation",
            ),
        )

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "n_layers": self.n_layers,
            "evaporation_coefficient_kg_m2_s": self.evaporation_coefficient,
            "latent_heat_J_kg": self.latent_heat,
            "initial_relative_humidity": self.initial_relative_humidity,
            "condensate_fate": "precipitates_immediately",
            "vapour_transport": "none",
        }

    def checkpoint(self) -> dict[str, Any]:
        return {VAPOUR_STATE: np.array(self.vapour)}

    def restore(self, checkpoint: dict[str, Any]) -> None:
        self.vapour = np.array(checkpoint[VAPOUR_STATE], dtype=float)


class PrescribedHumidityModule(Module):
    """A humidity profile imposed at fixed relative humidity (WAT-002, HYP-005).

    **This is a control, not a better water module.** It exists because two
    untested things were entangled: whether the vapour absorption relation gives
    the right feedback when the humidity distribution is realistic, and whether
    this model can produce a realistic humidity distribution. ``WaterModule``
    has no vertical transport, so it cannot -- all its vapour ends up in the
    lowest layer -- and a feedback measured through it tells you nothing about
    the radiation. Imposing the distribution pins the second question so the
    first can be answered on its own.

    So the *distribution* is prescribed and the *amount* is not. Relative
    humidity is held fixed and Clausius-Clapeyron does the rest, which is the
    idealisation Manabe and Wetherald used in 1967 and for the same reason. The
    feedback that results is therefore not a prediction of this configuration:
    its size follows largely from the assumption that relative humidity is
    conserved, and a claim registered against it would be confirming what was
    put in. What it *is* good for is a reference number -- the feedback the
    radiation produces when the humidity is right by construction -- against
    which a real transport scheme can later be judged.

    ``respond_to_temperature=False`` freezes the vapour at the profile the first
    step implies, which removes the feedback while changing nothing else. That
    is the null control for the control. Start it from the responding run's own
    equilibrium: freezing at a cold initial guess makes the control carry a
    spin-up transient the experiment does not, and their difference then is not
    the feedback.

    **The reference numbers this module exists to produce**, measured in an
    eight-layer convecting column at 0.7 relative humidity::

        vapour free      dT_2x 4.394 K   lambda 0.842 W m^-2 K^-1
        vapour frozen    dT_2x 1.768 K   lambda 2.093
        amplification    2.49x           feedback +1.251 W m^-2 K^-1

    The amplification agrees with the 2.51x a fixed-profile radiative
    calculation gives, which is the check that says the coupled column is doing
    what the absorption relation implies. The feedback is smaller than that
    calculation's +2.29 because this column convects: convective adjustment
    supplies something close to a lapse-rate feedback, and published estimates
    put water vapour and lapse rate together near +1.2 W m^-2 K^-1.

    Two caveats on using these as a target. The base state is 293.0 K against an
    observed 288, and sensitivity rises with base temperature because
    Clausius-Clapeyron is exponential, so a transport scheme that lands at a
    different temperature is not being compared like for like. And dT_2x of
    4.39 K is high for a column with no cloud and no surface-albedo feedback --
    higher than the real system's, which has both. Neither is tuned away here.

    Owns no state. It publishes ``atmosphere_water_vapour`` the way the
    prescribed cloud and greenhouse modules publish fields they do not own, so a
    configuration runs this module or ``WaterModule``, never both.
    """

    def __init__(
        self,
        *,
        module_id: str = "water",
        radiation_module_id: str = "radiation",
        atmosphere_module_id: str = "atmosphere",
        n_layers: int = 8,
        relative_humidity: float = 0.7,
        respond_to_temperature: bool = True,
        status: ModuleStatus = ModuleStatus.PRESCRIBED,
    ) -> None:
        contract = ModuleContract(
            module_id=module_id,
            module_name="prescribed fixed-relative-humidity vapour profile",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.PRESCRIBED,
            physical_basis=(
                "Relative humidity held fixed, so the vapour a layer carries is a fixed "
                "fraction of what it can carry. The saturation value comes from "
                "Clausius-Clapeyron, so the vapour amount rises with temperature at "
                "about seven percent per kelvin without that rate being imposed."
            ),
            governing_equation_ids=(
                EQUATION_IDS["clausius_clapeyron"],
                EQUATION_IDS["vapour_band_absorption"],
            ),
            assumptions=(
                "relative humidity is uniform in height and constant in time",
                "no water mass budget: this module moves no water and closes none",
                "no latent heat, so the energy a real phase change carries is absent",
            ),
            validity_domain=(
                "A reference distribution for measuring the radiative water-vapour "
                "feedback. Nothing about the water cycle, precipitation, latent heating "
                "or the humidity distribution itself is in scope, because all of those "
                "are imposed here rather than computed."
            ),
            required_inputs=(
                ATMOSPHERE_TEMPERATURE_STATE,
                "atmosphere_layer_mass",
                "atmosphere_layer_pressure",
            ),
            produced_outputs=(
                VAPOUR_STATE,
                "atmosphere_specific_humidity",
                "atmosphere_relative_humidity",
                "atmosphere_precipitable_water",
            ),
            conserved_quantities=(),
            uncertainty_outputs=("prescribed_relative_humidity_boundary_condition",),
            correctness_tests=(
                "test_the_prescribed_profile_follows_clausius_clapeyron",
                "test_holding_the_vapour_removes_the_feedback_and_nothing_else",
                "test_the_prescribed_profile_is_republished_every_step",
            ),
            known_limitations=(
                "Prescribed, so the humidity distribution is an input. The feedback that "
                "follows is not a prediction of this configuration: its magnitude comes "
                "largely from the assumption that relative humidity is conserved, so no "
                "claim about the water-vapour feedback can be registered against it.",
                "No water or latent-energy budget at all. Vapour appears and disappears "
                "with the temperature and nothing accounts for it, which is why this "
                "module conserves nothing and says so.",
                "Relative humidity is uniform in height. The real profile falls through "
                "the mid-troposphere, so this overstates the vapour aloft, which is "
                "where the feedback is most sensitive.",
            ),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(n_layers)",
                memory_complexity="O(n_layers)",
                scaling_variables=("n_layers",),
                work_units_per_step=0.3,
                benchmark_configuration=f"{n_layers} layers, global scalar column",
                peak_memory_bytes=8 * n_layers + 4096,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
                lower_fidelity_alternatives=("dry_column_no_water",),
            ),
        )
        super().__init__(contract, status=status)
        self.radiation_module_id = radiation_module_id
        self.atmosphere_module_id = atmosphere_module_id
        self.n_layers = n_layers
        self.relative_humidity = relative_humidity
        self.respond_to_temperature = respond_to_temperature
        self.vapour = np.zeros(n_layers)
        self._specific_humidity = np.zeros(n_layers)
        self._held: np.ndarray | None = None

    def validate_configuration(self) -> None:
        super().validate_configuration()
        if not 0.0 <= self.relative_humidity <= 1.0:
            raise ConfigurationError(
                f"relative humidity {self.relative_humidity} is not a fraction between zero and one"
            )
        if self.n_layers < 1:
            raise ConfigurationError("a humidity profile needs at least one layer")

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        self.vapour = np.zeros(self.n_layers)
        self._specific_humidity = np.zeros(self.n_layers)
        self._saturation = np.ones(self.n_layers)
        self._relative_humidity_profile = np.full(self.n_layers, float(self.relative_humidity))
        self._held = None
        self.mark_initialized()

    def compute_tendencies(self) -> list[Any]:
        temperature = self.input_field(ATMOSPHERE_TEMPERATURE_STATE).require_units("K").values
        layer_mass = self.input_field("atmosphere_layer_mass").require_units("kg m^-2").values
        pressure = self.input_field("atmosphere_layer_pressure").require_units("Pa").values

        saturation = saturation_specific_humidity(
            np.asarray(temperature, dtype=float), np.asarray(pressure, dtype=float)
        )
        self._saturation = saturation
        applied = self._applied_relative_humidity(
            np.asarray(temperature, dtype=float),
            np.asarray(pressure, dtype=float),
            np.asarray(layer_mass, dtype=float),
        )
        self._relative_humidity_profile = np.broadcast_to(
            np.asarray(applied, dtype=float), (self.n_layers,)
        ).copy()
        self._specific_humidity = self._relative_humidity_profile * saturation
        vapour = self._specific_humidity * np.asarray(layer_mass, dtype=float)

        if self.respond_to_temperature:
            self.vapour = vapour
        else:
            # The null control: frozen at what the first step implied, so the
            # column warms without its vapour being allowed to answer.
            if self._held is None:
                self._held = np.array(vapour)
            self.vapour = np.array(self._held)
            self._specific_humidity = self.vapour / np.asarray(layer_mass, dtype=float)
        return []

    def _applied_relative_humidity(
        self, temperature: np.ndarray, pressure: np.ndarray, layer_mass: np.ndarray
    ) -> np.ndarray | float:
        """The relative humidity to impose. A single number here, by definition.

        The hook exists so a subclass can *predict* a profile without
        reimplementing the Clausius-Clapeyron arithmetic that turns one into
        vapour.
        """
        return self.relative_humidity

    def precipitable_water(self) -> float:
        return float(np.sum(self.vapour))

    #: What the published fields say they came from. Overridden by subclasses
    #: that derive the humidity instead of imposing it, so the provenance never
    #: describes a different module than the one that produced the numbers.
    derivation_equation_id = EQUATION_IDS["clausius_clapeyron"]
    provenance_note = "prescribed distribution; relative humidity is imposed, not computed"

    def _fields(self, timestamp: float) -> dict[str, Field]:
        provenance = Provenance(
            producer=self.module_id,
            derivation=self.derivation_equation_id,
            notes=self.provenance_note,
        )
        profile = {
            "quantity_kind": QuantityKind.BOUNDARY_VALUE,
            "time_semantics": TimeSemantics.INSTANTANEOUS,
            "timestamp": timestamp,
            "spatial_support": SpatialSupport.LAYER_AVERAGE,
            "dimensions": ("layer",),
            "coordinate_system": "pressure_layers",
            "sign_convention": SignConvention.POSITIVE_DEFINITE,
            "provenance": provenance,
        }
        return {
            VAPOUR_STATE: Field(
                field_id=VAPOUR_STATE,
                physical_name="water vapour mass per unit area in each layer",
                symbol="w",
                values=np.array(self.vapour),
                units="kg m^-2",
                phase_or_species="water_vapour",
                **profile,  # type: ignore[arg-type]
            ),
            "atmosphere_specific_humidity": Field(
                field_id="atmosphere_specific_humidity",
                physical_name="specific humidity",
                symbol="q",
                values=np.array(self._specific_humidity),
                units="1",
                **profile,  # type: ignore[arg-type]
            ),
            "atmosphere_relative_humidity": Field(
                field_id="atmosphere_relative_humidity",
                physical_name="relative humidity",
                symbol="RH",
                # Frozen, the relative humidity is no longer the imposed one:
                # the vapour is held while the temperature moves, so it has to
                # be recovered from the vapour that is actually there. This
                # field used to carry the specific humidity in that branch,
                # which is a different quantity of the same nominal units.
                values=np.array(self._relative_humidity_profile)
                if self.respond_to_temperature
                else np.divide(
                    self._specific_humidity,
                    self._saturation,
                    out=np.zeros(self.n_layers),
                    where=self._saturation > 0.0,
                ),
                units="1",
                **profile,  # type: ignore[arg-type]
            ),
            "atmosphere_precipitable_water": Field(
                field_id="atmosphere_precipitable_water",
                physical_name="column integrated water vapour",
                symbol="W",
                values=np.array(self.precipitable_water()),
                units="kg m^-2",
                quantity_kind=QuantityKind.DIAGNOSTIC,
                time_semantics=TimeSemantics.END_OF_STEP,
                timestamp=timestamp,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                provenance=provenance,
            ),
        }

    def emit_state(self) -> dict[str, Field]:
        return self._fields(self._context.time if self._context else 0.0)

    def emit_outputs(self) -> dict[str, Field]:
        # Republished every step, because the profile follows a temperature that
        # changes. A PRESCRIBED module whose emit_state ran only at
        # initialisation would silently freeze the feedback it exists to supply.
        return self._fields(self.context.time)

    def report_uncertainty(self) -> UncertaintyReport:
        return UncertaintyReport(
            module_id=self.module_id,
            method="none; the humidity distribution is imposed",
            distributional_assumptions="unresolved",
            correlations="unresolved",
            unresolved=(
                "structural_uncertainty_of_a_prescribed_uniform_relative_humidity",
                "the observed profile falls through the mid-troposphere and this does not",
            ),
        )

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "relative_humidity": self.relative_humidity,
            "respond_to_temperature": self.respond_to_temperature,
            "n_layers": self.n_layers,
            "water_budget": "not_represented",
            "latent_energy": "not_represented",
        }


class EntrainingPlumeHumidityModule(PrescribedHumidityModule):
    """Relative humidity predicted by a bulk plume rather than chosen (WAT-002, CLM-003).

    The difference from ``PrescribedHumidityModule`` is one number and its
    provenance. That module holds relative humidity at 0.7 because Manabe and
    Wetherald did, and says plainly that no feedback claim may be registered
    against it, since the feedback follows from the assumption. This one
    computes the same quantity from the environment's water budget::

        RH = (a + alpha) / (1 + a)

    where ``a`` is the dimensionless bulk-plume parameter and ``alpha`` the
    fraction of condensate that evaporates instead of falling out. At the
    measured values -- ``a = 0.2`` for Earth's tropics, and ``alpha = 0.67`` from
    precipitation efficiencies near a third -- that is **0.725**, and nothing
    about the temperature record went into it.

    Why this is worth the trouble over prescribing 0.4 directly: the same ``a``
    also fixes the lapse rate an entraining plume leaves behind
    (``EQ-ENTRAINING-LAPSE-RATE``). Two observables, one parameter, so the pair
    is over-determined and can be wrong in a way a chosen humidity cannot be.
    That is the whole reason for preferring this to a second prescribed number.

    **What it still does not do.** It moves no water and closes no water budget,
    exactly like the module it derives from -- it predicts the *distribution*,
    not the cycle. It gives one humidity for the whole column, so the observed C
    shape is absent. And it describes a convecting column while this framework's
    is a global mean, so it should be expected to run *moist*: the global-mean
    free troposphere is far drier than a radiative-convective one, for reasons
    -- subtropical subsidence, last-saturation dynamics -- that no bulk plume
    contains.

    Set ``respond_to_temperature=False`` for the null control, which freezes the
    vapour and removes the feedback while changing nothing else.
    """

    derivation_equation_id = EQUATION_IDS["bulk_plume_humidity"]
    provenance_note = (
        "relative humidity predicted from the bulk-plume water budget; "
        "no temperature record enters the closure"
    )

    def __init__(
        self,
        *,
        module_id: str = "water",
        radiation_module_id: str = "radiation",
        atmosphere_module_id: str = "atmosphere",
        n_layers: int = 8,
        bulk_plume_parameter: float = BULK_PLUME_PARAMETER,
        evaporating_fraction: float = EVAPORATING_CONDENSATE_FRACTION,
        vertical_structure: bool = False,
        respond_to_temperature: bool = True,
        status: ModuleStatus = ModuleStatus.ACTIVE,
    ) -> None:
        relative_humidity = equilibrium_relative_humidity(
            bulk_plume_parameter, evaporating_fraction
        )
        contract = ModuleContract(
            module_id=module_id,
            module_name="bulk-plume predicted relative humidity profile",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.REDUCED,
            physical_basis=(
                "Subsidence dries the air between clouds and detrainment of saturated "
                "plume air moistens it; the balance fixes the relative humidity. The "
                "vapour amount then follows from Clausius-Clapeyron, so it rises with "
                "temperature at about seven percent per kelvin without that rate, or "
                "the humidity itself, being imposed."
            ),
            governing_equation_ids=(
                EQUATION_IDS["bulk_plume_humidity"],
                EQUATION_IDS["clausius_clapeyron"],
                EQUATION_IDS["vapour_band_absorption"],
            ),
            assumptions=(
                "radiative-convective equilibrium, so the water budget is steady",
                "one bulk plume stands for the whole cloud ensemble",
                (
                    "the convective mass flux is set by clear-sky radiative cooling, "
                    "and the detrainment rate by its vertical gradient"
                    if vertical_structure
                    else "constant mass flux with height, which makes the humidity uniform"
                ),
                "no water mass budget: this module moves no water and closes none",
                "no latent heat, so the energy a real phase change carries is absent",
            ),
            validity_domain=(
                "Convecting columns where cumulus mass flux rather than large-scale "
                "advection sets the free-tropospheric humidity. The water cycle itself "
                "-- precipitation, latent heating, the surface reservoir -- is out of "
                "scope, because this module predicts a distribution and not a flux."
            ),
            required_inputs=(
                ATMOSPHERE_TEMPERATURE_STATE,
                "atmosphere_layer_mass",
                "atmosphere_layer_pressure",
            )
            + (
                # The mass-flux constraint is what turns the closure from one
                # number into a profile, and the cooling that sets it is the
                # column's own.
                ("layer_radiative_heating", "atmosphere_layer_heat_capacity")
                if vertical_structure
                else ()
            ),
            produced_outputs=(
                VAPOUR_STATE,
                "atmosphere_specific_humidity",
                "atmosphere_relative_humidity",
                "atmosphere_precipitable_water",
            ),
            conserved_quantities=(),
            equation_variables=(
                (
                    EQUATION_IDS["bulk_plume_humidity"],
                    (
                        ("a", "bulk_plume_parameter", "1"),
                        ("alpha", "evaporating_condensate_fraction", "1"),
                    ),
                ),
                (
                    EQUATION_IDS["clausius_clapeyron"],
                    (("T", ATMOSPHERE_TEMPERATURE_STATE, "K"),),
                ),
            ),
            uncertainty_outputs=(
                "bulk_plume_parameter_uncertainty",
                "evaporating_condensate_fraction_uncertainty",
            ),
            correctness_tests=(
                "test_the_predicted_humidity_lands_in_the_observed_range",
                "test_the_plume_humidity_follows_clausius_clapeyron",
                "test_the_plume_module_derives_its_humidity_rather_than_taking_one",
            ),
            known_limitations=(
                "The evaporating condensate fraction comes from measured precipitation "
                "efficiencies spanning 19 to 38 percent, so the predicted humidity runs "
                "from 0.68 to 0.84. Any result that depends on it should be reported "
                "across that range and not at the midpoint.",
                "A convecting-column relation applied to a global-mean column. It "
                "reproduces the humidity of a radiative-convective column while the "
                "global-mean free troposphere is far drier, being dominated by "
                "subtropical subsidence no bulk plume represents, so this runs too "
                "moist by an amount the framework cannot estimate from within.",
                (
                    "The closure is a free-tropospheric relation and is applied to the "
                    "boundary layer as well, where cumulus mass flux is not what sets "
                    "the humidity; the lowest layers come out nearer saturation than "
                    "observed."
                    if vertical_structure
                    else "One humidity for the whole column, so the observed C-shaped "
                    "profile is absent and the vapour aloft, where the feedback is most "
                    "sensitive, is the least trustworthy part of it."
                ),
                "No water or latent-energy budget. The distribution is predicted; the "
                "cycle that would sustain it is not, which is why this module still "
                "conserves nothing and says so.",
            ),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(n_layers)",
                memory_complexity="O(n_layers)",
                scaling_variables=("n_layers",),
                work_units_per_step=0.3,
                benchmark_configuration=f"{n_layers} layers, global scalar column",
                peak_memory_bytes=8 * n_layers + 4096,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
                lower_fidelity_alternatives=("prescribed_fixed_relative_humidity",),
            ),
        )
        # Deliberately the grandparent's constructor: this module reuses every
        # bit of the prescribed module's machinery and none of its contract,
        # because the contract is the part that says where the number came from.
        Module.__init__(self, contract, status=status)
        self.radiation_module_id = radiation_module_id
        self.atmosphere_module_id = atmosphere_module_id
        self.n_layers = n_layers
        self.bulk_plume_parameter = bulk_plume_parameter
        self.evaporating_fraction = evaporating_fraction
        self.vertical_structure = vertical_structure
        self.relative_humidity = relative_humidity
        self.respond_to_temperature = respond_to_temperature
        self.vapour = np.zeros(n_layers)
        self._specific_humidity = np.zeros(n_layers)
        self._held: np.ndarray | None = None
        self._mass_flux = np.zeros(n_layers)

    def _applied_relative_humidity(
        self, temperature: np.ndarray, pressure: np.ndarray, layer_mass: np.ndarray
    ) -> np.ndarray | float:
        """The humidity profile the bulk plume predicts, height by height.

        The constant-mass-flux form gives one number for the column because it
        assumes ``epsilon = delta``. Dropping that assumption needs a second
        relation, and radiative-convective equilibrium supplies it: the
        environment subsides at whatever rate balances its own radiative
        cooling, so the convective mass flux is determined rather than chosen,
        and ``d_z M = (epsilon - delta) M`` then fixes the detrainment rate.

        What comes out is the observed C shape, and it comes out of the same
        parameter the constant-mass-flux form already used -- no new knob, and
        nothing calibrated against a temperature record.
        """
        if not self.vertical_structure:
            return self.relative_humidity

        heating = self.input_field("layer_radiative_heating").require_units("W m^-2").values
        capacity = (
            self.input_field("atmosphere_layer_heat_capacity").require_units("J m^-2 K^-1").values
        )
        cooling = -np.asarray(heating, dtype=float) / np.asarray(capacity, dtype=float)

        # Layer thickness from hydrostatic balance, so no surface pressure is
        # needed: layer_mass = rho dz with rho = p / (R T).
        thickness = layer_mass * _GAS_CONSTANT_DRY_AIR * temperature / pressure
        height = np.cumsum(thickness) - 0.5 * thickness

        lapse_rate = -np.gradient(temperature, height)
        self._mass_flux = convective_mass_flux(cooling, temperature, pressure, lapse_rate)
        gamma = water_vapour_lapse_rate(temperature, lapse_rate)

        # Where the profile is not decreasing with height the vapour lapse rate
        # goes negative and the plume equations have no solution -- there is
        # nothing there for a plume to do. Common during spin-up, and permanent
        # in a stratosphere. Those layers are excluded rather than clipped, so
        # the closure is never evaluated outside its own domain.
        convecting = (self._mass_flux > 0.0) & (gamma > 0.0)
        usable = np.where(convecting, gamma, 1.0)

        entrainment = fractional_entrainment_rate(
            usable, self.bulk_plume_parameter, self.evaporating_fraction
        )
        detrainment = detrainment_rate(self._mass_flux, entrainment, height)
        humidity = bulk_plume_relative_humidity(
            entrainment=entrainment,
            detrainment=detrainment,
            vapour_lapse_rate=usable,
            evaporating_fraction=self.evaporating_fraction,
        )
        # Outside the convecting layer the closure has nothing to say and the
        # constant-mass-flux value stands in. That is an admission, not a
        # result: see the module's known limitations.
        return np.clip(np.where(convecting, humidity, self.relative_humidity), 0.0, 1.0)

    def implied_entrainment_rate(self, temperature: float, lapse_rate: float) -> float:
        """The dimensional entrainment rate this closure stands for, in m^-1.

        Reported rather than used, because it is what the cloud-resolving
        literature quotes and therefore the number a reader can check the
        parameter choice against. At present it comes out several times smaller
        than the rates those simulations diagnose, which is recorded in the
        equation registry as a known failure mode rather than adjusted away.
        """
        gamma = water_vapour_lapse_rate(np.array(temperature), lapse_rate)
        return float(
            fractional_entrainment_rate(gamma, self.bulk_plume_parameter, self.evaporating_fraction)
        )

    def report_uncertainty(self) -> UncertaintyReport:
        low, high = EVAPORATING_FRACTION_RANGE
        return UncertaintyReport(
            module_id=self.module_id,
            quantified={
                "relative_humidity": self.relative_humidity,
                "relative_humidity_low_alpha": equilibrium_relative_humidity(
                    self.bulk_plume_parameter, low
                ),
                "relative_humidity_high_alpha": equilibrium_relative_humidity(
                    self.bulk_plume_parameter, high
                ),
            },
            method="range over the declared evaporating-fraction assumption",
            distributional_assumptions="unresolved; a range is asserted, not a distribution",
            correlations="unresolved",
            unresolved=(
                "the evaporating condensate fraction is assumed, not measured",
                "a tropical closure is applied to a global-mean column",
                "the humidity is uniform in height and the observed profile is not",
            ),
        )

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "bulk_plume_parameter": self.bulk_plume_parameter,
            "evaporating_condensate_fraction": self.evaporating_fraction,
            "relative_humidity": self.relative_humidity,
            "relative_humidity_origin": "predicted_from_the_bulk_plume_water_budget",
            "respond_to_temperature": self.respond_to_temperature,
            "n_layers": self.n_layers,
            "water_budget": "not_represented",
            "latent_energy": "not_represented",
        }
