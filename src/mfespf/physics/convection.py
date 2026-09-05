"""Dry convective adjustment.

Charter 9.4 and section 19 stage 3; requirements ATM-002, MOD-011, STA-006,
CON-009.

This module exists because the Stage 2 column was wrong in a way it declared:
`EQ-LAYER-ENERGY-BALANCE` lists a super-adiabatic lapse rate among its known
failure modes. Convective adjustment removes such instability the only way a
model without resolved motion can -- by asserting that where the profile would be
unstable, mixing has already restored it.

A finding worth recording, because it was not what was expected. Measured in the
Stage 2 configurations, radiative equilibrium alone is *already* dry-stable:
adding convective adjustment to it changed nothing at all, to the last digit.
The instability the textbooks describe lives in the 27 K jump between the
surface and the air above it, and a column whose only surface-air coupling is
radiative cannot act on that jump. The sensible heat flux below is therefore not
an embellishment; without it convective adjustment is inert, and with it the
surface cools by about 5 K while the lowest layers warm by nearly 20 K.

It is the framework's first case of **two modules contributing to one state**.
Radiation and convection both submit labeled tendencies for the atmospheric
layer temperatures; the column owns the state and commits their sum (A1.3). The
labels survive into the per-state budget, so a reader can see how much of a
layer's warming was radiative and how much was mixing.

Convection moves energy between layers and across no module boundary, so its
contribution to the column's *total* energy is exactly zero. That is checked,
not assumed: a convection scheme that quietly created energy would otherwise be
indistinguishable from one that did not.
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
from mfespf.core.errors import ConfigurationError
from mfespf.core.quantities import (
    Field,
    Provenance,
    QuantityKind,
    SignConvention,
    SpatialSupport,
    TimeSemantics,
)
from mfespf.core.versioning import InterfaceContract, Version, VersionRange
from mfespf.physics.atmosphere import ATMOSPHERE_TEMPERATURE_STATE
from mfespf.physics.equations import EQUATION_IDS
from mfespf.physics.surface import SURFACE_TEMPERATURE_STATE
from mfespf.physics.thermodynamics import (
    convective_adjustment_with_groups,
    potential_temperature,
)
from mfespf.physics.water import VAPOUR_STATE
from mfespf.registries.states import TendencyContribution

__all__ = ["CONVECTION_INTERFACE", "ConvectiveAdjustmentModule"]

V1 = Version(1, 0, 0)

CONVECTION_INTERFACE = InterfaceContract(
    interface_id="convective_adjustment",
    version=V1,
    physical_semantics={
        "tendency_units": "K s^-1",
        "layer_ordering": "index_0_is_the_lowest_layer",
        "stability_variable": "dry_potential_temperature",
        "energy_role": "redistributes_within_the_column_and_adds_nothing",
    },
)


class ConvectiveAdjustmentModule(Module):
    """Restores a dry-stable profile, conserving the column's enthalpy."""

    def __init__(
        self,
        *,
        module_id: str = "convection",
        atmosphere_module_id: str = "atmosphere",
        surface_module_id: str = "surface",
        sensible_heat_coefficient: float = 10.0,
        transport_vapour: bool = False,
        reference_pressure: float | None = None,
        status: ModuleStatus = ModuleStatus.ACTIVE,
    ) -> None:
        contract = ModuleContract(
            module_id=module_id,
            module_name="surface sensible heat exchange and dry convective adjustment",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.REDUCED,
            physical_basis=(
                "Turbulent exchange carries sensible heat from the surface into the air "
                "above it; where the resulting dry potential temperature decreases with "
                "height the column is unstable, and convection is assumed to have mixed "
                "such layers to neutrality."
            ),
            governing_equation_ids=(
                EQUATION_IDS["potential_temperature"],
                EQUATION_IDS["convective_adjustment"],
                EQUATION_IDS["sensible_heat"],
            ),
            assumptions=(
                "the surface exchange coefficient absorbs wind speed and stability",
                "convective mixing is instantaneous on the coupling timestep",
                "mixing is complete within an unstable group and absent outside it",
                "dry adiabatic neutrality is the target, with no moist correction",
                "no momentum, moisture or tracer transport accompanies the heat mixing",
            ),
            validity_domain=(
                "Dry columns, and moist columns only where the answer does not depend on "
                "condensation inside the mixed layer. It targets the dry adiabat, so it "
                "over-stabilises a saturated troposphere, whose real lapse rate is the "
                "smaller moist one."
            ),
            required_inputs=(
                ATMOSPHERE_TEMPERATURE_STATE,
                "atmosphere_layer_heat_capacity",
                "atmosphere_layer_pressure",
                "atmosphere_layer_mass",
                "layer_radiative_heating",
                SURFACE_TEMPERATURE_STATE,
                "surface_heat_capacity",
            ),
            produced_outputs=(
                "convective_temperature_tendency",
                "convective_layer_count",
                "atmosphere_potential_temperature",
                "surface_sensible_heat_flux",
            ),
            conserved_quantities=("energy",),
            internal_sources_and_sinks=(),
            provided_interfaces=(CONVECTION_INTERFACE,),
            required_interfaces=(("atmospheric_state", VersionRange.parse(">=1.0.0,<2.0.0")),),
            uncertainty_outputs=(
                "instantaneous_mixing_structural",
                "dry_versus_moist_target_structural",
                "operator_splitting_order",
            ),
            correctness_tests=(
                "test_convective_adjustment_conserves_column_enthalpy_exactly",
                "test_convective_adjustment_leaves_a_stable_profile",
                "test_an_already_stable_profile_is_left_untouched",
                "test_convection_adds_no_energy_to_the_column",
                "test_radiative_convective_equilibrium_is_stable_and_closes",
                "test_sensible_heat_flux_vanishes_when_isothermal",
                "test_sensible_heat_is_transferred_not_created",
                "test_convection_changes_the_profile_it_is_given",
            ),
            known_limitations=(
                "Targets the dry adiabat, so a saturated troposphere is left more stable "
                "than the moist adiabat it should relax to.",
                "Instantaneous mixing means the result is independent of the timestep, "
                "which is convenient and physically wrong at short timesteps.",
                "The surface exchange coefficient is a free parameter with no wind or "
                "stability dependence, so the partition between the radiative and the "
                "turbulent surface fluxes is set by that choice and not derived.",
                "Vapour transport, when enabled, mixes only the layers the adjustment "
                "mixes -- and in a column whose vapour is trapped in the lowest layer "
                "that is the wrong set. Trapped vapour makes that layer an efficient "
                "radiator, it cools 45 K below the layer above it, and the resulting "
                "inversion decouples it from the neutral region higher up. So the "
                "adjustment is inactive exactly where the water is, and the transport "
                "moves nothing: the missing transport sustains the stratification that "
                "blocks it. Measured at 20 layers, the potential-temperature step across "
                "the lowest interface is 50 K with the vapour radiatively active and 0 K "
                "with it inert.",
            ),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(n_layers^2) worst case, O(n_layers) typical",
                memory_complexity="O(n_layers)",
                scaling_variables=("n_layers",),
                work_units_per_step=0.5,
                benchmark_configuration="8 layers, global scalar column",
                peak_memory_bytes=4096,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
                lower_fidelity_alternatives=("no_convection_pure_radiative_equilibrium",),
            ),
        )
        super().__init__(contract, status=status)
        self.atmosphere_module_id = atmosphere_module_id
        self.surface_module_id = surface_module_id
        self.sensible_heat_coefficient = sensible_heat_coefficient
        self.transport_vapour = transport_vapour
        self._mixing_groups: tuple[tuple[int, int], ...] = ()
        self._vapour_moved = 0.0
        self.reference_pressure = reference_pressure
        self._tendency = np.zeros(0)
        self._adjusted_layers = 0
        self._potential_temperature = np.zeros(0)
        self._sensible_flux = 0.0

    def validate_configuration(self) -> None:
        super().validate_configuration()
        if self.sensible_heat_coefficient < 0.0:
            raise ConfigurationError(
                f"sensible heat coefficient {self.sensible_heat_coefficient} cannot be negative"
            )

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        self.mark_initialized()

    def _vapour_transport(
        self, groups: tuple[tuple[int, int], ...], layer_mass: np.ndarray
    ) -> np.ndarray | None:
        """Mix the vapour over the layers convection just mixed the heat over.

        The closure homogenises potential temperature over a group because that
        is what an overturning circulation does to a conserved scalar it
        carries. Water is another such scalar, so it is mixed to a uniform
        *specific* humidity over exactly the same layers, with the group's total
        vapour mass unchanged. No new parameter enters: the group boundaries
        were decided by the temperature adjustment, which had none.

        What this does not represent is entrainment. Real convection dilutes as
        it rises and detrains air that is not saturated, so the observed
        mid-troposphere is far drier than a fully mixed column. Mixing
        completely is therefore an upper bound on the moistening, and the water
        module's saturation adjustment will leave the layers above at saturation
        rather than at the observed thirty to fifty percent.
        """
        vapour = self._inputs.get(VAPOUR_STATE)
        if vapour is None:
            return None
        w = np.array(vapour.require_units("kg m^-2").values, dtype=float)
        mixed = np.array(w)
        for low, high in groups:
            span = slice(low, high + 1)
            total_water = float(np.sum(w[span]))
            total_mass = float(np.sum(layer_mass[span]))
            if total_mass <= 0.0:
                continue
            # Uniform specific humidity, and the group's water is unchanged by
            # construction rather than to a tolerance -- the same argument the
            # enthalpy-conserving mixed temperature rests on.
            mixed[span] = total_water * layer_mass[span] / total_mass
        return mixed - w

    def compute_tendencies(self) -> list[TendencyContribution]:
        """Adjust the profile the radiation module has already heated.

        Convection is applied to the *provisional* post-radiative profile rather
        than to the state at the start of the step, so the two tendencies compose
        as sequential operator splitting: the committed temperature is exactly the
        adjusted profile, not the adjusted profile plus a stale heating term.
        """
        ctx = self.context
        temperature = self.input_field(ATMOSPHERE_TEMPERATURE_STATE).require_units("K").values
        heat_capacity = (
            self.input_field("atmosphere_layer_heat_capacity").require_units("J m^-2 K^-1").values
        )
        pressure = self.input_field("atmosphere_layer_pressure").require_units("Pa").values
        heating = self.input_field("layer_radiative_heating").require_units("W m^-2").values

        if ctx.dt <= 0.0:
            raise ConfigurationError("convective adjustment needs a positive timestep")

        # In a moist configuration the water module has already released latent
        # heat this step. Convection must adjust the profile *after* that, or it
        # would stabilise a column that condensation is about to destabilise. The
        # field is optional because a dry configuration has no water module, and
        # its absence is a configuration fact rather than a missing input.
        latent = self._inputs.get("latent_heating_rate")
        if latent is not None:
            heating = heating + latent.require_units("W m^-2").values

        # Sensible heat out of the surface and into the lowest layer. Without it
        # the surface and the air exchange energy only by radiation, the
        # surface-to-air temperature jump cannot drive anything, and convective
        # adjustment finds nothing to do -- which is what a column lacking this
        # flux actually behaves like, not an argument for leaving it out.
        surface_temperature = (
            self.input_field(SURFACE_TEMPERATURE_STATE).require_units("K").scalar()
        )
        surface_heat_capacity = (
            self.input_field("surface_heat_capacity").require_units("J m^-2 K^-1").scalar()
        )
        self._sensible_flux = self.sensible_heat_coefficient * (
            surface_temperature - float(temperature[0])
        )
        heating = np.array(heating, dtype=float)
        heating[0] += self._sensible_flux

        provisional = temperature + (heating / heat_capacity) * ctx.dt
        adjusted, groups = convective_adjustment_with_groups(
            provisional,
            heat_capacity,
            pressure,
            reference_pressure=self.reference_pressure,
        )
        self._mixing_groups = groups

        self._tendency = (adjusted - provisional) / ctx.dt
        self._adjusted_layers = int(np.count_nonzero(~np.isclose(adjusted, provisional)))
        self._potential_temperature = potential_temperature(
            adjusted, pressure, self.reference_pressure
        )

        layer_mass = np.asarray(
            self.input_field("atmosphere_layer_mass").require_units("kg m^-2").values, dtype=float
        )
        vapour_change = (
            self._vapour_transport(groups, layer_mass) if self.transport_vapour else None
        )
        if vapour_change is not None:
            # A redistribution, so it sums to zero by construction. Anything
            # else would be convection creating or destroying water.
            self._vapour_moved = float(np.sum(np.abs(vapour_change)) / 2.0)

        # The sensible flux is submitted separately from the mixing so the two
        # appear as distinct lines in the layer budget: one moves energy across a
        # module boundary, the other only within the column.
        sensible_tendency = np.zeros_like(self._tendency)
        sensible_tendency[0] = self._sensible_flux / heat_capacity[0]

        return [
            TendencyContribution(
                state_id=ATMOSPHERE_TEMPERATURE_STATE,
                provider_module_id=self.module_id,
                label="surface_sensible_heat",
                value=sensible_tendency,
                units="K s^-1",
                equation_id=EQUATION_IDS["sensible_heat"],
            ),
            TendencyContribution(
                state_id=ATMOSPHERE_TEMPERATURE_STATE,
                provider_module_id=self.module_id,
                label="dry_convective_adjustment",
                value=self._tendency,
                units="K s^-1",
                equation_id=EQUATION_IDS["convective_adjustment"],
            ),
            TendencyContribution(
                state_id=SURFACE_TEMPERATURE_STATE,
                provider_module_id=self.module_id,
                label="surface_sensible_heat",
                value=-self._sensible_flux / surface_heat_capacity,
                units="K s^-1",
                equation_id=EQUATION_IDS["sensible_heat"],
            ),
            *(
                (
                    TendencyContribution(
                        state_id=VAPOUR_STATE,
                        provider_module_id=self.module_id,
                        label="convective_vapour_transport",
                        value=vapour_change / ctx.dt,
                        units="kg m^-2",
                        equation_id=EQUATION_IDS["convective_adjustment"],
                    ),
                )
                if vapour_change is not None
                else ()
            ),
        ]

    def emit_outputs(self) -> dict[str, Field]:
        ctx = self.context
        provenance = Provenance(
            producer=self.module_id, derivation=EQUATION_IDS["convective_adjustment"]
        )
        return {
            "convective_temperature_tendency": Field(
                field_id="convective_temperature_tendency",
                physical_name="temperature tendency from convective mixing",
                symbol="dT/dt|conv",
                values=np.array(self._tendency),
                units="K s^-1",
                quantity_kind=QuantityKind.TENDENCY,
                time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                time_interval=ctx.dt,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.LAYER_AVERAGE,
                dimensions=("layer",),
                coordinate_system="pressure_layers",
                sign_convention=SignConvention.POSITIVE_INTO_RECEIVER,
                provenance=provenance,
            ),
            "atmosphere_potential_temperature": Field(
                field_id="atmosphere_potential_temperature",
                physical_name="dry potential temperature after adjustment",
                symbol="theta",
                values=np.array(self._potential_temperature),
                units="K",
                quantity_kind=QuantityKind.DIAGNOSTIC,
                time_semantics=TimeSemantics.END_OF_STEP,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.LAYER_AVERAGE,
                dimensions=("layer",),
                coordinate_system="pressure_layers",
                provenance=provenance,
            ),
            "surface_sensible_heat_flux": Field(
                field_id="surface_sensible_heat_flux",
                physical_name="sensible heat flux from the surface into the lowest layer",
                symbol="H",
                values=np.array(self._sensible_flux),
                units="W m^-2",
                quantity_kind=QuantityKind.FLUX,
                time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                time_interval=ctx.dt,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                sign_convention=SignConvention.POSITIVE_UPWARD,
                provenance=Provenance(
                    producer=self.module_id, derivation=EQUATION_IDS["sensible_heat"]
                ),
            ),
            "convective_layer_count": Field(
                field_id="convective_layer_count",
                physical_name="layers the adjustment altered this step",
                values=np.array(float(self._adjusted_layers)),
                units="1",
                quantity_kind=QuantityKind.DIAGNOSTIC,
                time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                time_interval=ctx.dt,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.DOMAIN_INTEGRAL,
                provenance=provenance,
            ),
        }

    def report_uncertainty(self) -> UncertaintyReport:
        return UncertaintyReport(
            module_id=self.module_id,
            quantified={
                "layers_adjusted": float(self._adjusted_layers),
                "sensible_heat_flux_W_m2": self._sensible_flux,
            },
            method="none propagated; the closure's error is structural, not parametric",
            distributional_assumptions="unresolved",
            correlations="unresolved",
            unresolved=(
                "structural_uncertainty_of_instantaneous_complete_mixing",
                "structural_uncertainty_of_a_dry_rather_than_moist_target",
                "sensible_heat_coefficient_parameter_uncertainty",
            ),
        )

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "target": "dry_adiabat",
            "sensible_heat_coefficient_W_m2_K": self.sensible_heat_coefficient,
            "reference_pressure_Pa": self.reference_pressure,
            "mixing": "instantaneous_complete_within_unstable_groups",
        }
