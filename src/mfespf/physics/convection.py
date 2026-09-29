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

The surface-based option qualifies that last sentence, and had to declare itself
to keep it true. Comparing each layer with the one above it makes an inversion at
the bottom of the column invisible, however buoyant the ground underneath is, and
that is the state a column with radiatively active vapour trapped in its lowest
layer settles into. Admitting the surface as element zero of the profile removes
the blind spot for the price of an exchange that *does* cross a module boundary,
so the flux is published and both the surface and the atmosphere declare the
transfer. The framework found the omission before the physics did: the first
version balanced its own arithmetic perfectly and failed reconciliation at
3e-4.
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
from mfespf.contracts.uncertainty import DeclaredRange, RangeBasis
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
)
from mfespf.core.versioning import (
    InterfaceContract,
    RequiredInterface,
    Version,
    VersionRange,
)
from mfespf.physics.atmosphere import ATMOSPHERE_TEMPERATURE_STATE
from mfespf.physics.equations import EQUATION_IDS
from mfespf.physics.plume import (
    BULK_PLUME_PARAMETER,
    bulk_plume_declaration,
    entraining_lapse_rate,
)
from mfespf.physics.surface import SURFACE_TEMPERATURE_STATE
from mfespf.physics.thermodynamics import (
    convective_adjustment_with_groups,
    hydrostatic_heights,
    lapse_rate_offsets,
    potential_temperature,
    surface_coupled_adjustment,
)
from mfespf.physics.water import VAPOUR_STATE
from mfespf.registries.states import TendencyContribution

__all__ = ["CONVECTION_INTERFACE", "ConvectiveAdjustmentModule"]

V1 = Version(1, 0, 0)

#: Profiles the adjustment can relax a column to.
_TARGETS = frozenset({"dry_adiabat", "entraining_plume"})

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


#: Declared range of the bulk sensible-heat coefficient ``k_H = rho c_p C_H |U|``,
#: in W m^-2 K^-1. The bounds are the arithmetic of typical bulk values --
#: ``C_H`` of 1.0e-3 to 1.5e-3 and a global-mean near-surface wind of 5 to
#: 8 m s^-1 at ``rho c_p`` of 1205 J m^-3 K^-1 -- which is a judgement about
#: which values are typical, not a measured global mean, and is labelled so.
SENSIBLE_HEAT_COEFFICIENT_RANGE = (6.0, 14.5)


def sensible_heat_declaration(nominal: float) -> DeclaredRange:
    low, high = SENSIBLE_HEAT_COEFFICIENT_RANGE
    return DeclaredRange(
        parameter="sensible_heat_coefficient_W_m2_K",
        keyword="sensible_heat_coefficient",
        nominal=nominal,
        low=low,
        high=high,
        basis=RangeBasis.ASSERTED,
        source=(
            "rho c_p C_H |U| over C_H of 1.0e-3 to 1.5e-3 and |U| of 5 to 8 m/s; "
            "typical bulk values, not a measured global mean"
        ),
    )


#: The declared range of the convective relaxation timescale, in seconds.
#: Published: Frierson (2007) relaxes over 2 h and tested 1 h to 16 h, finding
#: an aquaplanet remarkably insensitive across that range; Betts (1986).
RELAXATION_TIMESCALE_RANGE_S = (3600.0, 16 * 3600.0)


def relaxation_declaration(nominal: float) -> DeclaredRange:
    low, high = RELAXATION_TIMESCALE_RANGE_S
    return DeclaredRange(
        parameter="relaxation_timescale_s",
        keyword="convective_relaxation_timescale",
        nominal=nominal,
        low=low,
        high=high,
        basis=RangeBasis.PUBLISHED,
        source=(
            "Frierson (2007) uses 2 h and tested 1 h to 16 h; Betts (1986). Selected "
            "under docs/protocols/convective-relaxation-timescale.md"
        ),
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
        surface_based: bool = False,
        surface_pressure: float | None = None,
        target: str = "dry_adiabat",
        bulk_plume_parameter: float = BULK_PLUME_PARAMETER,
        target_responds_to_temperature: bool = True,
        reference_pressure: float | None = None,
        relaxation_timescale: float = 0.0,
        status: ModuleStatus = ModuleStatus.ACTIVE,
    ) -> None:
        relaxed = relaxation_timescale > 0.0
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
                *((EQUATION_IDS["relaxed_convective_adjustment"],) if relaxed else ()),
                EQUATION_IDS["sensible_heat"],
                *((EQUATION_IDS["entraining_lapse_rate"],) if target == "entraining_plume" else ()),
            ),
            assumptions=(
                "the surface exchange coefficient absorbs wind speed and stability",
                (
                    "convective mixing relaxes toward neutrality over a declared timescale"
                    if relaxed
                    else "convective mixing is instantaneous on the coupling timestep"
                ),
                "mixing is complete within an unstable group and absent outside it",
                "dry adiabatic neutrality is the target, with no moist correction"
                if target == "dry_adiabat"
                else "the neutral profile is the one an entraining plume leaves, which "
                "is more stable than a moist adiabat and less stable than a dry one",
                "no momentum, moisture or tracer transport accompanies the heat mixing",
                *(
                    (
                        "the surface is a column element the mixing can reach, so a "
                        "buoyant surface is coupled to the air above it at whatever rate "
                        "instantaneous mixing implies rather than at a fitted one",
                    )
                    if surface_based
                    else ()
                ),
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
                # Optional in general -- a dry column has no water module -- but
                # not optional once vapour transport is switched on, and a module
                # that reads a field it declares nowhere is a contract defect
                # whatever its ``.get`` says.
                *(
                    (VAPOUR_STATE, "vapour_phase_change_rate", "latent_heating_rate")
                    if transport_vapour
                    else ()
                ),
            ),
            produced_outputs=(
                "convective_temperature_tendency",
                "convective_layer_count",
                "convectively_mixed_layer_temperature",
                "atmosphere_potential_temperature",
                "surface_sensible_heat_flux",
                *(("surface_convective_heat_flux",) if surface_based else ()),
            ),
            conserved_quantities=("energy",),
            internal_sources_and_sinks=(),
            provided_interfaces=(CONVECTION_INTERFACE,),
            required_interfaces=(
                RequiredInterface(
                    interface_id="atmospheric_state",
                    version_range=VersionRange.parse(">=1.0.0,<2.0.0"),
                    expects={
                        "layer_ordering": "index_0_is_the_lowest_layer",
                        "vertical_coordinate": "pressure",
                        "temperature_units": "K",
                    },
                ),
            ),
            equation_variables=(
                (
                    EQUATION_IDS["potential_temperature"],
                    (
                        ("T", ATMOSPHERE_TEMPERATURE_STATE, "K"),
                        ("p", "atmosphere_layer_pressure", "Pa"),
                    ),
                ),
                # Bound to the masked profile, not the whole column. The
                # adjustment mixes only inside its unstable groups, so a
                # stratospheric temperature judged against these ranges flags a
                # relation whose output never reaches that layer. The potential
                # temperature above is different and stays bound to the full
                # state: it is diagnosed for every layer.
                (
                    EQUATION_IDS["convective_adjustment"],
                    (("T", "convectively_mixed_layer_temperature", "K"),),
                ),
                *(
                    (
                        (
                            EQUATION_IDS["relaxed_convective_adjustment"],
                            (("T", "convectively_mixed_layer_temperature", "K"),),
                        ),
                    )
                    if relaxed
                    else ()
                ),
                (
                    EQUATION_IDS["entraining_lapse_rate"],
                    (
                        ("T", "convectively_mixed_layer_temperature", "K"),
                        ("a", "bulk_plume_parameter", "1"),
                    ),
                ),
            ),
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
                "test_the_transport_mixes_the_profile_condensation_leaves_not_the_one_it_found",
                "test_the_adjustment_publishes_the_layers_it_mixed_not_the_whole_profile",
            ),
            known_limitations=(
                "Targets the dry adiabat by default, so a saturated troposphere is left "
                "more stable than the moist adiabat it should relax to. The "
                "entraining-plume target removes that, and a column carrying the "
                "bulk-plume humidity should use it: measured at 20 layers, the dry "
                "target leaves 9.68 K/km against the 6.33 the same plume parameter "
                "predicts, so the two closures contradict each other.",
                "The entraining target makes the stratification depend on the "
                "saturation humidity, so any experiment holding the vapour fixed must "
                "set target_responds_to_temperature=False as well. A control that "
                "freezes only the radiative half of the vapour is not a control.",
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
                "with it inert. Removed by the surface-based option below, which is "
                "why that option exists; the finding stands for the layer-only form.",
                "Surface-based mixing is complete rather than entraining. Real plumes "
                "dilute as they rise and detrain air that is not saturated; mixing "
                "completely instead drives the top of the group past saturation, "
                "condensation clips it there, and nothing moistens the column above it "
                "at all. The vapour feedback that results is only about +0.5 W m^-2 "
                "K^-1 against +1.251 for a fixed-relative-humidity control, and it is "
                "not safely converged in the layer grid -- amplification runs 1.10x, "
                "1.22x, 1.88x, 1.28x and 1.30x at 10, 20, 30, 40 and 60 layers, and the "
                "excursion at 30 is unexplained.",
                "The surface joins the mixing at a temperature this module forms from "
                "the fluxes it can see -- its own sensible flux and the surface net "
                "radiative flux. Any other tendency the surface carries, deep-ocean "
                "exchange in particular, is left to the splitting, so in a coupled-ocean "
                "configuration the surface convects from a slightly stale temperature.",
            ),
            fidelity_justification=(
                "Pure radiative equilibrium is convectively unstable in the lower "
                "troposphere, so it produces a lapse rate steeper than any observed one "
                "and a surface far too warm. The comparison is the point of Stage 3 and "
                "is kept as the null control."
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
        self.surface_based = surface_based
        self.surface_pressure = surface_pressure
        self.target = target
        self.bulk_plume_parameter = bulk_plume_parameter
        self.target_responds_to_temperature = target_responds_to_temperature
        # Zero is the instantaneous closure, exactly; anything positive relaxes
        # toward the same adjusted profile over that many seconds.
        self.relaxation_timescale = relaxation_timescale
        self._held_rate: np.ndarray | None = None
        self._mixing_groups: tuple[tuple[int, int], ...] = ()
        self._adjusted = np.zeros(0)
        self._vapour_dt = 0.0
        self._vapour_moved = 0.0
        self._surface_adjustment = 0.0
        self._surface_convective_flux = 0.0
        self._surface_mixed = False
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
        if self.relaxation_timescale < 0.0:
            raise ConfigurationError(
                f"convective relaxation timescale {self.relaxation_timescale} cannot be "
                "negative; zero is the instantaneous adjustment"
            )
        if self.target not in _TARGETS:
            raise ConfigurationError(
                f"unknown convective target {self.target!r}; expected one of "
                f"{', '.join(sorted(_TARGETS))}"
            )
        if self.target == "entraining_plume" and not self.surface_based:
            # The plume target is a lapse rate measured upward from the ground,
            # so a scheme that cannot see the ground has nothing to measure it
            # from. Refused rather than silently anchored to the lowest layer.
            raise ConfigurationError(
                "the entraining-plume target needs surface_based convection, because "
                "its lapse rate is anchored at the surface"
            )

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        self.mark_initialized()

    def _reference_surface_pressure(self) -> float:
        """The pressure the surface element sits at, a boundary condition.

        Defaults to the standard surface pressure, which is what the atmosphere
        module builds its layers from; a configuration that moves one has to
        move both, so the column keeps a single bottom.
        """
        return (
            self.surface_pressure
            if self.surface_pressure is not None
            else value_of("standard_surface_pressure")
        )

    def _plume_offsets(self, temperature: np.ndarray, pressure: np.ndarray) -> np.ndarray:
        """The temperature drop the entraining plume implies, surface upward.

        Built from the column's own temperature at every step rather than once,
        because the target lapse rate depends on the saturation humidity and so
        on the profile it is adjusting. That makes the adjustment a fixed point
        rather than a relaxation to a fixed profile, which is the honest form:
        a warmer column has a smaller entraining lapse rate, and pretending
        otherwise would freeze the lapse-rate feedback this is meant to supply.
        """
        surface_pressure = self._reference_surface_pressure()
        height = hydrostatic_heights(temperature, pressure, surface_pressure)
        if self.target_responds_to_temperature:
            rate = entraining_lapse_rate(temperature, pressure, self.bulk_plume_parameter)
        else:
            # The null control. Water vapour acts on this column twice -- once
            # radiatively through the absorption, and once dynamically because
            # the entraining lapse rate depends on the saturation humidity -- so
            # a control that holds the vapour has to hold both, or it is not
            # holding the vapour. Freezing only the radiative half left the
            # frozen run with the entire lapse-rate response and a feedback
            # parameter of 12.4 W m^-2 K^-1 against a Planck response of 3.2.
            #
            # The *rate* is frozen and not the offsets: the column still expands
            # as it warms, which is hypsometry and not the vapour's doing.
            if self._held_rate is None:
                self._held_rate = entraining_lapse_rate(
                    temperature, pressure, self.bulk_plume_parameter
                )
            rate = self._held_rate
        # The surface joins as the bottom of the column, at zero height and
        # carrying the lowest layer's lapse rate.
        return lapse_rate_offsets(
            np.concatenate(([rate[0]], rate)), np.concatenate(([0.0], height))
        )

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

        The profile mixed here is the *provisional* one -- the committed vapour
        plus what the water module has already published as its own phase-change
        rate for this step -- for the same reason the temperature mixed here is
        the post-radiative, post-latent one. Mixing the start-of-step vapour
        instead made the two tendencies additive rather than sequential, so both
        modules charged the same layer for water only one of them could have:
        each affordable alone, jointly drawing 3.79 kg m^-2 from a layer holding
        3.17. Composed this way the committed vapour is exactly ``mixed``, which
        is non-negative by construction because the water module's provisional
        profile is, and the group total is untouched.
        """
        vapour = self._inputs.get(VAPOUR_STATE)
        if vapour is None:
            return None
        w = np.array(vapour.require_units("kg m^-2").values, dtype=float)
        self._vapour_dt = self.context.dt
        phase = self._inputs.get("vapour_phase_change_rate")
        if phase is not None:
            # Over the interval the *water* module published the rate for, and
            # the contribution below is divided by the same interval, because
            # that is the one the owner of the vapour will multiply back when it
            # commits. Using this module's step instead would make the
            # composition exact only where the two happen to share a schedule,
            # and silently inexact where they do not.
            self._vapour_dt = float(phase.time_interval or self.context.dt)
            w = (
                w
                + np.asarray(phase.require_units("kg m^-2 s^-1").values, dtype=float)
                * self._vapour_dt
            )
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

    def relaxation_fraction(self, dt: float) -> float:
        """The share of the instantaneous adjustment applied over ``dt``.

        ``1 - exp(-dt / tau)``: the exact solution of relaxing toward a fixed
        target at rate ``1 / tau``, so a steady target is reached the same way
        whatever the step. Exactly one at ``tau = 0``.
        """
        if self.relaxation_timescale <= 0.0:
            return 1.0
        return float(-np.expm1(-dt / self.relaxation_timescale))

    def _mixed_mask(self) -> np.ndarray:
        """True on the layers this run's adjustment has mixed.

        The closure evaluates its stability criterion everywhere and mixes only
        within the unstable groups, so these are the layers its cited relations
        were applied to -- and the only ones whose temperatures a declared
        validity range has any business being judged against.

        This step's groups, which is the honest statement and was briefly not
        available. While the domain check judged a run by its final state, a
        per-step mask made the identical configuration resolvable at 2000 steps
        and unresolvable at 400 -- this closure mixes on 79 percent of steps and
        on none of the first four hundred, the column being dry-stable while it
        spins up -- so the mask had to be the union over the run to have any
        footprint left at the end. That union was the right answer to "which
        layers does this relation govern" and the wrong one to "was it in range
        when it ran": a layer mixed at 250 K during spin-up was judged at
        whatever it held when the run stopped. Now that the check accumulates
        each relation's extent step by step, the union is no longer needed and
        the mask means what it says.
        """
        mask = np.zeros(len(self._adjusted), dtype=bool)
        for low, high in self._mixing_groups:
            mask[low : high + 1] = True
        return mask

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
        if self.surface_based:
            # The surface enters the profile the adjustment acts on, so its own
            # provisional temperature has to be formed the way the layers' is:
            # the state at the start of the step plus every flux this module can
            # see. The radiative flux is folded in here and subtracted again in
            # the tendency, exactly as ``layer_radiative_heating`` is for the
            # layers, so the surface module still commits its own radiation once
            # and only once.
            forcing = -self._sensible_flux
            net = self._inputs.get("surface_net_radiative_flux")
            if net is not None:
                forcing += net.require_units("W m^-2").scalar()
            provisional_surface = surface_temperature + (forcing / surface_heat_capacity) * ctx.dt
            (
                adjusted_surface,
                adjusted,
                groups,
                self._surface_mixed,
            ) = surface_coupled_adjustment(
                provisional,
                heat_capacity,
                pressure,
                surface_temperature=provisional_surface,
                surface_heat_capacity=surface_heat_capacity,
                surface_pressure=self._reference_surface_pressure(),
                target_offsets=(
                    self._plume_offsets(provisional, pressure)
                    if self.target == "entraining_plume"
                    else None
                ),
                reference_pressure=self.reference_pressure,
            )
            self._surface_adjustment = (adjusted_surface - provisional_surface) / ctx.dt
            # Positive upward: what the surface loses, the air gains, and the
            # mixing conserves their sum by construction.
            self._surface_convective_flux = -self._surface_adjustment * surface_heat_capacity
        else:
            adjusted, groups = convective_adjustment_with_groups(
                provisional,
                heat_capacity,
                pressure,
                reference_pressure=self.reference_pressure,
            )
            self._surface_adjustment = 0.0
            self._surface_convective_flux = 0.0
            self._surface_mixed = False
        self._mixing_groups = groups

        # The relaxed closure applies the exactly integrated fraction of the
        # same adjustment. One scalar across the whole column, so every group
        # still sums to zero and the enthalpy closes exactly, not on average.
        fraction = self.relaxation_fraction(ctx.dt)
        adjusted = provisional + fraction * (np.asarray(adjusted, dtype=float) - provisional)
        self._surface_adjustment *= fraction
        self._surface_convective_flux *= fraction

        self._tendency = (adjusted - provisional) / ctx.dt
        self._adjusted = np.asarray(adjusted, dtype=float)
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
            # A convex combination of the provisional and the mixed vapour, both
            # non-negative, so this cannot over-draw a layer either.
            vapour_change = fraction * vapour_change
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
                        state_id=SURFACE_TEMPERATURE_STATE,
                        provider_module_id=self.module_id,
                        label="surface_convective_adjustment",
                        value=self._surface_adjustment,
                        units="K s^-1",
                        equation_id=EQUATION_IDS["convective_adjustment"],
                    ),
                )
                if self.surface_based
                else ()
            ),
            *(
                (
                    TendencyContribution(
                        state_id=VAPOUR_STATE,
                        provider_module_id=self.module_id,
                        label="convective_vapour_transport",
                        # Divided by the interval the provisional profile was
                        # formed over, so the owner's commit reproduces the
                        # mixed profile exactly rather than to a schedule.
                        value=vapour_change / self._vapour_dt,
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
            # The layers the closure actually ran on. It computes a target for
            # every layer and mixes only inside the unstable groups, so judging
            # its declared temperature range against the whole profile lets a
            # stratospheric layer flag a relation whose output is never used up
            # there. NaN says "not evaluated here" rather than "missing".
            "convectively_mixed_layer_temperature": Field(
                field_id="convectively_mixed_layer_temperature",
                physical_name="layer temperature where the adjustment has mixed",
                symbol="T_conv",
                values=np.where(
                    self._mixed_mask(), np.asarray(self._adjusted, dtype=float), np.nan
                ),
                units="K",
                quantity_kind=QuantityKind.DIAGNOSTIC,
                time_semantics=TimeSemantics.END_OF_STEP,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.LAYER_AVERAGE,
                dimensions=("layer",),
                coordinate_system="pressure_layers",
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                # The mask is the content of this field, so it declares the
                # policy DAT-007 provides rather than slipping NaNs past a
                # default that exists to catch them.
                missing_value_policy=MissingValuePolicy.MASKED,
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
            **(
                {
                    # Mixing the surface into the column moves energy across a
                    # module boundary, which the sensible flux also does and
                    # which the column's energy ledger has to see. Published as
                    # a flux so the atmosphere can declare the import, rather
                    # than left implicit in a temperature tendency: an exchange
                    # no ledger records is indistinguishable from an invention.
                    "surface_convective_heat_flux": Field(
                        field_id="surface_convective_heat_flux",
                        physical_name=(
                            "heat carried from the surface into the air by convective mixing"
                        ),
                        symbol="H_conv",
                        values=np.array(self._surface_convective_flux),
                        units="W m^-2",
                        quantity_kind=QuantityKind.FLUX,
                        time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                        time_interval=ctx.dt,
                        timestamp=ctx.time,
                        spatial_support=SpatialSupport.GLOBAL_SCALAR,
                        sign_convention=SignConvention.POSITIVE_UPWARD,
                        provenance=provenance,
                    )
                }
                if self.surface_based
                else {}
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
                "surface_took_part_in_the_mixing": float(self._surface_mixed),
                "surface_convective_flux_W_m2": self._surface_convective_flux,
            },
            method=(
                "declared ranges, propagated only by an ensemble built from them; the "
                "closure's own error is structural and is not among them"
            ),
            distributional_assumptions=(
                "a range, not a distribution: no probability is assigned within the "
                "declared ranges, so a sweep across them states no coverage"
            ),
            correlations=(
                "the sensible heat coefficient fixes the evaporation coefficient exactly "
                "(k_E = k_H / c_p), so the two surface fluxes move together"
                + (
                    "; the bulk-plume parameter is shared with the humidity closure"
                    if self.target == "entraining_plume"
                    else ""
                )
            ),
            unresolved=(
                "structural_uncertainty_of_instantaneous_complete_mixing",
                "structural_uncertainty_of_a_dry_rather_than_moist_target",
            ),
            declared_ranges=(
                sensible_heat_declaration(self.sensible_heat_coefficient),
                *(
                    (relaxation_declaration(self.relaxation_timescale),)
                    if self.relaxation_timescale > 0.0
                    else ()
                ),
                *(
                    (bulk_plume_declaration(self.bulk_plume_parameter),)
                    if self.target == "entraining_plume"
                    else ()
                ),
            ),
        )

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "target": self.target,
            "bulk_plume_parameter": (
                self.bulk_plume_parameter if self.target == "entraining_plume" else None
            ),
            "target_responds_to_temperature": self.target_responds_to_temperature,
            "surface_based": self.surface_based,
            "transport_vapour": self.transport_vapour,
            "surface_pressure_Pa": (
                self._reference_surface_pressure() if self.surface_based else None
            ),
            "sensible_heat_coefficient_W_m2_K": self.sensible_heat_coefficient,
            "reference_pressure_Pa": self.reference_pressure,
            "relaxation_timescale_s": self.relaxation_timescale,
            "mixing": (
                "relaxed_toward_complete_within_unstable_groups"
                if self.relaxation_timescale > 0.0
                else "instantaneous_complete_within_unstable_groups"
            ),
        }
