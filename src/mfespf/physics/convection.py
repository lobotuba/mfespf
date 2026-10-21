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
from mfespf.physics.advection import advection_field, composed_advection
from mfespf.physics.atmosphere import ATMOSPHERE_TEMPERATURE_STATE
from mfespf.physics.equations import EQUATION_IDS
from mfespf.physics.mass_flux import (
    CAPE_DEFINITIONS,
    CAPE_TRIGGER_J_KG,
    CLOUD_WATER_CONVERSION_PER_M,
    CONDENSATE_PATHWAYS,
    CONVECTIVE_ADJUSTMENT_TIME_S,
    DETRAINED_STATES,
    PLUME_ENTRAINMENT_PER_M,
    PRECIPITATION_EVAPORATION_EFFICIENCY,
    PRECIPITATION_EVAPORATIONS,
    TRIGGER_CROSSINGS,
    WINDOW_CONTINUATIONS,
    convect,
)
from mfespf.physics.plume import (
    BULK_PLUME_PARAMETER,
    EVAPORATING_CONDENSATE_FRACTION,
    bulk_plume_declaration,
    entraining_lapse_rate,
)
from mfespf.physics.surface import SURFACE_TEMPERATURE_STATE
from mfespf.physics.thermodynamics import (
    KAPPA,
    convective_adjustment_with_groups,
    hydrostatic_heights,
    lapse_rate_offsets,
    lifting_condensation_pressure,
    potential_temperature,
    saturation_specific_humidity,
    surface_coupled_adjustment,
)
from mfespf.physics.water import VAPOUR_STATE
from mfespf.registries.states import TendencyContribution

__all__ = ["CONVECTION_INTERFACE", "ConvectiveAdjustmentModule"]

V1 = Version(1, 0, 0)

#: Profiles the adjustment can relax a column to.
_TARGETS = frozenset({"dry_adiabat", "entraining_plume"})
#: How far the surface-rooted group mixes vapour: over every layer the heat
#: adjustment mixed (the declared upper bound), or only up to the surface
#: parcel's lifting condensation level (the subcloud lower bound;
#: rising-mixing-variants.md).
VAPOUR_MIXING_TOPS = ("group", "lifting_condensation_level")

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


#: The mass-flux operator's decomposition, published per layer as rates
#: (detrainment-thermodynamic-state.md): field -> (MassFluxStep attribute,
#: physical name, units). The evaporation is added beside them.
MASS_FLUX_DECOMPOSITION: dict[str, tuple[str, str, str]] = {
    "mass_flux_subsidence_vapour_rate": (
        "subsidence_vapour",
        "vapour moved into each layer by compensating subsidence",
        "kg m^-2 s^-1",
    ),
    "mass_flux_detrainment_vapour_rate": (
        "detrainment_vapour",
        "net vapour added to each layer by detrainment, D (q_det - q)",
        "kg m^-2 s^-1",
    ),
    "mass_flux_subsidence_heating": (
        "subsidence_energy",
        "dry static energy moved into each layer by compensating subsidence",
        "W m^-2",
    ),
    "mass_flux_detrainment_heating": (
        "detrainment_energy",
        "net dry static energy added to each layer by detrainment, D (s_det - s)",
        "W m^-2",
    ),
    "detrained_vapour_rate": (
        "detrained_vapour",
        "vapour detrained into each layer, gross",
        "kg m^-2 s^-1",
    ),
    "detrained_condensate_rate": (
        "detrained_condensate",
        "condensate detrained into each layer",
        "kg m^-2 s^-1",
    ),
    "post_detrainment_excess_rate": (
        "post_detrainment_excess",
        "supersaturation detrainment leaves in a receiving layer, at constant h",
        "kg m^-2 s^-1",
    ),
    "convective_evaporation_rate": (
        "evaporation",
        "convective condensate evaporated in each layer",
        "kg m^-2 s^-1",
    ),
    "convective_precipitation_entering_rate": (
        "precipitation_entering",
        "convective condensate falling into each layer from above",
        "kg m^-2 s^-1",
    ),
    "convective_evaporation_limit": (
        "evaporation_limit_time",
        "the most each layer could evaporate, step-mean of each substep's limit",
        "kg m^-2",
    ),
    "convective_condensate_formation_rate": (
        "formation_gross",
        "convective condensate formed in each layer, gross (C_form)",
        "kg m^-2 s^-1",
    ),
    "convective_cloud_evaporation_rate": (
        "cloud_evaporation",
        "convective cloud condensate evaporated before conversion (E_cloud)",
        "kg m^-2 s^-1",
    ),
    "convective_rain_evaporation_rate": (
        "rain_evaporation",
        "convective precipitation evaporated (E_precip)",
        "kg m^-2 s^-1",
    ),
    "convective_precipitation_formation_rate": (
        "precipitation_formed",
        "convective precipitation formed in each layer (C_precip)",
        "kg m^-2 s^-1",
    ),
    "detrained_cloud_water_rate": (
        "detrained_cloud",
        "cloud water detrained into each layer",
        "kg m^-2 s^-1",
    ),
    "detrained_mass_rate": (
        "detrained_mass",
        "air detrained into each layer by the plume",
        "kg m^-2 s^-1",
    ),
    "plume_mass_flux_rate": (
        "plume_mass",
        "the plume's mass flux at each layer's lower interface, the subsidence above it",
        "kg m^-2 s^-1",
    ),
}


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
        fractional_top: bool = False,
        receives_advection: bool = False,
        advection_sources: tuple[str, ...] = ("",),
        diagnose_mixing_groups: bool = False,
        vapour_mixing_top: str = "group",
        condensation_follows: bool = False,
        mass_flux: bool = False,
        plume_entrainment_rate: float = PLUME_ENTRAINMENT_PER_M,
        convective_adjustment_time: float = CONVECTIVE_ADJUSTMENT_TIME_S,
        cape_trigger: float = CAPE_TRIGGER_J_KG,
        evaporating_fraction: float = EVAPORATING_CONDENSATE_FRACTION,
        detrained_state: str = "plume",
        precipitation_evaporation: str = "cascade",
        condensate_pathway: str = "immediate",
        cloud_water_conversion_rate: float = CLOUD_WATER_CONVERSION_PER_M,
        track_water_provenance: bool = False,
        convective_substep: float | None = None,
        convective_window_continuation: str = "stop",
        convective_trigger_crossing: str = "explicit",
        convective_cape_definition: str = "first_lnb",
        status: ModuleStatus = ModuleStatus.ACTIVE,
    ) -> None:
        relaxed = relaxation_timescale > 0.0
        if detrained_state not in DETRAINED_STATES:
            raise ConfigurationError(
                f"detrained_state {detrained_state!r} is not one of {DETRAINED_STATES}"
            )
        if detrained_state != "plume" and not mass_flux:
            raise ConfigurationError("detrained_state applies only to mass-flux convection")
        if precipitation_evaporation not in PRECIPITATION_EVAPORATIONS:
            raise ConfigurationError(
                f"precipitation_evaporation {precipitation_evaporation!r} is not one of "
                f"{PRECIPITATION_EVAPORATIONS}"
            )
        if precipitation_evaporation != "cascade" and not mass_flux:
            raise ConfigurationError(
                "precipitation_evaporation applies only to mass-flux convection"
            )
        if convective_cape_definition not in CAPE_DEFINITIONS:
            raise ConfigurationError(
                f"convective_cape_definition {convective_cape_definition!r} is not one of "
                f"{CAPE_DEFINITIONS}"
            )
        if convective_trigger_crossing not in TRIGGER_CROSSINGS:
            raise ConfigurationError(
                f"convective_trigger_crossing {convective_trigger_crossing!r} is not one of "
                f"{TRIGGER_CROSSINGS}"
            )
        if convective_window_continuation not in WINDOW_CONTINUATIONS:
            raise ConfigurationError(
                f"convective_window_continuation {convective_window_continuation!r} is not one "
                f"of {WINDOW_CONTINUATIONS}"
            )
        if condensate_pathway not in CONDENSATE_PATHWAYS:
            raise ConfigurationError(
                f"condensate_pathway {condensate_pathway!r} is not one of {CONDENSATE_PATHWAYS}"
            )
        if condensate_pathway == "carried" and not (
            mass_flux
            and detrained_state == "plume"
            and precipitation_evaporation in ("none", "rate_law_only")
        ):
            raise ConfigurationError(
                "carried cloud water needs mass-flux convection, the plume's own detrained "
                "state, and precipitation_evaporation 'none' or 'rate_law_only'"
            )
        if track_water_provenance and condensate_pathway != "carried":
            raise ConfigurationError(
                "water provenance is carried with the cloud-water pathway "
                "(fall-path-closure-candidate.md)"
            )
        # Mass-flux convection (mass-flux-moisture-transport.md) replaces the moist
        # adjustment, so it composes only with the dry one, and only where the
        # water module condenses after it, because the plume's condensation is
        # handed to the water module to book.
        if mass_flux and not (
            condensation_follows and transport_vapour and target == "dry_adiabat" and not relaxed
        ):
            raise ConfigurationError(
                "mass-flux convection replaces the entraining-plume adjustment: it needs "
                "target='dry_adiabat', vapour transport, condensation after convection, "
                "and no relaxation timescale"
            )
        if vapour_mixing_top not in VAPOUR_MIXING_TOPS:
            raise ConfigurationError(
                f"vapour_mixing_top {vapour_mixing_top!r} is not one of {VAPOUR_MIXING_TOPS}"
            )
        subcloud = vapour_mixing_top == "lifting_condensation_level"
        if subcloud and not (transport_vapour and surface_based and not mass_flux):
            raise ConfigurationError(
                "vapour mixing limited to the lifting condensation level needs vapour "
                "transport, surface-based convection and no mass-flux convection"
            )
        if condensation_follows and not transport_vapour:
            raise ConfigurationError(
                "condensation can only follow convection's vapour transport when "
                "convection transports vapour"
            )
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
                *(
                    (EQUATION_IDS["fractional_top_convective_adjustment"],)
                    if fractional_top
                    else ()
                ),
                EQUATION_IDS["sensible_heat"],
                *((EQUATION_IDS["entraining_lapse_rate"],) if target == "entraining_plume" else ()),
                *((EQUATION_IDS["lifting_condensation_level"],) if subcloud else ()),
                *(
                    (EQUATION_IDS["mass_flux_convection"], EQUATION_IDS["clausius_clapeyron"])
                    if mass_flux
                    else ()
                ),
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
                # When the water module runs after convection it has formed
                # neither yet, so convection mixes the committed vapour and the
                # water module condenses what the mixing lifts
                # (docs/protocols/convective-condensation-order.md).
                *(
                    (VAPOUR_STATE,)
                    if condensation_follows
                    else (VAPOUR_STATE, "vapour_phase_change_rate", "latent_heating_rate")
                    if transport_vapour
                    else ()
                ),
                # A circulation's heating and moistening this step, which the
                # adjustment must see so that it adjusts the profile the
                # circulation leaves rather than the one it found.
                *(
                    tuple(
                        advection_field(s, kind)
                        for s in advection_sources
                        for kind in ("heating", "vapour")
                    )
                    if receives_advection
                    else ()
                ),
                # The plume's layer thicknesses.
                *(("atmosphere_interface_pressure",) if mass_flux else ()),
                *(("water_provenance",) if track_water_provenance else ()),
            ),
            produced_outputs=(
                "convective_temperature_tendency",
                "convective_layer_count",
                *(("surface_mixed_layer_top",) if surface_based else ()),
                # The parcel the subcloud cut lifted, NaN on steps it was not
                # lifted, so the level's validity is judged where it ran.
                *(
                    (
                        "subcloud_parcel_temperature",
                        "subcloud_parcel_specific_humidity",
                        "subcloud_parcel_pressure",
                        "subcloud_condensation_pressure",
                    )
                    if subcloud
                    else ()
                ),
                # Off by default: two layer fields a step, wanted only to replay
                # the mixing offline (shallow-return-flow.md).
                *(
                    ("convective_mixing_group", "convective_partial_top_fraction")
                    if diagnose_mixing_groups
                    else ()
                ),
                "convectively_mixed_layer_temperature",
                "atmosphere_potential_temperature",
                "surface_sensible_heat_flux",
                *(("surface_convective_heat_flux",) if surface_based else ()),
                *(("convective_vapour_transport_rate",) if condensation_follows else ()),
                *(
                    (
                        "convective_condensation_rate",
                        "cloud_base_mass_flux",
                        "dilute_cape",
                        "cape_removal_per_unit_mass_flux",
                        "convective_surface_precipitation",
                        "plume_temperature",
                        "plume_buoyancy",
                        "plume_top_pressure",
                        "mass_flux_substeps",
                        "cape_host_start",
                        "cape_time_mean",
                        "cape_pre_convection",
                        "cape_post_convection",
                        "cape_host_end",
                        "cape_window_start",
                        "cape_window_end",
                        "convective_window_evaluated",
                        "convective_window_status",
                        "convective_window_mass",
                        "convective_window_fraction",
                        "convective_trigger_unresolved_crossings",
                        "cape_window_start_positive_area",
                        "cape_window_start_maximum_work",
                        "convective_active_time_fraction",
                        "requested_cloud_base_mass_flux",
                        "plume_top_pressure_continuous",
                        *MASS_FLUX_DECOMPOSITION,
                        *(
                            ("water_provenance_change_rate", "convective_provenance_milestones")
                            if track_water_provenance
                            else ()
                        ),
                    )
                    if mass_flux
                    else ()
                ),
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
                *(
                    (
                        (
                            EQUATION_IDS["fractional_top_convective_adjustment"],
                            (("T", "convectively_mixed_layer_temperature", "K"),),
                        ),
                    )
                    if fractional_top
                    else ()
                ),
                (
                    EQUATION_IDS["entraining_lapse_rate"],
                    (
                        ("T", "convectively_mixed_layer_temperature", "K"),
                        ("a", "bulk_plume_parameter", "1"),
                    ),
                ),
                # Bound to the parcel actually lifted, not to every layer's
                # pressure: the level is evaluated for one parcel a step.
                *(
                    (
                        (
                            EQUATION_IDS["lifting_condensation_level"],
                            (
                                ("T_0", "subcloud_parcel_temperature", "K"),
                                ("q_0", "subcloud_parcel_specific_humidity", "kg kg^-1"),
                                ("p_0", "subcloud_parcel_pressure", "Pa"),
                                ("p_L", "subcloud_condensation_pressure", "Pa"),
                            ),
                        ),
                    )
                    if subcloud
                    else ()
                ),
                # The plume saturates at its own temperature, where it exists.
                *(
                    (
                        (EQUATION_IDS["clausius_clapeyron"], (("T", "plume_temperature", "K"),)),
                        (EQUATION_IDS["mass_flux_convection"], (("T", "plume_temperature", "K"),)),
                    )
                    if mass_flux
                    else ()
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
        self.diagnose_mixing_groups = bool(diagnose_mixing_groups)
        self.vapour_mixing_top = vapour_mixing_top
        self.surface_pressure = surface_pressure
        self.target = target
        self.bulk_plume_parameter = bulk_plume_parameter
        self.target_responds_to_temperature = target_responds_to_temperature
        # Zero is the instantaneous closure, exactly; anything positive relaxes
        # toward the same adjusted profile over that many seconds.
        self.relaxation_timescale = relaxation_timescale
        # Whether a mixed group may stop part-way through a layer. Off, the
        # top of the mixing is always a layer edge and the column's answer
        # follows where the grid puts those edges.
        self.fractional_top = fractional_top
        self.receives_advection = receives_advection
        self.advection_sources = tuple(advection_sources)
        self.condensation_follows = condensation_follows
        self.mass_flux = mass_flux
        self.plume_entrainment_rate = plume_entrainment_rate
        self.convective_adjustment_time = convective_adjustment_time
        self.cape_trigger = cape_trigger
        self.evaporating_fraction = evaporating_fraction
        self.detrained_state = detrained_state
        self.precipitation_evaporation = precipitation_evaporation
        self.condensate_pathway = condensate_pathway
        self.cloud_water_conversion_rate = cloud_water_conversion_rate
        self.track_water_provenance = track_water_provenance
        self.convective_substep = convective_substep
        self.convective_window_continuation = convective_window_continuation
        self.convective_trigger_crossing = convective_trigger_crossing
        self.convective_cape_definition = convective_cape_definition
        self._unresolved_crossings = 0
        self._window_arrays: dict[str, np.ndarray] = {}
        self._cape_after = 0.0
        self._cape_time_mean = 0.0
        self._cape_host_start = 0.0
        self._cape_pre = 0.0
        self._cape_post = 0.0
        self._cape_window_start = np.zeros(1)
        self._active_fraction = 0.0
        self._requested_flux = 0.0
        self._top_continuous = float("nan")
        self._provenance_change: np.ndarray | None = None
        self._provenance_milestones: np.ndarray | None = None
        self._mass_flux_parts: dict[str, np.ndarray] = {}
        self._mass_flux_tendency: np.ndarray | None = None
        self._convective_condensation: np.ndarray | None = None
        self._plume: Any = None
        self._cloud_base_flux = 0.0
        self._cape_removal = 0.0
        self._surface_precipitation = 0.0
        self._mass_flux_substeps = 0
        self._vapour_base: np.ndarray | None = None
        self._vapour_change: np.ndarray | None = None
        self._partial_tops: dict[int, float] = {}
        self._held_rate: np.ndarray | None = None
        self._mixing_groups: tuple[tuple[int, int], ...] = ()
        # The groups the vapour is mixed over: the heat groups themselves, or
        # the surface group cut at its condensation level.
        self._vapour_groups: tuple[tuple[int, int], ...] = ()
        self._vapour_partials: dict[int, float] = {}
        self._condensation_level_layer = -1
        self._subcloud_parcel = (float("nan"),) * 4
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
        phase = None if self.condensation_follows else self._inputs.get("vapour_phase_change_rate")
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
        advection, first = composed_advection(self._inputs, self.advection_sources, "vapour")
        if advection is not None:
            w = w + advection * float(first.time_interval or self._vapour_dt)
        self._vapour_base = np.array(w)
        mixed = np.array(w)
        for low, high in groups:
            span = slice(low, high + 1)
            # The part of the layer above that the group rose into, by mass. It
            # carries its share of that layer's vapour into the mixing, and the
            # rest of the layer keeps the rest -- the same fraction the heat used.
            top = high + 1
            f = self._vapour_partials.get(top, 0.0) if top < w.size else 0.0
            total_water = float(np.sum(w[span])) + f * float(w[top] if f else 0.0)
            total_mass = float(np.sum(layer_mass[span])) + f * float(layer_mass[top] if f else 0.0)
            if total_mass <= 0.0:
                continue
            # Uniform specific humidity, and the group's water is unchanged by
            # construction rather than to a tolerance -- the same argument the
            # enthalpy-conserving mixed temperature rests on.
            # In this order so the whole-layer case is bit-identical to what it
            # always was: a limit-cycling column amplifies round-off.
            mixed[span] = total_water * layer_mass[span] / total_mass
            if f:
                mixed[top] = (1.0 - f) * w[top] + f * total_water * layer_mass[top] / total_mass
        return mixed - w

    def _cut_at_condensation_level(
        self, adjusted: np.ndarray, pressure: np.ndarray, layer_mass: np.ndarray
    ) -> None:
        """Limit the surface group's vapour mixing to the subcloud layer.

        The subcloud layer is found as a fixed point (rising-mixing-variants.md,
        erratum 1). Starting from the top of the heat adjustment's surface group,
        the parcel is layer 0's adjusted temperature and pressure carrying the
        mass-weighted mean of the pre-mixing specific humidity over the layers
        kept. If it is still unsaturated at the top layer kept, those layers are
        the subcloud layer; otherwise the top moves down to the layer holding the
        parcel's lifting condensation level (EQ-LIFTING-CONDENSATION-LEVEL) and
        the parcel is formed again. The top only descends, so this ends.

        The mean, not layer 0 alone, because layer 0 receives each step's
        evaporation after the mixing: a layer-0 parcel condenses inside layer 0,
        the cut stops the mixing that would dilute it, and the column traps its
        evaporation at the surface within about thirty steps.

        The surface group's vapour is mixed only over that subcloud layer; the
        layers above keep theirs, and a group cut this way gives its partial top
        no vapour. The heat groups are untouched. This is the lower bound on the
        transport, as complete mixing is the upper one.
        """
        vapour = self._inputs.get(VAPOUR_STATE)
        if vapour is None:
            return
        base = np.asarray(vapour.require_units("kg m^-2").values, dtype=float)
        phase = None if self.condensation_follows else self._inputs.get("vapour_phase_change_rate")
        if phase is not None:
            base = base + np.asarray(
                phase.require_units("kg m^-2 s^-1").values, dtype=float
            ) * float(phase.time_interval or self.context.dt)
        advection, first = composed_advection(self._inputs, self.advection_sources, "vapour")
        if advection is not None:
            base = base + advection * float(first.time_interval or self.context.dt)
        pressure = np.asarray(pressure, dtype=float)
        top = max((high for low, high in self._mixing_groups if low == 0), default=-1)
        T0, p0 = float(adjusted[0]), float(pressure[0])
        self._subcloud_parcel = (float("nan"),) * 4
        while top >= 1:
            q = float(np.sum(base[: top + 1]) / np.sum(layer_mass[: top + 1]))
            T_top = T0 * (float(pressure[top]) / p0) ** KAPPA
            if q < float(
                saturation_specific_humidity(np.array([T_top]), np.array([pressure[top]]))[0]
            ):
                break
            # Saturated by the top layer kept, so the level lies within the
            # layers already spanned and inside the relation's domain.
            level = lifting_condensation_pressure(T0, q, p0)
            self._subcloud_parcel = (T0, q, p0, level)
            lower = int(np.flatnonzero(pressure >= level)[-1])
            if lower >= top:
                break
            top = lower
        cut = max(top, 0)
        self._condensation_level_layer = cut
        groups, partials = [], dict(self._partial_tops)
        for low, high in self._mixing_groups:
            if low == 0 and high > cut:
                partials.pop(high + 1, None)
                if cut >= 1:
                    groups.append((0, cut))
                continue
            groups.append((low, high))
        self._vapour_groups, self._vapour_partials = tuple(groups), partials

    def _mix_classes(self, classes: np.ndarray, layer_mass: np.ndarray) -> np.ndarray:
        """Each provenance class mixed over this step's groups exactly as
        ``_vapour_transport`` mixes the vapour: uniform specific amount within a
        group, its total unchanged, the partial top taking its share."""
        mixed = np.array(classes, dtype=float)
        n = mixed.shape[1]
        for low, high in self._vapour_groups:
            span = slice(low, high + 1)
            top = high + 1
            f = self._vapour_partials.get(top, 0.0) if top < n else 0.0
            total = np.sum(classes[:, span], axis=1) + (f * classes[:, top] if f else 0.0)
            total_mass = float(np.sum(layer_mass[span])) + f * float(layer_mass[top] if f else 0.0)
            if total_mass <= 0.0:
                continue
            mixed[:, span] = np.outer(total, layer_mass[span]) / total_mass
            if f:
                mixed[:, top] = (1.0 - f) * classes[:, top] + f * total * layer_mass[top] / (
                    total_mass
                )
        return mixed

    def _group_fields(self, n: int) -> tuple[np.ndarray, np.ndarray]:
        """This step's groups as layer fields: each layer's group index or -1,
        and the fraction of each layer a group rose into. Output only: exactly
        what ``_vapour_transport`` mixed over, so the mixing can be replayed."""
        group = np.full(n, -1.0)
        for index, (low, high) in enumerate(self._vapour_groups):
            group[low : high + 1] = float(index)
        partial = np.zeros(n)
        for layer, fraction in self._vapour_partials.items():
            if 0 <= layer < n:
                partial[layer] = float(fraction)
        return group, partial

    def _surface_group_top(self) -> int:
        """The top layer of this step's surface-containing group, or 0 when the
        surface took no part in the mixing."""
        if not self._surface_mixed:
            return 0
        for low, high in self._mixing_groups:
            if low == 0:
                return int(high)
        return 0

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
        latent = None if self.condensation_follows else self._inputs.get("latent_heating_rate")
        if latent is not None:
            heating = heating + latent.require_units("W m^-2").values
        advective = composed_advection(self._inputs, self.advection_sources, "heating")[0]
        if advective is not None:
            heating = heating + advective

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
        partials: dict[int, float] | None = {} if self.fractional_top else None
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
                partial_tops=partials,
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
                partial_tops=partials,
            )
            self._surface_adjustment = 0.0
            self._surface_convective_flux = 0.0
            self._surface_mixed = False
        self._mixing_groups = groups
        self._partial_tops = dict(partials) if partials is not None else {}
        self._vapour_groups, self._vapour_partials = groups, self._partial_tops

        # The relaxed closure applies the exactly integrated fraction of the
        # same adjustment. One scalar across the whole column, so every group
        # still sums to zero and the enthalpy closes exactly, not on average.
        fraction = self.relaxation_fraction(ctx.dt)
        if fraction < 1.0:
            # Skipped, not multiplied by one: provisional + 1.0 * (adjusted -
            # provisional) is not bit-identical to adjusted, and a limit-cycling
            # column amplifies that round-off into a visible change.
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
        if self.vapour_mixing_top == "lifting_condensation_level":
            self._cut_at_condensation_level(np.asarray(adjusted, dtype=float), pressure, layer_mass)
        vapour_change = (
            self._vapour_transport(self._vapour_groups, layer_mass)
            if self.transport_vapour
            else None
        )
        if vapour_change is not None and fraction < 1.0:
            # A convex combination of the provisional and the mixed vapour, both
            # non-negative, so this cannot over-draw a layer either.
            vapour_change = fraction * vapour_change
        if vapour_change is not None:
            # A redistribution, so it sums to zero by construction. Anything
            # else would be convection creating or destroying water.
            self._vapour_moved = float(np.sum(np.abs(vapour_change)) / 2.0)
        self._mass_flux_tendency = None
        self._convective_condensation = None
        if self.mass_flux:
            vapour_change = self._apply_mass_flux(
                np.asarray(adjusted, dtype=float),
                heat_capacity,
                pressure,
                layer_mass,
                vapour_change,
            )
        self._vapour_change = vapour_change

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
                        state_id=ATMOSPHERE_TEMPERATURE_STATE,
                        provider_module_id=self.module_id,
                        label="mass_flux_transport",
                        value=self._mass_flux_tendency,
                        units="K s^-1",
                        equation_id=EQUATION_IDS["mass_flux_convection"],
                    ),
                )
                if self._mass_flux_tendency is not None
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

    def _apply_mass_flux(
        self,
        adjusted: np.ndarray,
        heat_capacity: np.ndarray,
        pressure: np.ndarray,
        layer_mass: np.ndarray,
        vapour_change: np.ndarray | None,
    ) -> np.ndarray:
        """The plume, its closure and its transport on the dry-adjusted profile.

        Returns the step's whole vapour redistribution (the dry groups' mixing
        plus the plume's transport), and keeps the plume's temperature tendency
        and its net condensation for the water module (the protocol's erratum 2).
        """
        ctx = self.context
        base = self._vapour_base
        if base is None or vapour_change is None:
            raise ConfigurationError("mass-flux convection needs the vapour")
        if abs(self._vapour_dt - ctx.dt) > 1.0e-9 * ctx.dt:
            raise ConfigurationError(
                "mass-flux convection composes with the water module over one interval"
            )
        vapour = base + vapour_change
        # The provenance classes start from the committed classes, which sum to
        # the committed vapour this module's base is, and are mixed over the
        # same groups (fall-path-closure-candidate.md).
        classes = None
        start = None
        if self.track_water_provenance:
            start = np.asarray(
                self.input_field("water_provenance").require_units("kg m^-2").values, dtype=float
            )
            classes = self._mix_classes(start, layer_mass)
        interfaces = np.asarray(
            self.input_field("atmosphere_interface_pressure").require_units("Pa").values,
            dtype=float,
        )
        common = dict(
            temperature=adjusted,
            vapour=vapour,
            layer_mass=layer_mass,
            heat_capacity=np.asarray(heat_capacity, dtype=float),
            pressure=np.asarray(pressure, dtype=float),
            interface_pressure=interfaces,
            latent_heat=value_of("water_latent_heat_vaporization_273"),
            gravity=value_of("earth_standard_gravity"),
            gas_constant=value_of("gas_constant_dry_air"),
        )
        result = convect(
            dt=ctx.dt,
            entrainment_rate=self.plume_entrainment_rate,
            evaporating_fraction=self.evaporating_fraction,
            adjustment_time=self.convective_adjustment_time,
            trigger=self.cape_trigger,
            detrained_state=self.detrained_state,
            precipitation_evaporation=self.precipitation_evaporation,
            condensate_pathway=self.condensate_pathway,
            conversion_rate=self.cloud_water_conversion_rate,
            provenance=classes,
            substep=self.convective_substep,
            continuation=self.convective_window_continuation,
            crossing=self.convective_trigger_crossing,
            cape_definition=self.convective_cape_definition,
            **common,
        )
        step = result.step
        self._plume = result.plume
        self._cloud_base_flux = result.mean_cloud_base_flux
        self._cape_removal = result.cape_removal
        self._surface_precipitation = step.surface_precipitation / ctx.dt
        self._mass_flux_substeps = step.substeps
        self._cape_after = result.cape_after
        self._cape_time_mean = result.cape_time_mean
        self._cape_host_start = result.cape
        self._cape_pre = result.cape_pre_convection
        self._cape_post = result.cape_post_convection
        if result.cape_window_start is not None:
            self._cape_window_start = np.array(result.cape_window_start)
        self._window_arrays = {
            name: np.array(values, dtype=float)
            for name, values in (
                ("cape_window_end", result.cape_window_end),
                ("convective_window_evaluated", result.window_evaluated),
                ("convective_window_status", result.window_status),
                ("convective_window_mass", result.window_mass),
                ("convective_window_fraction", result.window_fraction),
                ("cape_window_start_positive_area", result.window_start_positive_area),
                ("cape_window_start_maximum_work", result.window_start_maximum_work),
            )
            if values is not None
        }
        self._active_fraction = result.active_time_fraction
        self._unresolved_crossings = result.unresolved_crossings
        self._requested_flux = result.first_cloud_base_flux
        self._top_continuous = (
            result.plume.top_pressure_continuous if result.plume.active else float("nan")
        )
        if start is not None:
            if result.provenance is None or result.milestones is None:
                raise ConfigurationError("mass-flux convection lost the provenance classes")
            self._provenance_change = (result.provenance - start) / ctx.dt
            self._provenance_milestones = result.milestones.as_array()
        self._mass_flux_parts = {
            name: np.asarray(getattr(step, attribute), dtype=float) / ctx.dt
            for name, (attribute, _, _) in MASS_FLUX_DECOMPOSITION.items()
        }
        self._mass_flux_tendency = step.energy / np.asarray(heat_capacity, dtype=float) / ctx.dt
        self._convective_condensation = step.condensation / ctx.dt
        return vapour_change + step.transport

    def _mass_flux_outputs(self, ctx: Any) -> dict[str, Field]:
        if not self.mass_flux:
            return {}
        provenance = Provenance(
            producer=self.module_id, derivation=EQUATION_IDS["mass_flux_convection"]
        )
        n = np.asarray(self._adjusted).size
        plume = self._plume
        pressure = np.asarray(self._inputs["atmosphere_layer_pressure"].values, dtype=float)

        def layer(name: str, physical: str, values: np.ndarray, units: str) -> Field:
            return Field(
                field_id=name,
                physical_name=physical,
                values=np.asarray(values, dtype=float),
                units=units,
                quantity_kind=QuantityKind.DIAGNOSTIC,
                time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                time_interval=ctx.dt,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.LAYER_AVERAGE,
                dimensions=("layer",),
                missing_value_policy=MissingValuePolicy.MASKED,
                provenance=provenance,
            )

        def scalar(name: str, physical: str, value: float, units: str) -> Field:
            return Field(
                field_id=name,
                physical_name=physical,
                values=np.array(float(value)),
                units=units,
                quantity_kind=QuantityKind.DIAGNOSTIC,
                time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                time_interval=ctx.dt,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                # NaN only where the plume does not exist this step, and then it
                # says so: "no plume", not a missing number.
                missing_value_policy=(
                    MissingValuePolicy.MASKED
                    if not np.isfinite(value)
                    else MissingValuePolicy.NOT_PERMITTED
                ),
                provenance=provenance,
            )

        nan = np.full(n, np.nan)
        active = plume is not None and plume.active
        condensation = (
            self._convective_condensation
            if self._convective_condensation is not None
            else np.zeros(n)
        )
        return {
            "convective_condensation_rate": Field(
                field_id="convective_condensation_rate",
                physical_name="net convective condensation, formed less evaporated",
                values=np.asarray(condensation, dtype=float),
                units="kg m^-2 s^-1",
                quantity_kind=QuantityKind.TENDENCY,
                time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                time_interval=ctx.dt,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.LAYER_AVERAGE,
                dimensions=("layer",),
                provenance=provenance,
            ),
            "cloud_base_mass_flux": scalar(
                "cloud_base_mass_flux",
                "plume mass flux at cloud base",
                self._cloud_base_flux,
                "kg m^-2 s^-1",
            ),
            "dilute_cape": scalar(
                "dilute_cape",
                "the entraining plume's CAPE",
                plume.cape if plume else 0.0,
                "J kg^-1",
            ),
            "cape_removal_per_unit_mass_flux": scalar(
                "cape_removal_per_unit_mass_flux",
                "dilute CAPE removed per unit of cloud-base mass moved",
                self._cape_removal,
                "J m^2 kg^-2",
            ),
            "convective_surface_precipitation": scalar(
                "convective_surface_precipitation",
                "convective condensate reaching the surface",
                self._surface_precipitation,
                "kg m^-2 s^-1",
            ),
            "plume_temperature": layer(
                "plume_temperature",
                "the plume's temperature where it exists",
                plume.temperature if active else nan,
                "K",
            ),
            "plume_buoyancy": layer(
                "plume_buoyancy",
                "the plume's buoyancy with condensate loading",
                plume.buoyancy if active else nan,
                "m s^-2",
            ),
            "plume_top_pressure": scalar(
                "plume_top_pressure",
                "pressure of the plume's termination layer",
                float(pressure[plume.top]) if active else float("nan"),
                "Pa",
            ),
            "mass_flux_substeps": scalar(
                "mass_flux_substeps",
                "substeps the transport took",
                self._mass_flux_substeps,
                "1",
            ),
            # Labelled CAPE (convective-cape-sampling.md, user disposition):
            # the host-start sample carries a first-order splitting bias and is
            # disclosed as such; convergence uses the time mean.
            "cape_host_start": scalar(
                "cape_host_start",
                "the plume's CAPE at the step's first window, an operator-splitting diagnostic",
                self._cape_host_start,
                "J kg^-1",
            ),
            "cape_time_mean": scalar(
                "cape_time_mean",
                "the plume's CAPE averaged over the step, convection's idle remainder included",
                self._cape_time_mean,
                "J kg^-1",
            ),
            "cape_pre_convection": scalar(
                "cape_pre_convection",
                "the plume's CAPE at each convective window's start, averaged over the windows",
                self._cape_pre,
                "J kg^-1",
            ),
            "cape_post_convection": scalar(
                "cape_post_convection",
                "the plume's CAPE at each convective window's end, averaged over the windows",
                self._cape_post,
                "J kg^-1",
            ),
            "cape_host_end": scalar(
                "cape_host_end",
                "the plume's CAPE on the state convection leaves at the step's end",
                self._cape_after,
                "J kg^-1",
            ),
            "cape_window_start": Field(
                field_id="cape_window_start",
                physical_name="the plume's CAPE at the start of each convective window",
                values=np.asarray(self._cape_window_start, dtype=float),
                units="J kg^-1",
                quantity_kind=QuantityKind.DIAGNOSTIC,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                time_interval=ctx.dt,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                dimensions=("convective_window",),
                provenance=provenance,
            ),
            **{
                name: Field(
                    field_id=name,
                    physical_name=physical,
                    values=np.asarray(
                        self._window_arrays.get(name, np.zeros_like(self._cape_window_start)),
                        dtype=float,
                    ),
                    units=units,
                    quantity_kind=QuantityKind.DIAGNOSTIC,
                    time_semantics=semantics,
                    time_interval=ctx.dt,
                    timestamp=ctx.time,
                    spatial_support=SpatialSupport.GLOBAL_SCALAR,
                    dimensions=("convective_window",),
                    provenance=provenance,
                )
                for name, physical, units, semantics in (
                    (
                        "cape_window_end",
                        "the plume's CAPE at the end of each convective window",
                        "J kg^-1",
                        TimeSemantics.INSTANTANEOUS,
                    ),
                    (
                        "convective_window_evaluated",
                        "1 where the closure was evaluated at the window's start, 0 where skipped",
                        "1",
                        TimeSemantics.INSTANTANEOUS,
                    ),
                    (
                        "convective_window_status",
                        "each window's outcome: 0 active, 1 stopped within, 2 below the "
                        "trigger, 3 no plume, 4 no CAPE removal, 5 skipped",
                        "1",
                        TimeSemantics.INSTANTANEOUS,
                    ),
                    (
                        "convective_window_mass",
                        "the cloud-base mass each window moved",
                        "kg m^-2",
                        TimeSemantics.INTERVAL_ACCUMULATION,
                    ),
                    (
                        "cape_window_start_positive_area",
                        "the plume's positive buoyant area over its whole extent at each "
                        "window's start (D1)",
                        "J kg^-1",
                        TimeSemantics.INSTANTANEOUS,
                    ),
                    (
                        "cape_window_start_maximum_work",
                        "the largest integrated buoyancy work the plume reaches, at each "
                        "window's start (D2)",
                        "J kg^-1",
                        TimeSemantics.INSTANTANEOUS,
                    ),
                    (
                        "convective_window_fraction",
                        "the fraction of its substep a window convected for: f where limited to "
                        "the trigger, 1 where unlimited, 0 where it did not convect",
                        "1",
                        TimeSemantics.INSTANTANEOUS,
                    ),
                )
            },
            "convective_trigger_unresolved_crossings": scalar(
                "convective_trigger_unresolved_crossings",
                "fractional substeps whose CAPE could not be brought to the trigger within "
                "tolerance, the crossing being discontinuous",
                float(self._unresolved_crossings),
                "1",
            ),
            "convective_active_time_fraction": scalar(
                "convective_active_time_fraction",
                "the share of the step convection integrated with a positive mass flux",
                self._active_fraction,
                "1",
            ),
            "requested_cloud_base_mass_flux": scalar(
                "requested_cloud_base_mass_flux",
                "the closure's cloud-base mass flux at the step's first substep",
                self._requested_flux,
                "kg m^-2 s^-1",
            ),
            "plume_top_pressure_continuous": scalar(
                "plume_top_pressure_continuous",
                "where the plume's buoyancy work crosses zero, interpolated between layers",
                self._top_continuous,
                "Pa",
            ),
            **{
                name: Field(
                    field_id=name,
                    physical_name=physical,
                    values=np.asarray(self._mass_flux_parts.get(name, np.zeros(n)), dtype=float),
                    units=units,
                    quantity_kind=QuantityKind.DIAGNOSTIC,
                    time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                    time_interval=ctx.dt,
                    timestamp=ctx.time,
                    spatial_support=SpatialSupport.LAYER_AVERAGE,
                    dimensions=("layer",),
                    provenance=provenance,
                )
                for name, (_, physical, units) in MASS_FLUX_DECOMPOSITION.items()
            },
            **(self._provenance_outputs(ctx, provenance) if self.track_water_provenance else {}),
        }

    def _provenance_outputs(self, ctx: Any, provenance: Any) -> dict[str, Field]:
        n = np.asarray(self._adjusted).size
        change = (
            self._provenance_change
            if self._provenance_change is not None
            else np.zeros((2 * n + 1, n))
        )
        crossed = (
            self._provenance_milestones
            if self._provenance_milestones is not None
            else np.zeros((5, n))
        )
        common = dict(
            quantity_kind=QuantityKind.DIAGNOSTIC,
            time_semantics=TimeSemantics.INTERVAL_AVERAGE,
            time_interval=ctx.dt,
            timestamp=ctx.time,
            spatial_support=SpatialSupport.LAYER_AVERAGE,
            provenance=provenance,
        )
        return {
            "water_provenance_change_rate": Field(
                field_id="water_provenance_change_rate",
                physical_name="convection's change to each provenance class of the vapour",
                values=np.asarray(change, dtype=float),
                units="kg m^-2 s^-1",
                dimensions=("provenance_class", "layer"),
                **common,
            ),
            "convective_provenance_milestones": Field(
                field_id="convective_provenance_milestones",
                physical_name=(
                    "first condensation, first precipitation by origin and by location, and "
                    "first surface arrival, crossed in convection this step"
                ),
                values=np.asarray(crossed, dtype=float),
                units="kg m^-2",
                dimensions=("milestone", "layer"),
                **common,
            ),
        }

    def emit_outputs(self) -> dict[str, Field]:
        ctx = self.context
        provenance = Provenance(
            producer=self.module_id, derivation=EQUATION_IDS["convective_adjustment"]
        )
        return {
            **self._mass_flux_outputs(ctx),
            "convective_temperature_tendency": Field(
                field_id="convective_temperature_tendency",
                physical_name="temperature tendency from convective mixing",
                symbol="dT/dt|conv",
                values=np.array(
                    self._tendency
                    if self._mass_flux_tendency is None
                    else self._tendency + self._mass_flux_tendency
                ),
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
            # The mixing as a rate, for the water module to condense what it
            # lifts in the same step. Over the interval the vapour was formed
            # for, so the water module's composition reproduces it exactly.
            **(
                {
                    "convective_vapour_transport_rate": Field(
                        field_id="convective_vapour_transport_rate",
                        physical_name="vapour redistribution by the convective mixing",
                        values=np.array(
                            self._vapour_change / self._vapour_dt
                            if self._vapour_change is not None and self._vapour_dt > 0.0
                            else np.zeros(np.asarray(self._adjusted).size)
                        ),
                        units="kg m^-2 s^-1",
                        quantity_kind=QuantityKind.TENDENCY,
                        time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                        time_interval=self._vapour_dt or ctx.dt,
                        timestamp=ctx.time,
                        spatial_support=SpatialSupport.LAYER_AVERAGE,
                        dimensions=("layer",),
                        provenance=provenance,
                    )
                }
                if self.condensation_follows
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
            **(
                {
                    # The top of the group the surface mixed with this step, as a
                    # layer index; the lowest layer when the surface took no part
                    # (descending-boundary-layer-ventilation.md). Output only.
                    "surface_mixed_layer_top": Field(
                        field_id="surface_mixed_layer_top",
                        physical_name="top layer of the adjustment group containing the surface",
                        values=np.array(float(self._surface_group_top())),
                        units="1",
                        quantity_kind=QuantityKind.DIAGNOSTIC,
                        time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                        time_interval=ctx.dt,
                        timestamp=ctx.time,
                        spatial_support=SpatialSupport.GLOBAL_SCALAR,
                        provenance=provenance,
                    )
                }
                if self.surface_based
                else {}
            ),
            **(
                {
                    name: Field(
                        field_id=name,
                        physical_name=label,
                        values=np.array(value),
                        units=units,
                        quantity_kind=QuantityKind.DIAGNOSTIC,
                        time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                        time_interval=ctx.dt,
                        timestamp=ctx.time,
                        spatial_support=SpatialSupport.GLOBAL_SCALAR,
                        # NaN where the parcel was not lifted this step.
                        missing_value_policy=MissingValuePolicy.MASKED,
                        provenance=provenance,
                    )
                    for name, label, units, value in zip(
                        (
                            "subcloud_parcel_temperature",
                            "subcloud_parcel_specific_humidity",
                            "subcloud_parcel_pressure",
                            "subcloud_condensation_pressure",
                        ),
                        (
                            "temperature of the parcel the subcloud cut lifted",
                            "mean specific humidity of the layers that parcel stands for",
                            "pressure the parcel was lifted from",
                            "lifting condensation pressure of that parcel",
                        ),
                        ("K", "kg kg^-1", "Pa", "Pa"),
                        self._subcloud_parcel,
                        strict=True,
                    )
                }
                if self.vapour_mixing_top == "lifting_condensation_level"
                else {}
            ),
            **(
                {
                    name: Field(
                        field_id=name,
                        physical_name=label,
                        values=values,
                        units="1",
                        quantity_kind=QuantityKind.DIAGNOSTIC,
                        time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                        time_interval=ctx.dt,
                        timestamp=ctx.time,
                        spatial_support=SpatialSupport.LAYER_AVERAGE,
                        dimensions=("layer",),
                        provenance=provenance,
                    )
                    for name, label, values in zip(
                        ("convective_mixing_group", "convective_partial_top_fraction"),
                        (
                            "index of the adjustment group each layer was mixed in, or -1",
                            "fraction of each layer an adjustment group rose into",
                        ),
                        self._group_fields(len(self._adjusted)),
                        strict=True,
                    )
                }
                if self.diagnose_mixing_groups
                else {}
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

    def advection_fields(self, kind: str) -> tuple[str, ...]:
        """The fields this module composes with for one kind of rate."""
        if not self.receives_advection:
            return ()
        return tuple(advection_field(s, kind) for s in self.advection_sources)

    def advection_rate(self, kind: str) -> np.ndarray | None:
        """This step's advected rate of one kind, every declared transport summed."""
        if not self.receives_advection:
            return None
        return composed_advection(self._inputs, self.advection_sources, kind)[0]

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "target": self.target,
            "advection_sources": self.advection_fields("vapour"),
            "bulk_plume_parameter": (
                self.bulk_plume_parameter if self.target == "entraining_plume" else None
            ),
            "target_responds_to_temperature": self.target_responds_to_temperature,
            "surface_based": self.surface_based,
            "transport_vapour": self.transport_vapour,
            **({"condensation_follows": True} if self.condensation_follows else {}),
            "surface_pressure_Pa": (
                self._reference_surface_pressure() if self.surface_based else None
            ),
            "sensible_heat_coefficient_W_m2_K": self.sensible_heat_coefficient,
            "reference_pressure_Pa": self.reference_pressure,
            "relaxation_timescale_s": self.relaxation_timescale,
            "fractional_top": self.fractional_top,
            **(
                {"vapour_mixing_top": self.vapour_mixing_top}
                if self.vapour_mixing_top != "group"
                else {}
            ),
            **(
                {
                    "mass_flux": True,
                    "plume_entrainment_rate_per_m": self.plume_entrainment_rate,
                    "convective_adjustment_time_s": self.convective_adjustment_time,
                    "cape_trigger_J_kg": self.cape_trigger,
                    "evaporating_fraction": self.evaporating_fraction,
                    "detrained_state": self.detrained_state,
                    "precipitation_evaporation": self.precipitation_evaporation,
                    "condensate_pathway": self.condensate_pathway,
                    **({"tracks_water_provenance": True} if self.track_water_provenance else {}),
                    **(
                        {"convective_substep_s": self.convective_substep}
                        if self.convective_substep is not None
                        else {}
                    ),
                    **(
                        {"convective_window_continuation": self.convective_window_continuation}
                        if self.convective_window_continuation != "stop"
                        else {}
                    ),
                    **(
                        {"convective_trigger_crossing": self.convective_trigger_crossing}
                        if self.convective_trigger_crossing != "explicit"
                        else {}
                    ),
                    **(
                        {"convective_cape_definition": self.convective_cape_definition}
                        if self.convective_cape_definition != "first_lnb"
                        else {}
                    ),
                    **(
                        {"cloud_water_conversion_rate_per_m": self.cloud_water_conversion_rate}
                        if self.condensate_pathway == "carried"
                        else {}
                    ),
                    **(
                        {
                            "precipitation_evaporation_efficiency": (
                                PRECIPITATION_EVAPORATION_EFFICIENCY
                            )
                        }
                        if self.precipitation_evaporation in ("falling_rate_law", "rate_law_only")
                        else {}
                    ),
                }
                if self.mass_flux
                else {}
            ),
            "mixing": (
                "relaxed_toward_complete_within_unstable_groups"
                if self.relaxation_timescale > 0.0
                else "instantaneous_complete_within_unstable_groups"
            ),
        }
