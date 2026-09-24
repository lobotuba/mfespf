"""Reduced ocean heat storage: the deep reservoir beneath the mixed layer.

Charter 9.10 and section 19 stage 5; requirements OCN-001, OCN-002, STA-001,
CON-002, CON-009, TIM-008, KER-012.

Every configuration up to Stage 4 stored heat in exactly one place, and every one
of them carried the same declared limitation: *one temperature and one timescale
by construction*.  A single reservoir relaxes as a single exponential, so it
cannot exhibit the fast surface adjustment and slow continuing drift that the
charter's timescale philosophy (11.1) expects a real system to show.  This module
is the smallest addition that removes that limitation honestly.

**It is a second module, not a second number inside the surface.**  The mixed
layer is still owned by the surface module; this one owns a deep reservoir and
nothing else.  The heat that moves between them is a declared transfer, recorded
independently by both parties and reconciled by the ledger like any other
(charter 8.4).  The alternative -- giving the surface module two temperatures --
would have hidden the exchange inside one module where no audit could see it, and
would have made the deep reservoir impossible to disable as a single edge
(CAU-007).

**What this is not.**  The fidelity path in charter 9.10 runs from a fixed
boundary temperature through mixed-layer and deep reservoirs to a
three-dimensional dynamic ocean.  This module is the third rung.  It has no
salinity, no circulation, no meridional transport and no sea ice, and it does not
pretend that its exchange coefficient is a measured property of seawater: it is a
fitted number, classified accordingly, and any claim resting on it inherits that.
"""

from __future__ import annotations

import math
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
from mfespf.registries.states import (
    StateBudget,
    StateEntry,
    StateKind,
    TendencyContribution,
    UpdateAuthority,
)

__all__ = [
    "OCEAN_DEEP_STATE",
    "OCEAN_INTERFACE",
    "DeepOceanModule",
    "areal_heat_capacity_of_water",
    "deep_ocean_heat_flux",
    "ocean_deep_state_entry",
    "two_reservoir_modes",
    "two_reservoir_transfer_function",
]

V1 = Version(1, 0, 0)

OCEAN_DEEP_STATE = "ocean_deep_temperature"

OCEAN_INTERFACE = InterfaceContract(
    interface_id="ocean",
    version=V1,
    physical_semantics={
        "temperature_units": "K",
        "flux_units": "W m^-2",
        "flux_sign_convention": "positive_downward_from_mixed_layer_into_the_deep_reservoir",
        "heat_content_units": "J m^-2",
        "spatial_support": "global_scalar",
        "mixed_layer_owner": "surface",
        "deep_reservoir_owner": "ocean",
    },
    optional_capabilities=frozenset({"salinity", "circulation_transport", "sea_ice_exchange"}),
)


def ocean_deep_state_entry(
    *, owner: str = "ocean", readers: frozenset[str] = frozenset()
) -> StateEntry:
    """The registry row for the deep-reservoir temperature (STA-002, STA-003)."""
    return StateEntry(
        record_id=OCEAN_DEEP_STATE,
        version=V1,
        physical_name="temperature of the reduced deep-ocean reservoir",
        symbol="T_d",
        owner_module_id=owner,
        units="K",
        state_kind=StateKind.PROGNOSTIC,
        conservation_participation=("energy",),
        authorized_readers=readers,
        authorized_tendency_providers=frozenset({owner}),
        update_authority=UpdateAuthority.OWNER_ONLY,
        initialization_source="manifest",
        checkpoint_required=True,
        uncertainty_representation="initial_condition_and_parameter_uncertainty",
    )


# --- physics ---------------------------------------------------------------


def areal_heat_capacity_of_water(depth_m: float) -> float:
    """Areal heat capacity ``rho c_p h`` of a well-mixed water column, J m^-2 K^-1.

    Depth is the parameter a reader can check against a measurement; the product
    is not.  Configurations therefore declare a depth and let this derive the
    capacity, so two runs quoting different capacities can be compared without
    reverse-engineering which density and specific heat each assumed.
    """
    if depth_m <= 0.0:
        raise ConfigurationError(f"water depth {depth_m} m must be positive")
    return value_of("density_seawater") * value_of("specific_heat_seawater") * depth_m


def deep_ocean_heat_flux(
    mixed_layer_temperature: float, deep_temperature: float, exchange_coefficient: float
) -> float:
    """``H = gamma (T_s - T_d)``, positive downward, W m^-2.

    A linear exchange law is the reduced representation of everything that
    actually carries heat downward -- subduction, diffusion, overturning -- none
    of which is resolved here.  Its one defensible property is the limiting case
    that matters: when the two reservoirs are at the same temperature the
    exchange is exactly zero, so a spun-up column cannot drift because of the
    ocean's presence.
    """
    return exchange_coefficient * (mixed_layer_temperature - deep_temperature)


def two_reservoir_modes(
    *,
    mixed_layer_capacity: float,
    deep_capacity: float,
    feedback_parameter: float,
    exchange_coefficient: float,
) -> tuple[float, float]:
    """Analytic relaxation timescales of the linearised two-reservoir system.

    For anomalies about a radiative equilibrium,

        C_s dT_s/dt = F - lambda T_s - gamma (T_s - T_d)
        C_d dT_d/dt =                  gamma (T_s - T_d)

    which is linear, so its response is a sum of two decaying exponentials whose
    timescales are minus the reciprocal eigenvalues of the coefficient matrix.
    Returned fast first.

    This exists to be *compared against*, not used: the modules integrate the
    nonlinear system, and VAL-002 asks for an analytic comparison wherever a
    limiting case admits one.  A measured timescale that disagrees with this by
    more than the discretisation error is a defect in the coupling, not a
    discovery about the ocean.
    """
    for name, value in (
        ("mixed_layer_capacity", mixed_layer_capacity),
        ("deep_capacity", deep_capacity),
        ("feedback_parameter", feedback_parameter),
        ("exchange_coefficient", exchange_coefficient),
    ):
        if value <= 0.0:
            raise ConfigurationError(f"{name} must be positive, got {value}")

    a = (feedback_parameter + exchange_coefficient) / mixed_layer_capacity
    b = exchange_coefficient / mixed_layer_capacity
    c = exchange_coefficient / deep_capacity

    # Eigenvalues of [[-a, b], [c, -c]]; both roots are real and negative because
    # the discriminant below is a sum of squares and a positive cross term.
    trace = -(a + c)
    determinant = a * c - b * c
    discriminant = math.sqrt(trace * trace - 4.0 * determinant)
    slow_eigenvalue = 0.5 * (trace + discriminant)
    fast_eigenvalue = 0.5 * (trace - discriminant)
    return (-1.0 / fast_eigenvalue, -1.0 / slow_eigenvalue)


def two_reservoir_transfer_function(
    angular_frequency: float,
    *,
    mixed_layer_capacity: float,
    deep_capacity: float,
    feedback_parameter: float,
    exchange_coefficient: float,
) -> complex:
    """Complex gain of the mixed-layer temperature per unit absorbed forcing.

    Charter 11.3 asks for frequency-dependent gain and phase. Eliminating the
    deep temperature from the linearised pair gives, at angular frequency w,

        T_s(w) / F(w) = 1 / (i w C_s + lambda + gamma i w C_d / (i w C_d + gamma))

    whose magnitude is the gain in K per W m^-2 and whose negated argument is the
    phase lag in radians.

    The last term is what distinguishes a two-reservoir system from a slab, and
    it does two things at once: it raises the effective damping at high frequency
    -- the deep store is another place for heat to go -- while contributing a
    zero that *reduces* the phase lag. A measured lag smaller than a slab's is
    therefore a prediction of this model, not a sign that the coupling ran
    backwards, and it is worth asserting for exactly that reason.

    As ``exchange_coefficient`` approaches zero this reduces to the single-pole
    slab response ``1 / (i w C_s + lambda)``.
    """
    if angular_frequency <= 0.0:
        raise ConfigurationError("angular frequency must be positive")
    for name, value in (
        ("mixed_layer_capacity", mixed_layer_capacity),
        ("deep_capacity", deep_capacity),
        ("feedback_parameter", feedback_parameter),
    ):
        if value <= 0.0:
            raise ConfigurationError(f"{name} must be positive, got {value}")
    if exchange_coefficient < 0.0:
        raise ConfigurationError("exchange coefficient must be non-negative")

    s = 1j * angular_frequency
    deep_branch = (
        exchange_coefficient * s * deep_capacity / (s * deep_capacity + exchange_coefficient)
        if exchange_coefficient > 0.0
        else 0.0
    )
    return 1.0 / (s * mixed_layer_capacity + feedback_parameter + deep_branch)


class DeepOceanModule(Module):
    """A single deep reservoir exchanging heat with the mixed layer above it.

    The second half of the "reduced ocean heat storage" configuration KER-012
    names.  The first half is whichever surface module owns the mixed layer: this
    one adds a reservoir, not a surface.
    """

    def __init__(
        self,
        *,
        module_id: str = "ocean",
        surface_module_id: str = "surface",
        depth: float = 1200.0,
        areal_heat_capacity: float | None = None,
        exchange_coefficient: float = 0.7,
        initial_temperature: float = 250.0,
        status: ModuleStatus = ModuleStatus.ACTIVE,
    ) -> None:
        contract = ModuleContract(
            module_id=module_id,
            module_name="reduced two-reservoir ocean: deep store below the mixed layer",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.REDUCED,
            physical_basis=(
                "First law applied to a second well-mixed water reservoir that exchanges "
                "heat with the mixed layer at a rate linear in their temperature "
                "difference. The mixed layer itself is owned by the surface module; this "
                "module adds storage and a downward heat pathway, not a surface."
            ),
            governing_equation_ids=(
                EQUATION_IDS["deep_ocean_exchange"],
                EQUATION_IDS["slab"],
                EQUATION_IDS["forward_euler"],
            ),
            equation_variables=(
                (
                    EQUATION_IDS["slab"],
                    (
                        # Both of these were wrong before they were stated. The
                        # temperature was refused as ambiguous, which was the
                        # check working. The capacity was not: the mixed layer's
                        # surface_heat_capacity is a declared input of this
                        # module and carries the symbol C in the right units, so
                        # it was the *only* candidate and inference reported
                        # confidently on the wrong reservoir. A unique match is
                        # not evidence of a correct one -- the ambiguity refusal
                        # only fires when a second candidate exists to notice.
                        ("T", OCEAN_DEEP_STATE, "K"),
                        ("C", "areal_heat_capacity_J_m2_K", "J m^-2 K^-1"),
                    ),
                ),
            ),
            assumptions=(
                "the deep reservoir is isothermal and well mixed",
                "heat exchange is linear in the mixed-layer to deep temperature difference",
                "the exchange coefficient is constant in time and independent of state",
                "no salinity, circulation, meridional transport or sea-ice exchange",
                "no freshwater, momentum or chemical coupling at this fidelity",
            ),
            validity_domain=(
                "Global-mean heat uptake on timescales from the mixed-layer adjustment time "
                "to a few times the deep-reservoir time. Nothing about spatial patterns of "
                "heat uptake, about circulation change, or about sea level is in scope."
            ),
            owned_state_ids=frozenset({OCEAN_DEEP_STATE}),
            required_inputs=("surface_temperature", "surface_heat_capacity"),
            produced_outputs=(
                OCEAN_DEEP_STATE,
                "deep_ocean_heat_flux",
                "ocean_heat_uptake",
                "ocean_deep_stored_heat",
                "ocean_total_stored_heat",
                "ocean_exchange_coefficient",
            ),
            imported_transfers=("deep_ocean_heat",),
            exported_transfers=(),
            conserved_quantities=("energy",),
            provided_interfaces=(OCEAN_INTERFACE,),
            required_interfaces=(("surface", VersionRange.parse(">=1.0.0,<2.0.0")),),
            uncertainty_outputs=(
                "initial_condition",
                "deep_reservoir_depth_parameter",
                "exchange_coefficient_calibrated_parameter",
                "time_discretization",
                "linear_exchange_law_structural",
            ),
            correctness_tests=(
                "test_equal_temperatures_exchange_no_heat",
                "test_two_reservoir_timescales_match_the_analytic_modes",
                "test_deep_reservoir_conserves_energy_with_the_surface",
                "test_disabling_the_exchange_recovers_the_single_reservoir_result",
            ),
            known_limitations=(
                "The exchange coefficient is a calibrated parameter, not a measured "
                "property; a claim about heat-uptake magnitude inherits that provenance.",
                "One deep temperature, so the vertical structure of ocean heat uptake is "
                "not represented and no claim about it is admissible.",
                "No salinity, circulation, meridional transport, sea ice, freshwater, "
                "momentum or chemical exchange, all of which charter 9.10 lists for the "
                "full-fidelity module.",
                "Heat uptake is proportional to the temperature difference at all times, "
                "so the efficiency cannot change as a real ocean's does under warming.",
            ),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(1)",
                memory_complexity="O(1)",
                scaling_variables=(),
                work_units_per_step=1.0,
                benchmark_configuration="global scalar, one reservoir",
                peak_memory_bytes=4096,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
                lower_fidelity_alternatives=("fixed_deep_boundary_temperature",),
                checkpoint_cost_bytes=8,
            ),
        )
        super().__init__(contract, status=status)
        self.surface_module_id = surface_module_id
        self.depth = depth
        self.areal_heat_capacity = (
            areal_heat_capacity
            if areal_heat_capacity is not None
            else areal_heat_capacity_of_water(depth)
        )
        self.exchange_coefficient = exchange_coefficient
        self.initial_temperature = initial_temperature

        self.temperature = initial_temperature
        self._mixed_layer_temperature = initial_temperature
        self._mixed_layer_capacity = 0.0
        self._flux = 0.0
        self._interval = (0.0, 0.0)
        self._stored_initial = 0.0
        self._budget: StateBudget | None = None

    # --- configuration ----------------------------------------------------
    def validate_configuration(self) -> None:
        super().validate_configuration()
        if self.areal_heat_capacity <= 0.0:
            raise ConfigurationError(
                f"deep-reservoir areal heat capacity {self.areal_heat_capacity} must be positive"
            )
        if self.exchange_coefficient < 0.0:
            raise ConfigurationError(
                f"exchange coefficient {self.exchange_coefficient} must be non-negative; a "
                f"negative value would carry heat up the temperature gradient"
            )
        if self.initial_temperature <= 0.0:
            raise ConfigurationError("initial deep temperature must be positive")

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        self.temperature = float(initial_state.get(OCEAN_DEEP_STATE, self.initial_temperature))
        self.mark_initialized()

    # --- physics ----------------------------------------------------------
    def stored_energy(self, temperature: float | None = None) -> float:
        T = self.temperature if temperature is None else temperature
        return self.areal_heat_capacity * T

    def compute_tendencies(self) -> list[TendencyContribution]:
        ctx = self.context
        self._mixed_layer_temperature = (
            self.input_field("surface_temperature").require_units("K").scalar()
        )
        self._mixed_layer_capacity = (
            self.input_field("surface_heat_capacity").require_units("J m^-2 K^-1").scalar()
        )
        self._flux = deep_ocean_heat_flux(
            self._mixed_layer_temperature, self.temperature, self.exchange_coefficient
        )
        self._interval = (ctx.time, ctx.time + ctx.dt)
        self._stored_initial = self.stored_energy()

        return [
            TendencyContribution(
                state_id=OCEAN_DEEP_STATE,
                provider_module_id=self.module_id,
                label="mixed_layer_exchange",
                value=self._flux / self.areal_heat_capacity,
                units="K s^-1",
                equation_id=EQUATION_IDS["deep_ocean_exchange"],
            )
        ]

    def advance_state(self, dt: float, contributions: list[TendencyContribution]) -> None:
        initial = np.asarray(self.temperature, dtype=float)
        total = sum(float(np.asarray(c.value)) for c in contributions)
        self.temperature = float(initial + total * dt)
        self._budget = StateBudget(
            state_id=OCEAN_DEEP_STATE,
            initial_value=initial,
            final_value=np.asarray(self.temperature, dtype=float),
            dt=dt,
            contributions=list(contributions),
        )

    def state_budgets(self) -> dict[str, Any]:
        return {OCEAN_DEEP_STATE: self._budget} if self._budget else {}

    # --- outputs ----------------------------------------------------------
    def emit_state(self) -> dict[str, Field]:
        """State and static properties, published before the first step.

        The deep temperature is committed state, so it belongs here rather than
        in ``emit_outputs``: it has an instant, not an averaging interval.
        """
        timestamp = self._context.time if self._context else 0.0
        provenance = Provenance(
            producer=self.module_id, derivation=EQUATION_IDS["deep_ocean_exchange"]
        )
        return {
            OCEAN_DEEP_STATE: Field(
                field_id=OCEAN_DEEP_STATE,
                physical_name="temperature of the reduced deep-ocean reservoir",
                symbol="T_d",
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
                        "deep_reservoir_depth_parameter",
                        "time_discretization",
                    ),
                ),
            ),
            "ocean_deep_stored_heat": Field(
                field_id="ocean_deep_stored_heat",
                physical_name="heat stored in the deep reservoir",
                symbol="C_d T_d",
                values=np.array(self.stored_energy()),
                units="J m^-2",
                quantity_kind=QuantityKind.EXTENSIVE_STATE,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                timestamp=timestamp,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                provenance=provenance,
            ),
            "ocean_exchange_coefficient": Field(
                field_id="ocean_exchange_coefficient",
                physical_name="mixed-layer to deep-reservoir heat exchange coefficient",
                symbol="gamma",
                values=np.array(self.exchange_coefficient),
                units="W m^-2 K^-1",
                quantity_kind=QuantityKind.BOUNDARY_VALUE,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                timestamp=timestamp,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                provenance=Provenance(
                    producer=self.module_id,
                    notes=(
                        "calibrated parameter, not a measured property of seawater "
                        "(charter 4.13, EQR-004)"
                    ),
                ),
            ),
        }

    def emit_outputs(self) -> dict[str, Field]:
        ctx = self.context
        provenance = Provenance(
            producer=self.module_id, derivation=EQUATION_IDS["deep_ocean_exchange"]
        )
        flux_common = {
            "units": "W m^-2",
            "quantity_kind": QuantityKind.FLUX,
            "time_semantics": TimeSemantics.INTERVAL_AVERAGE,
            "time_interval": ctx.dt,
            "timestamp": ctx.time,
            "spatial_support": SpatialSupport.GLOBAL_SCALAR,
            "sign_convention": SignConvention.POSITIVE_DOWNWARD,
            "provenance": provenance,
        }
        return {
            "deep_ocean_heat_flux": Field(
                field_id="deep_ocean_heat_flux",
                physical_name="heat flux from the mixed layer into the deep reservoir",
                symbol="H",
                values=np.array(self._flux),
                **flux_common,  # type: ignore[arg-type]
            ),
            # Ocean heat uptake is the same number seen as an observable rather
            # than as a coupling term: it is what an observing system estimates
            # from the change in ocean heat content, and OCN-002 asks for it by
            # name. Published separately so a reader comparing against an
            # observational product is not reading a coupling diagnostic.
            "ocean_heat_uptake": Field(
                field_id="ocean_heat_uptake",
                physical_name="rate of heat accumulation in the deep reservoir",
                symbol="N_deep",
                values=np.array(self._flux),
                **flux_common,  # type: ignore[arg-type]
            ),
            "ocean_total_stored_heat": Field(
                field_id="ocean_total_stored_heat",
                physical_name="heat stored in the mixed layer and the deep reservoir together",
                symbol="OHC",
                values=np.array(
                    self._mixed_layer_capacity * self._mixed_layer_temperature
                    + self.stored_energy()
                ),
                units="J m^-2",
                quantity_kind=QuantityKind.DIAGNOSTIC,
                time_semantics=TimeSemantics.START_OF_STEP,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                provenance=Provenance(
                    producer=self.module_id,
                    derivation=EQUATION_IDS["deep_ocean_exchange"],
                    notes=(
                        "a diagnostic sum across two owners: the mixed-layer term is read "
                        "from the surface module's published state and is not stored here "
                        "(STA-001)"
                    ),
                ),
            ),
        }

    def report_transfers(self) -> tuple[list[Transfer], list[Transfer]]:
        """The deep reservoir imports what the mixed layer exports.

        Recorded as one signed transaction rather than as separate downward and
        upward channels. The sign convention is fixed at positive into the
        receiver (CON-003), so a cooling mixed layer produces a negative import
        here and a negative export there, and the pair still reconciles. Splitting
        it by direction would create two channels of which one is always zero and
        would make a sign error look like a missing record.
        """
        start, end = self._interval
        dt = end - start
        return (
            [],
            [
                Transfer(
                    transfer_id=f"deep_ocean_heat@{start:g}",
                    sender_module=self.surface_module_id,
                    receiver_module=self.module_id,
                    quantity=ConservedQuantity.ENERGY,
                    amount=self._flux * dt,
                    units="J m^-2",
                    start_time=start,
                    end_time=end,
                    spatial_mapping="global_mean_over_planetary_surface",
                    channel="deep_ocean_heat",
                    equation_id=EQUATION_IDS["deep_ocean_exchange"],
                )
            ],
        )

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
            quantified={
                "exchange_coefficient_W_m2_K": self.exchange_coefficient,
                "deep_reservoir_depth_m": self.depth,
            },
            method="none propagated at this fidelity; the parameters are declared, not sampled",
            distributional_assumptions="unresolved",
            correlations="unresolved",
            unresolved=(
                "initial_deep_temperature_uncertainty",
                "deep_reservoir_depth_parameter_uncertainty",
                "exchange_coefficient_is_calibrated_not_measured",
                "forward_euler_time_discretization_error",
                "structural_uncertainty_of_a_linear_state_independent_exchange_law",
            ),
        )

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "deep_reservoir_depth_m": self.depth,
            "areal_heat_capacity_J_m2_K": self.areal_heat_capacity,
            "exchange_coefficient_W_m2_K": self.exchange_coefficient,
            "exchange_coefficient_parameter_class": "calibrated",
            "initial_temperature_K": self.initial_temperature,
            "time_integration": "forward_euler",
        }

    def checkpoint(self) -> dict[str, Any]:
        return {OCEAN_DEEP_STATE: self.temperature}

    def restore(self, checkpoint: dict[str, Any]) -> None:
        self.temperature = float(checkpoint[OCEAN_DEEP_STATE])
