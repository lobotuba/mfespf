"""Hydrostatic atmospheric column.

Charter 9.4 and section 19 stage 2; requirements ATM-001, ATM-002, ATM-003,
STA-001, STA-006, CON-009.

The column owns the atmospheric thermodynamic state and nothing else.  It does
not compute radiation: the radiation module submits labeled per-layer heating
tendencies and the column commits them, which is the A1.3 pattern for a state
several processes act on.  At this fidelity there is no dynamics, so the layer
temperatures evolve under radiative heating alone and the module's known
limitations say so.

Mass is distributed hydrostatically in pressure coordinates: a layer spanning
``dp`` holds ``dp / g`` kilograms per square metre, independent of temperature.
That is exact for a hydrostatic atmosphere and is what makes the column's heat
capacity a fixed, declarable quantity rather than something that drifts with the
state.
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
    Provenance,
    QuantityKind,
    SignConvention,
    SpatialSupport,
    TimeSemantics,
    Uncertainty,
)
from mfespf.core.versioning import InterfaceContract, Version
from mfespf.coupling.conservation import ModuleBudget
from mfespf.coupling.transfers import ConservedQuantity, Transfer
from mfespf.physics.equations import EQUATION_IDS
from mfespf.physics.gravity import GravityField, as_gravity_field
from mfespf.physics.thermodynamics import hydrostatic_heights
from mfespf.registries.states import (
    StateBudget,
    StateEntry,
    StateKind,
    TendencyContribution,
    UpdateAuthority,
)

__all__ = [
    "ATMOSPHERE_INTERFACE",
    "ATMOSPHERE_TEMPERATURE_STATE",
    "AtmosphericColumnModule",
    "atmosphere_state_entry",
    "hydrostatic_layer_masses",
    "layer_interface_pressures",
]

V1 = Version(1, 0, 0)

ATMOSPHERE_TEMPERATURE_STATE = "atmosphere_layer_temperature"

ATMOSPHERE_INTERFACE = InterfaceContract(
    interface_id="atmospheric_state",
    version=V1,
    physical_semantics={
        "temperature_units": "K",
        "pressure_units": "Pa",
        "layer_mass_units": "kg m^-2",
        "layer_ordering": "index_0_is_the_lowest_layer",
        "vertical_coordinate": "pressure",
        "state_owner": "atmosphere",
    },
)


#: US Standard Atmosphere 1976 (NOAA/NASA/USAF): layer bases as (base pressure
#: Pa, base temperature K, lapse rate K per geopotential metre), to 51 km.
_US_STANDARD_1976 = (
    (101325.0, 288.15, -0.0065),
    (22632.06, 216.65, 0.0),
    (5474.889, 216.65, 0.001),
    (868.0187, 228.65, 0.0028),
    (110.9063, 270.65, 0.0),
    (66.93887, 270.65, -0.0028),
)


def standard_atmosphere_temperature(pressure: np.ndarray) -> np.ndarray:
    """The US Standard Atmosphere 1976 temperature at each pressure, K.

    An initial state, not a prediction: a column started from it spins up from
    a physical profile rather than from an isothermal one, whose stratosphere
    overshoots its equilibrium while it is still dry
    (``docs/protocols/tropical-column.md``). Within each of the standard's
    layers the temperature follows its constant lapse rate hydrostatically:
    ``T = T_b (p / p_b)^(-R L / g)``, or ``T_b`` where it is isothermal.
    Above 51 km it holds the 51 km value.
    """
    p = np.asarray(pressure, dtype=float)
    r_over_g = value_of("gas_constant_dry_air") / value_of("earth_standard_gravity")
    result = np.empty_like(p)
    for index, value in np.ndenumerate(p):
        base_p, base_t, lapse = _US_STANDARD_1976[0]
        for layer in _US_STANDARD_1976:
            if value <= layer[0]:
                base_p, base_t, lapse = layer
        if value < _US_STANDARD_1976[-1][0]:
            result[index] = _US_STANDARD_1976[-1][1]
            continue
        result[index] = base_t * (value / base_p) ** (-r_over_g * lapse) if lapse else base_t
    return result


def layer_interface_pressures(
    *,
    n_layers: int,
    surface_pressure: float,
    top_pressure: float,
    max_log_step: float | None = None,
) -> np.ndarray:
    """Interface pressures from the surface (index 0) to the model top.

    Without ``max_log_step``: ``n_layers`` layers equal in pressure. With it,
    the step stays ``(p_s - p_top) / n_layers`` wherever that is thinner than
    ``max_log_step`` in log pressure, and is exactly that depth above --
    ``dp_i = min(dp_0, p_i (1 - exp(-max_log_step)))``, built upward -- with
    a remainder under half a step merged into the top layer. So ``n_layers``
    is then the reference count that fixes the lower troposphere's step, and
    the layer count is whatever the cap makes it
    (``docs/protocols/stratospheric-resolution.md``). Below the first capped
    interface the grid is today's to rounding (it is built by repeated steps,
    not by ``linspace``).
    """
    if n_layers < 1:
        raise ConfigurationError("an atmospheric column needs at least one layer")
    if not 0.0 <= top_pressure < surface_pressure:
        raise ConfigurationError(
            f"top pressure {top_pressure} must be non-negative and below the surface "
            f"pressure {surface_pressure}"
        )
    if max_log_step is None:
        return np.linspace(surface_pressure, top_pressure, n_layers + 1)
    if not max_log_step > 0.0:
        raise ConfigurationError(f"the maximum log-pressure step {max_log_step} must be positive")
    if top_pressure <= 0.0:
        raise ConfigurationError("a log-pressure cap needs a model top above zero pressure")
    uniform = (surface_pressure - top_pressure) / n_layers
    fraction = 1.0 - np.exp(-max_log_step)
    interfaces = [float(surface_pressure)]
    while True:
        current = interfaces[-1]
        step = min(uniform, current * fraction)
        following = current - step
        # The last layer takes whatever is left, unless that would leave a
        # sliver under half a step, which joins the layer below it instead.
        if following - top_pressure < 0.5 * min(uniform, following * fraction):
            interfaces.append(float(top_pressure))
            break
        interfaces.append(following)
    return np.asarray(interfaces)


def hydrostatic_layer_masses(
    *,
    n_layers: int,
    surface_pressure: float,
    top_pressure: float,
    gravity: float | GravityField,
    layer_heights: np.ndarray | None = None,
    max_log_step: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Layer interface pressures and masses for a hydrostatic column.

    Returns ``(interface_pressures, layer_masses)`` with interface index 0 at the
    surface and the last at the model top, so layer ``i`` lies between
    interfaces ``i`` and ``i+1``. Layers are equal in pressure unless
    ``max_log_step`` caps their depth aloft (``layer_interface_pressures``);
    either way the total is exactly ``(p_s - p_top) / g`` with constant gravity,
    so the column's mass budget is checkable by hand.

    ``gravity`` may be a number -- which is the first rung of charter 9.15's
    fidelity path, ``ConstantGravity``, and keeps the total exact -- or a
    ``GravityField``. A field that varies with height needs ``layer_heights`` and
    **refuses to be evaluated without them** rather than quietly using the
    surface value, since the point of declaring the approximation is that it
    cannot be applied by accident.
    """
    field = as_gravity_field(gravity)
    interfaces = layer_interface_pressures(
        n_layers=n_layers,
        surface_pressure=surface_pressure,
        top_pressure=top_pressure,
        max_log_step=max_log_step,
    )
    n_layers = interfaces.size - 1
    if field.varies_with_height:
        if layer_heights is None:
            raise ConfigurationError(
                f"{type(field).__name__} varies with height, so the layer heights are "
                "needed to evaluate it; pass layer_heights or use ConstantGravity"
            )
        heights = np.asarray(layer_heights, dtype=float)
        if heights.size != n_layers:
            raise ConfigurationError(
                f"layer_heights has {heights.size} entries for {n_layers} layers"
            )
    else:
        heights = np.zeros(n_layers)
    masses = -np.diff(interfaces) / field.at_height(heights)
    return interfaces, masses


def atmosphere_state_entry(
    *,
    owner: str = "atmosphere",
    n_layers: int,
    tendency_providers: frozenset[str] = frozenset(),
    readers: frozenset[str] = frozenset(),
) -> StateEntry:
    """Registry row for the layer temperatures (STA-002, STA-003)."""
    return StateEntry(
        record_id=ATMOSPHERE_TEMPERATURE_STATE,
        version=V1,
        physical_name="atmospheric layer temperature",
        symbol="T_a",
        owner_module_id=owner,
        units="K",
        state_kind=StateKind.PROGNOSTIC,
        dimensions=("layer",),
        conservation_participation=("energy",),
        authorized_readers=readers,
        authorized_tendency_providers=frozenset({owner}) | tendency_providers,
        update_authority=UpdateAuthority.OWNER_ONLY,
        initialization_source="manifest",
        checkpoint_required=True,
        uncertainty_representation="initial_condition_and_parameter_uncertainty",
        validity_flags=(f"n_layers={n_layers}",),
    )


class AtmosphericColumnModule(Module):
    """A hydrostatic column of isothermal layers, heated only by radiation."""

    def __init__(
        self,
        *,
        module_id: str = "atmosphere",
        surface_module_id: str = "surface",
        space_module_id: str = "space",
        source_module_id: str = "orbit",
        n_layers: int = 8,
        surface_pressure: float | None = None,
        top_pressure: float = 1.0e3,
        initial_temperature: float | np.ndarray = 250.0,
        specific_heat: float | None = None,
        gravity: float | GravityField | None = None,
        couple_water: bool = False,
        water_module_id: str = "water",
        max_log_step: float | None = None,
        status: ModuleStatus = ModuleStatus.ACTIVE,
    ) -> None:
        # With a log-pressure cap, ``n_layers`` is the reference count that sets
        # the lower troposphere's step and the column has however many layers
        # the cap makes; everything below is sized by that count.
        reference_layers = n_layers
        n_layers = (
            layer_interface_pressures(
                n_layers=n_layers,
                surface_pressure=(
                    surface_pressure
                    if surface_pressure is not None
                    else value_of("standard_surface_pressure")
                ),
                top_pressure=top_pressure,
                max_log_step=max_log_step,
            ).size
            - 1
        )
        contract = ModuleContract(
            module_id=module_id,
            module_name=f"hydrostatic radiative column ({n_layers} layers)",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.REDUCED,
            physical_basis=(
                "First law applied to hydrostatically distributed isothermal layers "
                "exchanging energy by radiation only."
            ),
            governing_equation_ids=(
                EQUATION_IDS["hydrostatic"],
                EQUATION_IDS["layer_energy_balance"],
                EQUATION_IDS["forward_euler"],
            ),
            assumptions=(
                "hydrostatic balance with a prescribed, time-invariant pressure structure",
                "each layer is isothermal and well mixed",
                "dry air with constant specific heat",
                "no advection, convection, latent heat or kinetic energy",
            ),
            validity_domain=(
                "Radiative-only column experiments. With no convective adjustment the "
                "column reaches pure radiative equilibrium, which is known to be "
                "convectively unstable in the real lower troposphere; a claim about "
                "lapse rate is outside this module's domain until Stage 3 adds it."
            ),
            owned_state_ids=frozenset({ATMOSPHERE_TEMPERATURE_STATE}),
            required_inputs=(
                "layer_radiative_heating",
                "toa_downward_radiation",
                "toa_upward_radiation",
                "surface_downward_radiation",
                "surface_upward_radiation",
                *(("latent_heating_rate",) if couple_water else ()),
            ),
            produced_outputs=(
                ATMOSPHERE_TEMPERATURE_STATE,
                "atmosphere_layer_mass",
                "atmosphere_layer_heat_capacity",
                "atmosphere_layer_pressure",
                "atmosphere_interface_pressure",
                "atmosphere_column_energy",
            ),
            imported_transfers=(
                "radiation_down_toa",
                "radiation_up_surface",
                # Latent heat arrives as an ordinary transfer, not as a quietly
                # applied tendency: condensation moves energy across a module
                # boundary and CON-002 requires both parties to record it.
                *(("latent_heat_release",) if couple_water else ()),
            ),
            exported_transfers=("radiation_down_surface", "radiation_up_toa"),
            conserved_quantities=("energy", "mass"),
            internal_sources_and_sinks=(),
            provided_interfaces=(ATMOSPHERE_INTERFACE,),
            uncertainty_outputs=(
                "initial_condition",
                "specific_heat_parameter",
                "vertical_resolution",
                "missing_convection_structural",
            ),
            correctness_tests=(
                "test_transparent_column_recovers_the_airless_result",
                "test_single_grey_layer_matches_the_analytic_greenhouse",
                "test_column_energy_budget_closes",
                "test_column_mass_is_hydrostatic",
                "test_zero_insolation_column_only_cools",
            ),
            known_limitations=(
                "No dynamics, convection or water. Pure radiative equilibrium only.",
                "Layer temperatures may become convectively unstable; nothing corrects it.",
            ),
            fidelity_justification=(
                "A single layer has no lapse rate, and the lapse rate is the quantity "
                "convection acts on and the one the entraining-plume closure predicts. "
                "It would also put the whole column's emission at one temperature, "
                "which removes the height dependence the greenhouse effect is made of."
            ),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(n_layers)",
                memory_complexity="O(n_layers)",
                scaling_variables=("n_layers",),
                work_units_per_step=1.0,
                benchmark_configuration=f"{n_layers} layers, global scalar column",
                peak_memory_bytes=8 * n_layers + 4096,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
                lower_fidelity_alternatives=("single_layer_column",),
                checkpoint_cost_bytes=8 * n_layers,
            ),
        )
        super().__init__(contract, status=status)
        self.couple_water = couple_water
        self.water_module_id = water_module_id
        self.surface_module_id = surface_module_id
        self.space_module_id = space_module_id
        self.source_module_id = source_module_id
        self.n_layers = n_layers
        self.reference_layers = reference_layers
        self.max_log_step = max_log_step
        self.surface_pressure = (
            surface_pressure
            if surface_pressure is not None
            else value_of("standard_surface_pressure")
        )
        self.top_pressure = top_pressure
        self.specific_heat = (
            specific_heat if specific_heat is not None else value_of("specific_heat_dry_air")
        )
        self.gravity_field = as_gravity_field(gravity)
        self.initial_temperature = initial_temperature

        # The pressure grid does not depend on gravity, so the interfaces come
        # first and the masses follow, which is what lets a height-dependent
        # field be evaluated at all.
        interfaces = layer_interface_pressures(
            n_layers=reference_layers,
            surface_pressure=self.surface_pressure,
            top_pressure=top_pressure,
            max_log_step=max_log_step,
        )
        # Layer midpoint pressure, the level a layer's bulk properties are
        # attributed to: the arithmetic mean of its interfaces, which is its
        # mass-weighted mean pressure whatever its depth.
        self.layer_pressure = 0.5 * (interfaces[:-1] + interfaces[1:])
        self.layer_height = self._initial_heights()
        self.interface_pressure, self.layer_mass = hydrostatic_layer_masses(
            n_layers=reference_layers,
            surface_pressure=self.surface_pressure,
            top_pressure=top_pressure,
            gravity=self.gravity_field,
            layer_heights=self.layer_height if self.gravity_field.varies_with_height else None,
            max_log_step=max_log_step,
        )
        self.layer_heat_capacity = self.specific_heat * self.layer_mass
        self.temperature = np.full(n_layers, float(np.mean(np.asarray(initial_temperature))))

        self._interval = (0.0, 0.0)
        self._stored_initial = 0.0
        self._budget: StateBudget | None = None
        self._boundary_energy: dict[str, float] = {}
        self._latent_received = 0.0
        self._sensible_received: float | None = None
        self._convective_received: float | None = None

    def _initial_heights(self) -> np.ndarray:
        """Layer midpoint heights, from the initial profile, in metres.

        Only a height-dependent gravity field needs these, and it needs them
        before the column has taken a step, so they come from the initial
        temperature. That is an approximation with a name: the layer masses are
        fixed for the run either way -- the pressure grid is time-invariant by
        declaration -- so what is fixed here is which height each *pressure* is
        at, and that moves by a few percent over the run as the column warms. A
        few percent of a 0.2 percent correction is not a term anything here
        resolves, and declaring the mass field time-invariant is the assumption
        that was already being made.
        """
        initial = np.asarray(self.initial_temperature, dtype=float)
        profile = np.full(self.n_layers, float(np.mean(initial))) if initial.ndim == 0 else initial
        if profile.size != self.n_layers:
            # Left to validate_configuration to report properly; a placeholder
            # keeps construction from failing before it can.
            profile = np.full(self.n_layers, float(np.mean(profile)))
        return hydrostatic_heights(profile, self.layer_pressure, self.surface_pressure)

    @property
    def gravity(self) -> float:
        """The surface value, for callers that want one number (GRV-001)."""
        return self.gravity_field.surface_value()

    # --- configuration ----------------------------------------------------
    def validate_configuration(self) -> None:
        super().validate_configuration()
        if self.specific_heat <= 0.0:
            raise ConfigurationError(f"specific heat {self.specific_heat} must be positive")
        if self.gravity <= 0.0:
            raise ConfigurationError(f"gravity {self.gravity} must be positive")
        initial = np.asarray(self.initial_temperature, dtype=float)
        if initial.ndim == 1 and initial.size != self.n_layers:
            raise ConfigurationError(
                f"initial temperature has {initial.size} entries for {self.n_layers} layers"
            )
        if np.any(initial <= 0.0):
            raise ConfigurationError("initial layer temperatures must be positive")

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        supplied = initial_state.get(ATMOSPHERE_TEMPERATURE_STATE, self.initial_temperature)
        values = np.asarray(supplied, dtype=float)
        self.temperature = (
            np.full(self.n_layers, float(values))
            if values.ndim == 0
            else np.array(values, dtype=float)
        )
        if self.temperature.size != self.n_layers:
            raise ConfigurationError(
                f"initial atmospheric state has {self.temperature.size} layers, "
                f"expected {self.n_layers}"
            )
        self.mark_initialized()

    # --- physics ----------------------------------------------------------
    @property
    def column_mass(self) -> float:
        """Total dry mass per unit area, exactly ``(p_s - p_top) / g``."""
        return float(np.sum(self.layer_mass))

    def stored_energy(self, temperature: np.ndarray | None = None) -> float:
        """Column internal energy per unit area, ``sum_i c_p m_i T_i`` in J m^-2."""
        T = self.temperature if temperature is None else temperature
        return float(np.sum(self.layer_heat_capacity * T))

    def compute_tendencies(self) -> list[TendencyContribution]:
        """The column contributes nothing of its own; radiation supplies the heating.

        Returning an empty list is the honest statement of this fidelity level: at
        Stage 2 the only process acting on layer temperature is radiative flux
        divergence, and the radiation module owns that computation and labels it.
        """
        ctx = self.context
        self._interval = (ctx.time, ctx.time + ctx.dt)
        self._stored_initial = self.stored_energy()
        self._boundary_energy = {
            channel: self.input_field(field_id).require_units("W m^-2").scalar() * ctx.dt
            for channel, field_id in (
                ("radiation_down_toa", "toa_downward_radiation"),
                ("radiation_up_toa", "toa_upward_radiation"),
                ("radiation_down_surface", "surface_downward_radiation"),
                ("radiation_up_surface", "surface_upward_radiation"),
            )
        }
        sensible = self._inputs.get("surface_sensible_heat_flux")
        self._sensible_received = (
            sensible.require_units("W m^-2").scalar() * ctx.dt if sensible is not None else None
        )
        # Present only when convection mixes the surface into the column. The
        # energy is real and crosses a module boundary, so it is declared the
        # same way the sensible flux is rather than folded into the mixing.
        convective = self._inputs.get("surface_convective_heat_flux")
        self._convective_received = (
            convective.require_units("W m^-2").scalar() * ctx.dt if convective is not None else None
        )
        self._latent_received = (
            float(np.sum(self.input_field("latent_heating_rate").require_units("W m^-2").values))
            * ctx.dt
            if self.couple_water
            else 0.0
        )
        return []

    def advance_state(self, dt: float, contributions: list[TendencyContribution]) -> None:
        initial = np.array(self.temperature, dtype=float)
        total = np.zeros(self.n_layers)
        for contribution in contributions:
            total = total + np.broadcast_to(contribution.as_array(), (self.n_layers,))
        self.temperature = initial + total * dt
        # Charter 15.2: refuse an impossible state rather than integrating on
        # from it. Forward Euler beyond its stability limit produces negative
        # temperatures, and a column that keeps going produces numbers whose
        # budget still closes -- the accounting is consistent about nonsense --
        # so nothing downstream would have caught it.
        if np.any(self.temperature <= 0.0):
            worst = float(np.min(self.temperature))
            raise ConfigurationError(
                f"layer temperature reached {worst:.4g} K, which is unphysical. Forward "
                f"Euler is unstable when the step approaches a layer's radiative "
                f"relaxation time; reduce the master step or the opacity, and note that "
                f"conservation closing is not evidence the answer is physical"
            )
        self._budget = StateBudget(
            state_id=ATMOSPHERE_TEMPERATURE_STATE,
            initial_value=initial,
            final_value=np.array(self.temperature, dtype=float),
            dt=dt,
            contributions=list(contributions),
        )

    def state_budgets(self) -> dict[str, Any]:
        return {ATMOSPHERE_TEMPERATURE_STATE: self._budget} if self._budget else {}

    # --- outputs ----------------------------------------------------------
    def emit_state(self) -> dict[str, Field]:
        provenance = Provenance(producer=self.module_id, derivation=EQUATION_IDS["hydrostatic"])
        timestamp = self._context.time if self._context else 0.0
        return {
            ATMOSPHERE_TEMPERATURE_STATE: Field(
                field_id=ATMOSPHERE_TEMPERATURE_STATE,
                physical_name="atmospheric layer temperature",
                symbol="T_a",
                values=np.array(self.temperature),
                units="K",
                quantity_kind=QuantityKind.INTENSIVE_STATE,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                timestamp=timestamp,
                spatial_support=SpatialSupport.LAYER_AVERAGE,
                dimensions=("layer",),
                coordinate_system="pressure_layers",
                provenance=provenance,
                uncertainty=Uncertainty(
                    method="none",
                    unresolved=("initial_condition", "vertical_resolution"),
                ),
            ),
            "atmosphere_layer_mass": Field(
                field_id="atmosphere_layer_mass",
                physical_name="dry air mass per unit area in each layer",
                symbol="m_i",
                values=np.array(self.layer_mass),
                units="kg m^-2",
                quantity_kind=QuantityKind.EXTENSIVE_STATE,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                timestamp=timestamp,
                spatial_support=SpatialSupport.LAYER_AVERAGE,
                dimensions=("layer",),
                coordinate_system="pressure_layers",
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                provenance=provenance,
            ),
            "atmosphere_layer_pressure": Field(
                field_id="atmosphere_layer_pressure",
                physical_name="midpoint pressure of each layer",
                symbol="p_i",
                values=np.array(self.layer_pressure),
                units="Pa",
                quantity_kind=QuantityKind.INTENSIVE_STATE,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                timestamp=timestamp,
                spatial_support=SpatialSupport.LAYER_AVERAGE,
                dimensions=("layer",),
                coordinate_system="pressure_layers",
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                provenance=provenance,
            ),
            "atmosphere_layer_heat_capacity": Field(
                field_id="atmosphere_layer_heat_capacity",
                physical_name="areal heat capacity of each layer",
                symbol="C_i",
                values=np.array(self.layer_heat_capacity),
                units="J m^-2 K^-1",
                quantity_kind=QuantityKind.EXTENSIVE_STATE,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                timestamp=timestamp,
                spatial_support=SpatialSupport.LAYER_AVERAGE,
                dimensions=("layer",),
                coordinate_system="pressure_layers",
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                provenance=provenance,
            ),
            "atmosphere_interface_pressure": Field(
                field_id="atmosphere_interface_pressure",
                physical_name="pressure at layer interfaces, surface first",
                symbol="p",
                values=np.array(self.interface_pressure),
                units="Pa",
                quantity_kind=QuantityKind.INTENSIVE_STATE,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                timestamp=timestamp,
                spatial_support=SpatialSupport.FACE_CENTERED,
                dimensions=("interface",),
                coordinate_system="pressure_layers",
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                provenance=provenance,
            ),
        }

    def emit_outputs(self) -> dict[str, Field]:
        ctx = self.context
        return {
            "atmosphere_column_energy": Field(
                field_id="atmosphere_column_energy",
                physical_name="column internal energy per unit area",
                symbol="E_a",
                values=np.array(self.stored_energy()),
                units="J m^-2",
                quantity_kind=QuantityKind.EXTENSIVE_STATE,
                time_semantics=TimeSemantics.START_OF_STEP,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.DOMAIN_INTEGRAL,
                provenance=Provenance(
                    producer=self.module_id, derivation=EQUATION_IDS["layer_energy_balance"]
                ),
            )
        }

    def report_transfers(self) -> tuple[list[Transfer], list[Transfer]]:
        """The column's own record of what crossed its two interfaces.

        Computed from the flux fields the radiation module published, not from
        the radiation module's internals, so the reconciliation against the
        surface's and the boundary sink's records is an independent check.
        """
        start, end = self._interval
        common = {
            "quantity": ConservedQuantity.ENERGY,
            "units": "J m^-2",
            "start_time": start,
            "end_time": end,
            "spatial_mapping": "global_mean_over_planetary_surface",
        }
        imports = [
            Transfer(
                transfer_id=f"radiation_down_toa@{start:g}",
                sender_module=self.source_module_id,
                receiver_module=self.module_id,
                amount=self._boundary_energy["radiation_down_toa"],
                channel="radiation_down_toa",
                **common,  # type: ignore[arg-type]
            ),
            Transfer(
                transfer_id=f"radiation_up_surface@{start:g}",
                sender_module=self.surface_module_id,
                receiver_module=self.module_id,
                amount=self._boundary_energy["radiation_up_surface"],
                channel="radiation_up_surface",
                **common,  # type: ignore[arg-type]
            ),
        ]
        exports = [
            Transfer(
                transfer_id=f"radiation_down_surface@{start:g}",
                sender_module=self.module_id,
                receiver_module=self.surface_module_id,
                amount=self._boundary_energy["radiation_down_surface"],
                channel="radiation_down_surface",
                **common,  # type: ignore[arg-type]
            ),
            Transfer(
                transfer_id=f"radiation_up_toa@{start:g}",
                sender_module=self.module_id,
                receiver_module=self.space_module_id,
                amount=self._boundary_energy["radiation_up_toa"],
                channel="radiation_up_toa",
                **common,  # type: ignore[arg-type]
            ),
        ]
        if self._sensible_received is not None:
            imports.append(
                Transfer(
                    transfer_id=f"sensible_heat@{start:g}",
                    sender_module=self.surface_module_id,
                    receiver_module=self.module_id,
                    amount=self._sensible_received,
                    channel="sensible_heat",
                    **common,  # type: ignore[arg-type]
                )
            )
        if self._convective_received is not None:
            imports.append(
                Transfer(
                    transfer_id=f"surface_convective_heat@{start:g}",
                    sender_module=self.surface_module_id,
                    receiver_module=self.module_id,
                    amount=self._convective_received,
                    channel="surface_convective_heat",
                    **common,  # type: ignore[arg-type]
                )
            )
        if self.couple_water:
            imports.append(
                Transfer(
                    transfer_id=f"latent_heat_release@{start:g}",
                    sender_module=self.water_module_id,
                    receiver_module=self.module_id,
                    amount=self._latent_received,
                    channel="latent_heat_release",
                    **common,  # type: ignore[arg-type]
                )
            )
        return (exports, imports)

    def report_conservation(self) -> list[ModuleBudget]:
        return [
            ModuleBudget(
                module_id=self.module_id,
                quantity=ConservedQuantity.ENERGY,
                units="J m^-2",
                stored_initial=self._stored_initial,
                stored_final=self.stored_energy(),
            ),
            # Mass is conserved trivially at this fidelity -- the column neither
            # gains nor loses air -- but the budget is published rather than
            # assumed, so a later module that does move mass has a line to use.
            ModuleBudget(
                module_id=self.module_id,
                quantity=ConservedQuantity.MASS,
                units="kg m^-2",
                stored_initial=self.column_mass,
                stored_final=self.column_mass,
            ),
        ]

    def report_uncertainty(self) -> UncertaintyReport:
        return UncertaintyReport(
            module_id=self.module_id,
            quantified={"column_mass_kg_m2": self.column_mass},
            method="none propagated at this fidelity",
            distributional_assumptions="unresolved",
            correlations="unresolved",
            unresolved=(
                "initial_condition_uncertainty",
                "specific_heat_parameter_uncertainty",
                "vertical_resolution_uncertainty",
                "structural_uncertainty_of_omitting_convection_and_dynamics",
            ),
        )

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "n_layers": self.n_layers,
            "surface_pressure_Pa": self.surface_pressure,
            "top_pressure_Pa": self.top_pressure,
            # Only when capped, so today's grid hashes as it always has.
            **(
                {"max_log_step": self.max_log_step, "reference_layers": self.reference_layers}
                if self.max_log_step is not None
                else {}
            ),
            "specific_heat_J_kg_K": self.specific_heat,
            **self.gravity_field.parameters(),
            "column_mass_kg_m2": self.column_mass,
            "couple_water": self.couple_water,
            "initial_temperature_K": np.asarray(self.initial_temperature, dtype=float)
            .reshape(-1)
            .tolist(),
            "time_integration": "forward_euler",
        }

    def checkpoint(self) -> dict[str, Any]:
        return {ATMOSPHERE_TEMPERATURE_STATE: np.array(self.temperature)}

    def restore(self, checkpoint: dict[str, Any]) -> None:
        self.temperature = np.array(checkpoint[ATMOSPHERE_TEMPERATURE_STATE], dtype=float)
