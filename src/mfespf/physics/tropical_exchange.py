"""A tropical column's exchange with the rest of the planet.

Why this exists
---------------

A global-mean column's stratosphere is saturated at its own cold point: 205.7 K
at 100 hPa gives 37.6 ppmv of water, against the observed 4, and that water
cools its middle stratosphere by 20 to 30 K (``docs/protocols/model-top.md``).
The observed 4 ppmv is made where air actually enters the stratosphere --
rising through the cold tropical tropopause in the Brewer-Dobson circulation.
A column that stands for the tropics has to exchange two things with the rest
of the planet that a global column does not
(``docs/protocols/tropical-column.md``):

- **its heat export.** The tropics absorb more sunlight than they emit and
  carry the surplus poleward. Here it is a declared lateral boundary flux,
  taken from the atmosphere below 100 hPa in proportion to layer mass.
- **its stratospheric upwelling.** Air rises through the tropical tropopause
  at a measured mass flux (Seviour et al. 2012) and leaves poleward above. It
  is carried in flux form, upwind: dry static energy ``c_p T + g z`` and
  vapour move upward, and air enters and leaves laterally with the local
  properties. So each layer changes at ``F_k (phi_(k-1) - phi_k) / m_k``, and
  the cooling of rising air is the transport, not an added term.

Whatever leaves the column goes to ``LateralBoundaryModule``, which records the
receiving side from the fields this module publishes. The budgets close as
they do at the top of the atmosphere.
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
)
from mfespf.core.versioning import InterfaceContract, Version
from mfespf.coupling.conservation import ModuleBudget
from mfespf.coupling.transfers import ConservedQuantity, Transfer
from mfespf.physics.equations import EQUATION_IDS
from mfespf.physics.thermodynamics import hydrostatic_heights
from mfespf.registries.states import (
    StateBudget,
    StateEntry,
    StateKind,
    TendencyContribution,
    UpdateAuthority,
)

__all__ = [
    "EXPORT_TOP_PA",
    "LATERAL_ENERGY_STATE",
    "LATERAL_WATER_STATE",
    "OBSERVED_UPWELLING_MASS_FLUX_KG_S",
    "UPWELLING_BASE_PA",
    "UPWELLING_REFERENCE_PA",
    "LateralBoundaryModule",
    "TropicalExchangeModule",
    "export_heating",
    "lateral_boundary_state_entries",
    "upward_advection_rate",
    "upwelling_mass_flux",
]

V1 = Version(1, 0, 0)

#: Annual-mean upward mass flux at 70 hPa between the turnaround latitudes,
#: ERA-Interim 1989-2009 (Seviour, Butchart and Hardiman 2012, QJRMS); 6.2e9
#: in UKMO analyses and 5.9e9 across chemistry-climate models.
OBSERVED_UPWELLING_MASS_FLUX_KG_S = 5.9e9
#: The pressure that flux is measured at.
UPWELLING_REFERENCE_PA = 7000.0
#: Below this the upwelling is zero: the troposphere's ascent is convection's.
UPWELLING_BASE_PA = 10000.0
#: The heat export is taken from the atmosphere below this pressure.
EXPORT_TOP_PA = 10000.0

LATERAL_ENERGY_STATE = "lateral_boundary_cumulative_energy"
LATERAL_WATER_STATE = "lateral_boundary_cumulative_water"

_TEMPERATURE = "atmosphere_layer_temperature"
_VAPOUR = "atmosphere_water_vapour"

EXCHANGE_INTERFACE = InterfaceContract(
    interface_id="lateral_exchange",
    version=V1,
    physical_semantics={
        "upwelling_sign": "positive_upward",
        "upwelling_basis": "per_unit_area_of_the_column",
        "transport": "flux_form_upwind_dry_static_energy_and_vapour",
        "export_placement": "atmosphere_below_export_top_by_layer_mass",
    },
)


def upwelling_mass_flux(
    interface_pressure: np.ndarray,
    *,
    flux_at_reference: float,
    reference_pressure: float = UPWELLING_REFERENCE_PA,
    base_pressure: float = UPWELLING_BASE_PA,
) -> np.ndarray:
    """Upward mass flux at each interface, kg m^-2 s^-1, surface (0) to top.

    Zero below ``base_pressure``; ``flux_at_reference`` from there up to
    ``reference_pressure``; falling linearly in pressure to zero at the top
    interface above that. The shape is declared, not measured
    (``docs/protocols/tropical-column.md``).
    """
    p = np.asarray(interface_pressure, dtype=float)
    if flux_at_reference < 0.0:
        raise ConfigurationError("the upwelling mass flux cannot be negative")
    top = float(p[-1])
    if not top < reference_pressure < base_pressure:
        raise ConfigurationError(
            "the upwelling needs the model top above its reference pressure, and the "
            "reference above its base"
        )
    flux = np.where(p <= base_pressure, flux_at_reference, 0.0)
    above = p < reference_pressure
    flux[above] = flux_at_reference * (p[above] - top) / (reference_pressure - top)
    flux[0] = 0.0
    flux[-1] = 0.0
    return flux


def upward_advection_rate(flux: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Rate of change of ``phi m`` per layer under an upward mass flux.

    Upwind: interface ``k`` carries the value of layer ``k-1`` below it, and the
    lateral exchange ``F_(k+1) - F_k`` enters or leaves with the layer's own
    value. So each layer changes at ``F_k (phi_(k-1) - phi_k)``, and the column
    total changes only by what the lateral exchange carries.
    """
    f = np.asarray(flux, dtype=float)
    q = np.asarray(phi, dtype=float)
    if f.size != q.size + 1:
        raise ConfigurationError("the fluxes must sit on the layers' interfaces")
    rate = np.zeros_like(q)
    rate[1:] = f[1:-1] * (q[:-1] - q[1:])
    return rate


def export_heating(
    layer_mass: np.ndarray, layer_pressure: np.ndarray, export: float, top: float = EXPORT_TOP_PA
) -> np.ndarray:
    """The heat export as a heating per layer, W m^-2 (negative), mass-weighted."""
    m = np.asarray(layer_mass, dtype=float)
    below = np.asarray(layer_pressure, dtype=float) > top
    if not np.any(below):
        raise ConfigurationError("no layer lies below the export's top")
    weight = np.where(below, m, 0.0)
    return -float(export) * weight / float(np.sum(weight))


class TropicalExchangeModule(Module):
    """Heat export and Brewer-Dobson upwelling for a column standing for the tropics."""

    def __init__(
        self,
        *,
        module_id: str = "tropical_exchange",
        atmosphere_module_id: str = "atmosphere",
        water_module_id: str = "water",
        boundary_module_id: str = "lateral",
        upwelling_flux: float,
        heat_export: float,
        export_from_surface_fraction: float = 0.0,
        status: ModuleStatus = ModuleStatus.ACTIVE,
    ) -> None:
        contract = ModuleContract(
            module_id=module_id,
            module_name="tropical lateral exchange: heat export and Brewer-Dobson upwelling",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.REDUCED,
            physical_basis=(
                "The tropics export their top-of-atmosphere surplus poleward, and air "
                "rises through the tropical tropopause in the Brewer-Dobson circulation "
                "and leaves poleward above it. Both are exchanges with the rest of the "
                "planet that a column standing for the tropics has to make."
            ),
            governing_equation_ids=(
                EQUATION_IDS["brewer_dobson_upwelling"],
                EQUATION_IDS["lateral_heat_export"],
            ),
            assumptions=(
                "the upwelling mass flux is prescribed from observation, with a declared "
                "vertical shape",
                "air enters and leaves the column laterally with the local properties",
                "the heat export is prescribed from observation and taken from the "
                "atmosphere below 100 hPa in proportion to layer mass",
            ),
            validity_domain=(
                "Annual-mean experiments with a column standing for the tropical band "
                "between the turnaround latitudes. Nothing here responds to the column: "
                "the export and the upwelling are boundary conditions."
            ),
            required_inputs=(
                _TEMPERATURE,
                _VAPOUR,
                "atmosphere_layer_mass",
                "atmosphere_layer_pressure",
                "atmosphere_interface_pressure",
                "atmosphere_layer_heat_capacity",
            ),
            produced_outputs=(
                "vapour_advection_rate",
                "advective_heating_rate",
                "lateral_energy_export",
                "lateral_water_export",
                "lateral_latent_energy_export",
                "brewer_dobson_mass_flux",
            ),
            exported_transfers=(
                "lateral_energy",
                "lateral_water",
                "lateral_latent_energy",
            ),
            conserved_quantities=("energy", "water"),
            internal_sources_and_sinks=(),
            provided_interfaces=(EXCHANGE_INTERFACE,),
            uncertainty_outputs=(
                "upwelling_mass_flux_observational",
                "upwelling_shape_structural",
                "heat_export_observational",
                "export_placement_structural",
            ),
            correctness_tests=(
                "test_the_upwelling_conserves_mass_layer_by_layer",
                "test_rising_air_cools_at_the_dry_adiabatic_rate_in_an_isothermal_layer",
                "test_the_exchange_budget_closes_through_the_lateral_boundary",
                "test_no_exchange_leaves_the_column_unchanged",
            ),
            known_limitations=(
                "The export and the upwelling do not respond to the column.",
                "The export is all atmospheric by default; part of the real transport is oceanic.",
            ),
            fidelity_justification=(
                "A global-mean column cannot represent where air enters the "
                "stratosphere, and its stratospheric water is saturated at a cold point "
                "5 to 15 K warmer than the tropical one."
            ),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(n_layers)",
                memory_complexity="O(n_layers)",
                scaling_variables=("n_layers",),
                work_units_per_step=0.2,
                benchmark_configuration="one 51-layer column",
                peak_memory_bytes=65536,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
                lower_fidelity_alternatives=("global_mean_column",),
            ),
        )
        super().__init__(contract, status=status)
        self.atmosphere_module_id = atmosphere_module_id
        self.water_module_id = water_module_id
        self.boundary_module_id = boundary_module_id
        self.upwelling_flux = float(upwelling_flux)
        self.heat_export = float(heat_export)
        self.export_from_surface_fraction = float(export_from_surface_fraction)
        self.latent_heat = value_of("water_latent_heat_vaporization_273")
        self.gravity = value_of("earth_standard_gravity")
        self._flux = np.zeros(0)
        self._rates: dict[str, np.ndarray] = {}
        self._energy_out = 0.0
        self._water_out = 0.0
        self._transfers: list[Transfer] = []

    def validate_configuration(self) -> None:
        super().validate_configuration()
        if self.upwelling_flux < 0.0:
            raise ConfigurationError("the upwelling mass flux cannot be negative")
        if not 0.0 <= self.export_from_surface_fraction <= 1.0:
            raise ConfigurationError("the surface's share of the export must lie in [0, 1]")
        if self.export_from_surface_fraction > 0.0:
            raise ConfigurationError(
                "an export taken from the surface needs the surface module to accept it, "
                "which it does not yet"
            )

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        self.mark_initialized()

    def compute_tendencies(self) -> list[TendencyContribution]:
        ctx = self.context

        def get(name: str, units: str) -> np.ndarray:
            return np.asarray(self.input_field(name).require_units(units).values, dtype=float)

        T = get(_TEMPERATURE, "K")
        w = get(_VAPOUR, "kg m^-2")
        m = get("atmosphere_layer_mass", "kg m^-2")
        p = get("atmosphere_layer_pressure", "Pa")
        p_i = get("atmosphere_interface_pressure", "Pa")
        C = get("atmosphere_layer_heat_capacity", "J m^-2 K^-1")

        # Gravity as the grid has it, the same way the dry static energy was
        # built for the circulation; heights from the column's own profile.
        z = hydrostatic_heights(T, p, float(p_i[0]))
        s = C / m * T + self.gravity * z
        q = w / m

        self._flux = upwelling_mass_flux(p_i, flux_at_reference=self.upwelling_flux)
        courant = self._flux[1:-1] * ctx.dt / np.minimum(m[:-1], m[1:])
        if np.any(courant > 1.0):
            raise ConfigurationError(
                f"the upwelling moves {float(courant.max()):.3g} of a layer per step; "
                "shorten the step"
            )
        energy_rate = upward_advection_rate(self._flux, s)
        vapour_rate = upward_advection_rate(self._flux, q)
        export = export_heating(m, p, self.heat_export)
        heating = energy_rate + export

        self._rates = {"vapour_advection_rate": vapour_rate, "advective_heating_rate": heating}
        # Out of the column over the step, per unit of its area: positive leaves.
        self._energy_out = -float(np.sum(heating)) * ctx.dt
        self._water_out = -float(np.sum(vapour_rate)) * ctx.dt
        self._transfers = [
            self._transfer(
                "energy",
                self.atmosphere_module_id,
                ConservedQuantity.ENERGY,
                "J m^-2",
                self._energy_out,
            ),
            self._transfer(
                "water", self.water_module_id, ConservedQuantity.WATER, "kg m^-2", self._water_out
            ),
            self._transfer(
                "latent_energy",
                self.water_module_id,
                ConservedQuantity.ENERGY,
                "J m^-2",
                self._water_out * self.latent_heat,
            ),
        ]
        return [
            TendencyContribution(
                state_id=_TEMPERATURE,
                provider_module_id=self.module_id,
                label="tropical_exchange",
                value=heating / C,
                units="K s^-1",
                equation_id=EQUATION_IDS["brewer_dobson_upwelling"],
            ),
            TendencyContribution(
                state_id=_VAPOUR,
                provider_module_id=self.module_id,
                label="brewer_dobson_upwelling",
                value=vapour_rate,
                units="kg m^-2 s^-1",
                equation_id=EQUATION_IDS["brewer_dobson_upwelling"],
            ),
        ]

    def _transfer(
        self, name: str, sender: str, quantity: ConservedQuantity, units: str, amount: float
    ) -> Transfer:
        """One exchange, from the column's owner to the boundary, or back if negative."""
        ctx = self.context
        outward = amount >= 0.0
        return Transfer(
            transfer_id=f"lateral_{name}@{ctx.time:g}",
            sender_module=sender if outward else self.boundary_module_id,
            receiver_module=self.boundary_module_id if outward else sender,
            quantity=quantity,
            amount=abs(amount),
            units=units,
            start_time=ctx.time,
            end_time=ctx.time + ctx.dt,
            channel=f"lateral_{name}",
        )

    def report_transfers(self) -> tuple[list[Transfer], list[Transfer]]:
        # The column's side only; the boundary records its own from the
        # published fields, so the reconciliation compares two computations.
        return list(self._transfers), []

    def report_conservation(self) -> list[ModuleBudget]:
        return []

    def emit_outputs(self) -> dict[str, Field]:
        ctx = self.context
        provenance = Provenance(
            producer=self.module_id, derivation=EQUATION_IDS["brewer_dobson_upwelling"]
        )
        common: dict[str, Any] = {
            "time_semantics": TimeSemantics.INTERVAL_AVERAGE,
            "time_interval": ctx.dt,
            "timestamp": ctx.time,
            "provenance": provenance,
        }
        layer: dict[str, Any] = {
            "spatial_support": SpatialSupport.LAYER_AVERAGE,
            "dimensions": ("layer",),
        }
        scalar: dict[str, Any] = {
            "spatial_support": SpatialSupport.GLOBAL_SCALAR,
            # Into the lateral boundary, which is what receives it.
            "sign_convention": SignConvention.POSITIVE_INTO_RECEIVER,
            "quantity_kind": QuantityKind.FLUX,
        }
        return {
            "vapour_advection_rate": Field(
                field_id="vapour_advection_rate",
                physical_name="vapour tendency from the lateral exchange",
                values=np.array(self._rates["vapour_advection_rate"]),
                units="kg m^-2 s^-1",
                quantity_kind=QuantityKind.TENDENCY,
                **layer,
                **common,
            ),
            "advective_heating_rate": Field(
                field_id="advective_heating_rate",
                physical_name="heating from the lateral exchange",
                values=np.array(self._rates["advective_heating_rate"]),
                units="W m^-2",
                quantity_kind=QuantityKind.TENDENCY,
                **layer,
                **common,
            ),
            "brewer_dobson_mass_flux": Field(
                field_id="brewer_dobson_mass_flux",
                physical_name="upward mass flux at each interface",
                values=np.array(self._flux),
                units="kg m^-2 s^-1",
                quantity_kind=QuantityKind.FLUX,
                spatial_support=SpatialSupport.LAYER_AVERAGE,
                dimensions=("interface",),
                **common,
            ),
            "lateral_energy_export": Field(
                field_id="lateral_energy_export",
                physical_name="dry static energy leaving the column laterally",
                values=np.array(self._energy_out / ctx.dt),
                units="W m^-2",
                **scalar,
                **common,
            ),
            "lateral_water_export": Field(
                field_id="lateral_water_export",
                physical_name="water leaving the column laterally",
                values=np.array(self._water_out / ctx.dt),
                units="kg m^-2 s^-1",
                **scalar,
                **common,
            ),
            "lateral_latent_energy_export": Field(
                field_id="lateral_latent_energy_export",
                physical_name="latent energy leaving the column with its water",
                values=np.array(self._water_out * self.latent_heat / ctx.dt),
                units="W m^-2",
                **scalar,
                **common,
            ),
        }

    def report_uncertainty(self) -> UncertaintyReport:
        return UncertaintyReport(
            module_id=self.module_id,
            quantified={
                "upwelling_flux_kg_m2_s": self.upwelling_flux,
                "heat_export_W_m2": self.heat_export,
            },
            method="declared; the latitude range 20-30 degrees is swept by the protocol",
            unresolved=(
                "the upwelling's vertical shape is declared, not measured",
                "the heat export's placement in the column",
            ),
        )

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "upwelling_flux_kg_m2_s": self.upwelling_flux,
            "upwelling_source": "Seviour, Butchart and Hardiman 2012, 5.9e9 kg s^-1 at 70 hPa",
            "heat_export_W_m2": self.heat_export,
            "heat_export_source": (
                "NCEP/NCAR R1 1991-2020 band-mean net TOA minus its global mean"
            ),
            "export_from_surface_fraction": self.export_from_surface_fraction,
        }


def lateral_boundary_state_entries(*, owner: str = "lateral") -> list[StateEntry]:
    """Registry rows for what the lateral boundary has received."""

    def entry(record_id: str, name: str, units: str, participation: str) -> StateEntry:
        return StateEntry(
            record_id=record_id,
            version=V1,
            physical_name=name,
            symbol=record_id,
            owner_module_id=owner,
            units=units,
            state_kind=StateKind.PROGNOSTIC,
            conservation_participation=(participation,),
            authorized_tendency_providers=frozenset({owner}),
            update_authority=UpdateAuthority.OWNER_ONLY,
            initialization_source="zero_at_run_start",
            checkpoint_required=True,
            uncertainty_representation="inherits_from_contributing_fluxes",
        )

    return [
        entry(LATERAL_ENERGY_STATE, "cumulative energy exported laterally", "J m^-2", "energy"),
        entry(LATERAL_WATER_STATE, "cumulative water exported laterally", "kg m^-2", "water"),
    ]


class LateralBoundaryModule(Module):
    """Accumulates what a column exports to the rest of the planet.

    The lateral counterpart of the space boundary: it records the receiving side
    of each exchange from the fields the exchange module publishes, so the
    reconciliation compares two computations rather than one restated.
    """

    def __init__(
        self,
        *,
        module_id: str = "lateral",
        atmosphere_module_id: str = "atmosphere",
        water_module_id: str = "water",
        status: ModuleStatus = ModuleStatus.ACTIVE,
    ) -> None:
        contract = ModuleContract(
            module_id=module_id,
            module_name="lateral boundary energy and water sink",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.NULL_CONTROL,
            physical_basis=(
                "Accounting boundary. What a regional column exports to the rest of the "
                "planet is accumulated so its budgets can be closed and reported."
            ),
            governing_equation_ids=("EQ-BOUNDARY-ACCUMULATION",),
            assumptions=("the rest of the planet is not represented; it only receives",),
            validity_domain=(
                "Bookkeeping only. It makes no physical prediction and cannot support a "
                "claim about anything other than closure of the modeled budget."
            ),
            owned_state_ids=frozenset({LATERAL_ENERGY_STATE, LATERAL_WATER_STATE}),
            required_inputs=(
                "lateral_energy_export",
                "lateral_water_export",
                "lateral_latent_energy_export",
            ),
            produced_outputs=(LATERAL_ENERGY_STATE, LATERAL_WATER_STATE),
            imported_transfers=("lateral_energy", "lateral_water", "lateral_latent_energy"),
            conserved_quantities=("energy", "water"),
            uncertainty_outputs=("inherited_from_contributing_fluxes",),
            correctness_tests=("test_the_exchange_budget_closes_through_the_lateral_boundary",),
            known_limitations=("Represents no physical process; it is an accounting sink.",),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(1)",
                memory_complexity="O(1)",
                work_units_per_step=0.1,
                benchmark_configuration="global scalar",
                peak_memory_bytes=1024,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
                checkpoint_cost_bytes=16,
            ),
        )
        super().__init__(contract, status=status)
        self.atmosphere_module_id = atmosphere_module_id
        self.water_module_id = water_module_id
        self.energy = 0.0
        self.water = 0.0
        self._received: dict[str, float] = {}
        self._interval = (0.0, 0.0)
        self._initial = (0.0, 0.0)
        self._budgets: dict[str, StateBudget] = {}

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        self.energy = float(initial_state.get(LATERAL_ENERGY_STATE, 0.0))
        self.water = float(initial_state.get(LATERAL_WATER_STATE, 0.0))
        self.mark_initialized()

    def compute_tendencies(self) -> list[TendencyContribution]:
        ctx = self.context
        self._interval = (ctx.time, ctx.time + ctx.dt)
        self._initial = (self.energy, self.water)
        energy = self.input_field("lateral_energy_export").require_units("W m^-2").scalar()
        water = self.input_field("lateral_water_export").require_units("kg m^-2 s^-1").scalar()
        latent = self.input_field("lateral_latent_energy_export").require_units("W m^-2").scalar()
        self._received = {
            "lateral_energy": energy * ctx.dt,
            "lateral_water": water * ctx.dt,
            "lateral_latent_energy": latent * ctx.dt,
        }
        return [
            TendencyContribution(
                state_id=LATERAL_ENERGY_STATE,
                provider_module_id=self.module_id,
                label="lateral_energy_accumulation",
                value=energy + latent,
                units="W m^-2",
                equation_id="EQ-BOUNDARY-ACCUMULATION",
            ),
            TendencyContribution(
                state_id=LATERAL_WATER_STATE,
                provider_module_id=self.module_id,
                label="lateral_water_accumulation",
                value=water,
                units="kg m^-2 s^-1",
                equation_id="EQ-BOUNDARY-ACCUMULATION",
            ),
        ]

    def advance_state(self, dt: float, contributions: list[TendencyContribution]) -> None:
        start = {LATERAL_ENERGY_STATE: self.energy, LATERAL_WATER_STATE: self.water}
        totals = dict.fromkeys(start, 0.0)
        for c in contributions:
            totals[c.state_id] += float(np.asarray(c.value))
        self.energy = start[LATERAL_ENERGY_STATE] + totals[LATERAL_ENERGY_STATE] * dt
        self.water = start[LATERAL_WATER_STATE] + totals[LATERAL_WATER_STATE] * dt
        final = {LATERAL_ENERGY_STATE: self.energy, LATERAL_WATER_STATE: self.water}
        self._budgets = {
            sid: StateBudget(
                state_id=sid,
                initial_value=np.asarray(start[sid], dtype=float),
                final_value=np.asarray(final[sid], dtype=float),
                dt=dt,
                contributions=[c for c in contributions if c.state_id == sid],
            )
            for sid in start
        }

    def state_budgets(self) -> dict[str, Any]:
        return dict(self._budgets)

    def emit_outputs(self) -> dict[str, Field]:
        ctx = self.context
        provenance = Provenance(producer=self.module_id, derivation="EQ-BOUNDARY-ACCUMULATION")
        return {
            LATERAL_ENERGY_STATE: Field(
                field_id=LATERAL_ENERGY_STATE,
                physical_name="cumulative energy exported laterally",
                values=np.array(self.energy),
                units="J m^-2",
                quantity_kind=QuantityKind.EXTENSIVE_STATE,
                time_semantics=TimeSemantics.START_OF_STEP,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                provenance=provenance,
            ),
            LATERAL_WATER_STATE: Field(
                field_id=LATERAL_WATER_STATE,
                physical_name="cumulative water exported laterally",
                values=np.array(self.water),
                units="kg m^-2",
                quantity_kind=QuantityKind.EXTENSIVE_STATE,
                time_semantics=TimeSemantics.START_OF_STEP,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                provenance=provenance,
            ),
        }

    def report_transfers(self) -> tuple[list[Transfer], list[Transfer]]:
        start, end = self._interval
        senders = {
            "lateral_energy": (self.atmosphere_module_id, ConservedQuantity.ENERGY, "J m^-2"),
            "lateral_water": (self.water_module_id, ConservedQuantity.WATER, "kg m^-2"),
            "lateral_latent_energy": (self.water_module_id, ConservedQuantity.ENERGY, "J m^-2"),
        }
        imports = []
        for channel, amount in self._received.items():
            sender, quantity, units = senders[channel]
            outward = amount >= 0.0
            imports.append(
                Transfer(
                    transfer_id=f"{channel}@{start:g}",
                    sender_module=sender if outward else self.module_id,
                    receiver_module=self.module_id if outward else sender,
                    quantity=quantity,
                    amount=abs(amount),
                    units=units,
                    start_time=start,
                    end_time=end,
                    channel=channel,
                )
            )
        return ([], imports)

    def report_conservation(self) -> list[ModuleBudget]:
        e0, w0 = self._initial
        return [
            ModuleBudget(
                module_id=self.module_id,
                quantity=ConservedQuantity.ENERGY,
                units="J m^-2",
                stored_initial=e0,
                stored_final=self.energy,
            ),
            ModuleBudget(
                module_id=self.module_id,
                quantity=ConservedQuantity.WATER,
                units="kg m^-2",
                stored_initial=w0,
                stored_final=self.water,
            ),
        ]

    def report_uncertainty(self) -> UncertaintyReport:
        return UncertaintyReport(
            module_id=self.module_id,
            method="none; the sink inherits the uncertainty of the fluxes it receives",
            unresolved=("inherited_flux_uncertainty",),
        )

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "atmosphere_module_id": self.atmosphere_module_id,
            "water_module_id": self.water_module_id,
        }

    def checkpoint(self) -> dict[str, Any]:
        return {LATERAL_ENERGY_STATE: self.energy, LATERAL_WATER_STATE: self.water}

    def restore(self, checkpoint: dict[str, Any]) -> None:
        self.energy = float(checkpoint[LATERAL_ENERGY_STATE])
        self.water = float(checkpoint[LATERAL_WATER_STATE])
