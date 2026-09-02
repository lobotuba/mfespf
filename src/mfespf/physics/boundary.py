"""System boundary sink: space.

Charter 4.21 and 8.1; requirements CON-001, CON-009, STA-001.

Energy leaving the modeled system has to go somewhere accountable, or the
top-level budget cannot be checked.  This module is that somewhere.  It owns a
cumulative reservoir of everything radiated or reflected to space and records its
own side of each transfer, computed from the fluxes the surface published rather
than from the surface's internal variables.  That independence is what makes the
CON-004 reconciliation a real check instead of a restatement.

The module claims no physics beyond bookkeeping, and its fidelity level says so.
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
from mfespf.registries.states import (
    StateBudget,
    StateEntry,
    StateKind,
    TendencyContribution,
    UpdateAuthority,
)

__all__ = ["SPACE_ENERGY_STATE", "SpaceBoundaryModule", "space_state_entry"]

V1 = Version(1, 0, 0)
SPACE_ENERGY_STATE = "space_cumulative_energy_received"

SPACE_BOUNDARY_INTERFACE = InterfaceContract(
    interface_id="system_boundary_sink",
    version=V1,
    physical_semantics={
        "energy_units": "J m^-2",
        "accounting_basis": "cumulative_per_planetary_surface_area",
        "sign_convention": "positive_into_receiver",
    },
)


def space_state_entry(*, owner: str = "space") -> StateEntry:
    return StateEntry(
        record_id=SPACE_ENERGY_STATE,
        version=V1,
        physical_name="cumulative energy radiated and reflected to space",
        symbol="E_space",
        owner_module_id=owner,
        units="J m^-2",
        state_kind=StateKind.PROGNOSTIC,
        conservation_participation=("energy",),
        authorized_tendency_providers=frozenset({owner}),
        update_authority=UpdateAuthority.OWNER_ONLY,
        initialization_source="zero_at_run_start",
        checkpoint_required=True,
        uncertainty_representation="inherits_from_contributing_fluxes",
    )


class SpaceBoundaryModule(Module):
    """Accumulates outgoing energy at the system boundary."""

    def __init__(
        self,
        *,
        module_id: str = "space",
        sender_module_id: str = "surface",
        incoming_flux_fields: tuple[tuple[str, str], ...] = (
            ("surface_reflected_shortwave", "reflected_shortwave"),
            ("surface_emitted_longwave", "emitted_longwave"),
        ),
        status: ModuleStatus = ModuleStatus.ACTIVE,
    ) -> None:
        contract = ModuleContract(
            module_id=module_id,
            module_name="system boundary energy sink",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.NULL_CONTROL,
            physical_basis=(
                "Accounting boundary. Energy leaving the modeled system is accumulated so "
                "that the top-level budget can be closed and reported."
            ),
            governing_equation_ids=("EQ-BOUNDARY-ACCUMULATION",),
            assumptions=(
                "space returns no energy to the system at this fidelity",
                "no cosmic microwave background or interplanetary source is represented",
            ),
            validity_domain=(
                "Bookkeeping only. It makes no physical prediction and cannot support a "
                "claim about anything other than closure of the modeled budget."
            ),
            owned_state_ids=frozenset({SPACE_ENERGY_STATE}),
            required_inputs=tuple(f for f, _ in incoming_flux_fields),
            produced_outputs=("space_cumulative_energy_received", "outgoing_total_flux"),
            imported_transfers=tuple(c for _, c in incoming_flux_fields),
            conserved_quantities=("energy",),
            provided_interfaces=(SPACE_BOUNDARY_INTERFACE,),
            uncertainty_outputs=("inherited_from_contributing_fluxes",),
            correctness_tests=(
                "test_space_accumulates_exactly_what_surface_exports",
                "test_zero_insolation_null_control",
            ),
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
                checkpoint_cost_bytes=8,
            ),
        )
        super().__init__(contract, status=status)
        self.sender_module_id = sender_module_id
        self.incoming_flux_fields = incoming_flux_fields
        self.cumulative_energy = 0.0
        self._received: dict[str, float] = {}
        self._interval = (0.0, 0.0)
        self._stored_initial = 0.0
        self._budget: StateBudget | None = None

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        self.cumulative_energy = float(initial_state.get(SPACE_ENERGY_STATE, 0.0))
        self.mark_initialized()

    def compute_tendencies(self) -> list[TendencyContribution]:
        ctx = self.context
        self._interval = (ctx.time, ctx.time + ctx.dt)
        self._stored_initial = self.cumulative_energy
        self._received = {
            channel: self.input_field(field_id).require_units("W m^-2").scalar() * ctx.dt
            for field_id, channel in self.incoming_flux_fields
        }
        total_rate = sum(self._received.values()) / ctx.dt if ctx.dt else 0.0
        return [
            TendencyContribution(
                state_id=SPACE_ENERGY_STATE,
                provider_module_id=self.module_id,
                label="outgoing_radiation_accumulation",
                value=total_rate,
                units="W m^-2",
                equation_id="EQ-BOUNDARY-ACCUMULATION",
            )
        ]

    def advance_state(self, dt: float, contributions: list[TendencyContribution]) -> None:
        initial = np.asarray(self.cumulative_energy, dtype=float)
        total = sum(float(np.asarray(c.value)) for c in contributions)
        self.cumulative_energy = float(initial + total * dt)
        self._budget = StateBudget(
            state_id=SPACE_ENERGY_STATE,
            initial_value=initial,
            final_value=np.asarray(self.cumulative_energy, dtype=float),
            dt=dt,
            contributions=list(contributions),
        )

    def state_budgets(self) -> dict[str, Any]:
        return {SPACE_ENERGY_STATE: self._budget} if self._budget else {}

    def emit_outputs(self) -> dict[str, Field]:
        ctx = self.context
        provenance = Provenance(producer=self.module_id, derivation="EQ-BOUNDARY-ACCUMULATION")
        outgoing = sum(self._received.values()) / ctx.dt if ctx.dt else 0.0
        return {
            "space_cumulative_energy_received": Field(
                field_id="space_cumulative_energy_received",
                physical_name="cumulative energy received at the system boundary",
                symbol="E_space",
                values=np.array(self.cumulative_energy),
                units="J m^-2",
                quantity_kind=QuantityKind.EXTENSIVE_STATE,
                time_semantics=TimeSemantics.START_OF_STEP,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                provenance=provenance,
            ),
            "outgoing_total_flux": Field(
                field_id="outgoing_total_flux",
                physical_name="total outgoing flux density at the top of the system",
                symbol="F_out",
                values=np.array(outgoing),
                units="W m^-2",
                quantity_kind=QuantityKind.FLUX,
                time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                time_interval=ctx.dt,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                sign_convention=SignConvention.POSITIVE_UPWARD,
                provenance=provenance,
            ),
        }

    def report_transfers(self) -> tuple[list[Transfer], list[Transfer]]:
        start, end = self._interval
        imports = [
            Transfer(
                transfer_id=f"{channel}_received@{start:g}",
                sender_module=self.sender_module_id,
                receiver_module=self.module_id,
                quantity=ConservedQuantity.ENERGY,
                amount=amount,
                units="J m^-2",
                start_time=start,
                end_time=end,
                channel=channel,
            )
            for channel, amount in self._received.items()
        ]
        return ([], imports)

    def report_conservation(self) -> list[ModuleBudget]:
        return [
            ModuleBudget(
                module_id=self.module_id,
                quantity=ConservedQuantity.ENERGY,
                units="J m^-2",
                stored_initial=self._stored_initial,
                stored_final=self.cumulative_energy,
            )
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
            "sender_module_id": self.sender_module_id,
            "incoming_channels": [channel for _f, channel in self.incoming_flux_fields],
        }

    def checkpoint(self) -> dict[str, Any]:
        return {SPACE_ENERGY_STATE: self.cumulative_energy}

    def restore(self, checkpoint: dict[str, Any]) -> None:
        self.cumulative_energy = float(checkpoint[SPACE_ENERGY_STATE])
