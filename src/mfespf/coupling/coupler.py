"""Multirate coupler.

Charter sections 4.7 and 11, requirements TIM-001 through TIM-004, CAU-007 and
KER-011.

The coupler schedules exchanges; it does not add physics (charter 4.7).  Every
number it moves came from a module and is recorded as a transfer.  The one
transformation it performs is temporal: modules advancing at different rates see
inputs held at their last exchanged value, which is a declared interpolation
choice (``zero_order_hold``) with a declared consequence -- it is first-order
accurate in the coupling interval, and TIM-004 requires that be reported rather
than assumed harmless.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from dataclasses import field as dc_field

from mfespf.contracts.module import Module, StepContext, UncertaintyReport
from mfespf.core.errors import ConfigurationError
from mfespf.core.quantities import Field
from mfespf.coupling.conservation import (
    ConservationAuditor,
    ConservationReport,
    ModuleBudget,
    apply_ledger_to_budgets,
)
from mfespf.coupling.graph import DependencyGraph
from mfespf.coupling.transfers import ConservedQuantity, TransferLedger
from mfespf.registries.states import StateBudget, StateRegistry, TendencyContribution

__all__ = ["Coupler", "CouplingSchedule", "ModuleSchedule", "StepReport"]


@dataclass
class ModuleSchedule:
    """One module's place in the coupling schedule (TIM-002, TIM-003)."""

    module: Module
    step_multiple: int = 1
    interpolation: str = "zero_order_hold"
    accumulation: str = "interval_sum"

    def __post_init__(self) -> None:
        if self.step_multiple < 1:
            raise ConfigurationError(
                f"module {self.module.module_id!r} has step_multiple "
                f"{self.step_multiple}; it must be a positive integer number of master steps"
            )

    def runs_at(self, step_index: int) -> bool:
        return step_index % self.step_multiple == 0

    @property
    def module_id(self) -> str:
        return self.module.module_id


@dataclass
class CouplingSchedule:
    """The declaration TIM-003 requires of a multirate configuration."""

    master_dt: float
    schedules: list[ModuleSchedule]
    splitting_order: tuple[str, ...] = ()
    synchronization_points: tuple[str, ...] = ("end_of_master_step",)
    convergence_criterion: str = "explicit_single_pass"
    master_time_units: str = "s"

    def __post_init__(self) -> None:
        if self.master_dt <= 0:
            raise ConfigurationError("master_dt must be positive")
        ids = [s.module_id for s in self.schedules]
        if len(set(ids)) != len(ids):
            raise ConfigurationError(f"duplicate module IDs in coupling schedule: {ids}")
        if not self.splitting_order:
            self.splitting_order = tuple(ids)
        missing = set(ids) - set(self.splitting_order)
        unknown = set(self.splitting_order) - set(ids)
        if missing or unknown:
            raise ConfigurationError(
                f"splitting_order must list every scheduled module exactly once "
                f"(missing {sorted(missing)}, unknown {sorted(unknown)})"
            )

    @property
    def coupling_period(self) -> int:
        """Master steps between synchronization points (charter 8.4, TIM-003).

        The least common multiple of the module step multiples: the shortest
        interval over which every scheduled module has advanced a whole number of
        its own steps, and therefore the shortest interval over which the two
        sides of every exchange can be compared.
        """
        period = 1
        for s in self.schedules:
            period = math.lcm(period, s.step_multiple)
        return period

    def ordered(self) -> list[ModuleSchedule]:
        by_id = {s.module_id: s for s in self.schedules}
        return [by_id[mid] for mid in self.splitting_order]

    def module_dt(self, module_id: str) -> float:
        for s in self.schedules:
            if s.module_id == module_id:
                return self.master_dt * s.step_multiple
        raise ConfigurationError(f"module {module_id!r} is not scheduled")

    def describe(self) -> str:
        lines = [
            f"master dt: {self.master_dt:g} {self.master_time_units}",
            f"splitting order: {' -> '.join(self.splitting_order)}",
            f"convergence: {self.convergence_criterion}",
            f"synchronization every {self.coupling_period} master steps",
        ]
        for s in self.ordered():
            lines.append(
                f"  {s.module_id:<24} dt={self.module_dt(s.module_id):g} "
                f"(every {s.step_multiple} master steps), {s.interpolation}"
            )
        return "\n".join(lines)


@dataclass
class StepReport:
    """What one master step produced."""

    step_index: int
    time: float
    dt: float
    ran: tuple[str, ...]
    held: tuple[str, ...]
    conservation: ConservationReport | None
    is_synchronization_point: bool = False
    state_budgets: dict[str, StateBudget] = dc_field(default_factory=dict)
    outputs: dict[str, Field] = dc_field(default_factory=dict)
    warnings: list[str] = dc_field(default_factory=list)


class Coupler:
    """Drives the charter 6.3 lifecycle across a set of scheduled modules."""

    def __init__(
        self,
        *,
        schedule: CouplingSchedule,
        state_registry: StateRegistry,
        graph: DependencyGraph | None = None,
        auditor: ConservationAuditor | None = None,
        start_time: float = 0.0,
    ) -> None:
        self.schedule = schedule
        self.states = state_registry
        self.graph = graph or DependencyGraph()
        self.auditor = auditor or ConservationAuditor()
        self.time = start_time
        self.step_index = 0
        self.bus: dict[str, Field] = {}
        self.uncertainty: dict[str, UncertaintyReport] = {}
        self.warnings: list[str] = []
        self._held_since: dict[str, float] = {}
        self._max_staleness: dict[str, float] = {}
        # Accumulated across the current coupling interval; reconciled and reset
        # at each synchronization point (charter 8.4).
        self._ledger = TransferLedger()
        self._budgets: dict[tuple[str, ConservedQuantity], ModuleBudget] = {}
        self._interval_start = start_time
        self.reports: list[ConservationReport] = []

    # --- setup ------------------------------------------------------------
    def validate_configuration(self) -> None:
        """Every pre-initialisation gate, in the order a failure is cheapest.

        Contract completeness first (MOD-001..005), then ownership (STA-001,
        STA-004, STA-005), then the dependency graph's undeclared cycles
        (CAU-005).  SYS-006 requires this to fail rather than proceed.
        """
        active_ids = {s.module_id for s in self.schedule.schedules if s.module.status.participates}
        for s in self.schedule.ordered():
            s.module.validate_configuration()

        declarations = {
            s.module_id: s.module.contract.owned_state_ids
            for s in self.schedule.ordered()
            if s.module.status.evolves_state
        }
        self.states.validate_declared_ownership(declarations)
        self.states.validate_configuration(active_ids)
        self.graph.validate()

    def initialize(
        self,
        initial_state: dict[str, object] | None = None,
        boundary_state: dict[str, object] | None = None,
    ) -> None:
        for s in self.schedule.ordered():
            if not s.module.status.participates:
                continue
            s.module.initialize(dict(initial_state or {}), dict(boundary_state or {}))
        # The bus stays empty until the first step. A field emitted outside a
        # step would have no interval to declare, and DAT-006 does not allow an
        # interval-averaged quantity without one; modules publish in splitting
        # order within the step instead.

    # --- stepping ---------------------------------------------------------
    def step(self) -> StepReport:
        """Advance the coupled system by one master step.

        Transfers and budgets accumulate into the current coupling interval and
        are reconciled when every module has completed a whole number of its own
        steps.  Reconciling every master step would compare a slow module's
        single daily record against a fast module's first hourly one and report a
        mismatch that is an artifact of the schedule, not of the physics.
        """
        t0 = self.time
        dt = self.schedule.master_dt
        contributions: list[TendencyContribution] = []
        ran: list[str] = []
        held: list[str] = []
        step_warnings: list[str] = []

        for s in self.schedule.ordered():
            module = s.module
            if not module.status.participates:
                continue
            if not s.runs_at(self.step_index):
                held.append(module.module_id)
                self._held_since[module.module_id] = (self.step_index % s.step_multiple) * dt
                continue

            ran.append(module.module_id)
            module_dt = self.schedule.module_dt(module.module_id)
            context = StepContext(
                time=t0,
                dt=module_dt,
                step_index=self.step_index,
                master_time_units=self.schedule.master_time_units,
            )
            module.begin_step(context)
            module.receive_inputs(self._inputs_for(module))
            new_contributions = module.compute_tendencies()
            for c in new_contributions:
                self.states.assert_may_contribute(c)
            contributions.extend(new_contributions)
            # Publish before the next module in the splitting order receives its
            # inputs: the order is declared, so downstream modules see this
            # step's diagnostics rather than last step's.
            self.bus.update(module.emit_outputs())

        # Every module that ran has published its tendencies before any state is
        # committed: the operator splitting order is declared, not incidental.
        state_budgets: dict[str, StateBudget] = {}
        for s in self.schedule.ordered():
            module = s.module
            if module.module_id not in ran or not module.status.evolves_state:
                continue
            owned = [c for c in contributions if c.state_id in module.contract.owned_state_ids]
            module.advance_state(self.schedule.module_dt(module.module_id), owned)

        for s in self.schedule.ordered():
            module = s.module
            if module.module_id not in ran:
                continue
            exports, imports = module.report_transfers()
            for t in exports:
                self._ledger.record_export(t)
            for t in imports:
                self._ledger.record_import(t)
            for budget in module.report_conservation():
                key = (budget.module_id, budget.quantity)
                existing = self._budgets.get(key)
                if existing is None:
                    self._budgets[key] = budget
                else:
                    existing.merge_later(budget)
            self.uncertainty[module.module_id] = module.report_uncertainty()
            state_budgets.update(module.state_budgets())

        self.time = t0 + dt
        self.step_index += 1

        report: ConservationReport | None = None
        is_sync = self.step_index % self.schedule.coupling_period == 0
        if is_sync:
            interval = (self._interval_start, self.time)
            apply_ledger_to_budgets(self._ledger, self._budgets)
            report = self.auditor.reconcile(self._ledger, list(self._budgets.values()), interval)
            for transfer in self._ledger.all_transfers():
                if not transfer.within(interval):
                    report.unmatched_transfers.append(
                        f"{transfer.transfer_id} covers [{transfer.start_time:g}, "
                        f"{transfer.end_time:g}], outside the coupling interval "
                        f"[{interval[0]:g}, {interval[1]:g}]"
                    )
            self.reports.append(report)
            self._ledger = TransferLedger()
            self._budgets = {}
            self._interval_start = self.time

        if held:
            for module_id in held:
                self._max_staleness[module_id] = max(
                    self._max_staleness.get(module_id, 0.0), self._held_since.get(module_id, 0.0)
                )
        if is_sync and self._max_staleness:
            # Reported once per coupling interval rather than once per step: the
            # fact TIM-004 wants recorded is the size of the interpolation
            # approximation, and repeating it every hour buries it.
            detail = ", ".join(
                f"{m} up to {v:g} {self.schedule.master_time_units}"
                for m, v in sorted(self._max_staleness.items())
            )
            step_warnings.append(
                f"zero-order hold applied between exchanges ({detail}); the coupling "
                f"interpolation is first-order in the exchange interval (TIM-004)"
            )
            self._max_staleness = {}
        self.warnings.extend(step_warnings)

        return StepReport(
            step_index=self.step_index - 1,
            time=t0,
            dt=dt,
            ran=tuple(ran),
            held=tuple(held),
            conservation=report,
            is_synchronization_point=is_sync,
            state_budgets=state_budgets,
            outputs=dict(self.bus),
            warnings=step_warnings,
        )

    def run(self, n_steps: int, *, assert_conservation: bool = True) -> list[StepReport]:
        reports = []
        for _ in range(n_steps):
            report = self.step()
            if assert_conservation and report.conservation is not None:
                self.auditor.assert_valid(report.conservation)
            reports.append(report)
        if reports and reports[-1].conservation is None:
            self.warnings.append(
                f"the run ended {self.step_index % self.schedule.coupling_period} master steps "
                f"into an unreconciled coupling interval; its transfers are unaudited. Run a "
                f"whole number of coupling periods ({self.schedule.coupling_period} steps) to "
                f"close the budget (CON-009)."
            )
        return reports

    def finalize(self) -> None:
        for s in self.schedule.ordered():
            if s.module.status.participates:
                s.module.finalize()

    # --- internals --------------------------------------------------------
    def _inputs_for(self, module: Module) -> dict[str, Field]:
        """Deliver declared inputs, honouring any disabled dependency edge.

        A disabled edge (CAU-007) withholds the field entirely rather than
        substituting zero, so a module that genuinely requires the input fails
        loudly instead of silently running a different experiment (DAT-007).
        """
        delivered: dict[str, Field] = {}
        for field_id in module.contract.required_inputs:
            field = self.bus.get(field_id)
            if field is None:
                continue
            producer = field.provenance.producer
            if producer and not self.graph.is_enabled(producer, module.module_id, field_id):
                continue
            delivered[field_id] = field
        # Optional inputs a module can use but does not require.
        for field_id, field in self.bus.items():
            if field_id not in delivered and field_id not in module.contract.required_inputs:
                delivered[field_id] = field
        return delivered
