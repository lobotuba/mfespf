"""Assembled configuration and run driver.

Charter section 20 steps 11-15; requirements SYS-003, SYS-006, SYS-008, KER-001
through KER-007, AUD-001, AUD-002, AUD-008.

``Kernel`` is the object a run is: the registries, the schedule, the graph, the
budget and the manifest, assembled and validated together.  The order in
``prepare`` is the charter's: validate the configuration, estimate the cost, then
initialise.  Nothing allocates state until every declared gate has passed, so a
run that cannot be valid never produces output that might later be mistaken for
evidence.
"""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field as dc_field
from typing import TYPE_CHECKING, Any

from mfespf.contracts.cost import CostEstimate, ResourceBudget
from mfespf.core.errors import MFESPFError
from mfespf.coupling.conservation import ConservationReport
from mfespf.coupling.coupler import Coupler, CouplingSchedule, StepReport
from mfespf.coupling.graph import DependencyGraph
from mfespf.coupling.validity import DomainReport, EvaluatedExtents, check_validity_domains
from mfespf.registries.claims import ClaimRegistry, ExperimentRegistry
from mfespf.registries.datasets import DatasetRegistry
from mfespf.registries.equations import EquationRegistry
from mfespf.registries.states import StateRegistry
from mfespf.run.manifest import RunManifest, RunStatus
from mfespf.run.modes import OperatingMode
from mfespf.run.preflight import assert_within_budget, estimate_cost

if TYPE_CHECKING:  # pragma: no cover - a type-only import, to keep run free of data
    from mfespf.data.adapter import HistoricalDataAdapter

__all__ = ["Kernel", "RunResult"]


@dataclass
class RunResult:
    """What a completed run hands back, including its own admissibility."""

    manifest: RunManifest
    steps: list[StepReport] = dc_field(default_factory=list)
    final_conservation: ConservationReport | None = None
    validity_domains: DomainReport | None = None
    """Which declared equation domains the run actually stayed inside (EQR-003).

    A **mandatory** gate, so a violation invalidates the run for inference. That
    is the harsher of the two readings and it is the right one: a run outside a
    declared domain is producing numbers from a relation nobody established for
    those conditions, and the framework has already been burnt five times by
    treating that as a caveat. Making it advisory would rebuild precisely the
    situation this was written to end.
    """
    cost_estimate: CostEstimate | None = None
    warnings: list[str] = dc_field(default_factory=list)

    @property
    def may_support_a_claim(self) -> bool:
        return self.manifest.may_support_a_claim

    @property
    def operating_mode(self) -> OperatingMode:
        return self.manifest.operating_mode

    def assert_may_be_reported_as(self, claimed: OperatingMode) -> None:
        """TIM-006 / VAL-007: refuse a description the run does not support."""
        self.manifest.assert_may_be_reported_as(claimed)

    @property
    def conservation_reports(self) -> list[ConservationReport]:
        return [s.conservation for s in self.steps if s.conservation is not None]

    def worst_relative_residual(self) -> float:
        return max((r.worst_relative_residual() for r in self.conservation_reports), default=0.0)

    def describe(self) -> str:
        lines = [self.manifest.describe(), f"  steps run:    {len(self.steps)}"]
        lines.append(
            f"  worst relative conservation residual: {self.worst_relative_residual():.3e}"
        )
        if self.cost_estimate is not None:
            lines.append("  cost estimate:")
            lines.extend(f"    {line}" for line in self.cost_estimate.describe().splitlines())
        for w in dict.fromkeys(self.warnings):
            lines.append(f"  warning: {w}")
        return "\n".join(lines)


@dataclass
class Kernel:
    """A validated, runnable configuration together with its registries (KER-001)."""

    run_id: str
    purpose: str
    schedule: CouplingSchedule
    states: StateRegistry
    equations: EquationRegistry
    datasets: DatasetRegistry = dc_field(default_factory=DatasetRegistry)
    claims: ClaimRegistry = dc_field(default_factory=ClaimRegistry)
    experiments: ExperimentRegistry = dc_field(default_factory=ExperimentRegistry)
    graph: DependencyGraph = dc_field(default_factory=DependencyGraph)
    budget: ResourceBudget = dc_field(default_factory=ResourceBudget)
    manifest: RunManifest | None = None
    initial_state: dict[str, Any] = dc_field(default_factory=dict)
    boundary_state: dict[str, Any] = dc_field(default_factory=dict)
    conservation_tolerance: float = 1e-10
    spectral_bins: int = 1
    adapter: HistoricalDataAdapter | None = None
    """The declared record driving this run, when there is one.

    A field rather than an attribute bolted on afterwards: a replay's record is
    part of what the run *is*, and a reader of the kernel should not have to
    know that some configurations attach one and others do not.
    """
    start_time: float = 0.0
    """Where the master clock begins (TIM-001).

    Zero for a controlled experiment, whose time origin is arbitrary. A replay
    of a dated record needs the run to start where the record does, or every
    lookup lands outside its coverage.
    """
    operating_mode: OperatingMode = OperatingMode.CONTROLLED_EXPERIMENT
    coupler: Coupler | None = None

    def declared_boundary_state(self) -> dict[str, Any]:
        """Boundary conditions as the modules receive them, clock origin included."""
        return {**self.boundary_state, "master_start_time": self.start_time}

    def build_manifest(self) -> RunManifest:
        """Assemble the manifest from what the configuration actually declares.

        Nothing here is typed in by hand: every entry is read back off the
        modules, registries and schedule, so a manifest cannot describe a run
        that was not configured that way.
        """
        modules = [s.module for s in self.schedule.ordered()]
        manifest = self.manifest or RunManifest(run_id=self.run_id, purpose=self.purpose)
        manifest.operating_mode = self.operating_mode
        manifest.module_status = {m.module_id: m.status.value for m in modules}
        manifest.module_implementations = {m.module_id: m.contract.module_name for m in modules}
        manifest.module_versions = {m.module_id: str(m.contract.module_version) for m in modules}
        manifest.interface_versions = {
            f"{m.module_id}:{iface.interface_id}": str(iface.version)
            for m in modules
            for iface in m.contract.provided_interfaces
        }
        manifest.equation_bindings = self.equations.bindings()
        manifest.parameters = {m.module_id: m.parameters() for m in modules}
        manifest.state_ownership = {
            entry.state_id: entry.owner_module_id for entry in self.states.current()
        }
        manifest.datasets = {
            f"{d.dataset_id}@{d.version}": d.checksum for d in self.datasets.current()
        }
        manifest.dataset_roles = {
            f"{u.dataset_id}@{u.dataset_version}:{u.claim_id}": u.role.value
            for u in self.datasets.uses()
        }
        manifest.coupling_schedule = {
            "master_dt": self.schedule.master_dt,
            "master_time_units": self.schedule.master_time_units,
            "splitting_order": list(self.schedule.splitting_order),
            "convergence_criterion": self.schedule.convergence_criterion,
            "module_steps": {
                s.module_id: {
                    "dt": self.schedule.module_dt(s.module_id),
                    "step_multiple": s.step_multiple,
                    "interpolation": s.interpolation,
                    "accumulation": s.accumulation,
                }
                for s in self.schedule.ordered()
            },
        }
        manifest.temporal_resolution = f"{self.schedule.master_dt:g} s master step"
        manifest.spectral_resolution = f"{self.spectral_bins} bins"
        manifest.tolerances = {"conservation_relative": self.conservation_tolerance}
        manifest.uncertainty_methods = {m.module_id: m.report_uncertainty().method for m in modules}
        manifest.resource_class = self.budget.resource_class
        manifest.resource_budget = {
            "max_memory_bytes": self.budget.max_memory_bytes,
            "max_work_units": self.budget.max_work_units,
            "max_output_bytes": self.budget.max_output_bytes,
            "on_exceed": self.budget.on_exceed,
        }
        manifest.initial_conditions = dict(self.initial_state)
        # The clock origin is a boundary condition of the run, not a detail of
        # the driver: two replays over different periods of one record are
        # different configurations and must not hash alike.
        manifest.boundary_conditions = self.declared_boundary_state()
        manifest.dependency_graph = self.graph.to_dict()
        self.manifest = manifest
        return manifest

    def prepare(self, *, n_steps: int) -> tuple[Coupler, CostEstimate]:
        """Run every pre-initialisation gate, then initialise (SYS-006, CMP-005)."""
        manifest = self.build_manifest()
        coupler = Coupler(
            schedule=self.schedule,
            state_registry=self.states,
            graph=self.graph,
            start_time=self.start_time,
        )
        coupler.auditor.tolerance = self.conservation_tolerance

        try:
            coupler.validate_configuration()
            manifest.record_gate("configuration_validation", True)
        except MFESPFError as exc:
            manifest.record_gate("configuration_validation", False)
            manifest.status = RunStatus.FAILED
            manifest.warnings.append(str(exc))
            raise

        estimate = estimate_cost(
            schedule=self.schedule,
            n_steps=n_steps,
            budget=self.budget,
            spectral_bins=self.spectral_bins,
        )
        try:
            assert_within_budget(estimate, self.budget)
            manifest.record_gate("resource_preflight", True)
        except MFESPFError:
            manifest.record_gate("resource_preflight", False)
            manifest.status = RunStatus.FAILED
            raise
        manifest.warnings.extend(estimate.warnings)

        coupler.initialize(self.initial_state, self.declared_boundary_state())
        self.coupler = coupler
        return coupler, estimate

    def run(self, n_steps: int, *, assert_conservation: bool = True) -> RunResult:
        """Execute the configuration and return a result that knows its own status."""
        coupler, estimate = self.prepare(n_steps=n_steps)
        manifest = self.manifest
        assert manifest is not None  # build_manifest ran inside prepare

        steps: list[StepReport] = []
        try:
            steps = coupler.run(n_steps, assert_conservation=assert_conservation)
        except MFESPFError as exc:
            manifest.record_gate("conservation_reconciliation", False)
            manifest.status = RunStatus.FAILED
            manifest.warnings.append(str(exc))
            coupler.finalize()
            # The steps completed before the failure are kept. They are not
            # evidence -- the manifest marks the run unusable -- but discarding
            # them would leave nothing to diagnose the failure from.
            steps = coupler.step_reports
            return RunResult(
                manifest=manifest,
                steps=steps,
                cost_estimate=estimate,
                warnings=[*coupler.warnings, str(exc)],
            )

        coupler.finalize()

        # Judged on the extent each relation was actually evaluated over, which
        # is a property of the whole trajectory rather than of its last step: a
        # closure that ran out of range while spinning up and recovered was
        # still evaluated out of range, and one that is merely idle at the end
        # leaves no footprint on a final state to be judged on at all.
        extents = EvaluatedExtents()
        for step in steps:
            extents.observe(step.outputs)
        domains = check_validity_domains(
            registry=self.equations,
            modules=[s.module for s in self.schedule.ordered()],
            fields=steps[-1].outputs if steps else {},
            extents=extents,
        )
        manifest.record_gate("validity_domains", domains.passed)
        manifest.warnings.extend(domains.describe_violations())

        reports = [s.conservation for s in steps if s.conservation is not None]
        final = reports[-1] if reports else None

        # ACC-005 and KER-014: a run that cannot be fingerprinted cannot be shown
        # to have been reproduced. The endpoint alone would not distinguish two
        # trajectories that met at the end, so the residual series goes in beside
        # it -- a scalar per step, which costs nothing and fingerprints the path.
        if steps:
            manifest.record_computed_output(
                "final_outputs", {k: f.values for k, f in steps[-1].outputs.items()}
            )
            manifest.record_computed_output(
                "conservation_trajectory", [r.worst_relative_residual() for r in reports]
            )
        manifest.record_computed_output(
            "final_state",
            {
                s.module.module_id: {k: f.values for k, f in s.module.emit_state().items()}
                for s in self.schedule.ordered()
            },
        )
        passed = bool(reports) and all(r.passed for r in reports)
        manifest.record_gate("conservation_reconciliation", passed)
        if manifest.status == RunStatus.PENDING:
            manifest.status = RunStatus.COMPLETE
        manifest.warnings.extend(coupler.warnings)

        return RunResult(
            manifest=manifest,
            steps=steps,
            final_conservation=final,
            validity_domains=domains,
            cost_estimate=estimate,
            warnings=list(coupler.warnings),
        )
