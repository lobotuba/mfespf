"""Module contract and lifecycle.

Charter sections 6.2-6.4 and requirements MOD-001 through MOD-014.

A module owns part of the model state, advances or diagnoses it, and touches
other modules only through this contract.  The lifecycle below is the charter's
6.3 sequence; a language may rename the steps but not reorder them, because the
coupler's conservation accounting depends on tendencies being computed from a
consistent set of inputs before any state advances.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from dataclasses import field as dc_field
from enum import StrEnum
from typing import Any

from mfespf.contracts.cost import CostDeclaration
from mfespf.core.errors import ConfigurationError
from mfespf.core.quantities import Field
from mfespf.core.versioning import InterfaceContract, Version, VersionRange
from mfespf.registries.states import TendencyContribution

__all__ = [
    "FidelityLevel",
    "Module",
    "ModuleContract",
    "ModuleStatus",
    "StepContext",
    "UncertaintyReport",
]


class ModuleStatus(StrEnum):
    """Charter 6.4 / SYS-003."""

    ACTIVE = "active"
    PRESCRIBED = "prescribed"
    DIAGNOSTIC = "diagnostic"
    INACTIVE = "inactive"
    STUB = "stub"

    @property
    def participates(self) -> bool:
        return self is not ModuleStatus.INACTIVE

    @property
    def evolves_state(self) -> bool:
        return self is ModuleStatus.ACTIVE

    @property
    def may_support_a_claim(self) -> bool:
        """SYS-004: stub output is test data and supports nothing."""
        return self is not ModuleStatus.STUB


class FidelityLevel(StrEnum):
    """Charter 4.28.  Module-specific and explicitly not a ranking across modules."""

    NULL_CONTROL = "null_control"
    PRESCRIBED = "prescribed"
    REDUCED = "reduced"
    STANDARD = "standard"
    HIGH = "high"
    REFERENCE = "reference"


@dataclass(frozen=True)
class ModuleContract:
    """The declaration every module implementation must publish (charter 6.2).

    Every field the charter lists is present.  Several default to an empty tuple
    rather than being optional in the scientific sense: an empty
    ``correctness_tests`` is caught by ``validate`` because MOD-014 requires at
    least one null, limiting-case or disabled-coupling test.
    """

    module_id: str
    module_name: str
    module_version: Version
    interface_version: Version
    fidelity_level: FidelityLevel
    physical_basis: str
    governing_equation_ids: tuple[str, ...]
    assumptions: tuple[str, ...]
    validity_domain: str
    owned_state_ids: frozenset[str] = frozenset()
    required_inputs: tuple[str, ...] = ()
    produced_outputs: tuple[str, ...] = ()
    exported_transfers: tuple[str, ...] = ()
    imported_transfers: tuple[str, ...] = ()
    conserved_quantities: tuple[str, ...] = ()
    internal_sources_and_sinks: tuple[str, ...] = ()
    supported_spatial_representations: tuple[str, ...] = ("global_scalar",)
    supported_temporal_representations: tuple[str, ...] = ("fixed_step",)
    supported_spectral_representations: tuple[str, ...] = ()
    provided_interfaces: tuple[InterfaceContract, ...] = ()
    required_interfaces: tuple[tuple[str, VersionRange], ...] = ()
    uncertainty_outputs: tuple[str, ...] = ()
    correctness_tests: tuple[str, ...] = ()
    known_limitations: tuple[str, ...] = ()
    equation_variables: tuple[tuple[str, tuple[tuple[str, str, str], ...]], ...] = ()
    """Which of this module's quantities plays each equation variable (EQR-003).

    ``(equation_id, ((symbol, name, units), ...))``, where ``name`` is one of
    the module's published field ids or a key of its ``parameters()``.

    The units are stated here rather than taken on trust. A parameter carries no
    units of its own, so binding ``concentration_ppm`` to a variable whose range
    is declared in ``1`` would otherwise compare parts per million against a
    fraction and report confidently -- the silent mismatch this check refuses
    everywhere else. Where the binding names a field, both are compared, so a
    module that mis-states its own units is caught too.

    This exists because an equation writes ``T`` and a column publishes ``T_a``
    and ``T_s``, and only the module knows which of the two the Schwarzschild
    equation meant. The run-time domain check falls back to matching by symbol
    root and units, which is enough for the unambiguous cases and correctly
    refuses the rest; a binding here turns a refusal into a check.

    Tuples rather than a dict because the contract is frozen and hashed into the
    manifest, and two runs that bind differently must not hash alike.
    """
    cost: CostDeclaration | None = None
    fidelity_justification: str = ""
    """Why this fidelity rather than the cheaper one this module offers (CMP-007).

    Required exactly when ``cost.lower_fidelity_alternatives`` is non-empty: a
    module that names a cheaper implementation and then does not use it has made
    a scientific choice, and CMP-007 asks that the reason travel in the manifest
    rather than in whoever configured the run.

    It must say what the cheaper form would *fail to produce*, because that is
    the falsifiable half. "More accurate" is not a justification; "a single layer
    cannot have a lapse rate, and the lapse rate is what convection acts on" is.
    """
    supports_checkpoint: bool = True

    def validate(self) -> None:
        """Reject a contract that cannot satisfy the mandatory module requirements."""
        problems = []
        if not self.governing_equation_ids:
            problems.append("declares no governing equation IDs (MOD-002, EQR-005)")
        if not self.physical_basis.strip():
            problems.append("declares no physical basis (MOD-002)")
        if not self.validity_domain.strip():
            problems.append("declares no validity domain (MOD-002)")
        if not self.correctness_tests:
            problems.append(
                "names no null, limiting-case or disabled-coupling test (MOD-014, VAL-001)"
            )
        if self.cost is None:
            problems.append("provides no computational-cost declaration (MOD-013, CMP-001)")
        elif self.cost.module_id != self.module_id:
            problems.append(
                f"cost declaration names module {self.cost.module_id!r}, not {self.module_id!r}"
            )
        if (
            self.cost is not None
            and self.cost.lower_fidelity_alternatives
            and not self.fidelity_justification.strip()
        ):
            problems.append(
                f"names the cheaper alternatives "
                f"{list(self.cost.lower_fidelity_alternatives)} and gives no reason for "
                f"running at {self.fidelity_level.value} instead (CMP-007)"
            )
        if self.owned_state_ids and not self.conserved_quantities:
            problems.append(
                "owns state but declares no conserved quantity; declare the applicable "
                "quantities or state explicitly that none apply (MOD-004)"
            )
        if problems:
            raise ConfigurationError(
                f"module contract {self.module_id!r} is incomplete:\n  - " + "\n  - ".join(problems)
            )

    def cite(self) -> str:
        return f"{self.module_id}@{self.module_version} (iface {self.interface_version})"


@dataclass(frozen=True, slots=True)
class StepContext:
    """What a module is told about the step it is about to take."""

    time: float
    dt: float
    step_index: int
    master_time_units: str = "s"

    @property
    def end_time(self) -> float:
        return self.time + self.dt


@dataclass
class UncertaintyReport:
    """MOD-009, UNC-001, UNC-004.

    ``unresolved`` is not decoration: a module that cannot quantify a
    contribution it knows applies must say so, and the run report prints it.
    """

    module_id: str
    quantified: dict[str, float] = dc_field(default_factory=dict)
    method: str = "unresolved"
    distributional_assumptions: str = "unresolved"
    correlations: str = "unresolved"
    unresolved: tuple[str, ...] = ()

    def describe(self) -> str:
        lines = [f"{self.module_id}: method={self.method}"]
        for name, value in sorted(self.quantified.items()):
            lines.append(f"  {name} = {value:.6g}")
        for name in self.unresolved:
            lines.append(f"  {name} = UNRESOLVED")
        return "\n".join(lines)


class Module(ABC):
    """Base class implementing the charter 6.3 lifecycle (MOD-007).

    Subclasses override the hooks they need.  The default implementations are
    inert rather than clever: a module that does not export transfers reports an
    empty list, and one that owns no state has nothing to advance.  What they are
    not allowed to do is write another module's state -- that path does not exist
    on this class, only ``compute_tendencies`` does (MOD-010, MOD-011).
    """

    def __init__(
        self, contract: ModuleContract, status: ModuleStatus = ModuleStatus.ACTIVE
    ) -> None:
        self.contract = contract
        self.status = status
        self.configuration: dict[str, Any] = {}
        self._inputs: dict[str, Field] = {}
        self._context: StepContext | None = None
        self._initialized = False

    # --- identity ---------------------------------------------------------
    @property
    def module_id(self) -> str:
        return self.contract.module_id

    def __repr__(self) -> str:  # pragma: no cover - convenience
        return f"<{type(self).__name__} {self.module_id} status={self.status.value}>"

    # --- lifecycle --------------------------------------------------------
    def configure(self, configuration: dict[str, Any]) -> None:
        self.configuration = dict(configuration)

    def validate_configuration(self) -> None:
        """MOD-005: validate before any state exists.

        The base implementation validates the contract itself; subclasses should
        call ``super().validate_configuration()`` and then check their own
        parameters against their declared validity domain.
        """
        self.contract.validate()

    @abstractmethod
    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        """MOD-006: deterministic initialisation from a recorded manifest."""

    def expose_capabilities(self) -> tuple[InterfaceContract, ...]:
        return self.contract.provided_interfaces

    def begin_step(self, context: StepContext) -> None:
        self._context = context
        self._inputs = {}

    def receive_inputs(self, packets: dict[str, Field]) -> None:
        missing = [f for f in self.contract.required_inputs if f not in packets]
        if missing and self.status.evolves_state:
            raise ConfigurationError(
                f"module {self.module_id!r} did not receive required inputs {missing}; "
                f"DAT-007 forbids substituting zero for a missing physical input"
            )
        self._inputs = dict(packets)

    def compute_tendencies(self) -> list[TendencyContribution]:
        """Labeled contributions to owned or foreign state (MOD-011)."""
        return []

    def advance_state(self, dt: float, contributions: list[TendencyContribution]) -> None:
        """Commit owned state from the aggregated labeled contributions.

        Only the owner is ever called here, and it receives every contribution to
        its states including those from other modules (A1.3).
        """
        return None

    def emit_outputs(self) -> dict[str, Field]:
        """Diagnostics and fluxes for the step just computed.

        Called after ``compute_tendencies`` and before the next module in the
        splitting order receives its inputs.
        """
        return {}

    def emit_state(self) -> dict[str, Field]:
        """The module's currently committed state, and any static properties.

        Called once at initialisation and again after every ``advance_state``, so
        the bus always carries the state as of the *start* of the next step. A
        module that needs another's state -- a radiation solver needing
        temperatures -- reads it from here rather than being ordered after the
        owner, which would make its inputs a step stale.

        Fields published here are instantaneous and carry no averaging interval,
        which is what lets them exist before any step has been taken (DAT-006).
        """
        return {}

    def parameters(self) -> dict[str, Any]:
        """Every value that determines this module's numbers (AUD-002).

        The run manifest hashes these, so a module that omits a parameter here
        makes two materially different runs indistinguishable in the audit trail.
        Subclasses override; the base returns whatever was passed to ``configure``.
        """
        return dict(self.configuration)

    def state_budgets(self) -> dict[str, Any]:
        """Per-state update budgets for the step just taken (STA-006).

        Keyed by state_id, valued by ``StateBudget``.  An owner that commits a
        state must be able to show which labeled contributions produced it.
        """
        return {}

    def report_conservation(self) -> list[Any]:
        """Return this module's ``ModuleBudget`` records for the coupling interval.

        A module reports only what it owns: stored quantities and its own labeled
        internal sources and sinks.  The coupler fills in the import and export
        lines from the reconciled ledger, so a module cannot certify its own
        boundary exchange.
        """
        return []

    def report_transfers(self) -> tuple[list[Any], list[Any]]:
        """Return ``(exports, imports)`` as this module independently recorded them.

        Both parties to a transfer record it, and the auditor reconciles the two
        (CON-004).  A module that reports its counterpart's side of an exchange
        defeats the check, so each side reports only its own.
        """
        return ([], [])

    def report_uncertainty(self) -> UncertaintyReport:
        return UncertaintyReport(
            module_id=self.module_id,
            unresolved=self.contract.uncertainty_outputs or ("all_contributions",),
        )

    def checkpoint(self) -> dict[str, Any]:
        """MOD-008.  Must round-trip through ``restore`` bit-for-bit."""
        return {}

    def restore(self, checkpoint: dict[str, Any]) -> None:
        return None

    def finalize(self) -> None:
        return None

    # --- helpers for subclasses ------------------------------------------
    def input_field(self, field_id: str) -> Field:
        try:
            return self._inputs[field_id]
        except KeyError:
            raise ConfigurationError(
                f"module {self.module_id!r} requested input {field_id!r}, which was not delivered "
                f"this step"
            ) from None

    @property
    def context(self) -> StepContext:
        if self._context is None:
            raise ConfigurationError(f"module {self.module_id!r} is not inside a step")
        return self._context

    def mark_initialized(self) -> None:
        self._initialized = True

    @property
    def is_initialized(self) -> bool:
        return self._initialized
