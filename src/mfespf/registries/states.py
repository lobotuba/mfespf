"""State ownership registry.

Amendment A1 and requirements STA-001 through STA-007.

The single rule this module exists to enforce: every prognostic state has
exactly one authoritative owner, and no other module writes it.  Other modules
contribute *labeled* tendencies which the owner aggregates and commits, so the
update budget of every state can be reconstructed after the fact (STA-006).
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from dataclasses import field as dc_field
from enum import StrEnum

import numpy as np

from mfespf.core.errors import ConfigurationError, OwnershipError
from mfespf.core.quantities import SpatialSupport
from mfespf.core.units import parse_unit
from mfespf.core.versioning import Version
from mfespf.registries.base import Registry

__all__ = [
    "StateBudget",
    "StateEntry",
    "StateKind",
    "StateRegistry",
    "TendencyContribution",
    "UpdateAuthority",
]


class UpdateAuthority(StrEnum):
    """Who is permitted to commit a new value for a state (A1.2)."""

    OWNER_ONLY = "owner_only"
    OWNER_WITH_COUPLER_ARBITRATION = "owner_with_coupler_arbitration"
    PRESCRIBED_SOURCE = "prescribed_source"


class StateKind(StrEnum):
    """Charter 4.9-4.12."""

    PROGNOSTIC = "prognostic"
    DIAGNOSTIC = "diagnostic"
    PRESCRIBED = "prescribed"
    CONTROL = "control"


@dataclass(frozen=True)
class StateEntry:
    """One row of the state-ownership registry (A1.2, STA-002, STA-003)."""

    record_id: str
    version: Version
    physical_name: str
    owner_module_id: str
    units: str
    state_kind: StateKind = StateKind.PROGNOSTIC
    symbol: str = ""
    dimensions: tuple[str, ...] = ()
    coordinate_support: SpatialSupport = SpatialSupport.GLOBAL_SCALAR
    phase_or_species: str | None = None
    conservation_participation: tuple[str, ...] = ()
    authorized_readers: frozenset[str] = frozenset()
    authorized_tendency_providers: frozenset[str] = frozenset()
    authorized_boundary_providers: frozenset[str] = frozenset()
    update_authority: UpdateAuthority = UpdateAuthority.OWNER_ONLY
    initialization_source: str = "manifest"
    checkpoint_required: bool = True
    uncertainty_representation: str = "unresolved"
    validity_flags: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        parse_unit(self.units)  # unknown units fail at registration, not at run time
        if self.state_kind is StateKind.PROGNOSTIC and not self.owner_module_id:
            raise ConfigurationError(
                f"prognostic state {self.record_id!r} declares no owner (STA-001)"
            )

    @property
    def state_id(self) -> str:
        return self.record_id

    def may_read(self, module_id: str) -> bool:
        return module_id == self.owner_module_id or module_id in self.authorized_readers

    def may_provide_tendency(self, module_id: str) -> bool:
        return module_id == self.owner_module_id or module_id in self.authorized_tendency_providers


@dataclass(frozen=True, slots=True)
class TendencyContribution:
    """A labeled contribution to another module's state (MOD-011, A1.3).

    ``label`` is the physical process name -- ``radiative_heating``,
    ``surface_evaporation``, ``numerical_limiter`` -- and it survives into the
    state budget, which is what makes STA-006 checkable and CON-008 honest.
    """

    state_id: str
    provider_module_id: str
    label: str
    value: np.ndarray | float
    units: str
    is_numerical: bool = False
    equation_id: str | None = None

    def as_array(self) -> np.ndarray:
        return np.asarray(self.value, dtype=float)


@dataclass
class StateBudget:
    """The reconstructable update record for one state over one step (STA-006)."""

    state_id: str
    initial_value: np.ndarray
    final_value: np.ndarray
    dt: float
    contributions: list[TendencyContribution] = dc_field(default_factory=list)

    @property
    def stored_change(self) -> np.ndarray:
        return self.final_value - self.initial_value

    @property
    def integrated_contributions(self) -> np.ndarray:
        total = np.zeros_like(np.asarray(self.initial_value, dtype=float))
        for c in self.contributions:
            total = total + c.as_array() * self.dt
        return total

    def residual(self) -> np.ndarray:
        """Change not explained by the labeled contributions.

        A non-zero residual means the owner committed a value the labels do not
        account for.  For an explicit forward step this must be zero to round-off;
        for a higher-order or implicit scheme it is the truncation term, and the
        owner is expected to declare it as such.
        """
        return self.stored_change - self.integrated_contributions

    def by_label(self) -> dict[str, np.ndarray]:
        out: dict[str, np.ndarray] = {}
        for c in self.contributions:
            out[c.label] = out.get(c.label, 0.0) + c.as_array() * self.dt
        return out

    def numerical_part(self) -> np.ndarray:
        """The portion attributable to stabilisation, limiting or filtering (CON-008)."""
        total = np.zeros_like(np.asarray(self.initial_value, dtype=float))
        for c in self.contributions:
            if c.is_numerical:
                total = total + c.as_array() * self.dt
        return total


class StateRegistry(Registry[StateEntry]):
    """Registry plus the ownership gates of STA-001, STA-004, STA-005 and MOD-010."""

    def __init__(self) -> None:
        super().__init__("state")

    def owner_of(self, state_id: str) -> str:
        return self.get(state_id).owner_module_id

    def assert_may_write(self, module_id: str, state_id: str) -> None:
        """MOD-010: a module shall not directly write state it does not own."""
        entry = self.get(state_id)
        if module_id != entry.owner_module_id:
            raise OwnershipError(
                f"module {module_id!r} attempted to write state {state_id!r}, which is owned by "
                f"{entry.owner_module_id!r}. Submit a labeled tendency through the coupler instead "
                f"(MOD-011)."
            )

    def assert_may_contribute(self, contribution: TendencyContribution) -> None:
        entry = self.get(contribution.state_id)
        if not entry.may_provide_tendency(contribution.provider_module_id):
            raise OwnershipError(
                f"module {contribution.provider_module_id!r} is not an authorized tendency "
                f"provider for state {contribution.state_id!r} (authorized: "
                f"{sorted(entry.authorized_tendency_providers) or 'none'})"
            )
        expected = parse_unit(entry.units).dimension
        actual = parse_unit(contribution.units).dimension
        # A tendency carries the state's dimension per unit time.
        per_time = expected / parse_unit("s").dimension
        if actual not in (expected, per_time):
            raise OwnershipError(
                f"tendency {contribution.label!r} for {contribution.state_id!r} is in "
                f"{contribution.units!r}; expected {entry.units!r} or its per-second rate"
            )

    def validate_configuration(self, active_module_ids: set[str]) -> None:
        """Reject duplicate or absent ownership before initialisation.

        STA-004, STA-005 and A1.5.  Called by the run driver ahead of
        ``initialize`` so that an unresolvable configuration never produces
        output that might later be mistaken for evidence.
        """
        owners: dict[str, list[str]] = defaultdict(list)
        problems: list[str] = []

        for entry in self.current():
            if entry.state_kind is not StateKind.PROGNOSTIC:
                continue
            if entry.owner_module_id in active_module_ids:
                owners[entry.state_id].append(entry.owner_module_id)
            else:
                problems.append(
                    f"prognostic state {entry.state_id!r} names owner "
                    f"{entry.owner_module_id!r}, which is not active in this "
                    f"configuration (STA-005)"
                )
            unknown = (
                (
                    entry.authorized_readers
                    | entry.authorized_tendency_providers
                    | entry.authorized_boundary_providers
                )
                - active_module_ids
                - {entry.owner_module_id}
            )
            if unknown:
                problems.append(
                    f"state {entry.state_id!r} authorises inactive modules {sorted(unknown)}"
                )

        for state_id, claimants in owners.items():
            if len(claimants) > 1:
                problems.append(
                    f"state {state_id!r} is claimed by {sorted(claimants)}; exactly one owner is "
                    f"permitted (STA-001, STA-004)"
                )

        if problems:
            raise OwnershipError(
                "state ownership validation failed:\n  - " + "\n  - ".join(problems)
            )

    def validate_declared_ownership(self, declarations: dict[str, frozenset[str]]) -> None:
        """Cross-check what modules claim to own against the registry.

        ``declarations`` maps module_id to the state_ids that module declares it
        owns.  Two modules claiming the same state is the A1.5 conflict case, and
        it can only be seen here: the registry itself is keyed by state_id and so
        cannot represent the conflict.  A module claiming a state the registry
        assigns elsewhere is the same failure seen from the other side.
        """
        claims: dict[str, list[str]] = defaultdict(list)
        for module_id, state_ids in declarations.items():
            for state_id in state_ids:
                claims[state_id].append(module_id)

        problems: list[str] = []
        for state_id, claimants in sorted(claims.items()):
            if len(claimants) > 1:
                problems.append(
                    f"state {state_id!r} is claimed by {sorted(claimants)}; exactly one owner is "
                    f"permitted (STA-001, STA-004)"
                )
                continue
            if state_id not in self:
                problems.append(
                    f"module {claimants[0]!r} claims unregistered state {state_id!r}; "
                    f"register it with units, coordinates and authorisations first (STA-002)"
                )
                continue
            registered_owner = self.get(state_id).owner_module_id
            if registered_owner != claimants[0]:
                problems.append(
                    f"module {claimants[0]!r} claims state {state_id!r}, but the registry assigns "
                    f"it to {registered_owner!r}"
                )

        registered_prognostic = {
            e.state_id for e in self.current() if e.state_kind is StateKind.PROGNOSTIC
        }
        orphaned = registered_prognostic - set(claims)
        for state_id in sorted(orphaned):
            problems.append(f"prognostic state {state_id!r} has no claiming module (STA-005)")

        if problems:
            raise OwnershipError(
                "declared state ownership validation failed:\n  - " + "\n  - ".join(problems)
            )
