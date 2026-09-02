"""Claim and experiment registries.

Charter section 14 and requirements CLM-001 through CLM-003, EXP-001, EXP-002.

The registry enforces one ordering rule that no amount of care enforces on its
own: failure and support criteria are frozen before the withheld evaluation runs
(HYP-009).  After freezing, changing a threshold produces a new claim version
with a link to the old one, so the record shows the change instead of absorbing
it.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import StrEnum

from mfespf.core.errors import ConfigurationError, ProvenanceError
from mfespf.core.hashing import digest
from mfespf.core.versioning import Version
from mfespf.registries.base import Registry

__all__ = [
    "ClaimRecord",
    "ClaimRegistry",
    "ClaimStatus",
    "ExpectedSignature",
    "ExperimentRecord",
    "ExperimentRegistry",
]


class ClaimStatus(StrEnum):
    """Charter 14.3 / CLM-002."""

    PROPOSED = "proposed"
    AWAITING_IMPLEMENTATION = "awaiting_implementation"
    TESTING = "testing"
    SUPPORTED_WITHIN_SCOPE = "supported_within_scope"
    QUALIFIED = "qualified"
    REFUTED = "refuted"
    MATHEMATICALLY_UNRESOLVED = "mathematically_unresolved"
    NUMERICALLY_UNRESOLVED = "numerically_unresolved"
    EMPIRICALLY_UNRESOLVED = "empirically_unresolved"
    UNABLE_TO_DISTINGUISH = "unable_to_distinguish"
    SUPERSEDED = "superseded"

    @property
    def is_terminal(self) -> bool:
        return self not in (
            ClaimStatus.PROPOSED,
            ClaimStatus.AWAITING_IMPLEMENTATION,
            ClaimStatus.TESTING,
        )


@dataclass(frozen=True, slots=True)
class ExpectedSignature:
    """HYP-003: sign, magnitude, and the temporal/spatial/spectral behaviour.

    Declared before evaluation.  A mechanism that produces the right global mean
    with the wrong spectral signature has failed 15.4, and that is only checkable
    if the expected signature was written down first.
    """

    sign: str = "unspecified"
    magnitude_bound: str = "unspecified"
    temporal: str = "unspecified"
    spatial: str = "unspecified"
    spectral: str = "not_applicable"


@dataclass(frozen=True)
class ClaimRecord:
    """Charter 14.1.  The fields the charter lists as mandatory are required here."""

    record_id: str
    version: Version
    statement: str
    mechanism: str
    physical_basis: str
    governing_equation_ids: tuple[str, ...]
    failure_criterion: str
    support_criterion: str
    assumptions: tuple[str, ...] = ()
    validity_domain: str = ""
    active_modules: tuple[str, ...] = ()
    inactive_modules: tuple[str, ...] = ()
    prescribed_quantities: tuple[str, ...] = ()
    independent_variables: tuple[str, ...] = ()
    dependent_variables: tuple[str, ...] = ()
    expected: ExpectedSignature = ExpectedSignature()
    calibration_data: tuple[str, ...] = ()
    validation_data: tuple[str, ...] = ()
    withheld_data: tuple[str, ...] = ()
    uncertainty_method: str = "unresolved"
    control_experiment: str = ""
    null_hypothesis: str = ""
    alternative_hypotheses: tuple[str, ...] = ()
    numerical_correctness_gates: tuple[str, ...] = ()
    reproduction_configuration: str = ""
    status: ClaimStatus = ClaimStatus.PROPOSED
    result: str = ""
    limitations: tuple[str, ...] = ()
    frozen_criteria_hash: str | None = None
    supersedes: str | None = None

    def __post_init__(self) -> None:
        if not self.failure_criterion.strip():
            raise ConfigurationError(
                f"claim {self.record_id!r} declares no failure criterion; an unfalsifiable claim "
                f"cannot be registered (charter 4.32)"
            )
        if not self.governing_equation_ids:
            raise ConfigurationError(
                f"claim {self.record_id!r} names no governing equations (CLM-001, EQR-005)"
            )
        overlap = set(self.calibration_data) & set(self.withheld_data)
        if overlap:
            raise ConfigurationError(
                f"claim {self.record_id!r} lists {sorted(overlap)} as both calibration and "
                f"withheld data (INF-003)"
            )

    @property
    def claim_id(self) -> str:
        return self.record_id

    def criteria_hash(self) -> str:
        """Hash over exactly the things HYP-009 forbids changing after freeze."""
        return digest(
            {
                "failure_criterion": self.failure_criterion,
                "support_criterion": self.support_criterion,
                "expected": self.expected,
                "uncertainty_method": self.uncertainty_method,
                "withheld_data": sorted(self.withheld_data),
                "numerical_correctness_gates": sorted(self.numerical_correctness_gates),
            },
            length=16,
        )

    def freeze(self) -> ClaimRecord:
        """Return a copy with its criteria hash pinned (HYP-009, A4.5)."""
        return replace(self, frozen_criteria_hash=self.criteria_hash())

    @property
    def is_frozen(self) -> bool:
        return self.frozen_criteria_hash is not None

    def assert_criteria_intact(self) -> None:
        """Verify that nothing meaning-bearing changed since the freeze."""
        if not self.is_frozen:
            raise ProvenanceError(
                f"claim {self.record_id!r} has not been frozen; a withheld evaluation requires "
                f"frozen failure thresholds (HYP-009)"
            )
        if self.criteria_hash() != self.frozen_criteria_hash:
            raise ProvenanceError(
                f"claim {self.record_id!r} criteria changed after freeze "
                f"({self.frozen_criteria_hash} -> {self.criteria_hash()}). Register a new claim "
                f"version instead (HYP-009)."
            )

    def conclude(
        self, status: ClaimStatus, result: str, limitations: tuple[str, ...] = ()
    ) -> ClaimRecord:
        """Attach an outcome.  CLM-003 bars an unqualified confirmation status."""
        if not status.is_terminal:
            raise ConfigurationError(f"{status.value} is not an outcome status")
        return replace(self, status=status, result=result, limitations=limitations)


@dataclass(frozen=True)
class ExperimentRecord:
    """Charter 14.2 / EXP-001, EXP-002."""

    record_id: str
    version: Version
    purpose: str
    related_claim_ids: tuple[str, ...]
    controlled_variables: tuple[str, ...] = ()
    manipulated_variables: tuple[str, ...] = ()
    measured_outputs: tuple[str, ...] = ()
    predeclared_analysis: str = ""
    expected_result: str = ""
    failure_thresholds: tuple[str, ...] = ()
    stopping_criteria: tuple[str, ...] = ()
    run_ids: tuple[str, ...] = ()
    deviations: tuple[str, ...] = ()
    outcome: str = ""

    def __post_init__(self) -> None:
        if not self.related_claim_ids:
            raise ConfigurationError(f"experiment {self.record_id!r} references no claim (CLM-001)")
        if not self.predeclared_analysis.strip():
            raise ConfigurationError(
                f"experiment {self.record_id!r} declares no analysis in advance (EXP-001)"
            )

    @property
    def experiment_id(self) -> str:
        return self.record_id


class ClaimRegistry(Registry[ClaimRecord]):
    def __init__(self) -> None:
        super().__init__("claim")

    def open_claims(self) -> list[ClaimRecord]:
        return [c for c in self.current() if not c.status.is_terminal]

    def by_status(self, status: ClaimStatus) -> list[ClaimRecord]:
        return [c for c in self.current() if c.status is status]


class ExperimentRegistry(Registry[ExperimentRecord]):
    def __init__(self) -> None:
        super().__init__("experiment")
