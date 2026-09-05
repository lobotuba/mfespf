"""Hypothesis registry.

Charter section 14 and requirements HYP-001 through HYP-010.

A claim is a statement that can fail. A *hypothesis* is the mechanism proposed to
explain something, and the charter asks more of it than of a claim: not only what
would refute it, but where it enters the model, what it does to the conserved
quantities, what the control configuration is, and which competing explanations
the proposed experiment can and cannot tell it apart from.

That last one (HYP-008) is the field most likely to be left vague, and the one
most worth insisting on. An experiment that cannot distinguish a mechanism from
its alternatives can still support it -- and a record that never said so lets a
reader mistake "consistent with" for "evidence for".

HYP-010 is enforced by omission rather than by a check: ``mechanism_class`` is a
label with no effect on validation, so a mechanism declared unconventional
passes exactly the tests a conventional one does. There is no lenient path and no
strict one.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import StrEnum

from mfespf.core.errors import ConfigurationError, OwnershipError
from mfespf.core.hashing import digest
from mfespf.core.versioning import Version
from mfespf.registries.base import Registry
from mfespf.registries.claims import ExpectedSignature

__all__ = [
    "EntryPoint",
    "HypothesisRecord",
    "HypothesisRegistry",
    "MechanismClass",
]


class EntryPoint(StrEnum):
    """HYP-005.  The interfaces through which a mechanism may enter the model.

    The list is the charter's and is deliberately closed. A mechanism that fits
    none of these is not a mechanism the framework can test, because there is
    nowhere for it to act that the accounting can see.
    """

    MODULE = "module"
    TENDENCY = "tendency"
    BOUNDARY = "boundary"
    SOURCE = "source"
    SINK = "sink"
    CONSTITUTIVE_RELATION = "constitutive_relation"
    DIAGNOSTIC = "diagnostic"


class MechanismClass(StrEnum):
    """A label, and nothing more (HYP-010).

    It exists so a reader can see what kind of mechanism is being tested. It has
    no effect on validation, on the thresholds, or on what the evidence has to
    show: a conventional mechanism and an unconventional one face identical
    requirements, and the way to guarantee that is to give the field nothing to
    switch on.
    """

    CONVENTIONAL = "conventional"
    UNCONVENTIONAL = "unconventional"


@dataclass(frozen=True)
class HypothesisRecord:
    """Charter 14 / HYP-001 through HYP-008, with every field required."""

    record_id: str
    version: Version
    statement: str
    claim_ids: tuple[str, ...]
    physical_basis: str
    governing_equation_ids: tuple[str, ...]
    required_inputs: tuple[str, ...]
    produced_outputs: tuple[str, ...]
    expected: ExpectedSignature
    validity_domain: str
    failure_criterion: str

    conserved_quantities: tuple[str, ...]
    transfer_accounting: str

    entry_point: EntryPoint
    entry_point_detail: str
    implementing_modules: tuple[str, ...]

    null_configuration: str
    distinguishes_from: tuple[str, ...]
    cannot_distinguish_from: tuple[str, ...]

    mechanism_class: MechanismClass = MechanismClass.CONVENTIONAL
    assumptions: tuple[str, ...] = ()
    notes: str = ""
    frozen_criteria_hash: str | None = None

    def __post_init__(self) -> None:
        problems: list[str] = []
        if not self.claim_ids:
            problems.append("names no claim_id (HYP-001)")
        for name in ("statement", "physical_basis", "validity_domain", "failure_criterion"):
            if not str(getattr(self, name)).strip():
                problems.append(f"declares no {name} (HYP-002)")
        if not self.governing_equation_ids:
            problems.append("names no governing equations (HYP-002, EQR-005)")
        if not self.required_inputs or not self.produced_outputs:
            problems.append("declares no inputs or no outputs (HYP-002)")

        # HYP-003. "Where applicable" applies to the spectral entry only; the
        # other four are always applicable, and an unspecified sign is the
        # commonest way a prediction ends up unfalsifiable.
        for aspect in ("sign", "magnitude_bound", "temporal", "spatial"):
            if getattr(self.expected, aspect) == "unspecified":
                problems.append(f"leaves the expected {aspect} unspecified (HYP-003)")

        if not self.conserved_quantities:
            problems.append(
                "declares no affected conserved quantities; state 'none' explicitly if the "
                "mechanism moves nothing (HYP-004)"
            )
        if not self.transfer_accounting.strip():
            problems.append("does not say how its transfers are accounted (HYP-004)")
        if not self.entry_point_detail.strip():
            problems.append("does not say where it enters the model (HYP-005)")
        if not self.implementing_modules:
            problems.append("names no implementing module (HYP-005, HYP-006)")
        if not self.null_configuration.strip():
            problems.append("declares no null or control configuration (HYP-007)")
        if not self.distinguishes_from or not self.cannot_distinguish_from:
            problems.append(
                "does not say which competing explanations the experiment can and cannot "
                "distinguish it from; an experiment that separates nothing should say so "
                "explicitly (HYP-008)"
            )
        if problems:
            raise ConfigurationError(
                f"hypothesis {self.record_id!r} is not registrable: " + "; ".join(problems)
            )

    # --- freezing ---------------------------------------------------------
    def criteria_hash(self) -> str:
        """Exactly what HYP-009 forbids changing after the freeze."""
        return digest(
            {
                "failure_criterion": self.failure_criterion,
                "expected": self.expected,
                "validity_domain": self.validity_domain,
                "null_configuration": self.null_configuration,
                "cannot_distinguish_from": sorted(self.cannot_distinguish_from),
            },
            length=16,
        )

    def freeze(self) -> HypothesisRecord:
        return replace(self, frozen_criteria_hash=self.criteria_hash())

    @property
    def is_frozen(self) -> bool:
        return self.frozen_criteria_hash is not None

    def assert_criteria_intact(self) -> None:
        if not self.is_frozen:
            raise ConfigurationError(
                f"hypothesis {self.record_id!r} has not been frozen; HYP-009 requires the "
                f"failure threshold to be fixed before the withheld evaluation"
            )
        if self.frozen_criteria_hash != self.criteria_hash():
            raise ConfigurationError(
                f"hypothesis {self.record_id!r} changed after freezing "
                f"({self.frozen_criteria_hash} -> {self.criteria_hash()}); register a new "
                f"version that supersedes this one rather than editing it (HYP-009)"
            )

    def describe(self) -> str:
        lines = [
            f"hypothesis {self.record_id}@{self.version} [{self.mechanism_class.value}]",
            f"  {self.statement}",
            f"  claims:        {', '.join(self.claim_ids)}",
            f"  enters as:     {self.entry_point.value} -- {self.entry_point_detail}",
            f"  implemented by:{' ' + ', '.join(self.implementing_modules)}",
            f"  conserves:     {', '.join(self.conserved_quantities)}",
            f"  control:       {self.null_configuration}",
            f"  expected:      sign {self.expected.sign}, {self.expected.magnitude_bound}",
            f"  fails if:      {self.failure_criterion}",
            f"  separates it from:     {'; '.join(self.distinguishes_from)}",
            f"  does NOT separate it from: {'; '.join(self.cannot_distinguish_from)}",
        ]
        if self.frozen_criteria_hash:
            lines.append(f"  frozen:        {self.frozen_criteria_hash}")
        return "\n".join(lines)


class HypothesisRegistry(Registry[HypothesisRecord]):
    """Registered mechanisms, and the checks that are cheaper here than later."""

    def __init__(self) -> None:
        super().__init__("hypothesis")

    def assert_modules_present(self, record_id: str, module_ids: set[str]) -> None:
        """HYP-005: the mechanism must act through modules the run actually has.

        Presence, not ownership. A prescribed boundary module owns no state and
        is still where a boundary-entry mechanism enters, so requiring every
        implementing module to own something would refuse the commonest honest
        case.

        HYP-006 -- that an implementation must not write unowned state -- is not
        checked here and does not need to be. There is no method on ``Module``
        that writes foreign state, and every tendency contribution passes
        ``StateRegistry.assert_may_contribute`` before the coupler accepts it.
        The guarantee is structural rather than declared, which is stronger.
        """
        record = self.get(record_id)
        missing = sorted(set(record.implementing_modules) - module_ids)
        if missing:
            raise OwnershipError(
                f"hypothesis {record_id!r} names {missing} as implementing modules, but the "
                f"configuration has {sorted(module_ids)}; a mechanism that acts through a "
                f"module the run does not contain cannot be tested by it (HYP-005)"
            )

    def for_claim(self, claim_id: str) -> list[HypothesisRecord]:
        return [h for h in self.current() if claim_id in h.claim_ids]

    def undistinguished(self, record_id: str) -> tuple[str, ...]:
        """The competing explanations this experiment cannot rule out (HYP-008).

        Returned as data rather than prose so a report can be required to carry
        it. A mechanism supported by an experiment that could not separate it
        from three alternatives is supported jointly with those alternatives,
        and the list is what makes that statable.
        """
        return self.get(record_id).cannot_distinguish_from
