"""Declared experiments, frozen before they are run.

Charter 14.2 and requirements EXP-001 and EXP-002.

EXP-001 lists six things every experiment must declare -- controls, manipulated
variables, measured outputs, expected result, failure threshold, and analysis --
and adds the words that make the requirement bite: *before final evaluation*.
A declaration that can be edited after the numbers are in is not a declaration.

So the protocol is frozen: ``freeze()`` records a hash of the whole declaration,
and ``evaluate()`` refuses to run against a protocol that was never frozen or
whose declaration has changed since. The refusal is an
``InvalidForScientificInferenceError`` rather than a warning, because a result
compared against a threshold chosen after seeing it is not a weaker result but a
different kind of thing.

EXP-002 asks for deviations to be recorded. ``record_deviation`` appends to the
outcome rather than amending the protocol, so the frozen declaration and what
actually happened stay separately readable.
"""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field as dc_field
from typing import Any

from mfespf.core.errors import ConfigurationError, InvalidForScientificInferenceError
from mfespf.core.hashing import digest

__all__ = ["ExperimentOutcome", "ExperimentProtocol"]


@dataclass
class ExperimentProtocol:
    """The EXP-001 declaration, made before the run and hashed."""

    experiment_id: str
    question: str
    control_configuration: str
    manipulated_variables: tuple[str, ...]
    measured_outputs: tuple[str, ...]
    expected_result: str
    failure_threshold: str
    analysis_method: str
    forcing_regime: str = "none"
    assumptions: tuple[str, ...] = ()
    notes: str = ""
    _frozen_hash: str | None = dc_field(default=None, init=False, repr=False)

    def __post_init__(self) -> None:
        missing = [
            name
            for name in (
                "question",
                "control_configuration",
                "expected_result",
                "failure_threshold",
                "analysis_method",
            )
            if not str(getattr(self, name)).strip()
        ]
        if not self.manipulated_variables:
            missing.append("manipulated_variables")
        if not self.measured_outputs:
            missing.append("measured_outputs")
        if missing:
            raise ConfigurationError(
                f"experiment {self.experiment_id!r} declares no {', '.join(sorted(missing))}; "
                f"EXP-001 requires every field before the experiment is run"
            )

    # --- freezing ---------------------------------------------------------
    def declaration(self) -> dict[str, Any]:
        """Everything the hash covers: the declaration and nothing else."""
        return {
            "experiment_id": self.experiment_id,
            "question": self.question,
            "control_configuration": self.control_configuration,
            "manipulated_variables": list(self.manipulated_variables),
            "measured_outputs": list(self.measured_outputs),
            "expected_result": self.expected_result,
            "failure_threshold": self.failure_threshold,
            "analysis_method": self.analysis_method,
            "forcing_regime": self.forcing_regime,
            "assumptions": list(self.assumptions),
            "notes": self.notes,
        }

    def freeze(self) -> str:
        """Record the declaration hash. Called before the run, not after."""
        self._frozen_hash = digest(self.declaration(), length=32)
        return self._frozen_hash

    @property
    def frozen_hash(self) -> str | None:
        return self._frozen_hash

    @property
    def is_frozen(self) -> bool:
        return self._frozen_hash is not None

    def assert_unchanged(self) -> None:
        """Refuse to proceed if the declaration was never frozen or has moved."""
        if self._frozen_hash is None:
            raise InvalidForScientificInferenceError(
                f"experiment {self.experiment_id!r} was never frozen; EXP-001 requires the "
                f"expected result, failure threshold and analysis to be declared before "
                f"evaluation, and an unfrozen protocol cannot show that they were"
            )
        current = digest(self.declaration(), length=32)
        if current != self._frozen_hash:
            raise InvalidForScientificInferenceError(
                f"experiment {self.experiment_id!r} was edited after it was frozen "
                f"({self._frozen_hash} -> {current}); a threshold chosen after seeing the "
                f"result is not a test of the prediction (EXP-001)"
            )

    # --- evaluation -------------------------------------------------------
    def evaluate(
        self, *, passed: bool, measured: dict[str, float], summary: str
    ) -> ExperimentOutcome:
        """Compare against the frozen declaration and return the outcome.

        The caller supplies the verdict because the threshold is expressed in
        physical terms this class cannot interpret. What this class enforces is
        that the threshold existed, unedited, before the verdict was formed.
        """
        self.assert_unchanged()
        unmeasured = [name for name in self.measured_outputs if name not in measured]
        outcome = ExperimentOutcome(
            experiment_id=self.experiment_id,
            protocol_hash=self._frozen_hash or "",
            passed=passed,
            measured=dict(measured),
            summary=summary,
        )
        if unmeasured:
            outcome.record_deviation(
                f"declared outputs not measured: {', '.join(sorted(unmeasured))}"
            )
        return outcome

    def describe(self) -> str:
        lines = [
            f"experiment {self.experiment_id}: {self.question}",
            f"  control:      {self.control_configuration}",
            f"  manipulated:  {', '.join(self.manipulated_variables)}",
            f"  measured:     {', '.join(self.measured_outputs)}",
            f"  forcing:      {self.forcing_regime}",
            f"  expected:     {self.expected_result}",
            f"  fails if:     {self.failure_threshold}",
            f"  analysis:     {self.analysis_method}",
        ]
        if self.assumptions:
            lines.append(f"  assumes:      {'; '.join(self.assumptions)}")
        lines.append(
            f"  declaration hash: {self._frozen_hash}"
            if self._frozen_hash
            else "  declaration hash: NOT FROZEN"
        )
        return "\n".join(lines)


@dataclass
class ExperimentOutcome:
    """What the experiment produced, kept separate from what it declared."""

    experiment_id: str
    protocol_hash: str
    passed: bool
    measured: dict[str, float] = dc_field(default_factory=dict)
    summary: str = ""
    deviations: list[str] = dc_field(default_factory=list)

    def record_deviation(self, description: str) -> None:
        """EXP-002. Appended here, never folded back into the declaration."""
        self.deviations.append(description)

    @property
    def as_declared(self) -> bool:
        """True when the experiment ran as written, with nothing to explain."""
        return not self.deviations

    def describe(self) -> str:
        verdict = "PASSED" if self.passed else "FAILED"
        lines = [f"experiment {self.experiment_id}: {verdict} ({self.protocol_hash})"]
        if self.summary:
            lines.append(f"  {self.summary}")
        for name, value in sorted(self.measured.items()):
            lines.append(f"  {name} = {value:.6g}")
        for deviation in self.deviations:
            lines.append(f"  DEVIATION: {deviation}")
        return "\n".join(lines)
