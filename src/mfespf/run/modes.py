"""Operating modes, and the language each one is entitled to.

Charter 11 and 14; requirements TIM-005, TIM-006, VAL-007.

TIM-005 asks the framework to support seven operating modes. TIM-006 adds the
rule that makes the list worth keeping: *reports shall not describe one
operating mode using the evidentiary language of another*. VAL-007 names the
case that matters most -- a historical reconstruction reported as a prospective
prediction -- and it is the easiest mistake in the field to make honestly,
because the arithmetic is identical and only the information available at the
time differs.

The modes are therefore ordered by how much information about the compared
period was available to the model when it ran. A run may always be described as
a mode of equal or *lower* standing, since claiming less than you are entitled
to is not a misrepresentation. Describing it as a higher one is refused.

The ordering is not about quality. A controlled experiment is a better thing to
run than a bad hindcast; what the ordering ranks is the strength of the
evidentiary claim the mode licenses.
"""

from __future__ import annotations

from enum import StrEnum

from mfespf.core.errors import InvalidForScientificInferenceError
from mfespf.registries.claims import ClaimStatus

__all__ = ["OperatingMode"]


class OperatingMode(StrEnum):
    """TIM-005.  What kind of run this is, in the sense evidence cares about."""

    EQUILIBRIUM = "equilibrium"
    CONTROLLED_EXPERIMENT = "controlled_experiment"
    PRESCRIBED_HISTORY = "prescribed_history"
    FREE_EVOLUTION = "free_evolution"
    HISTORICAL_RECONSTRUCTION = "historical_reconstruction"
    HINDCAST = "hindcast"
    OUT_OF_SAMPLE_PREDICTION = "out_of_sample_prediction"

    @property
    def standing(self) -> int:
        """How strong a claim the mode licenses, on an ordinal scale.

        Prescribed history is deliberately the lowest of the comparison modes.
        A run whose inputs were taken from the period being compared against has
        been told part of the answer, and agreement is correspondingly cheap.
        """
        return {
            OperatingMode.EQUILIBRIUM: 0,
            OperatingMode.CONTROLLED_EXPERIMENT: 0,
            OperatingMode.FREE_EVOLUTION: 0,
            OperatingMode.PRESCRIBED_HISTORY: 1,
            OperatingMode.HISTORICAL_RECONSTRUCTION: 2,
            OperatingMode.HINDCAST: 3,
            OperatingMode.OUT_OF_SAMPLE_PREDICTION: 4,
        }[self]

    @property
    def compares_against_a_record(self) -> bool:
        """Whether the mode involves comparison with a recorded period at all."""
        return self.standing > 0

    @property
    def is_prospective(self) -> bool:
        """Whether the compared values were unavailable when the model was fixed."""
        return self in (OperatingMode.HINDCAST, OperatingMode.OUT_OF_SAMPLE_PREDICTION)

    @property
    def uses_target_period_inputs(self) -> bool:
        """Whether inputs from the compared period were supplied to the run.

        The distinguishing fact between a reconstruction and a prediction, and
        the reason the first cannot be reported as the second: a reconstruction
        was given the boundary conditions of the period it reproduces.
        """
        return self in (
            OperatingMode.PRESCRIBED_HISTORY,
            OperatingMode.HISTORICAL_RECONSTRUCTION,
        )

    def assert_may_be_reported_as(self, claimed: OperatingMode) -> None:
        """TIM-006 and VAL-007, as a check rather than as a convention.

        Refuses when a run would be described in the language of a mode that
        licenses a stronger claim than the run actually supports.
        """
        if claimed.standing <= self.standing:
            return
        detail = ""
        if self.uses_target_period_inputs and claimed.is_prospective:
            detail = (
                " The run was given inputs from the period it is being compared against, so "
                "its agreement is not evidence of prospective skill (VAL-007)."
            )
        raise InvalidForScientificInferenceError(
            f"a run executed in {self.value!r} mode may not be reported as {claimed.value!r}: "
            f"that would describe it with the evidentiary language of a mode it does not "
            f"support (TIM-006).{detail}"
        )

    def status_for_result(self, passed: bool) -> ClaimStatus:
        """What a pass or a fail in this mode is worth (CLM-002, VAL-007).

        The asymmetry is the point. A failure is decisive in any mode: a
        mechanism given every advantage the experiment could offer and still
        unable to account for the observations has been refuted, and being handed
        the boundary conditions of the period only makes that stronger.

        A pass is worth what the mode's standing says. Reproducing observations
        while holding the target period's own inputs is a qualified success, not
        a supported claim -- the run was told part of the answer. Only a
        prospective comparison earns SUPPORTED_WITHIN_SCOPE, and a mode that
        compared against nothing earns nothing either way.
        """
        if not self.compares_against_a_record:
            return ClaimStatus.EMPIRICALLY_UNRESOLVED
        if not passed:
            return ClaimStatus.REFUTED
        if self.uses_target_period_inputs:
            return ClaimStatus.QUALIFIED
        return ClaimStatus.SUPPORTED_WITHIN_SCOPE

    def describe(self) -> str:
        parts = [f"operating mode: {self.value}"]
        if self.uses_target_period_inputs:
            parts.append("inputs taken from the compared period")
        if self.is_prospective:
            parts.append("compared values unavailable when the model was fixed")
        return "; ".join(parts)
