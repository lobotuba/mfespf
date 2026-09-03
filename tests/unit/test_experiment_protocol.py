"""Declared experiments: EXP-001 and EXP-002.

EXP-001's six required declarations are the easy half. The half that decides
whether the requirement means anything is the phrase *before final evaluation*,
and that is what most of this file is about: a threshold that can be edited after
the numbers arrive is not a threshold.
"""

from __future__ import annotations

import pytest

from mfespf.core.errors import ConfigurationError, InvalidForScientificInferenceError
from mfespf.experiments.protocol import ExperimentProtocol

pytestmark = pytest.mark.stage(5)


def _protocol(**overrides) -> ExperimentProtocol:
    fields = dict(
        experiment_id="EXP-STAGE5-001",
        question="Does adding a deep reservoir produce a second response timescale?",
        control_configuration="surface_thermal_storage, identical but for the reservoir",
        manipulated_variables=("deep_ocean",),
        measured_outputs=("fast_timescale_s", "slow_timescale_s"),
        expected_result="two separated modes in the two-reservoir case, one in the control",
        failure_threshold="a two-mode fit that does not beat a one-mode fit by 10x",
        analysis_method="two-mode exponential fit to the step response",
        forcing_regime="step",
    )
    fields.update(overrides)
    return ExperimentProtocol(**fields)  # type: ignore[arg-type]


# --- the declaration -------------------------------------------------------


@pytest.mark.req("EXP-001")
def test_every_declared_field_is_required():
    for missing in ("expected_result", "failure_threshold", "analysis_method", "question"):
        with pytest.raises(ConfigurationError, match=missing):
            _protocol(**{missing: "  "})
    with pytest.raises(ConfigurationError, match="measured_outputs"):
        _protocol(measured_outputs=())
    with pytest.raises(ConfigurationError, match="manipulated_variables"):
        _protocol(manipulated_variables=())


@pytest.mark.req("EXP-001", "AUD-002")
def test_freezing_records_a_hash_of_the_whole_declaration():
    protocol = _protocol()
    assert not protocol.is_frozen
    first = protocol.freeze()
    assert protocol.is_frozen
    assert protocol.frozen_hash == first
    # A second protocol with the same declaration hashes the same; one with a
    # different threshold does not.
    assert _protocol().freeze() == first
    assert _protocol(failure_threshold="anything at all").freeze() != first


# --- before final evaluation -----------------------------------------------


@pytest.mark.req("EXP-001", "INF-006")
def test_an_unfrozen_protocol_cannot_be_evaluated():
    """Nothing shows that the threshold predated the result, so nothing may."""
    with pytest.raises(InvalidForScientificInferenceError, match="never frozen"):
        _protocol().evaluate(passed=True, measured={}, summary="")


@pytest.mark.req("EXP-001", "INF-006")
def test_a_threshold_edited_after_freezing_invalidates_the_evaluation():
    """The failure this requirement exists to catch, made mechanical.

    Loosening a threshold once the numbers are in converts a test of a
    prediction into a description of an outcome. The two are not different in
    degree.
    """
    protocol = _protocol()
    protocol.freeze()
    protocol.failure_threshold = "a two-mode fit that does not beat a one-mode fit by 1.1x"
    with pytest.raises(InvalidForScientificInferenceError, match="after it was frozen"):
        protocol.evaluate(passed=True, measured={}, summary="")


@pytest.mark.req("EXP-001")
def test_an_unedited_protocol_evaluates_normally():
    protocol = _protocol()
    protocol.freeze()
    outcome = protocol.evaluate(
        passed=True,
        measured={"fast_timescale_s": 4.5e7, "slow_timescale_s": 8.3e9},
        summary="two modes, separated by a factor of 185",
    )
    assert outcome.passed
    assert outcome.as_declared
    assert outcome.protocol_hash == protocol.frozen_hash
    assert "PASSED" in outcome.describe()


# --- deviations ------------------------------------------------------------


@pytest.mark.req("EXP-002")
def test_an_output_declared_but_not_measured_is_recorded_as_a_deviation():
    """Silently dropping a declared output would narrow the experiment.

    The run still produces an outcome -- refusing would be worse, since the
    measurements that were taken are still informative -- but it is no longer
    the experiment that was declared, and it says so.
    """
    protocol = _protocol()
    protocol.freeze()
    outcome = protocol.evaluate(
        passed=True, measured={"fast_timescale_s": 4.5e7}, summary="only the fast mode resolved"
    )
    assert not outcome.as_declared
    assert any("slow_timescale_s" in d for d in outcome.deviations)
    assert "DEVIATION" in outcome.describe()


@pytest.mark.req("EXP-002")
def test_deviations_are_appended_to_the_outcome_not_the_declaration():
    """What was declared and what happened stay separately readable."""
    protocol = _protocol()
    frozen = protocol.freeze()
    outcome = protocol.evaluate(
        passed=False,
        measured={"fast_timescale_s": 4.5e7, "slow_timescale_s": 8.3e9},
        summary="",
    )
    outcome.record_deviation("the run was stopped early after a power failure")
    assert protocol.frozen_hash == frozen  # the declaration is untouched
    protocol.assert_unchanged()
    assert len(outcome.deviations) == 1
    assert "FAILED" in outcome.describe()


@pytest.mark.req("EXP-001")
def test_the_protocol_describes_itself_for_a_reader():
    protocol = _protocol()
    assert "NOT FROZEN" in protocol.describe()
    protocol.freeze()
    text = protocol.describe()
    assert "fails if:" in text
    assert protocol.frozen_hash in text
