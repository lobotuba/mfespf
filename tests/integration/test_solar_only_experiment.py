"""The solar-only attribution experiment, end to end.

Requirements HYP-001 through HYP-010, EXP-001, EXP-002, INF-003, CLM-002,
VAL-003, VAL-005, VAL-007, TIM-005, TIM-006, CON-009.

The framework's first experiment against a real observational record, and one
designed to fail: the kernel has no carbon dioxide, no aerosols and no volcanic
forcing, so an irradiance record alone cannot account for the observed warming.

These tests check that it fails *for the right reasons and with the right
paperwork* -- the threshold taken from the record rather than the model, the
control producing exactly nothing, the status set by what the operating mode
licenses, and the list of alternatives the experiment could not separate carried
into the result.

Skipped when the records are absent; `python tools/fetch_data.py` fetches them.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.data.loaders import available, load_declared
from mfespf.data.records import ProvenanceClass
from mfespf.experiments.solar_only import (
    BASELINE_END,
    BASELINE_START,
    CLAIM_ID,
    compared_period,
    observational_threshold,
    run_solar_only_experiment,
    solar_only_hypothesis,
)
from mfespf.registries.claims import ClaimStatus
from mfespf.registries.datasets import DatasetRole
from mfespf.run.modes import OperatingMode

pytestmark = pytest.mark.stage(6)

present = available()
needs_data = pytest.mark.skipif(
    not {"NRLTSI2-MONTHLY", "HADCRUT5-GLOBAL-MONTHLY"} <= set(present),
    reason="the declared records are not in data/raw; run `python tools/fetch_data.py`",
)


@pytest.fixture(scope="module")
def experiment():
    return run_solar_only_experiment()


# --- the declaration, before anything runs ----------------------------------


@pytest.mark.req("HYP-002", "HYP-003", "HYP-007", "HYP-008")
def test_the_hypothesis_is_complete_and_frozen_without_touching_any_data():
    """It has to be registrable before the run, or the run shapes it."""
    hypothesis = solar_only_hypothesis()
    assert hypothesis.is_frozen
    hypothesis.assert_criteria_intact()
    assert hypothesis.claim_ids == (CLAIM_ID,)
    assert hypothesis.expected.sign.startswith("positive")
    assert "constant irradiance" in hypothesis.null_configuration
    assert len(hypothesis.cannot_distinguish_from) >= 4
    assert any("greenhouse" in a for a in hypothesis.cannot_distinguish_from)


@needs_data
@pytest.mark.req("EXP-001", "INF-006")
def test_the_threshold_comes_from_the_record_and_not_from_the_model():
    """The number has to come from somewhere that is not the model's error.

    It is the record's own mean 95 percent half-width over the compared period,
    and the compared period is the overlap of the two records' coverage -- a
    property of the data, knowable before anything runs.
    """
    irradiance = load_declared("NRLTSI2-MONTHLY").series("total_solar_irradiance")
    observed = next(iter(load_declared("HADCRUT5-GLOBAL-MONTHLY")._series.values()))
    period = compared_period(irradiance, observed)
    threshold = observational_threshold(observed, over=period)

    inside = np.logical_and(observed.times >= period[0], observed.times <= period[1])
    assert threshold == pytest.approx(float(np.mean(observed.uncertainty[inside])))
    # Not the whole record: the early years are more uncertain and the
    # comparison never reaches them.
    assert threshold != pytest.approx(float(np.mean(observed.uncertainty)))
    assert 0.05 < threshold < 0.3


# --- the control ------------------------------------------------------------


@needs_data
@pytest.mark.req("HYP-007", "VAL-001", "MOD-014")
def test_the_constant_irradiance_control_produces_exactly_nothing(experiment):
    """The strongest statement a null control can make.

    Not "a small range" -- zero range, and zero variance explained. Anything
    else would mean something other than the solar variation was moving the
    modelled temperature, and every number in the forced run would carry that
    contamination too. It caught a real one: the run had been starting away from
    equilibrium, and its first decade of spin-up was appearing as a modelled
    range eight times larger than the solar response.
    """
    control = experiment.control_comparison
    assert control.modelled_range == pytest.approx(0.0, abs=1e-12)
    assert control.explained_fraction() == pytest.approx(0.0, abs=1e-12)
    assert control.n == experiment.comparison.n


@needs_data
@pytest.mark.req("HYP-007", "CLM-002")
def test_the_solar_variation_does_something_and_it_is_far_too_little(experiment):
    """Both halves matter.

    The mechanism is not inert -- it explains a few percent more than the
    control, so the experiment did test something. And a few percent is not the
    hundred the hypothesis needs.
    """
    forced = experiment.comparison
    control = experiment.control_comparison

    assert forced.modelled_range > 0.0
    assert forced.explained_fraction() > control.explained_fraction()
    assert forced.explained_fraction() < 0.2
    # Two orders of magnitude between what the mechanism produces and what
    # happened.
    assert forced.observed_range > 20.0 * forced.modelled_range


# --- the verdict ------------------------------------------------------------


@needs_data
@pytest.mark.req("CLM-002", "EXP-001", "VAL-007")
def test_the_claim_is_refuted_against_a_threshold_declared_beforehand(experiment):
    assert experiment.comparison.rms_error > experiment.threshold
    assert not experiment.outcome.passed
    assert experiment.claim_status is ClaimStatus.REFUTED
    assert experiment.refuted
    # The protocol was frozen before the run and is unedited, which is what
    # makes the verdict a test of a prediction rather than a description.
    experiment.protocol.assert_unchanged()
    assert experiment.outcome.protocol_hash == experiment.protocol.frozen_hash
    assert experiment.outcome.as_declared


@needs_data
@pytest.mark.req("CLM-002", "VAL-007", "TIM-006")
def test_a_refutation_from_a_reconstruction_is_decisive_where_a_pass_would_not_be(
    experiment,
):
    """The asymmetry that makes this experiment worth running at all.

    The run was handed the irradiance of the period it is judged on. Succeeding
    under those conditions would have been a qualified result at best. Failing
    under them is not qualified: the mechanism was given every advantage the
    experiment could offer.
    """
    mode = experiment.operating_mode
    assert mode is OperatingMode.PRESCRIBED_HISTORY
    assert mode.uses_target_period_inputs
    assert mode.status_for_result(False) is ClaimStatus.REFUTED
    assert mode.status_for_result(True) is ClaimStatus.QUALIFIED
    assert mode.status_for_result(True) is not ClaimStatus.SUPPORTED_WITHIN_SCOPE


@needs_data
@pytest.mark.req("HYP-008", "VAL-005")
def test_the_result_carries_what_the_experiment_could_not_distinguish(experiment):
    """A refutation of *this configuration*, not of solar physics.

    The experiment observes one global-mean scalar. It cannot separate the
    mechanism from greenhouse, aerosol or volcanic forcing, or from internal
    variability -- all of which are absent from the configuration entirely --
    and the result says so rather than leaving a reader to supply the caveat.
    """
    text = experiment.describe()
    assert "cannot distinguish" in text
    for alternative in ("greenhouse", "aerosol", "volcanic", "internal variability"):
        assert alternative in text
    assert any(
        "not on the physics of solar forcing itself" in limitation or "sufficiency" in limitation
        for limitation in experiment.claim.limitations
    )


# --- the paperwork ----------------------------------------------------------


@needs_data
@pytest.mark.req("INF-001", "INF-003", "INF-007")
def test_the_records_are_declared_in_their_roles_and_are_not_contaminated(experiment):
    """The irradiance was supplied to the run; the temperature was compared with it.

    A variable prescribed to a run can never afterwards be independent evidence
    about it, and declaring the roles is what lets the registry say so.
    """
    assert DatasetRole.BOUNDARY_CONDITION.shapes_the_model
    assert DatasetRole.HISTORICAL_RECONSTRUCTION_EVALUATION.is_evidence
    assert not DatasetRole.HISTORICAL_RECONSTRUCTION_EVALUATION.claims_independence

    # Every step between the model's output and the comparison is recorded.
    operations = {t.operation for t in experiment.transformations}
    assert any(op.startswith("anomaly:") for op in operations)
    assert any(op.startswith("resample:") for op in operations)
    assert any(op.startswith("align:") for op in operations)


@needs_data
@pytest.mark.req("VAL-003")
def test_the_model_is_differenced_against_the_records_own_baseline(experiment):
    """1961-1990, because that is HadCRUT5's, not because it is convenient."""
    from datetime import UTC, datetime

    from mfespf.data.loaders import seconds_since_epoch

    assert BASELINE_START == seconds_since_epoch(datetime(1961, 1, 1, tzinfo=UTC))
    assert BASELINE_END == seconds_since_epoch(datetime(1991, 1, 1, tzinfo=UTC))
    anomaly = next(t for t in experiment.transformations if t.operation.startswith("anomaly:"))
    assert f"baseline_window=[{BASELINE_START:g},{BASELINE_END:g}]" in anomaly.parameters
    assert "subtracted=" in anomaly.parameters
    # The model's own mean over the window is what was removed, so the modelled
    # anomaly averages to zero there by construction.
    assert abs(experiment.comparison.bias) < 0.05


@needs_data
@pytest.mark.req("CON-009", "ACC-004")
def test_the_run_behind_the_experiment_is_admissible(experiment):
    """A refutation from a run whose budget did not close would be worthless."""
    assert experiment.worst_relative_residual < 1e-12


@needs_data
@pytest.mark.req("INF-002")
def test_the_model_output_is_classified_as_model_output(experiment):
    """It is compared against a record; it is not one."""
    assert not ProvenanceClass.MODEL_OUTPUT.constrains_the_world
    assert ProvenanceClass.PROXY.constrains_the_world
