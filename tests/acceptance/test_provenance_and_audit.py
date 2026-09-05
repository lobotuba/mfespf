"""ACC-008, ACC-009 and the audit trail.

The two acceptance criteria here are about evidence rather than physics: a
calibration/validation contamination must be rejected by the provenance
controls, and a parent/child run pair must produce an accurate difference
record.  Both guard against the same class of error -- a conclusion that looks
supported because the reader cannot see what the model was told.
"""

from __future__ import annotations

import numpy as np
import pytest

import mfespf
from mfespf.core.errors import ProvenanceError
from mfespf.core.versioning import Version
from mfespf.evaluation.withheld import FreezeRecord, WithheldDataset, WithheldEvaluator
from mfespf.registries.claims import ClaimRecord, ClaimRegistry, ClaimStatus, ExpectedSignature
from mfespf.registries.datasets import (
    DatasetRecord,
    DatasetRegistry,
    DatasetRole,
    DatasetUse,
    ProvenanceClass,
)
from mfespf.run.configurations import airless_radiating_body
from mfespf.run.manifest import FRAMEWORK_VERSION, RunManifest, diff_manifests

pytestmark = pytest.mark.stage(0)

V1 = Version(1, 0, 0)


def _dataset(
    dataset_id: str,
    variables: tuple[str, ...] = ("temperature",),
    provenance_class: ProvenanceClass = ProvenanceClass.OBSERVATIONAL,
) -> DatasetRecord:
    # These records stand in for measured ones, so they declare themselves
    # observational: a record offered as evidence must say what kind of thing it
    # is, and the registry refuses an undeclared one in an evidence role.
    return DatasetRecord(
        record_id=dataset_id,
        version=V1,
        title=f"observations {dataset_id}",
        source="stand-in for a measured record",
        checksum="0" * 16,
        variables=variables,
        uncertainty_record="declared in the source record",
        provenance_class=provenance_class,
    )


def _claim(claim_id: str = "CLM-TEST-001", **kwargs) -> ClaimRecord:
    defaults = {
        "statement": "The slab reaches the analytic radiative equilibrium temperature.",
        "mechanism": "Radiative balance of an isothermal layer.",
        "physical_basis": "First law plus grey-body emission.",
        "governing_equation_ids": ("EQ-SLAB-ENERGY-BALANCE",),
        "failure_criterion": "RMSE above 0.5 K against withheld observations.",
        "support_criterion": "RMSE at or below 0.5 K.",
        "expected": ExpectedSignature(sign="positive", magnitude_bound="within 0.5 K"),
        "uncertainty_method": "none propagated; discretization error reported separately",
    }
    defaults.update(kwargs)
    return ClaimRecord(record_id=claim_id, version=V1, **defaults)


# --- ACC-008: contamination ------------------------------------------------


@pytest.mark.req("ACC-008", "INF-003", "VAL-006")
def test_calibration_data_cannot_be_reported_as_independent_validation():
    registry = DatasetRegistry()
    registry.register(_dataset("OBS-A"))
    registry.declare_use(
        DatasetUse("OBS-A", "1.0.0", DatasetRole.PARAMETER_CALIBRATION, "CLM-1", ("temperature",))
    )
    registry.declare_use(
        DatasetUse(
            "OBS-A",
            "1.0.0",
            DatasetRole.OUT_OF_SAMPLE_PREDICTION_EVALUATION,
            "CLM-1",
            ("temperature",),
        )
    )

    findings = registry.check_contamination("CLM-1")
    assert len(findings) == 1
    assert "calibration data cannot be reported as independent validation" in findings[0].describe()
    with pytest.raises(ProvenanceError, match="information-provenance validation failed"):
        registry.assert_clean("CLM-1")


@pytest.mark.req("INF-003")
def test_prescribing_a_variable_and_then_evaluating_against_it_is_contamination():
    registry = DatasetRegistry()
    registry.register(_dataset("REANALYSIS", ("wind", "temperature")))
    registry.declare_use(
        DatasetUse("REANALYSIS", "1.0.0", DatasetRole.BOUNDARY_CONDITION, "CLM-2", ("temperature",))
    )
    registry.declare_use(
        DatasetUse(
            "REANALYSIS",
            "1.0.0",
            DatasetRole.HISTORICAL_RECONSTRUCTION_EVALUATION,
            "CLM-2",
            ("temperature",),
        )
    )
    findings = registry.check_contamination("CLM-2")
    assert findings and "partly definitional" in findings[0].describe()


@pytest.mark.req("INF-001", "INF-003")
def test_disjoint_variables_from_one_dataset_are_not_contamination():
    """Granularity matters: using the winds does not spend the radiation record."""
    registry = DatasetRegistry()
    registry.register(_dataset("REANALYSIS", ("wind", "olr")))
    registry.declare_use(
        DatasetUse("REANALYSIS", "1.0.0", DatasetRole.BOUNDARY_CONDITION, "CLM-3", ("wind",))
    )
    registry.declare_use(
        DatasetUse("REANALYSIS", "1.0.0", DatasetRole.PROCESS_VALIDATION, "CLM-3", ("olr",))
    )
    assert registry.check_contamination("CLM-3") == []
    registry.assert_clean("CLM-3")


@pytest.mark.req("INF-004", "A4.6" if False else "INF-004")
def test_independence_statement_reports_what_was_available_at_each_stage():
    registry = DatasetRegistry()
    registry.register(_dataset("OBS-CAL"))
    registry.register(_dataset("OBS-HELD"))
    registry.declare_use(DatasetUse("OBS-CAL", "1.0.0", DatasetRole.PARAMETER_CALIBRATION, "CLM-4"))
    registry.declare_use(
        DatasetUse("OBS-HELD", "1.0.0", DatasetRole.WITHHELD_HINDCAST_EVALUATION, "CLM-4")
    )

    statement = registry.independence_statement("CLM-4")
    assert any("OBS-CAL" in e for e in statement["parameter_determination"])
    assert any("OBS-HELD" in e for e in statement["explicitly_withheld"])
    assert not any("OBS-HELD" in e for e in statement["parameter_determination"])


@pytest.mark.req("DAT-009")
def test_a_dataset_without_a_checksum_cannot_be_registered():
    with pytest.raises(Exception, match="checksum"):
        DatasetRecord(record_id="NO-SUM", version=V1, title="x", source="y", checksum="")


# --- Withheld data (KER-008, INF-008) --------------------------------------


@pytest.mark.req("INF-008", "KER-008")
def test_withheld_data_cannot_be_read_before_the_claim_is_frozen():
    dataset = WithheldDataset("OBS-HELD", "1.0.0", lambda: np.array([1.0, 2.0]), claim_id="CLM-5")
    with pytest.raises(ProvenanceError, match="accessed before its claim was frozen"):
        dataset.values()


@pytest.mark.req("HYP-009", "INF-005", "INF-006")
def test_thresholds_cannot_move_after_the_freeze():
    claim = _claim("CLM-6").freeze()
    evaluator = WithheldEvaluator()
    freeze = FreezeRecord(
        claim_id="CLM-6",
        information_cutoff="2026-01-01",
        source_revision="abc123",
        configuration_hash="cfg",
        parameter_set_hash="par",
        criteria_hash=claim.criteria_hash(),
    )
    evaluator.freeze(claim, freeze)

    relaxed = claim.__class__(**{**claim.__dict__, "failure_criterion": "RMSE above 50 K."})
    with pytest.raises(ProvenanceError, match="criteria changed after freeze"):
        relaxed.assert_criteria_intact()


@pytest.mark.req("KER-008", "CLM-003")
def test_a_frozen_claim_evaluates_against_withheld_data_exactly_once():
    registry = DatasetRegistry()
    registry.register(_dataset("OBS-HELD"))
    registry.declare_use(
        DatasetUse("OBS-HELD", "1.0.0", DatasetRole.OUT_OF_SAMPLE_PREDICTION_EVALUATION, "CLM-7")
    )
    claim = _claim("CLM-7", withheld_data=("OBS-HELD",)).freeze()
    evaluator = WithheldEvaluator(datasets=registry)
    evaluator.freeze(
        claim,
        FreezeRecord(
            claim_id="CLM-7",
            information_cutoff="2026-01-01",
            source_revision="abc123",
            configuration_hash="cfg",
            parameter_set_hash="par",
            criteria_hash=claim.criteria_hash(),
        ),
    )

    observations = np.array([254.5, 254.6, 254.7])
    withheld = WithheldDataset("OBS-HELD", "1.0.0", lambda: observations, claim_id="CLM-7")
    result = evaluator.evaluate(
        claim=claim,
        prediction=observations + 0.1,
        withheld=withheld,
        threshold=0.5,
    )

    assert result.status is ClaimStatus.SUPPORTED_WITHIN_SCOPE
    assert result.metric_value == pytest.approx(0.1, abs=1e-9)
    # CLM-003: never an unqualified confirmation.
    assert any("within the tested configuration" in n for n in result.notes)


@pytest.mark.req("KER-008", "CLM-002")
def test_a_prediction_outside_the_declared_threshold_is_refuted():
    registry = DatasetRegistry()
    registry.register(_dataset("OBS-HELD"))
    registry.declare_use(
        DatasetUse("OBS-HELD", "1.0.0", DatasetRole.WITHHELD_HINDCAST_EVALUATION, "CLM-8")
    )
    claim = _claim("CLM-8").freeze()
    evaluator = WithheldEvaluator(datasets=registry)
    evaluator.freeze(
        claim,
        FreezeRecord(
            claim_id="CLM-8",
            information_cutoff="2026-01-01",
            source_revision="abc123",
            configuration_hash="cfg",
            parameter_set_hash="par",
            criteria_hash=claim.criteria_hash(),
        ),
    )
    observations = np.array([254.5, 254.6])
    result = evaluator.evaluate(
        claim=claim,
        prediction=observations + 5.0,
        withheld=WithheldDataset("OBS-HELD", "1.0.0", lambda: observations, claim_id="CLM-8"),
        threshold=0.5,
    )
    assert result.status is ClaimStatus.REFUTED


@pytest.mark.req("INF-008")
def test_evaluation_is_refused_when_the_dataset_is_contaminated():
    registry = DatasetRegistry()
    registry.register(_dataset("OBS-BOTH"))
    for role in (
        DatasetRole.PARAMETER_CALIBRATION,
        DatasetRole.OUT_OF_SAMPLE_PREDICTION_EVALUATION,
    ):
        registry.declare_use(DatasetUse("OBS-BOTH", "1.0.0", role, "CLM-9"))

    claim = _claim("CLM-9").freeze()
    evaluator = WithheldEvaluator(datasets=registry)
    evaluator.freeze(
        claim,
        FreezeRecord(
            claim_id="CLM-9",
            information_cutoff="2026-01-01",
            source_revision="abc",
            configuration_hash="cfg",
            parameter_set_hash="par",
            criteria_hash=claim.criteria_hash(),
        ),
    )
    observations = np.array([1.0])
    with pytest.raises(ProvenanceError, match="information-provenance validation failed"):
        evaluator.evaluate(
            claim=claim,
            prediction=observations,
            withheld=WithheldDataset("OBS-BOTH", "1.0.0", lambda: observations, claim_id="CLM-9"),
            threshold=1.0,
        )


# --- ACC-009: lineage ------------------------------------------------------


@pytest.mark.req("ACC-009", "AUD-003", "AUD-004")
def test_a_parent_child_run_pair_produces_an_accurate_difference_record():
    parent = airless_radiating_body(run_id="PARENT-001", albedo=0.3).build_manifest()
    child = airless_radiating_body(
        run_id="CHILD-001", albedo=0.5, master_dt=1800.0
    ).build_manifest()
    child.parent_run_ids = ("PARENT-001",)
    child.reason_for_change = "raise the albedo and halve the master step"

    report = diff_manifests(parent, child)

    assert not report.is_empty
    assert any("master_dt" in e for e in report.changes["coupling"])
    assert any("temporal_resolution" in e for e in report.changes["resolution"])
    described = report.describe()
    assert "PARENT-001" in described and "CHILD-001" in described


@pytest.mark.req("AUD-003")
def test_an_identical_configuration_produces_an_empty_difference_record():
    a = airless_radiating_body(run_id="A").build_manifest()
    b = airless_radiating_body(run_id="B").build_manifest()
    report = diff_manifests(a, b)
    assert report.is_empty
    assert "configurationally identical" in report.describe()
    assert a.configuration_hash() == b.configuration_hash()


@pytest.mark.req("AUD-002", "AUD-001")
def test_the_configuration_hash_covers_physics_and_not_wall_clock():
    a = airless_radiating_body(run_id="A").build_manifest()
    b = airless_radiating_body(run_id="B").build_manifest()
    b.execution_time = "1999-01-01T00:00:00+00:00"
    assert a.configuration_hash() == b.configuration_hash()

    c = airless_radiating_body(run_id="C", albedo=0.31).build_manifest()
    assert c.configuration_hash() != a.configuration_hash()


@pytest.mark.req("AUD-005", "AUD-006")
def test_retained_evidence_cannot_be_overwritten(tmp_path):
    first = RunManifest(run_id="RUN-1", purpose="original")
    path = tmp_path / "manifest.json"
    first.write(path)

    # Rewriting the identical manifest is a no-op and permitted.
    first.write(path)

    second = RunManifest(run_id="RUN-1", purpose="quietly corrected")
    with pytest.raises(ProvenanceError, match="Retained evidence is immutable"):
        second.write(path)


@pytest.mark.req("AUD-008", "SYS-008")
def test_a_failed_gate_invalidates_the_run_for_inference():
    manifest = RunManifest(run_id="RUN-2", purpose="test")
    manifest.record_gate("conservation_reconciliation", False)
    assert not manifest.may_support_a_claim
    assert manifest.invalidation_reasons


@pytest.mark.req("CLM-001")
def test_a_claim_without_a_failure_criterion_cannot_be_registered():
    registry = ClaimRegistry()
    with pytest.raises(Exception, match="failure criterion"):
        registry.register(_claim("CLM-BAD", failure_criterion="   "))


@pytest.mark.req("INF-003", "CLM-001")
def test_a_claim_cannot_list_one_dataset_as_both_calibration_and_withheld():
    with pytest.raises(Exception, match="INF-003"):
        _claim("CLM-BAD2", calibration_data=("OBS-A",), withheld_data=("OBS-A",))


@pytest.mark.req("AUD-008", "SYS-008")
def test_the_manifest_records_the_version_of_the_code_that_produced_the_run():
    """A manifest whose framework version is a separate literal will go stale.

    This one did. It was a third copy of the version string beside the one in
    ``pyproject.toml`` and the one in ``mfespf.__version__``, and it stopped
    being updated an entire release ago, so every run written after 0.11.0
    recorded 0.10.0. Nothing failed, which is the problem: a manifest exists to
    say which code produced a number, and one that says the wrong thing is worse
    than one that says nothing.
    """
    assert FRAMEWORK_VERSION == mfespf.__version__
    assert RunManifest(run_id="R", purpose="p").framework_version == mfespf.__version__
