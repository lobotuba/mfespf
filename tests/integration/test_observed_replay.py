"""Replaying a real record, and the cap lifting because it is one.

Requirements DAT-009, INF-002, KER-010, TIM-001, TIM-005, TIM-006, VAL-006,
VAL-007, CON-009.

Stage 6 made a design claim with a synthetic record: that nothing in the replay
pathway is specific to the record being generated, and that swapping in a real
one would change the provenance class and nothing else. This is the test of that
claim -- the same builder, the same adapter, the same resampling, the same
evaluator, and a status that is allowed to reach SUPPORTED_WITHIN_SCOPE because
the record's lineage is a measurement rather than a computation.

Every test here skips when the files are absent. The bytes are not in the
repository; `python tools/fetch_data.py` puts them in `data/raw/`.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import InvalidForScientificInferenceError, ProvenanceError
from mfespf.core.versioning import Version
from mfespf.data.loaders import EPOCH_LABEL, available, load_declared
from mfespf.data.records import ProvenanceClass, UncertaintyKind
from mfespf.evaluation.reconstruction import reconstruct_budget
from mfespf.evaluation.withheld import FreezeRecord, WithheldDataset, WithheldEvaluator
from mfespf.registries.claims import ClaimRecord, ClaimStatus
from mfespf.registries.datasets import DatasetRole, DatasetUse
from mfespf.run.historical import observed_replay, replay_steps
from mfespf.run.modes import OperatingMode

pytestmark = pytest.mark.stage(6)

MONTH = 30.0 * 86400.0
VARIABLE = "total_solar_irradiance"
IRRADIANCE = "NRLTSI2-MONTHLY"

present = available()
needs_data = pytest.mark.skipif(
    IRRADIANCE not in present,
    reason=f"{IRRADIANCE} is not in data/raw; run `python tools/fetch_data.py`",
)


def _replay(master_dt: float = MONTH, **kwargs):
    kernel = observed_replay(master_dt=master_dt, **kwargs)
    return kernel, kernel.run(replay_steps(kernel.adapter, VARIABLE, master_dt))


def _claim() -> ClaimRecord:
    return ClaimRecord(
        record_id="CLM-STAGE6-REPLAY-001",
        version=Version(1, 0, 0),
        statement=(
            "The airless body driven by the NRLTSI2 irradiance record reproduces the "
            "equilibrium temperature that record implies, to within 0.05 K."
        ),
        mechanism="Radiative balance of an isothermal layer under a prescribed input.",
        physical_basis="First law plus grey-body emission.",
        governing_equation_ids=("EQ-SLAB-ENERGY-BALANCE", "EQ-GREYBODY-EMISSION"),
        failure_criterion="root-mean-square surface temperature error above 0.05 K",
        support_criterion="root-mean-square surface temperature error at or below 0.05 K",
        withheld_data=(IRRADIANCE,),
    ).freeze()


# --- the record ------------------------------------------------------------


@needs_data
@pytest.mark.req("DAT-009", "INF-002", "INF-007")
def test_the_loaded_record_carries_its_checksum_provenance_and_transformation():
    """The registry entry is built from the file, not typed alongside it."""
    adapter = load_declared(IRRADIANCE)
    record = adapter.record

    assert record.provenance_class is ProvenanceClass.PROXY
    assert record.checksum  # from the bytes on disk, via DatasetRecord.from_file
    assert "lasp.colorado.edu" in record.source
    assert "licence" in record.acquisition_record
    assert "citation" in record.acquisition_record

    step = record.transformations[0]
    assert step.operation == "parse:nrltsi2"
    assert EPOCH_LABEL in step.parameters
    assert step.software.endswith("load_nrltsi2")
    assert "standard_uncertainty" in record.uncertainty_record


@needs_data
@pytest.mark.req("DAT-009", "SYS-006")
def test_a_file_that_does_not_match_its_recorded_checksum_is_refused(tmp_path, monkeypatch):
    """Checked before the file is parsed, not after.

    Reading first and noticing second would mean the values had already reached
    the caller by the time anyone knew they were the wrong ones.
    """
    from mfespf.data import loaders

    corrupted = tmp_path / "raw"
    corrupted.mkdir()
    original = loaders.RAW / "nrl2_tsi_monthly.csv"
    (corrupted / "nrl2_tsi_monthly.csv").write_text(
        original.read_text(encoding="utf-8").replace("1360.5418701171875", "9999.0"),
        encoding="utf-8",
    )
    monkeypatch.setattr(loaders, "RAW", corrupted)

    assert IRRADIANCE not in loaders.available()
    with pytest.raises(ProvenanceError, match="does not match the checksum"):
        loaders.load_declared(IRRADIANCE)


@needs_data
@pytest.mark.req("UNC-002")
def test_the_published_uncertainty_is_carried_with_its_meaning():
    """NRLTSI2 publishes a standard uncertainty; HadCRUT5 publishes an interval.

    Both are loaded, and each says which it is. Combining them as though they
    were the same quantity would understate one by about a factor of two.
    """
    tsi = next(iter(load_declared(IRRADIANCE)._series.values()))
    assert tsi.uncertainty_kind is UncertaintyKind.STANDARD_UNCERTAINTY
    assert tsi.uncertainty_kind.is_one_sigma
    assert np.all(tsi.uncertainty > 0.0)

    if "HADCRUT5-GLOBAL-MONTHLY" in present:
        temperature = next(iter(load_declared("HADCRUT5-GLOBAL-MONTHLY")._series.values()))
        assert temperature.uncertainty_kind is UncertaintyKind.HALF_WIDTH_95
        assert not temperature.uncertainty_kind.is_one_sigma
        assert temperature.uncertainty_kind is not tsi.uncertainty_kind


# --- the run ---------------------------------------------------------------


@needs_data
@pytest.mark.req("TIM-001", "KER-010")
def test_the_run_starts_where_the_record_does():
    """A dated replay whose clock began at zero would fail on its first step.

    The record covers 1882 onward, which is not zero on a coordinate whose epoch
    is 1850, so the clock origin is a declared boundary condition of the run.
    """
    kernel, result = _replay()
    coverage = kernel.adapter.prepared(VARIABLE).coverage

    assert kernel.start_time == pytest.approx(coverage[0])
    assert result.steps[0].time == pytest.approx(coverage[0])
    assert result.manifest.boundary_conditions["master_start_time"] == pytest.approx(coverage[0])
    # Two replays over different periods of one record are different
    # configurations, so the clock origin has to reach the hash.
    later = observed_replay(master_dt=MONTH, run_id="B").build_manifest()
    assert "master_start_time" in later.configuration_payload()["boundary_conditions"]


@needs_data
@pytest.mark.req("KER-010", "MOD-001")
def test_the_run_plays_back_the_real_record_exactly():
    kernel, result = _replay()
    prepared = kernel.adapter.prepared(VARIABLE)
    played = np.array([s.outputs["solar_total_irradiance_1au"].scalar() for s in result.steps])

    assert played.size > 1500  # monthly from 1882
    assert np.array_equal(played, prepared.values[: played.size])
    assert 1358.0 < played.min() < played.max() < 1365.0  # a plausible solar constant


@needs_data
@pytest.mark.req("CON-009", "KER-013", "ACC-004")
def test_the_budget_closes_over_a_century_of_replay():
    _, result = _replay()
    assert result.may_support_a_claim
    for report in result.conservation_reports:
        assert report.passed, report.describe()

    reconstruction = reconstruct_budget(result.conservation_reports)
    assert reconstruction.duration > 100.0 * 365.25 * 86400.0
    assert reconstruction.passed(1e-12)
    assert abs(reconstruction.mean_imbalance) < 1e-9


# --- the claim the synthetic record could not support ------------------------


@needs_data
@pytest.mark.req("INF-002", "VAL-006", "CLM-002")
def test_a_measured_lineage_record_lifts_the_cap_a_generated_one_imposes():
    """The design claim Stage 6 made, tested by changing only the record.

    With the synthetic record the evaluator caps the status at
    EMPIRICALLY_UNRESOLVED however well the metric passes. Nothing else about
    the pathway differs here -- same builder, same adapter, same resampling,
    same evaluator -- and the status reaches SUPPORTED_WITHIN_SCOPE because the
    record's provenance class says it constrains the world.
    """
    kernel, result = _replay()
    predicted = np.array([s.outputs["surface_temperature"].scalar() for s in result.steps])
    withheld = WithheldDataset(
        IRRADIANCE,
        "1.0.0",
        lambda: predicted + 0.01,
        claim_id="CLM-STAGE6-REPLAY-001",
        provenance_class=ProvenanceClass.PROXY,
    )

    claim = _claim()
    kernel.datasets.declare_use(
        DatasetUse(
            dataset_id=IRRADIANCE,
            dataset_version="1.0.0",
            role=DatasetRole.WITHHELD_HINDCAST_EVALUATION,
            claim_id=claim.record_id,
            variables=("surface_temperature",),
        )
    )
    evaluator = WithheldEvaluator(datasets=kernel.datasets, equations=kernel.equations)
    evaluator.freeze(
        claim,
        FreezeRecord(
            claim_id=claim.record_id,
            information_cutoff="the record's final sample",
            source_revision="test",
            configuration_hash=result.manifest.configuration_hash(),
            parameter_set_hash="unversioned",
            criteria_hash=claim.criteria_hash(),
            equation_set_hash=kernel.equations.content_hash(),
        ),
    )
    evaluation = evaluator.evaluate(
        claim=claim, prediction=predicted, withheld=withheld, threshold=0.05
    )

    assert evaluation.provenance_class.constrains_the_world
    assert evaluation.status is ClaimStatus.SUPPORTED_WITHIN_SCOPE
    assert evaluation.passed
    assert not any("does not constrain the world" in note for note in evaluation.notes)
    # Still scoped: supported within the tested configuration, not in general.
    assert any("within the tested configuration" in note for note in evaluation.notes)


@needs_data
@pytest.mark.req("TIM-005", "TIM-006", "VAL-007")
def test_a_real_record_does_not_make_a_replay_into_a_prediction():
    """The provenance cap lifted; the operating-mode rule did not.

    These are two separate guards against two separate mistakes, and using a
    measured record relaxes only the first. The run was still handed the
    irradiance of the period it covers.
    """
    _, result = _replay()
    assert result.operating_mode is OperatingMode.PRESCRIBED_HISTORY
    for overstated in (OperatingMode.HINDCAST, OperatingMode.OUT_OF_SAMPLE_PREDICTION):
        with pytest.raises(InvalidForScientificInferenceError, match="evidentiary language"):
            result.assert_may_be_reported_as(overstated)
