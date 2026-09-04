"""Stage 6 of the validation ladder: historical replay.

Charter section 19 stage 6 -- declared historical inputs, withheld output
records, frozen equations, energy-budget reconstruction -- and requirements
DAT-008, INF-001, INF-002, INF-007, KER-006, TIM-004, TIM-005, TIM-006, VAL-006,
VAL-007.

The physics here is Stage 1's. What is under test is everything around it: where
the boundary condition came from, what was done to it on the way in, whether the
equations were still the frozen ones when the evaluation ran, whether the budget
closes over the whole replay rather than interval by interval, and what the run
is allowed to be called afterwards.

The record shipped with the framework is synthetic. That is stated in the record,
enforced by the registry, and enforced again by the evaluator, which caps what
this stage may conclude at EMPIRICALLY_UNRESOLVED however well the comparison
goes. So Stage 6 demonstrates the pathway and makes no claim about the Earth --
and the last test in this file is the one that checks the framework refuses to
let it pretend otherwise.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import InvalidForScientificInferenceError, ProvenanceError
from mfespf.core.versioning import Version
from mfespf.data.records import ProvenanceClass, SamplingConvention
from mfespf.data.synthetic import SYNTHETIC_RECORD_ID, synthetic_irradiance_adapter
from mfespf.evaluation.reconstruction import reconstruct_budget
from mfespf.evaluation.withheld import FreezeRecord, WithheldDataset, WithheldEvaluator
from mfespf.physics.equations import build_equation_registry
from mfespf.registries.claims import ClaimRecord, ClaimStatus, ExpectedSignature
from mfespf.registries.datasets import DatasetRole, DatasetUse
from mfespf.run.historical import REPLAY_CLAIM_ID, historical_replay, replay_steps
from mfespf.run.modes import OperatingMode

pytestmark = pytest.mark.stage(6)

DAY = 86400.0
MONTH = 30.0 * DAY
VARIABLE = "total_solar_irradiance"


def _replay(master_dt: float = MONTH, **kwargs):
    kernel = historical_replay(master_dt=master_dt, **kwargs)
    steps = replay_steps(kernel.adapter, VARIABLE, master_dt)
    return kernel, kernel.run(steps)


def _claim(**overrides) -> ClaimRecord:
    fields = {
        "record_id": REPLAY_CLAIM_ID,
        "version": Version(1, 0, 0),
        "statement": (
            "Replaying a declared irradiance record through the airless body reproduces "
            "the surface temperature history that record implies."
        ),
        "mechanism": "Radiative balance of an isothermal layer under a prescribed input.",
        "physical_basis": "First law plus grey-body emission.",
        "governing_equation_ids": ("EQ-SLAB-ENERGY-BALANCE", "EQ-GREYBODY-EMISSION"),
        "failure_criterion": "root-mean-square surface temperature error above 0.05 K",
        "support_criterion": "root-mean-square surface temperature error at or below 0.05 K",
        "expected": ExpectedSignature(
            sign="positive",
            magnitude_bound="within 0.05 K",
            temporal="follows the record's eleven-year cycle with a lag of order a year",
        ),
        "withheld_data": (SYNTHETIC_RECORD_ID,),
        "uncertainty_method": "none propagated; the record's uncertainty is declared not used",
    }
    fields.update(overrides)
    return ClaimRecord(**fields)  # type: ignore[arg-type]


# --- declared historical inputs --------------------------------------------


@pytest.mark.req("INF-001", "INF-002", "DAT-009")
def test_the_record_is_registered_with_its_provenance_and_its_role():
    """The two facts a reader needs, and they are different facts.

    The role says the irradiance was supplied to the run. The provenance class
    says the numbers were generated. Neither implies the other, and a reader
    given only the first would have no way to tell that agreement with this
    record means nothing about the Sun.
    """
    kernel = historical_replay()
    record = kernel.datasets.get(SYNTHETIC_RECORD_ID)

    assert record.provenance_class is ProvenanceClass.SYNTHETIC
    assert not record.provenance_class.constrains_the_world
    assert "GENERATED" in record.source
    assert record.checksum

    uses = kernel.datasets.uses(claim_id=REPLAY_CLAIM_ID)
    assert len(uses) == 1
    assert uses[0].role is DatasetRole.BOUNDARY_CONDITION
    assert uses[0].role.shapes_the_model
    assert "solar" in uses[0].influenced


@pytest.mark.req("DAT-008", "INF-007", "TIM-004")
def test_the_resampling_is_recorded_on_the_dataset_it_transformed():
    """DAT-008 wants preprocessing reproducible from a versioned record.

    The transformation is written onto the registered dataset, so the registry
    alone says what was done: the operation, its parameters, the software, and
    the effect on uncertainty and on frequency content.
    """
    kernel = historical_replay()
    record = kernel.datasets.get(SYNTHETIC_RECORD_ID)

    assert len(record.transformations) == 1
    step = record.transformations[0]
    assert step.operation == f"resample:{VARIABLE}"
    assert "aggregation=" in step.parameters and "interpolation=" in step.parameters
    assert "source_coverage=" in step.parameters
    assert step.software
    assert "UNRESOLVED" in step.uncertainty_effect
    assert "frequency content" in step.uncertainty_effect


@pytest.mark.req("TIM-004")
def test_a_coarse_model_clock_reports_the_variability_it_removed():
    """A model asked to reproduce a smoothed record has been asked for less."""
    fine = historical_replay(master_dt=MONTH)
    coarse = historical_replay(master_dt=12.0 * MONTH, run_id="STAGE6-COARSE")

    assert not fine.adapter.reports[0].is_coarsening
    assert coarse.adapter.reports[0].is_coarsening
    assert coarse.adapter.reports[0].variance_retained < fine.adapter.reports[0].variance_retained
    assert "not represented after resampling" in coarse.adapter.reports[0].lost_frequency_note


@pytest.mark.req("KER-010", "MOD-001")
def test_the_recorded_source_reproduces_its_record_exactly():
    """A third implementation of the solar interface, and it plays back the record.

    Exactly, not closely: the run's boundary condition must be the resampled
    record itself, or the transformation report describes something other than
    what the model actually saw.
    """
    kernel, result = _replay()
    prepared = kernel.adapter.prepared(VARIABLE)
    played = np.array([s.outputs["solar_total_irradiance_1au"].scalar() for s in result.steps])
    assert np.array_equal(played, prepared.values[: played.size])

    solar = next(s.module for s in kernel.schedule.ordered() if s.module_id == "solar")
    assumptions = " ".join(solar.contract.assumptions)
    assert SYNTHETIC_RECORD_ID in solar.contract.module_name
    assert "synthetic" in assumptions
    assert "no claim about solar variability" in solar.contract.validity_domain


@pytest.mark.req("DAT-007", "SYS-006")
def test_running_past_the_record_coverage_fails():
    """A boundary condition invented past the end of a record is not a record.

    The failure is the feature. A framework that held the last value would
    produce a longer run that looks the same and is partly fabricated.
    """
    kernel = historical_replay()
    available = replay_steps(kernel.adapter, VARIABLE, MONTH)
    result = kernel.run(available + 5)

    # Recorded rather than merely raised: the run stops, the manifest says why,
    # and the partial output it produced is marked unusable as evidence.
    assert not result.may_support_a_claim
    assert len(result.steps) == available
    assert any("beyond its coverage" in w for w in result.manifest.warnings)
    assert result.manifest.status == "failed"


# --- energy-budget reconstruction ------------------------------------------


@pytest.mark.req("CON-009", "KER-013", "ACC-004")
def test_the_budget_closes_over_the_whole_replay_not_just_per_interval():
    """The question the per-interval checks cannot ask.

    Every interval closes to round-off, and has since Stage 0. Whether four
    hundred of them still close when summed depends on whether the errors are
    biased, and that is only visible after the sum.
    """
    _, result = _replay()
    assert result.may_support_a_claim
    for report in result.conservation_reports:
        assert report.passed, report.describe()

    reconstruction = reconstruct_budget(result.conservation_reports)
    assert reconstruction.n_intervals > 400
    assert reconstruction.passed(1e-12)
    assert reconstruction.relative_residual < 1e-14
    # Expressed as a flux, because that is the form a reader can weigh against
    # something physical: a drift of 0.01 W m^-2 would matter, 1e-13 does not.
    assert abs(reconstruction.mean_imbalance) < 1e-9
    assert "solar_input_across_system_boundary" in reconstruction.source_labels


@pytest.mark.req("CON-009")
def test_the_reconstruction_names_the_system_and_the_boundary():
    """Which stores are inside is a modelling decision, and it is recorded."""
    _, result = _replay()
    reconstruction = reconstruct_budget(result.conservation_reports)
    assert "surface" in reconstruction.system_modules
    assert reconstruction.boundary_modules == ("space",)
    assert "space" not in reconstruction.system_modules
    assert "boundary:   space" in reconstruction.describe()


# --- frozen equations ------------------------------------------------------


@pytest.mark.req("EQR-001", "INF-005", "INF-006")
def test_the_equation_set_is_pinned_by_one_hash_over_the_whole_registry():
    """Charter stage 6's *frozen equations*, and why naming three is not enough.

    A claim frozen against the equations it cites says nothing about a fourth
    that was added or edited before the evaluation ran. The hash covers the set,
    so any change to any registered equation invalidates the freeze.
    """
    equations = build_equation_registry()
    before = equations.content_hash()

    claim = _claim().freeze()
    evaluator = WithheldEvaluator(equations=equations)
    freeze = FreezeRecord(
        claim_id=claim.record_id,
        information_cutoff="the record's final sample",
        source_revision="test",
        configuration_hash="c",
        parameter_set_hash="p",
        criteria_hash=claim.criteria_hash(),
        equation_set_hash=before,
    )
    evaluator.freeze(claim, freeze)
    assert "equation set hash" in freeze.describe()

    # Supersede one equation with a different validity domain: nothing the claim
    # cites by name has changed, and the freeze must still be invalidated.
    from mfespf.registries.equations import ValidityDomain

    record = equations.get("EQ-BOUNDARY-ACCUMULATION")
    from dataclasses import replace

    equations.supersede(
        record.record_id,
        record.version,
        replace(
            record,
            version=Version(1, 1, 0),
            validity_domain=ValidityDomain(conditions=("edited after the freeze",)),
        ),
    )
    assert equations.content_hash() != before
    with pytest.raises(ProvenanceError, match="not the model that was frozen"):
        evaluator.evaluate(
            claim=claim,
            prediction=np.zeros(4),
            withheld=WithheldDataset(
                SYNTHETIC_RECORD_ID,
                "1.0.0",
                lambda: np.zeros(4),
                claim_id=claim.record_id,
                provenance_class=ProvenanceClass.SYNTHETIC,
            ),
        )


# --- withheld output records ------------------------------------------------


@pytest.mark.req("INF-008", "KER-008", "VAL-006")
def test_the_replay_evaluates_against_a_withheld_record_through_the_frozen_path():
    """The whole Stage 6 pathway, run end to end.

    A withheld temperature record, a frozen claim, a frozen equation set, and an
    evaluation that will not touch the data until every gate has passed.
    """
    kernel, result = _replay()
    predicted = np.array([s.outputs["surface_temperature"].scalar() for s in result.steps])

    # The "record" is the same trajectory with a small offset: this exercises the
    # pathway, and the next test is about what it is therefore allowed to mean.
    observations = predicted + 0.01
    withheld = WithheldDataset(
        SYNTHETIC_RECORD_ID,
        "1.0.0",
        lambda: observations,
        claim_id=REPLAY_CLAIM_ID,
        variables=("surface_temperature",),
        provenance_class=ProvenanceClass.SYNTHETIC,
    )
    assert not withheld.is_released
    with pytest.raises(ProvenanceError, match="before its claim was frozen"):
        withheld.values()

    claim = _claim().freeze()
    kernel.datasets.declare_use(
        DatasetUse(
            dataset_id=SYNTHETIC_RECORD_ID,
            dataset_version="1.0.0",
            role=DatasetRole.WITHHELD_HINDCAST_EVALUATION,
            claim_id=REPLAY_CLAIM_ID,
            variables=("surface_temperature",),
        )
    )
    equations = kernel.equations
    evaluator = WithheldEvaluator(datasets=kernel.datasets, equations=equations)
    evaluator.freeze(
        claim,
        FreezeRecord(
            claim_id=claim.record_id,
            information_cutoff="the record's final sample",
            source_revision="test",
            configuration_hash=result.manifest.configuration_hash(),
            parameter_set_hash=result.manifest.parameter_set_version,
            criteria_hash=claim.criteria_hash(),
            equation_set_hash=equations.content_hash(),
        ),
    )
    evaluation = evaluator.evaluate(
        claim=claim, prediction=predicted, withheld=withheld, threshold=0.05
    )

    assert withheld.is_released
    assert evaluation.metric_value == pytest.approx(0.01, rel=1e-6)
    assert evaluation.metric_value <= evaluation.threshold
    assert evaluation.n_observations == predicted.size


@pytest.mark.req("INF-002", "VAL-006", "CLM-002")
def test_a_generated_record_cannot_support_an_empirical_claim():
    """The refusal this stage exists to demonstrate.

    The comparison passed its threshold. The status is still
    EMPIRICALLY_UNRESOLVED, because agreement with numbers a computation
    produced is a test of the evaluation pathway and not of the physics against
    the world. A framework that returned SUPPORTED here would have recorded a
    self-consistency check as observational support, which is the specific
    failure the provenance machinery exists to prevent.
    """
    kernel, result = _replay()
    predicted = np.array([s.outputs["surface_temperature"].scalar() for s in result.steps])
    withheld = WithheldDataset(
        SYNTHETIC_RECORD_ID,
        "1.0.0",
        lambda: predicted + 0.001,
        claim_id=REPLAY_CLAIM_ID,
        provenance_class=ProvenanceClass.SYNTHETIC,
    )
    claim = _claim().freeze()
    kernel.datasets.declare_use(
        DatasetUse(
            dataset_id=SYNTHETIC_RECORD_ID,
            dataset_version="1.0.0",
            role=DatasetRole.WITHHELD_HINDCAST_EVALUATION,
            claim_id=REPLAY_CLAIM_ID,
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

    assert evaluation.metric_value < evaluation.threshold  # it passed
    assert evaluation.status is ClaimStatus.EMPIRICALLY_UNRESOLVED  # and still concludes nothing
    assert not evaluation.passed
    assert any("does not constrain the world" in note for note in evaluation.notes)
    assert "tests the evaluation pathway, not the physics" in evaluation.describe()


@pytest.mark.req("INF-002")
def test_an_undeclared_provenance_cannot_be_offered_as_evidence():
    """A record used as evidence has to say what kind of thing it is."""
    kernel = historical_replay()
    from dataclasses import replace

    record = kernel.datasets.get(SYNTHETIC_RECORD_ID)
    kernel.datasets.supersede(
        record.record_id,
        record.version,
        replace(
            record,
            version=Version(2, 0, 0),
            provenance_class=ProvenanceClass.UNDECLARED,
        ),
    )
    with pytest.raises(ProvenanceError, match="declares no provenance class"):
        kernel.datasets.declare_use(
            DatasetUse(
                dataset_id=SYNTHETIC_RECORD_ID,
                dataset_version="2.0.0",
                role=DatasetRole.PROCESS_VALIDATION,
                claim_id="CLM-ANY",
            )
        )


# --- what the run may be called ---------------------------------------------


@pytest.mark.req("TIM-005", "TIM-006", "VAL-007")
def test_a_prescribed_history_run_may_not_be_reported_as_a_prediction():
    """The run was handed the boundary conditions of the period it replays.

    Its agreement is not evidence of prospective skill, and the framework
    refuses the description rather than relying on the report's author to
    remember. This is the same arithmetic a genuine hindcast would produce,
    which is exactly why it needs a check rather than a convention.
    """
    _, result = _replay()
    assert result.operating_mode is OperatingMode.PRESCRIBED_HISTORY
    assert result.manifest.to_dict()["operating_mode"] == "prescribed_history"

    for overstated in (
        OperatingMode.HISTORICAL_RECONSTRUCTION,
        OperatingMode.HINDCAST,
        OperatingMode.OUT_OF_SAMPLE_PREDICTION,
    ):
        with pytest.raises(InvalidForScientificInferenceError, match="evidentiary language"):
            result.assert_may_be_reported_as(overstated)

    # Describing it as what it is, or as something weaker, is fine.
    result.assert_may_be_reported_as(OperatingMode.PRESCRIBED_HISTORY)
    result.assert_may_be_reported_as(OperatingMode.CONTROLLED_EXPERIMENT)


@pytest.mark.req("AUD-002", "CMP-009")
def test_the_record_reaches_the_configuration_hash_and_the_run_reproduces():
    """Two replays of different records must not hash alike, and one must repeat.

    The dataset's checksum is in the solar module's parameters, so a run driven
    by a different record is a different configuration even though every line of
    code is the same.
    """
    baseline = historical_replay(run_id="A").build_manifest()
    repeat = historical_replay(run_id="A").build_manifest()
    assert baseline.configuration_hash() == repeat.configuration_hash()

    other = historical_replay(
        adapter=synthetic_irradiance_adapter(cycle_amplitude=2.0), run_id="B"
    ).build_manifest()
    assert other.configuration_hash() != baseline.configuration_hash()
    assert (
        baseline.parameters["solar"]["dataset_checksum"]
        != other.parameters["solar"]["dataset_checksum"]
    )
    assert baseline.parameters["solar"]["dataset_provenance_class"] == "synthetic"
    assert baseline.parameters["solar"]["record_sampling"] == SamplingConvention.INTERVAL_MEAN
