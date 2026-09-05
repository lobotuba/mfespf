"""The hypothesis registry.

Requirements HYP-001 through HYP-010.

A claim says what would refute it. A hypothesis has to say more: where the
mechanism enters the model, what it does to the conserved quantities, what the
control is, and which competing explanations the proposed experiment can and
cannot separate it from. Every one of those is required rather than encouraged,
because each is a field that reads fine when left vague and destroys the
inference when it is.
"""

from __future__ import annotations

from dataclasses import replace

import pytest

from mfespf.core.errors import ConfigurationError, OwnershipError
from mfespf.core.versioning import Version
from mfespf.registries.claims import ExpectedSignature
from mfespf.registries.hypotheses import (
    EntryPoint,
    HypothesisRecord,
    HypothesisRegistry,
    MechanismClass,
)

pytestmark = pytest.mark.stage(6)

SIGNATURE = ExpectedSignature(
    sign="positive",
    magnitude_bound="within 0.1 K",
    temporal="follows the forcing with a one-year lag",
    spatial="global mean only",
    spectral="not applicable",
)


def _hypothesis(**overrides) -> HypothesisRecord:
    fields = {
        "record_id": "HYP-TEST-001",
        "version": Version(1, 0, 0),
        "statement": "A mechanism that does something measurable.",
        "claim_ids": ("CLM-TEST-001",),
        "physical_basis": "First law applied to a slab.",
        "governing_equation_ids": ("EQ-SLAB-ENERGY-BALANCE",),
        "required_inputs": ("total_solar_irradiance",),
        "produced_outputs": ("surface_temperature",),
        "expected": SIGNATURE,
        "validity_domain": "Global mean, annual to centennial.",
        "failure_criterion": "rms error above 0.1 K",
        "conserved_quantities": ("energy",),
        "transfer_accounting": "ordinary transfers through the TOA boundary",
        "entry_point": EntryPoint.BOUNDARY,
        "entry_point_detail": "the solar-source interface",
        "implementing_modules": ("solar", "surface"),
        "null_configuration": "the same configuration with constant irradiance",
        "distinguishes_from": ("no forcing at all",),
        "cannot_distinguish_from": ("any other mechanism with the same global-mean signature",),
    }
    fields.update(overrides)
    return HypothesisRecord(**fields)  # type: ignore[arg-type]


# --- the required declarations ---------------------------------------------


@pytest.mark.req("HYP-001", "HYP-002")
def test_a_hypothesis_names_a_claim_and_declares_its_basis():
    record = _hypothesis()
    assert record.claim_ids == ("CLM-TEST-001",)
    for missing, requirement in (
        ("claim_ids", "HYP-001"),
        ("physical_basis", "HYP-002"),
        ("validity_domain", "HYP-002"),
        ("failure_criterion", "HYP-002"),
    ):
        empty = () if missing == "claim_ids" else "   "
        with pytest.raises(ConfigurationError, match=requirement):
            _hypothesis(**{missing: empty})
    with pytest.raises(ConfigurationError, match="EQR-005"):
        _hypothesis(governing_equation_ids=())
    with pytest.raises(ConfigurationError, match="no inputs or no outputs"):
        _hypothesis(produced_outputs=())


@pytest.mark.req("HYP-003")
def test_an_unspecified_expected_signature_is_refused():
    """An unspecified sign is the commonest way a prediction becomes unfalsifiable.

    Spectral behaviour is the only one the charter qualifies with "where
    applicable", so it alone may say "not applicable"; the other four have to be
    stated even when the answer is that the configuration cannot resolve them.
    """
    for aspect in ("sign", "magnitude_bound", "temporal", "spatial"):
        with pytest.raises(ConfigurationError, match=f"expected {aspect} unspecified"):
            _hypothesis(expected=replace(SIGNATURE, **{aspect: "unspecified"}))
    # An explicit "not applicable" is a statement and is accepted.
    _hypothesis(
        expected=ExpectedSignature(
            sign="positive",
            magnitude_bound="within 0.1 K",
            temporal="follows the forcing",
            spatial="not resolvable in a zero-dimensional configuration",
            spectral="not applicable",
        )
    )


@pytest.mark.req("HYP-004")
def test_a_hypothesis_declares_what_it_does_to_the_conserved_quantities():
    """Silence is not the same as 'nothing', so silence is refused."""
    with pytest.raises(ConfigurationError, match="state 'none' explicitly"):
        _hypothesis(conserved_quantities=())
    with pytest.raises(ConfigurationError, match="how its transfers are accounted"):
        _hypothesis(transfer_accounting="  ")
    assert _hypothesis(conserved_quantities=("none",)).conserved_quantities == ("none",)


@pytest.mark.req("HYP-005")
def test_a_mechanism_enters_through_one_of_the_declared_interfaces():
    """The charter's list is closed, and that is the point.

    A mechanism that fits none of these has nowhere to act that the accounting
    can see, and is therefore not a mechanism this framework can test.
    """
    assert {e.value for e in EntryPoint} == {
        "module",
        "tendency",
        "boundary",
        "source",
        "sink",
        "constitutive_relation",
        "diagnostic",
    }
    with pytest.raises(ConfigurationError, match="where it enters the model"):
        _hypothesis(entry_point_detail="")
    with pytest.raises(ConfigurationError, match="no implementing module"):
        _hypothesis(implementing_modules=())


@pytest.mark.req("HYP-005")
def test_the_implementing_modules_must_be_in_the_configuration():
    registry = HypothesisRegistry()
    registry.register(_hypothesis())
    registry.assert_modules_present("HYP-TEST-001", {"solar", "orbit", "surface", "space"})
    with pytest.raises(OwnershipError, match="cannot be tested by it"):
        registry.assert_modules_present("HYP-TEST-001", {"solar", "space"})


@pytest.mark.req("HYP-007")
def test_a_hypothesis_declares_its_null_control():
    with pytest.raises(ConfigurationError, match="null or control configuration"):
        _hypothesis(null_configuration="")


@pytest.mark.req("HYP-008")
def test_a_hypothesis_says_what_its_experiment_cannot_separate_it_from():
    """The field most likely to be left vague, and the one that decides the
    inference.

    A mechanism supported by an experiment that could not distinguish it from
    four alternatives is supported jointly with those four, and a record that
    never said so lets a reader mistake "consistent with" for "evidence for".
    """
    with pytest.raises(ConfigurationError, match="HYP-008"):
        _hypothesis(cannot_distinguish_from=())
    with pytest.raises(ConfigurationError, match="HYP-008"):
        _hypothesis(distinguishes_from=())

    registry = HypothesisRegistry()
    registry.register(_hypothesis())
    assert registry.undistinguished("HYP-TEST-001") == (
        "any other mechanism with the same global-mean signature",
    )


@pytest.mark.req("HYP-009")
def test_the_failure_threshold_freezes_and_an_edit_is_caught():
    record = _hypothesis().freeze()
    record.assert_criteria_intact()

    loosened = replace(record, failure_criterion="rms error above 10 K")
    with pytest.raises(ConfigurationError, match="changed after freezing"):
        loosened.assert_criteria_intact()
    with pytest.raises(ConfigurationError, match="has not been frozen"):
        _hypothesis().assert_criteria_intact()


@pytest.mark.req("HYP-009")
def test_the_frozen_hash_covers_what_the_experiment_could_not_separate():
    """Quietly shortening that list after the fact would strengthen the claim.

    It is therefore inside the frozen hash, alongside the threshold: the two are
    equally capable of turning a weak result into a strong-looking one.
    """
    record = _hypothesis().freeze()
    shortened = replace(record, cannot_distinguish_from=("almost nothing",))
    assert shortened.criteria_hash() != record.frozen_criteria_hash
    with pytest.raises(ConfigurationError, match="changed after freezing"):
        shortened.assert_criteria_intact()


@pytest.mark.req("HYP-010")
def test_an_unconventional_mechanism_faces_exactly_the_same_requirements():
    """There is no lenient path and no strict one.

    ``mechanism_class`` is a label with nothing to switch on, which is how the
    requirement is met: a validator that treated the two differently would need
    a branch, and there is none to write.
    """
    conventional = _hypothesis(mechanism_class=MechanismClass.CONVENTIONAL)
    unconventional = _hypothesis(mechanism_class=MechanismClass.UNCONVENTIONAL)

    assert conventional.criteria_hash() == unconventional.criteria_hash()
    for mechanism_class in MechanismClass:
        with pytest.raises(ConfigurationError, match="HYP-008"):
            _hypothesis(mechanism_class=mechanism_class, cannot_distinguish_from=())
        with pytest.raises(ConfigurationError, match="expected sign unspecified"):
            _hypothesis(
                mechanism_class=mechanism_class,
                expected=replace(SIGNATURE, sign="unspecified"),
            )


@pytest.mark.req("HYP-001")
def test_the_registry_finds_the_hypotheses_behind_a_claim():
    registry = HypothesisRegistry()
    registry.register(_hypothesis())
    registry.register(_hypothesis(record_id="HYP-TEST-002", claim_ids=("CLM-TEST-001", "CLM-2")))
    assert len(registry.for_claim("CLM-TEST-001")) == 2
    assert len(registry.for_claim("CLM-2")) == 1
    assert registry.for_claim("CLM-ABSENT") == []
