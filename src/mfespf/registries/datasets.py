"""Dataset provenance and information-separation control.

Amendment A4 and requirements INF-001 through INF-008, DAT-008, DAT-009.

This registry tracks evidentiary information flow, not physical transfer.  Its
job is to make calibration/validation contamination a machine-detectable
condition instead of a matter of authorial discipline, because the failure it
guards against -- using a record to shape the model and then reporting the same
record as independent evidence -- is invisible in the physics and fatal to the
claim.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

from mfespf.core.errors import ConfigurationError, ProvenanceError
from mfespf.core.hashing import file_digest
from mfespf.core.versioning import Version
from mfespf.registries.base import Registry

__all__ = [
    "ContaminationFinding",
    "DatasetRecord",
    "DatasetRegistry",
    "DatasetRole",
    "DatasetUse",
    "TransformationStep",
]


class DatasetRole(StrEnum):
    """A4.2 / INF-002."""

    FUNDAMENTAL_OR_REFERENCE_CONSTANT = "fundamental_or_reference_constant"
    MATERIAL_PROPERTY = "material_property"
    INITIAL_CONDITION = "initial_condition"
    BOUNDARY_CONDITION = "boundary_condition"
    PARAMETER_CALIBRATION = "parameter_calibration"
    PROCESS_VALIDATION = "process_validation"
    HISTORICAL_RECONSTRUCTION_EVALUATION = "historical_reconstruction_evaluation"
    WITHHELD_HINDCAST_EVALUATION = "withheld_hindcast_evaluation"
    OUT_OF_SAMPLE_PREDICTION_EVALUATION = "out_of_sample_prediction_evaluation"
    DIAGNOSTIC_ONLY = "diagnostic_only"

    @property
    def shapes_the_model(self) -> bool:
        """Roles in which the data influenced what the model *is*."""
        return self in (
            DatasetRole.PARAMETER_CALIBRATION,
            DatasetRole.INITIAL_CONDITION,
            DatasetRole.BOUNDARY_CONDITION,
        )

    @property
    def is_evidence(self) -> bool:
        """Roles in which the data is offered as a test of the model."""
        return self in (
            DatasetRole.PROCESS_VALIDATION,
            DatasetRole.HISTORICAL_RECONSTRUCTION_EVALUATION,
            DatasetRole.WITHHELD_HINDCAST_EVALUATION,
            DatasetRole.OUT_OF_SAMPLE_PREDICTION_EVALUATION,
        )

    @property
    def claims_independence(self) -> bool:
        """Roles that assert the data did not touch model construction."""
        return self in (
            DatasetRole.WITHHELD_HINDCAST_EVALUATION,
            DatasetRole.OUT_OF_SAMPLE_PREDICTION_EVALUATION,
        )


@dataclass(frozen=True, slots=True)
class TransformationStep:
    """One reproducible preprocessing operation (DAT-008, INF-007)."""

    operation: str
    parameters: str = ""
    software: str = ""
    uncertainty_effect: str = "unresolved"


@dataclass(frozen=True)
class DatasetRecord:
    """A4.3 / DAT-009.  Identity and provenance of an input dataset."""

    record_id: str
    version: Version
    title: str
    source: str
    checksum: str
    variables: tuple[str, ...] = ()
    units: tuple[str, ...] = ()
    temporal_coverage: tuple[str, str] | None = None
    spatial_coverage: str = ""
    uncertainty_record: str = "unresolved"
    quality_flags: tuple[str, ...] = ()
    acquisition_record: str = ""
    transformations: tuple[TransformationStep, ...] = ()
    path: str | None = None

    def __post_init__(self) -> None:
        if not self.checksum:
            raise ConfigurationError(
                f"dataset {self.record_id!r} has no checksum; a dataset without one cannot "
                f"support a reproducible result (DAT-009)"
            )
        if self.units and len(self.units) != len(self.variables):
            raise ConfigurationError(
                f"dataset {self.record_id!r} declares {len(self.variables)} variables but "
                f"{len(self.units)} units"
            )

    @property
    def dataset_id(self) -> str:
        return self.record_id

    @classmethod
    def from_file(
        cls, record_id: str, version: Version, path: str | Path, **kwargs: object
    ) -> DatasetRecord:
        """Build a record with the checksum computed from the file on disk."""
        return cls(
            record_id=record_id,
            version=version,
            checksum=file_digest(path),
            path=str(Path(path).as_posix()),
            **kwargs,  # type: ignore[arg-type]
        )


@dataclass(frozen=True, slots=True)
class DatasetUse:
    """One assignment of a dataset to a role, for a claim (INF-001).

    Recorded per use rather than per dataset: the same record may legitimately
    serve as a boundary condition for one claim and as withheld evidence for
    another.  What it may not do is both, for the same claim.
    """

    dataset_id: str
    dataset_version: str
    role: DatasetRole
    claim_id: str
    variables: tuple[str, ...] = ()
    run_id: str | None = None
    influenced: tuple[str, ...] = ()  # modules, equations or parameters affected
    note: str = ""

    @property
    def key(self) -> tuple[str, str]:
        return (self.dataset_id, self.dataset_version)


@dataclass(frozen=True, slots=True)
class ContaminationFinding:
    claim_id: str
    dataset_id: str
    variables: tuple[str, ...]
    shaping_role: DatasetRole
    evidence_role: DatasetRole
    detail: str

    def describe(self) -> str:
        return (
            f"claim {self.claim_id}: dataset {self.dataset_id} variables "
            f"{list(self.variables)} were used as {self.shaping_role.value} and are also offered "
            f"as {self.evidence_role.value} -- {self.detail}"
        )


class DatasetRegistry(Registry[DatasetRecord]):
    """Dataset records plus the use log that INF-003 and INF-004 are checked against."""

    def __init__(self) -> None:
        super().__init__("dataset")
        self._uses: list[DatasetUse] = []

    def declare_use(self, use: DatasetUse) -> DatasetUse:
        """Record an evidentiary use.  The dataset must already be registered."""
        self.entry(use.dataset_id, use.dataset_version)
        self._uses.append(use)
        return use

    def uses(
        self, *, claim_id: str | None = None, dataset_id: str | None = None
    ) -> list[DatasetUse]:
        return [
            u
            for u in self._uses
            if (claim_id is None or u.claim_id == claim_id)
            and (dataset_id is None or u.dataset_id == dataset_id)
        ]

    def check_contamination(self, claim_id: str) -> list[ContaminationFinding]:
        """Find data offered as evidence for a claim it also helped construct.

        A4.4 and INF-003.  Overlap is evaluated at variable granularity: using a
        reanalysis wind field as a boundary condition does not contaminate a
        claim evaluated against that same dataset's independent radiation
        variables, but reusing the *same* variables does.
        """
        findings: list[ContaminationFinding] = []
        by_dataset: dict[str, list[DatasetUse]] = defaultdict(list)
        for use in self.uses(claim_id=claim_id):
            by_dataset[use.dataset_id].append(use)

        for dataset_id, uses in sorted(by_dataset.items()):
            shaping = [u for u in uses if u.role.shapes_the_model]
            evidence = [u for u in uses if u.role.is_evidence]
            for s in shaping:
                for e in evidence:
                    overlap = self._overlap(s, e, dataset_id)
                    if overlap is None:
                        continue
                    if s.role is DatasetRole.PARAMETER_CALIBRATION:
                        detail = (
                            "calibration data cannot be reported as independent validation "
                            "(VAL-006)"
                        )
                    elif e.role.claims_independence:
                        detail = (
                            "a prescribed input cannot simultaneously be withheld out-of-sample "
                            "evidence (INF-008)"
                        )
                    else:
                        detail = (
                            "the same variables were prescribed to the run and then evaluated "
                            "against it, so agreement is partly definitional (A4.4)"
                        )
                    findings.append(
                        ContaminationFinding(
                            claim_id=claim_id,
                            dataset_id=dataset_id,
                            variables=overlap,
                            shaping_role=s.role,
                            evidence_role=e.role,
                            detail=detail,
                        )
                    )
        return findings

    @staticmethod
    def _overlap(a: DatasetUse, b: DatasetUse, dataset_id: str) -> tuple[str, ...] | None:
        """Variables common to two uses, or all of them when neither narrows scope."""
        if not a.variables or not b.variables:
            return ("<whole dataset>",)
        common = tuple(sorted(set(a.variables) & set(b.variables)))
        return common or None

    def assert_clean(self, claim_id: str) -> None:
        """Raise on any contamination finding for a claim (ACC-008)."""
        findings = self.check_contamination(claim_id)
        if findings:
            raise ProvenanceError(
                "information-provenance validation failed:\n  - "
                + "\n  - ".join(f.describe() for f in findings)
            )

    def independence_statement(self, claim_id: str) -> dict[str, list[str]]:
        """The INF-004 information-independence statement for a claim.

        Answers A4.6's question -- what was available at each stage -- from the
        recorded uses rather than from recollection.
        """
        statement: dict[str, list[str]] = {
            "mechanism_and_equation_selection": [],
            "parameter_determination": [],
            "configuration_and_boundary_setup": [],
            "final_evaluation": [],
            "explicitly_withheld": [],
        }
        for use in sorted(self.uses(claim_id=claim_id), key=lambda u: (u.dataset_id, u.role.value)):
            label = f"{use.dataset_id}@{use.dataset_version} {list(use.variables) or '[all]'}"
            if use.role is DatasetRole.PARAMETER_CALIBRATION:
                statement["parameter_determination"].append(label)
            elif use.role in (DatasetRole.INITIAL_CONDITION, DatasetRole.BOUNDARY_CONDITION):
                statement["configuration_and_boundary_setup"].append(label)
            elif use.role in (
                DatasetRole.MATERIAL_PROPERTY,
                DatasetRole.FUNDAMENTAL_OR_REFERENCE_CONSTANT,
            ):
                statement["mechanism_and_equation_selection"].append(label)
            elif use.role.claims_independence:
                statement["explicitly_withheld"].append(label)
                statement["final_evaluation"].append(label)
            elif use.role.is_evidence:
                statement["final_evaluation"].append(label)
        return statement
