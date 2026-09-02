"""Withheld-data evaluation path.

Amendment A4.5, requirements INF-005, INF-006, INF-008, HYP-009, and KER-008.

The evaluation path is separated from the evolution path by construction, not by
convention.  A ``WithheldDataset`` will not hand out its values at all until a
``FreezeRecord`` exists pinning the source revision, configuration and claim
criteria.  A module running inside a step therefore cannot read it, whatever the
author intended -- which is the only way INF-008 can mean anything, since the
person most likely to peek is the person who wrote the model.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from dataclasses import field as dc_field
from datetime import UTC, datetime

import numpy as np

from mfespf.core.errors import ProvenanceError
from mfespf.core.hashing import digest
from mfespf.registries.claims import ClaimRecord, ClaimStatus
from mfespf.registries.datasets import DatasetRegistry

__all__ = ["EvaluationResult", "FreezeRecord", "WithheldDataset", "WithheldEvaluator"]


@dataclass(frozen=True)
class FreezeRecord:
    """A4.5 / INF-005, INF-006: what was fixed, and when, before evaluation."""

    claim_id: str
    information_cutoff: str
    source_revision: str
    configuration_hash: str
    parameter_set_hash: str
    criteria_hash: str
    frozen_at: str = ""
    note: str = ""

    def __post_init__(self) -> None:
        if not self.frozen_at:
            object.__setattr__(self, "frozen_at", datetime.now(UTC).isoformat(timespec="seconds"))

    def fingerprint(self) -> str:
        return digest(self, length=16)

    def describe(self) -> str:
        return (
            f"claim {self.claim_id} frozen at {self.frozen_at}\n"
            f"  information cutoff:  {self.information_cutoff}\n"
            f"  source revision:     {self.source_revision}\n"
            f"  configuration hash:  {self.configuration_hash}\n"
            f"  parameter set hash:  {self.parameter_set_hash}\n"
            f"  criteria hash:       {self.criteria_hash}\n"
            f"  fingerprint:         {self.fingerprint()}"
        )


class WithheldDataset:
    """Observations sealed against the evolution and calibration path (INF-008).

    The values live behind a loader that is not called until ``release`` is
    invoked with a matching freeze record.  Attempting to read them earlier
    raises rather than returning anything -- including anything empty, which a
    caller might otherwise treat as a legitimate absence of data.
    """

    def __init__(
        self,
        dataset_id: str,
        version: str,
        loader: Callable[[], np.ndarray],
        *,
        claim_id: str,
        variables: Sequence[str] = (),
        description: str = "",
    ) -> None:
        self.dataset_id = dataset_id
        self.version = version
        self.claim_id = claim_id
        self.variables = tuple(variables)
        self.description = description
        self._loader = loader
        self._released_under: FreezeRecord | None = None
        self._access_log: list[str] = []

    @property
    def is_released(self) -> bool:
        return self._released_under is not None

    def release(self, freeze: FreezeRecord) -> None:
        if freeze.claim_id != self.claim_id:
            raise ProvenanceError(
                f"freeze record is for claim {freeze.claim_id!r}, but dataset "
                f"{self.dataset_id!r} is withheld for claim {self.claim_id!r}"
            )
        self._released_under = freeze
        self._access_log.append(
            f"released under freeze {freeze.fingerprint()} at {freeze.frozen_at}"
        )

    def values(self) -> np.ndarray:
        if self._released_under is None:
            raise ProvenanceError(
                f"withheld dataset {self.dataset_id!r} was accessed before its claim was frozen. "
                f"INF-008 keeps withheld data out of the evolution and calibration path; it is "
                f"readable only through the post-run evaluator."
            )
        self._access_log.append("values read")
        return np.asarray(self._loader(), dtype=float)

    @property
    def access_log(self) -> tuple[str, ...]:
        return tuple(self._access_log)


@dataclass
class EvaluationResult:
    """The outcome of comparing a frozen prediction with withheld observations."""

    claim_id: str
    freeze: FreezeRecord
    metric_name: str
    metric_value: float
    threshold: float
    n_observations: int
    status: ClaimStatus
    predictive_interval_coverage: float | None = None
    nominal_coverage: float | None = None
    notes: list[str] = dc_field(default_factory=list)

    @property
    def passed(self) -> bool:
        return self.status is ClaimStatus.SUPPORTED_WITHIN_SCOPE

    def describe(self) -> str:
        lines = [
            f"claim {self.claim_id}: {self.status.value}",
            f"  {self.metric_name} = {self.metric_value:.6g} (threshold {self.threshold:.6g})",
            f"  observations: {self.n_observations}",
            f"  freeze fingerprint: {self.freeze.fingerprint()}",
        ]
        if self.predictive_interval_coverage is not None:
            lines.append(
                f"  predictive interval coverage: {self.predictive_interval_coverage:.3f} "
                f"(nominal {self.nominal_coverage:.3f})"
                if self.nominal_coverage is not None
                else f"  predictive interval coverage: {self.predictive_interval_coverage:.3f}"
            )
        lines.extend(f"  note: {n}" for n in self.notes)
        return "\n".join(lines)


class WithheldEvaluator:
    """The only path by which withheld data may reach a claim (KER-008)."""

    def __init__(self, datasets: DatasetRegistry | None = None) -> None:
        self.datasets = datasets
        self.freezes: dict[str, FreezeRecord] = {}

    def freeze(self, claim: ClaimRecord, freeze: FreezeRecord) -> FreezeRecord:
        """Pin a claim's criteria and configuration before any evaluation."""
        claim.assert_criteria_intact()
        if freeze.criteria_hash != claim.criteria_hash():
            raise ProvenanceError(
                f"freeze record criteria hash {freeze.criteria_hash} does not match claim "
                f"{claim.record_id!r} ({claim.criteria_hash()}); the thresholds moved between "
                f"registration and freeze (HYP-009)"
            )
        self.freezes[claim.record_id] = freeze
        return freeze

    def evaluate(
        self,
        *,
        claim: ClaimRecord,
        prediction: np.ndarray,
        withheld: WithheldDataset,
        metric_name: str = "root_mean_square_error",
        threshold: float | None = None,
        predictive_interval: tuple[np.ndarray, np.ndarray] | None = None,
        nominal_coverage: float | None = None,
    ) -> EvaluationResult:
        """Compare a frozen prediction with withheld observations exactly once.

        Every guard runs before the data is touched: the claim must be frozen,
        the freeze must match the claim's current criteria, and the provenance
        registry must show no contamination for this claim.  Only then is the
        dataset released.
        """
        freeze = self.freezes.get(claim.record_id)
        if freeze is None:
            raise ProvenanceError(
                f"claim {claim.record_id!r} has no freeze record; INF-006 requires the "
                f"configuration, parameter and criteria hashes be recorded before evaluation "
                f"data are introduced"
            )
        claim.assert_criteria_intact()
        if freeze.criteria_hash != claim.criteria_hash():
            raise ProvenanceError(
                f"claim {claim.record_id!r} criteria changed after freeze; register a new claim "
                f"version (HYP-009)"
            )
        if self.datasets is not None:
            self.datasets.assert_clean(claim.record_id)
            roles = {
                u.role
                for u in self.datasets.uses(
                    claim_id=claim.record_id, dataset_id=withheld.dataset_id
                )
            }
            if not any(r.claims_independence for r in roles):
                raise ProvenanceError(
                    f"dataset {withheld.dataset_id!r} is not declared as withheld evidence for "
                    f"claim {claim.record_id!r}; declare its role before evaluating (INF-001)"
                )

        withheld.release(freeze)
        observations = withheld.values()
        predicted = np.asarray(prediction, dtype=float)
        if predicted.shape != observations.shape:
            raise ProvenanceError(
                f"prediction shape {predicted.shape} does not match withheld observations "
                f"{observations.shape}; reshaping here would be an undeclared transformation"
            )

        residual = predicted - observations
        metric = float(np.sqrt(np.mean(residual**2)))
        limit = threshold if threshold is not None else float("inf")

        coverage = None
        notes: list[str] = []
        if predictive_interval is not None:
            lower, upper = (np.asarray(a, dtype=float) for a in predictive_interval)
            inside = np.logical_and(observations >= lower, observations <= upper)
            coverage = float(np.mean(inside))
            if nominal_coverage is not None and observations.size < 20:
                notes.append(
                    f"coverage estimated from {observations.size} observations; too few for a "
                    "meaningful test of interval calibration (UNC-008)"
                )

        if threshold is None:
            status = ClaimStatus.EMPIRICALLY_UNRESOLVED
            notes.append("no failure threshold supplied, so the comparison cannot decide the claim")
        elif metric <= limit:
            status = ClaimStatus.SUPPORTED_WITHIN_SCOPE
            notes.append(
                "supported within the tested configuration and validity domain only (CLM-003)"
            )
        else:
            status = ClaimStatus.REFUTED

        return EvaluationResult(
            claim_id=claim.record_id,
            freeze=freeze,
            metric_name=metric_name,
            metric_value=metric,
            threshold=limit,
            n_observations=int(observations.size),
            status=status,
            predictive_interval_coverage=coverage,
            nominal_coverage=nominal_coverage,
            notes=notes,
        )
