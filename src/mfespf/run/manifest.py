"""Run manifest and lineage.

Charter section 16, amendment A8, requirements AUD-001 through AUD-008, and
KER-006.

A manifest is the object that makes a result reproducible and a claim auditable.
It records what was active, at which versions, from which data, under which
criteria -- and, once a run has produced output, the checksums of that output.
A8.4 forbids overwriting retained evidence, so a manifest written to disk refuses
to replace a different manifest already at that path.
"""

from __future__ import annotations

import getpass
import platform
import socket
import sys
from dataclasses import dataclass
from dataclasses import field as dc_field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from mfespf.contracts.cost import ResourceClass
from mfespf.core.errors import ProvenanceError
from mfespf.core.hashing import canonical_json, digest, file_digest
from mfespf.run.modes import OperatingMode

__all__ = [
    "DifferenceReport",
    "EnvironmentRecord",
    "RunManifest",
    "RunStatus",
    "diff_manifests",
    "lineage",
]

FRAMEWORK_VERSION = "0.10.0"


class RunStatus:
    """AUD-008 / SYS-008.  Deliberately a small set of strings, not an enum with
    room to grow: a run either may support inference or it may not."""

    PENDING = "pending"
    COMPLETE = "complete"
    FAILED = "failed"
    INVALID_FOR_SCIENTIFIC_INFERENCE = "invalid_for_scientific_inference"


@dataclass(frozen=True)
class EnvironmentRecord:
    """AUD-002: enough of the environment to explain a numerical difference."""

    python_version: str
    platform: str
    machine: str
    hostname: str
    user: str
    packages: dict[str, str] = dc_field(default_factory=dict)

    @classmethod
    def capture(cls, packages: dict[str, str] | None = None) -> EnvironmentRecord:
        if packages is None:
            packages = {}
            for name in ("numpy", "yaml"):
                try:
                    module = __import__(name)
                    packages[name] = getattr(module, "__version__", "unknown")
                except ImportError:  # pragma: no cover - optional dependency
                    packages[name] = "absent"
        try:
            user = getpass.getuser()
        except Exception:  # pragma: no cover - unusual environments
            user = "unknown"
        return cls(
            python_version=sys.version.split()[0],
            platform=platform.platform(),
            machine=platform.machine(),
            hostname=socket.gethostname(),
            user=user,
            packages=packages,
        )


@dataclass
class RunManifest:
    """Charter section 16 plus A8.2.

    ``configuration_hash`` covers everything that determines the run's numbers and
    excludes everything that does not (timestamps, host, output checksums).  Two
    runs with the same configuration hash and different results indicate
    non-determinism, which is a reportable failure rather than noise.
    """

    run_id: str
    purpose: str
    parent_run_ids: tuple[str, ...] = ()
    experiment_id: str | None = None
    claim_ids: tuple[str, ...] = ()
    framework_version: str = FRAMEWORK_VERSION
    source_revision: str = "unknown"
    execution_time: str = ""
    reason_for_change: str = ""
    operating_mode: OperatingMode = OperatingMode.CONTROLLED_EXPERIMENT
    """TIM-005.  What kind of run this is, in the sense evidence cares about.

    Deliberately outside ``configuration_payload``: two runs differing only in
    their declared mode produce identical numbers, and the configuration hash is
    for numbers. What the mode governs is which language a report may use about
    them, which ``assert_may_be_reported_as`` checks.
    """

    module_status: dict[str, str] = dc_field(default_factory=dict)
    module_implementations: dict[str, str] = dc_field(default_factory=dict)
    module_versions: dict[str, str] = dc_field(default_factory=dict)
    interface_versions: dict[str, str] = dc_field(default_factory=dict)
    equation_bindings: dict[str, str] = dc_field(default_factory=dict)
    state_ownership: dict[str, str] = dc_field(default_factory=dict)

    parameters: dict[str, Any] = dc_field(default_factory=dict)
    parameter_set_version: str = "unversioned"
    initial_conditions: dict[str, Any] = dc_field(default_factory=dict)
    boundary_conditions: dict[str, Any] = dc_field(default_factory=dict)
    datasets: dict[str, str] = dc_field(default_factory=dict)
    dataset_roles: dict[str, str] = dc_field(default_factory=dict)

    spatial_representation: str = "global_scalar"
    spatial_resolution: str = "0d"
    temporal_representation: str = "fixed_master_step"
    temporal_resolution: str = ""
    spectral_representation: str = ""
    spectral_resolution: str = ""
    coupling_schedule: dict[str, Any] = dc_field(default_factory=dict)
    solver_settings: dict[str, Any] = dc_field(default_factory=dict)
    tolerances: dict[str, float] = dc_field(default_factory=dict)
    random_seeds: dict[str, int] = dc_field(default_factory=dict)
    uncertainty_methods: dict[str, str] = dc_field(default_factory=dict)
    support_criterion: str = ""
    failure_criterion: str = ""

    resource_class: ResourceClass = ResourceClass.ENGINEERING_WORKSTATION
    resource_budget: dict[str, Any] = dc_field(default_factory=dict)
    environment: EnvironmentRecord | None = None

    dependency_graph: dict[str, Any] = dc_field(default_factory=dict)
    output_checksums: dict[str, str] = dc_field(default_factory=dict)
    warnings: list[str] = dc_field(default_factory=list)
    gate_results: dict[str, bool] = dc_field(default_factory=dict)
    status: str = RunStatus.PENDING
    invalidation_reasons: list[str] = dc_field(default_factory=list)

    def __post_init__(self) -> None:
        if not self.execution_time:
            self.execution_time = datetime.now(UTC).isoformat(timespec="seconds")
        if self.environment is None:
            self.environment = EnvironmentRecord.capture()

    # --- hashing ----------------------------------------------------------
    def configuration_payload(self) -> dict[str, Any]:
        """Exactly the fields that determine the run's numerical result."""
        return {
            "framework_version": self.framework_version,
            "source_revision": self.source_revision,
            "module_status": self.module_status,
            "module_implementations": self.module_implementations,
            "module_versions": self.module_versions,
            "interface_versions": self.interface_versions,
            "equation_bindings": self.equation_bindings,
            "state_ownership": self.state_ownership,
            "parameters": self.parameters,
            "parameter_set_version": self.parameter_set_version,
            "initial_conditions": self.initial_conditions,
            "boundary_conditions": self.boundary_conditions,
            "datasets": self.datasets,
            "spatial_resolution": self.spatial_resolution,
            "temporal_resolution": self.temporal_resolution,
            "spectral_resolution": self.spectral_resolution,
            "coupling_schedule": self.coupling_schedule,
            "solver_settings": self.solver_settings,
            "tolerances": self.tolerances,
            "random_seeds": self.random_seeds,
        }

    def configuration_hash(self) -> str:
        return digest(self.configuration_payload(), length=32)

    def criteria_hash(self) -> str:
        """Hash of the thresholds INF-006 requires be pinned before evaluation."""
        return digest(
            {
                "support_criterion": self.support_criterion,
                "failure_criterion": self.failure_criterion,
                "tolerances": self.tolerances,
                "uncertainty_methods": self.uncertainty_methods,
            },
            length=16,
        )

    def assert_may_be_reported_as(self, claimed: OperatingMode) -> None:
        """TIM-006 / VAL-007, checked against what the run actually was."""
        self.operating_mode.assert_may_be_reported_as(claimed)

    # --- status -----------------------------------------------------------
    def invalidate(self, reason: str) -> None:
        """SYS-008 / AUD-007: mark a run unusable for supporting or refuting a claim."""
        self.status = RunStatus.INVALID_FOR_SCIENTIFIC_INFERENCE
        self.invalidation_reasons.append(reason)

    @property
    def may_support_a_claim(self) -> bool:
        return self.status == RunStatus.COMPLETE and not self.invalidation_reasons

    def record_gate(self, name: str, passed: bool, *, mandatory: bool = True) -> None:
        self.gate_results[name] = passed
        if not passed and mandatory:
            self.invalidate(f"mandatory gate {name!r} failed")

    def record_output(self, path: str | Path) -> None:
        self.output_checksums[str(Path(path).name)] = file_digest(path)

    # --- serialisation ----------------------------------------------------
    def to_dict(self) -> dict[str, Any]:
        payload = {
            "run_id": self.run_id,
            "purpose": self.purpose,
            "parent_run_ids": list(self.parent_run_ids),
            "experiment_id": self.experiment_id,
            "claim_ids": list(self.claim_ids),
            "execution_time": self.execution_time,
            "reason_for_change": self.reason_for_change,
            "operating_mode": self.operating_mode.value,
            "status": self.status,
            "invalidation_reasons": list(self.invalidation_reasons),
            "configuration_hash": self.configuration_hash(),
            "criteria_hash": self.criteria_hash(),
            "configuration": self.configuration_payload(),
            "resource_class": self.resource_class.value,
            "resource_budget": self.resource_budget,
            "environment": self.environment,
            "dependency_graph": self.dependency_graph,
            "uncertainty_methods": self.uncertainty_methods,
            "support_criterion": self.support_criterion,
            "failure_criterion": self.failure_criterion,
            "dataset_roles": self.dataset_roles,
            "gate_results": self.gate_results,
            "output_checksums": self.output_checksums,
            "warnings": self.warnings,
        }
        return payload

    def to_json(self) -> str:
        return canonical_json(self.to_dict())

    def write(self, path: str | Path, *, allow_overwrite: bool = False) -> Path:
        """Write the manifest, refusing to overwrite different retained evidence.

        A8.4: a correction creates a new versioned record rather than replacing
        the old one.  Rewriting an identical manifest is permitted because it
        changes nothing; rewriting a different one is not.
        """
        target = Path(path)
        payload = self.to_json()
        if target.exists() and not allow_overwrite:
            existing = target.read_text(encoding="utf-8")
            if existing != payload:
                raise ProvenanceError(
                    f"{target} already holds a different manifest. Retained evidence is "
                    f"immutable (AUD-005, A8.4); write the correction to a new run_id and link "
                    f"it as a child instead."
                )
            return target
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(payload, encoding="utf-8")
        return target

    def describe(self) -> str:
        lines = [
            f"run {self.run_id} [{self.status}]",
            f"  purpose:      {self.purpose}",
            f"  framework:    {self.framework_version} @ {self.source_revision}",
            f"  config hash:  {self.configuration_hash()}",
            f"  claims:       {list(self.claim_ids) or 'none'}",
            f"  parents:      {list(self.parent_run_ids) or 'none'}",
            "  modules:",
        ]
        for module_id in sorted(self.module_status):
            lines.append(
                f"    {module_id:<12} {self.module_status[module_id]:<12} "
                f"{self.module_implementations.get(module_id, '')}"
            )
        if self.gate_results:
            lines.append("  gates:")
            for name, passed in sorted(self.gate_results.items()):
                lines.append(f"    {'PASS' if passed else 'FAIL'}  {name}")
        for reason in self.invalidation_reasons:
            lines.append(f"  INVALID: {reason}")
        return "\n".join(lines)


@dataclass
class DifferenceReport:
    """AUD-003, AUD-004, A8.3: what changed between a parent and a child run."""

    parent_run_id: str
    child_run_id: str
    changes: dict[str, list[str]] = dc_field(default_factory=dict)

    @property
    def is_empty(self) -> bool:
        return not any(self.changes.values())

    def to_dict(self) -> dict[str, Any]:
        return {
            "parent_run_id": self.parent_run_id,
            "child_run_id": self.child_run_id,
            "changes": {k: v for k, v in sorted(self.changes.items()) if v},
        }

    def describe(self) -> str:
        if self.is_empty:
            return (
                f"{self.child_run_id} is configurationally identical to {self.parent_run_id}; "
                f"any difference in results is non-determinism, not a change of setup"
            )
        lines = [f"differences from {self.parent_run_id} to {self.child_run_id}:"]
        for category, entries in sorted(self.changes.items()):
            if not entries:
                continue
            lines.append(f"  {category}:")
            lines.extend(f"    - {e}" for e in entries)
        return "\n".join(lines)


_DIFF_CATEGORIES = {
    "modules": ("module_status", "module_implementations", "module_versions"),
    "equations": ("equation_bindings",),
    "interfaces": ("interface_versions",),
    "state_ownership": ("state_ownership",),
    "parameters": ("parameters", "parameter_set_version"),
    "data": ("datasets", "dataset_roles"),
    "initial_and_boundary_conditions": ("initial_conditions", "boundary_conditions"),
    "resolution": (
        "spatial_resolution",
        "temporal_resolution",
        "spectral_resolution",
        "spatial_representation",
        "temporal_representation",
        "spectral_representation",
    ),
    "thresholds_and_criteria": (
        "tolerances",
        "support_criterion",
        "failure_criterion",
        "uncertainty_methods",
    ),
    "coupling": ("coupling_schedule", "solver_settings", "random_seeds"),
    "software_and_environment": ("framework_version", "source_revision", "environment"),
    "resources": ("resource_class", "resource_budget"),
}


def diff_manifests(parent: RunManifest, child: RunManifest) -> DifferenceReport:
    """Produce the machine-readable difference record AUD-003 requires."""
    report = DifferenceReport(parent_run_id=parent.run_id, child_run_id=child.run_id)
    for category, attributes in _DIFF_CATEGORIES.items():
        entries: list[str] = []
        for attribute in attributes:
            old = getattr(parent, attribute)
            new = getattr(child, attribute)
            if isinstance(old, dict) and isinstance(new, dict):
                entries.extend(_diff_mappings(attribute, old, new))
            elif old != new:
                entries.append(f"{attribute}: {_short(old)} -> {_short(new)}")
        report.changes[category] = entries
    return report


def _diff_mappings(label: str, old: dict[str, Any], new: dict[str, Any]) -> list[str]:
    entries = []
    for key in sorted(set(old) | set(new)):
        if key not in new:
            entries.append(f"{label}.{key} removed (was {_short(old[key])})")
        elif key not in old:
            entries.append(f"{label}.{key} added ({_short(new[key])})")
        elif old[key] != new[key]:
            entries.append(f"{label}.{key}: {_short(old[key])} -> {_short(new[key])}")
    return entries


def _short(value: Any, limit: int = 80) -> str:
    text = str(value)
    return text if len(text) <= limit else text[: limit - 3] + "..."


def lineage(manifest: RunManifest, known: dict[str, RunManifest]) -> list[RunManifest]:
    """The ancestry of a run, oldest first (KER-006, AUD-005).

    A manifest names its parents; this walks them into the chain a reader needs
    to answer "what did this result descend from". A parent that is not in
    ``known`` stops that branch rather than being invented, and a cycle raises:
    a run cannot be its own ancestor, and a lineage that says otherwise is
    corrupt rather than merely odd.
    """
    order: list[RunManifest] = []
    seen: set[str] = set()

    def walk(current: RunManifest, path: tuple[str, ...]) -> None:
        if current.run_id in path:
            raise ProvenanceError(
                f"run {current.run_id!r} appears in its own ancestry "
                f"({' -> '.join([*path, current.run_id])}); a lineage cannot contain a cycle"
            )
        for parent_id in current.parent_run_ids:
            parent = known.get(parent_id)
            if parent is not None:
                walk(parent, (*path, current.run_id))
        if current.run_id not in seen:
            seen.add(current.run_id)
            order.append(current)

    walk(manifest, ())
    return order
