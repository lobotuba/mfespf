"""Computational cost declarations and resource budgets.

Amendment A5 and requirements CMP-001 through CMP-009, plus SYS-005.

The charter treats architectural completeness and executable complexity as
separate axes.  That only stays true if the cost of a configuration can be
estimated before it runs, which is what these declarations are for: a module
says how it scales, the preflight estimator adds the declarations up, and a
configuration that would not fit the declared resource class fails or warns
before consuming the workstation (CMP-005).
"""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field as dc_field
from enum import StrEnum

__all__ = [
    "CostDeclaration",
    "CostEstimate",
    "ParallelismModel",
    "ResourceBudget",
    "ResourceClass",
]


class ResourceClass(StrEnum):
    """A5.3 / charter 17.2, ordered from smallest to largest."""

    LAPTOP = "laptop"
    ENGINEERING_WORKSTATION = "engineering_workstation"
    MULTICORE_SERVER = "multicore_server"
    DISTRIBUTED_ENVIRONMENT = "distributed_environment"
    SPECIALIZED_HIGH_PERFORMANCE_ENVIRONMENT = "specialized_high_performance_environment"

    @property
    def rank(self) -> int:
        return list(ResourceClass).index(self)

    def satisfies(self, required: ResourceClass) -> bool:
        return self.rank >= required.rank


class ParallelismModel(StrEnum):
    SERIAL = "serial"
    THREADED = "threaded"
    PROCESS_PARALLEL = "process_parallel"
    DISTRIBUTED = "distributed"
    ACCELERATOR = "accelerator"


@dataclass(frozen=True)
class CostDeclaration:
    """A5.2 / MOD-013, CMP-001, CMP-002.

    ``work_units_per_step`` is a benchmark-anchored estimate in an arbitrary but
    consistent unit; what matters is that it is measured once for a declared
    configuration (CMP-002) and scaled by the declared scaling variables, so
    estimates stay comparable across modules.
    """

    module_id: str
    time_complexity: str
    memory_complexity: str
    scaling_variables: tuple[str, ...] = ()
    work_units_per_step: float = 1.0
    benchmark_configuration: str = ""
    peak_memory_bytes: int = 0
    temporary_storage_bytes: int = 0
    persistent_output_bytes_per_step: int = 0
    parallelism: ParallelismModel = ParallelismModel.SERIAL
    accelerator_required: bool = False
    supported_resource_classes: tuple[ResourceClass, ...] = (ResourceClass.ENGINEERING_WORKSTATION,)
    lower_fidelity_alternatives: tuple[str, ...] = ()
    checkpoint_cost_bytes: int = 0
    ensemble_cost_multiplier: float = 1.0

    def scaled(self, factors: dict[str, float]) -> CostDeclaration:
        """Scale the benchmark estimate by the declared scaling variables.

        Only variables the module declared in ``scaling_variables`` are applied;
        an undeclared factor is ignored rather than guessed at, which keeps the
        estimate traceable to the declaration.
        """
        multiplier = 1.0
        for name in self.scaling_variables:
            multiplier *= float(factors.get(name, 1.0))
        return CostDeclaration(
            **{
                **self.__dict__,
                "work_units_per_step": self.work_units_per_step * multiplier,
                "peak_memory_bytes": int(self.peak_memory_bytes * multiplier),
            }
        )


@dataclass(frozen=True)
class ResourceBudget:
    """The per-run ceiling declared by CMP-003 and A5.4."""

    resource_class: ResourceClass = ResourceClass.ENGINEERING_WORKSTATION
    max_memory_bytes: int = 16 * 1024**3
    max_work_units: float = 1.0e9
    max_output_bytes: int = 4 * 1024**3
    max_spectral_bins: int = 100_000
    max_spatial_cells: int = 1_000_000
    max_ensemble_members: int = 1000
    on_exceed: str = "fail"  # "fail" or "warn" (CMP-005 leaves this to policy)

    def __post_init__(self) -> None:
        if self.on_exceed not in ("fail", "warn"):
            raise ValueError("on_exceed must be 'fail' or 'warn'")


@dataclass
class CostEstimate:
    """The preflight result (CMP-004)."""

    total_work_units: float = 0.0
    peak_memory_bytes: int = 0
    retained_trajectory_bytes: int = 0
    """The part of ``peak_memory_bytes`` that grows with the step count.

    Reported separately because it behaves differently from the rest: the
    modules' working set is fixed and this is linear in ``n_steps``, so a
    configuration that fits at 2000 steps and not at 200000 differs only here.
    Kept out of the budget checks as a separate line -- it is already inside
    ``peak_memory_bytes`` -- and printed on its own by ``describe``.
    """
    output_bytes: int = 0
    spectral_bins: int = 0
    spatial_cells: int = 0
    ensemble_members: int = 1
    per_module: dict[str, float] = dc_field(default_factory=dict)
    warnings: list[str] = dc_field(default_factory=list)
    violations: list[str] = dc_field(default_factory=list)

    @property
    def within_budget(self) -> bool:
        return not self.violations

    def describe(self) -> str:
        lines = [
            f"work units:    {self.total_work_units:.4g}",
            f"peak memory:   {self.peak_memory_bytes / 1024**3:.3f} GiB",
            f"  of which retained trajectory: {self.retained_trajectory_bytes / 1024**3:.3f} GiB",
            f"output:        {self.output_bytes / 1024**2:.1f} MiB",
            f"spectral bins: {self.spectral_bins}",
            f"spatial cells: {self.spatial_cells}",
            f"ensemble:      {self.ensemble_members}",
        ]
        if self.per_module:
            lines.append("per module:")
            for name, work in sorted(self.per_module.items(), key=lambda kv: -kv[1]):
                share = work / self.total_work_units if self.total_work_units else 0.0
                lines.append(f"  {name:<28} {work:12.4g}  ({share:5.1%})")
        for w in self.warnings:
            lines.append(f"WARNING: {w}")
        for v in self.violations:
            lines.append(f"OVER BUDGET: {v}")
        return "\n".join(lines)
