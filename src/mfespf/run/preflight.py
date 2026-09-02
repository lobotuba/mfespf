"""Resource preflight estimator.

Amendment A5.4 and requirements CMP-003 through CMP-005, plus KER-007.

The estimator sums module cost declarations, scales them by the configuration's
actual dimensions, and compares the total against the declared budget.  It runs
before initialisation so that a configuration which cannot fit the declared
resource class is refused rather than discovered halfway through.

The estimate is only as good as the declarations behind it, which is why
``CostDeclaration.benchmark_configuration`` is mandatory: a module that has never
been measured contributes a guess, and the report says which modules those are.
"""

from __future__ import annotations

from mfespf.contracts.cost import CostEstimate, ResourceBudget
from mfespf.core.errors import ResourceBudgetError
from mfespf.coupling.coupler import CouplingSchedule

__all__ = ["assert_within_budget", "estimate_cost"]


def estimate_cost(
    *,
    schedule: CouplingSchedule,
    n_steps: int,
    budget: ResourceBudget,
    spectral_bins: int = 1,
    spatial_cells: int = 1,
    ensemble_members: int = 1,
    output_bytes_per_step: int = 0,
) -> CostEstimate:
    """CMP-004: estimate memory, work, storage and problem size for a configuration."""
    estimate = CostEstimate(
        spectral_bins=spectral_bins,
        spatial_cells=spatial_cells,
        ensemble_members=ensemble_members,
    )
    scaling = {
        "n_bands": float(spectral_bins),
        "n_cells": float(spatial_cells),
        "n_layers": 1.0,
    }

    for module_schedule in schedule.ordered():
        module = module_schedule.module
        if not module.status.participates:
            continue
        cost = module.contract.cost
        if cost is None:  # pragma: no cover - blocked by ModuleContract.validate
            estimate.warnings.append(
                f"{module.module_id} has no cost declaration; its contribution is unknown"
            )
            continue
        if not cost.benchmark_configuration:
            estimate.warnings.append(
                f"{module.module_id} cost is an estimate, not a benchmark (CMP-002)"
            )
        scaled = cost.scaled(scaling)
        invocations = n_steps / module_schedule.step_multiple
        work = scaled.work_units_per_step * invocations * ensemble_members
        estimate.per_module[module.module_id] = work
        estimate.total_work_units += work
        estimate.peak_memory_bytes += scaled.peak_memory_bytes * ensemble_members
        estimate.output_bytes += scaled.persistent_output_bytes_per_step * int(invocations)

        if not any(
            supported.satisfies(budget.resource_class) or budget.resource_class.satisfies(supported)
            for supported in cost.supported_resource_classes
        ):  # pragma: no cover - defensive
            estimate.warnings.append(
                f"{module.module_id} declares no support for resource class "
                f"{budget.resource_class.value}"
            )
        if cost.accelerator_required and budget.resource_class.rank < 4:
            estimate.violations.append(
                f"{module.module_id} requires an accelerator, which resource class "
                f"{budget.resource_class.value} does not provide (CMP-008)"
            )

    estimate.output_bytes += output_bytes_per_step * n_steps

    checks = (
        ("work units", estimate.total_work_units, budget.max_work_units),
        ("peak memory", float(estimate.peak_memory_bytes), float(budget.max_memory_bytes)),
        ("output bytes", float(estimate.output_bytes), float(budget.max_output_bytes)),
        ("spectral bins", float(spectral_bins), float(budget.max_spectral_bins)),
        ("spatial cells", float(spatial_cells), float(budget.max_spatial_cells)),
        ("ensemble members", float(ensemble_members), float(budget.max_ensemble_members)),
    )
    for name, value, limit in checks:
        if value > limit:
            estimate.violations.append(f"{name}: {value:.6g} exceeds budget {limit:.6g}")

    return estimate


def assert_within_budget(estimate: CostEstimate, budget: ResourceBudget) -> CostEstimate:
    """Apply the CMP-005 policy: fail or warn according to the budget's declaration."""
    if estimate.within_budget:
        return estimate
    message = "resource preflight exceeded the declared budget:\n  - " + "\n  - ".join(
        estimate.violations
    )
    if budget.on_exceed == "fail":
        raise ResourceBudgetError(message)
    estimate.warnings.append(message)
    return estimate
