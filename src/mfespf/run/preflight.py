"""Resource preflight estimator.

Amendment A5.4 and requirements CMP-003 through CMP-005, plus KER-007.

The estimator sums module cost declarations, scales them by the configuration's
actual dimensions, and compares the total against the declared budget.  It runs
before initialisation so that a configuration which cannot fit the declared
resource class is refused rather than discovered halfway through.

The estimate is only as good as the declarations behind it, which is why
``CostDeclaration.benchmark_configuration`` is mandatory: a module that has never
been measured contributes a guess, and the report says which modules those are.

The memory that is actually spent
---------------------------------

For a long time this estimator's memory figure was wrong by a factor of twenty,
in the direction that admits a run which then exhausts memory -- and that run
happened: a 120000-step sweep reached 4.3 GB and had to be killed. Summing the
modules' declared working memory answered a question nobody was asking. What
grows is the *retained trajectory*: the driver keeps every step's published
fields so the manifest can be audited, and that is linear in the step count while
the working set is not.

Measuring it produced a more interesting number than expected. Across six
configurations spanning 4 to 24 layers and 8 to 32 bands, the retained cost fits

    bytes per step = 714 * n_fields + 17.6 * n_values

to within 8.8 percent. **The per-field term is forty times the per-value term.**
What is being stored is not the physics -- 49 fields of the moist column carry
1.7 kB of floating-point numbers per step out of 37.7 kB retained. The rest is
the metadata each ``Field`` carries: its units, quantity kind, time semantics,
spatial support, provenance, uncertainty, sign convention and missing-value
policy. That metadata is what this framework is for, and it costs about 700 bytes
per field per step, so it is worth knowing that it is the memory profile.

The field count needs no declaration: it is ``len(produced_outputs)`` summed over
the participating modules, which matches the retained count exactly on every
configuration tested. The value count does not follow from the contracts, so the
preflight omits it and says so -- which leaves the estimate low by 9 to 19
percent rather than by a factor of twenty, and the driver records the exact
figure once the fields exist.
"""

from __future__ import annotations

from mfespf.contracts.cost import CostEstimate, ResourceBudget
from mfespf.core.errors import ResourceBudgetError
from mfespf.coupling.coupler import CouplingSchedule

__all__ = [
    "RETAINED_BYTES_PER_FIELD",
    "RETAINED_BYTES_PER_VALUE",
    "assert_within_budget",
    "estimate_cost",
    "retained_trajectory_bytes",
]

#: Measured, by least squares over six configurations from 4 to 24 layers and 8
#: to 32 spectral bands, against ``tracemalloc`` peak divided by step count. Not
#: a guess at ``sys.getsizeof``: a ``Field`` is a frozen dataclass holding two
#: further dataclasses, several enums and several strings, and the interned and
#: shared parts of that are exactly what a structural estimate would get wrong.
RETAINED_BYTES_PER_FIELD = 714.0

#: The per-element term, forty times smaller. A float64 is 8 bytes; the rest is
#: the amortised numpy array header.
RETAINED_BYTES_PER_VALUE = 17.6


def retained_trajectory_bytes(*, n_fields: int, n_steps: int, n_values: int = 0) -> int:
    """Memory the driver will hold to keep every step auditable.

    ``n_values`` defaults to zero because the preflight cannot know it, and the
    term is 4 percent of the total at 8 layers. Pass it once the fields exist to
    get the exact figure.
    """
    per_step = RETAINED_BYTES_PER_FIELD * n_fields + RETAINED_BYTES_PER_VALUE * n_values
    return int(per_step * n_steps)


def estimate_cost(
    *,
    schedule: CouplingSchedule,
    n_steps: int,
    budget: ResourceBudget,
    spectral_bins: int = 1,
    spatial_cells: int = 1,
    ensemble_members: int = 1,
    output_bytes_per_step: int = 0,
    retained_values_per_step: int = 0,
) -> CostEstimate:
    """CMP-004: estimate memory, work, storage and problem size for a configuration.

    ``retained_values_per_step`` refines the retained-trajectory term when the
    caller knows how many numbers the run publishes per step. The preflight does
    not, so it defaults to zero and the estimate comes out 9 to 19 percent low
    rather than the factor of twenty it was low by before the term existed.
    """
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

    retained_fields = 0

    for module_schedule in schedule.ordered():
        module = module_schedule.module
        if not module.status.participates:
            continue
        retained_fields += len(module.contract.produced_outputs)
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

    # The retained trajectory, which is the memory that actually grows. Counted
    # after the module loop because it is a property of the run rather than of
    # any module: every step's published fields are kept so the manifest can be
    # audited, and that is linear in the step count while the working set is not.
    estimate.retained_trajectory_bytes = retained_trajectory_bytes(
        n_fields=retained_fields, n_steps=n_steps, n_values=retained_values_per_step
    )
    estimate.peak_memory_bytes += estimate.retained_trajectory_bytes * ensemble_members
    if retained_values_per_step == 0 and retained_fields:
        estimate.warnings.append(
            "the retained-trajectory estimate omits its per-value term, which is "
            "about 4 percent of the total at 8 layers and 19 percent at 24; the "
            "driver records the exact figure once the fields exist"
        )

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
