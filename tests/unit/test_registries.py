"""Stage 0: registry immutability, equation categories, cost preflight."""

from __future__ import annotations

import pytest

from mfespf.contracts.cost import (
    CostDeclaration,
    ResourceBudget,
    ResourceClass,
)
from mfespf.core.errors import ConfigurationError, RegistryError, ResourceBudgetError
from mfespf.core.versioning import Version
from mfespf.physics.equations import EQUATION_IDS, build_equation_registry
from mfespf.registries.equations import EquationCategory, EquationRecord, ValidityDomain
from mfespf.run.configurations import airless_radiating_body
from mfespf.run.preflight import assert_within_budget, estimate_cost

pytestmark = pytest.mark.stage(0)

V1 = Version(1, 0, 0)
V2 = Version(2, 0, 0)


def _equation(record_id: str = "EQ-TEST", version: Version = V1, **kwargs) -> EquationRecord:
    defaults = {
        "name": "test equation",
        "category": EquationCategory.REDUCED_PHYSICAL_EQUATION,
        "mathematical_form": "y = a x",
        "physical_basis": "test",
        "correctness_tests": ("test_something",),
    }
    defaults.update(kwargs)
    return EquationRecord(record_id=record_id, version=version, **defaults)


@pytest.mark.req("EQR-001", "KER-001")
def test_the_stage1_equation_set_registers_and_binds():
    registry = build_equation_registry()
    for equation_id in EQUATION_IDS.values():
        record = registry.get(equation_id)
        assert record.mathematical_form
        assert record.physical_basis
    bindings = registry.bind({equation_id: "1.0.0" for equation_id in EQUATION_IDS.values()})
    assert len(bindings) == len(EQUATION_IDS)


@pytest.mark.req("EQR-006", "AUD-006")
def test_a_registered_version_cannot_be_silently_changed():
    registry = build_equation_registry()
    with pytest.raises(RegistryError, match="already registered with different content"):
        registry.register(_equation(EQUATION_IDS["slab"], V1))


@pytest.mark.req("EQR-007")
def test_a_superseded_equation_stays_retrievable():
    registry = build_equation_registry()
    registry.supersede(
        EQUATION_IDS["greybody"],
        V1,
        _equation(
            EQUATION_IDS["greybody"],
            V2,
            name="spectral emission replacing the grey approximation",
            mathematical_form="M = integral of epsilon(lambda) B_lambda(T) d lambda",
        ),
    )
    old = registry.entry(EQUATION_IDS["greybody"], V1)
    assert not old.is_current
    assert old.superseded_by == (EQUATION_IDS["greybody"], "2.0.0")
    # The published result that used v1 can still be reconstructed.
    assert (
        registry.get(EQUATION_IDS["greybody"], V1).mathematical_form == "M = epsilon * sigma * T^4"
    )


@pytest.mark.req("EQR-004")
def test_a_calibrated_parameter_cannot_be_declared_a_fundamental_law():
    with pytest.raises(ConfigurationError, match="EQR-004"):
        _equation(
            "EQ-DISGUISED",
            category=EquationCategory.FUNDAMENTAL_LAW,
            parameters=(("tuning_factor", "calibrated", "fit to the target record"),),
        )


@pytest.mark.req("EQR-003", "INF-003")
def test_fitted_equations_are_identifiable_for_the_independence_statement():
    registry = build_equation_registry()
    registry.register(
        _equation(
            "EQ-FITTED",
            category=EquationCategory.CALIBRATED_CORRECTION,
            parameters=(("offset", "calibrated", "fit to observations"),),
        )
    )
    fitted = {e.record_id for e in registry.fitted_equations()}
    assert "EQ-FITTED" in fitted
    assert EQUATION_IDS["slab"] not in fitted


@pytest.mark.req("EQR-002", "VAL-001")
def test_an_equation_without_a_correctness_test_is_refused():
    with pytest.raises(ConfigurationError, match="VAL-001"):
        _equation("EQ-UNTESTED", correctness_tests=())


@pytest.mark.req("MOD-002")
def test_validity_domain_reports_which_variable_left_its_range():
    domain = ValidityDomain(ranges=(("T", 100.0, 400.0, "K"),))
    ok, violations = domain.contains(T=300.0)
    assert ok and violations == ()
    ok, violations = domain.contains(T=5000.0)
    assert not ok and "T=5000" in violations[0]


# --- resource preflight ----------------------------------------------------


@pytest.mark.req("CMP-004", "KER-007")
def test_preflight_estimates_cost_per_module():
    kernel = airless_radiating_body(master_dt=3600.0, solar_step_multiple=24)
    estimate = estimate_cost(
        schedule=kernel.schedule, n_steps=240, budget=kernel.budget, spectral_bins=1
    )
    assert set(estimate.per_module) == {"solar", "orbit", "surface", "space"}
    # The hourly surface dominates the daily solar and orbital modules.
    assert estimate.per_module["surface"] > estimate.per_module["solar"] * 10
    assert estimate.within_budget


@pytest.mark.req("CMP-005")
def test_exceeding_the_declared_budget_fails_before_the_run():
    kernel = airless_radiating_body()
    tight = ResourceBudget(max_work_units=1.0, on_exceed="fail")
    estimate = estimate_cost(
        schedule=kernel.schedule, n_steps=10_000, budget=tight, spectral_bins=1
    )
    assert not estimate.within_budget
    with pytest.raises(ResourceBudgetError, match="exceeded the declared budget"):
        assert_within_budget(estimate, tight)


@pytest.mark.req("CMP-005")
def test_the_warn_policy_records_the_overrun_instead_of_stopping():
    kernel = airless_radiating_body()
    lenient = ResourceBudget(max_work_units=1.0, on_exceed="warn")
    estimate = estimate_cost(
        schedule=kernel.schedule, n_steps=10_000, budget=lenient, spectral_bins=1
    )
    assert_within_budget(estimate, lenient)
    assert any("exceeded the declared budget" in w for w in estimate.warnings)


@pytest.mark.req("CMP-001", "MOD-013")
def test_cost_scales_by_the_declared_scaling_variables_only():
    cost = CostDeclaration(
        module_id="m",
        time_complexity="O(n_bands)",
        memory_complexity="O(n_bands)",
        scaling_variables=("n_bands",),
        work_units_per_step=2.0,
        peak_memory_bytes=1000,
        benchmark_configuration="1 band",
    )
    scaled = cost.scaled({"n_bands": 8.0, "n_cells": 1000.0})
    assert scaled.work_units_per_step == pytest.approx(16.0)
    assert scaled.peak_memory_bytes == 8000  # n_cells is not declared, so it is ignored


@pytest.mark.req("CMP-008", "SYS-005")
def test_the_reference_kernel_fits_a_workstation_budget():
    kernel = airless_radiating_body()
    assert kernel.budget.resource_class is ResourceClass.ENGINEERING_WORKSTATION
    for schedule in kernel.schedule.ordered():
        cost = schedule.module.contract.cost
        assert cost is not None
        assert not cost.accelerator_required
        assert ResourceClass.ENGINEERING_WORKSTATION in cost.supported_resource_classes
