"""Stage 0: state ownership and the labeled-tendency discipline."""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import OwnershipError
from mfespf.core.versioning import Version
from mfespf.registries.states import (
    StateBudget,
    StateEntry,
    StateKind,
    StateRegistry,
    TendencyContribution,
    UpdateAuthority,
)

pytestmark = pytest.mark.stage(0)

V1 = Version(1, 0, 0)


def _entry(state_id: str, owner: str, **kwargs) -> StateEntry:
    defaults = {
        "physical_name": state_id,
        "units": "K",
        "state_kind": StateKind.PROGNOSTIC,
        "update_authority": UpdateAuthority.OWNER_ONLY,
    }
    defaults.update(kwargs)
    return StateEntry(record_id=state_id, version=V1, owner_module_id=owner, **defaults)


@pytest.mark.req("STA-001", "MOD-010")
def test_a_module_cannot_write_state_it_does_not_own():
    registry = StateRegistry()
    registry.register(_entry("air_temperature", "atmosphere"))
    registry.assert_may_write("atmosphere", "air_temperature")
    with pytest.raises(OwnershipError, match="owned by 'atmosphere'"):
        registry.assert_may_write("ocean", "air_temperature")


@pytest.mark.req("STA-004", "STA-001")
def test_duplicate_ownership_is_rejected_before_initialization():
    registry = StateRegistry()
    registry.register(_entry("air_temperature", "atmosphere"))
    with pytest.raises(OwnershipError, match="exactly one owner"):
        registry.validate_declared_ownership(
            {
                "atmosphere": frozenset({"air_temperature"}),
                "radiation": frozenset({"air_temperature"}),
            }
        )


@pytest.mark.req("STA-005")
def test_an_unclaimed_prognostic_state_is_rejected():
    registry = StateRegistry()
    registry.register(_entry("air_temperature", "atmosphere"))
    with pytest.raises(OwnershipError, match="no claiming module"):
        registry.validate_declared_ownership({"atmosphere": frozenset()})


@pytest.mark.req("STA-003", "MOD-011")
def test_unauthorized_tendency_provider_is_refused():
    registry = StateRegistry()
    registry.register(
        _entry(
            "air_temperature",
            "atmosphere",
            authorized_tendency_providers=frozenset({"radiation"}),
        )
    )
    registry.assert_may_contribute(
        TendencyContribution(
            state_id="air_temperature",
            provider_module_id="radiation",
            label="radiative_heating",
            value=1.0e-5,
            units="K s^-1",
        )
    )
    with pytest.raises(OwnershipError, match="not an authorized tendency provider"):
        registry.assert_may_contribute(
            TendencyContribution(
                state_id="air_temperature",
                provider_module_id="biosphere",
                label="mystery_heating",
                value=1.0,
                units="K s^-1",
            )
        )


@pytest.mark.req("STA-002")
def test_tendency_units_must_match_the_state_or_its_rate():
    registry = StateRegistry()
    registry.register(
        _entry("air_temperature", "atmosphere", authorized_tendency_providers=frozenset({"rad"}))
    )
    with pytest.raises(OwnershipError, match="expected"):
        registry.assert_may_contribute(
            TendencyContribution(
                state_id="air_temperature",
                provider_module_id="rad",
                label="wrong_units",
                value=1.0,
                units="W m^-2",
            )
        )


@pytest.mark.req("STA-006")
def test_state_budget_reconstructs_the_update_from_labeled_contributions():
    contributions = [
        TendencyContribution("T", "atm", "radiative_heating", 2.0, "K s^-1"),
        TendencyContribution("T", "atm", "convective_cooling", -0.5, "K s^-1"),
        TendencyContribution("T", "atm", "limiter", 0.1, "K s^-1", is_numerical=True),
    ]
    dt = 10.0
    initial = np.array(300.0)
    final = initial + sum(c.value for c in contributions) * dt
    budget = StateBudget("T", initial, np.asarray(final), dt, contributions)

    assert budget.stored_change == pytest.approx(16.0)
    assert float(budget.residual()) == pytest.approx(0.0, abs=1e-12)
    assert budget.by_label()["radiative_heating"] == pytest.approx(20.0)
    # CON-008: the stabilisation term is separable from the physical ones.
    assert float(budget.numerical_part()) == pytest.approx(1.0)


@pytest.mark.req("STA-006")
def test_state_budget_residual_exposes_an_unlabeled_update():
    contributions = [TendencyContribution("T", "atm", "radiative_heating", 1.0, "K s^-1")]
    budget = StateBudget("T", np.array(300.0), np.array(315.0), 10.0, contributions)
    # 5 K of the change is not accounted for by any label.
    assert float(budget.residual()) == pytest.approx(5.0)


@pytest.mark.req("STA-007")
def test_registry_reports_the_owner_for_a_derived_view():
    registry = StateRegistry()
    registry.register(_entry("air_temperature", "atmosphere"))
    assert registry.owner_of("air_temperature") == "atmosphere"
