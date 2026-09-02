"""Stage 0: transfer reconciliation and conservation accounting."""

from __future__ import annotations

import pytest

from mfespf.core.errors import ConservationError, CouplingError
from mfespf.coupling.conservation import (
    ConservationAuditor,
    ModuleBudget,
    apply_ledger_to_budgets,
)
from mfespf.coupling.transfers import (
    ConservedQuantity,
    Transfer,
    TransferKind,
    TransferLedger,
)

pytestmark = pytest.mark.stage(0)

ENERGY = ConservedQuantity.ENERGY


def _transfer(sender: str, receiver: str, amount: float, **kwargs) -> Transfer:
    defaults = {
        "transfer_id": f"{sender}->{receiver}@{kwargs.get('start_time', 0.0)}",
        "quantity": ENERGY,
        "units": "J m^-2",
        "start_time": 0.0,
        "end_time": 100.0,
    }
    defaults.update(kwargs)
    return Transfer(sender_module=sender, receiver_module=receiver, amount=amount, **defaults)


@pytest.mark.req("CON-002", "CON-003")
def test_transfer_declares_both_parties_and_its_interval():
    t = _transfer("a", "b", 5.0)
    assert t.interval == 100.0
    assert t.mean_rate == pytest.approx(0.05)
    assert t.pairing_key == ("a", "b", "energy", "", "")


@pytest.mark.req("CON-002")
def test_a_module_cannot_transfer_to_itself():
    with pytest.raises(CouplingError, match="internal source or sink is not a transfer"):
        _transfer("a", "a", 1.0)


@pytest.mark.req("CON-004", "KER-005")
def test_matched_transfers_reconcile_to_zero():
    ledger = TransferLedger()
    ledger.record_export(_transfer("a", "b", 10.0))
    ledger.record_import(_transfer("a", "b", 10.0))
    report = ConservationAuditor(tolerance=1e-12).reconcile(ledger, [], (0.0, 100.0))
    assert not report.unmatched_transfers
    assert report.worst_coupling_mismatch() == pytest.approx(0.0)


@pytest.mark.req("CON-004", "CON-005")
def test_sender_receiver_mismatch_is_reported_as_a_coupling_residual():
    ledger = TransferLedger()
    ledger.record_export(_transfer("a", "b", 10.0))
    ledger.record_import(_transfer("a", "b", 9.0))
    report = ConservationAuditor(tolerance=1e-12).reconcile(ledger, [], (0.0, 100.0))
    assert report.coupling_residuals[0].mismatch == pytest.approx(-1.0)
    assert not report.passed
    with pytest.raises(ConservationError):
        ConservationAuditor(tolerance=1e-12).assert_valid(report)


@pytest.mark.req("CON-004")
def test_a_transfer_recorded_by_only_one_side_is_unmatched():
    ledger = TransferLedger()
    ledger.record_export(_transfer("a", "b", 10.0))
    report = ConservationAuditor().reconcile(ledger, [], (0.0, 100.0))
    assert any("no matching import" in m for m in report.unmatched_transfers)


@pytest.mark.req("CON-006")
def test_overlapping_records_of_one_transaction_are_a_double_count():
    ledger = TransferLedger()
    ledger.record_import(_transfer("a", "b", 5.0, start_time=0.0, end_time=100.0))
    ledger.record_import(
        Transfer(
            transfer_id="second",
            sender_module="a",
            receiver_module="b",
            quantity=ENERGY,
            amount=5.0,
            units="J m^-2",
            start_time=50.0,
            end_time=150.0,
        )
    )
    assert ledger.duplicates()


@pytest.mark.req("CON-006")
def test_tiled_substep_records_are_not_a_double_count():
    """A fast module reporting one exchange per substep is not double counting."""
    ledger = TransferLedger()
    for i in range(4):
        ledger.record_import(_transfer("a", "b", 2.5, start_time=i * 25.0, end_time=(i + 1) * 25.0))
    assert ledger.duplicates() == []


@pytest.mark.req("CON-008")
def test_numerical_sources_are_accounted_separately_from_physical_ones():
    budget = ModuleBudget("m", ENERGY, "J m^-2", stored_initial=0.0, stored_final=15.0)
    budget.add_source("phase_change", 10.0)
    budget.add_source("flux_limiter", 5.0, numerical=True)
    assert budget.net_physical == pytest.approx(10.0)
    assert budget.net_numerical == pytest.approx(5.0)
    assert budget.residual == pytest.approx(0.0)


@pytest.mark.req("CON-008")
def test_a_numerical_transfer_must_name_its_operator():
    with pytest.raises(CouplingError, match="CON-008"):
        _transfer("a", "b", 1.0, kind=TransferKind.NUMERICAL)


@pytest.mark.req("CON-009")
def test_budget_publishes_every_line_of_the_accounting():
    budget = ModuleBudget("m", ENERGY, "J m^-2", stored_initial=100.0, stored_final=120.0)
    ledger = TransferLedger()
    ledger.record_import(_transfer("src", "m", 30.0))
    ledger.record_export(_transfer("m", "sink", 10.0))
    apply_ledger_to_budgets(ledger, {("m", ENERGY): budget})
    assert budget.imports == pytest.approx(30.0)
    assert budget.exports == pytest.approx(10.0)
    assert budget.residual == pytest.approx(0.0)
    assert "residual" in budget.describe()


@pytest.mark.req("CON-010")
def test_a_residual_as_large_as_the_claimed_effect_invalidates_it():
    report = ConservationAuditor().reconcile(TransferLedger(), [], (0.0, 1.0))
    leaky = ModuleBudget("m", ENERGY, "J m^-2", stored_initial=0.0, stored_final=1.0)
    report.budgets = [leaky]
    # Residual of 1.0 against a claimed effect of 0.5: the run cannot decide it.
    assert not report.can_support_effect(ENERGY, 0.5)
    # Against a much larger claimed effect it can.
    assert report.can_support_effect(ENERGY, 100.0)


@pytest.mark.req("CON-009")
def test_residual_is_scaled_against_the_largest_budget_term_not_the_net():
    """A budget whose flows cancel must still be judged against those flows."""
    budget = ModuleBudget("m", ENERGY, "J m^-2", stored_initial=0.0, stored_final=0.0)
    budget.imports = 1.0e6
    budget.exports = 1.0e6 - 1.0
    assert budget.stored_change == 0.0
    assert budget.relative_residual == pytest.approx(1.0e-6)
