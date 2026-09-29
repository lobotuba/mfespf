"""Conservation reconciliation and budget reporting.

Charter sections 4.25 and 8.4-8.6, requirements CON-004 through CON-010, and
KER-005.

For a conserved quantity Q over an accounting interval, the residual is

    R_Q = dQ_stored - (imports - exports + physical_sources - physical_sinks
                       + numerical_sources - numerical_sinks)

Numerical sources and sinks are tracked separately from physical ones
throughout (CON-008).  That separation is the point: a budget that closes only
because a limiter injected energy is not a closed budget, and CON-010 needs to
compare the unexplained residual against the size of the effect being claimed.
"""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field as dc_field

from mfespf.core.errors import ConservationError
from mfespf.coupling.transfers import ConservedQuantity, Transfer, TransferKind, TransferLedger

__all__ = [
    "ConservationAuditor",
    "ConservationReport",
    "CouplingResidual",
    "ModuleBudget",
]


@dataclass
class ModuleBudget:
    """One module's accounting for one quantity over one interval (CON-009)."""

    module_id: str
    quantity: ConservedQuantity
    units: str
    stored_initial: float = 0.0
    stored_final: float = 0.0
    imports: float = 0.0
    exports: float = 0.0
    physical_sources: float = 0.0
    physical_sinks: float = 0.0
    numerical_sources: float = 0.0
    numerical_sinks: float = 0.0
    source_labels: dict[str, float] = dc_field(default_factory=dict)

    @property
    def stored_change(self) -> float:
        return self.stored_final - self.stored_initial

    @property
    def net_transfer(self) -> float:
        return self.imports - self.exports

    @property
    def net_physical(self) -> float:
        return self.physical_sources - self.physical_sinks

    @property
    def net_numerical(self) -> float:
        return self.numerical_sources - self.numerical_sinks

    @property
    def residual(self) -> float:
        """Charter 4.25.  Zero to round-off for a correct explicit budget."""
        return self.stored_change - (self.net_transfer + self.net_physical + self.net_numerical)

    @property
    def scale(self) -> float:
        """Magnitude the residual is judged against.

        The largest term in the budget, not the stored change: a budget whose
        stored change is near zero because large imports and exports cancel is
        still wrong if it is wrong by a percent of those flows.
        """
        return max(
            abs(self.stored_change),
            abs(self.imports),
            abs(self.exports),
            abs(self.net_physical),
            abs(self.stored_initial),
            1e-300,
        )

    @property
    def relative_residual(self) -> float:
        return abs(self.residual) / self.scale

    def add_source(self, label: str, amount: float, *, numerical: bool = False) -> None:
        """Add a labeled internal source (positive) or sink (negative)."""
        self.source_labels[label] = self.source_labels.get(label, 0.0) + amount
        if numerical:
            if amount >= 0:
                self.numerical_sources += amount
            else:
                self.numerical_sinks += -amount
        else:
            if amount >= 0:
                self.physical_sources += amount
            else:
                self.physical_sinks += -amount

    def merge_later(self, later: ModuleBudget) -> None:
        """Absorb a subsequent sub-interval's budget for the same module.

        A module stepping faster than the coupling interval reports once per
        substep.  Merging keeps the stored value from the start of the interval
        and the one from its end, and sums the labeled sources between, so the
        reconciled budget describes the whole interval rather than its last slice.
        """
        if (later.module_id, later.quantity) != (self.module_id, self.quantity):
            raise ConservationError(
                f"cannot merge budgets for different module/quantity pairs: "
                f"{(self.module_id, self.quantity.value)} and "
                f"{(later.module_id, later.quantity.value)}"
            )
        self.stored_final = later.stored_final
        self.imports += later.imports
        self.exports += later.exports
        self.physical_sources += later.physical_sources
        self.physical_sinks += later.physical_sinks
        self.numerical_sources += later.numerical_sources
        self.numerical_sinks += later.numerical_sinks
        for label, amount in later.source_labels.items():
            self.source_labels[label] = self.source_labels.get(label, 0.0) + amount

    def describe(self) -> str:
        return (
            f"{self.module_id} [{self.quantity.value}, {self.units}]\n"
            f"  stored      {self.stored_initial:+.6g} -> {self.stored_final:+.6g} "
            f"(d = {self.stored_change:+.6g})\n"
            f"  imports     {self.imports:+.6g}\n"
            f"  exports     {self.exports:+.6g}\n"
            f"  phys src    {self.physical_sources:+.6g}\n"
            f"  phys sink   {self.physical_sinks:+.6g}\n"
            f"  num  src    {self.numerical_sources:+.6g}\n"
            f"  num  sink   {self.numerical_sinks:+.6g}\n"
            f"  residual    {self.residual:+.6g} (relative {self.relative_residual:.3e})"
        )


@dataclass(frozen=True, slots=True)
class CouplingResidual:
    """A sender/receiver mismatch for one transaction (CON-004, CON-005)."""

    sender: str
    receiver: str
    quantity: ConservedQuantity
    exported: float
    imported: float
    units: str
    start_time: float
    end_time: float

    @property
    def mismatch(self) -> float:
        return self.imported - self.exported

    @property
    def relative_mismatch(self) -> float:
        scale = max(abs(self.exported), abs(self.imported), 1e-300)
        return abs(self.mismatch) / scale

    def describe(self) -> str:
        return (
            f"{self.sender} -> {self.receiver} [{self.quantity.value}]: exported "
            f"{self.exported:+.6g} {self.units}, imported {self.imported:+.6g} "
            f"(mismatch {self.mismatch:+.3g}, relative {self.relative_mismatch:.3e})"
        )


@dataclass
class ConservationReport:
    """The per-run conservation publication required by CON-009."""

    interval: tuple[float, float]
    budgets: list[ModuleBudget] = dc_field(default_factory=list)
    coupling_residuals: list[CouplingResidual] = dc_field(default_factory=list)
    unmatched_transfers: list[str] = dc_field(default_factory=list)
    double_counts: list[str] = dc_field(default_factory=list)
    tolerance: float = 1e-10

    def budgets_for(self, quantity: ConservedQuantity) -> list[ModuleBudget]:
        return [b for b in self.budgets if b.quantity is quantity]

    def system_residual(self, quantity: ConservedQuantity) -> float:
        """Summed residual across modules for one quantity.

        Internal transfers cancel in the sum, so this is the residual of the
        modeled system against its boundary -- the number a Level 1 correctness
        claim (charter 2.2) rests on.
        """
        return sum(b.residual for b in self.budgets_for(quantity))

    def worst_relative_residual(self) -> float:
        return max((b.relative_residual for b in self.budgets), default=0.0)

    def worst_coupling_mismatch(self) -> float:
        return max((r.relative_mismatch for r in self.coupling_residuals), default=0.0)

    @property
    def passed(self) -> bool:
        return (
            not self.unmatched_transfers
            and not self.double_counts
            and self.worst_relative_residual() <= self.tolerance
            and self.worst_coupling_mismatch() <= self.tolerance
        )

    def failures(self) -> list[str]:
        out: list[str] = []
        for t in self.unmatched_transfers:
            out.append(f"unmatched transfer: {t}")
        for d in self.double_counts:
            out.append(f"double-counted transfer: {d}")
        for b in self.budgets:
            if b.relative_residual > self.tolerance:
                out.append(
                    f"{b.module_id} {b.quantity.value}: relative residual "
                    f"{b.relative_residual:.3e} exceeds tolerance {self.tolerance:.3e}"
                )
        for r in self.coupling_residuals:
            if r.relative_mismatch > self.tolerance:
                out.append(f"coupling residual: {r.describe()}")
        return out

    def can_support_effect(self, quantity: ConservedQuantity, effect_magnitude: float) -> bool:
        """CON-010: is the unexplained residual small compared with the claimed effect?

        Returns False when the residual is large enough to account for the effect,
        in which case the run must not be used to support or refute the claim.
        """
        residual = abs(self.system_residual(quantity))
        if effect_magnitude <= 0:
            return False
        return residual < 0.1 * abs(effect_magnitude)

    def describe(self) -> str:
        lines = [f"conservation report for interval [{self.interval[0]:g}, {self.interval[1]:g}]"]
        for b in self.budgets:
            lines.append(b.describe())
        for r in self.coupling_residuals:
            lines.append("  " + r.describe())
        if self.passed:
            lines.append(f"PASS (tolerance {self.tolerance:.3e})")
        else:
            lines.append("FAIL:")
            lines.extend(f"  - {f}" for f in self.failures())
        return "\n".join(lines)


class ConservationAuditor:
    """Reconciles transfers and assembles the conservation report (KER-005)."""

    def __init__(self, tolerance: float = 1e-10) -> None:
        self.tolerance = tolerance

    def reconcile(
        self,
        ledger: TransferLedger,
        budgets: list[ModuleBudget],
        interval: tuple[float, float],
    ) -> ConservationReport:
        report = ConservationReport(interval=interval, tolerance=self.tolerance)
        report.budgets = list(budgets)

        exports_by_key: dict[tuple, list[Transfer]] = {}
        for t in ledger.exports:
            exports_by_key.setdefault(t.pairing_key, []).append(t)
        imports_by_key: dict[tuple, list[Transfer]] = {}
        for t in ledger.imports:
            imports_by_key.setdefault(t.pairing_key, []).append(t)

        for key in sorted(set(exports_by_key) | set(imports_by_key), key=str):
            exported = sum(t.amount for t in exports_by_key.get(key, []))
            imported = sum(t.amount for t in imports_by_key.get(key, []))
            sample = (exports_by_key.get(key) or imports_by_key[key])[0]

            if key not in exports_by_key or key not in imports_by_key:
                side = "import" if key not in imports_by_key else "export"
                report.unmatched_transfers.append(
                    f"{sample.sender_module} -> {sample.receiver_module} "
                    f"[{sample.quantity.value}] has no matching {side} record"
                )

            report.coupling_residuals.append(
                CouplingResidual(
                    sender=sample.sender_module,
                    receiver=sample.receiver_module,
                    quantity=sample.quantity,
                    exported=exported,
                    imported=imported,
                    units=sample.units,
                    start_time=sample.start_time,
                    end_time=sample.end_time,
                )
            )

        for label, ids in ledger.duplicates():
            report.double_counts.append(f"{label}: {ids}")

        return report

    def assert_valid(self, report: ConservationReport) -> None:
        """Raise when a report fails its declared tolerances (ACC-004)."""
        if not report.passed:
            raise ConservationError(
                "conservation reconciliation failed:\n  - " + "\n  - ".join(report.failures())
            )


def apply_ledger_to_budgets(
    ledger: TransferLedger, budgets: dict[tuple[str, ConservedQuantity], ModuleBudget]
) -> None:
    """Post every transfer into the sender's and receiver's budgets.

    Numerical transfers land in the numerical source/sink lines so that CON-008
    stays satisfied without the caller having to remember which kind it recorded.
    """
    # Each budget is per unit of its own module's area; a transfer between
    # regions is recorded per unit of global area and converted here, once.
    for t in ledger.imports:
        receiver = budgets.get((t.receiver_module, t.quantity))
        if receiver is not None:
            if t.kind is TransferKind.NUMERICAL:
                receiver.add_source(
                    t.source_or_sink_reference or t.transfer_id,
                    t.receiver_amount,
                    numerical=True,
                )
            else:
                receiver.imports += t.receiver_amount
    for t in ledger.exports:
        sender = budgets.get((t.sender_module, t.quantity))
        if sender is not None:
            if t.kind is TransferKind.NUMERICAL:
                sender.add_source(
                    t.source_or_sink_reference or t.transfer_id,
                    -t.sender_amount,
                    numerical=True,
                )
            else:
                sender.exports += t.sender_amount
