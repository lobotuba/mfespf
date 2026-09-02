"""Conserved-quantity transfers between modules.

Charter section 8 and requirements CON-001 through CON-008.

A transfer is a transaction, not a number: it names both parties, the quantity,
the interval and the spatial mapping, and it is recorded once by the sender and
once by the receiver so the two records can be reconciled (CON-004).  The sign
convention is fixed at positive into the receiver (CON-003), which is what makes
double counting (CON-006) detectable rather than a matter of reading two modules'
prose carefully.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from dataclasses import field as dc_field
from enum import StrEnum

from mfespf.core.errors import CouplingError
from mfespf.core.units import parse_unit

__all__ = [
    "ConservedQuantity",
    "Transfer",
    "TransferKind",
    "TransferLedger",
]


class ConservedQuantity(StrEnum):
    """CON-001 / charter 8.1."""

    ENERGY = "energy"
    MASS = "mass"
    WATER = "water"
    MOMENTUM = "momentum"
    ANGULAR_MOMENTUM = "angular_momentum"
    CHARGE = "charge"
    ELEMENT = "element"
    SPECIES = "species"

    @property
    def si_units(self) -> str:
        """Coherent SI units for a *stock* of this quantity, per unit area.

        The framework's zero-dimensional and column configurations account per
        unit area, so energy is J m^-2 rather than J.  A three-dimensional
        configuration declares its own areal or volumetric basis on the interface;
        the units travel with the transfer either way.
        """
        return {
            ConservedQuantity.ENERGY: "J m^-2",
            ConservedQuantity.MASS: "kg m^-2",
            ConservedQuantity.WATER: "kg m^-2",
            ConservedQuantity.MOMENTUM: "kg m s^-1",
            ConservedQuantity.ANGULAR_MOMENTUM: "kg m s^-1",
            ConservedQuantity.CHARGE: "C",
            ConservedQuantity.ELEMENT: "mol",
            ConservedQuantity.SPECIES: "mol",
        }[self]


class TransferKind(StrEnum):
    PHYSICAL = "physical"
    NUMERICAL = "numerical"  # CON-008: stabilisation, clipping, filtering


@dataclass(frozen=True)
class Transfer:
    """One conserved-quantity transaction (charter 8.3 / CON-002).

    ``amount`` is the quantity transferred over ``[start_time, end_time]``,
    positive into ``receiver_module``.  A flux is expressed by giving the amount
    over its interval, so no consumer has to guess whether a number was a rate or
    a total -- the ambiguity that CON-006 double counting usually hides in.
    """

    transfer_id: str
    sender_module: str
    receiver_module: str
    quantity: ConservedQuantity
    amount: float
    units: str
    start_time: float
    end_time: float
    spatial_mapping: str = "global_scalar"
    kind: TransferKind = TransferKind.PHYSICAL
    numerical_transfer_error: float = 0.0
    source_or_sink_reference: str = ""
    equation_id: str | None = None
    species: str | None = None
    channel: str = ""
    """Distinguishes concurrent transfers of the same quantity between the same
    pair, such as reflected shortwave and emitted longwave leaving a surface for
    space. Without it the two would look like one transaction recorded twice."""

    def __post_init__(self) -> None:
        parse_unit(self.units)
        if self.end_time < self.start_time:
            raise CouplingError(
                f"transfer {self.transfer_id!r} ends before it starts "
                f"({self.end_time} < {self.start_time})"
            )
        if self.sender_module == self.receiver_module:
            raise CouplingError(
                f"transfer {self.transfer_id!r} names {self.sender_module!r} as both sender and "
                f"receiver; an internal source or sink is not a transfer (charter 4.24)"
            )
        if self.kind is TransferKind.NUMERICAL and not self.source_or_sink_reference:
            raise CouplingError(
                f"numerical transfer {self.transfer_id!r} must name the operator responsible "
                f"(CON-008)"
            )

    @property
    def interval(self) -> float:
        return self.end_time - self.start_time

    @property
    def mean_rate(self) -> float:
        return self.amount / self.interval if self.interval else 0.0

    @property
    def pairing_key(self) -> tuple[str, str, str, str, str]:
        """What must match between a sender's and a receiver's record.

        Deliberately excludes the time bounds.  Under multirate coupling the two
        sides of one exchange are recorded over different intervals -- a source
        stepping daily writes one record while the receiver stepping hourly
        writes twenty-four -- and charter 8.4 reconciles at the end of the
        coupling interval, not per step.  The interval is checked separately by
        ``within``; what must pair is the transaction's identity.
        """
        return (
            self.sender_module,
            self.receiver_module,
            self.quantity.value,
            self.species or "",
            self.channel,
        )

    def within(self, interval: tuple[float, float], *, tolerance: float = 1e-9) -> bool:
        """True when this transfer lies inside a reconciliation interval."""
        lo, hi = interval
        span = max(abs(hi - lo), 1.0)
        return self.start_time >= lo - tolerance * span and self.end_time <= hi + tolerance * span

    def describe(self) -> str:
        return (
            f"{self.transfer_id}: {self.sender_module} -> {self.receiver_module} "
            f"{self.amount:+.6g} {self.units} of {self.quantity.value} "
            f"over [{self.start_time:g}, {self.end_time:g}]"
        )


@dataclass
class TransferLedger:
    """Paired sender and receiver records for one coupling interval.

    Both sides record the same physical transaction.  Keeping the two records
    separate and reconciling them (rather than having one side write for both) is
    what turns a sign-convention mistake or a missing reciprocal exchange into a
    reported residual instead of a silent energy leak.
    """

    exports: list[Transfer] = dc_field(default_factory=list)
    imports: list[Transfer] = dc_field(default_factory=list)

    def record_export(self, transfer: Transfer) -> None:
        self.exports.append(transfer)

    def record_import(self, transfer: Transfer) -> None:
        self.imports.append(transfer)

    def all_transfers(self) -> list[Transfer]:
        return [*self.exports, *self.imports]

    def clear(self) -> None:
        self.exports.clear()
        self.imports.clear()

    def net_for(self, module_id: str, quantity: ConservedQuantity) -> float:
        """Net amount of ``quantity`` received by ``module_id`` this interval.

        Counted from the receiver-side records, with the sender's own budget
        debited by the same magnitude.  Reconciliation is what checks the two
        agree; this method assumes nothing about that.
        """
        total = 0.0
        for t in self.imports:
            if t.quantity is not quantity:
                continue
            if t.receiver_module == module_id:
                total += t.amount
            elif t.sender_module == module_id:
                total -= t.amount
        return total

    def duplicates(self) -> list[tuple[str, list[str]]]:
        """Transfers claiming the same transaction twice (CON-006).

        Detection is by overlapping time intervals rather than by count.  Under
        multirate coupling a fast module legitimately writes many records for one
        exchange, tiling the coupling interval end to end; those are one
        transaction reported at finer granularity.  Two records of the same
        transaction whose intervals *overlap* are the real double count -- the
        same energy entering the coupled budget more than once.
        """
        found: list[tuple[str, list[str]]] = []
        for side in (self.exports, self.imports):
            grouped: dict[tuple, list[Transfer]] = {}
            for t in side:
                grouped.setdefault(t.pairing_key, []).append(t)
            for key, entries in grouped.items():
                ordered = sorted(entries, key=lambda t: (t.start_time, t.end_time))
                for earlier, later in itertools.pairwise(ordered):
                    span = max(earlier.interval, later.interval, 1.0)
                    if later.start_time < earlier.end_time - 1e-9 * span:
                        channel = f" [{key[4]}]" if key[4] else ""
                        found.append(
                            (
                                f"{key[0]} -> {key[1]} {key[2]}{channel}",
                                [earlier.transfer_id, later.transfer_id],
                            )
                        )
        return found
