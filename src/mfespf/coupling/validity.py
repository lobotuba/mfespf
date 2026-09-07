"""Checking, at run time, that equations are evaluated inside their declared domains.

Charter 4.26 and section 19; requirements EQR-002, EQR-003, SYS-006.

Every equation record in this framework carries a ``ValidityDomain`` saying what
conditions the relation was established for. Until this module existed,
``ValidityDomain.contains`` was called from the registry's own tests and from no
running code, so the framework declared the conditions each relation is good for
and then never checked whether a run was inside them.

That was not a hypothetical gap. Evaluating a relation outside its calibration
conditions and getting a plausible-looking answer had happened five times before
this was written -- the carbon dioxide band's background opacity, which stood
behind four registered claims; the vapour band, fitted on a column carrying
23.2 kg m^-2 and evaluated on one carrying 8.0; the evaporating condensate
fraction, asserted over a range that excluded the measured value; the bulk plume
applied to a boundary layer it does not describe; and the same closure asked
about a global mean when it describes a convecting column. Every one was caught
by hand, late, and by accident, and the second was already outside a domain this
project had written down.

How a variable is resolved
--------------------------

An equation declares ranges against physical symbols -- ``T``, ``p``, ``W``.
Fields published during a run carry a ``symbol`` too, so the two vocabularies
already meet and no separate lookup table has to be invented and then kept in
step.

The resolution is deliberately **scoped to the citing module**: only fields that
module declares as an input or an output are considered. A global match by
symbol would be actively dangerous, because symbols collide across subject
areas. ``alpha`` is the evaporating condensate fraction in the plume equations
and the surface albedo in the radiation fields; matching globally would check
one against the other's range and report confidently on nothing.

What "unchecked" means
----------------------

Three outcomes, not two. A variable can be inside its range, outside it, or
**not resolvable** -- because no field carries that symbol, or because the units
do not match the declaration. Unchecked is reported rather than folded into
success. A domain check that quietly passed when it could not find the variable
would be worse than no check at all, since it would licence exactly the
confidence it fails to earn.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np

__all__ = [
    "DomainCheck",
    "DomainReport",
    "check_validity_domains",
]


@dataclass(frozen=True, slots=True)
class DomainCheck:
    """One declared range, checked against one run's values."""

    module_id: str
    equation_id: str
    variable: str
    status: str  # "inside" | "violated" | "unchecked"
    detail: str

    @property
    def violated(self) -> bool:
        return self.status == "violated"

    @property
    def unchecked(self) -> bool:
        return self.status == "unchecked"


@dataclass(frozen=True, slots=True)
class DomainReport:
    """Every declared range a run could be judged against, and the verdict."""

    checks: tuple[DomainCheck, ...]

    @property
    def violations(self) -> tuple[DomainCheck, ...]:
        return tuple(c for c in self.checks if c.violated)

    @property
    def unchecked(self) -> tuple[DomainCheck, ...]:
        return tuple(c for c in self.checks if c.unchecked)

    @property
    def passed(self) -> bool:
        """True when nothing checkable was violated.

        Deliberately *not* affected by unchecked variables: a run is not invalid
        because a domain could not be resolved. The count is reported so the
        coverage is visible, because a report of no violations over two checks
        means something different from no violations over forty.

        A violation *does* invalidate: the driver records this as a mandatory
        gate. Coverage is presently 9 of 30 declared ranges on a column
        configuration, so a pass is a weaker statement than a failure.
        """
        return not self.violations

    def summary(self) -> str:
        inside = sum(1 for c in self.checks if c.status == "inside")
        return (
            f"validity domains: {inside} inside, {len(self.violations)} violated, "
            f"{len(self.unchecked)} unresolvable"
        )

    def describe_violations(self) -> tuple[str, ...]:
        return tuple(f"{c.module_id} cites {c.equation_id}: {c.detail}" for c in self.violations)


def _extent(values: Any) -> tuple[float, float] | None:
    """The range a field actually spans, or None if it is not numeric."""
    array = np.asarray(values, dtype=float)
    if array.size == 0 or not np.all(np.isfinite(array)):
        return None
    return (float(np.min(array)), float(np.max(array)))


def check_validity_domains(
    *,
    registry: Any,
    modules: Iterable[Any],
    fields: Mapping[str, Any],
) -> DomainReport:
    """Judge a run against the declared domains of the equations it cites.

    ``registry`` supplies equation records by id, ``modules`` are the modules
    that actually ran, and ``fields`` are the quantities they published.
    """
    checks: list[DomainCheck] = []

    for module in modules:
        contract = getattr(module, "contract", None)
        if contract is None:
            continue

        # Scope: this module's own declared inputs and outputs, nothing else.
        visible = set(getattr(contract, "produced_outputs", ()) or ())
        visible |= set(getattr(contract, "required_inputs", ()) or ())
        in_scope: list[tuple[str, Any, str]] = []
        for field_id in sorted(visible):
            field = fields.get(field_id)
            symbol = getattr(field, "symbol", None) if field is not None else None
            if symbol:
                in_scope.append((field_id, field, str(symbol)))

        for equation_id in getattr(contract, "governing_equation_ids", ()) or ():
            record = _lookup(registry, equation_id)
            if record is None:
                continue
            for name, low, high, units in record.validity_domain.ranges:
                resolved, ambiguity = _resolve(name, units, in_scope)
                checks.append(
                    _judge(
                        module_id=contract.module_id,
                        equation_id=equation_id,
                        name=name,
                        low=low,
                        high=high,
                        units=units,
                        resolved=resolved,
                        ambiguity=ambiguity,
                    )
                )

    return DomainReport(tuple(checks))


def _root(symbol: str) -> str:
    """The symbol without its subscript: ``T_a`` and ``T_i`` are both ``T``."""
    return symbol.split("_", 1)[0]


def _resolve(
    name: str, units: str, in_scope: list[tuple[str, Any, str]]
) -> tuple[tuple[str, Any] | None, str]:
    """Find the field a declared range refers to, or say why not.

    Exact symbol first. Failing that, a field whose symbol has the same root and
    the same units -- an equation writes ``T`` where a column publishes ``T_a``,
    and within one module's own declared fields that is the same quantity.

    Ambiguity is refused rather than broken by a rule. If a module can see both
    ``T_a`` and ``T_s`` in kelvin, nothing here knows which one the equation's
    ``T`` meant, and guessing would produce a confident check of the wrong
    quantity -- the exact failure this module was written to stop.
    """
    for field_id, field, symbol in in_scope:
        if symbol == name:
            return (field_id, field), ""

    root = _root(name)
    candidates = [
        (field_id, field, symbol)
        for field_id, field, symbol in in_scope
        if _root(symbol) == root and str(getattr(field, "units", "")) == units
    ]
    if len(candidates) == 1:
        field_id, field, _symbol = candidates[0]
        return (field_id, field), ""
    if len(candidates) > 1:
        listed = ", ".join(sorted(f"{s} ({f})" for f, _v, s in candidates))
        return None, f"{name!r} could be any of {listed}; the equation does not say which"
    return None, ""


def _judge(
    *,
    module_id: str,
    equation_id: str,
    name: str,
    low: float,
    high: float,
    units: str,
    resolved: tuple[str, Any] | None,
    ambiguity: str = "",
) -> DomainCheck:
    if resolved is None:
        return DomainCheck(
            module_id,
            equation_id,
            name,
            "unchecked",
            ambiguity or f"no field this module declares carries the symbol {name!r}",
        )

    field_id, field = resolved
    field_units = str(getattr(field, "units", ""))
    if field_units != units:
        # Not a silent conversion. Comparing a number against a range stated in
        # different units is the failure this module exists to prevent.
        return DomainCheck(
            module_id,
            equation_id,
            name,
            "unchecked",
            f"{field_id} is in {field_units!r}, the domain is stated in {units!r}",
        )

    extent = _extent(getattr(field, "values", None))
    if extent is None:
        return DomainCheck(
            module_id, equation_id, name, "unchecked", f"{field_id} holds no finite values"
        )

    lowest, highest = extent
    if lowest < low or highest > high:
        return DomainCheck(
            module_id,
            equation_id,
            name,
            "violated",
            f"{name} spans [{lowest:g}, {highest:g}] {units} in {field_id}, "
            f"outside the declared [{low:g}, {high:g}]",
        )
    return DomainCheck(
        module_id,
        equation_id,
        name,
        "inside",
        f"{name} spans [{lowest:g}, {highest:g}] {units} within [{low:g}, {high:g}]",
    )


def _lookup(registry: Any, equation_id: str) -> Any:
    """Fetch a record without caring which accessor the registry offers."""
    for accessor in ("current_record", "record", "get"):
        method = getattr(registry, accessor, None)
        if callable(method):
            try:
                found = method(equation_id)
            except Exception:
                continue
            if found is not None:
                return found
    records = getattr(registry, "current", None)
    if callable(records):
        for record in records():
            if getattr(record, "record_id", None) == equation_id:
                return record
    return None
