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

Over the run, not at the end of it
----------------------------------

The question a domain answers is whether a relation was *ever* evaluated
outside the conditions it was established for, and that is a property of the
whole trajectory. Judging the final state answers a different and much weaker
question, and it answers it wrongly in both directions: a closure that ran out
of range during spin-up and recovered is reported clean, and a closure that is
simply idle on the last step has no footprint to be judged on at all. The
convective adjustment is the case that forced this. It mixes on 79 percent of
steps and on none of the first four hundred, so the identical configuration was
resolvable at 2000 steps and unresolvable at 400 -- a fact about where the run
stopped rather than about the physics.

So ``EvaluatedExtents`` accumulates, per published field, the span it covered
while the run proceeded, and the judgement is made against that. **Resolution
still happens once, against the final state**, because which field a relation's
``T`` refers to is a property of the configuration's declared contracts and
does not vary from step to step; only the values do.

The immediate consequence is that a module's NaN mask can now say what it
should say -- *this step's* extent -- rather than the union over the run that
judging the final state forced it into. A layer mixed at 250 K during spin-up
now contributes 250 K, at the step it was mixed, instead of contributing
whatever it happens to hold when the run stops.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np

__all__ = [
    "INTEGRATED_OVER",
    "DomainCheck",
    "DomainReport",
    "EvaluatedExtents",
    "check_validity_domains",
]

#: Bind a symbol to this when the module integrates over the variable rather
#: than evaluating the relation at any value of it. The surface's thermal
#: emission is the Planck function integrated over every wavelength, so its
#: ``lambda`` has no value a range could be checked against -- and that is a
#: different thing from a variable the checker merely failed to find. Both stay
#: unchecked, because neither was checked; the distinction is that one is a gap
#: to close and the other is a category error that closing would introduce.
INTEGRATED_OVER = "<integrated over its full support>"


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
        gate. Coverage is presently 31 of 32 declared ranges on a cloudy column,
        the one exception being a variable a module integrates over rather than
        evaluates at, so a pass is no longer the much weaker statement it was.
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
    """The range a field actually spans, or None if there is nothing to judge.

    **A NaN means the relation was not evaluated there, and is skipped.** That
    is what lets a module narrow a check from the extent a quantity *spans* to
    the region the equation was *applied* to, which is otherwise this module's
    standing coarseness. The latent energy of vapour is the case that forced it:
    ``E = L_v sum_i w_i`` was flagged for reaching 225 K against a floor of 230,
    and the column holds exactly zero vapour below 230 K, so every temperature
    that broke the range multiplied a weight of zero.

    An infinity is *not* skipped. A masked-out entry is a statement the module
    is making on purpose; an infinity is a defect, and treating it as a mask
    would let the check pass by being broken.
    """
    array = np.asarray(values, dtype=float)
    if array.size == 0 or bool(np.any(np.isinf(array))):
        return None
    finite = array[~np.isnan(array)]
    if finite.size == 0:
        return None
    return (float(np.min(finite)), float(np.max(finite)))


class EvaluatedExtents:
    """The span every published field covered while a run proceeded.

    Fed one step's outputs at a time, it keeps a running minimum and maximum per
    field and nothing else, so a hundred thousand steps cost two floats each --
    which matters, because the alternative of keeping the fields themselves is
    how a 120000-step manifest reached four gigabytes.

    The per-step rules are ``_extent``'s, applied step by step: a NaN is a
    statement that the relation was not evaluated there and is skipped, an
    infinity is a defect that disqualifies the field permanently rather than
    being masked away, and a field that is NaN at every step never enters at all
    -- which is the correct report that its relation never ran.
    """

    __slots__ = ("_disqualified", "_spans")

    def __init__(self) -> None:
        self._spans: dict[str, tuple[float, float]] = {}
        self._disqualified: set[str] = set()

    def observe(self, fields: Mapping[str, Any]) -> None:
        """Fold one step's published quantities into the accumulated spans."""
        for field_id, field in fields.items():
            if field_id in self._disqualified:
                continue
            try:
                array = np.asarray(getattr(field, "values", None), dtype=float)
            except (TypeError, ValueError):
                # Not a numeric field. No range could have been stated against
                # it, so there is nothing here to accumulate.
                continue
            if array.size == 0:
                continue
            # This runs on every field on every step, so it takes the two-pass
            # route when it can: a plain min and max propagate a NaN, which is
            # the signal to pay for the mask rather than a reason to build one
            # every time. Most fields carry none.
            low = float(array.min())
            high = float(array.max())
            if low != low or high != high:
                finite = array[~np.isnan(array)]
                if finite.size == 0:
                    continue
                low = float(finite.min())
                high = float(finite.max())
            if low == -np.inf or high == np.inf:
                self._disqualified.add(field_id)
                self._spans.pop(field_id, None)
                continue
            previous = self._spans.get(field_id)
            if previous is not None:
                low = min(low, previous[0])
                high = max(high, previous[1])
            self._spans[field_id] = (low, high)

    def __contains__(self, field_id: str) -> bool:
        """Whether this field was seen at all, in either state.

        A disqualified field counts as seen: the judgement must use the
        accumulated verdict, not fall back to the last step's values and pass a
        field that carried an infinity earlier in the run.
        """
        return field_id in self._spans or field_id in self._disqualified

    def extent_of(self, field_id: str) -> tuple[float, float] | None:
        return self._spans.get(field_id)


def check_validity_domains(
    *,
    registry: Any,
    modules: Iterable[Any],
    fields: Mapping[str, Any],
    extents: EvaluatedExtents | None = None,
) -> DomainReport:
    """Judge a run against the declared domains of the equations it cites.

    ``registry`` supplies equation records by id, ``modules`` are the modules
    that actually ran, and ``fields`` are the quantities they published.

    ``extents`` carries what each field spanned over the whole run, and is what
    the ranges are judged against when it is given. ``fields`` still resolves
    the bindings -- symbols, units and ambiguity are properties of the
    contracts, not of any one step. Without it the judgement falls back to the
    values in ``fields`` alone, which is a single-state check and is what the
    unit tests of the resolution machinery want.
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

        # A symbol may be bound more than once. One equation's ``r_eff`` covers
        # both the droplet and the crystal radius, and checking whichever was
        # declared first would leave the other silently unjudged -- the same
        # half-answer this module refuses when it refuses an ambiguity.
        declared: dict[str, dict[str, list[tuple[str, str]]]] = {}
        for equation, bindings in getattr(contract, "equation_variables", ()) or ():
            per_symbol = declared.setdefault(equation, {})
            for symbol, name, stated in bindings:
                per_symbol.setdefault(symbol, []).append((name, stated))
        parameters = _parameters_of(module)

        for equation_id in getattr(contract, "governing_equation_ids", ()) or ():
            record = _lookup(registry, equation_id)
            if record is None:
                continue
            bound = declared.get(equation_id, {})
            for name, low, high, units in record.validity_domain.ranges:
                # A binding the module declared beats anything inferred: it is
                # the one statement here that cannot be a coincidence of naming.
                attempts = (
                    [_bound(binding, fields, parameters) for binding in bound[name]]
                    if name in bound
                    else [_resolve(name, units, in_scope)]
                )
                for resolved, ambiguity in attempts:
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
                            extents=extents,
                        )
                    )

    return DomainReport(tuple(checks))


def _parameters_of(module: Any) -> Mapping[str, Any]:
    """A module's own parameters, for variables that never become fields.

    Optical depths, orbital eccentricity and the bulk-plume parameter are all
    declared with ranges and none of them is published as a field, so without
    this they are permanently unresolvable.
    """
    getter = getattr(module, "parameters", None)
    if not callable(getter):
        return {}
    try:
        values = getter()
    except Exception:
        return {}
    return values if isinstance(values, Mapping) else {}


def _bound(
    binding: tuple[str, str],
    fields: Mapping[str, Any],
    parameters: Mapping[str, Any],
) -> tuple[tuple[str, Any] | None, str]:
    """Resolve a binding the module declared, as a field or as a parameter.

    The module states the units of what it bound, so a parameter -- which has
    none of its own -- is still judged on a stated basis rather than on the
    assumption that it matches whatever the equation declared.
    """
    name, stated = binding
    if name == INTEGRATED_OVER:
        return None, (
            "the module integrates over this variable rather than evaluating the "
            "relation at any value of it, so no range applies"
        )
    field = fields.get(name)
    if field is not None:
        return (name, _StatedField(field.values, stated)), ""
    if name in parameters:
        return (name, _StatedField(parameters[name], stated)), ""
    return None, f"the module binds this variable to {name!r}, which it did not publish"


@dataclass(frozen=True, slots=True)
class _StatedField:
    """A bound quantity carrying the units the module stated for it."""

    values: Any
    units: str
    symbol: str = ""


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
    extents: EvaluatedExtents | None = None,
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

    # A parameter is never accumulated -- it has no field to publish and does
    # not vary over the run -- so the fallback is not a special case for it.
    if extents is not None and field_id in extents:
        extent = extents.extent_of(field_id)
        over = " over the run"
    else:
        extent = _extent(getattr(field, "values", None))
        over = ""
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
            f"{name} spans [{lowest:g}, {highest:g}] {units} in {field_id}{over}, "
            f"outside the declared [{low:g}, {high:g}]",
        )
    return DomainCheck(
        module_id,
        equation_id,
        name,
        "inside",
        f"{name} spans [{lowest:g}, {highest:g}] {units} in {field_id}{over}, "
        f"within the declared [{low:g}, {high:g}]",
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
