"""What a module declares about its own uncertain parameters (UNC-001, UNC-002).

A module's ``report_uncertainty`` used to be able to say two things about a
parameter it knew was uncertain: a number it had quantified, or the word
UNRESOLVED. Neither is what the evaporating condensate fraction needed. It has a
nominal value, a range, and a *basis* for the range -- measured precipitation
efficiencies -- and it had been sitting in a code comment beside a constant,
where no ensemble could read it and nothing checked a sweep against it.

So there are three states a parameter's uncertainty can be in, and the report
says which:

``unresolved``
    Known to matter, not bounded. Still the honest answer for most things.
``DeclaredRange``
    Bounded, with the bounds' basis stated: measured, published, or asserted.
    Declared is not propagated -- a range says what to run, not what running it
    gives -- and an ensemble built from it is what propagates it.
``DerivedParameter``
    Not independently uncertain at all, because another parameter fixes it. The
    evaporation coefficient is ``k_H / c_p`` exactly; perturbing it on its own
    would reintroduce the fifty-fold inconsistency that was the whole Bowen-ratio
    error. That is a correlation of one, and it is the most useful correlation
    statement a module can make, because it forbids something.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Any

from mfespf.core.errors import ConfigurationError

__all__ = ["DeclaredRange", "DerivedParameter", "RangeBasis"]


class RangeBasis(StrEnum):
    """Where a range's bounds came from.

    Required because an asserted 0.1 to 0.5 and a measured 0.1 to 0.5 produce
    identical numbers. The evaporating fraction was asserted over exactly that
    range and it excluded the measured value; nothing in any sweep across it
    said the bounds were a judgement.
    """

    MEASURED = "measured"
    """The bounds come from observations or process-resolving studies."""
    PUBLISHED = "published"
    """A published model value and the spread published with it."""
    ASSERTED = "asserted"
    """Nobody measured it. The bounds are a judgement and the source says whose."""


@dataclass(frozen=True, slots=True)
class DeclaredRange:
    """A parameter's declared range and the basis for it.

    ``parameter`` is the key in the module's ``parameters()``, so the nominal can
    be checked against what the module actually runs with. ``keyword`` is the
    configuration keyword that sets it, which is what an ensemble overrides; it
    defaults to ``parameter`` and differs only where the two names do.

    The configured value is *not* required to lie inside the range. A run outside
    it is a legitimate exploration and the manifest says so; an *ensemble*
    reaching outside it is refused, because an ensemble claims a spread and a
    spread over a range nobody declared is a claim with no basis.
    """

    parameter: str
    nominal: float
    low: float
    high: float
    basis: RangeBasis
    source: str
    keyword: str = ""

    def __post_init__(self) -> None:
        if not self.high > self.low:
            raise ConfigurationError(
                f"declared range for {self.parameter!r} is [{self.low}, {self.high}]; a range "
                "must be increasing -- a parameter with no spread is not uncertain"
            )
        if not self.source.strip():
            raise ConfigurationError(
                f"declared range for {self.parameter!r} names no source. An asserted "
                "range still says who asserted it and why."
            )

    @property
    def setting(self) -> str:
        """The configuration keyword an ensemble sets this parameter through."""
        return self.keyword or self.parameter

    def contains(self, value: float) -> bool:
        return self.low <= value <= self.high

    @property
    def nominal_inside(self) -> bool:
        return self.contains(self.nominal)

    def describe(self) -> str:
        outside = "" if self.nominal_inside else " (CONFIGURED OUTSIDE ITS DECLARED RANGE)"
        return (
            f"{self.parameter} = {self.nominal:g} in [{self.low:g}, {self.high:g}], "
            f"{self.basis.value}: {self.source}{outside}"
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "parameter": self.parameter,
            "setting": self.setting,
            "nominal": self.nominal,
            "low": self.low,
            "high": self.high,
            "basis": self.basis.value,
            "source": self.source,
        }


@dataclass(frozen=True, slots=True)
class DerivedParameter:
    """A parameter fixed by another one, and so not to be perturbed alone."""

    parameter: str
    follows: str
    """The configuration keyword of the parameter that fixes this one."""
    relation: str
    keyword: str = ""

    def __post_init__(self) -> None:
        if not self.relation.strip():
            raise ConfigurationError(
                f"{self.parameter!r} is declared derived from {self.follows!r} with no "
                "relation; a dependency nobody can check is an assertion"
            )

    @property
    def setting(self) -> str:
        return self.keyword or self.parameter

    def describe(self) -> str:
        return f"{self.parameter} follows {self.follows}: {self.relation}"

    def to_dict(self) -> dict[str, Any]:
        return {
            "parameter": self.parameter,
            "setting": self.setting,
            "follows": self.follows,
            "relation": self.relation,
        }
