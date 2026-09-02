"""Equation and approximation registry.

Amendment A7 and requirements EQR-001 through EQR-007.

Every governing equation, constitutive relation, closure, empirical correction,
diagnostic identity and numerical operator used in a scientific run has a
versioned identity here.  EQR-004 is the reason ``category`` is mandatory: a
calibrated correction must be impossible to report as a fundamental law simply
by describing it differently in prose.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

from mfespf.core.errors import ConfigurationError
from mfespf.core.versioning import Version
from mfespf.registries.base import Registry

__all__ = ["EquationCategory", "EquationRecord", "EquationRegistry", "ValidityDomain"]


class EquationCategory(StrEnum):
    """A7.3 / EQR-003."""

    FUNDAMENTAL_LAW = "fundamental_law"
    CONSTITUTIVE_RELATION = "constitutive_relation"
    REDUCED_PHYSICAL_EQUATION = "reduced_physical_equation"
    CLOSURE_OR_PARAMETERIZATION = "closure_or_parameterization"
    NUMERICAL_OPERATOR = "numerical_operator"
    EMPIRICAL_RELATION = "empirical_relation"
    CALIBRATED_CORRECTION = "calibrated_correction"
    DIAGNOSTIC_IDENTITY = "diagnostic_identity"

    @property
    def is_fitted(self) -> bool:
        """Categories whose parameters came from data rather than from physics.

        A claim resting on one of these cannot be described as derived from first
        principles, and the data that fitted it cannot later be called
        independent validation (INF-003, VAL-006).
        """
        return self in (
            EquationCategory.EMPIRICAL_RELATION,
            EquationCategory.CALIBRATED_CORRECTION,
        )


@dataclass(frozen=True, slots=True)
class ValidityDomain:
    """Charter 4.26.  Ranges are inclusive and expressed in the equation's units."""

    ranges: tuple[tuple[str, float, float, str], ...] = ()
    conditions: tuple[str, ...] = ()

    def contains(self, **values: float) -> tuple[bool, tuple[str, ...]]:
        """Check named values against declared ranges.

        Returns ``(ok, violations)``.  A variable with no declared range is not a
        violation -- it is an undeclared dimension of the domain, which the
        equation record should say so in ``conditions``.
        """
        violations = []
        for name, low, high, _units in self.ranges:
            if name in values and not (low <= values[name] <= high):
                violations.append(f"{name}={values[name]:g} outside [{low:g}, {high:g}]")
        return (not violations, tuple(violations))

    def describe(self) -> str:
        parts = [f"{n} in [{lo:g}, {hi:g}] {u}" for n, lo, hi, u in self.ranges]
        return "; ".join(parts + list(self.conditions)) or "undeclared"


@dataclass(frozen=True)
class EquationRecord:
    """A7.2 / EQR-002."""

    record_id: str
    version: Version
    name: str
    category: EquationCategory
    mathematical_form: str
    physical_basis: str
    variables: tuple[tuple[str, str, str], ...] = ()  # (symbol, meaning, units)
    parameters: tuple[tuple[str, str, str], ...] = ()  # (name, parameter_class, source)
    assumptions: tuple[str, ...] = ()
    approximations: tuple[str, ...] = ()
    validity_domain: ValidityDomain = ValidityDomain()
    limiting_cases: tuple[str, ...] = ()
    conserved_quantities: tuple[str, ...] = ()
    symmetries: tuple[str, ...] = ()
    permitted_discretizations: tuple[str, ...] = ()
    correctness_tests: tuple[str, ...] = ()
    known_failure_modes: tuple[str, ...] = ()
    source: str = ""
    using_modules: tuple[str, ...] = ()
    supersedes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.mathematical_form.strip():
            raise ConfigurationError(
                f"equation {self.record_id!r} has no mathematical form (EQR-002)"
            )
        if self.category is EquationCategory.FUNDAMENTAL_LAW and any(
            pc == "calibrated" for _n, pc, _s in self.parameters
        ):
            raise ConfigurationError(
                f"equation {self.record_id!r} is declared FUNDAMENTAL_LAW but carries a calibrated "
                f"parameter; EQR-004 forbids labelling a calibrated correction as a fundamental law"
            )
        if not self.correctness_tests and self.category is not EquationCategory.DIAGNOSTIC_IDENTITY:
            raise ConfigurationError(
                f"equation {self.record_id!r} declares no correctness test (VAL-001); a "
                f"diagnostic identity is the only category exempt"
            )

    @property
    def equation_id(self) -> str:
        return self.record_id

    @property
    def calibrated_parameters(self) -> tuple[str, ...]:
        return tuple(n for n, pc, _s in self.parameters if pc == "calibrated")

    def cite(self) -> str:
        return f"{self.record_id}@{self.version}"


class EquationRegistry(Registry[EquationRecord]):
    def __init__(self) -> None:
        super().__init__("equation")

    def bind(self, equation_ids: dict[str, str]) -> dict[str, str]:
        """Resolve ``{equation_id: version}`` into manifest bindings (EQR-005, A7.4).

        Raises if any binding names an unregistered record, so a run cannot start
        while claiming to use an equation nobody wrote down.
        """
        out = {}
        for eq_id, version in sorted(equation_ids.items()):
            entry = self.entry(eq_id, version)
            out[f"{eq_id}@{version}"] = entry.content_hash
        return out

    def fitted_equations(self) -> list[EquationRecord]:
        """Every active equation whose parameters came from data (for INF-004)."""
        return [e for e in self.current() if e.category.is_fitted or e.calibrated_parameters]
