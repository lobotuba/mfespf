"""Exchanged-field data model.

Charter section 7 and requirements DAT-001, DAT-002, DAT-006 and DAT-007.  A
``Field`` is the only thing that crosses a module boundary.  It carries its full
metadata with it, so a receiver never has to infer units, sign convention, time
semantics or spatial support from context -- the failure mode the charter's
interface-failure taxonomy (section 15.3) is written against.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from dataclasses import field as dc_field
from enum import StrEnum
from typing import Any

import numpy as np

from mfespf.core.errors import ConfigurationError, DimensionalityError
from mfespf.core.units import Unit, parse_unit, require_compatible

__all__ = [
    "Field",
    "MissingValuePolicy",
    "Provenance",
    "QuantityKind",
    "SignConvention",
    "SpatialSupport",
    "TimeSemantics",
    "Uncertainty",
]


class QuantityKind(StrEnum):
    """Charter 7.2.  A receiver treats these differently under regridding and
    time aggregation, so the distinction is part of the contract, not a label."""

    INTENSIVE_STATE = "intensive_state"
    EXTENSIVE_STATE = "extensive_state"
    FLUX = "flux"
    SOURCE_RATE = "source_rate"
    SINK_RATE = "sink_rate"
    TENDENCY = "tendency"
    BOUNDARY_VALUE = "boundary_value"
    DIAGNOSTIC = "diagnostic"
    PROBABILITY_DISTRIBUTION = "probability_distribution"
    UNCERTAINTY_FIELD = "uncertainty_field"
    CATEGORICAL_STATE = "categorical_state"

    @property
    def is_extensive(self) -> bool:
        """Extensive quantities and fluxes must regrid conservatively (SPA-003)."""
        return self in (
            QuantityKind.EXTENSIVE_STATE,
            QuantityKind.FLUX,
            QuantityKind.SOURCE_RATE,
            QuantityKind.SINK_RATE,
        )


class TimeSemantics(StrEnum):
    """Charter 7.5 / DAT-006."""

    INSTANTANEOUS = "instantaneous"
    INTERVAL_AVERAGE = "interval_average"
    INTERVAL_ACCUMULATION = "interval_accumulation"
    START_OF_STEP = "start_of_step"
    MIDPOINT = "midpoint"
    END_OF_STEP = "end_of_step"
    PREDICTOR = "predictor"
    CORRECTED = "corrected"

    @property
    def requires_interval(self) -> bool:
        return self in (
            TimeSemantics.INTERVAL_AVERAGE,
            TimeSemantics.INTERVAL_ACCUMULATION,
        )


class SpatialSupport(StrEnum):
    """SPA-002."""

    GLOBAL_SCALAR = "global_scalar"
    POINT = "point"
    CELL_AVERAGE = "cell_average"
    FACE_CENTERED = "face_centered"
    LAYER_AVERAGE = "layer_average"
    DOMAIN_INTEGRAL = "domain_integral"


class MissingValuePolicy(StrEnum):
    """DAT-007.  ``ZERO_FILL`` is absent on purpose: silently replacing a missing
    physical input with zero is exactly what the requirement forbids."""

    NOT_PERMITTED = "not_permitted"
    PROPAGATE_NAN = "propagate_nan"
    EXPLICIT_SENTINEL = "explicit_sentinel"
    MASKED = "masked"


class SignConvention(StrEnum):
    """CON-003.  The default for transfers is positive into the receiver."""

    POSITIVE_INTO_RECEIVER = "positive_into_receiver"
    POSITIVE_UPWARD = "positive_upward"
    POSITIVE_DOWNWARD = "positive_downward"
    POSITIVE_DEFINITE = "positive_definite"
    NOT_APPLICABLE = "not_applicable"


@dataclass(frozen=True, slots=True)
class Provenance:
    """Where a value came from (DAT-001, DAT-009).

    ``producer`` is a module_id for computed fields and a dataset_id for
    prescribed ones; ``derivation`` names the equation_id when the value was
    computed, so a diagnostic can be traced back to the equation registry.
    """

    producer: str
    derivation: str | None = None
    dataset_id: str | None = None
    dataset_version: str | None = None
    notes: str = ""


@dataclass(frozen=True, slots=True)
class Uncertainty:
    """Charter 7.6 / UNC-002, UNC-004.

    ``unresolved`` names applicable uncertainty contributions the producer knows
    it has not quantified.  It is required rather than optional so that a silent
    omission and a declared gap are distinguishable in a report.
    """

    standard_uncertainty: np.ndarray | float | None = None
    lower_bound: np.ndarray | float | None = None
    upper_bound: np.ndarray | float | None = None
    method: str | None = None
    distribution: str | None = None
    covariance_ref: str | None = None
    unresolved: tuple[str, ...] = ()

    @property
    def is_quantified(self) -> bool:
        return self.standard_uncertainty is not None or self.lower_bound is not None


@dataclass(frozen=True)
class Field:
    """A quantity crossing a module interface, with its complete metadata.

    Every DAT-001 element is a required constructor argument or has a value whose
    default is itself a declaration (for example ``missing_value_policy`` defaults
    to ``NOT_PERMITTED``, the strictest choice).
    """

    field_id: str
    physical_name: str
    values: np.ndarray
    units: str
    quantity_kind: QuantityKind
    time_semantics: TimeSemantics
    spatial_support: SpatialSupport
    provenance: Provenance
    symbol: str = ""
    coordinate_system: str = "global"
    dimensions: tuple[str, ...] = ()
    timestamp: float | None = None
    time_interval: float | None = None
    averaging_method: str | None = None
    spectral_coordinate: str | None = None
    angular_coordinate: str | None = None
    phase_or_species: str | None = None
    sign_convention: SignConvention = SignConvention.NOT_APPLICABLE
    missing_value_policy: MissingValuePolicy = MissingValuePolicy.NOT_PERMITTED
    uncertainty: Uncertainty | None = None
    quality_flags: tuple[str, ...] = ()
    stub: bool = False
    metadata: dict[str, Any] = dc_field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "values", np.asarray(self.values, dtype=float))
        parse_unit(self.units)  # resolve now so an unknown unit fails at construction
        if self.time_semantics.requires_interval and self.time_interval is None:
            raise ConfigurationError(
                f"field {self.field_id!r} declares {self.time_semantics.value} "
                f"but no time_interval (DAT-006)"
            )
        if self.missing_value_policy is MissingValuePolicy.NOT_PERMITTED and not np.all(
            np.isfinite(self.values)
        ):
            raise ConfigurationError(
                f"field {self.field_id!r} contains non-finite values but declares "
                f"missing_value_policy=not_permitted (DAT-007)"
            )
        if self.dimensions and self.values.ndim != len(self.dimensions):
            raise ConfigurationError(
                f"field {self.field_id!r} declares dimensions {self.dimensions} "
                f"but its array has rank {self.values.ndim}"
            )

    @property
    def unit(self) -> Unit:
        return parse_unit(self.units)

    @property
    def is_scalar(self) -> bool:
        return self.values.ndim == 0 or self.values.size == 1

    def scalar(self) -> float:
        """Return the single value of a scalar field, or raise."""
        if not self.is_scalar:
            raise ConfigurationError(
                f"field {self.field_id!r} has {self.values.size} values; scalar() needs one"
            )
        return float(self.values.reshape(-1)[0])

    def require_units(self, expected: str) -> Field:
        """Assert a unit symbol exactly, for a receiver that will not convert.

        Dimensional compatibility alone is not enough here: a receiver that
        multiplies by a hard-coded constant needs the scale to match too.
        """
        if self.units != expected:
            require_compatible(
                self.unit,
                parse_unit(expected),
                context=f"field {self.field_id!r}",
            )
            raise DimensionalityError(
                f"field {self.field_id!r} is in {self.units!r} but {expected!r} was required; "
                f"convert at a declared boundary (DAT-004) rather than here"
            )
        return self

    def with_values(self, values: np.ndarray, **changes: Any) -> Field:
        """Return a copy carrying new values, preserving all declared metadata."""
        return replace(self, values=np.asarray(values, dtype=float), **changes)

    def describe(self) -> str:
        marker = " [STUB - NOT PHYSICAL]" if self.stub else ""
        return (
            f"{self.field_id} ({self.physical_name}) "
            f"[{self.units}] {self.quantity_kind.value}/{self.time_semantics.value}"
            f"{marker}"
        )
