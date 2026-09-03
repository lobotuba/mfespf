"""Dimension and unit system.

Implements KER-002 (unit and dimension validator) and requirements DAT-003
through DAT-005: internal interfaces carry SI units, conversion happens only at
controlled boundaries and is logged, and dimensionally incompatible fields fail
validation rather than being silently coerced.

The system is deliberately small and dependency-free.  A physical dimension is a
vector of seven rational exponents over the SI base quantities; a unit is a named
affine map onto the coherent SI unit of its dimension.  That is enough to catch
every error class the charter cares about -- wrong quantity, wrong scale, silent
coercion -- without importing a general-purpose units library whose behaviour we
would then have to characterise for the correctness record.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import ClassVar

from mfespf.core.errors import DimensionalityError, UnitError

__all__ = [
    "DIMENSIONLESS",
    "Dimension",
    "Unit",
    "UnitConversion",
    "convert",
    "parse_unit",
    "require_compatible",
]

# Order is fixed and normative: it appears in serialised manifests.
_BASE_SYMBOLS = ("kg", "m", "s", "K", "A", "mol", "cd")


@dataclass(frozen=True, slots=True)
class Dimension:
    """A vector of rational exponents over the seven SI base quantities."""

    exponents: tuple[Fraction, ...]

    def __post_init__(self) -> None:
        if len(self.exponents) != 7:
            raise DimensionalityError(
                f"a dimension needs exactly 7 exponents, got {len(self.exponents)}"
            )

    @classmethod
    def of(
        cls,
        *,
        mass: float | Fraction = 0,
        length: float | Fraction = 0,
        time: float | Fraction = 0,
        temperature: float | Fraction = 0,
        current: float | Fraction = 0,
        amount: float | Fraction = 0,
        luminous: float | Fraction = 0,
    ) -> Dimension:
        return cls(
            tuple(
                Fraction(v).limit_denominator(1000)
                for v in (mass, length, time, temperature, current, amount, luminous)
            )
        )

    def __mul__(self, other: Dimension) -> Dimension:
        return Dimension(tuple(a + b for a, b in zip(self.exponents, other.exponents, strict=True)))

    def __truediv__(self, other: Dimension) -> Dimension:
        return Dimension(tuple(a - b for a, b in zip(self.exponents, other.exponents, strict=True)))

    def __pow__(self, power: float | Fraction) -> Dimension:
        p = Fraction(power).limit_denominator(1000)
        return Dimension(tuple(a * p for a in self.exponents))

    @property
    def is_dimensionless(self) -> bool:
        return all(e == 0 for e in self.exponents)

    def symbol(self) -> str:
        """Canonical ``kg^a m^b ...`` string, stable enough to serialise."""
        if self.is_dimensionless:
            return "1"
        parts = []
        for sym, exp in zip(_BASE_SYMBOLS, self.exponents, strict=True):
            if exp == 0:
                continue
            parts.append(sym if exp == 1 else f"{sym}^{exp}")
        return " ".join(parts)

    def named(self) -> str:
        """Human-readable name when one is registered, else the symbol."""
        return _DIMENSION_NAMES.get(self, self.symbol())

    def __str__(self) -> str:  # pragma: no cover - convenience
        return self.symbol()


DIMENSIONLESS = Dimension.of()

# Named dimensions used by the physical interfaces of charter section 9.
MASS = Dimension.of(mass=1)
LENGTH = Dimension.of(length=1)
TIME = Dimension.of(time=1)
TEMPERATURE = Dimension.of(temperature=1)
AREA = Dimension.of(length=2)
VOLUME = Dimension.of(length=3)
VELOCITY = Dimension.of(length=1, time=-1)
ACCELERATION = Dimension.of(length=1, time=-2)
FORCE = Dimension.of(mass=1, length=1, time=-2)
PRESSURE = Dimension.of(mass=1, length=-1, time=-2)
ENERGY = Dimension.of(mass=1, length=2, time=-2)
POWER = Dimension.of(mass=1, length=2, time=-3)
IRRADIANCE = Dimension.of(mass=1, time=-3)
SPECTRAL_IRRADIANCE = Dimension.of(mass=1, length=-1, time=-3)
AREAL_ENERGY = Dimension.of(mass=1, time=-2)
AREAL_HEAT_CAPACITY = Dimension.of(mass=1, time=-2, temperature=-1)
MOMENTUM = Dimension.of(mass=1, length=1, time=-1)
CHARGE = Dimension.of(time=1, current=1)
SUBSTANCE = Dimension.of(amount=1)
MASS_FLUX = Dimension.of(mass=1, length=-2, time=-1)

_DIMENSION_NAMES: dict[Dimension, str] = {
    DIMENSIONLESS: "dimensionless",
    MASS: "mass",
    LENGTH: "length",
    TIME: "time",
    TEMPERATURE: "temperature",
    AREA: "area",
    VOLUME: "volume",
    VELOCITY: "velocity",
    ACCELERATION: "acceleration",
    FORCE: "force",
    PRESSURE: "pressure",
    ENERGY: "energy",
    POWER: "power",
    IRRADIANCE: "irradiance",
    SPECTRAL_IRRADIANCE: "spectral_irradiance",
    AREAL_ENERGY: "areal_energy",
    AREAL_HEAT_CAPACITY: "areal_heat_capacity",
    MOMENTUM: "momentum",
    CHARGE: "electric_charge",
    SUBSTANCE: "amount_of_substance",
    MASS_FLUX: "mass_flux",
}


@dataclass(frozen=True, slots=True)
class Unit:
    """An affine map from a unit onto the coherent SI unit of its dimension.

    ``si_value = value * scale + offset``.  ``offset`` is non-zero only for the
    conventional temperature scales.  Affine units are refused inside arithmetic
    because sums and products of them are not physically meaningful; they exist
    so that an input dataset in degrees Celsius can be converted once, at a
    logged boundary, and never again.
    """

    symbol: str
    dimension: Dimension
    scale: float = 1.0
    offset: float = 0.0

    registry: ClassVar[dict[str, Unit]] = {}

    @property
    def is_si_coherent(self) -> bool:
        return self.scale == 1.0 and self.offset == 0.0

    @property
    def is_affine(self) -> bool:
        return self.offset != 0.0

    def to_si(self, value: float) -> float:
        return value * self.scale + self.offset

    def from_si(self, si_value: float) -> float:
        return (si_value - self.offset) / self.scale

    def __str__(self) -> str:  # pragma: no cover - convenience
        return self.symbol


def _register(symbol: str, dimension: Dimension, scale: float = 1.0, offset: float = 0.0) -> Unit:
    unit = Unit(symbol=symbol, dimension=dimension, scale=scale, offset=offset)
    Unit.registry[symbol] = unit
    return unit


# --- Coherent SI units -------------------------------------------------------
DIMENSIONLESS_UNIT = _register("1", DIMENSIONLESS)
_register("kg", MASS)
_register("m", LENGTH)
_register("s", TIME)
_register("K", TEMPERATURE)
_register("A", Dimension.of(current=1))
_register("mol", SUBSTANCE)
_register("cd", Dimension.of(luminous=1))
_register("m^2", AREA)
_register("m^3", VOLUME)
_register("m s^-1", VELOCITY)
_register("m s^-2", ACCELERATION)
_register("N", FORCE)
_register("Pa", PRESSURE)
_register("J", ENERGY)
_register("W", POWER)
_register("W m^-2", IRRADIANCE)
_register("W m^-2 m^-1", SPECTRAL_IRRADIANCE)
_register("W m^-3", Dimension.of(mass=1, length=-1, time=-3))
_register("J m^-2", AREAL_ENERGY)
_register("J m^-2 K^-1", AREAL_HEAT_CAPACITY)
_register("J kg^-1 K^-1", Dimension.of(length=2, time=-2, temperature=-1))
_register("J kg^-1", Dimension.of(length=2, time=-2))
_register("kg m s^-1", MOMENTUM)
_register("C", CHARGE)
_register("kg m^-2 s^-1", MASS_FLUX)
_register("kg m^-2", Dimension.of(mass=1, length=-2))
_register("kg kg^-1", DIMENSIONLESS)
_register("K s^-1", Dimension.of(temperature=1, time=-1))
_register("W m^-2 K^-4", Dimension.of(mass=1, time=-3, temperature=-4))
_register("W m^-2 K^-1", Dimension.of(mass=1, time=-3, temperature=-1))
_register("J s", Dimension.of(mass=1, length=2, time=-1))
_register("J K^-1", Dimension.of(mass=1, length=2, time=-2, temperature=-1))
_register("kg m^-3", Dimension.of(mass=1, length=-3))
_register("m^3 kg^-1 s^-2", Dimension.of(mass=-1, length=3, time=-2))
_register("rad", DIMENSIONLESS)

# --- Accepted non-coherent units, converted only at boundaries (DAT-004) -----
_register("km", LENGTH, scale=1.0e3)
_register("nm", LENGTH, scale=1.0e-9)
_register("um", LENGTH, scale=1.0e-6)
_register("cm^-1", Dimension.of(length=-1), scale=1.0e2)
_register("W m^-2 nm^-1", SPECTRAL_IRRADIANCE, scale=1.0e9)
_register("W m^-2 um^-1", SPECTRAL_IRRADIANCE, scale=1.0e6)
_register("hPa", PRESSURE, scale=1.0e2)
_register("bar", PRESSURE, scale=1.0e5)
_register("minute", TIME, scale=60.0)
_register("hour", TIME, scale=3600.0)
_register("day", TIME, scale=86400.0)
_register("year", TIME, scale=365.25 * 86400.0)  # Julian year, declared explicitly
_register("au", LENGTH, scale=1.495978707e11)
_register("degC", TEMPERATURE, scale=1.0, offset=273.15)
_register("percent", DIMENSIONLESS, scale=1.0e-2)


def parse_unit(symbol: str) -> Unit:
    """Look up a registered unit symbol.

    An unknown symbol is an error rather than an assumption (DAT-005): a unit the
    framework cannot reason about cannot be validated, so it must not reach an
    interface.
    """
    try:
        return Unit.registry[symbol]
    except KeyError:
        raise UnitError(
            f"unit {symbol!r} is not registered; register it in mfespf.core.units "
            f"with an explicit dimension and scale before using it on an interface"
        ) from None


@dataclass(frozen=True, slots=True)
class UnitConversion:
    """A logged conversion event (DAT-004).

    Conversions are recorded rather than performed silently so that the run log
    can show every place a non-SI quantity crossed into the framework.
    """

    from_unit: str
    to_unit: str
    factor: float
    offset: float
    context: str

    def describe(self) -> str:
        tail = f" {self.offset:+g}" if self.offset else ""
        return f"{self.context}: {self.from_unit} -> {self.to_unit} (x{self.factor:g}{tail})"


def require_compatible(a: Unit, b: Unit, *, context: str) -> None:
    """Raise unless two units share a dimension.

    This is the DAT-005 gate.  It is a hard failure by design: there is no
    coercion path and no tolerance on dimensions.
    """
    if a.dimension != b.dimension:
        raise DimensionalityError(
            f"{context}: incompatible dimensions "
            f"{a.symbol} [{a.dimension.named()}] vs {b.symbol} [{b.dimension.named()}]"
        )


def convert(
    value: float, source: Unit, target: Unit, *, context: str
) -> tuple[float, UnitConversion]:
    """Convert a value, returning the result together with its conversion record."""
    require_compatible(source, target, context=context)
    si = source.to_si(value)
    out = target.from_si(si)
    record = UnitConversion(
        from_unit=source.symbol,
        to_unit=target.symbol,
        factor=source.scale / target.scale,
        offset=source.offset - target.offset,
        context=context,
    )
    return out, record
