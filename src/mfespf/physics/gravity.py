"""Gravitational acceleration, through a replaceable interface (GRV-001).

Charter 9.15 and requirement GRV-001: gravity shall be provided through a
replaceable interface beginning with a declared approximation. Until this module
existed it was provided by ``value_of("earth_standard_gravity")``, read directly
from the constants table at four call sites, and the constant's own note said
that the gravity module supersedes it.

That mattered for the reason the requirement exists. A constant read out of a
table is an approximation nobody has to declare, so nobody has to measure what it
costs -- and the hydrostatic equation's record asserted that constant gravity is
"about 1.5 percent at 100 km and negligible through the troposphere". The second
half is true and is about a column this framework does not run. **This column's
top is at 1000 Pa, which is 19 to 27 km up depending on the resolution**, and the
measured cost there is 0.21 percent of the column's mass and 0.6 to 0.85 percent
of the topmost layer's -- the layer that sets the outgoing longwave.

The fidelity path is charter 9.15's, and two of its four rungs are here:

1. ``ConstantGravity`` -- the conventional value, everywhere. What the framework
   has always done, now saying so.
2. ``NewtonianAltitudeGravity`` -- ``g(z) = g_0 (R / (R + z))^2``, the inverse
   square law about a point mass.
3. A spatial gravitational field, and
4. a relativistic correction where a defined claim requires it, are carve-outs
   in the charter's sense: no configuration here has a horizontal coordinate to
   vary gravity over, and nothing this framework predicts is sensitive at the
   scale of a relativistic correction.

**Latitude is deliberately not implemented**, although charter 9.15 pairs it with
altitude on the second rung. Latitudinal variation of gravity is about 0.5
percent from equator to pole, which is the same size as the altitude effect
measured above -- so it is not omitted for being small. It is omitted because
every configuration in this framework is a global mean or a single column, and
neither has a latitude to evaluate it at. Supplying one would be inventing a
coordinate the physics does not have, and the conventional value is already a
latitude-weighted compromise. This is a case where the honest implementation is
the missing one; see ``LATITUDE_NOT_REPRESENTED``.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

import numpy as np

from mfespf.contracts.module import FidelityLevel
from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError

__all__ = [
    "LATITUDE_NOT_REPRESENTED",
    "ConstantGravity",
    "GravityField",
    "NewtonianAltitudeGravity",
    "as_gravity_field",
]

#: Why the second rung of charter 9.15's fidelity path is implemented in altitude
#: and not in latitude. Quoted by the field's declared approximations so the
#: omission travels with the run rather than living only in this docstring.
LATITUDE_NOT_REPRESENTED = (
    "latitude is not represented: every configuration here is a global mean or a "
    "single column, which has no latitude to evaluate gravity at, and the "
    "conventional value is already a latitude-weighted compromise"
)

#: Mean Earth radius. Used only to place a column above the surface, so the mean
#: is the right choice rather than the equatorial or polar value -- which differ
#: by 0.3 percent and would be false precision beside a column that has no
#: latitude.
EARTH_MEAN_RADIUS_M = 6.371e6


class GravityField(ABC):
    """A replaceable source of gravitational acceleration.

    Deliberately *not* a ``Module``. A module owns state, advances it and takes
    part in the coupling; gravity does none of those things -- it is a property of
    the domain that other modules consume as a parameter, and giving it a
    lifecycle would put a participant in every dependency graph that contributes
    nothing to any budget. What GRV-001 asks for is that the approximation be
    declared and the implementation replaceable, and both are satisfied by an
    interface with more than one implementation behind it.

    Every implementation declares its fidelity level, its approximations and the
    conditions it holds under, and publishes its parameters so two runs that
    differ in gravity cannot hash alike (AUD-002).
    """

    fidelity_level: FidelityLevel
    approximations: tuple[str, ...]
    validity_domain: str

    @abstractmethod
    def at_height(self, height_m: Any) -> np.ndarray:
        """Gravitational acceleration at each height above the surface, in m s^-2."""

    @abstractmethod
    def parameters(self) -> dict[str, Any]:
        """Everything that determines this field's numbers."""

    @property
    def varies_with_height(self) -> bool:
        """Whether a caller must supply heights for this field to be evaluated.

        A caller that has no heights and a field that needs them is a refusal,
        not a silent fallback to the surface value: the whole point of declaring
        the approximation is that it cannot be applied by accident.
        """
        return True

    def surface_value(self) -> float:
        return float(np.asarray(self.at_height(0.0)).reshape(-1)[0])

    def describe(self) -> str:
        return (
            f"{type(self).__name__} at {self.fidelity_level.value}: "
            f"{self.surface_value():.5f} m s^-2 at the surface"
        )


class ConstantGravity(GravityField):
    """One value everywhere: charter 9.15's first rung.

    The conventional CGPM 1901 value by default, which is exact by definition and
    not a measurement of anywhere. Its error is not in the number but in the
    assumption that it applies at every height in the column.
    """

    fidelity_level = FidelityLevel.REDUCED
    validity_domain = (
        "Any column whose depth is small against the Earth's radius, where the "
        "z / R error is smaller than the accuracy claimed of anything computed "
        "from the layer masses."
    )

    def __init__(self, value: float | None = None) -> None:
        self.value = float(value if value is not None else value_of("earth_standard_gravity"))
        if self.value <= 0.0:
            raise ConfigurationError(f"gravity {self.value} must be positive")
        self.approximations = (
            "gravity constant over the depth of the column; the error is of order "
            "2 z / R_earth and is measured at 0.21 percent of this column's mass "
            "and 0.85 percent of its topmost layer at 32 layers over 1000 Pa",
            LATITUDE_NOT_REPRESENTED,
        )

    @property
    def varies_with_height(self) -> bool:
        return False

    def at_height(self, height_m: Any) -> np.ndarray:
        return np.full_like(np.asarray(height_m, dtype=float), self.value)

    def parameters(self) -> dict[str, Any]:
        return {"gravity_model": "constant", "surface_gravity_m_s2": self.value}


class NewtonianAltitudeGravity(GravityField):
    """``g(z) = g_0 (R / (R + z))^2``: charter 9.15's second rung, in altitude.

    The inverse square law outside a spherically symmetric mass. Two effects are
    folded into the surface value rather than modelled: the Earth's oblateness
    and its rotation, both of which are latitudinal and neither of which a column
    without a latitude can represent. So this is the altitude dependence alone,
    anchored at the conventional value -- which is what makes it a strictly
    better statement than a constant rather than a different one.
    """

    fidelity_level = FidelityLevel.STANDARD
    validity_domain = (
        "Heights above the surface small against the Earth's radius and outside "
        "the solid body, to a few hundred kilometres. Says nothing about "
        "latitude, so it is exact only for the latitude the surface value is "
        "conventional for."
    )

    def __init__(
        self,
        *,
        surface_gravity: float | None = None,
        earth_radius: float = EARTH_MEAN_RADIUS_M,
    ) -> None:
        self.surface_gravity = float(
            surface_gravity if surface_gravity is not None else value_of("earth_standard_gravity")
        )
        self.earth_radius = float(earth_radius)
        if self.surface_gravity <= 0.0:
            raise ConfigurationError(f"surface gravity {self.surface_gravity} must be positive")
        if self.earth_radius <= 0.0:
            raise ConfigurationError(f"earth radius {self.earth_radius} must be positive")
        self.approximations = (
            "the Earth is a point mass: oblateness and rotation are absorbed into "
            "the surface value rather than modelled, both being latitudinal",
            LATITUDE_NOT_REPRESENTED,
        )

    def at_height(self, height_m: Any) -> np.ndarray:
        z = np.asarray(height_m, dtype=float)
        if np.any(z <= -self.earth_radius):
            raise ConfigurationError("a height below the centre of the Earth has no gravity")
        ratio = self.earth_radius / (self.earth_radius + z)
        return self.surface_gravity * ratio * ratio

    def parameters(self) -> dict[str, Any]:
        return {
            "gravity_model": "newtonian_altitude",
            "surface_gravity_m_s2": self.surface_gravity,
            "earth_radius_m": self.earth_radius,
        }


def as_gravity_field(gravity: float | GravityField | None) -> GravityField:
    """Accept a field, a bare number, or nothing, and return a declared field.

    A bare number stays valid because it is what four call sites and a good many
    tests already pass, and because a constant *is* the first rung of the fidelity
    path rather than a way around it. What changes is that it now arrives wrapped
    in something that states its approximation.
    """
    if isinstance(gravity, GravityField):
        return gravity
    return ConstantGravity(gravity)
