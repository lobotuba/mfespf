"""Forcing protocols for timescale experiments.

Charter 11.3 and requirement TIM-007, which asks for step, pulse, impulse and
periodic forcing.

Two decisions shape this module.

**A forcing declares the energy it delivers, not the value it takes at an
instant.**  Every forcing implements ``mean_amplitude(t0, t1)``, the exact
average over the interval a module is about to integrate, obtained by
integrating the forcing analytically rather than by sampling it.  A pulse whose
edges fall inside a step therefore delivers exactly its declared energy at any
step size, and a periodic forcing delivers exactly its declared cycle mean.
Sampling at the step start instead would make the delivered energy a function of
the time step, which turns a convergence study into a forcing study.

**An impulse in discrete time is a pulse one step wide.**  There is no shorter
one.  ``ImpulseForcing`` therefore fixes the *energy* and lets the flux follow
from the step, so the impulse response is reproducible across time steps even
though the shape of the impulse is not resolved below one step.  That is stated
in the forcing's own declaration rather than left for a reader to discover.

The forcing is applied to total solar irradiance, so the energy it adds enters
through the top-of-atmosphere boundary and is accounted by the existing
transfers. Nothing here creates energy inside the system, and no configuration
using it needs a fiat source term.
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any

from mfespf.core.errors import ConfigurationError

__all__ = [
    "Forcing",
    "ImpulseForcing",
    "NoForcing",
    "PeriodicForcing",
    "PulseForcing",
    "StepForcing",
]


class Forcing(ABC):
    """A declared perturbation to a boundary condition.

    ``regime`` names the experiment class the charter lists, and travels with
    every timescale measured under it: TIM-008 requires timescales to be reported
    per forcing regime, which is only possible if the regime is carried rather
    than remembered.
    """

    regime: str = "undeclared"

    @abstractmethod
    def mean_amplitude(self, start: float, end: float) -> float:
        """Exact mean of the forcing over ``[start, end]``, in W m^-2.

        A zero-length interval returns the instantaneous value, which is what a
        module asks for outside a step.
        """

    @abstractmethod
    def parameters(self) -> dict[str, Any]:
        """The declaration that enters the run manifest and its hash."""

    @abstractmethod
    def describe(self) -> str:
        """One line naming the experiment, for a report a person reads."""

    def energy_delivered(self, start: float, end: float) -> float:
        """Integrated forcing over an interval, J m^-2 -- the invariant quantity."""
        return self.mean_amplitude(start, end) * (end - start)

    def _check_interval(self, start: float, end: float) -> None:
        if end < start:
            raise ConfigurationError(f"forcing interval [{start}, {end}] ends before it starts")


@dataclass(frozen=True)
class NoForcing(Forcing):
    """The null control: a forcing that is identically zero.

    Present as a class rather than as ``None`` so that an unforced run declares
    an experiment regime like any other. VAL-001 wants the null case run through
    the same path as the perturbed one, not around it.
    """

    regime: str = "none"

    def mean_amplitude(self, start: float, end: float) -> float:
        self._check_interval(start, end)
        return 0.0

    def parameters(self) -> dict[str, Any]:
        return {"regime": self.regime}

    def describe(self) -> str:
        return "no forcing (null control)"


@dataclass(frozen=True)
class StepForcing(Forcing):
    """A permanent change of ``amplitude`` applied from ``onset`` onward.

    The classical experiment for a relaxation timescale: the system leaves one
    equilibrium and approaches another, and the approach carries the timescales.
    """

    amplitude: float
    onset: float = 0.0
    regime: str = "step"

    def __post_init__(self) -> None:
        if self.onset < 0.0:
            raise ConfigurationError("step onset must not precede the start of the run")

    def mean_amplitude(self, start: float, end: float) -> float:
        self._check_interval(start, end)
        if end == start:
            return self.amplitude if start >= self.onset else 0.0
        active = max(0.0, end - max(start, self.onset))
        return self.amplitude * active / (end - start)

    def parameters(self) -> dict[str, Any]:
        return {"regime": self.regime, "amplitude_W_m2": self.amplitude, "onset_s": self.onset}

    def describe(self) -> str:
        return f"step of {self.amplitude:+g} W m^-2 from t = {self.onset:g} s"


@dataclass(frozen=True)
class PulseForcing(Forcing):
    """A rectangular pulse of declared amplitude and duration.

    The system leaves and then returns to the same equilibrium, which is a
    stronger test than a step: a step response can look right while energy is
    quietly lost, whereas a pulse must give back everything it delivered.
    """

    amplitude: float
    onset: float = 0.0
    duration: float = 0.0
    regime: str = "pulse"

    def __post_init__(self) -> None:
        if self.duration <= 0.0:
            raise ConfigurationError(
                "a pulse must declare a positive duration; for the shortest resolvable "
                "perturbation use ImpulseForcing, which fixes the energy instead"
            )
        if self.onset < 0.0:
            raise ConfigurationError("pulse onset must not precede the start of the run")

    @property
    def offset(self) -> float:
        return self.onset + self.duration

    def mean_amplitude(self, start: float, end: float) -> float:
        self._check_interval(start, end)
        if end == start:
            return self.amplitude if self.onset <= start < self.offset else 0.0
        active = max(0.0, min(end, self.offset) - max(start, self.onset))
        return self.amplitude * active / (end - start)

    def parameters(self) -> dict[str, Any]:
        return {
            "regime": self.regime,
            "amplitude_W_m2": self.amplitude,
            "onset_s": self.onset,
            "duration_s": self.duration,
            "energy_J_m2": self.amplitude * self.duration,
        }

    def describe(self) -> str:
        return (
            f"pulse of {self.amplitude:+g} W m^-2 for {self.duration:g} s "
            f"from t = {self.onset:g} s ({self.amplitude * self.duration:+.4g} J m^-2)"
        )


@dataclass(frozen=True)
class ImpulseForcing(Forcing):
    """A fixed quantity of energy delivered in whichever step contains ``onset``.

    The energy is exact and independent of the time step; the *shape* is not
    resolved below one step and this forcing says so rather than implying a
    delta function it cannot represent. An impulse response measured from it is
    therefore reproducible across time steps, while any structure in the response
    faster than one step is absent by construction.
    """

    energy: float
    onset: float = 0.0
    regime: str = "impulse"

    def __post_init__(self) -> None:
        if self.onset < 0.0:
            raise ConfigurationError("impulse onset must not precede the start of the run")

    def mean_amplitude(self, start: float, end: float) -> float:
        self._check_interval(start, end)
        if end == start:
            # An impulse has no instantaneous value: reporting one would invent a
            # flux whose magnitude depends on a step this call does not know.
            return 0.0
        return self.energy / (end - start) if start <= self.onset < end else 0.0

    def parameters(self) -> dict[str, Any]:
        return {
            "regime": self.regime,
            "energy_J_m2": self.energy,
            "onset_s": self.onset,
            "shape": "one master step wide; not resolved below the step",
        }

    def describe(self) -> str:
        return (
            f"impulse of {self.energy:+.4g} J m^-2 at t = {self.onset:g} s, "
            f"spread over the one step containing it"
        )


@dataclass(frozen=True)
class PeriodicForcing(Forcing):
    """``A sin(2 pi (t - onset) / period + phase)``, for gain and phase analysis.

    Charter 11.3 asks for frequency-dependent gain and phase; this is the input
    side of that measurement. The mean over a whole number of periods is exactly
    zero, so a periodic experiment cannot warm the system by accident -- and the
    exact interval mean below is what makes that true at any time step rather
    than only in the limit.
    """

    amplitude: float
    period: float
    onset: float = 0.0
    phase: float = 0.0
    regime: str = "periodic"

    def __post_init__(self) -> None:
        if self.period <= 0.0:
            raise ConfigurationError("periodic forcing must declare a positive period")
        if self.onset < 0.0:
            raise ConfigurationError("periodic onset must not precede the start of the run")

    @property
    def angular_frequency(self) -> float:
        return 2.0 * math.pi / self.period

    def _angle(self, time: float) -> float:
        return self.angular_frequency * (time - self.onset) + self.phase

    def mean_amplitude(self, start: float, end: float) -> float:
        self._check_interval(start, end)
        begin = max(start, self.onset)
        if end <= self.onset:
            return 0.0
        if end == start:
            return self.amplitude * math.sin(self._angle(start))
        # Exact integral of A sin(w t + p): the mean over the step, not its
        # midpoint value. At a step of a tenth of a period the two differ by
        # about 0.4 percent, which is larger than the tolerances the gain
        # measurement below is judged against.
        omega = self.angular_frequency
        integral = (self.amplitude / omega) * (
            math.cos(self._angle(begin)) - math.cos(self._angle(end))
        )
        return integral / (end - start)

    def parameters(self) -> dict[str, Any]:
        return {
            "regime": self.regime,
            "amplitude_W_m2": self.amplitude,
            "period_s": self.period,
            "onset_s": self.onset,
            "phase_rad": self.phase,
        }

    def describe(self) -> str:
        return (
            f"periodic forcing of {self.amplitude:g} W m^-2 amplitude at a period of "
            f"{self.period:g} s from t = {self.onset:g} s"
        )
