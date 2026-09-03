"""Measuring characteristic timescales from a run.

Charter 11.1 and 11.3, requirement TIM-008.

TIM-008 is a prohibition as much as a capability: timescales are to be measured
*per mechanism, subsystem, observable and forcing regime* rather than imposed as
one universal climate index. A function returning a single number would satisfy
the letter and destroy the point, so every measurement here carries those four
labels, and ``TimescaleReport.single_index`` exists only to refuse when the
measured timescales disagree.

Two ways of getting a number are offered and they are not interchangeable:

``efolding_time`` measures. It reports when the response first covered 1 - 1/e of
its total change, assuming nothing about the shape of the curve. For a single
exponential it recovers the timescale exactly; for anything else it reports a
property of the trajectory that is well defined but is not "the" timescale.

``fit_exponential_modes`` fits. It assumes the response is a sum of decaying
exponentials and returns the modes, with the residual that says how well that
assumption held. A fit whose residual is not small is not evidence of two
timescales -- it is evidence that the model used to extract them was wrong.

Both are numpy-only and deterministic: the same series gives the same answer on
every platform, which CMP-009 requires of anything a claim rests on.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from dataclasses import field as dc_field

import numpy as np

from mfespf.core.errors import ConfigurationError, InvalidForScientificInferenceError

__all__ = [
    "ExponentialFit",
    "PeriodicResponse",
    "TimescaleMeasurement",
    "TimescaleReport",
    "efolding_time",
    "fit_exponential_modes",
    "lag_correlation",
    "periodic_response",
]


# --- the labelled measurement ----------------------------------------------


@dataclass(frozen=True)
class TimescaleMeasurement:
    """One timescale, with the four labels TIM-008 requires it to carry."""

    observable: str
    subsystem: str
    mechanism: str
    forcing_regime: str
    timescale_s: float
    method: str
    amplitude: float = 0.0
    relative_residual: float = 0.0
    constrained: bool = True
    """False when the record is shorter than the timescale it reports.

    A fit will happily return a timescale several times longer than the data it
    saw, and the number looks exactly like a measured one. It is an
    extrapolation under the assumed model, and the charter's separation of
    operating modes (TIM-006) applies inside a single analysis as much as
    between runs: it must not be reported in the language of a measurement."""
    notes: str = ""

    @property
    def timescale_days(self) -> float:
        return self.timescale_s / 86400.0

    @property
    def timescale_years(self) -> float:
        return self.timescale_s / (365.25 * 86400.0)

    def describe(self) -> str:
        qualifier = (
            "" if self.constrained else " [EXTRAPOLATED: longer than the record it was fitted to]"
        )
        return (
            f"{self.observable} / {self.subsystem} / {self.mechanism} "
            f"under {self.forcing_regime} forcing: tau = {self.timescale_s:.6g} s "
            f"({self.timescale_years:.4g} yr) by {self.method}{qualifier}"
        )


@dataclass
class TimescaleReport:
    """A set of measurements that is deliberately not reducible to one number."""

    measurements: list[TimescaleMeasurement] = dc_field(default_factory=list)

    def add(self, measurement: TimescaleMeasurement) -> None:
        self.measurements.append(measurement)

    def select(
        self,
        *,
        observable: str | None = None,
        subsystem: str | None = None,
        mechanism: str | None = None,
        forcing_regime: str | None = None,
    ) -> list[TimescaleMeasurement]:
        def matches(m: TimescaleMeasurement) -> bool:
            return (
                (observable is None or m.observable == observable)
                and (subsystem is None or m.subsystem == subsystem)
                and (mechanism is None or m.mechanism == mechanism)
                and (forcing_regime is None or m.forcing_regime == forcing_regime)
            )

        return [m for m in self.measurements if matches(m)]

    @property
    def spread(self) -> float:
        """Ratio of the longest measured timescale to the shortest."""
        if not self.measurements:
            return 1.0
        values = [m.timescale_s for m in self.measurements if m.timescale_s > 0.0]
        if not values:
            return 1.0
        return max(values) / min(values)

    def single_index(self, *, tolerated_spread: float = 2.0) -> float:
        """The one number TIM-008 forbids, available only when it is honest.

        If every measured timescale agrees within ``tolerated_spread`` there is
        no harm in quoting one; if they do not, quoting one would assert
        something the measurements contradict, and this raises instead of
        averaging them into a number that describes no mechanism.
        """
        if not self.measurements:
            raise ConfigurationError("no timescales have been measured")
        if self.spread > tolerated_spread:
            detail = "; ".join(sorted(m.describe() for m in self.measurements))
            raise InvalidForScientificInferenceError(
                f"the measured timescales span a factor of {self.spread:.4g}, which exceeds "
                f"the tolerated {tolerated_spread:g}; TIM-008 forbids collapsing them into a "
                f"single climate index. Report them per mechanism instead: {detail}"
            )
        return float(np.mean([m.timescale_s for m in self.measurements]))

    def describe(self) -> str:
        lines = [f"{len(self.measurements)} measured timescales, spread {self.spread:.4g}x"]
        lines.extend(f"  {m.describe()}" for m in self.measurements)
        return "\n".join(lines)


# --- measurement without a model -------------------------------------------


def efolding_time(times: np.ndarray, values: np.ndarray, *, final: float | None = None) -> float:
    """Time at which the response first covers ``1 - 1/e`` of its total change.

    No functional form is assumed. ``final`` names the asymptote when the run did
    not reach it; without it the last sample is used, which understates the
    timescale of a trajectory still moving and is therefore reported as measured
    rather than as extrapolated.
    """
    t = np.asarray(times, dtype=float)
    y = np.asarray(values, dtype=float)
    if t.ndim != 1 or y.shape != t.shape:
        raise ConfigurationError("times and values must be one-dimensional and the same length")
    if t.size < 3:
        raise ConfigurationError("an e-folding time needs at least three samples")

    start = float(y[0])
    end = float(y[-1]) if final is None else float(final)
    total = end - start
    if total == 0.0:
        raise ConfigurationError(
            "the response does not change, so it has no e-folding time; a null control "
            "should be checked for being zero, not fitted"
        )

    target = start + total * (1.0 - 1.0 / math.e)
    progress = (y - start) / total
    threshold = 1.0 - 1.0 / math.e
    crossed = np.nonzero(progress >= threshold)[0]
    if crossed.size == 0:
        raise ConfigurationError(
            f"the response never reaches {threshold:.4f} of its change within the record; "
            f"run longer or supply the asymptote through `final`"
        )
    index = int(crossed[0])
    if index == 0:
        return float(t[0])
    # Linear interpolation between the bracketing samples: the crossing is
    # generally inside a step, and rounding it to a step would quantise every
    # measured timescale by the output interval.
    y0, y1 = float(y[index - 1]), float(y[index])
    if y1 == y0:
        return float(t[index])
    fraction = (target - y0) / (y1 - y0)
    return float(t[index - 1] + fraction * (t[index] - t[index - 1]))


# --- measurement under an assumed model ------------------------------------


@dataclass(frozen=True)
class ExponentialFit:
    """``y(t) = asymptote + sum_i amplitude_i exp(-t / tau_i)``."""

    asymptote: float
    amplitudes: tuple[float, ...]
    timescales: tuple[float, ...]
    rms_residual: float
    relative_residual: float
    n_modes: int
    record_span: float = 0.0

    @property
    def constrained(self) -> tuple[bool, ...]:
        """Which modes the record is long enough to have measured.

        A mode longer than the record was inferred from the first part of a
        curve the run never finished. The fit is not wrong to report it -- it is
        the best estimate under the assumed model -- but calling it measured
        would be, so the distinction is carried rather than left to whoever
        remembers how long the run was.
        """
        return tuple(tau <= self.record_span for tau in self.timescales)

    def evaluate(self, times: np.ndarray) -> np.ndarray:
        t = np.asarray(times, dtype=float)
        out = np.full(t.shape, self.asymptote)
        for amplitude, tau in zip(self.amplitudes, self.timescales, strict=True):
            out = out + amplitude * np.exp(-t / tau)
        return out

    def describe(self) -> str:
        modes = ", ".join(
            f"{a:+.5g} exp(-t/{tau:.5g} s{'' if ok else ', extrapolated'})"
            for a, tau, ok in zip(self.amplitudes, self.timescales, self.constrained, strict=True)
        )
        return f"{self.asymptote:.6g} + {modes}  (relative residual {self.relative_residual:.3g})"


def _linear_amplitudes(t: np.ndarray, y: np.ndarray, taus: np.ndarray) -> tuple[np.ndarray, float]:
    """Least-squares asymptote and amplitudes for fixed timescales."""
    design = np.empty((t.size, taus.size + 1))
    design[:, 0] = 1.0
    for column, tau in enumerate(taus, start=1):
        design[:, column] = np.exp(-t / tau)
    coefficients, *_ = np.linalg.lstsq(design, y, rcond=None)
    residual = float(np.sqrt(np.mean((design @ coefficients - y) ** 2)))
    return coefficients, residual


def fit_exponential_modes(
    times: np.ndarray,
    values: np.ndarray,
    *,
    n_modes: int = 1,
    tau_bounds: tuple[float, float] | None = None,
    max_samples: int = 400,
) -> ExponentialFit:
    """Fit a sum of decaying exponentials by search over timescales.

    For fixed timescales the model is linear in its amplitudes, so the fit is a
    grid search over log-spaced timescales with an exact least-squares solve
    inside it, refined twice around the best cell. That is slower than a gradient
    method and does not need one: it is deterministic, has no starting guess to
    tune, and cannot converge to a different answer on a different platform.

    The returned ``relative_residual`` is the quantity that decides whether the
    fit means anything. Comparing a one-mode and a two-mode fit of the same
    series through it is the honest way to ask whether a second timescale is
    present, rather than reading a second timescale off a fit that was granted
    one whether or not the data supported it.
    """
    t = np.asarray(times, dtype=float)
    y = np.asarray(values, dtype=float)
    if t.ndim != 1 or y.shape != t.shape:
        raise ConfigurationError("times and values must be one-dimensional and the same length")
    if n_modes < 1:
        raise ConfigurationError("a fit needs at least one mode")
    if t.size < 2 * n_modes + 2:
        raise ConfigurationError(
            f"a {n_modes}-mode fit needs at least {2 * n_modes + 2} samples, got {t.size}"
        )

    # Time is measured from the start of the record, so the amplitudes are the
    # mode amplitudes at t = 0 rather than at an arbitrary epoch.
    t = t - t[0]
    if t.size > max_samples:
        keep = np.unique(np.linspace(0, t.size - 1, max_samples).astype(int))
        t, y = t[keep], y[keep]

    span = float(t[-1])
    if span <= 0.0:
        raise ConfigurationError("the record has no duration")
    lower, upper = tau_bounds if tau_bounds is not None else (span / 500.0, span * 5.0)
    if not 0.0 < lower < upper:
        raise ConfigurationError(f"tau bounds {(lower, upper)} must satisfy 0 < lower < upper")

    scale = float(np.std(y))
    if scale == 0.0:
        raise ConfigurationError("the series is constant, so no timescale can be extracted from it")

    def search(bounds: list[tuple[float, float]], points: int) -> np.ndarray:
        grids = [np.logspace(math.log10(low), math.log10(high), points) for low, high in bounds]
        best_taus = np.array([])
        best_residual = math.inf
        for candidate in _increasing_products(grids):
            taus = np.asarray(candidate)
            # Modes closer than five percent are not separable by this data, and
            # admitting them lets the search report a spurious pair of nearly
            # equal timescales as though it had resolved two.
            if taus.size > 1 and np.any(np.diff(taus) < taus[:-1] * 0.05):
                continue
            _, residual = _linear_amplitudes(t, y, taus)
            if residual < best_residual:
                best_residual, best_taus = residual, taus
        if best_taus.size == 0:
            raise ConfigurationError("no separable set of timescales was found in the bounds")
        return best_taus

    # Coarse pass over the whole range, then refinements whose window is a few
    # cells of the previous grid. A window that does not shrink -- a fixed factor
    # around the current best -- leaves the answer quantised at the coarse grid's
    # resolution, which is a systematic error large enough to look like physics.
    points = 40 if n_modes == 1 else 26
    bounds = [(lower, upper)] * n_modes
    taus = search(bounds, points)
    for _ in range(3):
        new_bounds = []
        for index, tau in enumerate(taus):
            low, high = bounds[index]
            step = (high / low) ** (1.0 / (points - 1))
            width = step**2
            new_bounds.append((max(lower, tau / width), min(upper, tau * width)))
        bounds = new_bounds
        taus = search(bounds, points)

    taus = np.sort(taus)
    coefficients, residual = _linear_amplitudes(t, y, taus)
    return ExponentialFit(
        asymptote=float(coefficients[0]),
        amplitudes=tuple(float(c) for c in coefficients[1:]),
        timescales=tuple(float(v) for v in taus),
        rms_residual=residual,
        relative_residual=residual / scale,
        n_modes=n_modes,
        record_span=span,
    )


def _increasing_products(grids: list[np.ndarray]) -> list[tuple[float, ...]]:
    """Every strictly increasing selection, one value from each grid.

    Strictly increasing because a permutation of the same timescales is the same
    model: searching both orderings would double the cost and return one of two
    identical answers depending on floating-point ties.
    """
    if len(grids) == 1:
        return [(float(v),) for v in grids[0]]
    out: list[tuple[float, ...]] = []
    for value in grids[0]:
        for tail in _increasing_products(grids[1:]):
            if tail[0] > value:
                out.append((float(value), *tail))
    return out


# --- frequency response -----------------------------------------------------


@dataclass(frozen=True)
class PeriodicResponse:
    """Gain and phase of a response to a periodic forcing (charter 11.3)."""

    period: float
    forcing_amplitude: float
    response_amplitude: float
    phase_lag_rad: float

    @property
    def gain(self) -> float:
        """Response amplitude per unit forcing amplitude, in response units."""
        return self.response_amplitude / self.forcing_amplitude

    @property
    def phase_lag_s(self) -> float:
        return self.phase_lag_rad * self.period / (2.0 * math.pi)

    def describe(self) -> str:
        return (
            f"period {self.period:.6g} s: gain {self.gain:.6g} per W m^-2, "
            f"lag {self.phase_lag_rad:.4g} rad ({self.phase_lag_s:.6g} s)"
        )


def _single_frequency_amplitude_phase(
    t: np.ndarray, y: np.ndarray, omega: float
) -> tuple[float, float]:
    """Project onto sin and cos at one frequency: a one-bin discrete transform."""
    centred = y - float(np.mean(y))
    sine = np.sin(omega * t)
    cosine = np.cos(omega * t)
    # Least squares against the two basis functions rather than the textbook
    # 2/N sums, which assume the record is a whole number of periods sampled
    # uniformly. It usually is, but a fit does not need it to be.
    design = np.column_stack([sine, cosine])
    coefficients, *_ = np.linalg.lstsq(design, centred, rcond=None)
    a, b = float(coefficients[0]), float(coefficients[1])
    return math.hypot(a, b), math.atan2(b, a)


def periodic_response(
    times: np.ndarray,
    forcing: np.ndarray,
    response: np.ndarray,
    *,
    period: float,
    discard_fraction: float = 0.5,
) -> PeriodicResponse:
    """Gain and phase lag of ``response`` relative to ``forcing`` at ``period``.

    The first ``discard_fraction`` of the record is dropped: a periodic
    experiment starts from a transient that has nothing to do with the periodic
    response, and including it biases both the gain and the phase. Half is the
    default because a run of several periods has usually shed the transient by
    then -- but the fraction is a declared parameter, not a hidden choice, and a
    run whose transient is longer than half the record must say so.
    """
    t = np.asarray(times, dtype=float)
    f = np.asarray(forcing, dtype=float)
    r = np.asarray(response, dtype=float)
    if not (t.shape == f.shape == r.shape) or t.ndim != 1:
        raise ConfigurationError("times, forcing and response must be the same 1-D length")
    if period <= 0.0:
        raise ConfigurationError("period must be positive")
    if not 0.0 <= discard_fraction < 1.0:
        raise ConfigurationError("discard_fraction must lie in [0, 1)")

    start = int(t.size * discard_fraction)
    t, f, r = t[start:], f[start:], r[start:]
    if t.size < 8:
        raise ConfigurationError("too few samples remain after discarding the transient")
    if t[-1] - t[0] < period:
        raise ConfigurationError(
            f"the analysed record spans {t[-1] - t[0]:g} s, less than one {period:g} s period"
        )

    omega = 2.0 * math.pi / period
    forcing_amplitude, forcing_phase = _single_frequency_amplitude_phase(t, f, omega)
    response_amplitude, response_phase = _single_frequency_amplitude_phase(t, r, omega)
    if forcing_amplitude == 0.0:
        raise ConfigurationError("the forcing has no component at the analysed period")

    lag = forcing_phase - response_phase
    # Wrapped into [0, 2 pi): a lag is a delay, and a negative one here would
    # mean the response preceded the forcing, which a causal system cannot do.
    lag = lag % (2.0 * math.pi)
    return PeriodicResponse(
        period=period,
        forcing_amplitude=forcing_amplitude,
        response_amplitude=response_amplitude,
        phase_lag_rad=lag,
    )


def lag_correlation(
    times: np.ndarray, driver: np.ndarray, response: np.ndarray, *, max_lag_samples: int
) -> tuple[float, float]:
    """Lag of peak correlation between two series, and that correlation.

    Charter 11.3 lists lag and cross-correlation analysis alongside the transfer
    function; this is the coarse version, useful where the response is not
    sinusoidal. Returns ``(lag_in_seconds, correlation)`` with a positive lag
    meaning the response follows the driver.
    """
    t = np.asarray(times, dtype=float)
    a = np.asarray(driver, dtype=float)
    b = np.asarray(response, dtype=float)
    if not (t.shape == a.shape == b.shape) or t.ndim != 1:
        raise ConfigurationError("times, driver and response must be the same 1-D length")
    if not 0 < max_lag_samples < t.size:
        raise ConfigurationError("max_lag_samples must be between 1 and the record length")

    a = a - float(np.mean(a))
    b = b - float(np.mean(b))
    if np.all(a == 0.0) or np.all(b == 0.0):
        raise ConfigurationError("a constant series has no correlation lag")

    best_lag, best_correlation = 0, -math.inf
    for lag in range(max_lag_samples + 1):
        left = a[: a.size - lag] if lag else a
        right = b[lag:] if lag else b
        denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
        if denominator == 0.0:
            continue
        correlation = float(np.dot(left, right) / denominator)
        if correlation > best_correlation:
            best_lag, best_correlation = lag, correlation

    interval = float(np.mean(np.diff(t))) if t.size > 1 else 0.0
    return (best_lag * interval, best_correlation)
