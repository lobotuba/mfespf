"""How far the column is from the observed surface temperature, in kelvin.

The number everyone asks for first, which is why it needs the most guarding. It
has been quoted wrongly in this project in three distinct ways, and this module
exists to make each of them impossible rather than merely discouraged:

**From a run that had not equilibrated.** The dry radiative-convective column at
20000 steps reads -16.3 K, and it is not 16 K too cold -- it is unfinished,
carrying +43.8 W m^-2 at the top of the atmosphere and still warming. A deep
ocean under the same column looked 3 K cooler at eleven years for the same
reason. So ``surface_temperature_bias`` reports **no error at all** unless the
run is in balance, and says what the imbalance was instead.

**From an endpoint rather than a mean.** At a physically consistent exchange
coefficient the column limit-cycles with a surface span of about 1.75 K that does
*not* shrink with the timestep. A single step's value is a point on a cycle.
Everything here is a time-mean over the tail of the run.

**Without its process quantities.** The column has sat 3 K from the observed
temperature while its evaporation was twelve times too small and its Bowen ratio
thirty-seven times too large -- two large errors of opposite sign cancelling at
the one number anyone looks at. The report therefore carries the water cycle
beside the temperature, so a good temperature cannot be read on its own.

The reference value is 288.0 K, which is a global-mean surface air temperature
and not a measurement of anywhere in particular. It is stated to 0.1 K because
that is the precision the comparison supports; the modern observed value is
quoted between 287.6 and 288.2 K depending on the record and the period, so a
bias smaller than about 0.5 K is not a distinguishable result.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

__all__ = [
    "OBSERVED_SURFACE_TEMPERATURE_K",
    "OBSERVED_TEMPERATURE_SPREAD_K",
    "BiasReport",
    "surface_temperature_bias",
]

#: Global-mean surface air temperature, late twentieth century.
OBSERVED_SURFACE_TEMPERATURE_K = 288.0

#: How much the observational records disagree among themselves. A modelled bias
#: smaller than this is not a result, and reporting one to three decimals without
#: saying so would be false precision.
OBSERVED_TEMPERATURE_SPREAD_K = 0.5

#: A run is treated as equilibrated when the mean top-of-atmosphere imbalance
#: over the window is smaller than this. Chosen against what the imbalance *is*
#: rather than what would be convenient: an equilibrated column here sits at
#: 0.001 W m^-2 and an unfinished one at 43.8, so anything between those two
#: separates them, and 0.1 leaves three orders of magnitude of margin on the
#: side that matters.
EQUILIBRIUM_IMBALANCE_W_M2 = 0.1


@dataclass(frozen=True, slots=True)
class BiasReport:
    """The error in kelvin, or the reason there is no error to report."""

    label: str
    n_steps: int
    window: str
    surface_temperature_k: float
    surface_span_k: float
    imbalance_w_m2: float
    imbalance_span_w_m2: float
    evaporation_mm_day: float | None = None
    precipitable_water_kg_m2: float | None = None
    bowen_ratio: float | None = None

    @property
    def equilibrated(self) -> bool:
        return abs(self.imbalance_w_m2) < EQUILIBRIUM_IMBALANCE_W_M2

    @property
    def error_k(self) -> float | None:
        """Kelvin from the observed value, or ``None`` when the run has not settled.

        ``None`` rather than a number, because a bias from an unequilibrated run
        is not a small result or an uncertain one -- it is a different quantity,
        and returning it would let it be compared against a real one.
        """
        if not self.equilibrated:
            return None
        return self.surface_temperature_k - OBSERVED_SURFACE_TEMPERATURE_K

    @property
    def distinguishable(self) -> bool:
        """Whether the error is larger than the records disagree among themselves."""
        error = self.error_k
        return error is not None and abs(error) > OBSERVED_TEMPERATURE_SPREAD_K

    def describe(self) -> str:
        if not self.equilibrated:
            direction = "warming" if self.imbalance_w_m2 > 0 else "cooling"
            return (
                f"{self.label}: no error to report -- the run carries "
                f"{self.imbalance_w_m2:+.2f} W m^-2 and is still {direction} "
                f"(T_s {self.surface_temperature_k:.3f} K at {self.n_steps} steps)"
            )
        error = self.error_k
        assert error is not None
        verdict = "" if self.distinguishable else "  [inside the observational spread]"
        line = (
            f"{self.label}: {error:+.3f} K  "
            f"(T_s {self.surface_temperature_k:.3f} K, span {self.surface_span_k:.2f}, "
            f"imbalance {self.imbalance_w_m2:+.3f} W m^-2){verdict}"
        )
        process = []
        if self.evaporation_mm_day is not None:
            process.append(f"evaporation {self.evaporation_mm_day:.2f} mm/day (obs 2.7)")
        if self.precipitable_water_kg_m2 is not None:
            process.append(f"PW {self.precipitable_water_kg_m2:.1f} kg/m2 (obs 25)")
        if self.bowen_ratio is not None:
            process.append(f"Bowen {self.bowen_ratio:.3f} (obs 0.26)")
        if process:
            line += "\n    " + ", ".join(process)
        return line


def _series(result: Any, name: str) -> np.ndarray | None:
    steps = result.steps
    if not steps or name not in steps[-1].outputs:
        return None
    return np.array(
        [
            float(np.asarray(s.outputs[name].values).sum()) if name in s.outputs else np.nan
            for s in steps
        ]
    )


def surface_temperature_bias(result: Any, *, label: str, tail_fraction: float = 0.2) -> BiasReport:
    """Reduce a completed run to its error in kelvin, or to why it has none.

    ``tail_fraction`` is the part of the run averaged over. A fifth is long
    enough to cover several periods of the convective limit cycle at every
    timestep tried, and short enough to exclude the spin-up.
    """
    temperature = _series(result, "surface_temperature")
    if temperature is None:
        raise ValueError(f"{label}: the run published no surface temperature")
    tail = slice(int((1.0 - tail_fraction) * len(temperature)), None)

    absorbed = _series(result, "toa_downward_radiation")
    reflected = _series(result, "reflected_shortwave_radiation")
    outgoing = _series(result, "outgoing_longwave_radiation")
    if absorbed is None or reflected is None or outgoing is None:
        imbalance = np.zeros_like(temperature)
    else:
        imbalance = absorbed - reflected - outgoing

    def tail_mean(name: str, scale: float = 1.0) -> float | None:
        series = _series(result, name)
        return None if series is None else float(np.nanmean(series[tail])) * scale

    latent = _series(result, "surface_latent_heat_flux")
    sensible = _series(result, "surface_sensible_heat_flux")
    # Where convection reaches the surface, the ground also heats the air by
    # mixing, and that is sensible heat as much as the bulk flux is. Leaving it
    # out read the moist-convective column's Bowen ratio as 0.01: the bulk flux
    # is small there because the mixing has already carried the heat up.
    convective = _series(result, "surface_convective_heat_flux")
    if sensible is not None and convective is not None:
        sensible = sensible + convective
    bowen = None
    if latent is not None and sensible is not None:
        ratio = sensible[tail] / np.where(latent[tail] == 0.0, np.nan, latent[tail])
        if np.any(np.isfinite(ratio)):
            bowen = float(np.nanmean(ratio))

    return BiasReport(
        label=label,
        n_steps=len(temperature),
        window=f"last {tail_fraction:.0%} of the run",
        surface_temperature_k=float(np.mean(temperature[tail])),
        surface_span_k=float(np.max(temperature[tail]) - np.min(temperature[tail])),
        imbalance_w_m2=float(np.mean(imbalance[tail])),
        imbalance_span_w_m2=float(np.max(imbalance[tail]) - np.min(imbalance[tail])),
        evaporation_mm_day=tail_mean("surface_evaporation_flux", 86400.0),
        precipitable_water_kg_m2=tail_mean("atmosphere_precipitable_water"),
        bowen_ratio=bowen,
    )
