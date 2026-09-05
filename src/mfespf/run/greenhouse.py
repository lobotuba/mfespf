"""Stage 6 configuration: a radiating column under a recorded Sun and a recorded CO2.

Charter section 19 and 9.5; requirements CHM-001, HYP-005, INF-001, TIM-005.

The first configuration in this framework in which two declared records drive the
same run, and the first in which a greenhouse gas acts. Carbon dioxide enters
through the absorption coefficient, so the radiative transfer computes its effect
rather than being handed one.

**The base state is wrong, and knowing why matters.** This column equilibrates
near 315 K rather than 288 K. The cause is not the CO2 concentration but the
band's width: a reduced column with eight layers and no stratosphere has a
weaker contrast between the surface and the level it radiates from than the real
atmosphere does, so reproducing the accepted forcing per doubling requires a band
that is opaque over a wider range of wavelengths than the real 15 micron band.
The forcing is right and the width is not, and the surplus width closes part of
the atmospheric window.

Two things follow. The comparison below is of *anomalies*, differenced against
the model's own baseline, so a constant offset does not enter it. But the Planck
response at 315 K is about a third stronger than at 288 K, so the sensitivity is
suppressed relative to what this configuration would give at the right base
state -- and that is a declared limitation of every claim made with it, not a
detail.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from mfespf.data.adapter import HistoricalDataAdapter, RecordedSolarSource
from mfespf.data.loaders import load_declared
from mfespf.data.records import TimeSeries, splice
from mfespf.physics.greenhouse import PrescribedGreenhouseModule, greenhouse_grid
from mfespf.registries.datasets import DatasetRegistry, DatasetRole
from mfespf.run.column_configurations import radiating_column, window_absorption
from mfespf.run.driver import Kernel
from mfespf.run.historical import replay_clock
from mfespf.run.modes import OperatingMode

__all__ = [
    "CO2_JOIN_YEAR",
    "GREENHOUSE_BACKGROUND",
    "co2_record",
    "greenhouse_replay",
    "spin_up",
]

DAY = 86400.0
IRRADIANCE = "NRLTSI2-MONTHLY"
CO2_JOIN_YEAR = 1970
"""Where the ice-core record hands over to the instrumental one.

Inside the twenty-year overlap, and late enough that the instrumental record is
well established. The step between the two at the join is measured and reported
rather than smoothed away.
"""

GREENHOUSE_BACKGROUND = 3.0e-5
"""Infrared absorption of everything that is not CO2, m^2 per kg of air.

Set from the transparency of the atmospheric window rather than from any
temperature: over a column of about ten tonnes per square metre this gives an
optical depth near 0.3 in the window, which is the order the water-vapour
continuum produces. It is a free parameter of the reduced model, and it is held
fixed so that a difference between two runs is attributable to CO2 alone.
"""


def co2_record(join_year: int = CO2_JOIN_YEAR) -> tuple[TimeSeries, Any]:
    """The ice-core and instrumental records joined, with the join reported."""
    from datetime import UTC, datetime

    from mfespf.data.loaders import seconds_since_epoch

    ice = next(iter(load_declared("LAWDOME-CO2-ANNUAL")._series.values()))
    station = next(iter(load_declared("MLO-CO2-ANNUAL")._series.values()))
    return splice(ice, station, at=seconds_since_epoch(datetime(join_year, 1, 1, tzinfo=UTC)))


def spin_up(
    *, master_dt: float, n_layers: int, n_band: int, ppm: float, irradiance: float, steps: int
) -> tuple[float, np.ndarray]:
    """Equilibrate the column at fixed drivers and return the state it reaches.

    Every configuration in this framework that starts away from equilibrium
    spends its first years relaxing towards it, and that relaxation appears in
    the output as a trend. It caught the solar-only experiment once already --
    a constant-irradiance control with half a kelvin of range -- and it caught
    this one with seventeen.

    So the replay starts from the state the drivers imply rather than from a
    guess, and the constant-driver control then produces a flat line, which is
    the only thing that makes a difference between two runs attributable to the
    driver that changed.
    """
    from mfespf.physics.solar import ConstantSpectrumSolarSource

    grid = greenhouse_grid(n_background=24, n_band=n_band, band_lower_m=6.5e-6, band_upper_m=5.0e-5)
    kernel = radiating_column(
        run_id="GREENHOUSE-SPINUP",
        purpose="equilibrate the column at fixed drivers before the replay begins",
        absorption_coefficient=window_absorption(grid, kappa_opaque=GREENHOUSE_BACKGROUND),
        grid=grid,
        n_layers=n_layers,
        areal_heat_capacity=1.0e7,  # a thin slab: the equilibrium is the same, reached sooner
        master_dt=master_dt,
        convection=True,
        greenhouse=PrescribedGreenhouseModule(concentration=ppm),
        solar_source=ConstantSpectrumSolarSource(grid=grid, total_irradiance=irradiance),
        initial_surface_temperature=300.0,
        initial_air_temperature=270.0,
    )
    radiation = next(s.module for s in kernel.schedule.ordered() if s.module_id == "radiation")
    radiation.greenhouse_background = GREENHOUSE_BACKGROUND  # type: ignore[attr-defined]
    result = kernel.run(steps)
    outputs = result.steps[-1].outputs
    return (
        float(outputs["surface_temperature"].values),
        np.asarray(outputs["atmosphere_layer_temperature"].values, dtype=float).copy(),
    )


def greenhouse_replay(
    *,
    master_dt: float = 2.0 * DAY,
    n_layers: int = 8,
    areal_heat_capacity: float = 2.04e8,
    hold_co2: bool = False,
    hold_solar: bool = False,
    n_band: int = 48,
    spin_up_steps: int = 1500,
    **kwargs: Any,
) -> Kernel:
    """The column driven by the irradiance and CO2 records together.

    ``hold_co2`` and ``hold_solar`` freeze one driver at its record mean, which
    is how the two mechanisms are separated: three runs -- both, solar only, CO2
    only -- differ in exactly one thing each, so a difference between any pair is
    attributable to the driver that changed (HYP-007).
    """
    grid = greenhouse_grid(n_background=24, n_band=n_band, band_lower_m=6.5e-6, band_upper_m=5.0e-5)
    solar_adapter = load_declared(IRRADIANCE)
    clock = replay_clock(solar_adapter, "total_solar_irradiance", master_dt)
    if hold_solar:
        # Held at the record's *first* value, which is also what the spin-up
        # equilibrated to. Holding at the record mean instead would start the
        # control away from its own equilibrium, and the resulting transient
        # would belong to neither the control nor the experiment while showing
        # up in the difference between them.
        series = solar_adapter.series("total_solar_irradiance")
        first = float(series.value_at(clock[0]))
        flat = series.with_values(
            series.times,
            np.full(series.times.size, first),
            note="constant irradiance at the record's first value; the solar control",
        )
        solar_adapter = HistoricalDataAdapter(solar_adapter.record, {flat.variable: flat})
    solar_adapter.prepare("total_solar_irradiance", clock, target_interval=master_dt)

    co2, _ = co2_record()
    initial_co2 = float(co2.value_at(clock[0]))

    def concentration(time: float) -> float:
        return initial_co2 if hold_co2 else co2.value_at(time)

    datasets = DatasetRegistry()
    run_id = str(kwargs.setdefault("run_id", "STAGE6-GREENHOUSE-001"))
    solar_adapter.register(
        datasets,
        role=DatasetRole.BOUNDARY_CONDITION,
        claim_id=str(kwargs.pop("claim_id", "CLM-GREENHOUSE-001")),
        run_id=run_id,
        influenced=("solar",),
        note="prescribed top-of-atmosphere irradiance over the replay period",
    )

    kwargs.setdefault(
        "purpose",
        "Stage 6 greenhouse replay: a radiative-convective column driven by a recorded "
        "irradiance history and a recorded carbon dioxide history, with the greenhouse "
        "effect computed by the radiative transfer rather than imposed as a forcing.",
    )
    kwargs.setdefault("convection", True)
    if "initial_surface_temperature" not in kwargs:
        surface, layers = spin_up(
            master_dt=master_dt,
            n_layers=n_layers,
            n_band=n_band,
            ppm=float(co2.value_at(clock[0])),
            irradiance=float(solar_adapter.prepared("total_solar_irradiance").values[0]),
            steps=spin_up_steps,
        )
        kwargs["initial_surface_temperature"] = surface
        kwargs["initial_air_temperature"] = layers

    kernel = radiating_column(
        absorption_coefficient=window_absorption(grid, kappa_opaque=GREENHOUSE_BACKGROUND),
        grid=grid,
        n_layers=n_layers,
        areal_heat_capacity=areal_heat_capacity,
        master_dt=master_dt,
        start_time=float(clock[0]),
        operating_mode=OperatingMode.PRESCRIBED_HISTORY,
        datasets=datasets,
        greenhouse=PrescribedGreenhouseModule(concentration=concentration),
        solar_source=RecordedSolarSource(
            grid=grid, adapter=solar_adapter, variable="total_solar_irradiance"
        ),
        **kwargs,
    )
    radiation = next(s.module for s in kernel.schedule.ordered() if s.module_id == "radiation")
    radiation.greenhouse_background = GREENHOUSE_BACKGROUND  # type: ignore[attr-defined]
    kernel.adapter = solar_adapter
    kernel.manifest = None
    return kernel
