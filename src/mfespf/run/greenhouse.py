"""Stage 6 configuration: a radiating column under a recorded Sun and a recorded CO2.

Charter section 19 and 9.5; requirements CHM-001, HYP-005, INF-001, TIM-005.

The first configuration in this framework in which two declared records drive the
same run, and the first in which a greenhouse gas acts. Carbon dioxide enters
through the absorption coefficient, so the radiative transfer computes its effect
rather than being handed one.

**The base state is a prediction, and it is five kelvin low.** This column
equilibrates near 283 K against an observed 288, with a greenhouse effect of 28 K
against 33. Neither was a calibration target: the only number fitted is the CO2
band's wing decay, against a published forcing per doubling, and the base state
falls out.

Five kelvin is roughly the greenhouse warming water vapour would supply if the
vapour this framework has tracked since Stage 3 were not radiatively inert. The
same absence sets the sensitivity: this column's feedback parameter is about
2.5 W m^-2 K^-1 against a Planck-only 3.2 and an Earth-like net 1.2, so it will
*undershoot* the observed record rather than matching it. That is the honest
state, and it is a better place to stand than the accidental agreement an
over-sensitive column produced.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from mfespf.core.errors import ConfigurationError
from mfespf.data.adapter import HistoricalDataAdapter, RecordedSolarSource
from mfespf.data.loaders import load_declared
from mfespf.data.records import TimeSeries, splice
from mfespf.physics.greenhouse import (
    CO2_BAND_RESOLUTION,
    INFRARED_BACKGROUND_ABSORPTION,
    PrescribedGreenhouseModule,
    greenhouse_grid,
)
from mfespf.physics.volcanic import PrescribedVolcanicModule
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

GREENHOUSE_BACKGROUND = INFRARED_BACKGROUND_ABSORPTION
"""The infrared background, re-exported from where the calibration lives.

Not a second copy. The wing decay of the CO2 band is calibrated against this
value, and when the two were defined separately the configurations ran at one
and the calibration was fitted at another.
"""


def co2_record(join_year: int = CO2_JOIN_YEAR) -> tuple[TimeSeries, Any]:
    """The ice-core and instrumental records joined, with the join reported."""
    from datetime import UTC, datetime

    from mfespf.data.loaders import seconds_since_epoch

    ice = next(iter(load_declared("LAWDOME-CO2-ANNUAL")._series.values()))
    station = next(iter(load_declared("MLO-CO2-ANNUAL")._series.values()))
    return splice(ice, station, at=seconds_since_epoch(datetime(join_year, 1, 1, tzinfo=UTC)))


def spin_up(
    *,
    master_dt: float,
    n_layers: int,
    n_band: int,
    ppm: float,
    irradiance: float,
    steps: int,
    optical_depth: float = 0.0,
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

    ``optical_depth`` is here for the same reason. A stratosphere is never quite
    clear, so a control holding the aerosol at the record's background level
    would start slightly dimmed relative to a spin-up that had none, and drift
    by a few hundredths of a kelvin. Small, and it is the whole of what a control
    is supposed to show.
    """
    from mfespf.physics.solar import ConstantSpectrumSolarSource

    grid = greenhouse_grid(n_background=24, n_band=n_band)
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
        volcanic=(
            PrescribedVolcanicModule(optical_depth=optical_depth) if optical_depth > 0.0 else None
        ),
        solar_source=ConstantSpectrumSolarSource(grid=grid, total_irradiance=irradiance),
        initial_surface_temperature=300.0,
        initial_air_temperature=270.0,
    )
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
    n_band: int = CO2_BAND_RESOLUTION,
    spin_up_steps: int = 1500,
    deep_ocean: bool = False,
    volcanic: bool = False,
    hold_volcanic: bool = False,
    **kwargs: Any,
) -> Kernel:
    """The column driven by the irradiance and CO2 records together.

    ``deep_ocean`` puts the two-reservoir ocean under the column. It changes no
    equilibrium -- at equilibrium the two reservoirs are at one temperature and
    the exchange is exactly zero -- so the spin-up is unaffected and the control
    stays flat. What it changes is the path: heat that would have warmed the
    mixed layer goes into a reservoir that takes centuries to fill.

    ``hold_co2`` and ``hold_solar`` freeze one driver at its record mean, which
    is how the two mechanisms are separated: three runs -- both, solar only, CO2
    only -- differ in exactly one thing each, so a difference between any pair is
    attributable to the driver that changed (HYP-007).
    """
    grid = greenhouse_grid(n_background=24, n_band=n_band)
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
    aerosol = None
    if volcanic:
        aerosol = next(iter(load_declared("GISS-SAOD-MONTHLY")._series.values()))
        quiet = float(np.median(aerosol.values))
        # A replay driven by several records is bounded by the shortest of them.
        # The aerosol record ends in 2012, so a run using it covers less than one
        # using the solar and CO2 records alone -- which is why the volcanic
        # experiment reports its own compared period rather than borrowing one.
        clock = clock[clock <= aerosol.coverage[1]]
        if clock.size < 2:
            raise ConfigurationError("the aerosol record does not overlap the replay clock")

    solar_adapter.prepare("total_solar_irradiance", clock, target_interval=master_dt)

    co2, _ = co2_record()
    initial_co2 = float(co2.value_at(clock[0]))

    def concentration(time: float) -> float:
        return initial_co2 if hold_co2 else co2.value_at(time)

    def optical_depth(time: float) -> float:
        # Held at the record's *median*, not its first value: the first value
        # may sit inside an eruption's decay, and a control that started there
        # would be permanently dimmed relative to a quiet stratosphere.
        assert aerosol is not None
        return quiet if hold_volcanic else aerosol.value_at(time)

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
            optical_depth=(quiet if volcanic else 0.0),
        )
        kwargs["initial_surface_temperature"] = surface
        kwargs["initial_air_temperature"] = layers

    kernel = radiating_column(
        deep_ocean=deep_ocean,
        absorption_coefficient=window_absorption(grid, kappa_opaque=GREENHOUSE_BACKGROUND),
        grid=grid,
        n_layers=n_layers,
        areal_heat_capacity=areal_heat_capacity,
        master_dt=master_dt,
        start_time=float(clock[0]),
        operating_mode=OperatingMode.PRESCRIBED_HISTORY,
        datasets=datasets,
        greenhouse=PrescribedGreenhouseModule(concentration=concentration),
        volcanic=(PrescribedVolcanicModule(optical_depth=optical_depth) if volcanic else None),
        solar_source=RecordedSolarSource(
            grid=grid, adapter=solar_adapter, variable="total_solar_irradiance"
        ),
        **kwargs,
    )
    kernel.adapter = solar_adapter
    kernel.manifest = None
    return kernel
