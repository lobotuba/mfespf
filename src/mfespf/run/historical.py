"""Stage 6 configuration: replaying a declared record.

Charter section 19 stage 6 -- declared historical inputs, withheld output
records, frozen equations, energy-budget reconstruction -- and requirements
DAT-008, INF-001, INF-002, INF-007, KER-006, TIM-001, TIM-004, TIM-005, TIM-006,
VAL-007.

The configuration is the Stage 1 airless body with one thing changed: its solar
irradiance comes from a dataset instead of a constant. That is deliberate. The
physics of Stage 6 is the physics of Stage 1; what is new is everything around
it -- where the boundary condition came from, what was done to it on the way in,
what the run is allowed to be called afterwards, and whether the energy budget
still closes when summed over the whole replay rather than interval by interval.

The record shipped with the framework is synthetic and says so, so this
configuration exercises the pathway and makes no claim about the Earth. Nothing
in the pathway is specific to that: registering an observational record in its
place changes the provenance class and nothing else, and the evaluator's cap on
what may be concluded lifts on its own.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from mfespf.core.errors import ConfigurationError, ProvenanceError
from mfespf.data.adapter import HistoricalDataAdapter, RecordedSolarSource
from mfespf.data.synthetic import synthetic_irradiance_adapter
from mfespf.physics.solar import broadband_grid
from mfespf.registries.datasets import DatasetRegistry, DatasetRole
from mfespf.run.configurations import airless_radiating_body
from mfespf.run.driver import Kernel
from mfespf.run.modes import OperatingMode

__all__ = [
    "REPLAY_CLAIM_ID",
    "historical_replay",
    "observed_replay",
    "replay_clock",
    "replay_steps",
]

REPLAY_CLAIM_ID = "CLM-STAGE6-REPLAY-001"


def replay_clock(adapter: HistoricalDataAdapter, variable: str, master_dt: float) -> np.ndarray:
    """The master time coordinate the replay will run on (TIM-001).

    Built from the record's own coverage rather than from a step count, so the
    run cannot be asked for a boundary condition the record does not hold. A
    replay that wants to be longer needs a longer record, which is the right way
    round.
    """
    if master_dt <= 0.0:
        raise ConfigurationError("master_dt must be positive")
    low, high = adapter.series(variable).coverage
    n = int(np.floor((high - low) / master_dt)) + 1
    if n < 2:
        raise ConfigurationError(
            f"a {master_dt:g} s step does not fit twice inside the record's "
            f"{high - low:g} s coverage"
        )
    return low + np.arange(n, dtype=float) * master_dt


def historical_replay(
    *,
    adapter: HistoricalDataAdapter | None = None,
    variable: str = "total_solar_irradiance",
    master_dt: float = 30.0 * 86400.0,
    resampling_method: str | None = None,
    claim_id: str = REPLAY_CLAIM_ID,
    **kwargs: Any,
) -> Kernel:
    """KER-012-adjacent Stage 6 configuration: the airless body under a record.

    The dataset is registered as a *boundary condition* for the replay claim,
    which is what it is: the irradiance was supplied to the run, not produced by
    it. Declaring it that way is also what makes the contamination check
    meaningful later, since a variable prescribed to a run can never afterwards
    be offered as independent evidence about that run (INF-003).
    """
    working = adapter if adapter is not None else synthetic_irradiance_adapter()
    clock = replay_clock(working, variable, master_dt)
    working.prepare(variable, clock, method=resampling_method, target_interval=master_dt)

    datasets = DatasetRegistry()
    working.register(
        datasets,
        role=DatasetRole.BOUNDARY_CONDITION,
        claim_id=claim_id,
        run_id=str(kwargs.get("run_id", "STAGE6-REPLAY-001")),
        influenced=("solar",),
        note=(
            "prescribed top-of-atmosphere irradiance over the replay period; supplied to "
            "the run, not produced by it"
        ),
    )

    grid = kwargs.pop("spectral_grid", None) or broadband_grid()
    kwargs.setdefault("run_id", "STAGE6-REPLAY-001")
    kwargs.setdefault(
        "purpose",
        "Stage 6 historical replay: a declared irradiance record played through the "
        "airless body, with the transformation, the operating mode and the whole-run "
        "energy budget recorded.",
    )
    kwargs.setdefault("initial_temperature", 254.0)
    kwargs.setdefault("solar_step_multiple", 1)
    kwargs["spectral_grid"] = grid
    kwargs["master_dt"] = master_dt
    # The clock starts where the record does. A replay of a dated record whose
    # run began at zero would ask for a boundary condition outside the record's
    # coverage on its very first step (TIM-001).
    kwargs["start_time"] = float(clock[0])
    kwargs["solar_source"] = RecordedSolarSource(grid=grid, adapter=working, variable=variable)
    # Prescribed history, not reconstruction: the run is driven by a record and
    # is not being compared with one here. What it may be called afterwards is
    # checked rather than assumed (TIM-006, VAL-007).
    kwargs["operating_mode"] = OperatingMode.PRESCRIBED_HISTORY
    kwargs["datasets"] = datasets

    kernel = airless_radiating_body(**kwargs)
    kernel.manifest = None  # rebuilt from the configuration on the next prepare()
    kernel.adapter = working
    return kernel


def replay_steps(adapter: HistoricalDataAdapter | None, variable: str, master_dt: float) -> int:
    """How many master steps the record can support.

    Accepts ``None`` only to refuse it: a kernel with no adapter is not a replay,
    and answering "zero" would let a caller run a configuration that reads no
    record and believe it had replayed one.
    """
    if adapter is None:
        raise ConfigurationError(
            "this kernel has no declared record; replay_steps applies to a replay "
            "configuration, and a run without one has no record length to be bounded by"
        )
    try:
        prepared = adapter.prepared(variable)
    except ProvenanceError:
        return int(replay_clock(adapter, variable, master_dt).size) - 1
    # The prepared clock, not the raw record's: a replay driven by several
    # records is bounded by the shortest of them, and the configuration clips
    # the clock before preparing. Reading the raw coverage here would offer more
    # steps than the run can actually take.
    return int(prepared.times.size) - 1


def observed_replay(dataset_id: str = "NRLTSI2-MONTHLY", **kwargs: Any) -> Kernel:
    """The same replay, driven by a real irradiance record instead of a synthetic one.

    Nothing in the pathway changes. The record is loaded through the same
    adapter, resampled by the same function, registered with the same role --
    and because its provenance class is ``PROXY`` rather than ``SYNTHETIC``, the
    cap the evaluator places on what a comparison may conclude lifts by itself.
    That was the design claim Stage 6 made with the synthetic record, and this is
    the thing that tests it.

    NRLTSI2 is a reconstruction, not a measurement, before the satellite era. It
    is classified ``PROXY`` for that reason, and a claim resting on the early
    part of the record inherits the reconstruction's assumptions rather than
    only an instrument's uncertainty.
    """
    from mfespf.data.loaders import load_declared

    # An override exists for one purpose: a control configuration built from the
    # same record with its variation removed. Anything else passed here is not
    # the dataset it claims to be, and the run_id should say so.
    adapter = kwargs.pop("adapter", None) or load_declared(dataset_id)
    kwargs.setdefault("run_id", f"STAGE6-REPLAY-{dataset_id}")
    kwargs.setdefault(
        "purpose",
        f"Stage 6 historical replay of {dataset_id}: a measured-lineage irradiance record "
        f"played through the airless body, with the transformation, the operating mode and "
        f"the whole-run energy budget recorded.",
    )
    return historical_replay(adapter=adapter, **kwargs)
