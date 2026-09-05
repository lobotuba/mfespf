"""A synthetic historical record, and why the framework ships one.

Charter section 19 stage 6; requirements INF-002, DAT-009, VAL-006, TIM-006.

Stage 6 needs a declared historical input to replay. This framework has no
observational dataset in it, and inventing numbers that look like a solar
composite would be worse than having none: a plausible-looking record with a
real instrument's name on it is the one thing a provenance system cannot
protect a reader from.

So the record shipped here is generated, says so in every place it can, and is
classified ``SYNTHETIC``. Two consequences follow and both are enforced rather
than described. The dataset registry refuses to accept a record with no declared
provenance class in an evidence role, and the withheld evaluator caps a claim
evaluated against a generated record at ``EMPIRICALLY_UNRESOLVED`` however well
the metric passes.

What that buys is a Stage 6 that is honestly what it is: the replay pathway --
declared input, frozen equations, withheld output, budget reconstruction -- run
end to end and shown to work, with no claim about the Earth anywhere in it. When
an observational record is registered later, nothing in the pathway changes
except the provenance class, and the cap lifts on its own.

The generator is deterministic from its seed, so the record is reproducible from
the manifest and its checksum is stable across platforms (CMP-009).
"""

from __future__ import annotations

import numpy as np

from mfespf.core.errors import ConfigurationError
from mfespf.core.hashing import digest
from mfespf.core.versioning import Version
from mfespf.data.adapter import HistoricalDataAdapter
from mfespf.data.records import (
    ProvenanceClass,
    SamplingConvention,
    TimeSeries,
    UncertaintyKind,
)
from mfespf.registries.datasets import DatasetRecord

__all__ = ["SYNTHETIC_RECORD_ID", "synthetic_irradiance_adapter", "synthetic_irradiance_series"]

SYNTHETIC_RECORD_ID = "SYNTHETIC-TSI-001"
V1 = Version(1, 0, 0)

DAY = 86400.0
YEAR = 365.25 * DAY


def synthetic_irradiance_series(
    *,
    n_samples: int = 480,
    interval: float = 30.0 * DAY,
    baseline: float = 1361.0,
    cycle_amplitude: float = 0.55,
    cycle_period: float = 11.0 * YEAR,
    trend: float = 0.0,
    noise_amplitude: float = 0.08,
    seed: int = 20260904,
) -> TimeSeries:
    """A generated irradiance history with an eleven-year cycle and noise.

    The shape is chosen to exercise the pathway, not to imitate the Sun: a slow
    cycle the model can follow, and a fast component that resampling onto a
    coarse model clock will visibly remove, so the lost-frequency report has
    something to report.
    """
    if n_samples < 2:
        raise ConfigurationError("a record needs at least two samples")
    if interval <= 0.0:
        raise ConfigurationError("sample interval must be positive")

    times = np.arange(n_samples, dtype=float) * interval
    years = times / YEAR
    generator = np.random.default_rng(seed)
    values = (
        baseline
        + cycle_amplitude * np.sin(2.0 * np.pi * times / cycle_period)
        + trend * years
        + noise_amplitude * generator.standard_normal(n_samples)
    )
    return TimeSeries(
        variable="total_solar_irradiance",
        times=times,
        values=values,
        units="W m^-2",
        sampling=SamplingConvention.INTERVAL_MEAN,
        dataset_id=SYNTHETIC_RECORD_ID,
        dataset_version=str(V1),
        provenance_class=ProvenanceClass.SYNTHETIC,
        uncertainty=np.full(n_samples, noise_amplitude),
        # The generator's noise term is drawn from a standard normal scaled by
        # this amplitude, so the amplitude is a standard deviation. Declaring it
        # is not a formality: a half-width of a 95 percent interval would be
        # about twice this, and nothing in the array says which it is (UNC-002).
        uncertainty_kind=UncertaintyKind.STANDARD_UNCERTAINTY,
        note=(
            "GENERATED, NOT MEASURED. Produced by "
            f"mfespf.data.synthetic.synthetic_irradiance_series(seed={seed}). "
            "It resembles a solar composite in shape only and constrains nothing "
            "about the Sun."
        ),
    )


def synthetic_irradiance_adapter(**kwargs: object) -> HistoricalDataAdapter:
    """The record wrapped in an adapter, with a checksum over its own numbers.

    The checksum is taken over the generated values rather than over a file,
    because there is no file: the record is defined by its generator and seed,
    both of which are in the acquisition record. That is a stronger reproduction
    guarantee than a file hash, not a weaker one -- the numbers can be
    regenerated and checked, not merely compared.
    """
    series = synthetic_irradiance_series(**kwargs)  # type: ignore[arg-type]
    checksum = digest(
        {
            "generator": "mfespf.data.synthetic.synthetic_irradiance_series",
            "times": series.times.tolist(),
            "values": series.values.tolist(),
        }
    )
    record = DatasetRecord(
        record_id=SYNTHETIC_RECORD_ID,
        version=V1,
        title="Synthetic total solar irradiance history for validation-ladder stage 6",
        source=(
            "GENERATED by mfespf.data.synthetic. Not an observation, not a reanalysis, and "
            "not derived from any instrument record."
        ),
        checksum=checksum,
        variables=("total_solar_irradiance",),
        units=("W m^-2",),
        temporal_coverage=(f"{series.coverage[0]:g} s", f"{series.coverage[1]:g} s"),
        spatial_coverage="global scalar at 1 au",
        uncertainty_record=(
            "standard_uncertainty: the declared per-sample uncertainty is the standard "
            "deviation of the generator's own noise term. It is not a measurement "
            "uncertainty and must not be reported as one"
        ),
        quality_flags=("synthetic", "no_missing_values", "uniform_sampling"),
        acquisition_record=(
            "deterministic from the recorded seed; regenerate with "
            "synthetic_irradiance_series(**parameters) and compare the checksum"
        ),
        provenance_class=ProvenanceClass.SYNTHETIC,
    )
    return HistoricalDataAdapter(record, {series.variable: series})
