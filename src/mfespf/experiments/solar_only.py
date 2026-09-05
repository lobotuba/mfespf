"""The solar-only attribution experiment.

Charter sections 2.2, 3.5, 14 and 19; requirements HYP-001 through HYP-010,
EXP-001, EXP-002, INF-003, VAL-003, VAL-005, VAL-007, CLM-002, TIM-005, TIM-006.

This is the framework's first experiment against a real observational record,
and it is designed to fail. The kernel has no carbon dioxide, no aerosols and no
volcanic forcing; driven by a solar irradiance reconstruction alone it cannot
reproduce the observed temperature record, and the useful thing is to say so on
the record with the threshold declared beforehand.

**Why a failure is worth more here than a success.** The run is handed the
irradiance of the period it is compared against, so it is a reconstruction and
not a prediction. Reproducing observations under those conditions would be weak
evidence -- the model was given part of the answer. *Failing* to reproduce them
under those conditions is not weak at all: the mechanism was given every
advantage the experiment could offer and still could not account for what
happened. The asymmetry is enforced in ``OperatingMode.status_for_result``
rather than left to the reader.

**What this experiment cannot do.** It observes one global-mean scalar. It can
separate the solar hypothesis from a null control and from a response of the
opposite sign, and it cannot separate it from any other mechanism that would
produce the same global-mean history. That list is in the hypothesis record,
because a mechanism supported by an experiment that could not distinguish it
from four alternatives is supported jointly with those four.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import UTC, datetime

import numpy as np

from mfespf.core.versioning import Version
from mfespf.data.loaders import load_declared, seconds_since_epoch
from mfespf.data.records import TimeSeries
from mfespf.evaluation.comparison import (
    Comparison,
    align_to_record,
    anomaly_against_baseline,
    model_series,
)
from mfespf.experiments.protocol import ExperimentOutcome, ExperimentProtocol
from mfespf.physics.radiation import equilibrium_temperature
from mfespf.registries.claims import ClaimRecord, ClaimStatus, ExpectedSignature
from mfespf.registries.datasets import DatasetRole, DatasetUse, TransformationStep
from mfespf.registries.hypotheses import (
    EntryPoint,
    HypothesisRecord,
    HypothesisRegistry,
    MechanismClass,
)
from mfespf.run.driver import RunResult
from mfespf.run.historical import observed_replay, replay_steps
from mfespf.run.modes import OperatingMode

__all__ = [
    "BASELINE_END",
    "BASELINE_START",
    "CLAIM_ID",
    "EXPERIMENT_ID",
    "HYPOTHESIS_ID",
    "SolarOnlyResult",
    "compared_period",
    "observational_threshold",
    "run_solar_only_experiment",
    "solar_only_claim",
    "solar_only_hypothesis",
    "solar_only_protocol",
]

HYPOTHESIS_ID = "HYP-SOLAR-ONLY-001"
CLAIM_ID = "CLM-SOLAR-ONLY-001"
EXPERIMENT_ID = "EXP-SOLAR-ONLY-001"

IRRADIANCE = "NRLTSI2-MONTHLY"
TEMPERATURE = "HADCRUT5-GLOBAL-MONTHLY"
MONTH = 30.0 * 86400.0

# HadCRUT5's own stated baseline. Differencing the model against any other
# window would shift it by a constant that a reader could mistake for a bias.
BASELINE_START = seconds_since_epoch(datetime(1961, 1, 1, tzinfo=UTC))
BASELINE_END = seconds_since_epoch(datetime(1991, 1, 1, tzinfo=UTC))


def observational_threshold(record: TimeSeries, *, over: tuple[float, float]) -> float:
    """The failure threshold, taken from the record's own uncertainty.

    The number has to come from somewhere that is not the model's error, or the
    experiment tests nothing. The record publishes a 95 percent interval per
    sample; its mean half-width over the compared period is the level below which
    a disagreement is indistinguishable from the record's own noise, and that is
    the natural place to put the line.

    ``over`` is the compared period, which is knowable before the run because it
    is the overlap of the two records' coverage -- a property of the data alone.
    Taking the mean over the whole record instead would set the threshold partly
    from years the comparison never reaches.
    """
    if record.uncertainty is None:
        raise ValueError(f"{record.variable!r} publishes no uncertainty to set a threshold from")
    inside = np.logical_and(record.times >= over[0], record.times <= over[1])
    if not np.any(inside):
        raise ValueError(f"{record.variable!r} has no samples in [{over[0]:g}, {over[1]:g}]")
    return float(np.mean(record.uncertainty[inside]))


def compared_period(modelled: TimeSeries, observed: TimeSeries) -> tuple[float, float]:
    """The overlap of two records' coverage, known before either is used."""
    return (
        max(modelled.coverage[0], observed.coverage[0]),
        min(modelled.coverage[1], observed.coverage[1]),
    )


# --- the declarations -------------------------------------------------------


def solar_only_hypothesis() -> HypothesisRecord:
    """HYP-001..008, frozen."""
    return HypothesisRecord(
        record_id=HYPOTHESIS_ID,
        version=Version(1, 0, 0),
        statement=(
            "Variation in total solar irradiance is sufficient to account for the observed "
            "change in global-mean surface temperature since 1882."
        ),
        claim_ids=(CLAIM_ID,),
        physical_basis=(
            "The first law applied to an isothermal surface layer whose only energy input "
            "is absorbed solar radiation and whose only loss is thermal emission. A change "
            "in irradiance changes the equilibrium temperature through the Stefan-Boltzmann "
            "law, and the layer approaches that equilibrium on a timescale set by its heat "
            "capacity."
        ),
        governing_equation_ids=(
            "EQ-SLAB-ENERGY-BALANCE",
            "EQ-GREYBODY-EMISSION",
            "EQ-SPHERICAL-MEAN-INSOLATION",
            "EQ-RADIATIVE-EQUILIBRIUM-TEMPERATURE",
        ),
        required_inputs=("total_solar_irradiance",),
        produced_outputs=("surface_temperature",),
        expected=ExpectedSignature(
            sign=(
                "positive: a rise in irradiance raises the temperature, and the observed "
                "warming requires a rising irradiance to explain it"
            ),
            magnitude_bound=(
                "the modelled anomaly must track the observed one to within the "
                "observational uncertainty of the record, in root-mean-square"
            ),
            temporal=(
                "follows the record's eleven-year solar cycle with a lag of order the "
                "mixed layer's relaxation time, about one to two years, and follows any "
                "longer-term irradiance trend on the same timescale"
            ),
            spatial=(
                "not applicable at this fidelity beyond the global mean: the configuration "
                "is zero-dimensional, and no spatial signature is claimed or testable"
            ),
            spectral=(
                "not applicable: the irradiance record is broadband and the surface "
                "albedo is prescribed, so no spectral signature is claimed"
            ),
        ),
        validity_domain=(
            "Global-mean surface temperature on annual to centennial timescales, over the "
            "period the irradiance record covers. Nothing about regional temperature, "
            "seasonal cycle, or any observable other than the global mean is in scope."
        ),
        failure_criterion=(
            "the root-mean-square difference between the modelled and observed global-mean "
            "temperature anomaly exceeds the observational record's own mean 95 percent "
            "half-width over the compared period"
        ),
        conserved_quantities=("energy",),
        transfer_accounting=(
            "ordinary transfers: irradiance crosses the top-of-atmosphere boundary as an "
            "import recorded by the orbital module, and reflected and emitted radiation "
            "leave as exports recorded by both the surface and the space boundary. The "
            "mechanism adds no source or sink of its own."
        ),
        entry_point=EntryPoint.BOUNDARY,
        entry_point_detail=(
            "the solar-source interface: a recorded irradiance series replaces the constant "
            "one, which is a boundary condition and not a change to any equation"
        ),
        implementing_modules=("solar", "orbit", "surface"),
        null_configuration=(
            "the identical configuration with a constant irradiance equal to the record's "
            "mean, so that the solar variation is removed and nothing else is"
        ),
        distinguishes_from=(
            "no forcing at all, which the constant-irradiance control provides",
            "a solar response of the opposite sign",
            "a solar response larger than the record's variation could produce",
        ),
        cannot_distinguish_from=(
            "greenhouse forcing, which is absent from this configuration entirely",
            "aerosol forcing, likewise absent",
            "volcanic forcing, likewise absent",
            "internal variability of the coupled system, which this kernel cannot represent",
            "any other mechanism producing the same global-mean history, since the "
            "experiment observes one global-mean scalar and nothing else",
        ),
        mechanism_class=MechanismClass.CONVENTIONAL,
        assumptions=(
            "planetary albedo constant at 0.3 over the whole period",
            "emissivity unity, so the surface radiates as a black body",
            "no atmosphere, and therefore no greenhouse effect of any kind",
            "a single well-mixed reservoir with one response timescale",
        ),
        notes=(
            "Registered in order to be tested, not in order to be supported. The "
            "configuration is known in advance to omit the mechanisms that dominate the "
            "observed record, and the experiment's value is in recording that omission's "
            "consequence against a declared threshold rather than asserting it."
        ),
    ).freeze()


def solar_only_claim(threshold: float) -> ClaimRecord:
    """CLM-001..002, frozen. The threshold comes from the record, not the model."""
    return ClaimRecord(
        record_id=CLAIM_ID,
        version=Version(1, 0, 0),
        statement=(
            "The airless radiating body driven by the NRLTSI2 total solar irradiance "
            "reconstruction reproduces the HadCRUT5 global-mean surface temperature anomaly "
            f"to within {threshold:.4f} K, root-mean-square, from 1882 onward."
        ),
        mechanism="Radiative balance of an isothermal surface under a varying solar input.",
        physical_basis="First law plus grey-body emission; no atmosphere.",
        governing_equation_ids=("EQ-SLAB-ENERGY-BALANCE", "EQ-GREYBODY-EMISSION"),
        failure_criterion=(
            f"root-mean-square anomaly difference above {threshold:.6f} K, which is the "
            f"record's own mean 95 percent half-width over the compared period"
        ),
        support_criterion=f"root-mean-square anomaly difference at or below {threshold:.6f} K",
        expected=ExpectedSignature(
            sign="positive",
            magnitude_bound=f"within {threshold:.4f} K root-mean-square",
            temporal="tracks the irradiance record with a one-to-two-year lag",
            spatial="global mean only",
            spectral="not applicable",
        ),
        assumptions=(
            "albedo 0.3 and emissivity 1 throughout",
            "no greenhouse, aerosol or volcanic forcing in the configuration",
        ),
        validity_domain="Global-mean surface temperature, annual to centennial, 1882 onward.",
        active_modules=("solar", "orbit", "surface", "space"),
        inactive_modules=("radiation", "atmosphere", "water", "cloud", "convection", "ocean"),
        prescribed_quantities=("total_solar_irradiance",),
        independent_variables=("total_solar_irradiance",),
        dependent_variables=("surface_temperature",),
        calibration_data=(),
        validation_data=(TEMPERATURE,),
        uncertainty_method=(
            "the record's published 95 percent interval is used to set the threshold and is "
            "not propagated; the model carries no uncertainty of its own"
        ),
        control_experiment=("the same configuration with constant irradiance at the record's mean"),
        null_hypothesis=(
            "solar irradiance variation produces no detectable change in the modelled "
            "global-mean temperature relative to the constant-irradiance control"
        ),
        alternative_hypotheses=(
            "greenhouse forcing",
            "aerosol forcing",
            "volcanic forcing",
            "internal variability",
        ),
        numerical_correctness_gates=(
            "whole-run energy budget closes to better than 1e-12 relative",
            "every coupling interval reconciles",
        ),
        limitations=(
            "The configuration omits every mechanism except solar forcing, so a refutation "
            "bears on this configuration's sufficiency and not on the physics of solar "
            "forcing itself.",
        ),
    ).freeze()


def solar_only_protocol(threshold: float) -> ExperimentProtocol:
    """EXP-001, frozen before the run."""
    protocol = ExperimentProtocol(
        experiment_id=EXPERIMENT_ID,
        question=(
            "Can variation in total solar irradiance alone account for the observed change "
            "in global-mean surface temperature since 1882?"
        ),
        control_configuration=(
            "the identical airless configuration driven by a constant irradiance equal to "
            "the record's mean over the same period"
        ),
        manipulated_variables=("total_solar_irradiance",),
        measured_outputs=(
            "rms_anomaly_difference_K",
            "modelled_anomaly_range_K",
            "observed_anomaly_range_K",
            "variance_explained",
        ),
        expected_result=(
            "the modelled anomaly follows the irradiance record with a one-to-two-year lag "
            "and a range far smaller than the observed one, because the configuration omits "
            "every other forcing"
        ),
        failure_threshold=(
            f"the claim fails if the root-mean-square anomaly difference exceeds "
            f"{threshold:.6f} K, the record's own mean 95 percent half-width over the "
            f"compared period. The threshold is taken from the observations' uncertainty "
            f"and not from any model output"
        ),
        analysis_method=(
            "difference the modelled absolute temperature against its own mean over "
            "HadCRUT5's 1961-1990 baseline, resample the model onto the record's monthly "
            "times over their overlap, and take the root-mean-square residual"
        ),
        forcing_regime="prescribed_history",
        assumptions=(
            "the modelled temperature is an end-of-step instantaneous value compared with a "
            "monthly-mean record; at a one-month step and a one-to-two-year response time "
            "the difference between the two is small but is an approximation, not an "
            "identity",
            "the record's uncertainty is used to set the threshold and is not propagated",
        ),
        notes=(
            "Registered in the expectation that it fails. A reconstruction that reproduces "
            "observations it was given the boundary conditions for is weak evidence; one "
            "that fails to is not weak at all."
        ),
    )
    protocol.freeze()
    return protocol


# --- running it -------------------------------------------------------------


@dataclass
class SolarOnlyResult:
    """Everything the experiment produced, including what it could not settle."""

    hypothesis: HypothesisRecord
    claim: ClaimRecord
    protocol: ExperimentProtocol
    outcome: ExperimentOutcome
    comparison: Comparison
    control_comparison: Comparison
    threshold: float
    claim_status: ClaimStatus
    operating_mode: OperatingMode
    worst_relative_residual: float
    transformations: tuple[TransformationStep, ...] = ()
    """Every operation performed between the model's output and the comparison.

    Kept because the comparison is only meaningful given them: an anomaly
    differenced against the wrong window and a model resampled the wrong way
    both produce a number that looks like a result.
    """

    @property
    def refuted(self) -> bool:
        return self.claim_status is ClaimStatus.REFUTED

    def describe(self) -> str:
        lines = [
            self.hypothesis.describe(),
            "",
            self.protocol.describe(),
            "",
            "measured:",
            *(f"  {line}" for line in self.comparison.describe().splitlines()),
            "",
            "constant-irradiance control:",
            *(f"  {line}" for line in self.control_comparison.describe().splitlines()),
            "",
            f"threshold:    {self.threshold:.6f} K (the record's own mean 95% half-width)",
            f"measured rms: {self.comparison.rms_error:.6f} K",
            f"claim status: {self.claim_status.value}",
            f"operating mode: {self.operating_mode.value}",
            "",
            "this experiment cannot distinguish the mechanism from:",
            *(f"  - {alternative}" for alternative in self.hypothesis.cannot_distinguish_from),
        ]
        return "\n".join(lines)


def _anomaly_from_run(
    result: RunResult, baseline_label: str
) -> tuple[TimeSeries, tuple[TransformationStep, ...]]:
    absolute = model_series(
        result.steps,
        "surface_temperature",
        units="K",
        variable="surface_temperature_anomaly",
    )
    anomaly, step = anomaly_against_baseline(
        absolute, start=BASELINE_START, end=BASELINE_END, label=baseline_label
    )
    return anomaly, (step,)


def run_solar_only_experiment(*, master_dt: float = MONTH) -> SolarOnlyResult:
    """Declare, freeze, run, compare, and record the result -- in that order."""
    irradiance = load_declared(IRRADIANCE)
    temperature_adapter = load_declared(TEMPERATURE)
    observed = next(iter(temperature_adapter._series.values()))
    series = irradiance.series("total_solar_irradiance")

    # The compared period and the threshold are both properties of the two
    # records, fixed before anything runs. Deriving either from model output
    # would make the experiment a description of its own result.
    period = compared_period(series, observed)
    threshold = observational_threshold(observed, over=period)
    hypothesis_registry = HypothesisRegistry()
    hypothesis = hypothesis_registry.register(solar_only_hypothesis()).record
    claim = solar_only_claim(threshold)
    protocol = solar_only_protocol(threshold)

    # The run. Its irradiance is a boundary condition and its temperature record
    # is validation data: declaring both roles is what makes the contamination
    # check meaningful, since a variable prescribed to a run can never afterwards
    # be independent evidence about it (INF-003).
    # Start at the equilibrium the record's own mean irradiance implies. Left at
    # an arbitrary initial temperature the run spends its first decade relaxing,
    # and that transient dominates the modelled range -- a number that would then
    # describe the initial condition rather than the Sun. The control run below
    # is what exposed it: a constant-irradiance configuration should have no
    # range at all, and it had half a kelvin.
    mean_irradiance = float(np.mean(series.values))
    spun_up = equilibrium_temperature(mean_irradiance, 0.3, 1.0)

    kernel = observed_replay(
        IRRADIANCE,
        master_dt=master_dt,
        run_id="EXP-SOLAR-ONLY-001",
        initial_temperature=spun_up,
    )
    kernel.claims.register(claim)
    temperature_adapter.register(
        kernel.datasets,
        role=DatasetRole.HISTORICAL_RECONSTRUCTION_EVALUATION,
        claim_id=CLAIM_ID,
        run_id="EXP-SOLAR-ONLY-001",
        note="compared against, not supplied to, the run",
    )
    kernel.datasets.declare_use(
        DatasetUse(
            dataset_id=IRRADIANCE,
            dataset_version="1.0.0",
            role=DatasetRole.BOUNDARY_CONDITION,
            claim_id=CLAIM_ID,
            variables=("total_solar_irradiance",),
            influenced=("solar",),
        )
    )
    kernel.datasets.assert_clean(CLAIM_ID)
    hypothesis_registry.assert_modules_present(
        HYPOTHESIS_ID, {s.module_id for s in kernel.schedule.ordered()}
    )

    steps = replay_steps(kernel.adapter, "total_solar_irradiance", master_dt)
    run = kernel.run(steps)
    modelled, transformations = _anomaly_from_run(run, "HadCRUT5 1961-1990")
    comparison = align_to_record(modelled, observed)

    # The control (HYP-007): the same configuration with the variation removed
    # and nothing else changed.
    flat = series.with_values(
        series.times,
        np.full(series.times.size, mean_irradiance),
        note="constant irradiance at the record's own mean; the null control",
    )
    from mfespf.data.adapter import HistoricalDataAdapter

    control_kernel = observed_replay(
        IRRADIANCE,
        adapter=HistoricalDataAdapter(irradiance.record, {flat.variable: flat}),
        master_dt=master_dt,
        run_id="EXP-SOLAR-ONLY-CONTROL",
        initial_temperature=spun_up,
    )
    control_run = control_kernel.run(steps)
    control_modelled, _ = _anomaly_from_run(control_run, "HadCRUT5 1961-1990")
    control_comparison = align_to_record(control_modelled, observed)

    passed = comparison.rms_error <= threshold
    outcome = protocol.evaluate(
        passed=passed,
        measured={
            "rms_anomaly_difference_K": comparison.rms_error,
            "modelled_anomaly_range_K": comparison.modelled_range,
            "observed_anomaly_range_K": comparison.observed_range,
            "variance_explained": comparison.explained_fraction(),
        },
        summary=(
            f"modelled range {comparison.modelled_range:.4f} K against an observed "
            f"{comparison.observed_range:.4f} K; rms residual "
            f"{comparison.rms_error:.4f} K against a threshold of {threshold:.4f} K"
        ),
    )
    status = run.operating_mode.status_for_result(passed)
    return SolarOnlyResult(
        hypothesis=hypothesis,
        claim=claim,
        protocol=protocol,
        outcome=outcome,
        comparison=comparison,
        control_comparison=control_comparison,
        threshold=threshold,
        claim_status=status,
        operating_mode=run.operating_mode,
        worst_relative_residual=run.worst_relative_residual(),
        transformations=(*transformations, *comparison.transformations),
    )
