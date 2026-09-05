"""The solar-and-carbon-dioxide attribution experiment.

Charter sections 2.2, 3.5, 14 and 19; requirements HYP-001 through HYP-010,
EXP-001, EXP-002, INF-003, VAL-003, VAL-005, VAL-007, CLM-002, CHM-001.

The companion to ``solar_only``, and the first experiment in this framework
whose outcome was not obvious before it ran. Carbon dioxide enters through the
absorption coefficient, so the radiative transfer computes its effect rather
than being handed one; the concentration is a declared record; and the same
threshold, taken from the same observational record's own uncertainty, decides
the result.

**Four runs, differing in one thing each.** Both drivers varying, solar held,
CO2 held, and both held. The last is exactly flat -- zero range, zero variance
explained -- which is what makes a difference between any pair attributable to
the driver that changed rather than to a transient neither of them owns.

**What the configuration cannot do, stated in advance.** It has no aerosols, no
volcanic forcing, no ocean heat uptake beyond a fifty-metre slab, no
water-vapour feedback, and a base state about four kelvin too cool. Every one of
those bears on the result, and the ones that push in opposite directions are
named in the hypothesis rather than left for a reader to balance.

**The base state changed sign.** An earlier version of this experiment ran a
column that equilibrated near 315 K, twenty-seven kelvin too warm, and the
claims said so. That was not physics: the CO2 band was calibrated at one
background opacity and run at another, so the runs used a forcing 1.75 times the
calibrated one. Recalibrated at the opacity the runs actually use -- and on a
spectral grid where the answer has stopped moving -- the column settles at
283.6 K, which is four kelvin *low*, and the missing warmth is about the size of
the water-vapour greenhouse this configuration does not have. The claims below
are version 2.0.0 for that reason, and the version 1.0.0 results are superseded
rather than deleted (charter 3.5).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from mfespf.core.versioning import Version
from mfespf.data.loaders import load_declared
from mfespf.data.records import TimeSeries
from mfespf.evaluation.comparison import (
    Comparison,
    align_to_record,
    anomaly_against_baseline,
    model_series,
)
from mfespf.experiments.protocol import ExperimentOutcome, ExperimentProtocol
from mfespf.experiments.solar_only import (
    BASELINE_END,
    BASELINE_START,
    compared_period,
    observational_threshold,
)
from mfespf.registries.claims import ClaimRecord, ClaimStatus, ExpectedSignature
from mfespf.registries.datasets import DatasetRole
from mfespf.registries.hypotheses import (
    EntryPoint,
    HypothesisRecord,
    HypothesisRegistry,
    MechanismClass,
)
from mfespf.run.driver import Kernel
from mfespf.run.greenhouse import co2_record, greenhouse_replay
from mfespf.run.historical import replay_steps
from mfespf.run.modes import OperatingMode

__all__ = [
    "CLAIM_ID",
    "EXPERIMENT_ID",
    "HYPOTHESIS_ID",
    "OCEAN_CLAIM_ID",
    "OCEAN_EXPERIMENT_ID",
    "VOLCANIC_CLAIM_ID",
    "SolarAndCO2Result",
    "run_solar_and_co2_experiment",
    "run_solar_co2_ocean_experiment",
    "run_solar_co2_volcanic_experiment",
    "solar_and_co2_claim",
    "solar_and_co2_hypothesis",
    "solar_and_co2_protocol",
]

HYPOTHESIS_ID = "HYP-SOLAR-CO2-001"
CLAIM_ID = "CLM-SOLAR-CO2-001"
OCEAN_CLAIM_ID = "CLM-SOLAR-CO2-OCEAN-001"
EXPERIMENT_ID = "EXP-SOLAR-CO2-001"
OCEAN_EXPERIMENT_ID = "EXP-SOLAR-CO2-OCEAN-001"
VOLCANIC_CLAIM_ID = "CLM-SOLAR-CO2-VOLCANIC-001"
VOLCANIC_EXPERIMENT_ID = "EXP-SOLAR-CO2-VOLCANIC-001"


def _ids(deep_ocean: bool, volcanic: bool = False) -> tuple[str, str]:
    """One hypothesis, three claims.

    The mechanism is the same in all three -- carbon dioxide absorbing in a
    saturated band -- so it is one hypothesis. What differs is the configuration
    the claim is about, and a claim names a configuration. Reusing one claim id
    across them would make each silently supersede the last one's result.

    The volcanic configuration also compares over a shorter period, because the
    aerosol record ends in 2012 and a replay is bounded by the shortest record
    driving it. That alone would make its numbers a different claim.
    """
    if volcanic:
        return (VOLCANIC_CLAIM_ID, VOLCANIC_EXPERIMENT_ID)
    return (OCEAN_CLAIM_ID, OCEAN_EXPERIMENT_ID) if deep_ocean else (CLAIM_ID, EXPERIMENT_ID)


IRRADIANCE = "NRLTSI2-MONTHLY"
TEMPERATURE = "HADCRUT5-GLOBAL-MONTHLY"
DAY = 86400.0


def solar_and_co2_hypothesis(
    *, deep_ocean: bool = False, volcanic: bool = False
) -> HypothesisRecord:
    """One mechanism, two versions, differing in what the experiment can separate.

    The physics is identical. What differs is the configuration it is tested in,
    and that changes the answer to HYP-008: a column with no deep reservoir
    cannot distinguish its result from one in which ocean heat uptake mattered,
    and a column with one can. Since that list is inside the frozen criteria
    hash, the difference has to be a new version rather than an edit -- which is
    what HYP-009 is for.
    """
    return HypothesisRecord(
        record_id=HYPOTHESIS_ID,
        version=Version(1, 2, 0)
        if volcanic
        else (Version(1, 1, 0) if deep_ocean else Version(1, 0, 0)),
        statement=(
            "Variation in total solar irradiance together with the observed rise in "
            "carbon dioxide is sufficient to account for the observed change in "
            "global-mean surface temperature since 1882."
        ),
        claim_ids=(CLAIM_ID, OCEAN_CLAIM_ID, VOLCANIC_CLAIM_ID),
        physical_basis=(
            "Radiative-convective equilibrium of a hydrostatic column. Carbon dioxide "
            "absorbs in a band near 15 microns whose core is saturated, so raising the "
            "concentration widens the opaque region rather than deepening it, and the "
            "level from which the column radiates to space rises. The surface warms "
            "until the outgoing flux again balances the absorbed solar flux."
        ),
        governing_equation_ids=(
            "EQ-CO2-BAND-ABSORPTION",
            "EQ-SCHWARZSCHILD-TRANSFER",
            "EQ-TWO-STREAM-LAYER-COEFFICIENTS",
            "EQ-DRY-CONVECTIVE-ADJUSTMENT",
            "EQ-SLAB-ENERGY-BALANCE",
        ),
        required_inputs=("total_solar_irradiance", "co2_mole_fraction"),
        produced_outputs=("surface_temperature", "outgoing_longwave_radiation"),
        expected=ExpectedSignature(
            sign=(
                "positive: rising carbon dioxide raises the surface temperature, and the "
                "sign is a refutation rather than an inaccuracy if it comes out otherwise"
            ),
            magnitude_bound=(
                "the modelled anomaly tracks the observed one to within the record's own "
                "uncertainty, in root-mean-square, over the compared period"
            ),
            temporal=(
                "follows the concentration record's centennial rise with a lag of order "
                "the mixed layer's relaxation time, and carries the eleven-year solar "
                "cycle as a small superimposed oscillation"
            ),
            spatial="global mean only; the configuration is a single column",
            spectral=(
                "the outgoing radiance is depressed in the 15 micron band and unchanged "
                "in the window, which this configuration produces but which no available "
                "record here can test"
            ),
        ),
        validity_domain=(
            "Global-mean surface temperature on decadal to centennial timescales over the "
            "period both records cover. The column's base state is about 4 K too cool, "
            "so the Planck response is slightly weaker than at the real base state; the "
            "column's feedback parameter is 2.45 W m^-2 K^-1 against a no-feedback Planck "
            "response near 3.2, so it remains over-sensitive to any forcing."
        ),
        failure_criterion=(
            "the root-mean-square difference between the modelled and observed "
            "global-mean temperature anomaly exceeds the observational record's own mean "
            "95 percent half-width over the compared period"
        ),
        conserved_quantities=("energy",),
        transfer_accounting=(
            "ordinary transfers throughout: the mechanism changes an absorption "
            "coefficient and adds no source or sink of its own, so the column's energy "
            "budget closes to round-off with and without it"
        ),
        entry_point=EntryPoint.CONSTITUTIVE_RELATION,
        entry_point_detail=(
            "the absorption coefficient of the radiative transfer, through "
            "co2_band_absorption. The transfer solver is unchanged and unaware: it is "
            "handed a different kappa and computes the consequence, which is what makes "
            "the greenhouse effect a result of this configuration and not an input to it"
        ),
        implementing_modules=("greenhouse", "radiation", "atmosphere", "surface"),
        null_configuration=(
            "the identical configuration with both drivers held at their records' first "
            "values, which produces a flat line: zero range and zero variance explained"
        ),
        distinguishes_from=(
            "a configuration with no varying driver at all, which the null control gives",
            "solar variation alone, which the CO2-held run gives",
            "carbon dioxide alone, which the solar-held run gives",
            "a greenhouse response of the opposite sign",
        ),
        cannot_distinguish_from=(
            "aerosol forcing, which is absent and would have cooled",
            *(
                ()
                if volcanic
                else (
                    "volcanic forcing, which is absent from this configuration and "
                    "would have added decadal structure",
                )
            ),
            *(
                ()
                if deep_ocean
                else (
                    "ocean heat uptake below the mixed layer, which is absent from this "
                    "configuration and would have delayed the warming",
                )
            ),
            "water-vapour feedback, which is absent and would have amplified it",
            "internal variability of the coupled system, which this kernel cannot represent at all",
            "any combination of the above summing to the same global-mean history, since "
            "the experiment observes one global-mean scalar",
        ),
        mechanism_class=MechanismClass.CONVENTIONAL,
        assumptions=(
            "carbon dioxide uniformly mixed, and only its 15 micron band represented",
            "the background infrared absorption is fixed, so there is no water-vapour feedback",
            "planetary albedo fixed at 0.3; no cloud response of any kind",
            "a fifty-metre slab ocean, so heat uptake below the mixed layer is absent",
        ),
        notes=(
            "The first experiment here whose outcome was not obvious before it ran. Two "
            "of the omissions push the model warm -- no aerosols, no deep-ocean heat "
            "uptake -- and one pushes it cool: no water-vapour feedback. Which dominates "
            "was not predictable from the configuration alone."
        ),
    ).freeze()


def solar_and_co2_claim(
    threshold: float, *, deep_ocean: bool = False, volcanic: bool = False
) -> ClaimRecord:
    claim_id, _ = _ids(deep_ocean, volcanic)
    ocean = (
        " over a two-reservoir ocean, so that heat taken up below the mixed layer "
        "delays the response,"
        if deep_ocean
        else ""
    )
    return ClaimRecord(
        record_id=claim_id,
        version=Version(2, 0, 0),
        supersedes=f"{claim_id} 1.0.0",
        statement=(
            f"A radiative-convective column{ocean} driven by the NRLTSI2 irradiance "
            "record and the Law Dome and Mauna Loa carbon dioxide records reproduces the "
            f"HadCRUT5 global-mean temperature anomaly to within {threshold:.4f} K, "
            "root-mean-square, from 1882 onward."
        ),
        mechanism="Greenhouse absorption in a saturated band under a varying solar input.",
        physical_basis="Radiative-convective equilibrium with a spectrally resolved absorber.",
        governing_equation_ids=("EQ-CO2-BAND-ABSORPTION", "EQ-SCHWARZSCHILD-TRANSFER"),
        failure_criterion=(
            f"root-mean-square anomaly difference above {threshold:.6f} K, the record's "
            f"own mean 95 percent half-width over the compared period"
        ),
        support_criterion=f"root-mean-square anomaly difference at or below {threshold:.6f} K",
        expected=ExpectedSignature(
            sign="positive",
            magnitude_bound=f"within {threshold:.4f} K root-mean-square",
            temporal="a centennial rise following the concentration record",
            spatial="global mean only",
            spectral="not testable against any record registered here",
        ),
        assumptions=(
            "no aerosol, volcanic or land-use forcing",
            "no water-vapour or cloud feedback",
            (
                "a fifty-metre mixed layer over a 1200 m reservoir exchanging heat at a "
                "calibrated 0.7 W m^-2 K^-1, so the deep uptake is represented but its "
                "efficiency is a fitted number"
                if deep_ocean
                else "a fifty-metre slab ocean and no deeper heat uptake"
            ),
            "a base state about 4 K cooler than observed, the shortfall being roughly "
            "the size of the absent water-vapour greenhouse effect",
        ),
        validity_domain="Global-mean surface temperature, decadal to centennial, 1882 onward.",
        active_modules=("solar", "orbit", "greenhouse", "radiation", "atmosphere", "surface"),
        inactive_modules=("water", "cloud") + (() if deep_ocean else ("ocean",)),
        prescribed_quantities=("total_solar_irradiance", "co2_mole_fraction"),
        independent_variables=("total_solar_irradiance", "co2_mole_fraction"),
        dependent_variables=("surface_temperature",),
        calibration_data=("published_forcing_per_doubling_3.7_W_m2",),
        validation_data=(TEMPERATURE,),
        uncertainty_method=(
            "the record's published 95 percent interval sets the threshold and is not "
            "propagated; the model carries no uncertainty of its own"
        ),
        control_experiment=(
            "the same configuration with both drivers held at their records' first values"
        ),
        null_hypothesis=(
            "neither driver produces a detectable change relative to the held control"
        ),
        alternative_hypotheses=(
            "aerosol forcing",
            "volcanic forcing",
            "ocean heat uptake",
            "internal variability",
        ),
        numerical_correctness_gates=(
            "whole-run energy budget closes to better than 1e-11 relative",
            "the forcing per doubling is converged in spectral resolution",
        ),
        limitations=(
            "A refutation bears on this configuration's sufficiency, not on the physics "
            "of greenhouse absorption. The configuration omits every other forcing and "
            "every feedback.",
        ),
    ).freeze()


def solar_and_co2_protocol(
    threshold: float, *, deep_ocean: bool = False, volcanic: bool = False
) -> ExperimentProtocol:
    _, experiment_id = _ids(deep_ocean, volcanic)
    protocol = ExperimentProtocol(
        experiment_id=experiment_id,
        question=(
            "Can solar irradiance variation and the observed carbon dioxide rise together "
            "account for the observed change in global-mean surface temperature since 1882?"
        ),
        control_configuration=(
            "the identical column"
            + (" over the same ocean" if deep_ocean else "")
            + " with both drivers held at their records' first values"
        ),
        manipulated_variables=("total_solar_irradiance", "co2_mole_fraction"),
        measured_outputs=(
            "rms_anomaly_difference_K",
            "modelled_anomaly_range_K",
            "observed_anomaly_range_K",
            "variance_explained",
            "variance_explained_solar_only",
            "variance_explained_co2_only",
            *(("variance_explained_volcanic_held",) if volcanic else ()),
        ),
        expected_result=(
            "carbon dioxide accounts for most of the modelled response and solar variation "
            "for very little; the modelled range is of the same order as the observed one; "
            "and the residual is dominated by mechanisms the configuration omits"
            + (
                ". Adding the deep reservoir reduces the modelled range, because heat that "
                "would have warmed the mixed layer goes into it instead. Whether that "
                "improves the agreement depends on whether the amplitude or the timing was "
                "the larger error, and the configuration alone does not say which"
                if deep_ocean
                else ""
            )
        ),
        failure_threshold=(
            f"the claim fails if the root-mean-square anomaly difference exceeds "
            f"{threshold:.6f} K, the record's own mean 95 percent half-width over the "
            f"compared period. Taken from the observations' uncertainty and from no model "
            f"output, and identical to the threshold the solar-only experiment used"
        ),
        analysis_method=(
            "difference the modelled absolute temperature against its own mean over "
            "HadCRUT5's 1961-1990 baseline, resample onto the record's monthly times over "
            "the overlap, and take the root-mean-square residual. Repeat with each driver "
            "held to attribute the response between them"
        ),
        forcing_regime="prescribed_history",
        assumptions=(
            "the column is spun up to equilibrium at the records' first values before the "
            "replay begins, so the control is flat and a difference between runs is not a "
            "transient",
            "the modelled temperature is an end-of-step instantaneous value compared with "
            "a monthly-mean record",
            "the base state is about 4 K too cool, and the column's feedback parameter "
            "is 2.45 W m^-2 K^-1 against a no-feedback Planck response near 3.2",
        ),
        notes=(
            "The outcome was not predictable from the configuration: two omissions push "
            "the model warm and one pushes it cool."
        ),
    )
    protocol.freeze()
    return protocol


@dataclass
class SolarAndCO2Result:
    hypothesis: HypothesisRecord
    claim: ClaimRecord
    protocol: ExperimentProtocol
    outcome: ExperimentOutcome
    comparison: Comparison
    solar_only: Comparison
    co2_only: Comparison
    control: Comparison
    threshold: float
    claim_status: ClaimStatus
    operating_mode: OperatingMode
    worst_relative_residual: float
    co2_join_offset_ppm: float
    deep_ocean: bool = False
    volcanic: bool = False
    volcanic_held: Comparison | None = None

    @property
    def refuted(self) -> bool:
        return self.claim_status is ClaimStatus.REFUTED

    def describe(self) -> str:
        rows = (
            ("everything varying", self.comparison),
            *((("volcanic held", self.volcanic_held),) if self.volcanic_held is not None else ()),
            ("solar only", self.solar_only),
            ("CO2 only", self.co2_only),
            ("neither (control)", self.control),
        )
        lines = [
            self.hypothesis.describe(),
            "",
            self.protocol.describe(),
            "",
            f"{'run':<20}{'range K':>10}{'rms K':>10}{'var. expl.':>12}",
        ]
        for name, comparison in rows:
            lines.append(
                f"{name:<20}{comparison.modelled_range:>10.4f}"
                f"{comparison.rms_error:>10.4f}{comparison.explained_fraction():>12.4f}"
            )
        lines.extend(
            [
                f"{'observed':<20}{self.comparison.observed_range:>10.4f}",
                "",
                f"threshold:    {self.threshold:.6f} K (the record's own mean 95% half-width)",
                f"measured rms: {self.comparison.rms_error:.6f} K",
                f"claim status: {self.claim_status.value}",
                f"operating mode: {self.operating_mode.value}",
                f"deep ocean: {'present' if self.deep_ocean else 'absent'}",
                f"volcanic aerosol: {'present' if self.volcanic else 'absent'}",
                f"CO2 records disagree by {self.co2_join_offset_ppm:+.3f} ppm at the join",
                "",
                "this experiment cannot distinguish the mechanism from:",
                *(f"  - {a}" for a in self.hypothesis.cannot_distinguish_from),
            ]
        )
        return "\n".join(lines)


def _compare(kernel: Kernel, observed: TimeSeries, master_dt: float) -> tuple[Comparison, float]:
    steps = replay_steps(kernel.adapter, "total_solar_irradiance", master_dt)
    run = kernel.run(steps)
    modelled = model_series(
        run.steps, "surface_temperature", units="K", variable="surface_temperature_anomaly"
    )
    anomaly, _ = anomaly_against_baseline(
        modelled, start=BASELINE_START, end=BASELINE_END, label="HadCRUT5 1961-1990"
    )
    return align_to_record(anomaly, observed), run.worst_relative_residual()


def run_solar_and_co2_experiment(
    *,
    master_dt: float = 2.0 * DAY,
    deep_ocean: bool = False,
    volcanic: bool = False,
    **kwargs: Any,
) -> SolarAndCO2Result:
    """Declare, freeze, run four configurations, compare, and record."""
    claim_id, experiment_id = _ids(deep_ocean, volcanic)
    observed = next(iter(load_declared(TEMPERATURE)._series.values()))
    irradiance = load_declared(IRRADIANCE).series("total_solar_irradiance")
    threshold = observational_threshold(observed, over=compared_period(irradiance, observed))

    registry = HypothesisRegistry()
    hypothesis = registry.register(
        solar_and_co2_hypothesis(deep_ocean=deep_ocean, volcanic=volcanic)
    ).record
    claim = solar_and_co2_claim(threshold, deep_ocean=deep_ocean, volcanic=volcanic)
    protocol = solar_and_co2_protocol(threshold, deep_ocean=deep_ocean, volcanic=volcanic)
    _, join = co2_record()

    runs = {}
    residuals = []
    for name, options in (
        ("both", {}),
        *(
            (
                (
                    "volcanic_held",
                    {"hold_volcanic": True, "run_id": f"{experiment_id}-HOLD-VOLCANIC"},
                ),
            )
            if volcanic
            else ()
        ),
        (
            "solar_only",
            {
                "hold_co2": True,
                "hold_volcanic": True,
                "run_id": f"{experiment_id}-HOLD-CO2",
            },
        ),
        ("co2_only", {"hold_solar": True, "run_id": f"{experiment_id}-HOLD-SOLAR"}),
        (
            "neither",
            {
                "hold_co2": True,
                "hold_solar": True,
                "hold_volcanic": True,
                "run_id": f"{experiment_id}-CONTROL",
            },
        ),
    ):
        options.setdefault("run_id", experiment_id)
        kernel = greenhouse_replay(
            master_dt=master_dt,
            deep_ocean=deep_ocean,
            volcanic=volcanic,
            **{**kwargs, **options},
        )
        if name == "both":
            kernel.claims.register(claim)
            registry.assert_modules_present(
                HYPOTHESIS_ID, {s.module_id for s in kernel.schedule.ordered()}
            )
            # Registered through the adapter so the record's checksum, licence and
            # transformation reach the run's manifest, not just its identifier.
            load_declared(TEMPERATURE).register(
                kernel.datasets,
                role=DatasetRole.HISTORICAL_RECONSTRUCTION_EVALUATION,
                claim_id=claim_id,
                run_id=experiment_id,
                note="compared against, not supplied to, the run",
            )
            for dataset_id in ("LAWDOME-CO2-ANNUAL", "MLO-CO2-ANNUAL"):
                load_declared(dataset_id).register(
                    kernel.datasets,
                    role=DatasetRole.BOUNDARY_CONDITION,
                    claim_id=claim_id,
                    run_id=experiment_id,
                    influenced=("greenhouse",),
                    note="prescribed concentration; supplied to the run, not produced by it",
                )
            kernel.datasets.assert_clean(claim_id)
        comparison, residual = _compare(kernel, observed, master_dt)
        runs[name] = comparison
        residuals.append(residual)

    comparison = runs["both"]
    passed = comparison.rms_error <= threshold
    outcome = protocol.evaluate(
        passed=passed,
        measured={
            "rms_anomaly_difference_K": comparison.rms_error,
            "modelled_anomaly_range_K": comparison.modelled_range,
            "observed_anomaly_range_K": comparison.observed_range,
            "variance_explained": comparison.explained_fraction(),
            "variance_explained_solar_only": runs["solar_only"].explained_fraction(),
            "variance_explained_co2_only": runs["co2_only"].explained_fraction(),
            **(
                {"variance_explained_volcanic_held": runs["volcanic_held"].explained_fraction()}
                if volcanic
                else {}
            ),
        },
        summary=(
            f"modelled range {comparison.modelled_range:.4f} K against an observed "
            f"{comparison.observed_range:.4f} K; {comparison.explained_fraction():.1%} of the "
            f"observed variance explained, against "
            f"{runs['solar_only'].explained_fraction():.1%} with carbon dioxide held"
        ),
    )
    return SolarAndCO2Result(
        hypothesis=hypothesis,
        claim=claim,
        protocol=protocol,
        outcome=outcome,
        comparison=comparison,
        solar_only=runs["solar_only"],
        co2_only=runs["co2_only"],
        control=runs["neither"],
        volcanic_held=runs.get("volcanic_held"),
        threshold=threshold,
        claim_status=OperatingMode.PRESCRIBED_HISTORY.status_for_result(passed),
        operating_mode=OperatingMode.PRESCRIBED_HISTORY,
        worst_relative_residual=max(residuals),
        deep_ocean=deep_ocean,
        volcanic=volcanic,
        co2_join_offset_ppm=float(
            join.uncertainty_effect.split("disagree by ")[1].split(" ppm")[0]
        ),
    )


def run_solar_co2_ocean_experiment(**kwargs: Any) -> SolarAndCO2Result:
    """The same experiment with the two-reservoir ocean under the column."""
    return run_solar_and_co2_experiment(deep_ocean=True, **kwargs)


def run_solar_co2_volcanic_experiment(**kwargs: Any) -> SolarAndCO2Result:
    """The same experiment with volcanic aerosol as well.

    Its compared period ends in 2012 where the aerosol record does, so its
    numbers are not directly comparable with the other two. That is why it is a
    third claim rather than a revision of the second.
    """
    return run_solar_and_co2_experiment(deep_ocean=True, volcanic=True, **kwargs)
