"""Modules declare their uncertain parameters, and ensembles are held to it.

UNC-001 and UNC-002, joined to UNC-005's machinery. Before this, the ensemble
code existed and nothing reported through it: the evaporating fraction's
measured range sat in a comment beside a constant, the evaporation coefficient's
dependence on the sensible-heat one lived in a docstring, and five sweeps were
run by hand from scratch scripts whose runs did not know they were members of
anything.

The two refusals here are the old mistakes, replayed. A sweep of the evaporating
fraction over 0.1 to 0.5 is the one this project actually ran, over a range that
excluded the measured value; the declaration now refuses it. Perturbing the
evaporation coefficient on its own is how the fifty-fold Bowen-ratio
inconsistency would come back; the derivation now refuses it.
"""

from __future__ import annotations

from typing import Any

import pytest

from mfespf.core.errors import ConfigurationError
from mfespf.evaluation.bias import OBSERVED_TEMPERATURE_SPREAD_K
from mfespf.evaluation.uncertainty import (
    IntervalMeaning,
    across,
    declared_uncertainty,
    run_ensemble,
    sweep,
)
from mfespf.physics.greenhouse import PrescribedGreenhouseModule
from mfespf.physics.radiative_transfer import window_absorption
from mfespf.run import column_configurations as cc
from mfespf.run.driver import Kernel
from mfespf.run.greenhouse import greenhouse_grid

pytestmark = pytest.mark.stage(3)


def plume_column(**kwargs: Any) -> Kernel:
    grid = greenhouse_grid()
    return cc.radiating_column(
        run_id="UNC-PLUME",
        purpose="Declared-uncertainty checks on the entraining-plume humidity closure.",
        absorption_coefficient=window_absorption(grid, kappa_opaque=1.0e-4),
        grid=grid,
        greenhouse=PrescribedGreenhouseModule(concentration=280.0),
        convection=True,
        plume_humidity=True,
        **kwargs,
    )


def _declared(kernel: Kernel):
    return declared_uncertainty(
        slot.module.report_uncertainty()
        for slot in kernel.schedule.ordered()
        if slot.module.status.participates
    )


def final_surface_temperature(result) -> float:
    return float(result.steps[-1].outputs["surface_temperature"].values.sum())


@pytest.mark.req("UNC-001", "UNC-002", "UNC-004")
def test_the_plume_column_declares_its_ranges_with_their_basis():
    kernel = plume_column()
    declared = _declared(kernel)
    alpha = declared.ranges["evaporating_fraction"]
    assert (alpha.low, alpha.high, alpha.basis.value) == (0.62, 0.81, "measured")
    assert declared.ranges["bulk_plume_parameter"].basis.value == "asserted"

    manifest = kernel.build_manifest()
    humidity = next(
        entry
        for entry in manifest.declared_uncertainty.values()
        if any(r["setting"] == "evaporating_fraction" for r in entry["declared_ranges"])
    )
    assert humidity["distributional_assumptions"].startswith("a range, not a distribution")
    # The comment that said the fraction was "assumed, not measured" was stale;
    # what is still unresolved is structural and stays marked so (UNC-004).
    assert not any("evaporating" in u for u in humidity["unresolved"])
    assert humidity["unresolved"]


@pytest.mark.req("UNC-001", "UNC-005")
def test_the_sweep_this_project_actually_ran_is_now_refused():
    """0.1 to 0.5, which put the measured 0.67 outside the range being studied."""
    design = sweep(
        "evaporating_fraction",
        [0.1, 0.3, 0.5],
        basis="asserted: the original range, before anyone looked the number up",
    )
    built: list[str] = []

    def build(**overrides):
        built.append("built")
        return plume_column(**overrides)

    with pytest.raises(ConfigurationError, match="outside its declared range"):
        run_ensemble(design, build=build, observe=final_surface_temperature, n_steps=2)
    assert built == ["built"], "refused before the first member ran, not after"


@pytest.mark.req("UNC-002", "UNC-005")
def test_a_derived_parameter_cannot_be_perturbed_on_its_own():
    kernel = cc.moist_radiative_convective()
    declared = _declared(kernel)
    assert declared.derived["evaporation_coefficient"].follows == "sensible_heat_coefficient"
    design = sweep("evaporation_coefficient", [5e-3, 1.5e-2], basis="asserted: for the test")
    with pytest.raises(ConfigurationError, match="Perturb 'sensible_heat_coefficient' instead"):
        run_ensemble(
            design,
            build=lambda **o: cc.moist_radiative_convective(**o),
            observe=final_surface_temperature,
            n_steps=2,
        )


@pytest.mark.req("UNC-005", "UNC-002", "AUD-001")
def test_every_member_knows_the_ensemble_it_belongs_to():
    declared = _declared(plume_column()).ranges["evaporating_fraction"]
    design = across(declared)
    manifests = []

    def surface_temperature(result) -> float:
        manifests.append(result.manifest)
        return final_surface_temperature(result)

    result = run_ensemble(design, build=plume_column, observe=surface_temperature, n_steps=4)

    memberships = [m.ensemble_membership for m in manifests]
    assert {m["design_id"] for m in memberships} == {design.design_id()}
    assert [m["member"] for m in memberships] == list(result.labels)
    assert all(m["n_members"] == 2 and m["method"] == "bounds" for m in memberships)
    assert all(m["basis"].startswith("measured: alpha = 1 - PE") for m in memberships)

    assert result.member_runs == tuple(m.configuration_hash() for m in manifests)
    assert len(set(result.member_runs)) == 2
    # Membership is metadata about the run, not configuration of it: a member
    # hashes exactly like the same configuration run alone.
    alone = plume_column(evaporating_fraction=declared.low).run(4)
    assert alone.manifest.ensemble_membership == {}
    assert alone.manifest.configuration_hash() == result.member_runs[0]


@pytest.mark.req("UNC-001", "AUD-001")
def test_a_run_outside_its_declared_range_is_allowed_and_never_unremarked():
    manifest = plume_column(evaporating_fraction=0.3).build_manifest()
    assert any("CONFIGURED OUTSIDE ITS DECLARED RANGE" in w for w in manifest.warnings)
    assert not any("OUTSIDE" in w for w in plume_column().build_manifest().warnings), (
        "the nominal configuration must not warn"
    )


#: As in test_surface_temperature_error.py, plus the step count that lets the
#: upper end of the range settle with margin: at k_H = 14.5 the imbalance is
#: 0.096 W m^-2 at 16000 steps and 0.007 at 24000, with the time-mean moving
#: 0.04 K between them.
FULL_WATER_CYCLE = dict(
    grid=cc.column_grid(16),
    n_layers=8,
    areal_heat_capacity=1.0e6,
    master_dt=3600.0,
    initial_water_store=1.0e7,
    radiative_vapour=True,
    vapour_transport=True,
)
SETTLED_STEPS = 20_000


@pytest.mark.req("UNC-005", "UNC-007", "KER-012", "VAL-005")
def test_the_surface_exchange_coefficient_cannot_be_the_warm_bias(report_bias):
    """The error in kelvin across the declared range of the one surface parameter.

    The sensible heat coefficient fixes both surface fluxes, so it is the first
    suspect for a warm surface. Across its declared range the error never comes
    near zero, and the *low* end is the warmest: a weaker coupling traps more
    heat at the surface. Whatever the +6.7 K is, no value of this parameter
    anyone has declared plausible removes it.

    The interval is the full range with no coverage, because the members were
    chosen, not drawn.
    """
    declared = _declared(cc.moist_radiative_convective(**FULL_WATER_CYCLE)).ranges[
        "sensible_heat_coefficient"
    ]
    design = across(declared, points=3)
    assert [m.overrides["sensible_heat_coefficient"] for m in design.members] == [6.0, 10.0, 14.5]

    reports = []

    def surface_temperature_error_k(result) -> float:
        report = report_bias(
            result,
            label=f"full water cycle, {result.manifest.ensemble_membership['member']}",
        )
        reports.append(report)
        assert report.error_k is not None, report.describe()
        return report.error_k

    result = run_ensemble(
        design,
        build=lambda **o: cc.moist_radiative_convective(**FULL_WATER_CYCLE, **o),
        observe=surface_temperature_error_k,
        n_steps=SETTLED_STEPS,
    )

    interval = result.interval()
    assert interval.meaning is IntervalMeaning.FULL_RANGE
    assert interval.nominal_coverage is None
    assert interval.low > 6.0, interval.describe()
    assert interval.low > 10 * OBSERVED_TEMPERATURE_SPREAD_K
    assert interval.high < 9.0, interval.describe()
    weakest, nominal, _strongest = result.values
    assert weakest > nominal + 1.0, "the weakest coupling should be the warmest"
