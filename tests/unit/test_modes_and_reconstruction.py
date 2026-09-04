"""Stage 6: operating modes and the whole-run budget.

Requirements TIM-005, TIM-006, VAL-007, CON-009, CON-010, KER-006.

Two small pieces with one thing in common: both exist to catch something that a
per-item check cannot see. An operating mode governs what a run may be *called*,
which no numerical gate constrains; a whole-run reconstruction asks whether a
residual that is negligible in every interval is still negligible after fifty
thousand of them.
"""

from __future__ import annotations

import pytest

from mfespf.core.errors import (
    ConfigurationError,
    InvalidForScientificInferenceError,
    ProvenanceError,
)
from mfespf.coupling.conservation import ConservationReport, CouplingResidual, ModuleBudget
from mfespf.coupling.transfers import ConservedQuantity
from mfespf.evaluation.reconstruction import reconstruct_budget
from mfespf.run.manifest import RunManifest, lineage
from mfespf.run.modes import OperatingMode

pytestmark = pytest.mark.stage(6)


# --- operating modes -------------------------------------------------------


@pytest.mark.req("TIM-005")
def test_every_mode_the_requirement_names_exists():
    """TIM-005 lists seven modes; a framework that supports six supports none."""
    names = {m.value for m in OperatingMode}
    assert names == {
        "equilibrium",
        "controlled_experiment",
        "prescribed_history",
        "free_evolution",
        "historical_reconstruction",
        "hindcast",
        "out_of_sample_prediction",
    }


@pytest.mark.req("TIM-006", "VAL-007")
def test_a_reconstruction_may_not_be_reported_as_a_prediction():
    """VAL-007, and the reason it is the easiest mistake to make honestly.

    The arithmetic of a reconstruction and a prediction is identical. What
    differs is what was known when the model was fixed, and nothing in the
    numbers records that -- so it has to be carried and checked.
    """
    reconstruction = OperatingMode.HISTORICAL_RECONSTRUCTION
    with pytest.raises(InvalidForScientificInferenceError, match="VAL-007"):
        reconstruction.assert_may_be_reported_as(OperatingMode.OUT_OF_SAMPLE_PREDICTION)
    with pytest.raises(InvalidForScientificInferenceError, match="VAL-007"):
        reconstruction.assert_may_be_reported_as(OperatingMode.HINDCAST)


@pytest.mark.req("TIM-006")
def test_claiming_less_than_the_run_supports_is_allowed():
    """The rule is against overstating, not against modesty."""
    OperatingMode.HINDCAST.assert_may_be_reported_as(OperatingMode.HISTORICAL_RECONSTRUCTION)
    OperatingMode.HINDCAST.assert_may_be_reported_as(OperatingMode.CONTROLLED_EXPERIMENT)
    OperatingMode.HINDCAST.assert_may_be_reported_as(OperatingMode.HINDCAST)
    with pytest.raises(InvalidForScientificInferenceError):
        OperatingMode.HINDCAST.assert_may_be_reported_as(OperatingMode.OUT_OF_SAMPLE_PREDICTION)


@pytest.mark.req("TIM-005", "TIM-006")
def test_prescribed_history_is_the_weakest_of_the_comparison_modes():
    """A run told the boundary conditions of the period it reproduces.

    Its agreement is correspondingly cheap, and the ordering says so rather
    than leaving a reader to work it out.
    """
    prescribed = OperatingMode.PRESCRIBED_HISTORY
    assert prescribed.compares_against_a_record
    assert prescribed.uses_target_period_inputs
    assert not prescribed.is_prospective
    assert prescribed.standing < OperatingMode.HISTORICAL_RECONSTRUCTION.standing
    assert "inputs taken from the compared period" in prescribed.describe()


@pytest.mark.req("TIM-005", "AUD-002")
def test_the_mode_is_on_the_manifest_but_not_in_the_configuration_hash():
    """Two runs differing only in declared mode produce identical numbers.

    The configuration hash is for numbers, so the mode stays out of it -- and
    stays on the manifest, because what the mode governs is which language a
    report may use about those numbers.
    """
    first = RunManifest(run_id="A", purpose="p", operating_mode=OperatingMode.HINDCAST)
    second = RunManifest(
        run_id="A", purpose="p", operating_mode=OperatingMode.CONTROLLED_EXPERIMENT
    )
    assert first.configuration_hash() == second.configuration_hash()
    assert first.to_dict()["operating_mode"] == "hindcast"
    with pytest.raises(InvalidForScientificInferenceError):
        second.assert_may_be_reported_as(OperatingMode.HINDCAST)


# --- lineage ---------------------------------------------------------------


@pytest.mark.req("KER-006", "AUD-005")
def test_lineage_walks_a_run_back_to_its_ancestors():
    root = RunManifest(run_id="R0", purpose="baseline")
    middle = RunManifest(run_id="R1", purpose="varied", parent_run_ids=("R0",))
    leaf = RunManifest(run_id="R2", purpose="varied again", parent_run_ids=("R1",))
    known = {m.run_id: m for m in (root, middle, leaf)}

    chain = [m.run_id for m in lineage(leaf, known)]
    assert chain == ["R0", "R1", "R2"]
    # A parent that is not known stops that branch rather than being invented.
    orphan = RunManifest(run_id="R3", purpose="orphan", parent_run_ids=("MISSING",))
    assert [m.run_id for m in lineage(orphan, {"R3": orphan})] == ["R3"]


@pytest.mark.req("KER-006", "SYS-006")
def test_a_lineage_cycle_is_refused():
    a = RunManifest(run_id="A", purpose="a", parent_run_ids=("B",))
    b = RunManifest(run_id="B", purpose="b", parent_run_ids=("A",))
    with pytest.raises(ProvenanceError, match="its own ancestry"):
        lineage(a, {"A": a, "B": b})


# --- whole-run reconstruction ----------------------------------------------


def _report(
    interval: tuple[float, float],
    stored_initial: float,
    stored_final: float,
    entered: float,
    left: float,
) -> ConservationReport:
    return ConservationReport(
        interval=interval,
        budgets=[
            ModuleBudget(
                module_id="surface",
                quantity=ConservedQuantity.ENERGY,
                units="J m^-2",
                stored_initial=stored_initial,
                stored_final=stored_final,
                imports=entered,
                exports=left,
            ),
            ModuleBudget(
                module_id="space",
                quantity=ConservedQuantity.ENERGY,
                units="J m^-2",
                stored_initial=0.0,
                stored_final=left,
                imports=left,
            ),
        ],
        coupling_residuals=[
            CouplingResidual(
                sender="orbit",
                receiver="surface",
                quantity=ConservedQuantity.ENERGY,
                exported=entered,
                imported=entered,
                units="J m^-2",
                start_time=interval[0],
                end_time=interval[1],
            ),
            CouplingResidual(
                sender="surface",
                receiver="space",
                quantity=ConservedQuantity.ENERGY,
                exported=left,
                imported=left,
                units="J m^-2",
                start_time=interval[0],
                end_time=interval[1],
            ),
        ],
    )


@pytest.mark.req("CON-009", "KER-013")
def test_the_whole_run_budget_sums_the_intervals_it_was_given():
    reports = [
        _report((float(i), float(i + 1)), 100.0 * i, 100.0 * (i + 1), 300.0, 200.0)
        for i in range(10)
    ]
    reconstruction = reconstruct_budget(reports)

    assert reconstruction.n_intervals == 10
    assert reconstruction.interval == (0.0, 10.0)
    assert reconstruction.system_modules == ("surface",)
    assert reconstruction.entered == pytest.approx(3000.0)
    assert reconstruction.left == pytest.approx(2000.0)
    assert reconstruction.stored_change == pytest.approx(1000.0)
    assert reconstruction.residual == pytest.approx(0.0)
    assert reconstruction.passed()


@pytest.mark.req("CON-009", "CON-010")
def test_a_drift_too_small_to_see_per_interval_is_visible_over_the_run():
    """The whole point of reconstructing rather than trusting the interval checks.

    Each interval below is off by one part in ten thousand of its own flows,
    which many tolerances would pass. Biased the same way a thousand times, it
    becomes a percent of the run -- and a mean imbalance a reader can compare
    against something physical.
    """
    leak = 0.03  # per interval, against flows of 300
    reports = [
        _report((float(i), float(i + 1)), 100.0 * i, 100.0 * (i + 1) - leak * i, 300.0, 200.0)
        for i in range(1000)
    ]
    reconstruction = reconstruct_budget(reports)

    assert not reconstruction.passed(1e-10)
    assert reconstruction.relative_residual > 1e-4
    assert reconstruction.residual < 0.0  # the system stored less than it received
    assert reconstruction.mean_imbalance == pytest.approx(
        reconstruction.residual / reconstruction.duration
    )
    assert "as a mean imbalance" in reconstruction.describe()


@pytest.mark.req("CON-009")
def test_the_residual_is_scaled_against_the_largest_term_not_the_net_change():
    """A run whose flows nearly cancel is still wrong if it is wrong by a percent
    of those flows."""
    reports = [_report((0.0, 1.0), 0.0, 1.0, 1000.0, 999.0)]
    reconstruction = reconstruct_budget(reports)
    assert reconstruction.scale == pytest.approx(1000.0)
    assert reconstruction.stored_change == pytest.approx(1.0)


@pytest.mark.req("CON-009", "SYS-006")
def test_the_boundary_must_be_named_rather_than_inferred():
    """What counts as the system is a modelling decision, not something to guess.

    Counting the space boundary inside the system would make every run conserve
    trivially, so the reconstruction refuses when nothing is left outside.
    """
    reports = [_report((0.0, 1.0), 0.0, 100.0, 300.0, 200.0)]
    with pytest.raises(ConfigurationError, match="no system"):
        reconstruct_budget(reports, boundary_modules=("surface", "space"))
    with pytest.raises(ConfigurationError, match="no conservation reports"):
        reconstruct_budget([])
