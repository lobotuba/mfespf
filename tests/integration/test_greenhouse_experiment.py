"""The greenhouse attribution experiments, and what separates the two claims.

Requirements HYP-001, HYP-005, HYP-007, HYP-008, EXP-001, INF-003, CLM-002,
CON-009, OCN-002.

The declarations are checked without touching any data, because they have to be
registrable before a run or the run shapes them. The runs themselves are gated on
the records being present and use a coarse spectral grid, since what they are
checking is the structure of the comparison rather than a converged forcing --
the converged number has its own test in the unit suite.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.data.loaders import available
from mfespf.evaluation.comparison import (
    align_to_record,
    anomaly_against_baseline,
    model_series,
)
from mfespf.experiments.solar_and_co2 import (
    CLAIM_ID,
    OCEAN_CLAIM_ID,
    solar_and_co2_claim,
    solar_and_co2_hypothesis,
    solar_and_co2_protocol,
)
from mfespf.experiments.solar_only import BASELINE_END, BASELINE_START
from mfespf.registries.hypotheses import EntryPoint
from mfespf.run.greenhouse import GREENHOUSE_BACKGROUND, co2_record, greenhouse_replay
from mfespf.run.historical import replay_steps
from mfespf.run.modes import OperatingMode

pytestmark = pytest.mark.stage(6)

DAY = 86400.0
STEP = 2.0 * DAY
NEEDED = {"NRLTSI2-MONTHLY", "HADCRUT5-GLOBAL-MONTHLY", "LAWDOME-CO2-ANNUAL", "MLO-CO2-ANNUAL"}
needs_data = pytest.mark.skipif(
    not NEEDED <= set(available()),
    reason="the declared records are not in data/raw; run `python tools/fetch_data.py`",
)


def _anomaly(kernel):
    run = kernel.run(replay_steps(kernel.adapter, "total_solar_irradiance", STEP))
    modelled = model_series(
        run.steps, "surface_temperature", units="K", variable="surface_temperature_anomaly"
    )
    anomaly, _ = anomaly_against_baseline(modelled, start=BASELINE_START, end=BASELINE_END)
    from mfespf.data.loaders import load_declared

    observed = next(iter(load_declared("HADCRUT5-GLOBAL-MONTHLY")._series.values()))
    return align_to_record(anomaly, observed), run


# --- the declarations -------------------------------------------------------


@pytest.mark.req("HYP-001", "HYP-005", "CLM-001")
def test_one_hypothesis_carries_two_claims():
    """The mechanism is the same; the configurations are not.

    Carbon dioxide absorbing in a saturated band is one hypothesis. A claim
    names a configuration, and reusing one claim id for two of them would make
    the second silently supersede the first's result.
    """
    hypothesis = solar_and_co2_hypothesis()
    assert set(hypothesis.claim_ids) == {CLAIM_ID, OCEAN_CLAIM_ID}
    assert hypothesis.is_frozen
    assert hypothesis.entry_point is EntryPoint.CONSTITUTIVE_RELATION
    assert "absorption coefficient" in hypothesis.entry_point_detail
    assert "greenhouse" in hypothesis.implementing_modules

    plain = solar_and_co2_claim(0.1)
    ocean = solar_and_co2_claim(0.1, deep_ocean=True)
    assert plain.record_id != ocean.record_id
    assert "two-reservoir ocean" in ocean.statement
    assert "ocean" in plain.inactive_modules
    assert "ocean" not in ocean.inactive_modules


@pytest.mark.req("HYP-008", "HYP-009", "VAL-005")
def test_the_ocean_version_no_longer_lists_heat_uptake_as_unseparated():
    """What an experiment can separate changes when the mechanism is present.

    A column with no deep reservoir cannot distinguish its result from one in
    which ocean heat uptake mattered, because it has none. A column with one
    can, so that entry leaves the list. And because the list is inside the
    frozen criteria hash, the difference has to be a new hypothesis version
    rather than a quiet edit -- which is what HYP-009 is for.
    """
    plain = solar_and_co2_hypothesis()
    ocean = solar_and_co2_hypothesis(deep_ocean=True)

    assert any("ocean heat uptake" in a for a in plain.cannot_distinguish_from)
    assert not any("ocean heat uptake" in a for a in ocean.cannot_distinguish_from)
    assert ocean.version > plain.version
    assert ocean.criteria_hash() != plain.criteria_hash()

    # Everything the configuration still cannot separate is unchanged.
    for hypothesis in (plain, ocean):
        unseparated = " ".join(hypothesis.cannot_distinguish_from)
        assert "aerosol" in unseparated
        assert "water-vapour feedback" in unseparated
        assert "internal variability" in unseparated


@pytest.mark.req("EXP-001", "HYP-007")
def test_both_protocols_freeze_and_declare_the_same_threshold_basis():
    plain = solar_and_co2_protocol(0.113716)
    ocean = solar_and_co2_protocol(0.113716, deep_ocean=True)

    for protocol in (plain, ocean):
        assert protocol.is_frozen
        protocol.assert_unchanged()
        assert "record's own mean 95 percent half-width" in protocol.failure_threshold
        assert "from no model output" in protocol.failure_threshold
        assert "held at their records' first values" in protocol.control_configuration
    assert plain.experiment_id != ocean.experiment_id
    # The ocean protocol says in advance that the direction of the change is
    # not predictable from the configuration.
    assert "does not say which" in ocean.expected_result


# --- the runs ---------------------------------------------------------------


@needs_data
@pytest.mark.req("HYP-007", "VAL-001", "MOD-014")
def test_the_control_is_flat_with_the_ocean_as_well_as_without():
    """Adding a reservoir must not disturb a spun-up state.

    At equilibrium the two reservoirs are at one temperature and the exchange is
    exactly zero, so the control is unchanged by the ocean's presence. If it
    were not, every difference attributed to the ocean below would include a
    transient the ocean created rather than revealed.
    """
    comparison, run = _anomaly(
        greenhouse_replay(
            master_dt=STEP,
            deep_ocean=True,
            hold_co2=True,
            hold_solar=True,
            n_band=16,
            run_id="TEST-OCEAN-CONTROL",
        )
    )
    assert comparison.modelled_range == pytest.approx(0.0, abs=1e-9)
    assert comparison.explained_fraction() == pytest.approx(0.0, abs=1e-9)
    assert run.may_support_a_claim
    assert run.worst_relative_residual() < 1e-11


@needs_data
@pytest.mark.req("OCN-002", "CLM-003", "CON-009")
def test_the_deep_reservoir_reduces_the_modelled_range_and_the_budget_still_closes():
    """The declared direction: heat that would have warmed the surface goes below.

    Whether that *improves* the agreement is a separate question the protocol
    says the configuration cannot answer in advance. What is predictable is the
    sign of the change in amplitude, and that is what this asserts.
    """
    plain, plain_run = _anomaly(
        greenhouse_replay(master_dt=STEP, n_band=16, run_id="TEST-GH-PLAIN")
    )
    ocean, ocean_run = _anomaly(
        greenhouse_replay(master_dt=STEP, deep_ocean=True, n_band=16, run_id="TEST-GH-OCEAN")
    )

    assert ocean.modelled_range < plain.modelled_range
    assert ocean.modelled_range > 0.5 * plain.modelled_range  # delayed, not suppressed
    for run in (plain_run, ocean_run):
        assert run.may_support_a_claim
        assert run.worst_relative_residual() < 1e-11
        assert run.operating_mode is OperatingMode.PRESCRIBED_HISTORY


@needs_data
@pytest.mark.req("CHM-001", "INF-001")
def test_the_carbon_dioxide_the_column_saw_is_the_spliced_record():
    """The concentration reaching the radiation module is the declared record.

    Not a smooth curve fitted to it: the two source records disagree by nearly
    two parts per million where they meet, and the step is left in.
    """
    merged, report = co2_record()
    kernel = greenhouse_replay(master_dt=STEP, n_band=16, run_id="TEST-GH-CO2")
    run = kernel.run(400)

    played = np.array([float(s.outputs["co2_mole_fraction"].values) for s in run.steps])
    expected = np.array([merged.value_at(s.time) for s in run.steps])
    assert np.array_equal(played, expected)
    assert 280.0 < played.min() < played.max() < 300.0  # the 1880s, from the ice core
    assert "disagree by" in report.uncertainty_effect

    radiation = next(s.module for s in kernel.schedule.ordered() if s.module_id == "radiation")
    assert radiation.couple_greenhouse
    assert radiation.greenhouse_background == GREENHOUSE_BACKGROUND
