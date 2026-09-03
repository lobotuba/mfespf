"""Stage 5 of the validation ladder: thermal storage and response timescales.

Charter section 19 stage 5 -- surface heat capacity, mixed-layer ocean, deeper
reservoir, step, pulse and periodic forcing, measured response timescales -- and
requirements KER-012, OCN-001, OCN-002, TIM-007 and TIM-008.

These two configurations complete KER-012: all eight the specification names now
exist and run.

The whole stage rests on one comparison. ``surface_thermal_storage`` and
``reduced_ocean_heat_storage`` are the same configuration, built by the same
function, differing only in whether a deep reservoir is present. Anything that
appears in the second and not the first is attributable to the reservoir and to
nothing else, which is what makes this an experiment rather than two runs.

The reservoir is deliberately shallower and the exchange stronger here than the
reference values, so that both timescales fall inside a record a test can afford
to produce. The measured numbers are therefore properties of *this
configuration*, not of Earth; the CLI reference run uses the reference values.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from mfespf.core.errors import InvalidForScientificInferenceError
from mfespf.evaluation.timescales import (
    TimescaleMeasurement,
    TimescaleReport,
    efolding_time,
    fit_exponential_modes,
    periodic_response,
)
from mfespf.experiments.forcing import (
    ImpulseForcing,
    NoForcing,
    PeriodicForcing,
    PulseForcing,
    StepForcing,
)
from mfespf.experiments.protocol import ExperimentProtocol
from mfespf.physics.ocean import (
    OCEAN_DEEP_STATE,
    areal_heat_capacity_of_water,
    two_reservoir_modes,
    two_reservoir_transfer_function,
)
from mfespf.physics.radiation import equilibrium_temperature
from mfespf.run.configurations import reduced_ocean_heat_storage, surface_thermal_storage

pytestmark = pytest.mark.stage(5)

DAY = 86400.0
DT = 4.0 * DAY
STEPS = 4000

# A 200 m reservoir with a strong exchange: the slow mode is then about twenty
# times the fast one instead of two hundred, and a record long enough to
# constrain both is 4000 steps rather than 40000.
DEPTH = 200.0
GAMMA = 2.0
C_MIX = areal_heat_capacity_of_water(50.0)
C_DEEP = areal_heat_capacity_of_water(DEPTH)

_SIGMA = 5.670374419e-8
T_EQ = equilibrium_temperature(1361.0, 0.3, 1.0)
LAMBDA = 4.0 * _SIGMA * T_EQ**3

BASE = dict(master_dt=DT, solar_step_multiple=1)
OCEAN = dict(deep_ocean_depth=DEPTH, deep_ocean_exchange=GAMMA, **BASE)

ANALYTIC_FAST, ANALYTIC_SLOW = two_reservoir_modes(
    mixed_layer_capacity=C_MIX,
    deep_capacity=C_DEEP,
    feedback_parameter=LAMBDA,
    exchange_coefficient=GAMMA,
)


def _series(result, field: str = "surface_temperature"):
    times = np.array([s.time for s in result.steps])
    values = np.array([float(s.outputs[field].values) for s in result.steps])
    return times, values


# --- the null control ------------------------------------------------------


@pytest.mark.req("KER-012", "VAL-001", "VAL-004")
def test_the_unforced_configuration_does_not_move_at_all():
    """Started at the analytic equilibrium, an unforced run is exactly stationary.

    Not nearly stationary. Any drift here would be added to every step response
    measured below, so the strongest available statement is the one worth making:
    the temperature is bit-for-bit unchanged after four hundred steps.
    """
    result = surface_thermal_storage(forcing=NoForcing(), **BASE).run(400)
    _, values = _series(result)
    assert float(np.ptp(values)) == 0.0
    assert values[0] == pytest.approx(T_EQ, rel=1e-12)


@pytest.mark.req("KER-012", "VAL-001", "OCN-001")
def test_a_deep_reservoir_at_the_same_temperature_changes_nothing():
    """The disabled-coupling control for the ocean (MOD-014, CAU-007).

    Both reservoirs start at the equilibrium, so the exchange is exactly zero
    and the two configurations must agree to the last digit. A reservoir that
    perturbed a spun-up state would make every result below suspect.
    """
    control = surface_thermal_storage(forcing=NoForcing(), **BASE).run(200)
    with_ocean = reduced_ocean_heat_storage(forcing=NoForcing(), **OCEAN).run(200)
    _, plain = _series(control)
    _, coupled = _series(with_ocean)
    assert np.array_equal(plain, coupled)
    _, deep = _series(with_ocean, OCEAN_DEEP_STATE)
    assert float(np.ptp(deep)) == 0.0
    _, flux = _series(with_ocean, "deep_ocean_heat_flux")
    assert np.all(flux == 0.0)


# --- one reservoir, one timescale ------------------------------------------


@pytest.mark.req("KER-012", "TIM-007", "TIM-008", "VAL-002")
def test_the_single_reservoir_step_response_is_one_exponential():
    """The control, and the claim its declared limitation makes.

    A slab has one storage term, so its step response should be describable by
    one exponential -- and the measured e-folding time should be the linearised
    ``C / (4 sigma T^3)``. It is not exactly one exponential, because the
    emission is quartic in temperature rather than linear, and the residual
    below is that nonlinearity rather than solver error.
    """
    result = surface_thermal_storage(forcing=StepForcing(amplitude=20.0), **BASE).run(STEPS)
    times, values = _series(result)

    one = fit_exponential_modes(times, values, n_modes=1)
    assert one.relative_residual < 2e-3
    assert one.timescales[0] == pytest.approx(C_MIX / LAMBDA, rel=0.02)
    assert efolding_time(times, values) == pytest.approx(C_MIX / LAMBDA, rel=0.02)

    # Granting a second mode buys little, and the pair it returns is degenerate.
    two = fit_exponential_modes(times, values, n_modes=2)
    assert two.relative_residual > one.relative_residual / 20.0
    assert two.timescales[1] / two.timescales[0] < 3.0


# --- two reservoirs, two timescales ----------------------------------------


@pytest.mark.req("KER-012", "OCN-002", "TIM-008", "VAL-002")
def test_the_two_reservoir_step_response_has_two_measured_timescales():
    """The result Stage 5 exists for, checked against the analytic eigenvalues.

    A one-mode fit fails by two orders of magnitude where it succeeded for the
    slab; a two-mode fit succeeds, and the two timescales it finds match the
    eigenvalues of the linearised system to within the discretisation error.
    That agreement is what says the coupling is right rather than merely
    producing a curve with two bends in it.
    """
    result = reduced_ocean_heat_storage(forcing=StepForcing(amplitude=20.0), **OCEAN).run(STEPS)
    times, values = _series(result)

    one = fit_exponential_modes(times, values, n_modes=1)
    two = fit_exponential_modes(times, values, n_modes=2)
    assert one.relative_residual > 0.1
    assert two.relative_residual < one.relative_residual / 100.0

    assert two.timescales[0] == pytest.approx(ANALYTIC_FAST, rel=0.03)
    assert two.timescales[1] == pytest.approx(ANALYTIC_SLOW, rel=0.03)
    # Both inside the record, so both were measured rather than extrapolated.
    assert all(two.constrained)
    assert two.timescales[1] / two.timescales[0] > 10.0


@pytest.mark.req("TIM-008", "INF-006")
def test_the_measured_timescales_refuse_to_become_one_climate_index():
    """TIM-008 made mechanical on real output rather than on a fixture."""
    result = reduced_ocean_heat_storage(forcing=StepForcing(amplitude=20.0), **OCEAN).run(STEPS)
    times, values = _series(result)
    fit = fit_exponential_modes(times, values, n_modes=2)

    report = TimescaleReport()
    for tau, mechanism in zip(
        fit.timescales,
        ("radiative_relaxation_of_the_mixed_layer", "deep_reservoir_heat_uptake"),
        strict=True,
    ):
        report.add(
            TimescaleMeasurement(
                observable="surface_temperature",
                subsystem="mixed_layer_and_deep_reservoir",
                mechanism=mechanism,
                forcing_regime="step",
                timescale_s=tau,
                method="two-mode exponential fit",
                relative_residual=fit.relative_residual,
            )
        )
    assert report.spread > 10.0
    with pytest.raises(InvalidForScientificInferenceError, match="per mechanism"):
        report.single_index()
    assert "deep_reservoir_heat_uptake" in report.describe()


@pytest.mark.req("OCN-002", "TIM-008")
def test_the_surface_reaches_its_new_equilibrium_more_slowly_with_the_reservoir():
    """The observable consequence, stated without reference to any fit.

    At the same forcing and the same moment, the two-reservoir surface is cooler
    than the slab's -- the heat that would have warmed it is in the deep store --
    and the deep store's temperature is what says where it went.
    """
    slab = surface_thermal_storage(forcing=StepForcing(amplitude=20.0), **BASE).run(STEPS)
    ocean = reduced_ocean_heat_storage(forcing=StepForcing(amplitude=20.0), **OCEAN).run(STEPS)
    _, plain = _series(slab)
    _, coupled = _series(ocean)
    _, deep = _series(ocean, OCEAN_DEEP_STATE)

    assert coupled[-1] < plain[-1]
    assert deep[-1] > deep[0]  # the reservoir took up heat
    assert deep[-1] < coupled[-1]  # and is still colder than the layer feeding it
    # Both are heading for the same equilibrium; only the route differs.
    assert plain[-1] == pytest.approx(coupled[-1], rel=2e-3)


# --- pulse and impulse -----------------------------------------------------


@pytest.mark.req("TIM-007", "CON-009", "KER-013")
def test_a_pulse_gives_back_everything_it_delivered():
    """A stronger statement than a step response, and a different one.

    A step response can look right while energy leaks: the system settles
    wherever its own budget balances. A pulse must return the system to exactly
    the equilibrium it left, because the forcing that displaced it is gone.
    """
    pulse = PulseForcing(amplitude=40.0, onset=0.0, duration=200.0 * DT)
    result = surface_thermal_storage(forcing=pulse, **BASE).run(3000)
    _, values = _series(result)

    assert values.max() > T_EQ + 1.0  # the pulse did something
    assert values[-1] == pytest.approx(T_EQ, abs=1e-4)  # and it all came back
    for report in result.conservation_reports:
        assert report.passed, report.describe()


@pytest.mark.req("TIM-007", "CMP-009")
def test_an_impulse_delivers_the_same_energy_at_every_time_step():
    """The reason ``ImpulseForcing`` fixes energy rather than flux.

    The peak temperature anomaly is the delivered energy over the mixed layer's
    heat capacity, so if the energy were step-dependent the impulse response
    would be too, and a convergence study would be measuring the forcing.
    """
    energy = 2.0e8
    peaks = {}
    for multiplier in (2, 4, 8):
        step = multiplier * DAY
        result = surface_thermal_storage(
            forcing=ImpulseForcing(energy=energy, onset=0.0),
            master_dt=step,
            solar_step_multiple=1,
        ).run(int(2000 * 4 / multiplier))
        _, values = _series(result)
        peaks[multiplier] = values.max() - T_EQ

    reference = peaks[2]
    for multiplier, peak in peaks.items():
        assert peak == pytest.approx(reference, rel=1e-9), multiplier
    # And the peak is the absorbed energy over the heat capacity, a quarter of
    # the irradiance anomaly reaching the surface after a 0.3 albedo.
    assert reference == pytest.approx(0.7 * 0.25 * energy / C_MIX, rel=1e-3)


# --- periodic forcing ------------------------------------------------------


@pytest.mark.req("TIM-007", "TIM-008", "VAL-002")
def test_the_frequency_response_matches_the_analytic_transfer_function():
    """Charter 11.3's gain and phase, against a closed form rather than a fixture.

    Two results at once: the two-reservoir system has a lower gain than the slab
    at the same period, and it also has a *smaller* phase lag -- which looks like
    a sign error until the transfer function is written down. Both are checked
    against the analytic value, not merely against each other.
    """
    period = 4000.0 * DAY
    omega = 2.0 * math.pi / period
    steps = int(8.0 * period / DT)
    forcing = PeriodicForcing(amplitude=20.0, period=period)

    # The gain is per watt of irradiance anomaly at 1 au; the analytic transfer
    # function is per watt absorbed, and a quarter of the anomaly reaches unit
    # area of the sphere, of which the albedo reflects three tenths.
    absorbed_fraction = 0.7 * 0.25
    measurements = {}
    for name, kernel in (
        ("slab", surface_thermal_storage(forcing=forcing, **BASE)),
        ("ocean", reduced_ocean_heat_storage(forcing=forcing, **OCEAN)),
    ):
        result = kernel.run(steps)
        times, values = _series(result)
        _, anomaly = _series(result, "solar_irradiance_anomaly")
        measurements[name] = periodic_response(times, anomaly, values, period=period)

    common = dict(mixed_layer_capacity=C_MIX, deep_capacity=C_DEEP, feedback_parameter=LAMBDA)
    analytic = {
        "slab": two_reservoir_transfer_function(omega, exchange_coefficient=0.0, **common),
        "ocean": two_reservoir_transfer_function(omega, exchange_coefficient=GAMMA, **common),
    }
    for name, measured in measurements.items():
        expected_gain = abs(analytic[name]) * absorbed_fraction
        expected_lag = float(-np.angle(analytic[name]))
        assert measured.gain == pytest.approx(expected_gain, rel=0.02), name
        assert measured.phase_lag_rad == pytest.approx(expected_lag, rel=0.05), name

    assert measurements["ocean"].gain < measurements["slab"].gain
    assert measurements["ocean"].phase_lag_rad < measurements["slab"].phase_lag_rad


# --- accounting -------------------------------------------------------------


@pytest.mark.req("KER-012", "KER-013", "CON-004", "CON-009", "ACC-004")
def test_the_exchange_reconciles_and_the_budgets_close():
    """Both parties record the same transaction, and the ledger checks it.

    Putting the reservoir in its own module means the heat leaving the mixed
    layer crosses a module boundary, so it is a transfer the auditor can pair
    rather than an internal number nothing can see. Had the surface owned both
    temperatures, this test would have had nothing to check.
    """
    result = reduced_ocean_heat_storage(forcing=StepForcing(amplitude=20.0), **OCEAN).run(400)
    assert result.may_support_a_claim
    assert result.conservation_reports
    for report in result.conservation_reports:
        assert report.passed, report.describe()
    assert result.worst_relative_residual() < 1e-12

    text = "\n".join(r.describe() for r in result.conservation_reports)
    assert "surface -> ocean [energy]" in text


@pytest.mark.req("STA-001", "STA-006", "MOD-011")
def test_the_exchange_appears_as_a_labelled_contribution_to_the_surface_state():
    """The radiative and oceanic terms stay apart in the surface's own budget.

    Summed, a surface losing 3 W m^-2 to the deep ocean looks exactly like one
    emitting 3 W m^-2 more. STA-006 keeps the labels so a reader can tell which
    of the two a configuration is doing.
    """
    result = reduced_ocean_heat_storage(forcing=StepForcing(amplitude=20.0), **OCEAN).run(50)
    budget = result.steps[-1].state_budgets["surface_temperature"]
    labels = {c.label for c in budget.contributions}
    assert labels == {"net_radiative_flux", "deep_ocean_exchange"}
    exchange = next(c for c in budget.contributions if c.label == "deep_ocean_exchange")
    assert exchange.value < 0.0  # the warming surface is losing heat downward
    assert exchange.provider_module_id == "surface"


@pytest.mark.req("KER-012", "AUD-002", "AUD-004")
def test_the_forcing_and_the_reservoir_both_reach_the_configuration_hash():
    """Two runs that differ only in their experiment must not hash alike."""
    step = surface_thermal_storage(
        forcing=StepForcing(amplitude=20.0), run_id="A", **BASE
    ).build_manifest()
    pulse = surface_thermal_storage(
        forcing=PulseForcing(amplitude=20.0, duration=1.0e7), run_id="B", **BASE
    ).build_manifest()
    ocean = reduced_ocean_heat_storage(
        forcing=StepForcing(amplitude=20.0), run_id="C", **OCEAN
    ).build_manifest()

    hashes = {m.configuration_hash() for m in (step, pulse, ocean)}
    assert len(hashes) == 3
    assert step.parameters["solar"]["forcing"]["regime"] == "step"
    assert ocean.parameters["ocean"]["exchange_coefficient_W_m2_K"] == GAMMA


@pytest.mark.req("MOD-001", "TIM-006")
def test_a_forced_solar_source_does_not_claim_its_output_is_constant():
    """A contract that keeps an assumption its module abandoned is worse than none."""
    forced = surface_thermal_storage(forcing=StepForcing(amplitude=20.0), **BASE)
    solar = next(s.module for s in forced.schedule.ordered() if s.module_id == "solar")
    assumptions = " ".join(solar.contract.assumptions)
    assert "baseline total solar irradiance constant" in assumptions
    assert "step perturbation is applied" in assumptions
    # And the run is labelled a controlled experiment, not a solar record.
    assert "no statement about actual solar variability" in solar.contract.validity_domain


# --- the experiment, declared before it was run -----------------------------


@pytest.mark.req("EXP-001", "EXP-002", "HYP-003")
def test_stage_five_as_a_declared_experiment():
    """The whole stage run the way charter 14.2 asks an experiment to be run.

    The threshold below was written before the run and hashed; the evaluation
    refuses to proceed against a declaration that moved. That is the difference
    between testing a prediction and describing an outcome, and it is the reason
    this file ends here rather than with another assertion about a number.
    """
    protocol = ExperimentProtocol(
        experiment_id="EXP-STAGE5-TWO-TIMESCALES",
        question=(
            "Does a deep reservoir beneath the mixed layer produce a second, slower "
            "response timescale that a single reservoir does not have?"
        ),
        control_configuration=(
            "surface_thermal_storage: the identical configuration without the reservoir"
        ),
        manipulated_variables=("deep_ocean_reservoir_present",),
        measured_outputs=("one_mode_residual", "two_mode_residual", "timescale_ratio"),
        expected_result=(
            "for the control a one-mode fit suffices; with the reservoir it fails and a "
            "two-mode fit succeeds, with the two timescales separated by more than 10x"
        ),
        failure_threshold=(
            "the experiment fails if the two-mode fit does not improve on the one-mode "
            "fit by at least 100x, or if the two timescales differ by less than 10x"
        ),
        analysis_method="two-mode exponential fit to the step response of the surface",
        forcing_regime="step",
        assumptions=(
            "both configurations start at the analytic radiative equilibrium",
            "the reservoir is 200 m with a 2 W m^-2 K^-1 exchange, chosen so both "
            "timescales fall inside an affordable record",
        ),
    )
    protocol.freeze()

    result = reduced_ocean_heat_storage(forcing=StepForcing(amplitude=20.0), **OCEAN).run(STEPS)
    times, values = _series(result)
    one = fit_exponential_modes(times, values, n_modes=1)
    two = fit_exponential_modes(times, values, n_modes=2)
    ratio = two.timescales[1] / two.timescales[0]
    improvement = one.relative_residual / two.relative_residual

    outcome = protocol.evaluate(
        passed=improvement >= 100.0 and ratio >= 10.0,
        measured={
            "one_mode_residual": one.relative_residual,
            "two_mode_residual": two.relative_residual,
            "timescale_ratio": ratio,
        },
        summary=f"two-mode fit improved the residual {improvement:.0f}x; modes {ratio:.0f}x apart",
    )
    assert outcome.passed, outcome.describe()
    assert outcome.as_declared


# --- the same reservoir under a different surface ---------------------------


@pytest.mark.req("OCN-001", "OCN-002", "KER-010", "MOD-003")
def test_the_reservoir_couples_to_a_radiating_column_unchanged():
    """One ocean module, two surface implementations, no special casing.

    ``ColumnSurfaceModule`` and ``SlabSurfaceModule`` share nothing but an
    interface, and the deep reservoir was written against neither. Putting the
    reservoir in its own module rather than inside a surface is what makes this
    possible, and running it is the only way to know that the interface, and not
    a coincidence of implementation, is what it depends on.

    The feedback parameter here is not ``4 sigma T^3``: the atmosphere, the
    convective adjustment and the surface water are all inside the loop. So the
    analytic modes do not apply and nothing here checks against them -- what is
    checked is that the accounting still closes and the reservoir still takes up
    heat from a warming surface.
    """
    from mfespf.run.column_configurations import column_grid, column_over_reduced_ocean

    kernel = column_over_reduced_ocean(
        grid=column_grid(16),
        n_layers=8,
        areal_heat_capacity=1.0e7,
        deep_ocean_depth=100.0,
        deep_ocean_exchange=2.0,
        master_dt=3600.0,
    )
    result = kernel.run(400)
    _, surface = _series(result)
    _, deep = _series(result, OCEAN_DEEP_STATE)
    _, uptake = _series(result, "ocean_heat_uptake")

    assert result.may_support_a_claim
    for report in result.conservation_reports:
        assert report.passed, report.describe()
    assert result.worst_relative_residual() < 1e-11

    assert surface[-1] > surface[0]  # the column warms towards its equilibrium
    assert deep[-1] > deep[0]  # and the reservoir takes heat up
    assert uptake[0] == 0.0  # both started at the same temperature
    assert np.all(uptake[1:] > 0.0)
    assert deep[-1] < surface[-1]  # still colder than the layer feeding it

    # OCN-002: the outputs the ocean interface is required to provide.
    outputs = result.steps[-1].outputs
    for field_id in (
        OCEAN_DEEP_STATE,
        "ocean_heat_uptake",
        "ocean_deep_stored_heat",
        "ocean_total_stored_heat",
    ):
        assert field_id in outputs, field_id
    # The total is a sum across two owners and says so, rather than being a
    # second copy of state this module does not own (STA-001).
    total = float(outputs["ocean_total_stored_heat"].values)
    assert total > float(outputs["ocean_deep_stored_heat"].values)
    assert "not stored here" in outputs["ocean_total_stored_heat"].provenance.notes
