"""Stage 2 of the validation ladder: the radiating atmospheric column.

Charter section 19 stage 2 and requirement KER-012, which asks for a transparent
atmosphere and an absorbing-and-emitting one.

The regression test is the important one. The transparent configuration is the
same code as the absorbing one with the opacity set to zero, so when it
reproduces the Stage 1 airless equilibrium it is evidence that the column
machinery -- radiation module, layer state, four new transfers, a longer
splitting order -- introduced no bias. A separately written special case would
prove nothing.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.constants import value_of
from mfespf.coupling.transfers import ConservedQuantity
from mfespf.physics.atmosphere import ATMOSPHERE_TEMPERATURE_STATE
from mfespf.physics.radiation import equilibrium_temperature
from mfespf.run.column_configurations import (
    absorbing_emitting_atmosphere,
    column_grid,
    single_grey_layer,
    transparent_atmosphere,
    window_atmosphere,
)
from mfespf.run.configurations import airless_radiating_body

pytestmark = pytest.mark.stage(2)

TSI = value_of("total_solar_irradiance_reference")
C_SMALL = 1.0e6  # a thin slab so equilibrium is reached in a tractable run
GRID = column_grid(16)


def _final(kernel, steps: int) -> tuple[float, np.ndarray, object]:
    result = kernel.run(steps)
    last = result.steps[-1].outputs
    return (
        last["surface_temperature"].scalar(),
        np.array(last[ATMOSPHERE_TEMPERATURE_STATE].values),
        result,
    )


@pytest.mark.req("KER-012", "VAL-004", "SYS-001")
def test_transparent_column_recovers_the_airless_result():
    """The column null control, compared against Stage 1 under identical optics.

    Both surfaces are grey with albedo 0.3, so both have infrared emissivity 0.7
    by Kirchhoff's law. That equality is what makes the comparison meaningful:
    the airless configuration's usual pairing of a 0.3 albedo with unit
    emissivity is only self-consistent if the reflectance is confined to short
    wavelengths, and the column model -- which derives emissivity from albedo band
    by band -- will not let that inconsistency pass silently.
    """
    albedo = 0.3
    column_temperature, layers, result = _final(
        transparent_atmosphere(
            grid=GRID,
            albedo=albedo,
            areal_heat_capacity=C_SMALL,
            master_dt=600.0,
            initial_surface_temperature=100.0,
        ),
        9000,
    )
    airless = airless_radiating_body(
        albedo=albedo,
        emissivity=1.0 - albedo,
        areal_heat_capacity=C_SMALL,
        master_dt=600.0,
        initial_temperature=100.0,
        solar_step_multiple=1,
    )
    airless_temperature = airless.run(9000).steps[-1].outputs["surface_temperature"].scalar()

    assert result.may_support_a_claim
    assert column_temperature == pytest.approx(airless_temperature, rel=1e-6)
    # A transparent layer neither absorbs nor emits, so it never changes.
    assert np.allclose(layers, layers[0])
    assert layers[0] == pytest.approx(250.0, abs=1e-12)


@pytest.mark.req("KER-012", "SPC-009")
def test_the_shortwave_only_albedo_differs_from_a_grey_one_by_the_spectral_overlap():
    """Two readings of "albedo 0.3" give measurably different equilibria.

    A surface reflective only below 4 microns is black in the infrared and
    reaches the textbook 254.6 K; a genuinely grey one has emissivity 0.7 and
    reaches 278 K. Both are legitimate configurations and the framework records
    which was used -- SPC-009's point, in the smallest possible setting.
    """
    grey, _, _ = _final(
        transparent_atmosphere(
            grid=GRID,
            albedo=0.3,
            areal_heat_capacity=C_SMALL,
            master_dt=600.0,
            initial_surface_temperature=200.0,
        ),
        9000,
    )
    shortwave_only, _, _ = _final(
        transparent_atmosphere(
            grid=GRID,
            areal_heat_capacity=C_SMALL,
            master_dt=600.0,
            initial_surface_temperature=200.0,
        ),
        9000,
    )
    textbook = equilibrium_temperature(TSI, 0.3, 1.0)

    assert grey > shortwave_only + 20.0
    # The shortwave-only surface lands near the textbook value; the residual is
    # the solar spectrum's own tail above the 4 micron cutoff, about 0.1 percent.
    assert shortwave_only == pytest.approx(textbook, rel=2e-3)
    assert shortwave_only != pytest.approx(textbook, rel=1e-5)


@pytest.mark.req("KER-012", "CON-009", "ACC-004")
def test_column_energy_budget_closes():
    """Six modules, four transfers, and a budget that still closes to round-off."""
    kernel = absorbing_emitting_atmosphere(
        kappa=1.0e-3, grid=GRID, areal_heat_capacity=C_SMALL, master_dt=1800.0
    )
    result = kernel.run(200)

    reports = result.conservation_reports
    assert reports
    for report in reports:
        assert report.passed, report.describe()
        assert {b.module_id for b in report.budgets} == {
            "orbit",
            "atmosphere",
            "surface",
            "space",
        }
    assert result.worst_relative_residual() < 1e-11


@pytest.mark.req("CON-001", "CON-009")
def test_column_system_energy_is_conserved_across_the_boundary():
    """Storage gained by every store must equal what crossed the top."""
    kernel = absorbing_emitting_atmosphere(
        kappa=1.0e-3, grid=GRID, areal_heat_capacity=C_SMALL, master_dt=1800.0
    )
    result = kernel.run(100)

    supplied = 0.0
    stored = 0.0
    for report in result.conservation_reports:
        for budget in report.budgets_for(ConservedQuantity.ENERGY):
            if budget.module_id == "orbit":
                supplied += budget.net_physical
            else:
                stored += budget.stored_change
    assert stored == pytest.approx(supplied, rel=1e-10)


@pytest.mark.req("ATM-003", "CON-001")
def test_column_mass_is_hydrostatic():
    """Total column mass must be exactly (p_s - p_top) / g."""
    kernel = absorbing_emitting_atmosphere(grid=GRID, n_layers=12)
    atmosphere = next(s.module for s in kernel.schedule.ordered() if s.module_id == "atmosphere")
    expected = (atmosphere.surface_pressure - atmosphere.top_pressure) / atmosphere.gravity
    assert atmosphere.column_mass == pytest.approx(expected, rel=1e-14)
    assert atmosphere.layer_mass.size == 12
    assert np.all(atmosphere.layer_mass > 0.0)


@pytest.mark.req("KER-012", "CLM-003")
def test_the_greenhouse_response_has_the_declared_sign_and_magnitude():
    """Adding infrared opacity must warm the surface, and by roughly the analytic amount.

    The expected sign is declared before the run, as HYP-003 requires. The
    magnitude is compared against the closed form for a single grey layer, and
    the agreement is deliberately loose: the configuration violates two of the
    identity's assumptions -- the layer absorbs the 0.9 percent of sunlight that
    lies above 4 microns, and its emission below 4 microns is not quite zero --
    so exact agreement would indicate a coincidence rather than a correct model.
    The exact test lives at the solver level, free of those effects.
    """
    baseline, _, _ = _final(
        single_grey_layer(
            0.0001,
            grid=GRID,
            areal_heat_capacity=C_SMALL,
            master_dt=1200.0,
            initial_surface_temperature=250.0,
            initial_air_temperature=230.0,
        ),
        5000,
    )
    warmed, layer_temperatures, _ = _final(
        single_grey_layer(
            0.75,
            grid=GRID,
            areal_heat_capacity=C_SMALL,
            master_dt=1200.0,
            initial_surface_temperature=250.0,
            initial_air_temperature=230.0,
        ),
        5000,
    )

    assert warmed > baseline  # the declared sign
    reference = equilibrium_temperature(TSI, 0.3, 1.0)
    analytic = reference * (2.0 / (2.0 - 0.75)) ** 0.25
    assert warmed == pytest.approx(analytic, rel=0.01)
    # The layer sits colder than the surface, as the identity requires.
    assert float(layer_temperatures[0]) < warmed


@pytest.mark.req("VAL-004", "MOD-014", "KER-012")
def test_zero_insolation_column_only_cools():
    """The column null control: with the star off, nothing may warm."""
    kernel = absorbing_emitting_atmosphere(
        kappa=1.0e-3,
        grid=GRID,
        total_irradiance=0.0,
        areal_heat_capacity=C_SMALL,
        master_dt=1800.0,
        initial_surface_temperature=280.0,
        initial_air_temperature=280.0,
    )
    result = kernel.run(200)

    surface = [s.outputs["surface_temperature"].scalar() for s in result.steps]
    column = [float(np.sum(s.outputs[ATMOSPHERE_TEMPERATURE_STATE].values)) for s in result.steps]
    assert surface[-1] < surface[0]
    assert column[-1] < column[0]
    for report in result.conservation_reports:
        assert report.passed, report.describe()


@pytest.mark.req("STA-006", "MOD-011", "STA-003")
def test_radiation_supplies_the_labeled_tendency_the_column_commits():
    """A2.3: another module contributes, the owner commits, and the label survives."""
    kernel = absorbing_emitting_atmosphere(
        kappa=1.0e-3, grid=GRID, areal_heat_capacity=C_SMALL, master_dt=1800.0
    )
    result = kernel.run(5)

    for step in result.steps:
        budget = step.state_budgets[ATMOSPHERE_TEMPERATURE_STATE]
        assert [c.label for c in budget.contributions] == ["radiative_flux_divergence"]
        assert [c.provider_module_id for c in budget.contributions] == ["radiation"]
        assert np.allclose(budget.residual(), 0.0, atol=1e-12)


@pytest.mark.req("SPC-009", "KER-009")
def test_a_spectral_window_changes_the_equilibrium_at_equal_opacity():
    """Same integrated opacity, different distribution, different surface temperature."""
    windowed, _, _ = _final(
        window_atmosphere(
            kappa_opaque=1.0e-3,
            grid=GRID,
            areal_heat_capacity=C_SMALL,
            master_dt=1200.0,
            initial_surface_temperature=260.0,
            initial_air_temperature=250.0,
        ),
        4000,
    )
    kernel = window_atmosphere(kappa_opaque=1.0e-3, grid=GRID)
    radiation = next(s.module for s in kernel.schedule.ordered() if s.module_id == "radiation")
    equivalent_grey = float(np.mean(radiation.absorption_coefficient))

    grey, _, _ = _final(
        absorbing_emitting_atmosphere(
            kappa=equivalent_grey,
            grid=GRID,
            areal_heat_capacity=C_SMALL,
            master_dt=1200.0,
            initial_surface_temperature=260.0,
            initial_air_temperature=250.0,
        ),
        4000,
    )
    assert abs(windowed - grey) > 0.5  # K, far beyond integration noise


@pytest.mark.req("MOD-008")
def test_column_checkpoint_and_restore_round_trip():
    kernel = absorbing_emitting_atmosphere(
        kappa=1.0e-3, grid=GRID, areal_heat_capacity=C_SMALL, master_dt=1800.0
    )
    coupler, _ = kernel.prepare(n_steps=20)
    coupler.run(10)

    atmosphere = next(s.module for s in kernel.schedule.ordered() if s.module_id == "atmosphere")
    saved = atmosphere.checkpoint()
    coupler.run(5)
    assert not np.array_equal(atmosphere.temperature, saved[ATMOSPHERE_TEMPERATURE_STATE])

    atmosphere.restore(saved)
    assert np.array_equal(atmosphere.temperature, saved[ATMOSPHERE_TEMPERATURE_STATE])
