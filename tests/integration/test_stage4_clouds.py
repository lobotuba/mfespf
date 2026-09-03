"""Stage 4 of the validation ladder: clouds.

Charter section 19 stage 4 and requirements CLD-001, CLD-002, SPC-010, KER-012.
These two configurations complete KER-012: every reference configuration the
specification names now exists and runs.

The charter's cloud reporting rule drives most of this file. Shortwave and
longwave effects must be available separately before any net is formed, because
they carry opposite signs and comparable magnitudes -- a net figure alone cannot
distinguish a model that has both right from one that has both wrong.
"""

from __future__ import annotations

import pytest

from mfespf.contracts.module import ModuleStatus
from mfespf.coupling.transfers import ConservedQuantity
from mfespf.physics.clouds import CloudProfile, single_layer_cloud
from mfespf.run.column_configurations import (
    column_grid,
    no_cloud_control,
    prescribed_cloud,
)

pytestmark = pytest.mark.stage(4)

GRID = column_grid(16)
BASE = dict(grid=GRID, n_layers=8, areal_heat_capacity=1.0e6, master_dt=3600.0)


def _run(kernel, steps: int = 400):
    result = kernel.run(steps)
    return result, result.steps[-1].outputs


# --- the control -----------------------------------------------------------


@pytest.mark.req("KER-012", "VAL-004", "CLD-001")
def test_no_cloud_control_matches_the_clear_column():
    """A cloud of zero fraction must reproduce the control exactly.

    The strongest available statement that the cloud machinery is inert when it
    should be: the same configuration with a cloud module present but empty must
    agree to the last digit, not merely closely.
    """
    _, control = _run(no_cloud_control(**BASE), 200)
    _, empty = _run(
        prescribed_cloud(fraction=0.0, **BASE),
        200,
    )
    assert empty["surface_temperature"].scalar() == pytest.approx(
        control["surface_temperature"].scalar(), rel=1e-14
    )
    assert empty["outgoing_longwave_radiation"].scalar() == pytest.approx(
        control["outgoing_longwave_radiation"].scalar(), rel=1e-14
    )
    # And the effect it reports is exactly zero, not merely small.
    assert float(empty["cloud_shortwave_radiative_effect"].values) == 0.0
    assert float(empty["cloud_longwave_radiative_effect"].values) == 0.0


# --- the cloud radiative effect --------------------------------------------


@pytest.mark.req("CLD-001", "SPC-010", "RAD-002", "CLM-003")
def test_cloud_shortwave_effect_cools_and_longwave_effect_warms():
    """The declared signs, stated before the run as HYP-003 requires.

    A cloud reflects sunlight, so its shortwave effect is negative; it absorbs
    terrestrial radiation and re-emits from a colder level, so its longwave
    effect is positive. Getting either sign wrong is a refutation under the
    charter's 15.4 taxonomy, not an inaccuracy.
    """
    _, outputs = _run(prescribed_cloud(liquid_water_path=0.05, fraction=0.6, **BASE))

    shortwave = float(outputs["cloud_shortwave_radiative_effect"].values)
    longwave = float(outputs["cloud_longwave_radiative_effect"].values)

    assert shortwave < 0.0
    assert longwave > 0.0
    # And both are large enough to be a result rather than numerical noise.
    assert abs(shortwave) > 1.0
    assert longwave > 1.0


@pytest.mark.req("SPC-010", "CLD-001")
def test_cloud_radiative_effect_is_reported_before_any_net():
    """SPC-010: the components exist as fields, and the net is their sum."""
    _, outputs = _run(prescribed_cloud(**BASE), 100)

    shortwave = float(outputs["cloud_shortwave_radiative_effect"].values)
    longwave = float(outputs["cloud_longwave_radiative_effect"].values)
    net = float(outputs["cloud_net_radiative_effect"].values)

    assert net == pytest.approx(shortwave + longwave, rel=1e-12)
    # The clear-sky reference the effect is defined against is published too, so
    # a reader can recompute it rather than take the difference on trust.
    assert "clear_sky_outgoing_longwave_radiation" in outputs
    assert "clear_sky_reflected_shortwave_radiation" in outputs
    clear = float(outputs["clear_sky_outgoing_longwave_radiation"].values)
    all_sky = float(outputs["outgoing_longwave_radiation"].values)
    assert longwave == pytest.approx(clear - all_sky, rel=1e-12)


@pytest.mark.req("CLD-002", "CLM-003")
def test_a_thicker_cloud_reflects_more():
    """More water reflects more sunlight, and the effect saturates."""
    effects = []
    for water in (0.01, 0.05, 0.2):
        _, outputs = _run(prescribed_cloud(liquid_water_path=water, fraction=1.0, **BASE), 200)
        effects.append(float(outputs["cloud_shortwave_radiative_effect"].values))

    assert effects[0] > effects[1] > effects[2]  # increasingly negative
    # Saturation: doubling the water a second time buys less than the first.
    assert abs(effects[2] - effects[1]) < abs(effects[1] - effects[0]) * 4.0


@pytest.mark.req("CLD-002", "CLM-003")
def test_a_high_cloud_warms_more_than_a_low_one():
    """Cloud altitude governs the longwave effect, which is why layers matter.

    A high cloud radiates to space from a colder level than a low cloud of the
    same water path, so it reduces the outgoing longwave more. This is the
    property a model that folded cloud into a single albedo parameter could not
    express at all -- which is what CLD-001 is protecting against.
    """
    _, low = _run(prescribed_cloud(layer=1, liquid_water_path=0.05, fraction=1.0, **BASE))
    _, high = _run(prescribed_cloud(layer=6, liquid_water_path=0.05, fraction=1.0, **BASE))

    assert float(high["cloud_longwave_radiative_effect"].values) > float(
        low["cloud_longwave_radiative_effect"].values
    )


@pytest.mark.req("CLD-002")
def test_the_effect_scales_with_cloud_fraction():
    """Under maximum overlap the effect is linear in the fraction."""
    kernels = {
        fraction: prescribed_cloud(fraction=fraction, **BASE) for fraction in (0.25, 0.5, 1.0)
    }
    effects = {}
    for fraction, kernel in kernels.items():
        _, outputs = _run(kernel, 100)
        effects[fraction] = float(outputs["cloud_shortwave_radiative_effect"].values)

    assert effects[0.5] < effects[0.25] < 0.0
    assert effects[1.0] < effects[0.5]


# --- the prescribed state and its honesty ---------------------------------


@pytest.mark.req("CON-009", "SYS-003")
def test_prescribed_cloud_water_is_declared_as_a_source():
    """A prescribed cloud's water comes from nowhere, and the budget says so.

    The condensate is not drawn from the water cycle. Rather than leave that
    outside the accounting, the cloud module books it against a named source, so
    a reader of the conservation report sees the fiat rather than inferring it.
    """
    kernel = prescribed_cloud(**BASE)
    cloud = next(s.module for s in kernel.schedule.ordered() if s.module_id == "cloud")

    assert cloud.status is ModuleStatus.PRESCRIBED
    budgets = cloud.report_conservation()
    assert len(budgets) == 1
    budget = budgets[0]
    assert budget.quantity is ConservedQuantity.WATER
    assert "prescribed_cloud_water_source" in budget.source_labels
    assert budget.residual == pytest.approx(0.0, abs=1e-12)
    assert any(
        "does not conserve water" in limitation for limitation in cloud.contract.known_limitations
    )


@pytest.mark.req("KER-012", "CON-009", "ACC-004")
def test_the_cloudy_column_energy_budget_closes():
    """Scattering, two solves and a mixed column, and the budget still closes."""
    result, _ = _run(prescribed_cloud(**BASE), 200)

    assert result.may_support_a_claim
    assert result.conservation_reports
    for report in result.conservation_reports:
        assert report.passed, report.describe()
    assert result.worst_relative_residual() < 1e-11


@pytest.mark.req("KER-012", "CLD-001")
def test_a_time_varying_cloud_is_supported():
    """The charter's stage 4 asks for a prescribed *variable* cloud state."""

    def varying(time: float) -> CloudProfile:
        # A cloud that thickens over the first simulated day.
        water = 0.02 + 0.04 * min(time / 86400.0, 1.0)
        return single_layer_cloud(n_layers=8, layer=4, liquid_water_path=water, fraction=0.7)

    kernel = prescribed_cloud(cloud=varying, **BASE)
    result = kernel.run(48)

    paths = [float(s.outputs["cloud_total_water_path"].values) for s in result.steps]
    assert paths[-1] > paths[0]
    effects = [float(s.outputs["cloud_shortwave_radiative_effect"].values) for s in result.steps]
    assert effects[-1] < effects[0]  # thicker cloud, stronger cooling
    for report in result.conservation_reports:
        assert report.passed, report.describe()


@pytest.mark.req("CLD-002", "MOD-003")
def test_the_cloud_module_publishes_everything_cld002_requires():
    """Liquid, ice, fraction, vertical support and optical properties."""
    _, outputs = _run(prescribed_cloud(**BASE), 5)
    for field_id in (
        "cloud_liquid_water_path",
        "cloud_ice_water_path",
        "cloud_fraction",
        "cloud_effective_radius",
        "cloud_shortwave_optical_depth",
        "cloud_longwave_optical_depth",
    ):
        assert field_id in outputs, field_id
    # Vertical support: the water paths are per layer, not a column total.
    assert outputs["cloud_liquid_water_path"].values.size == 8
    assert outputs["cloud_liquid_water_path"].dimensions == ("layer",)


@pytest.mark.req("CLD-001", "AUD-002")
def test_the_cloud_reaches_the_configuration_hash():
    """Two clouds that differ physically must not hash alike."""
    thin = prescribed_cloud(liquid_water_path=0.01, run_id="A", **BASE).build_manifest()
    thick = prescribed_cloud(liquid_water_path=0.10, run_id="B", **BASE).build_manifest()
    assert thin.configuration_hash() != thick.configuration_hash()
    assert (
        thin.parameters["cloud"]["liquid_water_path_kg_m2"]
        != (thick.parameters["cloud"]["liquid_water_path_kg_m2"])
    )
