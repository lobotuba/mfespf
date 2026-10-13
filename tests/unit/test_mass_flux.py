"""Mass-flux convection: the plume, its operator, and the column's closure.

mass-flux-moisture-transport.md, gates A1-A5. Checked against the scheme's own
defining statements: the plume conserves mass, water and moist static energy;
the environment sinks at exactly the plume's flux; the operator redistributes
enthalpy and removes water only as precipitation; the subsiding environment warms
and dries; no CAPE, no convection.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError
from mfespf.physics.mass_flux import (
    CLOUD_WATER_CONVERSION_LAND_PER_M,
    CLOUD_WATER_CONVERSION_PER_M,
    DETRAINED_STATES,
    PRECIPITATION_EVAPORATION_EFFICIENCY,
    PRECIPITATION_EVAPORATIONS,
    apply_mass_flux,
    cloud_base_mass_flux,
    convect,
    deposit_evaporation,
    entraining_plume,
    evaporation_capacity,
    fall_and_evaporate,
)
from mfespf.physics.thermodynamics import saturation_specific_humidity
from mfespf.run.column_configurations import column_grid, reference_column

pytestmark = pytest.mark.stage(3)

G, R, CP = 9.80665, 287.04, 1004.0
LATENT = value_of("water_latent_heat_vaporization_273")


def _column(lapse=7.5e-3, rh=0.95, surface=300.0, n=30):
    """A conditionally unstable sounding: a fixed lapse rate to a 200 K floor."""
    pi = np.linspace(101325.0, 1.0e4, n + 1)
    p = 0.5 * (pi[:-1] + pi[1:])
    z = R * 255.0 / G * np.log(101325.0 / p)
    T = np.maximum(surface - lapse * z, 200.0)
    m = (pi[:-1] - pi[1:]) / G
    q = rh * saturation_specific_humidity(T, p)
    return dict(
        temperature=T,
        vapour=q * m,
        layer_mass=m,
        heat_capacity=CP * m,
        pressure=p,
        interface_pressure=pi,
        latent_heat=LATENT,
        gravity=G,
        gas_constant=R,
    )


def _h_and_q(c):
    T, m = c["temperature"], c["layer_mass"]
    dz = R * T / G * np.log(c["interface_pressure"][:-1] / c["interface_pressure"][1:])
    z = np.concatenate(([0.0], np.cumsum(dz)))[:-1] + 0.5 * dz
    q = c["vapour"] / m
    return CP * T + G * z + LATENT * q, q


@pytest.mark.req("CON-001", "EQR-002")
def test_the_plume_conserves_mass_water_and_moist_static_energy():
    """A1 and the plume's budgets: dM = E - D in every layer; the environment sinks
    at -M at every interface; sum(E h_e - D h_u) = 0 and sum(E q_e - D q_u) = sum(P)."""
    c = _column()
    plume = entraining_plume(entrainment_rate=1.0e-3, **c)
    assert plume.active and plume.cape > 0.0
    M = plume.flux
    np.testing.assert_allclose(np.diff(M), plume.entrainment - plume.detrainment, atol=1e-14)
    assert M[0] == 0.0 and M[-1] == 0.0
    h_e, q_e = _h_and_q(c)
    reached = ~np.isnan(plume.s_u)
    h_u = np.where(reached, plume.s_u + LATENT * np.nan_to_num(plume.q_u), 0.0)
    q_u = np.nan_to_num(plume.q_u)
    E, D = plume.entrainment, plume.detrainment
    scale = float(np.max(np.abs(h_e)))
    assert abs(float(np.sum(E * h_e - D * h_u))) < 1e-12 * scale * np.sum(E)
    assert float(np.sum(E * q_e - D * q_u)) == pytest.approx(
        float(np.sum(plume.condensate)), rel=1e-12
    )


@pytest.mark.req("CON-001", "CON-004", "EQR-002")
@pytest.mark.parametrize("alpha", [0.0, 0.67])
@pytest.mark.parametrize(
    ("state", "evaporation"),
    [(state, "cascade") for state in DETRAINED_STATES]
    + [("plume", evaporation) for evaporation in PRECIPITATION_EVAPORATIONS[1:4]],
)
def test_the_operator_redistributes_energy_and_removes_only_precipitation(
    alpha, state, evaporation
):
    """A2 and A3 in the operator, for every detrained state and evaporation
    placement: the temperature part sums to zero, the transport sums to zero, the
    net condensation sums to the surface precipitation, the column's water falls
    by exactly that, and no more than the evaporating fraction evaporates."""
    c = _column()
    result = convect(
        dt=3600.0,
        entrainment_rate=1.0e-3,
        evaporating_fraction=alpha,
        detrained_state=state,
        precipitation_evaporation=evaporation,
        **c,
    )
    assert result.mean_cloud_base_flux > 0.0 and result.step.substeps >= 1
    step = result.step
    typical = float(np.max(np.abs(step.energy)))
    assert abs(float(np.sum(step.energy))) < 1e-10 * typical
    assert abs(float(np.sum(step.transport))) < 1e-12 * float(np.sum(c["vapour"]))
    assert float(np.sum(step.condensation)) == pytest.approx(step.surface_precipitation, rel=1e-12)
    water_change = float(np.sum(step.transport - step.condensation))
    assert water_change == pytest.approx(-step.surface_precipitation, rel=1e-9)
    assert np.all(step.evaporation >= 0.0)
    assert float(np.sum(step.evaporation)) <= alpha * float(np.sum(step.formed)) + 1e-15


@pytest.mark.req("CON-004", "EQR-002", "VAL-001")
@pytest.mark.parametrize(
    ("state", "evaporation"),
    [(state, "cascade") for state in DETRAINED_STATES]
    + [("plume", evaporation) for evaporation in PRECIPITATION_EVAPORATIONS[1:4]],
)
def test_subsidence_dries_and_warms_a_column_without_destroying_water(state, evaporation):
    """A4: one hour of the closure's mass flux on a moist, conditionally unstable
    column with no surface fluxes. The environment between the level of free
    convection and the detrainment layer warms and its relative humidity falls,
    and water leaves only as precipitation."""
    c = _column()
    result = convect(
        dt=3600.0,
        entrainment_rate=1.0e-3,
        evaporating_fraction=0.67,
        detrained_state=state,
        precipitation_evaporation=evaporation,
        **c,
    )
    plume, step = result.plume, result.step
    C, m, p = c["heat_capacity"], c["layer_mass"], c["pressure"]
    T_new = c["temperature"] + (step.energy + LATENT * step.condensation) / C
    w_new = c["vapour"] + step.transport - step.condensation
    interior = slice(plume.lfc + 1, plume.top - 1)
    rh_old = c["vapour"] / m / saturation_specific_humidity(c["temperature"], p)
    rh_new = w_new / m / saturation_specific_humidity(T_new, p)
    assert np.mean(T_new[interior] - c["temperature"][interior]) > 0.0
    assert np.mean(rh_new[interior]) < np.mean(rh_old[interior])
    assert float(np.sum(w_new) - np.sum(c["vapour"])) == pytest.approx(
        -step.surface_precipitation, rel=1e-9
    )


@pytest.mark.req("EQR-002", "VAL-001")
def test_no_cape_means_no_convection():
    """A5: a stable, dry sounding has no dilute CAPE above the trigger, no mass flux
    and no change; and the evaporation cascade sends what a saturated layer cannot
    take to the layer below."""
    c = _column(lapse=6.5e-3, rh=0.6)
    plume = entraining_plume(entrainment_rate=1.0e-3, **c)
    flux, _ = cloud_base_mass_flux(plume, entrainment_rate=1.0e-3, evaporating_fraction=0.67, **c)
    assert flux == 0.0
    for state, evaporation in [(state, "cascade") for state in DETRAINED_STATES] + [
        ("plume", evaporation) for evaporation in PRECIPITATION_EVAPORATIONS[:4]
    ]:
        step = apply_mass_flux(
            plume,
            cloud_base_flux=flux,
            dt=3600.0,
            evaporating_fraction=0.67,
            detrained_state=state,
            precipitation_evaporation=evaporation,
            **c,
        )
        assert not np.any(step.energy) and not np.any(step.transport)
        assert not np.any(step.condensation) and not np.any(step.detrained_condensate)
    evaporated, reaching = deposit_evaporation(np.array([0.0, 0.0, 2.0]), np.array([0.5, 1.0, 0.0]))
    np.testing.assert_allclose(evaporated, [0.5, 1.0, 0.0])
    assert reaching == pytest.approx(0.5)


def _one_step(state, rh=0.95, flux=0.05, dt=600.0, alpha=0.67, evaporation="cascade"):
    """One operator step of ``state`` on the unstable sounding, and the state it leaves."""
    c = _column(rh=rh)
    plume = entraining_plume(entrainment_rate=1.0e-3, **c)
    step = apply_mass_flux(
        plume,
        cloud_base_flux=flux,
        dt=dt,
        evaporating_fraction=alpha,
        detrained_state=state,
        precipitation_evaporation=evaporation,
        **c,
    )
    T_new = c["temperature"] + (step.energy + LATENT * step.condensation) / c["heat_capacity"]
    w_new = c["vapour"] + step.transport - step.condensation
    return c, plume, step, T_new, w_new


@pytest.mark.req("CON-004", "EQR-002")
@pytest.mark.parametrize("state", DETRAINED_STATES)
def test_the_advective_split_is_the_operator(state):
    """detrainment-thermodynamic-state.md: subsidence M_(k+1)(phi_(k+1) - phi_k) and
    net detrainment D (phi_det - phi), computed independently, add up to the
    flux-form operator, because M_(k+1) - M_k = E - D."""
    _, plume, step, _, _ = _one_step(state)
    moved_vapour = step.transport - step.formed
    moved_energy = step.energy + LATENT * step.formed
    if state == "mixed_equilibrium":
        # The equilibration after the operator turns that much vapour into
        # condensate, with its latent heat into s.
        moved_vapour = moved_vapour + step.detrained_condensate
        moved_energy = moved_energy - LATENT * step.detrained_condensate
    np.testing.assert_allclose(
        step.subsidence_vapour + step.detrainment_vapour,
        moved_vapour,
        rtol=0.0,
        atol=1e-12 * float(np.max(np.abs(moved_vapour))),
    )
    np.testing.assert_allclose(
        step.subsidence_energy + step.detrainment_energy,
        moved_energy,
        rtol=0.0,
        atol=1e-10 * float(np.max(np.abs(moved_energy))),
    )
    interior = slice(plume.lfc + 1, plume.top)
    assert np.all(step.subsidence_vapour[interior] < 0.0)


@pytest.mark.req("EQR-002", "VAL-001")
def test_mixed_equilibrium_leaves_no_receiving_layer_supersaturated():
    """D1: on a saturated sounding the plume's detrainment leaves a receiving layer
    above saturation, and D1 brings each back to it at constant moist static
    energy, detraining that water as condensate. With no evaporating fraction,
    so the equilibration is seen alone: the evaporation is limited by the deficit
    before its own cooling, and that is the existing pathway, unchanged here."""
    _, _, before, _, _ = _one_step("plume", rh=1.0, alpha=0.0)
    assert float(np.sum(before.post_detrainment_excess)) > 0.0
    assert not np.any(before.detrained_condensate)
    c, plume, step, T_new, w_new = _one_step("mixed_equilibrium", rh=1.0, alpha=0.0)
    receiving = plume.detrainment > 0.0
    q_new = w_new / c["layer_mass"]
    q_s = saturation_specific_humidity(T_new, c["pressure"])
    assert np.all(q_new[receiving] <= q_s[receiving] * (1.0 + 1e-9))
    np.testing.assert_allclose(step.detrained_condensate, step.post_detrainment_excess)
    assert float(np.sum(step.detrained_condensate)) > 0.0


@pytest.mark.req("EQR-002", "VAL-001")
def test_environment_saturation_detrains_no_vapour_above_the_layers_saturation():
    """D2: detrained vapour is at most the receiving layer's saturation at its own
    temperature, and the buoyant plume's excess over it detrains as condensate."""
    c, plume, step, _, _ = _one_step("environment_saturation")
    receiving = plume.detrainment > 0.0
    scale = 0.05 * 600.0
    q_det = step.detrained_vapour[receiving] / (scale * plume.detrainment[receiving])
    q_s = saturation_specific_humidity(c["temperature"], c["pressure"])[receiving]
    assert np.all(q_det <= q_s * (1.0 + 1e-12))
    assert float(np.sum(step.detrained_condensate)) > 0.0
    _, _, plume_step, _, _ = _one_step("plume")
    assert float(np.sum(step.detrained_vapour)) < float(np.sum(plume_step.detrained_vapour))


@pytest.mark.req("EQR-002", "VAL-001")
def test_wet_bulb_capacity_leaves_no_layer_above_saturation():
    """E1: limited by the deficit before its own cooling, the cascade leaves a
    layer it evaporates into above saturation; limited by the wet-bulb capacity
    it leaves none, and still evaporates."""
    c, _, _, T_old, w_old = _one_step("plume", rh=0.9, evaporation="cascade")
    q_s = saturation_specific_humidity(T_old, c["pressure"])
    assert np.any(w_old / c["layer_mass"] > q_s * (1.0 + 1e-6))
    c, _, step, T_new, w_new = _one_step("plume", rh=0.9, evaporation="wet_bulb_cascade")
    q_s = saturation_specific_humidity(T_new, c["pressure"])
    assert np.all(w_new / c["layer_mass"] <= q_s * (1.0 + 1e-9))
    assert float(np.sum(step.evaporation)) > 0.0


@pytest.mark.req("EQR-002", "VAL-001")
def test_the_wet_bulb_capacity_is_bracketed_near_the_model_top():
    """Near 1 hPa a warm layer's q_s(T) exceeds what it could ever take: the
    capacity still lands on the wet-bulb condition, with no trial temperature
    below zero (the protocol's erratum)."""
    T = np.array([300.0, 270.0, 230.0])
    p = np.array([1.0e5, 107.0, 107.0])
    q = np.array([0.01, 1.0e-6, 1.0e-6])
    cp = np.full(3, CP)
    d = evaporation_capacity(T, q, p, cp, LATENT)
    assert np.all(d > 0.0)
    cooled = T - LATENT * d / cp
    assert np.all(cooled > 0.0)
    np.testing.assert_allclose(q + d, saturation_specific_humidity(cooled, p), rtol=1e-9)


@pytest.mark.req("EQR-002", "VAL-001")
def test_the_falling_rate_law_follows_its_closed_form():
    """E2: rain formed in a layer is not evaporated there; the layer below, with
    every other limit slack, evaporates ke (1 - q/q_s) sqrt(F) m dt; and the
    evaporating fraction caps the total."""
    m, dt = np.full(3, 1000.0), 600.0
    common = dict(
        capacity=np.full(3, 1.0e3),
        specific_humidity=np.full(3, 0.005),
        saturation=np.full(3, 0.010),
        layer_mass=m,
        dt=dt,
    )
    evaporated, reaching = fall_and_evaporate(np.array([0.0, 0.0, 1.0]), potential=1.0, **common)
    expected = PRECIPITATION_EVAPORATION_EFFICIENCY * 0.5 * np.sqrt(1.0 / dt) * m[1] * dt
    assert evaporated[2] == 0.0
    assert evaporated[1] == pytest.approx(expected, rel=1e-12)
    assert float(np.sum(evaporated)) + reaching == pytest.approx(1.0, rel=1e-12)
    capped, _ = fall_and_evaporate(np.array([0.0, 0.0, 1.0]), potential=1.0e-6, **common)
    assert float(np.sum(capped)) == pytest.approx(1.0e-6, rel=1e-12)


@pytest.mark.req("EQR-002", "VAL-001")
def test_below_cloud_base_evaporates_nothing_at_or_above_cloud_base():
    """E3: below the lowest condensing layer only."""
    _, plume, step, _, _ = _one_step("plume", rh=0.8, evaporation="below_cloud_base")
    base = int(np.flatnonzero(plume.condensate > 0.0)[0])
    assert not np.any(step.evaporation[base:])
    assert float(np.sum(step.formed)) > 0.0


@pytest.mark.req("CON-001", "EQR-002")
@pytest.mark.parametrize("c0", [CLOUD_WATER_CONVERSION_PER_M, CLOUD_WATER_CONVERSION_LAND_PER_M])
def test_the_carried_plume_conserves_its_cloud_water(c0):
    """convective-condensate-ledger.md: carrying cloud water, the plume's net
    condensation is what it converts plus what it detrains; its moist static
    energy and vapour budgets hold as before; and a slower conversion carries
    more water up."""
    c = _column()
    plume = entraining_plume(
        entrainment_rate=1.0e-3, condensate_pathway="carried", conversion_rate=c0, **c
    )
    assert plume.active and plume.liquid is not None and plume.conversion is not None
    assert plume.formation is not None and plume.cloud_evaporation is not None
    detrained = plume.detrainment * np.nan_to_num(plume.liquid)
    net = float(np.sum(plume.condensate))
    assert net == pytest.approx(float(np.sum(plume.conversion) + np.sum(detrained)), rel=1e-12)
    np.testing.assert_allclose(
        plume.condensate, plume.formation - plume.cloud_evaporation, atol=1e-15
    )
    h_e, q_e = _h_and_q(c)
    reached = ~np.isnan(plume.s_u)
    h_u = np.where(reached, plume.s_u + LATENT * np.nan_to_num(plume.q_u), 0.0)
    E, D = plume.entrainment, plume.detrainment
    scale = float(np.max(np.abs(h_e)))
    assert abs(float(np.sum(E * h_e - D * h_u))) < 1e-12 * scale * np.sum(E)
    assert float(np.sum(E * q_e - D * np.nan_to_num(plume.q_u))) == pytest.approx(net, rel=1e-12)


@pytest.mark.req("CON-001", "CON-004", "EQR-002")
@pytest.mark.parametrize(
    ("pathway", "evaporation"),
    [("immediate", evaporation) for evaporation in PRECIPITATION_EVAPORATIONS]
    + [("carried", "none"), ("carried", "rate_law_only")],
)
def test_the_condensate_ledger_closes(pathway, evaporation):
    """Every unit of condensate formed is evaporated as cloud, evaporated as
    rain, or reaches the surface: C_form = E_cloud + E_precip + P_surface. The
    redistribution still sums to zero and the water leaves only as
    precipitation."""
    c = _column()
    result = convect(
        dt=3600.0,
        entrainment_rate=1.0e-3,
        evaporating_fraction=0.67,
        precipitation_evaporation=evaporation,
        condensate_pathway=pathway,
        **c,
    )
    step = result.step
    assert result.mean_cloud_base_flux > 0.0
    formed = float(np.sum(step.formation_gross))
    destinations = (
        float(np.sum(step.cloud_evaporation))
        + float(np.sum(step.rain_evaporation))
        + step.surface_precipitation
    )
    assert destinations == pytest.approx(formed, rel=1e-12)
    assert float(np.sum(step.precipitation_formed)) == pytest.approx(
        formed - float(np.sum(step.cloud_evaporation)), rel=1e-12
    )
    typical = float(np.max(np.abs(step.energy)))
    assert abs(float(np.sum(step.energy))) < 1e-10 * typical
    assert abs(float(np.sum(step.transport))) < 1e-12 * float(np.sum(c["vapour"]))
    assert float(np.sum(step.condensation)) == pytest.approx(step.surface_precipitation, rel=1e-12)
    if evaporation == "none":
        assert not np.any(step.rain_evaporation)
    if pathway == "carried":
        assert float(np.sum(step.detrained_cloud)) > 0.0


@pytest.mark.req("CON-004", "KER-005")
@pytest.mark.parametrize("evaporation", ["none", "cascade"])
def test_the_recycling_counterfactual(evaporation):
    """convective-condensate-ledger.md: the water module's re-condensation of
    convection's evaporation, as a counterfactual. With no evaporation it is
    exactly zero. With the cascade's evaporation, limited by the deficit before
    its own cooling, some is re-condensed, and never more than was evaporated."""
    run = _warm_moist_column(
        precipitation_evaporation=evaporation, diagnose_convective_recycling=True
    ).run(4)
    assert run.may_support_a_claim, run.manifest.invalidation_reasons + run.manifest.warnings
    recycled = np.array(
        [s.outputs["convective_evaporation_recondensation_rate"].values for s in run.steps]
    )
    evaporated = np.array([s.outputs["convective_evaporation_rate"].values for s in run.steps])
    if evaporation == "none":
        assert not np.any(recycled)
    else:
        assert float(np.sum(recycled)) > 0.0
        assert np.all(recycled <= evaporated * (1.0 + 1e-9) + 1e-18)


@pytest.mark.req("CON-001", "CON-004", "KER-005", "VAL-001")
def test_the_provenance_classes_follow_the_vapour():
    """fall-path-closure-candidate.md, F5 and F6 through a coupled run: the
    classes sum to the vapour at every step, none goes negative, the milestones
    are ordered, PE = FE SE, both inventories close against their milestones,
    and the bulk state is bit-identical with the tracer off."""
    options = dict(condensate_pathway="carried", precipitation_evaporation="rate_law_only")
    run = _warm_moist_column(**options, track_water_provenance=True).run(6)
    assert run.may_support_a_claim, run.manifest.invalidation_reasons + run.manifest.warnings
    plain = _warm_moist_column(**options).run(6)
    initial = None
    convective = np.zeros((5, 30))
    large_scale = np.zeros((5, 30))
    evaporated = 0.0
    for tracked, untracked in zip(run.steps, plain.steps, strict=True):
        o = tracked.outputs
        for name in ("atmosphere_water_vapour", "atmosphere_layer_temperature"):
            assert np.array_equal(o[name].values, untracked.outputs[name].values)
        classes = np.asarray(o["water_provenance"].values, dtype=float)
        vapour = np.asarray(o["atmosphere_water_vapour"].values, dtype=float)
        assert np.all(np.abs(classes.sum(axis=0) - vapour) <= 1e-12 * vapour)
        assert np.all(classes >= -1e-12 * vapour)
        assert not np.any(classes[:, vapour == 0.0])
        convective += np.asarray(o["convective_provenance_milestones"].values) * 1.0
        large_scale += np.asarray(o["large_scale_provenance_milestones"].values) * 1.0
        evaporated += float(o["surface_evaporation_flux"].values) * 3600.0
        total = convective + large_scale
        c_unique, c_precip, p_surface = (float(np.sum(total[i])) for i in (0, 2, 4))
        assert p_surface <= c_precip * (1 + 1e-12) and c_precip <= c_unique * (1 + 1e-12)
    assert c_unique > 0.0 and c_precip > 0.0 and p_surface > 0.0
    FE, SE, PE = c_precip / c_unique, p_surface / c_precip, p_surface / c_unique
    assert PE == pytest.approx(FE * SE, rel=1e-14)
    initial = np.asarray(_warm_moist_column(**options).initial_state["atmosphere_water_vapour"])
    final = np.asarray(run.steps[-1].outputs["water_provenance"].values, dtype=float)
    throughput = c_unique + float(np.sum(large_scale[1])) + evaporated
    # Convective-born water enters only by first condensation and leaves only
    # by reaching the surface; never-condensed water enters only by surface
    # evaporation and leaves only by condensing.
    assert float(np.sum(final[1:])) == pytest.approx(c_unique - p_surface, abs=1e-12 * throughput)
    assert float(np.sum(final[0]) - np.sum(initial)) == pytest.approx(
        evaporated - c_unique - float(np.sum(large_scale[1])), abs=1e-12 * throughput
    )


def _warm_moist_column(**overrides):
    n = 30
    options = dict(
        grid=column_grid(16),
        n_layers=n,
        master_dt=3600.0,
        measured_solar_spectrum=True,
        ozone=True,
        solar_vapour_absorption=True,
        eddy_diffusivity=0.39,
        condense_after_convection=True,
        end_of_step_saturation=True,
        convective_target="dry_adiabat",
        mass_flux_convection=True,
    )
    options.update(overrides)
    probe = reference_column(
        **{
            **options,
            "mass_flux_convection": False,
            "convective_target": "entraining_plume",
            "precipitation_evaporation": "cascade",
            "diagnose_convective_recycling": False,
            "condensate_pathway": "immediate",
            "track_water_provenance": False,
        }
    )
    atmosphere = next(s.module for s in probe.schedule.ordered() if s.module_id == "atmosphere")
    pressure = np.asarray(atmosphere.layer_pressure, dtype=float)
    z = R * 255.0 / G * np.log(101325.0 / pressure)
    T = np.maximum(300.0 - 7.5e-3 * z, 200.0)
    kernel = reference_column(
        **options, initial_air_temperature=T, initial_surface_temperature=301.0
    )
    mass = np.asarray(atmosphere.layer_mass, dtype=float)
    vapour = 0.95 * saturation_specific_humidity(T, pressure) * mass
    vapour[pressure < 1.0e4] = 4.0e-6 * mass[pressure < 1.0e4]
    return replace(
        kernel, initial_state={**kernel.initial_state, "atmosphere_water_vapour": vapour}
    )


@pytest.mark.req("CON-001", "CON-004", "KER-005")
def test_a_convecting_column_closes_its_water_and_energy():
    """A2 and A3 coupled: a warm, moist, unstable column convects on its first step,
    and the run's water and energy budgets reconcile, the plume's condensation booked
    by the water module as precipitation and latent heat."""
    run = _warm_moist_column().run(4)
    assert run.may_support_a_claim, run.manifest.invalidation_reasons + run.manifest.warnings
    assert run.worst_relative_residual() < 1e-10
    first = run.steps[0].outputs
    assert float(first["cloud_base_mass_flux"].values) > 0.0
    assert float(first["convective_surface_precipitation"].values) > 0.0


@pytest.mark.req("SPC-001", "EQR-002")
def test_the_option_refuses_the_moist_adjustment_it_replaces():
    with pytest.raises(ConfigurationError, match="replaces the entraining-plume adjustment"):
        reference_column(
            grid=column_grid(16),
            n_layers=10,
            condense_after_convection=True,
            mass_flux_convection=True,
        )
