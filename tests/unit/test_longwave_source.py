"""The longwave source function: correlated-k emission and the linear-in-tau source.

longwave-source-function.md. T1-T7 are the frozen unit tests; P2 checks that
this model's solver, with both, reproduces RRTMG's own transfer on AER's four
standard atmospheres to P1's limits.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError
from mfespf.physics.radiation import band_exitance
from mfespf.physics.radiative_transfer import (
    PlanckFractions,
    interface_temperatures,
    linear_in_tau_weight,
    solve_column_radiation,
)
from mfespf.run.column_configurations import column_grid

pytestmark = pytest.mark.stage(6)


def _column(n_layers=6, *, scatter=False):
    grid = column_grid(16)
    rng = np.random.default_rng(3)
    tau = rng.uniform(0.0, 3.0, (n_layers, grid.n_bins))
    return dict(
        grid=grid,
        incident_spectral_irradiance=np.zeros(grid.n_bins),
        layer_temperature=np.linspace(290.0, 220.0, n_layers),
        layer_mass=np.full(n_layers, 1.0e3),
        absorption_coefficient=np.zeros(grid.n_bins),
        surface_temperature=295.0,
        surface_albedo=np.full(grid.n_bins, 0.05),
        absorption_optical_depth=tau,
        **({"scattering_optical_depth": 0.5 * tau, "asymmetry": 0.8} if scatter else {}),
    )


@pytest.mark.req("RAD-002", "EQR-002")
def test_both_options_off_leave_the_solver_unchanged():
    """T1: the defaults are the present solver, to the bit."""
    base = solve_column_radiation(**_column(scatter=True))
    explicit = solve_column_radiation(
        **_column(scatter=True), planck_fractions=None, interface_temperature=None
    )
    assert np.array_equal(base.upward, explicit.upward)
    assert np.array_equal(base.downward, explicit.downward)


@pytest.mark.req("RAD-002", "EQR-002")
def test_an_isothermal_column_is_unchanged_by_the_linear_source():
    """T2: with every interface at the layers' temperature, the linear source is the
    isothermal one."""
    case = _column(scatter=True)
    case["layer_temperature"] = np.full(6, 255.0)
    base = solve_column_radiation(**case)
    linear = solve_column_radiation(**case, interface_temperature=np.full(7, 255.0))
    np.testing.assert_allclose(linear.upward, base.upward, rtol=1e-12, atol=1e-300)
    np.testing.assert_allclose(linear.downward, base.downward, rtol=1e-12, atol=1e-300)


@pytest.mark.req("RAD-002", "EQR-002")
def test_the_linear_weight_has_its_limits():
    """T3: f(x)/x -> 1/6 as x -> 0, f -> 1 as x -> infinity, continuous at 0.06."""
    small = np.array([1.0e-8, 1.0e-5])
    np.testing.assert_allclose(linear_in_tau_weight(small) / small, 1.0 / 6.0, rtol=1e-6)
    assert linear_in_tau_weight(np.array([1.0e4]))[0] == pytest.approx(1.0, abs=1e-3)
    below, above = linear_in_tau_weight(np.array([0.06, 0.06 + 1e-12]))
    assert abs(above - below) < 1e-4
    x = np.linspace(0.0, 50.0, 2001)
    assert np.all(np.diff(linear_in_tau_weight(x)) >= 0.0)


@pytest.mark.req("RAD-002", "EQR-002")
def test_an_opaque_layer_emits_at_its_edges():
    """T4 (second erratum): a 280 K top and a 300 K bottom around a 290 K layer. At
    D tau = 50 each edge gets B_lay + f(50)(B_edge - B_lay), f(50) = 0.96 -- an opaque
    layer emits from about one optical depth inside its edge -- and at D tau = 1e9 the
    totals are the edges' Planck fluxes."""
    grid = column_grid(16)
    D = 1.66

    def solve(x):
        return solve_column_radiation(
            grid=grid,
            incident_spectral_irradiance=np.zeros(grid.n_bins),
            layer_temperature=np.array([290.0]),
            layer_mass=np.array([1.0e3]),
            absorption_coefficient=np.zeros(grid.n_bins),
            surface_temperature=300.0,
            surface_albedo=np.zeros(grid.n_bins),
            absorption_optical_depth=np.full((1, grid.n_bins), x / D),
            diffusivity_factor=D,
            interface_temperature=np.array([300.0, 280.0]),
        )

    lay, top, bottom = (band_exitance(grid, t) for t in (290.0, 280.0, 300.0))
    f = 1.0 - 2.0 * (1.0 / 50.0 - np.exp(-50.0) / (1.0 - np.exp(-50.0)))
    assert f == pytest.approx(0.96, abs=1e-12)
    thick = solve(50.0)
    np.testing.assert_allclose(thick.upward[1], lay + f * (top - lay), rtol=1e-6)
    # The surface's own emission is all absorbed; what reaches it is the layer's.
    np.testing.assert_allclose(thick.downward[0], lay + f * (bottom - lay), rtol=1e-6)
    limit = solve(1.0e9)
    assert np.sum(limit.upward[1]) == pytest.approx(np.sum(top), rel=1e-6)
    assert np.sum(limit.downward[0]) == pytest.approx(np.sum(bottom), rel=1e-6)


def _rrtmg_grid():
    from mfespf.physics.rrtmg_longwave import (
        _NGS,
        rrtmg_fixed_fractions,
        rrtmg_longwave_grid,
    )

    grid, gp = rrtmg_longwave_grid(column_grid(16))
    lw = gp >= 0
    band = np.full(gp.size, -1)
    band[lw] = np.searchsorted(_NGS, gp[lw], side="right") - 1
    frac = np.zeros(gp.size)
    frac[lw] = rrtmg_fixed_fractions()[gp[lw]]
    return grid, band, frac


@pytest.mark.req("RAD-002", "EQR-002")
def test_fixed_fractions_close_an_isothermal_column():
    """T5: layer-independent fractions, an isothermal column over a black surface at
    the same temperature: each band's OLR is that band's Planck flux."""
    grid, band, frac = _rrtmg_grid()
    n = 5
    T = 250.0
    fractions = PlanckFractions(band=band, layer=np.tile(frac, (n, 1)), surface=frac)
    sol = solve_column_radiation(
        grid=grid,
        incident_spectral_irradiance=np.zeros(grid.n_bins),
        layer_temperature=np.full(n, T),
        layer_mass=np.full(n, 1.0e3),
        absorption_coefficient=np.zeros(grid.n_bins),
        surface_temperature=T,
        surface_albedo=np.zeros(grid.n_bins),
        absorption_optical_depth=np.random.default_rng(1).uniform(0, 5, (n, grid.n_bins)),
        planck_fractions=fractions,
        interface_temperature=np.full(n + 1, T),
    )
    own = band_exitance(grid, T)
    for b in range(int(band.max()) + 1):
        members = band == b
        assert np.sum(sol.upward[-1][members]) == pytest.approx(np.sum(own[members]), rel=1e-9)


@pytest.mark.req("RAD-002", "EQR-002")
def test_the_source_options_refuse_what_they_cannot_honour():
    """T6: fractions that do not close, interfaces of the wrong length, and correlated-k
    emission without g-points."""
    grid, band, frac = _rrtmg_grid()
    case = dict(
        grid=grid,
        incident_spectral_irradiance=np.zeros(grid.n_bins),
        layer_temperature=np.full(3, 250.0),
        layer_mass=np.full(3, 1.0e3),
        absorption_coefficient=np.zeros(grid.n_bins),
        surface_temperature=260.0,
        surface_albedo=np.zeros(grid.n_bins),
    )
    open_ = PlanckFractions(band=band, layer=np.tile(frac, (3, 1)), surface=1.001 * frac)
    with pytest.raises(ConfigurationError, match="sum to one"):
        solve_column_radiation(**case, planck_fractions=open_)
    with pytest.raises(ConfigurationError, match="interfaces"):
        solve_column_radiation(**case, interface_temperature=np.full(3, 250.0))
    from mfespf.run.column_configurations import reference_column

    with pytest.raises(ConfigurationError, match="g-points"):
        reference_column(grid=column_grid(16), n_layers=10, k_distribution_source=True)


@pytest.mark.req("RAD-002", "CON-006")
def test_the_column_closes_its_energy_with_both_options():
    """T7: the module with both options on; its solve is the direct one with CAM's
    interfaces and RRTMG's fractions, the surface emits its greybody total, and the
    column's radiative energy closes to round-off."""
    from mfespf.physics.rrtmg_longwave import rrtmg_optical_depth, rrtmg_planck_fractions
    from mfespf.run.column_configurations import reference_column

    kernel = reference_column(
        grid=column_grid(16),
        n_layers=10,
        master_dt=3600.0,
        measured_solar_spectrum=True,
        ozone=True,
        solar_vapour_absorption=True,
        rrtmg_longwave=True,
        k_distribution_source=True,
        linear_in_tau_source=True,
    )
    kernel.run(2)
    rad = next(s.module for s in kernel.schedule.ordered() if s.module_id == "radiation")

    def value(name, units):
        return np.asarray(rad.input_field(name).require_units(units).values, dtype=float)

    temperature = value("atmosphere_layer_temperature", "K")
    mass = value("atmosphere_layer_mass", "kg m^-2")
    pressure = value("atmosphere_layer_pressure", "Pa")
    surface = float(value("surface_temperature", "K"))
    tau, fractions = rrtmg_optical_depth(
        rad.rrtmg,
        vapour=value("atmosphere_water_vapour", "kg m^-2"),
        layer_mass=mass,
        pressure=pressure,
        temperature=temperature,
        ozone_path=rad.ozone_path,
        with_fractions=True,
    )
    common = dict(
        grid=rad.grid,
        incident_spectral_irradiance=value("toa_spectral_irradiance", "W m^-2 m^-1"),
        layer_temperature=temperature,
        layer_mass=mass,
        absorption_coefficient=rad.absorption_coefficient,
        surface_temperature=surface,
        surface_albedo=value("surface_spectral_albedo", "1"),
        diffusivity_factor=rad.diffusivity_factor,
        absorption_optical_depth=tau,
    )
    direct = solve_column_radiation(
        **common,
        planck_fractions=rrtmg_planck_fractions(rad.rrtmg, fractions),
        interface_temperature=interface_temperatures(
            temperature, pressure, value("atmosphere_interface_pressure", "Pa"), surface
        ),
    )
    thermal = rad.rrtmg.gpoint >= 0
    np.testing.assert_allclose(
        np.asarray(rad.clear_sky.upward[-1])[thermal],
        np.asarray(direct.upward[-1])[thermal],
        rtol=1e-9,
    )
    plain = solve_column_radiation(**common)
    assert direct.surface_emitted == pytest.approx(plain.surface_emitted, rel=1e-12)
    assert abs(rad.clear_sky.energy_residual()) < 1e-9
    assert rad.parameters()["k_distribution_source"] is True
    assert rad.parameters()["linear_in_tau_source"] is True


@pytest.mark.req("RAD-002", "EQR-002")
def test_interface_temperatures_follow_cam_s_convention():
    """Interior interfaces lie linearly in ln p between layers; the top takes the top
    layer's temperature and the bottom the surface's."""
    t = np.array([280.0, 260.0, 240.0])
    p = np.array([9.0e4, 7.0e4, 5.0e4])
    pi = np.array([1.0e5, 8.0e4, 6.0e4, 4.0e4])
    out = interface_temperatures(t, p, pi, 290.0)
    assert out[0] == 290.0
    assert out[-1] == 240.0
    expected = 280.0 + (260.0 - 280.0) * np.log(8.0e4 / 9.0e4) / np.log(7.0e4 / 9.0e4)
    assert out[1] == pytest.approx(expected, rel=1e-14)


@pytest.mark.req("EQR-002", "VAL-001")
@pytest.mark.parametrize(
    ("atmosphere", "predicted"),
    [
        ("MLS", (-0.60, 0.24)),
        ("SAW", (-0.16, 0.49)),
        ("TROP", (-0.97, 1.12)),
        ("MLW", (-0.23, 0.07)),
    ],
)
def test_the_solver_reproduces_rrtmg_s_transfer(atmosphere, predicted):
    """P2: on AER's atmospheres, the grid to 2600 cm^-1, against AER's total less
    band 16, OLR within 1.5 and surface downward within 2.0 W m^-2; and the residual
    is P1's (d), the diffusivity alone, within 0.15, as predicted before the code."""
    from test_rrtmg_longwave import ROOT, read_input, read_output, rtrnmr

    from mfespf.physics.rrtmg_longwave import (
        RRTMGLongwave,
        rrtmg_gas_optics,
        rrtmg_longwave_grid,
        rrtmg_planck_fractions,
    )

    col, pz, tz, tbound, pwvcm = read_input(ROOT / f"input_rrtm_{atmosphere}-clr")
    total = read_output(ROOT / f"output_rrtm_{atmosphere}-clr")[(10.0, 3250.0)]
    taug, fracs, _ = rrtmg_gas_optics(col)
    up16, dn16 = rtrnmr(col, pz, tz, tbound, pwvcm, taug, fracs, [15], band16_total=True)
    olr_ref, down_ref = total[-1, 2] - up16[-1], total[0, 3] - dn16[0]

    grid, gp = rrtmg_longwave_grid(column_grid(16), upper_cm=2600.0)
    lw = gp >= 0
    spec = RRTMGLongwave(gpoint=gp, co2_ppmv=0.0)
    tau = np.zeros((col.pressure_hpa.size, grid.n_bins))
    tau[:, lw] = taug[:, gp[lw]]
    layer = np.zeros_like(tau)
    layer[:, lw] = fracs[:, gp[lw]]
    sol = solve_column_radiation(
        grid=grid,
        incident_spectral_irradiance=np.zeros(grid.n_bins),
        layer_temperature=col.temperature,
        layer_mass=(pz[:-1] - pz[1:]) * 100.0 / 9.80665,
        absorption_coefficient=np.zeros(grid.n_bins),
        surface_temperature=tbound,
        surface_albedo=np.zeros(grid.n_bins),
        absorption_optical_depth=tau,
        planck_fractions=rrtmg_planck_fractions(spec, layer),
        interface_temperature=tz,
    )
    d_olr = float(np.sum(sol.upward[-1][lw])) - olr_ref
    d_down = float(np.sum(sol.downward[0][lw])) - down_ref
    assert abs(d_olr) <= 1.5
    assert abs(d_down) <= 2.0
    assert d_olr == pytest.approx(predicted[0], abs=0.15)
    assert d_down == pytest.approx(predicted[1], abs=0.15)
