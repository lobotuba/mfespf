"""RRTMG's zero-floor, its tabulated domain, and how a refusal is reported.

rrtmg-domain-and-refusals.md. RRTMG's binary-species interpolation can give a
g-point a small negative optical depth, and its transfer floors it at zero; the
k-tables cover T_ref(p) +/- 30 K, which the validity gate now judges; and a
module refusing mid-run is named as a refusal, not as a conservation failure.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError, ConservationError
from mfespf.physics.rrtmg_longwave import (
    _AMD,
    _AMW,
    _AVOGADRO,
    _LOSCHMIDT_PER_CM3,
    RRTMGColumn,
    RRTMGLongwave,
    rrtmg_gas_optics,
    rrtmg_longwave_grid,
    rrtmg_optics,
    rrtmg_reference_departure,
    rrtmg_tables,
)
from mfespf.physics.thermodynamics import saturation_vapour_pressure
from mfespf.run.column_configurations import column_grid, reference_column

pytestmark = pytest.mark.stage(6)


def _reference_state(relative_humidity: float):
    """RRTMG's own reference profile, 1050 to 1 hPa, at a uniform relative humidity."""
    t = rrtmg_tables()
    edges = np.geomspace(105000.0, 100.0, 81)
    p = np.sqrt(edges[:-1] * edges[1:])
    T = np.interp(-np.log(p / 100.0), -t["preflog"], t["tref"])
    mass = (edges[:-1] - edges[1:]) / 9.80665
    e = relative_humidity * saturation_vapour_pressure(T)
    # Above 100 hPa saturation is meaningless at these temperatures; 4 ppm by mass.
    q = np.where(p > 1.0e4, 0.622 * e / (p - 0.378 * e), 4.0e-6)
    return p, T, mass, q * mass


@pytest.mark.req("EQR-002", "VAL-001")
def test_rrtmg_s_zero_floor_is_applied_and_counted():
    """T1: on the reference profile at RH 0.5 no optical depth is negative after the
    floor, and the count is exactly the negatives the raw optics give there."""
    p, T, mass, vapour = _reference_state(0.5)
    _, gp = rrtmg_longwave_grid(column_grid(16))
    spec = RRTMGLongwave(gpoint=gp, co2_ppmv=284.7)
    dry = mass - vapour
    coldry = dry * 0.1 / _AMD * _AVOGADRO
    ozone_path = 1.0e-6 * coldry / _LOSCHMIDT_PER_CM3
    optics = rrtmg_optics(
        spec, vapour=vapour, layer_mass=mass, pressure=p, temperature=T, ozone_path=ozone_path
    )
    assert np.all(optics.tau >= 0.0)
    raw, _, _ = rrtmg_gas_optics(
        RRTMGColumn.from_mixing_ratios(
            pressure_hpa=p / 100.0,
            temperature=T,
            coldry=coldry,
            vmr={
                "h2o": vapour / dry * _AMD / _AMW,
                "co2": 284.7e-6,
                "o3": np.full(p.size, 1.0e-6),
                "ch4": 808.2e-9,
                "n2o": 273.0e-9,
                "o2": 0.209,
            },
        )
    )
    on_grid = raw[:, gp[gp >= 0]]
    assert optics.floored == int(np.count_nonzero(on_grid < 0.0)) > 0
    assert optics.most_negative == pytest.approx(float(on_grid.min()), rel=1e-12)


def _column(**overrides):
    options = dict(
        grid=column_grid(16),
        n_layers=10,
        master_dt=3600.0,
        measured_solar_spectrum=True,
        ozone=True,
        solar_vapour_absorption=True,
        rrtmg_longwave=True,
    )
    options.update(overrides)
    return reference_column(**options)


def _radiation(kernel):
    return next(s.module for s in kernel.schedule.ordered() if s.module_id == "radiation")


@pytest.mark.req("EQR-002", "VAL-001")
def test_the_gate_judges_rrtmg_s_departure_from_its_tables():
    """T2: the published departure is T - T_ref(p); one layer 35 K above the
    reference fails the validity gate by name, and 25 K above passes it."""
    probe = _column()
    probe.run(1)
    rad = _radiation(probe)
    p = np.asarray(rad.input_field("atmosphere_layer_pressure").values, dtype=float)
    T = np.asarray(rad.input_field("atmosphere_layer_temperature").values, dtype=float)
    published = np.asarray(rad.emit_outputs()["rrtmg_reference_temperature_departure"].values)
    t = rrtmg_tables()
    by_hand = T - np.interp(-np.log(p / 100.0), -t["preflog"], t["tref"])
    np.testing.assert_allclose(published, by_hand, rtol=0, atol=1e-12)
    np.testing.assert_allclose(rrtmg_reference_departure(p, T), by_hand, rtol=0, atol=1e-12)

    reference = T - by_hand
    for offset, passes in ((25.0, True), (35.0, False)):
        start = reference.copy()
        start[3] += offset
        run = _column(initial_air_temperature=start).run(1)
        assert run.manifest.gate_results["validity_domains"] is passes, offset
        if not passes:
            assert any("dT_ref" in w and "RRTMG" in w for w in run.manifest.warnings)


def _raise_on_third_call(module, error):
    original = module.compute_tendencies
    calls = {"n": 0}

    def wrapped():
        calls["n"] += 1
        if calls["n"] == 3:
            raise error
        return original()

    module.compute_tendencies = wrapped


@pytest.mark.req("KER-005", "CON-004")
def test_a_refusal_is_named_and_not_filed_as_conservation():
    """T3: a module raising ConfigurationError on its third step fails run_completed,
    naming the module and the step; a ConservationError still fails conservation."""
    kernel = _column(rrtmg_longwave=False, ozone=False)
    _raise_on_third_call(_radiation(kernel), ConfigurationError("refused on purpose"))
    run = kernel.run(5)
    gates = run.manifest.gate_results
    assert gates["run_completed"] is False
    assert gates.get("conservation_reconciliation", True) is True
    assert any(
        "ConfigurationError raised by radiation at step 2: refused on purpose" in w
        for w in run.manifest.warnings
    )

    kernel = _column(rrtmg_longwave=False, ozone=False)
    _raise_on_third_call(_radiation(kernel), ConservationError("leak on purpose"))
    run = kernel.run(5)
    assert run.manifest.gate_results["conservation_reconciliation"] is False
    assert "run_completed" not in run.manifest.gate_results
