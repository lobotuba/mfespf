"""The overturning circulation between a rising and a sinking column.

Checked against its own defining statements rather than against the
observations it is later judged by: whatever leaves one column arrives in the
other, exactly; a temperature difference in the gravest mode relaxes on
``tau_1`` and one in mode ``n`` on ``n tau_1``; no difference, no circulation;
the columns carry equal and opposite mass; and the global mean of the regions
preserves the integral.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError
from mfespf.core.quantities import Field, Provenance, QuantityKind, SpatialSupport, TimeSemantics
from mfespf.physics.circulation import (
    AreaWeightedMeansModule,
    spectral_wtg_subsidence,
    upwind_tendency,
)
from mfespf.run.column_configurations import column_grid
from mfespf.run.two_column import two_column

pytestmark = pytest.mark.stage(3)

G = 9.80665
C_P = 1004.0
TROPOPAUSE = 10000.0
#: Bulk dS/dp, J kg^-1 Pa^-1: about 60 kJ/kg across the troposphere.
STRATIFICATION = -0.65
TAU = 1.3e4


def grid(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Interfaces from the surface to 50 hPa, and the layers' zeta."""
    interfaces = np.linspace(101325.0, 5000.0, n + 1)
    centres = 0.5 * (interfaces[:-1] + interfaces[1:])
    zeta = (interfaces[0] - centres) / (interfaces[0] - TROPOPAUSE)
    return interfaces, zeta


def flux_for(anomaly: np.ndarray, interfaces: np.ndarray) -> np.ndarray:
    return spectral_wtg_subsidence(
        temperature_anomaly=anomaly,
        specific_heat=np.full(anomaly.size, C_P),
        bulk_stratification=STRATIFICATION,
        relaxation_time=TAU,
        interface_pressure=interfaces,
        tropopause=TROPOPAUSE,
        gravity=G,
    )


def warming_rate(flux: np.ndarray) -> np.ndarray:
    """dT/dt from subsidence across the uniform stratification, per layer."""
    omega = G * 0.5 * (flux[:-1] + flux[1:])
    return -omega * STRATIFICATION / C_P


def mode_shape(n: int, zeta: np.ndarray) -> np.ndarray:
    return np.where(zeta < 1.0, np.sin(n * np.pi * np.clip(zeta, 0.0, 1.0)), 0.0)


@pytest.mark.req("EQR-002", "VAL-001")
def test_the_gravest_mode_relaxes_on_tau_1():
    """A sinking column 1 K colder in the gravest mode warms at 1 K per tau_1."""
    interfaces, zeta = grid(400)
    anomaly = -1.0 * mode_shape(1, zeta)
    rate = warming_rate(flux_for(anomaly, interfaces))
    inside = zeta < 1.0
    assert rate[inside] == pytest.approx(-anomaly[inside] / TAU, abs=2e-3 / TAU)
    # Layers wholly above the tropopause are untouched.
    assert np.all(rate[interfaces[:-1] <= TROPOPAUSE] == 0.0)


@pytest.mark.req("EQR-002", "VAL-001")
def test_each_mode_relaxes_n_times_more_slowly():
    """Mode n travels at c_1 / n, so its difference relaxes on n tau_1."""
    interfaces, zeta = grid(400)
    inside = zeta < 1.0
    for n in (2, 3, 5):
        anomaly = -1.0 * mode_shape(n, zeta)
        rate = warming_rate(flux_for(anomaly, interfaces))
        assert rate[inside] == pytest.approx(-anomaly[inside] / (n * TAU), abs=5e-3 / TAU)


@pytest.mark.req("EQR-002", "VAL-001")
def test_no_temperature_difference_means_no_circulation():
    interfaces, _ = grid(30)
    assert np.all(flux_for(np.zeros(30), interfaces) == 0.0)
    # Colder than the mean sinks, warmer rises, and nothing crosses the
    # surface or the tropopause.
    colder = flux_for(np.full(30, -0.5), interfaces)
    warmer = flux_for(np.full(30, 0.5), interfaces)
    free = (interfaces > TROPOPAUSE) & (interfaces < interfaces[0])
    assert np.all(colder[free] > 0.0)
    assert warmer == pytest.approx(-colder, rel=1e-12)
    assert colder[0] == 0.0 and np.all(colder[interfaces <= TROPOPAUSE] == 0.0)
    with pytest.raises(ConfigurationError, match="not stable"):
        spectral_wtg_subsidence(
            temperature_anomaly=np.full(30, -0.5),
            specific_heat=np.full(30, C_P),
            bulk_stratification=0.1,
            relaxation_time=TAU,
            interface_pressure=interfaces,
        )


@pytest.mark.req("CON-001", "CON-004", "EQR-002")
def test_the_circulation_conserves_energy_and_water_exactly():
    """Whatever leaves one column arrives in the other, and a run closes."""
    rng = np.random.default_rng(7)
    sigma = 0.476
    interfaces, _ = grid(12)
    flux_down = flux_for(rng.uniform(-1.0, 0.2, 12), interfaces)
    flux_up = -(1.0 - sigma) / sigma * flux_down
    for _ in range(50):
        phi_up, phi_down = rng.uniform(0.0, 0.02, 12), rng.uniform(0.0, 0.02, 12)
        rate_down, _ = upwind_tendency(flux_down, phi_down, phi_up)
        rate_up, _ = upwind_tendency(flux_up, phi_up, phi_down)
        total = (1.0 - sigma) * rate_down.sum() + sigma * rate_up.sum()
        assert total == pytest.approx(0.0, abs=1e-15 * np.abs(rate_down).sum())

    run = two_column(grid=column_grid(16), areal_heat_capacity=1.0e6).run(200)
    assert run.may_support_a_claim, run.manifest.invalidation_reasons
    assert run.worst_relative_residual() < 1e-10


@pytest.mark.req("CON-001", "EQR-002")
def test_the_two_columns_carry_equal_and_opposite_mass():
    """Layer by layer, the lateral mass one column loses the other gains."""
    sigma = 0.476
    interfaces, zeta = grid(12)
    flux_down = flux_for(-1.0 * mode_shape(1, zeta), interfaces)
    flux_up = -(1.0 - sigma) / sigma * flux_down
    lateral_down = flux_down[:-1] - flux_down[1:]
    lateral_up = flux_up[:-1] - flux_up[1:]
    assert (1.0 - sigma) * lateral_down + sigma * lateral_up == pytest.approx(
        np.zeros(12), abs=1e-18
    )
    # The sinking column takes air in aloft and sends it out low down.
    assert lateral_down[-4:].sum() > 0.0 and lateral_down[:3].sum() < 0.0


@pytest.mark.req("SPA-003", "CON-001")
def test_the_global_mean_preserves_the_integral():
    """Regions to a global scalar: the area-weighted sum, exactly (SPA-003)."""
    module = AreaWeightedMeansModule(
        regions=(("a", 0.3), ("b", 0.7)), fields=(("outgoing_longwave_radiation", "W m^-2"),)
    )

    def field(value: float) -> Field:
        return Field(
            field_id="outgoing_longwave_radiation",
            physical_name="olr",
            values=np.array(value),
            units="W m^-2",
            quantity_kind=QuantityKind.FLUX,
            time_semantics=TimeSemantics.INTERVAL_AVERAGE,
            time_interval=3600.0,
            timestamp=0.0,
            spatial_support=SpatialSupport.GLOBAL_SCALAR,
            provenance=Provenance(producer="test"),
        )

    class _Context:
        dt = 3600.0
        time = 0.0

    module._inputs = {
        "a:outgoing_longwave_radiation": field(200.0),
        "b:outgoing_longwave_radiation": field(280.0),
    }
    module._context = _Context()
    out = module.emit_outputs()["outgoing_longwave_radiation"]
    assert float(out.values) == pytest.approx(0.3 * 200.0 + 0.7 * 280.0, rel=1e-15)
