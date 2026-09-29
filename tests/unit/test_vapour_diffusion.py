"""Vertical eddy diffusion of water vapour.

Checked against closed forms rather than against the model: the column total
is unchanged, a uniform humidity does not move, and in an isothermal column
with humidity linear in ``ln p`` -- linear in height -- the flux at every
interface is exactly ``-rho K dq/dz``.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError
from mfespf.physics.water import eddy_diffusion_rate

pytestmark = pytest.mark.stage(3)

G = 9.80665
R = value_of("gas_constant_dry_air")
N = 20
INTERFACES = np.linspace(101325.0, 1000.0, N + 1)
PRESSURE = 0.5 * (INTERFACES[:-1] + INTERFACES[1:])
MASS = (INTERFACES[:-1] - INTERFACES[1:]) / G
T0 = 250.0
TEMPERATURE = np.full(N, T0)
K = 0.39
DT = 3600.0


def rate_for(vapour: np.ndarray, diffusivity: float = K) -> np.ndarray:
    return eddy_diffusion_rate(
        vapour=vapour,
        layer_mass=MASS,
        temperature=TEMPERATURE,
        pressure=PRESSURE,
        interface_pressure=INTERFACES,
        diffusivity=diffusivity,
        dt=DT,
    )


@pytest.mark.req("WAT-001", "CON-001")
def test_eddy_diffusion_conserves_the_column():
    rng = np.random.default_rng(3)
    for _ in range(20):
        vapour = rng.uniform(0.0, 0.02, N) * MASS
        rate = rate_for(vapour)
        assert abs(rate.sum()) <= 1e-15 * np.abs(rate).sum()


@pytest.mark.req("WAT-001", "VAL-001")
def test_a_uniform_humidity_does_not_diffuse():
    rate = rate_for(4.0e-6 * MASS)
    assert np.all(np.abs(rate) < 1e-20)
    assert np.all(rate_for(np.zeros(N)) == 0.0)


@pytest.mark.req("WAT-001", "VAL-001")
def test_the_flux_is_minus_rho_k_dq_dz():
    """Isothermal, so z = (R T / g) ln(p_s / p): q linear in ln p is linear in z.

    ``dq/dz = -b g / (R T)`` for ``q = a + b ln p``, so the upward flux at an
    interface of pressure ``p_i`` is ``rho_i K b g / (R T)`` with
    ``rho_i = p_i / (R T)``, exactly.
    """
    a, b = 1.0e-3, 2.0e-4
    vapour = (a + b * np.log(PRESSURE)) * MASS
    rate = rate_for(vapour)
    # Rebuild the interface fluxes from the divergences, surface upward.
    flux = np.concatenate(([0.0], -np.cumsum(rate)))
    expected = INTERFACES / (R * T0) * K * b * G / (R * T0)
    assert flux[1:N] == pytest.approx(expected[1:N], rel=1e-10)
    assert flux[N] == pytest.approx(0.0, abs=1e-12 * np.abs(flux).max())
    # Moister below, so the flux is upward everywhere inside.
    assert np.all(flux[1:N] > 0.0)


@pytest.mark.req("WAT-001", "VAL-001")
def test_an_unstable_step_is_refused():
    """Explicit diffusion keeps layers non-negative only below one mass per step."""
    with pytest.raises(ConfigurationError, match="unstable"):
        rate_for(0.01 * MASS, diffusivity=1.0e4)
    with pytest.raises(ConfigurationError, match="negative"):
        rate_for(0.01 * MASS, diffusivity=-1.0)
