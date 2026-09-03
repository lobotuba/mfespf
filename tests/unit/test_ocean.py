"""Stage 5: the reduced two-reservoir ocean.

Charter 9.10 and requirements OCN-001, OCN-002, EQR-002.

The exchange law is one line of algebra, so the tests that matter are about its
limits and about the linear system it produces with the mixed layer above it --
whose eigenvalues and frequency response are known in closed form and are
therefore something to check against rather than something to report.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from mfespf.core.constants import ParameterClass, get
from mfespf.core.errors import ConfigurationError
from mfespf.physics.ocean import (
    OCEAN_DEEP_STATE,
    DeepOceanModule,
    areal_heat_capacity_of_water,
    deep_ocean_heat_flux,
    two_reservoir_modes,
    two_reservoir_transfer_function,
)

pytestmark = pytest.mark.stage(5)

C_MIX = areal_heat_capacity_of_water(50.0)
C_DEEP = areal_heat_capacity_of_water(1200.0)
LAMBDA = 3.742
GAMMA = 0.7


# --- the exchange law ------------------------------------------------------


@pytest.mark.req("OCN-001", "VAL-001")
def test_equal_temperatures_exchange_no_heat():
    """The limiting case a spun-up configuration depends on.

    If this were only approximately zero, every configuration with an ocean
    under it would drift slowly for no physical reason, and the drift would be
    indistinguishable from a result.
    """
    assert deep_ocean_heat_flux(288.0, 288.0, GAMMA) == 0.0
    assert deep_ocean_heat_flux(288.0, 288.0, 0.0) == 0.0


@pytest.mark.req("OCN-001", "EQR-002")
def test_the_exchange_is_linear_and_signed():
    """Heat flows down a temperature gradient and up when the gradient reverses."""
    assert deep_ocean_heat_flux(290.0, 288.0, GAMMA) == pytest.approx(2.0 * GAMMA)
    assert deep_ocean_heat_flux(286.0, 288.0, GAMMA) == pytest.approx(-2.0 * GAMMA)
    doubled = deep_ocean_heat_flux(292.0, 288.0, GAMMA)
    assert doubled == pytest.approx(2.0 * deep_ocean_heat_flux(290.0, 288.0, GAMMA))


@pytest.mark.req("OCN-002")
def test_heat_capacity_follows_from_a_declared_depth():
    """Depth is checkable against a measurement; the product is not."""
    assert areal_heat_capacity_of_water(100.0) == pytest.approx(2.0 * C_MIX)
    with pytest.raises(ConfigurationError, match="must be positive"):
        areal_heat_capacity_of_water(0.0)


@pytest.mark.req("EQR-004", "SYS-003")
def test_the_exchange_coefficient_is_declared_calibrated():
    """It was fitted to model behaviour, and EQR-004 forbids dressing that up.

    The classification is not decoration: a claim about heat-uptake magnitude
    rests on this number, and a reader has to be able to see that the number
    came from a fit rather than from a measurement on seawater.
    """
    constant = get("ocean_deep_exchange_coefficient")
    assert constant.parameter_class is ParameterClass.CALIBRATED
    assert constant.standard_uncertainty is not None
    assert constant.relative_uncertainty > 0.3  # poorly constrained, and it says so


# --- the two-reservoir system ----------------------------------------------


@pytest.mark.req("OCN-002", "VAL-002", "TIM-008")
def test_the_analytic_modes_are_the_eigenvalues_of_the_linear_system():
    """Checked against a numerically computed eigendecomposition, not itself."""
    fast, slow = two_reservoir_modes(
        mixed_layer_capacity=C_MIX,
        deep_capacity=C_DEEP,
        feedback_parameter=LAMBDA,
        exchange_coefficient=GAMMA,
    )
    matrix = np.array(
        [
            [-(LAMBDA + GAMMA) / C_MIX, GAMMA / C_MIX],
            [GAMMA / C_DEEP, -GAMMA / C_DEEP],
        ]
    )
    eigenvalues = np.sort(np.linalg.eigvals(matrix).real)  # most negative first
    assert fast == pytest.approx(-1.0 / eigenvalues[0], rel=1e-12)
    assert slow == pytest.approx(-1.0 / eigenvalues[1], rel=1e-12)
    assert fast < slow


@pytest.mark.req("OCN-002", "TIM-008")
def test_a_deep_reservoir_speeds_the_surface_up_and_slows_the_system_down():
    """Both halves of the two-timescale claim, stated as inequalities.

    The fast mode is *shorter* than the slab's own relaxation time, because the
    deep store gives heat somewhere else to go; the slow mode is far longer,
    because that heat has to come back out before the system finishes adjusting.
    A model with only the second effect would be missing the first.
    """
    slab = C_MIX / LAMBDA
    fast, slow = two_reservoir_modes(
        mixed_layer_capacity=C_MIX,
        deep_capacity=C_DEEP,
        feedback_parameter=LAMBDA,
        exchange_coefficient=GAMMA,
    )
    assert fast < slab
    assert slow > 20.0 * slab


@pytest.mark.req("OCN-002", "VAL-001")
def test_a_vanishing_exchange_recovers_the_single_reservoir_response():
    """The disabled-coupling control for the whole module (MOD-014)."""
    omega = 2.0 * math.pi / (4000.0 * 86400.0)
    slab = 1.0 / (1j * omega * C_MIX + LAMBDA)
    with_ocean = two_reservoir_transfer_function(
        omega,
        mixed_layer_capacity=C_MIX,
        deep_capacity=C_DEEP,
        feedback_parameter=LAMBDA,
        exchange_coefficient=0.0,
    )
    assert with_ocean == pytest.approx(slab, rel=1e-14)


@pytest.mark.req("OCN-002", "TIM-008")
def test_the_deep_reservoir_lowers_the_gain_and_the_phase_lag_together():
    """A prediction worth stating because it looks like a bug.

    The deep store is another place for heat to go, so the gain falls. But the
    same branch contributes a zero to the transfer function, and the phase lag
    falls too -- a two-reservoir system lags *less* than a slab at the same
    frequency. A measured lag smaller than the slab's is therefore the expected
    result, not a sign that the coupling ran backwards.
    """
    omega = 2.0 * math.pi / (4000.0 * 86400.0)
    common = dict(mixed_layer_capacity=C_MIX, deep_capacity=C_DEEP, feedback_parameter=LAMBDA)
    slab = two_reservoir_transfer_function(omega, exchange_coefficient=0.0, **common)
    ocean = two_reservoir_transfer_function(omega, exchange_coefficient=GAMMA, **common)
    assert abs(ocean) < abs(slab)
    assert -np.angle(ocean) < -np.angle(slab)


@pytest.mark.req("OCN-002")
def test_the_gain_approaches_the_equilibrium_sensitivity_at_low_frequency():
    """As the forcing slows, both reservoirs equilibrate and the gain is 1/lambda."""
    common = dict(
        mixed_layer_capacity=C_MIX,
        deep_capacity=C_DEEP,
        feedback_parameter=LAMBDA,
        exchange_coefficient=GAMMA,
    )
    slow = two_reservoir_transfer_function(2.0 * math.pi / 1.0e14, **common)
    assert abs(slow) == pytest.approx(1.0 / LAMBDA, rel=1e-3)
    assert -np.angle(slow) == pytest.approx(0.0, abs=1e-3)


@pytest.mark.req("SYS-006")
def test_unphysical_parameters_are_refused():
    common = dict(mixed_layer_capacity=C_MIX, deep_capacity=C_DEEP, feedback_parameter=LAMBDA)
    with pytest.raises(ConfigurationError, match="exchange_coefficient must be positive"):
        two_reservoir_modes(exchange_coefficient=0.0, **common)
    with pytest.raises(ConfigurationError, match="angular frequency"):
        two_reservoir_transfer_function(0.0, exchange_coefficient=GAMMA, **common)
    with pytest.raises(ConfigurationError, match="non-negative"):
        DeepOceanModule(exchange_coefficient=-1.0).validate_configuration()
    with pytest.raises(ConfigurationError, match="initial deep temperature"):
        DeepOceanModule(initial_temperature=0.0).validate_configuration()


# --- the module ------------------------------------------------------------


@pytest.mark.req("MOD-001", "MOD-002", "OCN-002")
def test_the_module_contract_declares_what_it_does_not_do():
    """Charter 9.10 lists salinity, circulation and ice; this fidelity has none.

    An omission that is declared is a scope statement. The same omission left
    out of the contract would let a reader infer that the module represents
    something it does not.
    """
    module = DeepOceanModule()
    module.validate_configuration()
    limitations = " ".join(module.contract.known_limitations)
    for absent in ("salinity", "circulation", "sea ice", "calibrated"):
        assert absent in limitations
    assert module.contract.owned_state_ids == frozenset({OCEAN_DEEP_STATE})
    assert "surface_temperature" in module.contract.required_inputs
    assert "ocean_heat_uptake" in module.contract.produced_outputs
    # It owns a reservoir, not a surface: the sea-surface temperature stays with
    # whichever module owns the mixed layer (STA-001).
    assert "sea_surface_temperature" not in module.contract.produced_outputs


@pytest.mark.req("AUD-002")
def test_the_reservoir_reaches_the_configuration_hash():
    """Two oceans that differ physically must not record identical parameters."""
    shallow = DeepOceanModule(depth=400.0).parameters()
    deep = DeepOceanModule(depth=1200.0).parameters()
    assert shallow["areal_heat_capacity_J_m2_K"] != deep["areal_heat_capacity_J_m2_K"]
    assert deep["exchange_coefficient_parameter_class"] == "calibrated"


@pytest.mark.req("MOD-009")
def test_the_module_checkpoints_and_restores_its_only_state():
    module = DeepOceanModule(initial_temperature=277.0)
    module.initialize({}, {})
    checkpoint = module.checkpoint()
    module.temperature = 300.0
    module.restore(checkpoint)
    assert module.temperature == pytest.approx(277.0)
