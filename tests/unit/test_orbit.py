"""Stage 1: orbital geometry correctness gates (ORB-001, ORB-002)."""

from __future__ import annotations

import math

import pytest

from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError
from mfespf.physics.orbit import orbit_radius, solve_kepler, spherical_mean_factor

pytestmark = pytest.mark.stage(1)

AU = value_of("astronomical_unit")
YEAR = 365.256363004 * 86400.0


@pytest.mark.req("ORB-002", "VAL-001")
def test_intercepted_equals_distributed_power():
    """The one identity ORB-002 asks for, checked for arbitrary radius and flux.

    A sphere intercepts S over its projected disc and re-emits over its whole
    surface. If the geometry factor were wrong, the coupled energy budget would
    still close -- it would simply close around the wrong number, which is why
    this is checked against geometry rather than against the model.
    """
    for radius in (1.0, 6.371e6, 7.0e7):
        for irradiance in (1.0, 1361.0):
            intercepted = irradiance * math.pi * radius**2
            distributed = irradiance * spherical_mean_factor() * 4.0 * math.pi * radius**2
            assert distributed == pytest.approx(intercepted, rel=1e-15)


@pytest.mark.req("ORB-001")
def test_global_mean_factor_is_one_quarter():
    assert spherical_mean_factor() == 0.25


@pytest.mark.req("ORB-001", "VAL-001")
def test_kepler_circular_limit():
    """At zero eccentricity the distance is the semi-major axis at every phase."""
    for fraction in (0.0, 0.1, 0.25, 0.5, 0.75, 0.99):
        assert orbit_radius(fraction * YEAR, eccentricity=0.0) == pytest.approx(AU)


@pytest.mark.req("ORB-001")
def test_kepler_perihelion_and_aphelion():
    e = 0.0167086
    perihelion = orbit_radius(0.0, eccentricity=e)
    aphelion = orbit_radius(YEAR / 2.0, eccentricity=e)
    assert perihelion == pytest.approx(AU * (1.0 - e), rel=1e-12)
    assert aphelion == pytest.approx(AU * (1.0 + e), rel=1e-12)
    # Irradiance varies by roughly 6.8% between the two.
    ratio = (aphelion / perihelion) ** 2
    assert ratio == pytest.approx(1.068, abs=0.002)


@pytest.mark.req("ORB-001")
def test_kepler_radius_symmetry_about_perihelion():
    e = 0.2
    for offset in (0.1, 0.3, 0.45):
        forward = orbit_radius(offset * YEAR, eccentricity=e)
        backward = orbit_radius(-offset * YEAR, eccentricity=e)
        assert forward == pytest.approx(backward, rel=1e-12)


@pytest.mark.req("VAL-001")
def test_kepler_solver_residual_below_tolerance():
    for e in (0.0, 0.1, 0.5, 0.8):
        for M in (0.0, 0.7, math.pi, 5.0):
            E = solve_kepler(M, e)
            assert E - e * math.sin(E) - M == pytest.approx(0.0, abs=1e-12)


@pytest.mark.req("MOD-002", "SYS-006")
def test_eccentricity_outside_the_validity_domain_is_refused():
    with pytest.raises(ConfigurationError, match="validity domain"):
        solve_kepler(1.0, 0.95)


@pytest.mark.req("ORB-001", "VAL-001")
def test_inverse_square_scaling():
    """Doubling the distance quarters the irradiance."""
    near = orbit_radius(0.0, eccentricity=0.0, semi_major_axis_m=AU)
    far = orbit_radius(0.0, eccentricity=0.0, semi_major_axis_m=2.0 * AU)
    assert (AU / far) ** 2 == pytest.approx(0.25 * (AU / near) ** 2)


@pytest.mark.req("ORB-001")
def test_inverse_square_identity_at_the_reference_distance():
    assert (AU / orbit_radius(0.0, eccentricity=0.0)) ** 2 == pytest.approx(1.0)
