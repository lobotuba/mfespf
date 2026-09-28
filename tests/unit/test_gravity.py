"""Gravity through a replaceable interface, and what the constant one costs.

Charter 9.15 and requirement GRV-001. Until this existed gravity was
``value_of("earth_standard_gravity")`` read straight from the constants table at
four call sites, and the constant's own note said the gravity module supersedes
it. The requirement is not really about gravity: it is that an approximation
nobody has to declare is one nobody has to measure, and the hydrostatic
equation's record proved the point by asserting a bound that was wrong for the
column it was attached to.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError
from mfespf.physics.atmosphere import AtmosphericColumnModule, hydrostatic_layer_masses
from mfespf.physics.gravity import (
    EARTH_MEAN_RADIUS_M,
    ConstantGravity,
    NewtonianAltitudeGravity,
    as_gravity_field,
)
from mfespf.physics.thermodynamics import hydrostatic_heights

pytestmark = pytest.mark.stage(2)

G0 = value_of("earth_standard_gravity")


@pytest.mark.req("GRV-001")
def test_the_constant_field_is_the_declared_first_rung():
    """Charter 9.15's first rung, now saying what it is.

    The conventional value is exact by definition, so nothing here is measuring
    the number. What is being established is that it arrives as a field that
    declares an approximation, rather than as a float nobody has to account for.
    """
    field = ConstantGravity()
    assert field.surface_value() == pytest.approx(G0)
    assert np.all(field.at_height(np.array([0.0, 1.0e3, 1.0e5])) == G0)
    assert not field.varies_with_height
    assert field.approximations, "a field that declares no approximation declares nothing"
    assert any("z / R" in a or "z_/_R" in a for a in field.approximations)
    assert field.validity_domain.strip()

    with pytest.raises(ConfigurationError, match="must be positive"):
        ConstantGravity(-1.0)


@pytest.mark.req("GRV-001", "EQR-002")
def test_the_newtonian_field_is_the_inverse_square_law():
    """Checked against the closed form and against its own series expansion.

    The second check is the one worth having. ``g(z)/g_0 = (1 + z/R)^-2`` expands
    to ``1 - 2x + 3x^2 - 4x^3 + ...`` in ``x = z/R``. Agreeing with the closed
    form would only establish that the code computes the formula it was handed;
    agreeing with an independently derived expansion *term by term*, down to the
    size of the first term dropped, establishes that the formula is the right
    one.
    """
    field = NewtonianAltitudeGravity()
    assert field.surface_value() == pytest.approx(G0)

    z = np.array([0.0, 1.0e3, 1.0e4, 2.0e4, 1.0e5])
    ratio = np.asarray(field.at_height(z)) / G0
    closed_form = (EARTH_MEAN_RADIUS_M / (EARTH_MEAN_RADIUS_M + z)) ** 2
    assert ratio == pytest.approx(closed_form, rel=1e-12)

    x = z / EARTH_MEAN_RADIUS_M
    first_order = 1.0 - 2.0 * x
    second_order = first_order + 3.0 * x**2
    # Each truncation is wrong by the first term it dropped, which is the whole
    # content of a convergent expansion and is not something a wrong formula
    # would reproduce.
    assert ratio == pytest.approx(first_order, abs=4.0 * x.max() ** 2)
    assert (second_order - ratio)[1:] == pytest.approx(4.0 * x[1:] ** 3, rel=0.05)

    assert field.at_height(0.0) == pytest.approx(G0)  # the limiting case
    with pytest.raises(ConfigurationError, match="below the centre"):
        field.at_height(-2.0 * EARTH_MEAN_RADIUS_M)


@pytest.mark.req("GRV-001", "MOD-001")
def test_the_interface_is_replaceable_and_the_replacement_is_visible():
    """Replaceable is the requirement, and hashable is what makes it honest.

    Two runs that differ in gravity must not hash alike, so each field publishes
    its own parameters (AUD-002) rather than a single number that would make a
    constant field and an altitude-dependent one with the same surface value
    indistinguishable in the audit trail.
    """
    constant = ConstantGravity()
    newtonian = NewtonianAltitudeGravity()
    assert constant.surface_value() == pytest.approx(newtonian.surface_value())
    assert constant.parameters() != newtonian.parameters()
    assert constant.parameters()["gravity_model"] == "constant"
    assert newtonian.parameters()["gravity_model"] == "newtonian_altitude"

    # A bare number still works, because a constant is the first rung of the
    # fidelity path rather than a way around it -- it just arrives declared now.
    assert isinstance(as_gravity_field(9.5), ConstantGravity)
    assert as_gravity_field(None).surface_value() == pytest.approx(G0)
    assert as_gravity_field(newtonian) is newtonian


@pytest.mark.req("GRV-001", "DAT-007")
def test_a_height_dependent_field_refuses_to_be_evaluated_without_heights():
    """Declining is not falling back, which is the difference that matters.

    Silently using the surface value would give a plausible answer and would
    undo the whole point: the approximation would be back to being applied
    without anyone declaring it.
    """
    with pytest.raises(ConfigurationError, match="varies with height"):
        hydrostatic_layer_masses(
            n_layers=8,
            surface_pressure=1.0e5,
            top_pressure=1.0e3,
            gravity=NewtonianAltitudeGravity(),
        )
    with pytest.raises(ConfigurationError, match="entries for 8 layers"):
        hydrostatic_layer_masses(
            n_layers=8,
            surface_pressure=1.0e5,
            top_pressure=1.0e3,
            gravity=NewtonianAltitudeGravity(),
            layer_heights=np.zeros(4),
        )


@pytest.mark.req("GRV-001", "EQR-004", "CON-001")
def test_what_the_constant_approximation_actually_costs_this_column():
    """The measurement the equation record was asserting instead of making.

    It said constant gravity is "negligible through the troposphere". True, and
    about a column this framework does not run: the top at 1000 Pa is 19 km up at
    8 layers. Held against the inverse square law over the column's own heights,
    the constant puts 0.21 percent too little mass in the column and 0.6 percent
    too little in the topmost layer -- which is the layer that sets the outgoing
    longwave, so it is the least forgiving place for the error to sit.
    """
    n_layers = 8
    p_s, p_top = 1.0e5, 1.0e3
    interfaces = np.linspace(p_s, p_top, n_layers + 1)
    p_mid = 0.5 * (interfaces[:-1] + interfaces[1:])
    heights = hydrostatic_heights(np.full(n_layers, 250.0), p_mid, p_s)

    _, constant = hydrostatic_layer_masses(
        n_layers=n_layers, surface_pressure=p_s, top_pressure=p_top, gravity=ConstantGravity()
    )
    _, newtonian = hydrostatic_layer_masses(
        n_layers=n_layers,
        surface_pressure=p_s,
        top_pressure=p_top,
        gravity=NewtonianAltitudeGravity(),
        layer_heights=heights,
    )

    # The constant case stays exactly (p_s - p_top) / g, which is what makes the
    # column's mass budget checkable by hand.
    assert constant.sum() == pytest.approx((p_s - p_top) / G0, rel=1e-12)
    assert np.all(np.diff(constant) == pytest.approx(0.0))  # equal layers

    assert newtonian.sum() / constant.sum() - 1.0 == pytest.approx(0.00214, abs=2e-5)
    assert newtonian[-1] / constant[-1] - 1.0 == pytest.approx(0.0061, abs=2e-4)
    # Every layer gains, and the higher it is the more: gravity only weakens.
    gain = newtonian / constant - 1.0
    assert np.all(gain > 0.0)
    assert np.all(np.diff(gain) > 0.0)
    assert heights[-1] > 19.0e3, "the column's top is not in the troposphere"


@pytest.mark.req("GRV-001", "ATM-002")
def test_the_column_takes_a_field_and_the_default_is_unchanged():
    """The wiring, and the guarantee that adding the interface moved nothing.

    Every number this project has recorded was produced at constant gravity, so
    the default has to still be exactly that -- not approximately.
    """
    default = AtmosphericColumnModule(n_layers=8)
    explicit = AtmosphericColumnModule(n_layers=8, gravity=G0)
    declared = AtmosphericColumnModule(n_layers=8, gravity=ConstantGravity())
    assert default.layer_mass == pytest.approx(explicit.layer_mass, rel=0.0, abs=0.0)
    assert default.layer_mass == pytest.approx(declared.layer_mass, rel=0.0, abs=0.0)
    assert default.gravity == pytest.approx(G0)

    varying = AtmosphericColumnModule(n_layers=8, gravity=NewtonianAltitudeGravity())
    assert varying.layer_mass.sum() > default.layer_mass.sum()
    assert varying.parameters()["gravity_model"] == "newtonian_altitude"
    assert default.parameters()["gravity_model"] == "constant"
