"""Run-time enforcement of declared equation validity domains.

Charter 4.26; requirements EQR-002, EQR-003, SYS-006.

Every equation record carried a validity domain and nothing checked it, so a
relation could be evaluated far outside the conditions it was established for
and return a plausible number in silence. These tests are about the check that
closed that, and specifically about the two ways such a check can lie: by
matching the wrong quantity, and by counting a variable it could not find as a
pass.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pytest

from mfespf.coupling.validity import (
    INTEGRATED_OVER,
    EvaluatedExtents,
    check_validity_domains,
)
from mfespf.registries.equations import ValidityDomain

pytestmark = pytest.mark.stage(3)


@dataclass
class _Field:
    symbol: str
    units: str
    values: Any


@dataclass
class _Contract:
    module_id: str
    governing_equation_ids: tuple[str, ...]
    required_inputs: tuple[str, ...] = ()
    produced_outputs: tuple[str, ...] = ()
    equation_variables: tuple[tuple[str, tuple[tuple[str, str, str], ...]], ...] = ()


@dataclass
class _Module:
    contract: _Contract
    _parameters: dict[str, Any] | None = None

    def parameters(self) -> dict[str, Any]:
        return dict(self._parameters or {})


@dataclass
class _Record:
    record_id: str
    validity_domain: ValidityDomain


class _Registry:
    def __init__(self, *records: _Record) -> None:
        self._records = {r.record_id: r for r in records}

    def current(self) -> list[_Record]:
        return list(self._records.values())


def _report(record, contract, fields, parameters=None):
    return check_validity_domains(
        registry=_Registry(record), modules=[_Module(contract, parameters)], fields=fields
    )


@pytest.mark.req("EQR-003", "SYS-006")
def test_a_value_outside_its_declared_range_is_reported():
    record = _Record("EQ-X", ValidityDomain(ranges=(("T", 230.0, 340.0, "K"),)))
    contract = _Contract("water", ("EQ-X",), produced_outputs=("layer_temperature",))
    fields = {"layer_temperature": _Field("T_a", "K", np.array([172.0, 300.0]))}

    report = _report(record, contract, fields)

    assert not report.passed
    assert len(report.violations) == 1
    assert "outside the declared [230, 340]" in report.violations[0].detail
    # And the same range with an in-domain profile passes.
    inside = _report(
        record, contract, {"layer_temperature": _Field("T_a", "K", np.array([250.0, 300.0]))}
    )
    assert inside.passed and not inside.unchecked


@pytest.mark.req("EQR-003")
def test_a_symbol_is_never_matched_outside_the_citing_module():
    """``alpha`` is an evaporating fraction here and a surface albedo there.

    A global match by symbol would check one against the other's range. The
    module's own declared fields are the whole scope, so a field it does not
    declare is invisible however well the symbol fits.
    """
    record = _Record("EQ-PLUME", ValidityDomain(ranges=(("alpha", 0.6, 0.85, "1"),)))
    contract = _Contract("water", ("EQ-PLUME",), produced_outputs=("humidity",))
    fields = {
        "humidity": _Field("RH", "1", np.array([0.7])),
        # Squarely outside the plume's range, and none of the plume's business.
        "surface_spectral_albedo": _Field("alpha", "1", np.array([0.0, 0.3])),
    }

    report = _report(record, contract, fields)

    assert not report.violations, "matched a symbol the module never declared"
    assert len(report.unchecked) == 1


@pytest.mark.req("EQR-003")
def test_an_unresolvable_variable_is_never_counted_as_a_pass():
    """The failure mode that would make the whole check worthless."""
    record = _Record("EQ-Y", ValidityDomain(ranges=(("tau", 0.0, 10.0, "1"),)))
    contract = _Contract("radiation", ("EQ-Y",), produced_outputs=("heating",))
    fields = {"heating": _Field("H", "W m^-2", np.array([1.0]))}

    report = _report(record, contract, fields)

    assert len(report.unchecked) == 1
    assert not report.violations
    # passed() is about what was checkable; the coverage is reported separately
    # so that "nothing violated" cannot be read as "everything verified".
    assert report.passed
    assert "1 unresolvable" in report.summary()


@pytest.mark.req("EQR-003", "SYS-006")
def test_ambiguity_is_refused_rather_than_guessed():
    """Two candidates of the same root and units means the equation is silent."""
    record = _Record("EQ-Z", ValidityDomain(ranges=(("T", 230.0, 340.0, "K"),)))
    contract = _Contract(
        "radiation", ("EQ-Z",), required_inputs=("air", "ground"), produced_outputs=()
    )
    fields = {
        "air": _Field("T_a", "K", np.array([172.0])),  # would violate
        "ground": _Field("T_s", "K", np.array([288.0])),  # would pass
    }

    report = _report(record, contract, fields)

    assert not report.violations
    assert len(report.unchecked) == 1
    assert "could be any of" in report.unchecked[0].detail


@pytest.mark.req("EQR-003")
def test_a_range_in_other_units_is_not_silently_compared():
    """A degC range and a kelvin field differ by 273.15 and neither says so."""
    record = _Record("EQ-U", ValidityDomain(ranges=(("T", -40.0, 60.0, "degC"),)))
    contract = _Contract("surface", ("EQ-U",), produced_outputs=("temperature",))
    fields = {"temperature": _Field("T", "K", np.array([288.0]))}

    report = _report(record, contract, fields)

    assert not report.violations
    assert len(report.unchecked) == 1
    assert "degC" in report.unchecked[0].detail


@pytest.mark.req("EQR-003", "SYS-006")
def test_every_binding_on_one_symbol_is_judged_not_just_the_first():
    """One symbol can name several quantities, and checking one is checking none.

    The cloud optics declare a range on ``W`` and the module carries a liquid
    path and an ice path; the water module evaluates saturation at the layer
    temperature and again at the surface. Keeping only the first binding would
    leave the other quantity evaluated and silently unjudged, which is the
    half-answer this check refuses when it refuses an ambiguity.
    """
    record = _Record("EQ-W", ValidityDomain(ranges=(("W", 0.0, 10.0, "kg m^-2"),)))
    contract = _Contract(
        "cloud",
        ("EQ-W",),
        produced_outputs=("liquid", "ice"),
        equation_variables=(("EQ-W", (("W", "liquid", "kg m^-2"), ("W", "ice", "kg m^-2"))),),
    )
    fields = {
        "liquid": _Field("LWP", "kg m^-2", np.array([0.0, 0.05])),
        "ice": _Field("IWP", "kg m^-2", np.array([0.0, 40.0])),  # outside
    }

    report = _report(record, contract, fields)

    assert len(report.checks) == 2
    assert len(report.violations) == 1
    assert "ice" in report.violations[0].detail


@pytest.mark.req("EQR-003")
def test_a_binding_beats_a_unique_but_wrong_inferred_match():
    """A single candidate is not evidence of a correct one.

    The ambiguity refusal only fires when a second candidate exists to notice.
    The deep ocean's slab balance was matched against the *mixed layer's* heat
    capacity, which was the only field of that symbol and units in scope, and
    so was reported on confidently. Only the module can say otherwise.
    """
    record = _Record("EQ-C", ValidityDomain(ranges=(("C", 1.0e3, 1.0e10, "J m^-2 K^-1"),)))
    inferred = _Contract("ocean", ("EQ-C",), required_inputs=("surface_heat_capacity",))
    fields = {"surface_heat_capacity": _Field("C", "J m^-2 K^-1", np.array([1.0e6]))}

    # Left to inference it resolves -- to the wrong reservoir.
    assert "surface_heat_capacity" in _report(record, inferred, fields).checks[0].detail

    bound = _Contract(
        "ocean",
        ("EQ-C",),
        required_inputs=("surface_heat_capacity",),
        equation_variables=(("EQ-C", (("C", "deep_capacity", "J m^-2 K^-1"),)),),
    )
    check = _report(record, bound, fields, {"deep_capacity": 4.1e8}).checks[0]

    assert check.status == "inside"
    assert "deep_capacity" in check.detail
    assert "surface_heat_capacity" not in check.detail


@pytest.mark.req("EQR-003", "SYS-006")
def test_a_masked_value_narrows_the_check_to_where_the_relation_was_applied():
    """NaN means "not evaluated here", which is the standing coarseness's cure.

    ``E = L_v sum_i w_i`` was flagged for reaching 225 K against a floor of 230,
    and the column holds exactly zero vapour below 230 K -- every temperature
    that broke the range multiplied a weight of zero. Masking the layers that
    carry no vapour turns the extent a quantity *spans* into the region the
    equation was *applied* to.
    """
    record = _Record("EQ-L", ValidityDomain(ranges=(("T", 230.0, 340.0, "K"),)))
    contract = _Contract("water", ("EQ-L",), produced_outputs=("T_where_vapour",))
    unmasked = {"T_where_vapour": _Field("T", "K", np.array([259.6, 232.7, 225.2]))}
    masked = {"T_where_vapour": _Field("T", "K", np.array([259.6, np.nan, np.nan]))}

    assert _report(record, contract, unmasked).violations
    report = _report(record, contract, masked)
    assert not report.violations
    assert report.checks[0].status == "inside"
    assert "259.6" in report.checks[0].detail


@pytest.mark.req("EQR-003")
def test_an_infinity_is_a_defect_rather_than_a_mask():
    """Masking is a statement a module makes; an infinity is a bug.

    If a non-finite value simply meant "skip me", a field that had blown up
    could pass the check by being broken.
    """
    record = _Record("EQ-I", ValidityDomain(ranges=(("T", 230.0, 340.0, "K"),)))
    contract = _Contract("water", ("EQ-I",), produced_outputs=("temperature",))
    fields = {"temperature": _Field("T", "K", np.array([288.0, np.inf]))}

    report = _report(record, contract, fields)

    assert not report.violations
    assert len(report.unchecked) == 1
    assert "finite" in report.unchecked[0].detail


@pytest.mark.req("EQR-003")
def test_a_variable_integrated_over_is_declared_rather_than_left_unfound():
    """Two unresolvables that mean different things should not read alike.

    The surface's thermal emission is Planck integrated over every wavelength,
    so its ``lambda`` has no value a range could be checked against. That is not
    a gap to close -- closing it would state a wavelength support for a quantity
    that has none left.
    """
    record = _Record("EQ-P", ValidityDomain(ranges=(("lambda", 1e-9, 1e-2, "m"),)))
    contract = _Contract(
        "surface",
        ("EQ-P",),
        produced_outputs=("surface_temperature",),
        equation_variables=(("EQ-P", (("lambda", INTEGRATED_OVER, "m"),)),),
    )
    fields = {"surface_temperature": _Field("T_s", "K", np.array([288.0]))}

    report = _report(record, contract, fields)

    assert not report.violations
    assert len(report.unchecked) == 1
    assert "integrates over this variable" in report.unchecked[0].detail


def _over_run(record, contract, history, parameters=None):
    """Judge a whole trajectory the way the driver does."""
    extents = EvaluatedExtents()
    for fields in history:
        extents.observe(fields)
    return check_validity_domains(
        registry=_Registry(record),
        modules=[_Module(contract, parameters)],
        fields=history[-1],
        extents=extents,
    )


@pytest.mark.req("EQR-003", "SYS-006")
def test_a_relation_that_left_its_domain_and_came_back_is_still_caught():
    """The failure that judging a final state cannot see.

    A run's last step is not where its evidence lives. A closure that ran far
    outside its declared conditions during spin-up and settled back inside them
    was still evaluated outside them, and every number it produced on the way is
    in the answer.
    """
    record = _Record("EQ-X", ValidityDomain(ranges=(("T", 230.0, 340.0, "K"),)))
    contract = _Contract("water", ("EQ-X",), produced_outputs=("layer_temperature",))
    history = [
        {"layer_temperature": _Field("T_a", "K", np.array([288.0]))},
        {"layer_temperature": _Field("T_a", "K", np.array([201.4]))},  # out, mid-run
        {"layer_temperature": _Field("T_a", "K", np.array([291.0]))},
    ]

    assert _report(record, contract, history[-1]).passed, "the final state alone is innocent"

    report = _over_run(record, contract, history)

    assert not report.passed
    assert "over the run" in report.violations[0].detail
    assert "201.4" in report.violations[0].detail


@pytest.mark.req("EQR-003", "SYS-006", "DAT-007")
def test_a_closure_idle_at_the_end_is_judged_on_the_steps_it_ran():
    """Why the mask can now say what it means.

    The convective adjustment mixes on 79 percent of steps and on none of the
    first four hundred, so a per-step mask left nothing behind on a final state
    and the identical configuration was resolvable at 2000 steps and
    unresolvable at 400. Accumulating removes the incentive to widen the mask to
    the union over the run, which would report a spin-up layer at whatever it
    holds when the run stops.
    """
    record = _Record("EQ-M", ValidityDomain(ranges=(("T", 230.0, 340.0, "K"),)))
    contract = _Contract("convection", ("EQ-M",), produced_outputs=("mixed_layer_temperature",))
    ran = {"mixed_layer_temperature": _Field("T", "K", np.array([301.0, 288.0]))}
    idle = {"mixed_layer_temperature": _Field("T", "K", np.array([np.nan, np.nan]))}

    # The closure was active, then stopped. The final state says nothing at all.
    assert _report(record, contract, idle).unchecked

    report = _over_run(record, contract, [ran, idle, idle])

    assert report.checks[0].status == "inside"
    assert "301" in report.checks[0].detail

    # And a closure that never ran at any step stays unchecked, rather than
    # being credited with a pass it did not earn.
    assert _over_run(record, contract, [idle, idle]).unchecked


@pytest.mark.req("EQR-003")
def test_an_infinity_at_any_step_disqualifies_the_field_for_the_whole_run():
    """A field that blew up once must not be rescued by recovering.

    An infinity is a defect rather than a mask, and accumulation must not turn
    it into one: taking a running minimum and maximum over finite values alone
    would quietly drop the step that was broken.
    """
    record = _Record("EQ-I", ValidityDomain(ranges=(("T", 230.0, 340.0, "K"),)))
    contract = _Contract("water", ("EQ-I",), produced_outputs=("temperature",))
    blew_up = {"temperature": _Field("T", "K", np.array([288.0, np.inf]))}
    fine = {"temperature": _Field("T", "K", np.array([288.0, 290.0]))}

    report = _over_run(record, contract, [fine, blew_up, fine])

    assert not report.violations
    assert len(report.unchecked) == 1
    assert "finite" in report.unchecked[0].detail
