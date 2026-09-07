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

from mfespf.coupling.validity import check_validity_domains
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


@dataclass
class _Module:
    contract: _Contract


@dataclass
class _Record:
    record_id: str
    validity_domain: ValidityDomain


class _Registry:
    def __init__(self, *records: _Record) -> None:
        self._records = {r.record_id: r for r in records}

    def current(self) -> list[_Record]:
        return list(self._records.values())


def _report(record, contract, fields):
    return check_validity_domains(
        registry=_Registry(record), modules=[_Module(contract)], fields=fields
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
