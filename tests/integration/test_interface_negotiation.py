"""VER-006 and VER-008 in running configurations, not in isolation.

``negotiate`` had unit tests for years and no caller. These tests are about the
caller: that a real configuration negotiates every declared pair before any
module initialises, that a disagreement stops it there, and that the one way
across a gap -- a declared migration adapter -- says what it costs and leaves
that statement in the manifest.
"""

from __future__ import annotations

import dataclasses

import pytest

from mfespf.core.errors import ConfigurationError, InterfaceVersionError
from mfespf.core.versioning import MigrationAdapter, Version, VersionRange
from mfespf.run import column_configurations as cc
from mfespf.run.driver import Kernel

pytestmark = pytest.mark.stage(2)


def _atmosphere(kernel: Kernel):
    (module,) = [
        s.module
        for s in kernel.schedule.ordered()
        if any(i.interface_id == "atmospheric_state" for i in s.module.contract.provided_interfaces)
    ]
    return module


def _republish(kernel: Kernel, **changes) -> None:
    """Replace the atmosphere's provided interface, as a new release of it would."""
    module = _atmosphere(kernel)
    (old,) = [
        i for i in module.contract.provided_interfaces if i.interface_id == "atmospheric_state"
    ]
    semantics = {**old.physical_semantics, **changes.pop("semantics", {})}
    new = dataclasses.replace(old, physical_semantics=semantics, **changes)
    module.contract = dataclasses.replace(module.contract, provided_interfaces=(new,))


def _initialisations(kernel: Kernel) -> list[str]:
    """Record every initialise call, so a test can show none happened."""
    calls: list[str] = []
    for slot in kernel.schedule.ordered():
        module = slot.module
        original = module.initialize

        def wrapper(*args, _m=module, _f=original, **kwargs):
            calls.append(_m.module_id)
            return _f(*args, **kwargs)

        module.initialize = wrapper
    return calls


@pytest.mark.req("VER-006", "VER-005")
def test_every_declared_requirement_is_negotiated_and_recorded():
    kernel = cc.moist_radiative_convective()
    kernel.run(1)
    declared = sum(
        len(s.module.contract.required_interfaces)
        for s in kernel.schedule.ordered()
        if s.module.status.participates
    )
    assert declared > 0
    assert kernel.coupler is not None
    assert len(kernel.coupler.negotiations) == declared
    assert kernel.manifest is not None
    assert len(kernel.manifest.negotiated_interfaces) == declared
    assert all(v.startswith("direct v") for v in kernel.manifest.negotiated_interfaces.values())


@pytest.mark.req("VER-006", "VER-007")
def test_a_reversed_layer_ordering_stops_the_run_before_anything_initialises():
    """The failure the ``expects`` declarations exist for.

    Same interface, same version, and the column upside down. Before
    negotiation was wired in, this ran.
    """
    kernel = cc.radiative_convective()
    _republish(kernel, semantics={"layer_ordering": "index_0_is_the_top_layer"})
    calls = _initialisations(kernel)
    with pytest.raises(InterfaceVersionError, match="layer_ordering"):
        kernel.run(1)
    assert calls == []
    assert kernel.manifest is not None
    assert kernel.manifest.gate_results["configuration_validation"] is False


@pytest.mark.req("VER-006")
def test_a_requirement_nobody_provides_stops_the_run_before_anything_initialises():
    kernel = cc.radiative_convective()
    module = _atmosphere(kernel)
    module.contract = dataclasses.replace(module.contract, provided_interfaces=())
    calls = _initialisations(kernel)
    with pytest.raises(InterfaceVersionError, match="no participating module provides it"):
        kernel.run(1)
    assert calls == []


def _major_bump_adapter(**overrides) -> MigrationAdapter:
    declared = dict(
        interface_id="atmospheric_state",
        from_range=VersionRange.parse(">=2.0.0,<3.0.0"),
        to_version=Version(1, 0, 0),
        approximation="none: the 2.x release renamed fields and changed no values",
        information_loss="none",
        conservation_residual="none",
    )
    declared.update(overrides)
    return MigrationAdapter(**declared)


@pytest.mark.req("VER-008", "VER-006")
def test_a_declared_adapter_bridges_a_major_version_and_its_cost_reaches_the_manifest():
    refused = cc.radiative_convective()
    _republish(refused, version=Version(2, 0, 0))
    with pytest.raises(InterfaceVersionError, match="No migration adapter covers this gap"):
        refused.run(1)

    kernel = cc.radiative_convective()
    _republish(kernel, version=Version(2, 0, 0))
    kernel.migration_adapters = (_major_bump_adapter(),)
    result = kernel.run(2)
    assert result.manifest is kernel.manifest
    via = [v for v in kernel.manifest.negotiated_interfaces.values() if not v.startswith("direct")]
    assert via, "the adapted pairs must say so in the manifest"
    for entry in via:
        assert "approximation: none: the 2.x release renamed fields" in entry
        assert "information loss: none" in entry
        assert "conservation residual: none" in entry


@pytest.mark.req("VER-008", "VER-007")
def test_an_adapter_cannot_carry_a_semantic_change_it_does_not_name():
    """An adapter written for a units change must not paper over a flipped column."""
    kernel = cc.radiative_convective()
    _republish(
        kernel,
        version=Version(2, 0, 0),
        semantics={"layer_ordering": "index_0_is_the_top_layer", "pressure_units": "hPa"},
    )
    kernel.migration_adapters = (
        _major_bump_adapter(
            translates=frozenset({"pressure_units"}),
            approximation="none: hPa to Pa is an exact factor of 100",
        ),
    )
    calls = _initialisations(kernel)
    with pytest.raises(InterfaceVersionError, match=r"does not translate \['layer_ordering'\]"):
        kernel.run(1)
    assert calls == []


@pytest.mark.req("VER-008")
@pytest.mark.parametrize("omitted", ["approximation", "information_loss", "conservation_residual"])
def test_an_adapter_that_does_not_declare_its_cost_is_refused(omitted):
    with pytest.raises(ConfigurationError, match=omitted):
        _major_bump_adapter(**{omitted: "  "})


@pytest.mark.req("VER-008")
def test_a_nonzero_residual_cannot_be_described_as_none():
    with pytest.raises(ConfigurationError, match="describes it as 'none'"):
        _major_bump_adapter(residual_magnitude=1e-6)
    declared = _major_bump_adapter(
        residual_magnitude=1e-6,
        conservation_residual="column mass, from re-gridding 2.x's pressure levels",
    )
    assert "1e-06 (column mass" in declared.describe()
