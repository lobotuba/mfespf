"""Stage 0: semantic interface versions and negotiation."""

from __future__ import annotations

import pytest

from mfespf.core.errors import ConfigurationError, InterfaceVersionError
from mfespf.core.versioning import InterfaceContract, Version, VersionRange, negotiate

pytestmark = pytest.mark.stage(0)


def _contract(**semantics) -> InterfaceContract:
    return InterfaceContract(
        interface_id="radiation",
        version=Version(1, 2, 0),
        physical_semantics={
            "flux_units": "W m^-2",
            "sign_convention": "positive_downward",
            "averaging": "interval_average",
            **semantics,
        },
        optional_capabilities=frozenset({"spectral_output"}),
    )


@pytest.mark.req("VER-001", "VER-002")
def test_version_parsing_and_ordering():
    assert Version.parse("1.2.3") == Version(1, 2, 3)
    assert Version(1, 2, 3) < Version(1, 3, 0) < Version(2, 0, 0)
    with pytest.raises(ConfigurationError):
        Version.parse("1.2")


@pytest.mark.req("VER-002", "VER-003", "VER-004")
def test_compatibility_follows_semantic_version_rules():
    assert Version(1, 3, 0).compatible_with(Version(1, 2, 0))  # newer minor is fine
    assert Version(1, 2, 5).compatible_with(Version(1, 2, 0))  # newer patch is fine
    assert not Version(1, 1, 0).compatible_with(Version(1, 2, 0))  # older minor is not
    assert not Version(2, 0, 0).compatible_with(Version(1, 2, 0))  # major break


@pytest.mark.req("VER-005", "VER-006", "KER-003")
def test_negotiation_succeeds_within_the_declared_range():
    provider = _contract()
    resolved = negotiate(
        provider=provider,
        consumer_interface_id="radiation",
        required_range=VersionRange.parse(">=1.0.0,<2.0.0"),
    )
    assert resolved.version == Version(1, 2, 0)


@pytest.mark.req("VER-006")
def test_negotiation_fails_before_initialization_on_a_major_mismatch():
    with pytest.raises(InterfaceVersionError, match="does not satisfy required range"):
        negotiate(
            provider=_contract(),
            consumer_interface_id="radiation",
            required_range=VersionRange.parse(">=2.0.0,<3.0.0"),
        )


@pytest.mark.req("VER-007")
def test_a_changed_sign_convention_cannot_hide_inside_a_compatible_version():
    """The failure VER-007 exists for: same version, opposite meaning."""
    provider = _contract(sign_convention="positive_upward")
    with pytest.raises(InterfaceVersionError, match="disagree on physical meaning"):
        negotiate(
            provider=provider,
            consumer_interface_id="radiation",
            required_range=VersionRange.parse(">=1.0.0,<2.0.0"),
            consumer_semantics={"sign_convention": "positive_downward"},
        )


@pytest.mark.req("VER-007")
def test_physical_signature_changes_when_meaning_changes():
    assert _contract().physical_signature() != _contract(flux_units="W").physical_signature()


@pytest.mark.req("VER-005")
def test_a_missing_optional_capability_is_a_negotiation_failure():
    with pytest.raises(InterfaceVersionError, match="lacks required optional capabilities"):
        negotiate(
            provider=_contract(),
            consumer_interface_id="radiation",
            required_range=VersionRange.parse(">=1.0.0"),
            required_capabilities=frozenset({"line_by_line"}),
        )
