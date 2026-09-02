"""Semantic interface versions and negotiation.

Amendment A6 and requirements VER-001 through VER-008, plus KER-003.

The rule that gives this module its teeth is A6.4 / VER-007: a change in
physical meaning -- units, sign convention, averaging semantics, spatial support,
state ownership, or spectral-density basis -- may not hide inside a compatible
patch.  ``InterfaceContract.physical_signature`` makes that mechanical: the
signature is hashed from exactly those attributes, and a mismatch under a
compatible version pair is a hard error rather than a warning.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from dataclasses import field as dc_field
from functools import total_ordering

from mfespf.core.errors import ConfigurationError, InterfaceVersionError

__all__ = ["InterfaceContract", "Version", "VersionRange", "negotiate"]

_VERSION_RE = re.compile(r"^(\d+)\.(\d+)\.(\d+)$")


@total_ordering
@dataclass(frozen=True, slots=True)
class Version:
    """MAJOR.MINOR.PATCH with the meanings fixed by A6.2."""

    major: int
    minor: int
    patch: int

    @classmethod
    def parse(cls, text: str) -> Version:
        m = _VERSION_RE.match(text.strip())
        if not m:
            raise ConfigurationError(
                f"invalid semantic version {text!r}; expected MAJOR.MINOR.PATCH"
            )
        return cls(int(m.group(1)), int(m.group(2)), int(m.group(3)))

    def __str__(self) -> str:
        return f"{self.major}.{self.minor}.{self.patch}"

    def __lt__(self, other: Version) -> bool:
        return (self.major, self.minor, self.patch) < (other.major, other.minor, other.patch)

    def compatible_with(self, required: Version) -> bool:
        """True when this version can satisfy a requirement for ``required``.

        Same major, and at least the required minor.  A newer minor is
        backward-compatible by A6.2; a different major is not.
        """
        if self.major != required.major:
            return False
        if self.minor != required.minor:
            return self.minor > required.minor
        return self.patch >= required.patch


@dataclass(frozen=True, slots=True)
class VersionRange:
    """An inclusive minimum and exclusive maximum, as declared by VER-005."""

    minimum: Version
    below: Version | None = None

    @classmethod
    def parse(cls, text: str) -> VersionRange:
        """Accept ``1.2.0``, ``>=1.2.0``, or ``>=1.2.0,<2.0.0``."""
        parts = [p.strip() for p in text.split(",") if p.strip()]
        minimum: Version | None = None
        below: Version | None = None
        for part in parts:
            if part.startswith(">="):
                minimum = Version.parse(part[2:])
            elif part.startswith("<"):
                below = Version.parse(part[1:])
            else:
                minimum = Version.parse(part)
                below = Version(minimum.major + 1, 0, 0)
        if minimum is None:
            raise ConfigurationError(f"version range {text!r} declares no minimum")
        return cls(minimum, below)

    def contains(self, version: Version) -> bool:
        if version < self.minimum:
            return False
        if self.below is not None and not version < self.below:
            return False
        return True

    def __str__(self) -> str:
        return f">={self.minimum}" + (f",<{self.below}" if self.below else "")


@dataclass(frozen=True)
class InterfaceContract:
    """A versioned interface a module provides or requires.

    ``physical_semantics`` holds the attributes whose change constitutes a change
    in physical meaning (VER-007).  Anything placed here participates in the
    signature hash; anything outside it is free to vary across patch versions.
    """

    interface_id: str
    version: Version
    physical_semantics: dict[str, str] = dc_field(default_factory=dict)
    optional_capabilities: frozenset[str] = frozenset()
    deprecated_fields: frozenset[str] = frozenset()

    def physical_signature(self) -> str:
        """Stable hash of the meaning-bearing attributes of this contract."""
        payload = "|".join(
            f"{k}={self.physical_semantics[k]}" for k in sorted(self.physical_semantics)
        )
        return hashlib.sha256(f"{self.interface_id}|{payload}".encode()).hexdigest()[:16]


def negotiate(
    *,
    provider: InterfaceContract,
    consumer_interface_id: str,
    required_range: VersionRange,
    consumer_semantics: dict[str, str] | None = None,
    required_capabilities: frozenset[str] = frozenset(),
    context: str = "",
) -> InterfaceContract:
    """Resolve one provider/consumer interface pair, or fail before initialisation.

    Implements VER-006 and A6.5: a run must fail here rather than start with an
    interface whose meaning the two sides read differently.
    """
    where = f" ({context})" if context else ""

    if provider.interface_id != consumer_interface_id:
        raise InterfaceVersionError(
            f"interface identity mismatch{where}: provider offers "
            f"{provider.interface_id!r}, consumer requires {consumer_interface_id!r}"
        )

    if not required_range.contains(provider.version):
        raise InterfaceVersionError(
            f"interface {provider.interface_id!r}{where}: provider version "
            f"{provider.version} does not satisfy required range {required_range}"
        )

    missing = required_capabilities - provider.optional_capabilities
    if missing:
        raise InterfaceVersionError(
            f"interface {provider.interface_id!r}{where}: provider lacks required "
            f"optional capabilities {sorted(missing)}"
        )

    if consumer_semantics:
        disagreements = {
            key: (provider.physical_semantics.get(key), value)
            for key, value in consumer_semantics.items()
            if provider.physical_semantics.get(key) != value
        }
        if disagreements:
            detail = ", ".join(
                f"{k}: provider={p!r} consumer={c!r}" for k, (p, c) in sorted(disagreements.items())
            )
            raise InterfaceVersionError(
                f"interface {provider.interface_id!r} v{provider.version}{where}: the two sides "
                f"disagree on physical meaning ({detail}). VER-007 forbids resolving this by "
                f"version compatibility alone."
            )

    return provider
