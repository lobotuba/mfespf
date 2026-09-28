"""Semantic interface versions and negotiation.

Amendment A6 and requirements VER-001 through VER-008, plus KER-003.

The rule that gives this module its teeth is A6.4 / VER-007: a change in
physical meaning -- units, sign convention, averaging semantics, spatial support,
state ownership, or spectral-density basis -- may not hide inside a compatible
patch.  ``InterfaceContract.physical_signature`` makes that mechanical: the
signature is hashed from exactly those attributes, and a mismatch under a
compatible version pair is a hard error rather than a warning.

Declared and then never checked
-------------------------------

For a long time that was the whole of it, and it had the defect this project
keeps finding: ``negotiate`` was called from this module's own tests and from no
running code. Every module declared what it provides and what it requires -- in
the moist column, eight provided interfaces carrying forty-four semantic keys
against seven consumer declarations -- and the framework matched none of them up.
It was the validity-domain situation exactly: the conditions were written down
and nothing compared a run against them.

``Coupler.validate_configuration`` negotiates now, before initialisation, which
is what VER-006 asks for. And because a consumer that states only an interface id
and a version range can only be checked for identity and version, a
``RequiredInterface`` also carries the semantics the consumer *assumes* --
``layer_ordering=index_0_is_the_lowest_layer``,
``flux_sign_convention=positive_upward_away_from_surface``. Those are the
statements a silent disagreement would reverse a column or a flux by.

Crossing a gap
--------------

Refusing is the right default and cannot be the only option, because the way
round a refusal that cannot be satisfied is to loosen the refusal. So there is
one declared path across a version or semantics gap: a ``MigrationAdapter``,
which must state the approximation it introduces, the information it loses and
the conservation residual it leaves (VER-008). Each is required, and "none" is a
valid declaration that has to be *made* -- an empty string is refused, the same
way an unquantified uncertainty must be marked UNRESOLVED rather than omitted.
"""

from __future__ import annotations

import hashlib
import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from dataclasses import field as dc_field
from functools import total_ordering

from mfespf.core.errors import ConfigurationError, InterfaceVersionError

__all__ = [
    "InterfaceContract",
    "MigrationAdapter",
    "Negotiation",
    "RequiredInterface",
    "Version",
    "VersionRange",
    "negotiate",
    "negotiate_configuration",
]

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


@dataclass(frozen=True)
class RequiredInterface:
    """What a consumer needs, including what it assumes the numbers mean.

    An interface id and a version range can only be checked for identity and
    version, and neither of those changes when a provider reverses its layer
    ordering inside a compatible minor release. ``expects`` is where the consumer
    writes down the assumptions its code is built on, so VER-007's rule -- a
    physical change shall not hide inside a compatible patch -- has something to
    compare against.

    Only the keys a consumer actually depends on belong here. Copying the
    provider's whole semantics block would make every provider change a
    consumer change and would say nothing about which assumptions are load
    bearing.
    """

    interface_id: str
    version_range: VersionRange
    expects: dict[str, str] = dc_field(default_factory=dict)
    required_capabilities: frozenset[str] = frozenset()


@dataclass(frozen=True)
class MigrationAdapter:
    """The one declared way across a version or semantics gap (VER-008).

    Refusing an incompatible pair is right and cannot be the only option: the way
    round a refusal that cannot be satisfied is to loosen the refusal, and
    ``negotiate`` is not a thing to loosen. So an adapter may bridge a gap, and
    the price is that it says what the bridge costs.

    All three declarations are required, and "none" is a valid one that has to be
    *made*. An empty string is refused, on the same reasoning as UNC-004: a
    silent omission and a declared absence must not look alike.

    ``translates`` names the semantic keys this adapter converts. A key the
    adapter does not name is a key it cannot bridge, so an adapter written for a
    units change cannot quietly paper over a sign-convention change that arrived
    in the same release.
    """

    interface_id: str
    from_range: VersionRange
    to_version: Version
    translates: frozenset[str] = frozenset()
    approximation: str = ""
    information_loss: str = ""
    conservation_residual: str = ""
    residual_magnitude: float = 0.0

    def __post_init__(self) -> None:
        missing = [
            name
            for name in ("approximation", "information_loss", "conservation_residual")
            if not str(getattr(self, name)).strip()
        ]
        if missing:
            raise ConfigurationError(
                f"migration adapter for {self.interface_id!r} does not declare "
                f"{missing} (VER-008). Write 'none' if there is none; an empty "
                f"declaration and an absent one must not look alike."
            )
        if self.residual_magnitude < 0.0:
            raise ConfigurationError("a conservation residual cannot be negative")
        if self.residual_magnitude > 0.0 and self.conservation_residual.strip() == "none":
            raise ConfigurationError(
                f"migration adapter for {self.interface_id!r} declares a residual of "
                f"{self.residual_magnitude:g} and describes it as 'none'"
            )

    def covers(self, provider: InterfaceContract, required: RequiredInterface) -> bool:
        return (
            provider.interface_id == self.interface_id == required.interface_id
            and self.from_range.contains(provider.version)
            and required.version_range.contains(self.to_version)
        )

    def describe(self) -> str:
        residual = (
            f"{self.residual_magnitude:g} ({self.conservation_residual})"
            if self.residual_magnitude
            else self.conservation_residual
        )
        return (
            f"{self.interface_id}: {self.from_range} -> {self.to_version}; "
            f"approximation: {self.approximation}; information loss: "
            f"{self.information_loss}; conservation residual: {residual}"
        )


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


@dataclass(frozen=True)
class Negotiation:
    """One resolved provider/consumer pair, and what it cost to resolve it."""

    interface_id: str
    provider_module_id: str
    consumer_module_id: str
    version: str
    adapter: MigrationAdapter | None = None

    def describe(self) -> str:
        via = f" via {self.adapter.describe()}" if self.adapter is not None else ""
        return (
            f"{self.consumer_module_id} <- {self.provider_module_id}: "
            f"{self.interface_id} v{self.version}{via}"
        )


def negotiate_configuration(
    *,
    providers: Mapping[str, tuple[str, InterfaceContract]],
    consumers: Iterable[tuple[str, RequiredInterface]],
    adapters: Iterable[MigrationAdapter] = (),
) -> tuple[Negotiation, ...]:
    """Resolve every declared provider/consumer pair, or fail before initialisation.

    VER-006, applied to a whole configuration rather than to one pair. This is
    the part that was missing: ``negotiate`` existed and nothing called it, so
    forty-four declared semantic keys in the moist column were compared against
    nothing.

    A pair that cannot be resolved directly is tried against the declared
    adapters, and a pair that no adapter covers is a failure. Refusing is the
    default and an adapter is the one way past it, which is what keeps the way
    past it from being a change to this function.
    """
    available = list(adapters)
    resolved: list[Negotiation] = []

    for consumer_id, required in consumers:
        entry = providers.get(required.interface_id)
        if entry is None:
            raise InterfaceVersionError(
                f"{consumer_id} requires interface {required.interface_id!r} and no "
                f"participating module provides it"
            )
        provider_id, provider = entry
        context = f"{consumer_id} <- {provider_id}"
        try:
            negotiate(
                provider=provider,
                consumer_interface_id=required.interface_id,
                required_range=required.version_range,
                consumer_semantics=dict(required.expects),
                required_capabilities=required.required_capabilities,
                context=context,
            )
        except InterfaceVersionError as direct:
            adapter = next((a for a in available if a.covers(provider, required)), None)
            if adapter is None:
                raise InterfaceVersionError(
                    f"{direct}\nNo migration adapter covers this gap. An adapter may "
                    f"bridge it, and must declare the approximation, information loss "
                    f"and conservation residual it introduces (VER-008)."
                ) from direct
            unbridged = {
                key
                for key, value in required.expects.items()
                if provider.physical_semantics.get(key) != value and key not in adapter.translates
            }
            if unbridged:
                raise InterfaceVersionError(
                    f"{direct}\nThe migration adapter for {required.interface_id!r} does "
                    f"not translate {sorted(unbridged)}, so it cannot bridge this gap. "
                    f"An adapter written for one semantic change must not carry another "
                    f"across unexamined (VER-007, VER-008)."
                ) from direct
            resolved.append(
                Negotiation(
                    interface_id=required.interface_id,
                    provider_module_id=provider_id,
                    consumer_module_id=consumer_id,
                    version=str(provider.version),
                    adapter=adapter,
                )
            )
            continue

        resolved.append(
            Negotiation(
                interface_id=required.interface_id,
                provider_module_id=provider_id,
                consumer_module_id=consumer_id,
                version=str(provider.version),
            )
        )
    return tuple(resolved)
