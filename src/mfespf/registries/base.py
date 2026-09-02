"""Generic versioned registry.

KER-001 requires registries for modules, states, equations, datasets, claims,
experiments and runs.  They share one behaviour, implemented once here: records
are immutable, identified by ``(record_id, version)``, and a superseding record
never overwrites the one it supersedes (AUD-006, EQR-007, A8.4).
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, replace
from typing import Generic, Protocol, TypeVar

from mfespf.core.errors import RegistryError
from mfespf.core.hashing import digest
from mfespf.core.versioning import Version

__all__ = ["Record", "Registry", "RegistryEntry"]


class Record(Protocol):
    """Anything a registry can hold: an id and a semantic version.

    Declared read-only so that frozen dataclasses -- which every registry record
    is, because records are immutable -- satisfy the protocol.
    """

    @property
    def record_id(self) -> str: ...

    @property
    def version(self) -> Version: ...


R = TypeVar("R", bound=Record)


@dataclass(frozen=True, slots=True)
class RegistryEntry(Generic[R]):
    record: R
    content_hash: str
    superseded_by: tuple[str, str] | None = None

    @property
    def key(self) -> tuple[str, str]:
        return (self.record.record_id, str(self.record.version))

    @property
    def is_current(self) -> bool:
        return self.superseded_by is None


class Registry(Generic[R]):
    """An append-only store of immutable versioned records."""

    def __init__(self, kind: str) -> None:
        self.kind = kind
        self._entries: dict[tuple[str, str], RegistryEntry[R]] = {}

    def register(self, record: R) -> RegistryEntry[R]:
        key = (record.record_id, str(record.version))
        content_hash = digest(record, length=16)
        existing = self._entries.get(key)
        if existing is not None:
            if existing.content_hash == content_hash:
                return existing
            raise RegistryError(
                f"{self.kind} {key[0]!r} version {key[1]} is already registered with different "
                f"content (hash {existing.content_hash} vs {content_hash}). Registering a changed "
                f"record under an existing version would break run lineage; bump the version."
            )
        entry = RegistryEntry(record=record, content_hash=content_hash)
        self._entries[key] = entry
        return entry

    def supersede(
        self, old_id: str, old_version: Version | str, replacement: R
    ) -> RegistryEntry[R]:
        """Register ``replacement`` and mark the older record as superseded.

        The old record stays retrievable (EQR-007, AUD-006): a published result
        must remain reconstructible from the equation versions it actually used.
        """
        old_key = (old_id, str(old_version))
        if old_key not in self._entries:
            raise RegistryError(f"cannot supersede unknown {self.kind} {old_key}")
        new_entry = self.register(replacement)
        old = self._entries[old_key]
        self._entries[old_key] = replace(old, superseded_by=new_entry.key)
        return new_entry

    def get(self, record_id: str, version: Version | str | None = None) -> R:
        if version is not None:
            key = (record_id, str(version))
            try:
                return self._entries[key].record
            except KeyError:
                raise RegistryError(f"unknown {self.kind} {key}") from None
        versions = self.versions_of(record_id)
        if not versions:
            raise RegistryError(f"unknown {self.kind} {record_id!r}")
        return self._entries[(record_id, str(max(versions)))].record

    def entry(self, record_id: str, version: Version | str) -> RegistryEntry[R]:
        try:
            return self._entries[(record_id, str(version))]
        except KeyError:
            raise RegistryError(f"unknown {self.kind} {(record_id, str(version))}") from None

    def versions_of(self, record_id: str) -> list[Version]:
        return sorted(Version.parse(v) for (rid, v) in self._entries if rid == record_id)

    def ids(self) -> list[str]:
        return sorted({rid for rid, _ in self._entries})

    def current(self) -> list[R]:
        """Latest non-superseded record for each id."""
        out = []
        for rid in self.ids():
            versions = [v for v in self.versions_of(rid) if self._entries[(rid, str(v))].is_current]
            if versions:
                out.append(self._entries[(rid, str(max(versions)))].record)
        return out

    def bindings(self) -> dict[str, str]:
        """``{record_id: version}`` for every registered record, for a manifest."""
        return {
            f"{rid}@{ver}": self._entries[(rid, ver)].content_hash
            for rid, ver in sorted(self._entries)
        }

    def __contains__(self, record_id: object) -> bool:
        return any(rid == record_id for rid, _ in self._entries)

    def __iter__(self) -> Iterator[R]:
        return (self._entries[k].record for k in sorted(self._entries))

    def __len__(self) -> int:
        return len(self._entries)

    def __repr__(self) -> str:  # pragma: no cover - convenience
        return f"<Registry {self.kind}: {len(self)} records, {len(self.ids())} ids>"
