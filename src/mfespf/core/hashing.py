"""Canonical serialisation and hashing.

AUD-002 requires a configuration hash, INF-006 requires source, configuration,
parameter and claim hashes recorded *before* evaluation data are introduced, and
A4.5 requires a freeze record.  All of those depend on one thing: a canonical
form whose hash does not change when a dict is rebuilt in a different order or a
float is round-tripped through YAML.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, is_dataclass
from enum import Enum
from pathlib import Path
from typing import Any

import numpy as np

__all__ = ["canonical", "canonical_json", "digest", "file_digest"]


def canonical(obj: Any) -> Any:
    """Reduce an object to JSON-safe primitives with deterministic ordering.

    Floats are rendered at 17 significant digits (round-trip exact for float64)
    so that a hash is stable across serialisation formats.  Non-finite values are
    encoded as strings rather than dropped: a NaN in a configuration is a fact the
    audit trail should record, not hide.
    """
    if obj is None or isinstance(obj, (bool, int, str)):
        return obj
    if isinstance(obj, float):
        if math.isnan(obj):
            return "NaN"
        if math.isinf(obj):
            return "Infinity" if obj > 0 else "-Infinity"
        return f"{obj:.17g}"
    if isinstance(obj, Enum):
        return canonical(obj.value)
    if isinstance(obj, Path):
        return obj.as_posix()
    if isinstance(obj, np.ndarray):
        return [canonical(float(v)) for v in obj.reshape(-1)]
    if isinstance(obj, (np.floating, np.integer)):
        return canonical(obj.item())
    if is_dataclass(obj) and not isinstance(obj, type):
        return canonical(asdict(obj))
    if isinstance(obj, dict):
        return {str(k): canonical(obj[k]) for k in sorted(obj, key=str)}
    if isinstance(obj, (list, tuple)):
        return [canonical(v) for v in obj]
    if isinstance(obj, (set, frozenset)):
        return [canonical(v) for v in sorted(obj, key=str)]
    return str(obj)


def canonical_json(obj: Any) -> str:
    return json.dumps(canonical(obj), sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def digest(obj: Any, *, length: int = 64) -> str:
    """SHA-256 of the canonical form, truncated to ``length`` hex characters."""
    return hashlib.sha256(canonical_json(obj).encode("utf-8")).hexdigest()[:length]


def file_digest(path: str | Path, *, length: int = 64) -> str:
    """SHA-256 of a file's bytes, for dataset provenance (A4.3, DAT-009)."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:length]
