"""The transports a column composes with, and the names they publish under.

A column can be advected by more than one module at once: the tropical
column's Brewer-Dobson upwelling and band export, and a circulation linking
it to another column. Each publishes its own rate under its own name, and a
module that has to act on the profile the transports leave (the water
module's saturation adjustment, the convective adjustment) declares which
sources it composes with and sums them. Two modules once published the same
name into a column and the coupler kept whichever came later, so the
consumers never saw the circulation's transport in the step it was made
(``docs/protocols/advection-composition.md``).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

__all__ = [
    "ADVECTION_KINDS",
    "CIRCULATION_SOURCE",
    "advection_field",
    "composed_advection",
]

#: What each kind of advected rate is called, and in what units, for the
#: source that publishes under the plain name.
ADVECTION_KINDS = {
    "vapour": ("vapour_advection_rate", "kg m^-2 s^-1"),
    "heating": ("advective_heating_rate", "W m^-2"),
}

#: The circulation's source name: it publishes ``circulation_<plain name>``.
CIRCULATION_SOURCE = "circulation"


def advection_field(source: str, kind: str) -> str:
    """The field a source publishes one kind of rate under.

    The empty source is the plain name, which the tropical exchange uses.
    """
    plain = ADVECTION_KINDS[kind][0]
    return f"{source}_{plain}" if source else plain


def composed_advection(
    inputs: Mapping[str, Any], sources: Sequence[str], kind: str
) -> tuple[np.ndarray | None, Any]:
    """The sum of every declared source's rate that is present, and the first
    field (for its time interval), or ``(None, None)`` if none is.

    One source is returned as that source's array, not as a sum that could
    round differently, so a column advected by one module composes exactly as
    it did before there could be two.
    """
    units = ADVECTION_KINDS[kind][1]
    present = [
        inputs[name]
        for name in (advection_field(s, kind) for s in sources)
        if inputs.get(name) is not None
    ]
    if not present:
        return None, None
    total = np.asarray(present[0].require_units(units).values, dtype=float)
    for field in present[1:]:
        total = total + np.asarray(field.require_units(units).values, dtype=float)
    return total, present[0]
