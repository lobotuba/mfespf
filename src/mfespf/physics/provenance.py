"""Water provenance: which water has condensed, formed precipitation, or neither.

Why this exists
---------------

Precipitation efficiency as Langhans, Yeo & Romps (2015) define it is a
probability per molecule: of the water that condensed at least once, the
fraction that reaches the surface. Gross Eulerian totals cannot give it. Water
that evaporates and condenses again is counted again, and E0's
evaporation-recondensation loop inflated both sides of every ratio
(``docs/protocols/convective-condensate-ledger.md``). No combination of gross
counters reconstructs how many distinct molecules crossed a boundary, so the
history has to travel with the water
(``docs/protocols/fall-path-closure-candidate.md``).

The classes
-----------

Water carries an irreversible history, U -> C -> R -> surface:

- U, never condensed;
- C_i, condensed at least once, first in convective layer i, never formed
  precipitation;
- R_i, formed precipitation at least once, first condensed in layer i.

Phase changes inside a class; history never goes back. Reaching the surface
ends the cohort, and the surface's evaporation starts a new one as U. A pool
is a vector of ``2 n + 1`` masses, U first, then C_0..C_(n-1), then
R_0..R_(n-1). Removing water from a pool takes each class in proportion to
its share: the declared well-mixed assumption.

The milestones count each molecule once: ``first_condensed`` by the layer of
first condensation, ``first_precip`` and ``first_surface`` by that same
origin.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "Milestones",
    "class_count",
    "condense",
    "evaporate",
    "form_precipitation",
    "proportional",
    "surface_arrival",
]


def class_count(n_layers: int) -> int:
    """U, then C_i and R_i for every layer."""
    return 2 * n_layers + 1


@dataclass
class Milestones:
    """Molecule-once counters, kg m^-2, per layer of first condensation.

    ``first_condensed`` is convective first condensation and
    ``first_condensed_large_scale`` the water module's, by the layer it
    happened in. ``first_precip`` and ``first_surface`` are by origin layer,
    and count convective-born water only: large-scale-born water reaches the
    surface in the step it condenses, so its three milestones coincide.
    ``first_precip_location`` is where convective-born water first formed
    precipitation.
    """

    first_condensed: np.ndarray
    first_condensed_large_scale: np.ndarray
    first_precip: np.ndarray
    first_precip_location: np.ndarray
    first_surface: np.ndarray

    @classmethod
    def zeros(cls, n_layers: int) -> Milestones:
        return cls(*(np.zeros(n_layers) for _ in range(5)))

    def add(self, other: Milestones) -> None:
        for name in (
            "first_condensed",
            "first_condensed_large_scale",
            "first_precip",
            "first_precip_location",
            "first_surface",
        ):
            setattr(self, name, getattr(self, name) + getattr(other, name))

    def as_array(self) -> np.ndarray:
        return np.stack(
            (
                self.first_condensed,
                self.first_condensed_large_scale,
                self.first_precip,
                self.first_precip_location,
                self.first_surface,
            )
        )


def proportional(pool: np.ndarray, amount: float) -> np.ndarray:
    """The classes of ``amount`` taken from ``pool`` in proportion to their
    shares. Exactly ``pool`` when the whole pool is taken, so a pool emptied
    is emptied to zero rather than to a rounding."""
    total = float(np.sum(pool))
    if amount <= 0.0 or total <= 0.0:
        return np.zeros_like(pool)
    if amount >= total:
        return np.array(pool, dtype=float)
    return pool * (amount / total)


def condense(
    vapour: np.ndarray, amount: float, layer: int, milestones: Milestones, *, large_scale: bool
) -> tuple[np.ndarray, np.ndarray]:
    """``amount`` of ``vapour`` condenses in ``layer``. Returns the vapour left
    and the condensate by class: U becomes C_layer, the first condensation of
    that water; C and R keep their class."""
    n = (vapour.size - 1) // 2
    taken = proportional(vapour, amount)
    condensate = np.array(taken)
    condensate[0] = 0.0
    condensate[1 + layer] += taken[0]
    if large_scale:
        milestones.first_condensed_large_scale[layer] += taken[0]
    else:
        milestones.first_condensed[layer] += taken[0]
    assert condensate.size == 2 * n + 1
    return vapour - taken, condensate


def evaporate(pool: np.ndarray, amount: float) -> tuple[np.ndarray, np.ndarray]:
    """``amount`` of condensate or rain in ``pool`` evaporates. The vapour keeps
    its class: evaporation never undoes a history. Returns the pool left and the
    vapour by class."""
    taken = proportional(pool, amount)
    return pool - taken, taken


def form_precipitation(
    cloud: np.ndarray, amount: float, layer: int, milestones: Milestones
) -> tuple[np.ndarray, np.ndarray]:
    """``amount`` of ``cloud`` water becomes precipitation in ``layer``. C_i
    becomes R_i, the first precipitation of that water, counted by its origin
    i; R_i stays R_i and is not counted again. Returns the cloud left and the
    precipitation by class."""
    n = (cloud.size - 1) // 2
    taken = proportional(cloud, amount)
    rain = np.zeros_like(taken)
    first = taken[1 : 1 + n]
    rain[1 + n :] = taken[1 + n :] + first
    milestones.first_precip += first
    milestones.first_precip_location[layer] += float(np.sum(first))
    return cloud - taken, rain


def surface_arrival(rain: np.ndarray, milestones: Milestones) -> None:
    """Precipitation reaches the surface and its cohort ends. Every molecule
    arriving is R, and arrives once, so each is counted by its origin."""
    n = (rain.size - 1) // 2
    if rain[0] != 0.0 or np.any(rain[1 : 1 + n] != 0.0):
        raise ValueError("only water that has formed precipitation reaches the surface")
    milestones.first_surface += rain[1 + n :]
