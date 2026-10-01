"""A column standing for the tropics: stage 1 of the stratospheric-water build.

The reference column over the band between the tropical upwelling's turnaround
latitudes, with the band's own sunlight, its heat export and the Brewer-Dobson
upwelling through its tropopause (``docs/protocols/tropical-column.md``). What
it has to show is the cold tropical tropopause that freeze-dries the air
entering the stratosphere, which a global-mean column cannot have.
"""

from __future__ import annotations

import math
from dataclasses import replace
from typing import Any

import numpy as np

from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError
from mfespf.coupling.coupler import ModuleSchedule
from mfespf.coupling.graph import Edge, EdgeType, Node, NodeType
from mfespf.physics.orbit import band_insolation_factor
from mfespf.physics.tropical_exchange import (
    LATERAL_ENERGY_STATE,
    LATERAL_WATER_STATE,
    OBSERVED_UPWELLING_MASS_FLUX_KG_S,
    LateralBoundaryModule,
    TropicalExchangeModule,
    lateral_boundary_state_entries,
)
from mfespf.registries.states import StateRegistry
from mfespf.run.column_configurations import reference_column
from mfespf.run.driver import Kernel

__all__ = ["OBSERVED_BAND_HEAT_EXPORT_W_M2", "tropical_column", "upwelling_per_unit_area"]

#: A band's heat export, W m^-2 of the band: NCEP/NCAR R1 1991-2020 band-mean net
#: top-of-atmosphere radiation minus its global mean (-12.87), which removes a
#: bias taken as uniform. Keyed by the band's latitude limit in degrees.
OBSERVED_BAND_HEAT_EXPORT_W_M2: dict[float, float] = {20.0: 42.51, 25.0: 38.464, 30.0: 32.95}

_TEMPERATURE = "atmosphere_layer_temperature"
_VAPOUR = "atmosphere_water_vapour"


def upwelling_per_unit_area(latitude_limit_deg: float) -> float:
    """The observed upwelling mass flux spread over the band, kg m^-2 s^-1."""
    radius = value_of("earth_radius_mean")
    band_area = 4.0 * math.pi * radius**2 * math.sin(math.radians(latitude_limit_deg))
    return OBSERVED_UPWELLING_MASS_FLUX_KG_S / band_area


def tropical_column(
    *,
    latitude_limit: float = 25.0,
    upwelling: bool = True,
    run_id: str | None = None,
    **column: Any,
) -> Kernel:
    """The reference column made to stand for the band ``|lat| <= latitude_limit``."""
    if latitude_limit not in OBSERVED_BAND_HEAT_EXPORT_W_M2:
        raise ConfigurationError(
            f"no measured heat export for a {latitude_limit} degree band; measured for "
            f"{sorted(OBSERVED_BAND_HEAT_EXPORT_W_M2)}"
        )
    column.setdefault("run_id", run_id or f"TROPICAL-COLUMN-{latitude_limit:g}")
    column.setdefault(
        "purpose",
        "A column standing for the tropics: the band's sunlight, its heat export and "
        "the Brewer-Dobson upwelling through its tropopause.",
    )
    kernel = reference_column(
        receives_advection=True,
        geometry_factor=band_insolation_factor(latitude_limit),
        **column,
    )
    exchange = TropicalExchangeModule(
        upwelling_flux=upwelling_per_unit_area(latitude_limit) if upwelling else 0.0,
        heat_export=OBSERVED_BAND_HEAT_EXPORT_W_M2[latitude_limit],
    )
    boundary = LateralBoundaryModule()

    # After the radiation, so this step's profile is the one exchanged, and
    # before the water and convection, which compose with it; the boundary last,
    # reading what the exchange published.
    # In the column's splitting order, not the order its schedules were
    # appended in: the two differ when a configuration reorders modules
    # (condense_after_convection runs convection before the water).
    schedules = list(kernel.schedule.ordered())
    ids = [s.module_id for s in schedules]
    if "radiation" not in ids:
        raise ConfigurationError("the tropical exchange needs a radiation module to follow")
    cut = ids.index("radiation") + 1
    schedules = [
        *schedules[:cut],
        ModuleSchedule(exchange),
        *schedules[cut:],
        ModuleSchedule(boundary),
    ]

    states = StateRegistry()
    for entry in kernel.states.current():
        if entry.record_id in (_TEMPERATURE, _VAPOUR):
            entry = replace(
                entry,
                authorized_tendency_providers=entry.authorized_tendency_providers
                | {exchange.module_id},
            )
        states.register(entry)
    for entry in lateral_boundary_state_entries(owner=boundary.module_id):
        states.register(entry)

    graph = kernel.graph
    graph.add_node(Node(exchange.module_id, NodeType.MODULE, "tropical lateral exchange"))
    graph.add_node(Node(boundary.module_id, NodeType.BOUNDARY, "the rest of the planet"))
    for source, destination, exchanged in (
        ("atmosphere", exchange.module_id, _TEMPERATURE),
        ("water", exchange.module_id, _VAPOUR),
        (exchange.module_id, boundary.module_id, "lateral_energy_export"),
        (exchange.module_id, boundary.module_id, "lateral_water_export"),
    ):
        graph.add_edge(
            Edge(
                source=source,
                destination=destination,
                exchanged=exchanged,
                edge_type=EdgeType.PHYSICAL_TRANSFER,
                interface_version="1.0.0",
                spatial_semantics="column_to_lateral_boundary",
                equation_ids=("EQ-BREWER-DOBSON-UPWELLING", "EQ-LATERAL-HEAT-EXPORT"),
            )
        )

    return replace(
        kernel,
        schedule=replace(
            kernel.schedule,
            schedules=schedules,
            splitting_order=tuple(s.module_id for s in schedules),
        ),
        states=states,
        graph=graph,
        initial_state={
            **kernel.initial_state,
            LATERAL_ENERGY_STATE: np.asarray(0.0),
            LATERAL_WATER_STATE: np.asarray(0.0),
        },
    )
