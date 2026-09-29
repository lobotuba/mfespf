"""A rising and a sinking column, linked by an overturning circulation.

Two ``reference_column``s filed as regions ``ascent`` and ``descent``, the
weak-temperature-gradient circulation between them, and the area-weighted
global means every global diagnostic reads. See ``physics/circulation.py`` and
``docs/protocols/two-column-circulation.md``.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

from mfespf.coupling.coupler import ModuleSchedule
from mfespf.coupling.graph import Edge, EdgeType, Node, NodeType
from mfespf.physics.circulation import (
    OBSERVED_ASCENT_FRACTION,
    OBSERVED_RELAXATION_DISTANCE_M,
    AreaWeightedMeansModule,
    WeakTemperatureGradientCirculation,
)
from mfespf.registries.states import StateRegistry
from mfespf.run.column_configurations import reference_column
from mfespf.run.driver import Kernel
from mfespf.run.regions import combine_regions, regionalize

__all__ = ["GLOBAL_MEAN_FIELDS", "two_column"]

#: Scalar fields averaged over the two regions, with their units.
GLOBAL_MEAN_FIELDS: tuple[tuple[str, str], ...] = (
    ("surface_temperature", "K"),
    ("outgoing_longwave_radiation", "W m^-2"),
    ("reflected_shortwave_radiation", "W m^-2"),
    ("toa_downward_radiation", "W m^-2"),
    ("toa_incident_irradiance", "W m^-2"),
    ("surface_latent_heat_flux", "W m^-2"),
    ("surface_sensible_heat_flux", "W m^-2"),
    ("surface_convective_heat_flux", "W m^-2"),
    ("surface_evaporation_flux", "kg m^-2 s^-1"),
    ("surface_net_radiative_flux", "W m^-2"),
    ("atmosphere_precipitable_water", "kg m^-2"),
    ("atmosphere_absorbed_shortwave", "W m^-2"),
    ("surface_absorbed_shortwave", "W m^-2"),
    ("cloud_shortwave_radiative_effect", "W m^-2"),
    ("cloud_longwave_radiative_effect", "W m^-2"),
    ("clear_sky_outgoing_longwave_radiation", "W m^-2"),
    ("clear_sky_reflected_shortwave_radiation", "W m^-2"),
)

_ASCENT, _DESCENT = "ascent", "descent"
_STATES = ("atmosphere_layer_temperature", "atmosphere_water_vapour")


def _authorize(kernel: Kernel, provider: str) -> Kernel:
    """Let the circulation contribute tendencies to the region's air and water."""
    states = StateRegistry()
    for entry in kernel.states.current():
        if entry.record_id.split(":", 1)[-1] in _STATES:
            entry = replace(
                entry,
                authorized_tendency_providers=entry.authorized_tendency_providers | {provider},
            )
        states.register(entry)
    return replace(kernel, states=states)


def two_column(
    *,
    ascent_fraction: float = OBSERVED_ASCENT_FRACTION,
    relaxation_distance: float = OBSERVED_RELAXATION_DISTANCE_M,
    run_id: str = "TWO-COLUMN",
    **column: Any,
) -> Kernel:
    """Two reference columns and the circulation between them."""
    circulation = WeakTemperatureGradientCirculation(
        ascent_fraction=ascent_fraction, relaxation_distance=relaxation_distance
    )
    regions = [
        _authorize(
            regionalize(reference_column(receives_advection=True, **column), region),
            circulation.module_id,
        )
        for region in (_ASCENT, _DESCENT)
    ]
    means = AreaWeightedMeansModule(
        regions=((_ASCENT, ascent_fraction), (_DESCENT, 1.0 - ascent_fraction)),
        fields=GLOBAL_MEAN_FIELDS,
    )
    kernel = combine_regions(
        regions,
        run_id=run_id,
        purpose=(
            "Two-column overturning circulation: a rising and a sinking reference "
            "column under the relaxed, spectral weak-temperature-gradient balance."
        ),
        middle_modules=[ModuleSchedule(circulation)],
        global_modules=[ModuleSchedule(means)],
    )
    graph = kernel.graph
    for module_id, label in (
        (circulation.module_id, "weak-temperature-gradient circulation"),
        (means.module_id, "area-weighted global means"),
    ):
        graph.add_node(Node(module_id, NodeType.MODULE, label))
    for region in (_ASCENT, _DESCENT):
        for owner, exchanged in (
            ("atmosphere", "atmosphere_layer_temperature"),
            ("water", "atmosphere_water_vapour"),
        ):
            graph.add_edge(
                Edge(
                    source=f"{region}/{owner}",
                    destination=circulation.module_id,
                    exchanged=exchanged,
                    edge_type=EdgeType.PHYSICAL_TRANSFER,
                    interface_version="1.0.0",
                    spatial_semantics="area_weighted_regions",
                    equation_ids=("EQ-WTG-OVERTURNING-CIRCULATION",),
                )
            )
    return kernel
