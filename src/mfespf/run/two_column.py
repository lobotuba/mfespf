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
from mfespf.physics.ocean_transport import OceanHeatTransport
from mfespf.physics.surface import SURFACE_TEMPERATURE_STATE
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


def _authorize(kernel: Kernel, provider: str, names: tuple[str, ...] = _STATES) -> Kernel:
    """Let a module between the regions contribute tendencies to the named states."""
    states = StateRegistry()
    for entry in kernel.states.current():
        if entry.record_id.split(":", 1)[-1] in names:
            entry = replace(
                entry,
                authorized_tendency_providers=entry.authorized_tendency_providers | {provider},
            )
        states.register(entry)
    return replace(kernel, states=states)


def _column(tropical_latitude: float | None, column: dict[str, Any]) -> Kernel:
    """The column each region holds: the reference column, or the tropical one."""
    if tropical_latitude is None:
        return reference_column(receives_advection=True, **column)
    from mfespf.run.tropical_column import tropical_column

    return tropical_column(latitude_limit=tropical_latitude, **column)


def two_column(
    *,
    ascent_fraction: float = OBSERVED_ASCENT_FRACTION,
    relaxation_distance: float = OBSERVED_RELAXATION_DISTANCE_M,
    run_id: str = "TWO-COLUMN",
    tropical_latitude: float | None = None,
    ocean_heat_contrast: float | None = None,
    prescribed_surface_temperatures: tuple[float, float] | None = None,
    regional_settings: dict[str, dict[str, Any]] | None = None,
    diagnose_lateral_exchange: bool = False,
    **column: Any,
) -> Kernel:
    """Two columns and the circulation between them.

    ``tropical_latitude`` makes each region the tropical column of that band
    (its sunlight, export and upwelling); ``ocean_heat_contrast`` adds the
    ocean's heat transport between their surfaces, in W m^-2, rising minus
    sinking (ocean-heat-transport-circulation.md). ``prescribed_surface_temperatures``
    holds the rising and sinking surfaces at those temperatures, in K, the energy
    that holds them declared as a boundary transfer (coupled-ascent-subsidence.md).
    ``regional_settings`` overrides column settings in one region only, by region name
    (descending-boundary-layer-ventilation.md, V1). ``diagnose_lateral_exchange``
    publishes the circulation's lateral and vertical transport per layer, output
    only (boundary-layer-lateral-exchange.md).
    """
    unknown = set(regional_settings or {}) - {_ASCENT, _DESCENT}
    if unknown:
        raise ValueError(f"regional settings for unknown regions {sorted(unknown)}")
    circulation = WeakTemperatureGradientCirculation(
        ascent_fraction=ascent_fraction,
        relaxation_distance=relaxation_distance,
        diagnose_lateral_exchange=diagnose_lateral_exchange,
    )
    ocean = (
        OceanHeatTransport(ascent_fraction=ascent_fraction, heat_contrast=ocean_heat_contrast)
        if ocean_heat_contrast is not None
        else None
    )
    regions = []
    held = (
        dict(zip((_ASCENT, _DESCENT), prescribed_surface_temperatures, strict=True))
        if (prescribed_surface_temperatures is not None)
        else {}
    )
    for region in (_ASCENT, _DESCENT):
        settings = {**column, **(regional_settings or {}).get(region, {})}
        if region in held:
            settings["prescribed_surface_temperature"] = held[region]
        kernel = _authorize(
            regionalize(_column(tropical_latitude, settings), region),
            circulation.module_id,
        )
        if ocean is not None:
            kernel = _authorize(kernel, ocean.module_id, (SURFACE_TEMPERATURE_STATE,))
        regions.append(kernel)
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
        middle_modules=[
            ModuleSchedule(circulation),
            *((ModuleSchedule(ocean),) if ocean is not None else ()),
        ],
        global_modules=[ModuleSchedule(means)],
    )
    graph = kernel.graph
    for module_id, label in (
        (circulation.module_id, "weak-temperature-gradient circulation"),
        (means.module_id, "area-weighted global means"),
    ):
        graph.add_node(Node(module_id, NodeType.MODULE, label))
    if ocean is not None:
        graph.add_node(Node(ocean.module_id, NodeType.MODULE, "ocean heat transport"))
        for region in (_ASCENT, _DESCENT):
            graph.add_edge(
                Edge(
                    source=f"{region}/surface",
                    destination=ocean.module_id,
                    exchanged=SURFACE_TEMPERATURE_STATE,
                    edge_type=EdgeType.PHYSICAL_TRANSFER,
                    interface_version="1.0.0",
                    spatial_semantics="area_weighted_regions",
                    equation_ids=("EQ-OCEAN-HEAT-TRANSPORT",),
                )
            )
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
