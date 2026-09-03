"""Stage 2 reference configurations: the radiating atmospheric column.

Charter section 19 stage 2 and requirement KER-012, which asks for a transparent
atmosphere and an absorbing-and-emitting one alongside the airless body.

All three share one builder and differ only in the absorption coefficient, which
is the point: the transparent case is the same code with ``kappa = 0``, so when
it reproduces the Stage 1 airless result it is evidence about the column
machinery rather than about a separately written special case.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from mfespf.contracts.cost import ResourceBudget, ResourceClass
from mfespf.contracts.module import ModuleStatus
from mfespf.contracts.spectral import SpectralGrid
from mfespf.core.constants import value_of
from mfespf.coupling.coupler import CouplingSchedule, ModuleSchedule
from mfespf.coupling.graph import Cycle, DependencyGraph, Edge, EdgeType, Node, NodeType
from mfespf.physics.atmosphere import (
    ATMOSPHERE_TEMPERATURE_STATE,
    AtmosphericColumnModule,
    atmosphere_state_entry,
)
from mfespf.physics.boundary import SPACE_ENERGY_STATE, SpaceBoundaryModule, space_state_entry
from mfespf.physics.clouds import CloudProfile, PrescribedCloudModule, single_layer_cloud
from mfespf.physics.column_surface import ColumnSurfaceModule
from mfespf.physics.convection import ConvectiveAdjustmentModule
from mfespf.physics.equations import EQUATION_IDS, build_equation_registry
from mfespf.physics.ocean import (
    OCEAN_DEEP_STATE,
    DeepOceanModule,
    ocean_deep_state_entry,
)
from mfespf.physics.orbit import OrbitalGeometryModule
from mfespf.physics.radiative_transfer import (
    SpectralRadiationModule,
    grey_absorption,
    window_absorption,
)
from mfespf.physics.solar import ConstantSpectrumSolarSource, banded_grid
from mfespf.physics.surface import (
    SURFACE_TEMPERATURE_STATE,
    SURFACE_WATER_STATE,
    surface_state_entry,
    surface_water_state_entry,
)
from mfespf.physics.water import VAPOUR_STATE, WaterModule, vapour_state_entry
from mfespf.registries.states import StateRegistry
from mfespf.run.driver import Kernel

__all__ = [
    "absorbing_emitting_atmosphere",
    "column_grid",
    "column_over_reduced_ocean",
    "declare_explicit_cycles",
    "moist_radiative_convective",
    "no_cloud_control",
    "prescribed_cloud",
    "radiating_column",
    "radiative_convective",
    "shortwave_albedo",
    "single_grey_layer",
    "transparent_atmosphere",
    "window_atmosphere",
]

# A 50 m ocean mixed layer, as in the airless configuration.
_DEFAULT_AREAL_HEAT_CAPACITY = 1025.0 * 3985.0 * 50.0

_SPLITTING_ORDER = ("solar", "orbit", "radiation", "atmosphere", "surface", "space")


def column_grid(n_bands: int = 32) -> SpectralGrid:
    """The default column grid: 0.1 micron to 1 mm, logarithmically spaced.

    The upper edge is well beyond the terrestrial peak on purpose. At 255 K a
    grid stopping at 100 microns omits about 0.7 percent of the Planck integral,
    which is large compared with the tolerances the analytic gates use; at 1 mm
    the omission is order 1e-5 and the redistribution ``band_exitance`` performs
    is negligible rather than merely small.
    """
    return banded_grid(n_bands, 1.0e-7, 1.0e-3)


def shortwave_albedo(
    grid: SpectralGrid, albedo: float = 0.3, cutoff_m: float = 4.0e-6
) -> np.ndarray:
    """Reflective to sunlight, black in the thermal infrared.

    This is the idealisation the airless Stage 1 configuration assumed without
    saying so: it reflected 30 percent of the incoming solar flux while emitting
    as a black body. Under Kirchhoff's law those two statements are only
    compatible if the reflectance is confined to short wavelengths -- a surface
    that is genuinely grey at 0.3 has an infrared emissivity of 0.7, and reaches
    a *different* equilibrium.

    The column model makes emissivity follow from albedo band by band, so the
    inconsistency becomes visible rather than staying implicit. Passing a plain
    float to a column configuration gives the genuinely grey surface; this helper
    gives the shortwave-only one, which is what reproduces the Stage 1 result.
    """
    if not 0.0 <= albedo <= 1.0:
        raise ValueError(f"albedo {albedo} must lie in [0, 1]")
    return np.where(grid.centers < cutoff_m, albedo, 0.0)


def declare_explicit_cycles(
    graph: DependencyGraph,
    *,
    splitting_order: tuple[str, ...],
    scheme: str,
    criterion: str,
) -> list[Cycle]:
    """Declare every cycle the graph contains under one resolution scheme.

    CAU-005 requires each feedback loop to name its coupling scheme, iteration
    order, convergence condition and residual. A radiating column contains many
    loops -- temperature feeds emission feeds heating feeds temperature, through
    both the surface and every layer -- and the coupler resolves all of them the
    same way: one explicit pass per master step in the declared splitting order.

    Declaring them uniformly is therefore an accurate statement rather than a
    rubber stamp, and the cycles are enumerated into the run manifest so a reader
    can see exactly which loops that claim covers. A configuration that resolved
    different loops differently would have to declare them separately, and this
    helper would be the wrong tool for it.
    """
    declared = []
    for nodes in graph.find_cycles():
        order = tuple(name for name in splitting_order if name in nodes)
        declared.append(
            graph.declare_cycle(
                Cycle(
                    nodes=nodes,
                    coupling_scheme=scheme,
                    iteration_order=order or nodes,
                    convergence_criterion=criterion,
                )
            )
        )
    return declared


def _build_state_registry(
    *,
    n_layers: int,
    convection: bool = False,
    water: bool = False,
    deep_ocean: bool = False,
) -> StateRegistry:
    """Register the states this configuration actually contains.

    The authorisation lists grow with the configuration rather than being fixed:
    a module that is not present must not appear as an authorised provider, or
    STA-003 would be satisfied by a declaration nothing backs.
    """
    providers = {"radiation"}
    if convection:
        providers.add("convection")
    if water:
        providers.add("water")

    readers = {"radiation", "surface", *providers}
    states = StateRegistry()
    states.register(
        surface_state_entry(
            owner="surface",
            readers=frozenset(readers - {"surface"}),
            tendency_providers=frozenset({"convection"}) if convection else frozenset(),
        )
    )
    states.register(
        atmosphere_state_entry(
            owner="atmosphere",
            n_layers=n_layers,
            tendency_providers=frozenset(providers),
            readers=frozenset(readers),
        )
    )
    states.register(space_state_entry(owner="space"))
    if water:
        states.register(
            vapour_state_entry(
                owner="water", n_layers=n_layers, readers=frozenset({"radiation", "atmosphere"})
            )
        )
        states.register(surface_water_state_entry(owner="surface"))
    if deep_ocean:
        states.register(ocean_deep_state_entry(owner="ocean", readers=frozenset({"surface"})))
    return states


def _build_graph(
    *,
    configuration_id: str,
    coupling_frequency: str,
    convection: bool = False,
    water: bool = False,
    clouds: bool = False,
    deep_ocean: bool = False,
    splitting_order: tuple[str, ...] = _SPLITTING_ORDER,
) -> DependencyGraph:
    graph = DependencyGraph(configuration_id)
    modules = [
        ("solar", "solar source"),
        ("orbit", "orbital geometry and TOA boundary"),
        ("radiation", "spectral column radiative transfer"),
        ("atmosphere", "hydrostatic radiative column"),
        ("surface", "slab surface"),
        ("space", "system boundary sink"),
    ]
    if convection:
        modules.append(("convection", "dry convective adjustment"))
    if water:
        modules.append(("water", "atmospheric water vapour"))
    if clouds:
        modules.append(("cloud", "prescribed cloud state"))
    if deep_ocean:
        modules.append(("ocean", "reduced deep-ocean reservoir"))
    for module_id, label in modules:
        graph.add_node(Node(module_id, NodeType.MODULE, label))
    graph.add_node(Node(SURFACE_TEMPERATURE_STATE, NodeType.STATE, "slab temperature"))
    graph.add_node(Node(ATMOSPHERE_TEMPERATURE_STATE, NodeType.STATE, "layer temperatures"))
    graph.add_node(Node(SPACE_ENERGY_STATE, NodeType.STATE, "cumulative outgoing energy"))

    def edge(source: str, destination: str, exchanged: str, **kwargs: Any) -> None:
        defaults: dict[str, Any] = {
            "edge_type": EdgeType.BOUNDARY_DEPENDENCY,
            "interface_version": "1.0.0",
            "coupling_frequency": coupling_frequency,
            "spatial_semantics": "global_scalar",
        }
        defaults.update(kwargs)
        graph.add_edge(
            Edge(source=source, destination=destination, exchanged=exchanged, **defaults)
        )

    edge("solar", "orbit", "solar_total_irradiance_1au", units="W m^-2")
    edge(
        "solar",
        "orbit",
        "solar_spectral_irradiance_1au",
        units="W m^-2 m^-1",
        spectral_semantics="carried_with_its_own_grid",
    )
    edge(
        "orbit",
        "radiation",
        "toa_spectral_irradiance",
        units="W m^-2 m^-1",
        spectral_semantics="carried_with_its_own_grid",
        equation_ids=("EQ-SPHERICAL-MEAN-INSOLATION", "EQ-INVERSE-SQUARE-IRRADIANCE"),
    )
    edge(
        "orbit",
        "atmosphere",
        "toa_incident_energy",
        edge_type=EdgeType.PHYSICAL_TRANSFER,
        units="J m^-2",
        time_semantics="interval_accumulation",
        equation_ids=("EQ-SPHERICAL-MEAN-INSOLATION",),
    )

    # Radiation reads the state it needs; the owners publish it before the step.
    edge(
        ATMOSPHERE_TEMPERATURE_STATE,
        "radiation",
        ATMOSPHERE_TEMPERATURE_STATE,
        edge_type=EdgeType.STATE_DEPENDENCY,
        units="K",
        equation_ids=("EQ-SCHWARZSCHILD-TRANSFER", "EQ-PLANCK-SPECTRAL-RADIANCE"),
        uncertainty_rule="layer temperature uncertainty propagates into emission",
    )
    edge(
        SURFACE_TEMPERATURE_STATE,
        "radiation",
        SURFACE_TEMPERATURE_STATE,
        edge_type=EdgeType.STATE_DEPENDENCY,
        units="K",
        equation_ids=("EQ-PLANCK-SPECTRAL-RADIANCE",),
    )
    edge(
        "surface",
        "radiation",
        "surface_spectral_albedo",
        edge_type=EdgeType.PARAMETER_DEPENDENCY,
        units="1",
        spectral_semantics="carried_with_its_own_grid",
    )
    edge(
        "atmosphere",
        "radiation",
        "atmosphere_layer_mass",
        edge_type=EdgeType.PARAMETER_DEPENDENCY,
        units="kg m^-2",
        equation_ids=("EQ-HYDROSTATIC-LAYER-MASS", "EQ-BEER-LAMBERT-EXTINCTION"),
    )

    # Radiation determines the tendencies the owners commit.
    edge(
        "radiation",
        "atmosphere",
        "layer_radiative_heating",
        edge_type=EdgeType.STATE_DEPENDENCY,
        units="W m^-2",
        equation_ids=("EQ-SCHWARZSCHILD-TRANSFER", "EQ-LAYER-ENERGY-BALANCE"),
    )
    edge(
        "radiation",
        "surface",
        "surface_downward_radiation",
        edge_type=EdgeType.STATE_DEPENDENCY,
        units="W m^-2",
        equation_ids=("EQ-SCHWARZSCHILD-TRANSFER", "EQ-SLAB-ENERGY-BALANCE"),
    )

    # The conserved-energy transfers between the stores.
    for source, destination, exchanged, equation in (
        ("atmosphere", "surface", "radiation_down_surface", "EQ-SCHWARZSCHILD-TRANSFER"),
        ("surface", "atmosphere", "radiation_up_surface", "EQ-PLANCK-SPECTRAL-RADIANCE"),
        ("atmosphere", "space", "radiation_up_toa", "EQ-SCHWARZSCHILD-TRANSFER"),
    ):
        edge(
            source,
            destination,
            exchanged,
            edge_type=EdgeType.PHYSICAL_TRANSFER,
            units="J m^-2",
            time_semantics="interval_accumulation",
            equation_ids=(equation,),
        )

    for owner, state, equation in (
        ("atmosphere", ATMOSPHERE_TEMPERATURE_STATE, "EQ-LAYER-ENERGY-BALANCE"),
        ("surface", SURFACE_TEMPERATURE_STATE, "EQ-SLAB-ENERGY-BALANCE"),
        ("space", SPACE_ENERGY_STATE, "EQ-BOUNDARY-ACCUMULATION"),
    ):
        edge(
            owner,
            state,
            f"{state}_update",
            edge_type=EdgeType.NUMERICAL_DEPENDENCY,
            units="K" if state != SPACE_ENERGY_STATE else "J m^-2",
            equation_ids=(equation, "EQ-NUM-FORWARD-EULER"),
        )

    if convection:
        edge(
            "radiation",
            "convection",
            "layer_radiative_heating",
            edge_type=EdgeType.DIAGNOSTIC_DEPENDENCY,
            units="W m^-2",
        )
        edge(
            ATMOSPHERE_TEMPERATURE_STATE,
            "convection",
            ATMOSPHERE_TEMPERATURE_STATE,
            edge_type=EdgeType.STATE_DEPENDENCY,
            units="K",
            equation_ids=("EQ-POTENTIAL-TEMPERATURE",),
        )
        edge(
            "convection",
            "atmosphere",
            "convective_temperature_tendency",
            edge_type=EdgeType.STATE_DEPENDENCY,
            units="K s^-1",
            equation_ids=("EQ-DRY-CONVECTIVE-ADJUSTMENT",),
            uncertainty_rule="structural; the closure has no parametric uncertainty",
        )

    if water:
        graph.add_node(Node(VAPOUR_STATE, NodeType.STATE, "layer water vapour"))
        graph.add_node(Node(SURFACE_WATER_STATE, NodeType.STATE, "surface water store"))
        edge(
            SURFACE_TEMPERATURE_STATE,
            "water",
            SURFACE_TEMPERATURE_STATE,
            edge_type=EdgeType.STATE_DEPENDENCY,
            units="K",
            equation_ids=("EQ-CLAUSIUS-CLAPEYRON-SATURATION",),
        )
        edge(
            ATMOSPHERE_TEMPERATURE_STATE,
            "water",
            ATMOSPHERE_TEMPERATURE_STATE,
            edge_type=EdgeType.STATE_DEPENDENCY,
            units="K",
            equation_ids=("EQ-SATURATION-SPECIFIC-HUMIDITY",),
        )
        edge(
            "surface",
            "water",
            "evaporation_water",
            edge_type=EdgeType.PHYSICAL_TRANSFER,
            units="kg m^-2",
            time_semantics="interval_accumulation",
            equation_ids=("EQ-BULK-SURFACE-EVAPORATION",),
        )
        edge(
            "surface",
            "water",
            "evaporation_latent_energy",
            edge_type=EdgeType.PHYSICAL_TRANSFER,
            units="J m^-2",
            time_semantics="interval_accumulation",
            equation_ids=("EQ-LATENT-ENERGY-OF-VAPOUR",),
        )
        edge(
            "water",
            "surface",
            "precipitation_water",
            edge_type=EdgeType.PHYSICAL_TRANSFER,
            units="kg m^-2",
            time_semantics="interval_accumulation",
            equation_ids=("EQ-SATURATION-ADJUSTMENT-CONDENSATION",),
        )
        edge(
            "water",
            "atmosphere",
            "latent_heat_release",
            edge_type=EdgeType.PHYSICAL_TRANSFER,
            units="J m^-2",
            time_semantics="interval_accumulation",
            equation_ids=("EQ-SATURATION-ADJUSTMENT-CONDENSATION",),
        )
        edge(
            "water",
            VAPOUR_STATE,
            f"{VAPOUR_STATE}_update",
            edge_type=EdgeType.NUMERICAL_DEPENDENCY,
            units="kg m^-2",
            equation_ids=("EQ-NUM-FORWARD-EULER",),
        )
        edge(
            "surface",
            SURFACE_WATER_STATE,
            f"{SURFACE_WATER_STATE}_update",
            edge_type=EdgeType.NUMERICAL_DEPENDENCY,
            units="kg m^-2",
            equation_ids=("EQ-NUM-FORWARD-EULER",),
        )

    if clouds:
        for exchanged, units in (
            ("cloud_fraction", "1"),
            ("cloud_liquid_water_path", "kg m^-2"),
            ("cloud_ice_water_path", "kg m^-2"),
            ("cloud_effective_radius", "m"),
        ):
            edge(
                "cloud",
                "radiation",
                exchanged,
                edge_type=EdgeType.BOUNDARY_DEPENDENCY,
                units=units,
                equation_ids=("EQ-CLOUD-OPTICAL-PROPERTIES",),
                uncertainty_rule=(
                    "prescribed; the cloud responds to nothing, so no feedback path exists"
                ),
            )

    if deep_ocean:
        graph.add_node(Node(OCEAN_DEEP_STATE, NodeType.STATE, "deep reservoir temperature"))
        edge(
            SURFACE_TEMPERATURE_STATE,
            "ocean",
            SURFACE_TEMPERATURE_STATE,
            edge_type=EdgeType.STATE_DEPENDENCY,
            units="K",
            equation_ids=(EQUATION_IDS["deep_ocean_exchange"],),
        )
        edge(
            "surface",
            "ocean",
            "surface_heat_capacity",
            edge_type=EdgeType.PARAMETER_DEPENDENCY,
            units="J m^-2 K^-1",
        )
        edge(
            "ocean",
            "surface",
            "deep_ocean_heat_flux",
            edge_type=EdgeType.STATE_DEPENDENCY,
            units="W m^-2",
            equation_ids=(EQUATION_IDS["deep_ocean_exchange"],),
        )
        edge(
            "surface",
            "ocean",
            "deep_ocean_heat",
            edge_type=EdgeType.PHYSICAL_TRANSFER,
            units="J m^-2",
            time_semantics="interval_accumulation",
            equation_ids=(EQUATION_IDS["deep_ocean_exchange"],),
        )
        edge(
            "ocean",
            OCEAN_DEEP_STATE,
            f"{OCEAN_DEEP_STATE}_update",
            edge_type=EdgeType.NUMERICAL_DEPENDENCY,
            units="K",
            equation_ids=("EQ-NUM-FORWARD-EULER",),
        )

    declare_explicit_cycles(
        graph,
        splitting_order=splitting_order,
        scheme="explicit_single_pass_per_master_step_in_the_declared_splitting_order",
        criterion=(
            "no iteration; the residual of the lagged coupling is bounded by the "
            "time-step convergence study rather than by an inner convergence test"
        ),
    )
    return graph


def radiating_column(
    *,
    run_id: str,
    purpose: str,
    absorption_coefficient: np.ndarray,
    grid: SpectralGrid | None = None,
    n_layers: int = 8,
    albedo: float | np.ndarray | None = None,
    areal_heat_capacity: float = _DEFAULT_AREAL_HEAT_CAPACITY,
    initial_surface_temperature: float = 250.0,
    initial_air_temperature: float | np.ndarray = 250.0,
    total_irradiance: float | None = None,
    eccentricity: float = 0.0,
    surface_pressure: float | None = None,
    top_pressure: float = 1.0e3,
    master_dt: float = 3600.0,
    solar_step_multiple: int = 1,
    diffusivity_factor: float | None = None,
    convection: bool = False,
    water: bool = False,
    cloud: CloudProfile | None = None,
    deep_ocean: bool = False,
    deep_ocean_depth: float = 1200.0,
    deep_ocean_exchange: float | None = None,
    initial_deep_temperature: float | None = None,
    evaporation_coefficient: float = 2.0e-4,
    initial_water_store: float = 1.0e4,
    conservation_tolerance: float = 1e-11,
    budget: ResourceBudget | None = None,
) -> Kernel:
    """Assemble a Stage 2 column configuration.

    The absorption coefficient is the only thing that distinguishes the
    transparent, grey and spectrally structured cases from one another.
    """
    working_grid = grid if grid is not None else column_grid()
    kappa = np.asarray(absorption_coefficient, dtype=float)
    # Default to a surface that reflects sunlight and is black in the infrared.
    # A plain float would be grey in every band, which is a different -- and
    # physically distinguishable -- surface; see ``shortwave_albedo``.
    surface_reflectance = shortwave_albedo(working_grid) if albedo is None else albedo

    solar = ConstantSpectrumSolarSource(
        grid=working_grid, total_irradiance=total_irradiance, status=ModuleStatus.PRESCRIBED
    )
    orbit = OrbitalGeometryModule(eccentricity=eccentricity, receiver_module_id="atmosphere")
    radiation = SpectralRadiationModule(
        grid=working_grid,
        absorption_coefficient=kappa,
        diffusivity_factor=diffusivity_factor,
        couple_clouds=cloud is not None,
    )
    atmosphere = AtmosphericColumnModule(
        n_layers=n_layers,
        surface_pressure=surface_pressure,
        top_pressure=top_pressure,
        initial_temperature=initial_air_temperature,
        couple_water=water,
    )
    surface = ColumnSurfaceModule(
        grid=working_grid,
        areal_heat_capacity=areal_heat_capacity,
        albedo=surface_reflectance,
        initial_temperature=initial_surface_temperature,
        couple_water=water,
        initial_water_store=initial_water_store,
        couple_deep_ocean=deep_ocean,
    )
    space = SpaceBoundaryModule(
        sender_module_id="atmosphere",
        incoming_flux_fields=(("toa_upward_radiation", "radiation_up_toa"),),
    )

    schedules = [
        ModuleSchedule(solar, step_multiple=solar_step_multiple),
        ModuleSchedule(orbit, step_multiple=solar_step_multiple),
    ]
    order_head = ["solar", "orbit"]
    if cloud is not None:
        # The cloud is published before radiation reads it: a prescribed state is
        # a boundary condition on the radiation, not a response to it.
        schedules.append(
            ModuleSchedule(
                PrescribedCloudModule(grid=working_grid, n_layers=n_layers, profile=cloud),
                step_multiple=1,
            )
        )
        order_head.append("cloud")
    schedules.append(ModuleSchedule(radiation, step_multiple=1))
    order_head.append("radiation")
    # Splitting order is physics, not bookkeeping: water releases latent heat
    # before convection so that convection stabilises the profile condensation
    # has already produced, rather than one it is about to destabilise.
    order = order_head
    if water:
        schedules.append(
            ModuleSchedule(
                WaterModule(n_layers=n_layers, evaporation_coefficient=evaporation_coefficient),
                step_multiple=1,
            )
        )
        order.append("water")
    if convection:
        schedules.append(ModuleSchedule(ConvectiveAdjustmentModule(), step_multiple=1))
        order.append("convection")
    if deep_ocean:
        # Ahead of the surface, like every other module the surface reads from:
        # the exchange is computed from the temperature committed at the end of
        # the previous step and applied within this one.
        schedules.append(
            ModuleSchedule(
                DeepOceanModule(
                    depth=deep_ocean_depth,
                    exchange_coefficient=(
                        deep_ocean_exchange
                        if deep_ocean_exchange is not None
                        else value_of("ocean_deep_exchange_coefficient")
                    ),
                    initial_temperature=(
                        initial_deep_temperature
                        if initial_deep_temperature is not None
                        else initial_surface_temperature
                    ),
                ),
                step_multiple=1,
            )
        )
        order.append("ocean")
    schedules.extend(
        [
            ModuleSchedule(atmosphere, step_multiple=1),
            ModuleSchedule(surface, step_multiple=1),
            ModuleSchedule(space, step_multiple=1),
        ]
    )
    order.extend(["atmosphere", "surface", "space"])

    schedule = CouplingSchedule(
        master_dt=master_dt,
        schedules=schedules,
        splitting_order=tuple(order),
        convergence_criterion="explicit_single_pass",
    )

    coupling_frequency = (
        "every_master_step"
        if solar_step_multiple == 1
        else f"solar_every_{solar_step_multiple}_master_steps"
    )
    graph = _build_graph(
        configuration_id=run_id,
        coupling_frequency=coupling_frequency,
        convection=convection,
        water=water,
        clouds=cloud is not None,
        deep_ocean=deep_ocean,
        splitting_order=tuple(order),
    )

    initial_air = np.asarray(initial_air_temperature, dtype=float)
    return Kernel(
        run_id=run_id,
        purpose=purpose,
        schedule=schedule,
        states=_build_state_registry(
            n_layers=n_layers, convection=convection, water=water, deep_ocean=deep_ocean
        ),
        equations=build_equation_registry(),
        graph=graph,
        budget=budget
        or ResourceBudget(
            resource_class=ResourceClass.ENGINEERING_WORKSTATION,
            max_memory_bytes=2 * 1024**3,
            max_work_units=1.0e8,
        ),
        initial_state={
            SURFACE_TEMPERATURE_STATE: initial_surface_temperature,
            ATMOSPHERE_TEMPERATURE_STATE: (
                np.full(n_layers, float(initial_air)) if initial_air.ndim == 0 else initial_air
            ),
            SPACE_ENERGY_STATE: 0.0,
            **(
                {
                    OCEAN_DEEP_STATE: (
                        initial_deep_temperature
                        if initial_deep_temperature is not None
                        else initial_surface_temperature
                    )
                }
                if deep_ocean
                else {}
            ),
        },
        conservation_tolerance=conservation_tolerance,
        spectral_bins=working_grid.n_bins,
    )


def transparent_atmosphere(**kwargs: Any) -> Kernel:
    """KER-012 configuration 2: an atmosphere that absorbs and emits nothing.

    The null control for the whole column. With zero opacity the surface must
    reach exactly the airless equilibrium, the layers must never change
    temperature, and every joule entering the top must leave it.
    """
    grid = kwargs.pop("grid", None) or column_grid()
    kwargs.setdefault("run_id", "STAGE2-TRANSPARENT-001")
    kwargs.setdefault(
        "purpose",
        "Stage 2 transparent atmosphere: the column null control, which must "
        "reproduce the airless radiating body exactly.",
    )
    return radiating_column(
        absorption_coefficient=grey_absorption(grid, 0.0),
        grid=grid,
        **kwargs,
    )


def absorbing_emitting_atmosphere(kappa: float = 1.0e-3, **kwargs: Any) -> Kernel:
    """KER-012 configuration 3: a grey absorbing and emitting column.

    ``kappa`` is a mass absorption coefficient in m^2 kg^-1 applied in every
    band. It is a free parameter of the reduced model, not a measured property of
    any gas, and the equation registry classifies it accordingly.
    """
    grid = kwargs.pop("grid", None) or column_grid()
    kwargs.setdefault("run_id", f"STAGE2-GREY-{kappa:g}")
    kwargs.setdefault(
        "purpose",
        "Stage 2 absorbing and emitting grey atmosphere: multilayer energy "
        "closure and the sign and magnitude of the greenhouse response.",
    )
    return radiating_column(
        absorption_coefficient=grey_absorption(grid, kappa),
        grid=grid,
        **kwargs,
    )


def window_atmosphere(
    kappa_opaque: float = 1.0e-4, kappa_window: float = 0.0, **kwargs: Any
) -> Kernel:
    """A spectrally structured column: transparent to sunlight, open in the window.

    Exists so SPC-009 can be tested against something other than a grey
    coefficient -- the same integrated opacity distributed differently must be
    allowed to give a different answer.

    The default ``kappa_opaque`` gives a column optical depth of order unity,
    which is Earth-like in magnitude. It is a free parameter of the reduced
    model, not a measured property of any gas, and nothing here claims otherwise.
    """
    grid = kwargs.pop("grid", None) or column_grid()
    kwargs.setdefault("run_id", "STAGE2-WINDOW-001")
    kwargs.setdefault(
        "purpose",
        "Stage 2 spectrally structured atmosphere: an infrared window changes the "
        "response at equal integrated opacity.",
    )
    return radiating_column(
        absorption_coefficient=window_absorption(
            grid, kappa_opaque=kappa_opaque, kappa_window=kappa_window
        ),
        grid=grid,
        **kwargs,
    )


def single_grey_layer(emissivity: float = 0.9, **kwargs: Any) -> Kernel:
    """One grey layer opaque in the infrared and transparent to sunlight.

    The configuration the analytic greenhouse identity applies to. The requested
    layer emissivity is converted to the mass absorption coefficient that
    produces it, so the test compares the solver against the closed form rather
    than against a coefficient chosen to make them agree.
    """
    grid = kwargs.pop("grid", None) or column_grid()
    n_layers = 1
    surface_pressure = float(kwargs.get("surface_pressure") or 101325.0)
    top_pressure = float(kwargs.get("top_pressure", 1.0e3))
    gravity = 9.80665
    layer_mass = (surface_pressure - top_pressure) / gravity
    diffusivity = float(kwargs.get("diffusivity_factor") or 1.66)

    if not 0.0 <= emissivity < 1.0:
        raise ValueError(
            "a layer emissivity of exactly 1 needs infinite opacity; use a value below 1"
        )
    kappa_lw = -np.log1p(-emissivity) / (diffusivity * layer_mass)

    # Transparent to sunlight, opaque in the thermal infrared: the spectral
    # separation is what makes this a greenhouse rather than a uniform shade.
    kappa = np.where(grid.centers < 4.0e-6, 0.0, kappa_lw)

    kwargs.setdefault("run_id", f"STAGE2-GREY-LAYER-EPS{emissivity:g}")
    kwargs.setdefault(
        "purpose",
        "Stage 2 single grey layer: comparison against the closed-form greenhouse "
        "temperatures for a one-layer atmosphere.",
    )
    kwargs["n_layers"] = n_layers
    kwargs["surface_pressure"] = surface_pressure
    kwargs["top_pressure"] = top_pressure
    kwargs["diffusivity_factor"] = diffusivity
    return radiating_column(
        absorption_coefficient=kappa,
        grid=grid,
        **kwargs,
    )


def radiative_convective(kappa: float = 1.0e-4, **kwargs: Any) -> Kernel:
    """Stage 3: a radiating column with dry convective adjustment.

    The correction to Stage 2. Left to radiation alone the column reaches pure
    radiative equilibrium, which is convectively unstable in its lower layers;
    with adjustment the profile is stable everywhere and the surface is cooler,
    because convection carries heat upward that radiation alone had to leave at
    the bottom.
    """
    grid = kwargs.pop("grid", None) or column_grid()
    kwargs.setdefault("run_id", f"STAGE3-RADCONV-{kappa:g}")
    kwargs.setdefault(
        "purpose",
        "Stage 3 radiative-convective equilibrium: a dry column whose lapse rate "
        "is limited by convective adjustment.",
    )
    kwargs["convection"] = True
    return radiating_column(
        absorption_coefficient=window_absorption(grid, kappa_opaque=kappa),
        grid=grid,
        **kwargs,
    )


def moist_radiative_convective(kappa: float = 1.0e-4, **kwargs: Any) -> Kernel:
    """Stage 3: the same column with water vapour, evaporation and precipitation.

    Two conserved quantities move at once, which is the point. Water leaves the
    surface as vapour carrying latent energy, condenses aloft releasing that
    energy as sensible heat, and returns as precipitation carrying none. Both
    budgets have to close simultaneously, which is what CON-007 asks of a phase
    change.

    The vapour is radiatively inert: the absorption coefficient does not respond
    to it, so this configuration has no water-vapour feedback and no claim about
    one is in scope.
    """
    grid = kwargs.pop("grid", None) or column_grid()
    kwargs.setdefault("run_id", f"STAGE3-MOIST-{kappa:g}")
    kwargs.setdefault(
        "purpose",
        "Stage 3 moist radiative-convective column: simultaneous closure of the "
        "water and energy budgets across a phase change.",
    )
    kwargs["convection"] = True
    kwargs["water"] = True
    return radiating_column(
        absorption_coefficient=window_absorption(grid, kappa_opaque=kappa),
        grid=grid,
        **kwargs,
    )


def no_cloud_control(kappa: float = 1.0e-4, **kwargs: Any) -> Kernel:
    """KER-012 configuration: the no-cloud control.

    The reference every cloud result is measured against. It is the cloudy
    configuration with the cloud removed, not a separately written column, so a
    difference between the two is attributable to the cloud and to nothing else.
    """
    kwargs.setdefault("run_id", "STAGE4-NO-CLOUD-001")
    kwargs.setdefault(
        "purpose",
        "Stage 4 no-cloud control: the clear-sky reference for every cloud "
        "radiative effect reported by the prescribed-cloud configuration.",
    )
    kwargs["cloud"] = None
    return radiative_convective(kappa, **kwargs)


def prescribed_cloud(
    kappa: float = 1.0e-4,
    *,
    layer: int = 4,
    liquid_water_path: float = 0.05,
    fraction: float = 0.6,
    droplet_effective_radius: float = 1.0e-5,
    ice: bool = False,
    **kwargs: Any,
) -> Kernel:
    """KER-012 configuration: a prescribed cloud, static or time-varying.

    Pass a callable as ``cloud`` for the time-varying case; the defaults give a
    single mid-level liquid cloud covering 60 percent of the column.

    The cloud's water is supplied by fiat rather than drawn from the water cycle,
    so a run using this does not conserve water through its cloud. The cloud
    module's budget names that source rather than hiding it.
    """
    n_layers = int(kwargs.get("n_layers", 8))
    kwargs.setdefault("run_id", f"STAGE4-CLOUD-LWP{liquid_water_path:g}-F{fraction:g}")
    kwargs.setdefault(
        "purpose",
        "Stage 4 prescribed cloud: shortwave and longwave cloud radiative effects "
        "reported separately against the no-cloud control.",
    )
    kwargs.setdefault(
        "cloud",
        single_layer_cloud(
            n_layers=n_layers,
            layer=layer,
            liquid_water_path=liquid_water_path,
            fraction=fraction,
            ice=ice,
            droplet_effective_radius=droplet_effective_radius,
        ),
    )
    return radiative_convective(kappa, **kwargs)


def column_over_reduced_ocean(kappa: float = 1.0e-4, **kwargs: Any) -> Kernel:
    """A moist radiative-convective column standing over the two-reservoir ocean.

    The same ``DeepOceanModule`` that sits under the airless slab, under a very
    different surface. That is what putting the reservoir in its own module
    buys: it couples to whatever owns the mixed layer, and neither surface
    implementation had to learn anything about the other.

    Here the feedback parameter is not the bare ``4 sigma T^3`` of a radiating
    body -- the atmosphere, the convective adjustment and the surface water all
    sit inside the loop -- so the analytic two-mode expression does not apply
    and the timescales must be measured rather than predicted. That is the
    ordinary case. The airless configuration exists to check the measurement
    against the one case where a prediction is available.
    """
    kwargs.setdefault("run_id", "STAGE5-COLUMN-OCEAN-001")
    kwargs.setdefault(
        "purpose",
        "Stage 5 column over a reduced ocean: a moist radiative-convective column whose "
        "mixed layer exchanges heat with a deep reservoir, for the response timescales "
        "of a configuration whose feedback parameter is not known in closed form.",
    )
    kwargs["deep_ocean"] = True
    return moist_radiative_convective(kappa, **kwargs)
