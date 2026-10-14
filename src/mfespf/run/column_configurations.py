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
from mfespf.core.errors import ConfigurationError
from mfespf.coupling.coupler import CouplingSchedule, ModuleSchedule
from mfespf.coupling.graph import Cycle, DependencyGraph, Edge, EdgeType, Node, NodeType
from mfespf.physics.atmosphere import (
    ATMOSPHERE_TEMPERATURE_STATE,
    AtmosphericColumnModule,
    atmosphere_state_entry,
)
from mfespf.physics.boundary import SPACE_ENERGY_STATE, SpaceBoundaryModule, space_state_entry
from mfespf.physics.cloud_field import (
    ISCCP_SURFACE_REFLECTANCE,
    CloudType,
    ObservedCloudFieldModule,
)
from mfespf.physics.clouds import CloudProfile, PrescribedCloudModule, single_layer_cloud
from mfespf.physics.co2_band import co2_band_grid
from mfespf.physics.column_surface import ColumnSurfaceModule
from mfespf.physics.convection import ConvectiveAdjustmentModule
from mfespf.physics.equations import EQUATION_IDS, build_equation_registry
from mfespf.physics.gravity import GravityField
from mfespf.physics.greenhouse import CALIBRATION_N_LAYERS
from mfespf.physics.mass_flux import (
    CLOUD_WATER_CONVERSION_PER_M,
    CONVECTIVE_ADJUSTMENT_TIME_S,
    PLUME_ENTRAINMENT_PER_M,
)
from mfespf.physics.ocean import (
    OCEAN_DEEP_STATE,
    DeepOceanModule,
    ocean_deep_state_entry,
)
from mfespf.physics.orbit import OrbitalGeometryModule
from mfespf.physics.ozone_infrared import ozone_infrared_grid
from mfespf.physics.plume import (
    BULK_PLUME_PARAMETER,
    EVAPORATING_CONDENSATE_FRACTION,
)
from mfespf.physics.radiative_transfer import (
    SpectralRadiationModule,
    grey_absorption,
    window_absorption,
)
from mfespf.physics.rrtmg_longwave import RRTMGLongwave, rrtmg_longwave_grid
from mfespf.physics.solar import (
    ConstantSpectrumSolarSource,
    MeasuredSpectrumSolarSource,
    banded_grid,
)
from mfespf.physics.solar_absorption import ozone_bands, solar_vapour_bands
from mfespf.physics.surface import (
    SURFACE_TEMPERATURE_STATE,
    SURFACE_WATER_STATE,
    surface_state_entry,
    surface_water_state_entry,
)
from mfespf.physics.thermodynamics import evaporation_coefficient_for
from mfespf.physics.water import (
    VAPOUR_STATE,
    EntrainingPlumeHumidityModule,
    PrescribedHumidityModule,
    WaterModule,
    vapour_state_entry,
)
from mfespf.physics.water_vapour_longwave import ChouLongwave, chou_longwave_grid
from mfespf.registries.datasets import DatasetRegistry
from mfespf.registries.states import StateRegistry
from mfespf.run.driver import Kernel
from mfespf.run.modes import OperatingMode

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
    vapour_transport: bool = False,
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
                owner="water",
                n_layers=n_layers,
                readers=frozenset({"radiation", "atmosphere"}),
                # Convection redistributes the vapour it carries, so it submits a
                # tendency the water module commits -- the same arrangement the
                # surface temperature already has.
                tendency_providers=frozenset({"convection"}) if vapour_transport else frozenset(),
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
    deep_ocean: bool = False,
    clouds: bool = False,
    cloud_field: bool = False,
    greenhouse: bool = False,
    volcanic: bool = False,
    radiative_vapour: bool = False,
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
    if cloud_field:
        modules.append(("cloud_field", "observed global cloud field (ISCCP)"))
    if deep_ocean:
        modules.append(("ocean", "reduced deep-ocean reservoir"))
    if greenhouse:
        modules.append(("greenhouse", "prescribed carbon dioxide"))
    if volcanic:
        modules.append(("volcanic", "prescribed stratospheric aerosol"))
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
            "surface",
            "water",
            "deposition_fusion_energy",
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
        if radiative_vapour:
            # The feedback, as an edge. Without it the vapour is tracked and
            # radiatively inert, which is what every configuration before this
            # one did; with it the water the column carries changes what the
            # column radiates, and the loop through temperature closes.
            edge(
                VAPOUR_STATE,
                "radiation",
                VAPOUR_STATE,
                edge_type=EdgeType.STATE_DEPENDENCY,
                units="kg m^-2",
                equation_ids=("EQ-VAPOUR-BAND-ABSORPTION", "EQ-BEER-LAMBERT-EXTINCTION"),
                uncertainty_rule=(
                    "vapour path uncertainty propagates into the window's optical depth"
                ),
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

    if cloud_field:
        for exchanged, units in (
            ("cloud_type_amount", "1"),
            ("cloud_type_top_pressure", "Pa"),
            ("cloud_type_visible_optical_thickness", "1"),
            ("cloud_type_infrared_optical_thickness", "1"),
        ):
            edge(
                "cloud_field",
                "radiation",
                exchanged,
                edge_type=EdgeType.BOUNDARY_DEPENDENCY,
                units=units,
                equation_ids=("EQ-ISCCP-CLOUD-FIELD",),
                uncertainty_rule=(
                    "observed climatology; responds to nothing, so no feedback path exists"
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

    if greenhouse:
        edge(
            "greenhouse",
            "radiation",
            "co2_mole_fraction",
            edge_type=EdgeType.PARAMETER_DEPENDENCY,
            units="ppm",
            equation_ids=(EQUATION_IDS["co2_band_absorption"],),
            uncertainty_rule=(
                "prescribed; the concentration responds to nothing, so no carbon-cycle "
                "feedback path exists"
            ),
        )

    if volcanic:
        edge(
            "volcanic",
            "radiation",
            "stratospheric_aerosol_optical_depth",
            edge_type=EdgeType.PARAMETER_DEPENDENCY,
            units="1",
            equation_ids=(EQUATION_IDS["volcanic_aerosol_optics"],),
            uncertainty_rule=(
                "prescribed; eruptions are an input and nothing here explains their timing"
            ),
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
    n_layers: int = CALIBRATION_N_LAYERS,
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
    radiative_vapour: bool = False,
    vapour_transport: bool = False,
    surface_based_convection: bool = False,
    convective_target: str = "dry_adiabat",
    mass_flux_convection: bool = False,
    plume_entrainment_rate: float | None = None,
    convective_adjustment_time: float | None = None,
    detrained_state: str = "plume",
    precipitation_evaporation: str = "cascade",
    condensate_pathway: str = "immediate",
    cloud_water_conversion_rate: float | None = None,
    diagnose_convective_recycling: bool = False,
    track_water_provenance: bool = False,
    convective_substep: float | None = None,
    convective_window_continuation: str = "stop",
    convective_trigger_crossing: str = "explicit",
    prescribed_humidity: float | None = None,
    plume_humidity: bool = False,
    plume_vertical_structure: bool = False,
    bulk_plume_parameter: float = BULK_PLUME_PARAMETER,
    evaporating_fraction: float = EVAPORATING_CONDENSATE_FRACTION,
    humidity_responds_to_temperature: bool = True,
    cloud: CloudProfile | None = None,
    greenhouse: Any = None,
    volcanic: Any = None,
    solar_source: Any = None,
    start_time: float = 0.0,
    operating_mode: OperatingMode = OperatingMode.CONTROLLED_EXPERIMENT,
    datasets: DatasetRegistry | None = None,
    deep_ocean: bool = False,
    deep_ocean_depth: float = 1200.0,
    deep_ocean_exchange: float | None = None,
    initial_deep_temperature: float | None = None,
    sensible_heat_coefficient: float = 10.0,
    convective_relaxation_timescale: float = 0.0,
    convective_fractional_top: bool = False,
    solar_vapour_absorption: bool = False,
    measured_solar_spectrum: bool = False,
    ozone: bool = False,
    ozone_above_1hpa: str = "none",
    co2_band: bool = False,
    co2_ppmv: float = 284.7,
    ozone_infrared: bool = False,
    chou_water_vapour: bool = False,
    rrtmg_longwave: bool = False,
    k_distribution_source: bool = False,
    linear_in_tau_source: bool = False,
    observed_cloud_field: bool = False,
    cloud_types: tuple[CloudType, ...] | None = None,
    receives_advection: bool = False,
    eddy_diffusivity: float = 0.0,
    condense_after_convection: bool = False,
    end_of_step_saturation: bool = False,
    max_log_step: float | None = None,
    geometry_factor: float | None = None,
    evaporation_coefficient: float | None = None,
    gravity: float | GravityField | None = None,
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
    solar_vapour_k = None
    ozone_k = None
    # Ozone absorbs where a black-body sun is half again too bright, so it is
    # only ever run on the measured spectrum (ozone-solar-heating.md).
    measured = measured_solar_spectrum or ozone

    def carry(original: SpectralGrid, rebuilt: SpectralGrid) -> np.ndarray:
        """Each rebuilt band's home band on the original grid."""
        return np.asarray(
            np.clip(
                np.searchsorted(original.edges, rebuilt.centers, side="right") - 1,
                0,
                original.n_bins - 1,
            )
        )

    if solar_vapour_absorption:
        # Water vapour absorbs sunlight (Chou and Suarez 1999). The solar region
        # of the grid is rebuilt around the source's bands; everything given per
        # band on the original grid is carried onto the new one from the band
        # that contains each new band, which is exact for the thermal bands
        # (kept edge for edge) and for anything flat inside a solar band.
        original = working_grid
        working_grid, solar_vapour_k = solar_vapour_bands(original, measured_spectrum=measured)
        home = np.clip(
            np.searchsorted(original.edges, working_grid.centers, side="right") - 1,
            0,
            original.n_bins - 1,
        )
        if kappa.ndim >= 1 and kappa.shape[-1] == original.n_bins:
            kappa = kappa[..., home]
        if albedo is not None and np.ndim(albedo) >= 1 and np.shape(albedo)[-1] == original.n_bins:
            albedo = np.asarray(albedo, dtype=float)[..., home]
    if chou_water_vapour and not (co2_band and radiative_vapour):
        raise ConfigurationError(
            "Chou's water vapour shares band 3 with the CO2 band and absorbs only "
            "radiative vapour; switch on co2_band and radiative_vapour with it"
        )
    if ozone_infrared and not ozone:
        raise ConfigurationError(
            "ozone's infrared band uses the ozone profile; switch ozone on with it"
        )
    if rrtmg_longwave and (chou_water_vapour or co2_band or ozone_infrared):
        raise ConfigurationError(
            "RRTMG's longwave carries the vapour, CO2 and ozone itself; it replaces "
            "chou_water_vapour, co2_band and ozone_infrared rather than joining them"
        )
    if rrtmg_longwave and not (radiative_vapour and ozone):
        raise ConfigurationError(
            "RRTMG's longwave needs the radiative vapour and the ozone profile; switch on "
            "radiative_vapour and ozone with it"
        )
    rrtmg_spec: RRTMGLongwave | None = None
    if rrtmg_longwave:
        # RRTMG's gas optics (rrtmg-longwave.md): the longwave rebuilt on its
        # g-points; everything per band below 4 microns is carried by home band.
        original = working_grid
        working_grid, gpoint = rrtmg_longwave_grid(original)
        home = carry(original, working_grid)
        kappa = np.zeros(working_grid.n_bins)
        if albedo is not None and np.ndim(albedo) >= 1 and np.shape(albedo)[-1] == original.n_bins:
            albedo = np.asarray(albedo, dtype=float)[..., home]
        if solar_vapour_k is not None:
            solar_vapour_k = solar_vapour_k[home]
        rrtmg_spec = RRTMGLongwave(gpoint=gpoint, co2_ppmv=float(co2_ppmv))
    co2_band_spec = None
    chou_spec: ChouLongwave | None = None
    ozone_ir_terms = None
    if chou_water_vapour:
        # Chou et al.'s water vapour (global-column-chou-longwave.md): the whole
        # longwave is rebuilt on its bands, with the CO2 terms (and ozone's, if
        # on) split by the water terms they share a band with. It replaces both
        # co2_band_grid and ozone_infrared_grid, and the fitted envelope.
        original = working_grid
        working_grid, chou_spec = chou_longwave_grid(original, ozone=ozone_infrared)
        home = carry(original, working_grid)
        kappa = np.zeros(working_grid.n_bins)
        if albedo is not None and np.ndim(albedo) >= 1 and np.shape(albedo)[-1] == original.n_bins:
            albedo = np.asarray(albedo, dtype=float)[..., home]
        if solar_vapour_k is not None:
            solar_vapour_k = solar_vapour_k[home]
        co2_band_spec = (chou_spec.co2_k, chou_spec.co2_kind, float(co2_ppmv))
        if ozone_infrared:
            ozone_ir_terms = chou_spec.ozone_term
    elif co2_band:
        # Chou et al.'s CO2 band replaces the window caricature entirely: with a
        # real band, the caricature's stand-in for it would count CO2 twice
        # (co2-band.md). The region 12.5-18.5 microns is rebuilt around its
        # terms, and everything else per band is carried by home band.
        original = working_grid
        working_grid, co2_k, co2_kind = co2_band_grid(original)
        home = carry(original, working_grid)
        kappa = np.zeros(working_grid.n_bins)
        if albedo is not None and np.ndim(albedo) >= 1 and np.shape(albedo)[-1] == original.n_bins:
            albedo = np.asarray(albedo, dtype=float)[..., home]
        if solar_vapour_k is not None:
            solar_vapour_k = solar_vapour_k[home]
        co2_band_spec = (co2_k, co2_kind, float(co2_ppmv))
    if ozone_infrared and chou_spec is None:
        # Ozone's 9.6-micron band (ozone-infrared-band.md), on the SWOOSH
        # profile the solar ozone already places, so it needs ozone on.
        original = working_grid
        working_grid, ozone_ir_terms = ozone_infrared_grid(original)
        home = carry(original, working_grid)
        if kappa.ndim >= 1 and kappa.shape[-1] == original.n_bins:
            kappa = kappa[..., home]
        if albedo is not None and np.ndim(albedo) >= 1 and np.shape(albedo)[-1] == original.n_bins:
            albedo = np.asarray(albedo, dtype=float)[..., home]
        if solar_vapour_k is not None:
            solar_vapour_k = solar_vapour_k[home]
        if co2_band_spec is not None:
            co2_band_spec = (co2_band_spec[0][home], co2_band_spec[1][home], co2_band_spec[2])
    if ozone:
        # The same carrying, for ozone's rebuild of the region below 0.4 microns.
        # The vapour's sub-bands above it are kept edge for edge.
        original = working_grid
        working_grid, ozone_k = ozone_bands(original)
        home = carry(original, working_grid)
        if kappa.ndim >= 1 and kappa.shape[-1] == original.n_bins:
            kappa = kappa[..., home]
        if albedo is not None and np.ndim(albedo) >= 1 and np.shape(albedo)[-1] == original.n_bins:
            albedo = np.asarray(albedo, dtype=float)[..., home]
        if solar_vapour_k is not None:
            solar_vapour_k = solar_vapour_k[home]
        if co2_band_spec is not None:
            co2_band_spec = (co2_band_spec[0][home], co2_band_spec[1][home], co2_band_spec[2])
        if ozone_ir_terms is not None:
            ozone_ir_terms = ozone_ir_terms[home]
        if chou_spec is not None:
            chou_spec = chou_spec.carried(home)
        if rrtmg_spec is not None:
            rrtmg_spec = rrtmg_spec.carried(home)
    # Default to a surface that reflects sunlight and is black in the infrared.
    # A plain float would be grey in every band, which is a different -- and
    # physically distinguishable -- surface; see ``shortwave_albedo``.
    # With observed clouds the surface is only the surface: ISCCP's reflectance,
    # 0.17, not the planetary 0.30 that stood in for the clouds. The two go
    # together or the clouds' reflection is counted twice.
    surface_reflectance: float | np.ndarray
    if albedo is None:
        surface_reflectance = shortwave_albedo(
            working_grid, ISCCP_SURFACE_REFLECTANCE if observed_cloud_field else 0.3
        )
    else:
        surface_reflectance = albedo

    solar = solar_source or (
        MeasuredSpectrumSolarSource(
            grid=working_grid, total_irradiance=total_irradiance, status=ModuleStatus.PRESCRIBED
        )
        if measured
        else ConstantSpectrumSolarSource(
            grid=working_grid, total_irradiance=total_irradiance, status=ModuleStatus.PRESCRIBED
        )
    )
    orbit = OrbitalGeometryModule(
        eccentricity=eccentricity, receiver_module_id="atmosphere", geometry_factor=geometry_factor
    )
    radiation = SpectralRadiationModule(
        grid=working_grid,
        absorption_coefficient=kappa,
        diffusivity_factor=diffusivity_factor,
        couple_clouds=cloud is not None,
        couple_greenhouse=greenhouse is not None,
        couple_volcanic=volcanic is not None,
        couple_water=radiative_vapour,
        solar_vapour_k=solar_vapour_k,
        ozone_k=ozone_k,
        co2_band=co2_band_spec,
        ozone_infrared=ozone_ir_terms,
        water_longwave=chou_spec,
        rrtmg=rrtmg_spec,
        k_distribution_source=k_distribution_source,
        linear_in_tau_source=linear_in_tau_source,
        couple_cloud_field=observed_cloud_field,
        ozone_above_1hpa=ozone_above_1hpa,
    )
    atmosphere = AtmosphericColumnModule(
        n_layers=n_layers,
        surface_pressure=surface_pressure,
        top_pressure=top_pressure,
        initial_temperature=initial_air_temperature,
        couple_water=water,
        gravity=gravity,
        max_log_step=max_log_step,
    )
    # With a log-pressure cap the column has more layers than the reference
    # count that fixed its lower step; every module is sized by the real one.
    n_layers = atmosphere.n_layers
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
    if greenhouse is not None:
        # Published before radiation reads it: a prescribed concentration is a
        # boundary condition on the transfer, not a response to it.
        schedules.append(ModuleSchedule(greenhouse, step_multiple=1))
        order_head.append("greenhouse")
    if volcanic is not None:
        schedules.append(ModuleSchedule(volcanic, step_multiple=1))
        order_head.append("volcanic")
    if observed_cloud_field:
        if cloud is not None:
            raise ConfigurationError(
                "a column takes the observed cloud field or a prescribed cloud, not both: "
                "each would count the other's reflection"
            )
        schedules.append(
            ModuleSchedule(
                ObservedCloudFieldModule(cloud_types=cloud_types)
                if cloud_types is not None
                else ObservedCloudFieldModule(),
                step_multiple=1,
            )
        )
        order_head.append("cloud_field")
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
    humidity_sources = sum((bool(water), prescribed_humidity is not None, bool(plume_humidity)))
    if condense_after_convection and not (
        water
        and convection
        and vapour_transport
        and prescribed_humidity is None
        and not plume_humidity
    ):
        raise ConfigurationError(
            "condensing after convection needs the water module, convection and "
            "convection's vapour transport: it composes the one with the other"
        )
    if humidity_sources > 1:
        raise ConfigurationError(
            "a column runs the water module, a prescribed humidity profile or the "
            "bulk-plume closure, never more than one: they publish the same field, and "
            "two producers of one field is an ownership question the coupler cannot "
            "answer"
        )
    if plume_humidity:
        # The humidity is predicted here rather than chosen. It publishes the
        # same field the prescribed control does, so a run against that control
        # differs in one thing only: whether the number came from the plume's
        # water budget or from Manabe and Wetherald.
        schedules.append(
            ModuleSchedule(
                EntrainingPlumeHumidityModule(
                    n_layers=n_layers,
                    bulk_plume_parameter=bulk_plume_parameter,
                    evaporating_fraction=evaporating_fraction,
                    vertical_structure=plume_vertical_structure,
                    respond_to_temperature=humidity_responds_to_temperature,
                ),
                step_multiple=1,
            )
        )
        order.append("water")
    elif prescribed_humidity is not None:
        # The control for the feedback. It publishes the same field the water
        # module does, so radiation needs no second input and the only
        # difference between the two configurations is where the humidity
        # distribution came from -- which is the thing being controlled for.
        schedules.append(
            ModuleSchedule(
                PrescribedHumidityModule(
                    n_layers=n_layers,
                    relative_humidity=prescribed_humidity,
                    respond_to_temperature=humidity_responds_to_temperature,
                ),
                step_multiple=1,
            )
        )
        order.append("water")
    elif water:
        schedules.append(
            ModuleSchedule(
                WaterModule(
                    n_layers=n_layers,
                    evaporation_coefficient=evaporation_coefficient
                    if evaporation_coefficient is not None
                    else evaporation_coefficient_for(sensible_heat_coefficient),
                    evaporation_coefficient_follows=(
                        "" if evaporation_coefficient is not None else "sensible_heat_coefficient"
                    ),
                    receives_advection=receives_advection,
                    eddy_diffusivity=eddy_diffusivity,
                    receives_convective_transport=condense_after_convection,
                    end_of_step_saturation=end_of_step_saturation,
                    receives_convective_condensation=mass_flux_convection,
                    diagnose_convective_recycling=diagnose_convective_recycling,
                    track_water_provenance=track_water_provenance,
                ),
                step_multiple=1,
            )
        )
        order.append("water")
    if convection:
        schedules.append(
            ModuleSchedule(
                ConvectiveAdjustmentModule(
                    sensible_heat_coefficient=sensible_heat_coefficient,
                    transport_vapour=vapour_transport,
                    surface_based=surface_based_convection,
                    surface_pressure=surface_pressure,
                    target=convective_target,
                    bulk_plume_parameter=bulk_plume_parameter,
                    # One switch, because holding the vapour is one physical
                    # intervention even though it now has two consequences.
                    target_responds_to_temperature=humidity_responds_to_temperature,
                    relaxation_timescale=convective_relaxation_timescale,
                    fractional_top=convective_fractional_top,
                    receives_advection=receives_advection,
                    condensation_follows=condense_after_convection,
                    mass_flux=mass_flux_convection,
                    plume_entrainment_rate=(
                        plume_entrainment_rate
                        if plume_entrainment_rate is not None
                        else PLUME_ENTRAINMENT_PER_M
                    ),
                    convective_adjustment_time=(
                        convective_adjustment_time
                        if convective_adjustment_time is not None
                        else CONVECTIVE_ADJUSTMENT_TIME_S
                    ),
                    detrained_state=detrained_state,
                    precipitation_evaporation=precipitation_evaporation,
                    condensate_pathway=condensate_pathway,
                    track_water_provenance=track_water_provenance,
                    convective_substep=convective_substep,
                    convective_window_continuation=convective_window_continuation,
                    convective_trigger_crossing=convective_trigger_crossing,
                    cloud_water_conversion_rate=(
                        cloud_water_conversion_rate
                        if cloud_water_conversion_rate is not None
                        else CLOUD_WATER_CONVERSION_PER_M
                    ),
                ),
                step_multiple=1,
            )
        )
        if condense_after_convection:
            # The opposite of the default order, and for the opposite reason:
            # the vapour convection lifts is what condenses, so the saturation
            # adjustment follows the lifting
            # (docs/protocols/convective-condensation-order.md).
            order.insert(order.index("water"), "convection")
        else:
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
        water=water or prescribed_humidity is not None,
        clouds=cloud is not None,
        cloud_field=observed_cloud_field,
        deep_ocean=deep_ocean,
        greenhouse=greenhouse is not None,
        volcanic=volcanic is not None,
        radiative_vapour=radiative_vapour,
        splitting_order=tuple(order),
    )

    initial_air = np.asarray(initial_air_temperature, dtype=float)
    return Kernel(
        start_time=start_time,
        operating_mode=operating_mode,
        datasets=datasets or DatasetRegistry(),
        run_id=run_id,
        purpose=purpose,
        schedule=schedule,
        states=_build_state_registry(
            n_layers=n_layers,
            convection=convection,
            water=water,
            deep_ocean=deep_ocean,
            vapour_transport=vapour_transport,
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


def reference_column(kappa: float = 1.0e-4, **kwargs: Any) -> Kernel:
    """The column the headline error is quoted from: every mechanism it has.

    The Stage 3 moist column with its vapour radiatively active and carried by
    convection, and with the convection relaxing toward the profile an
    entraining plume leaves rather than toward the dry adiabat. That last
    choice is physics, not a setting: air rising through a moist column
    condenses its water, and the latent heat released warms the air aloft. A
    column that condenses but relaxes to the dry adiabat has the heat and not
    its effect on the profile -- measured at 8.5 K/km against an observed 6.5,
    worth about six kelvin of surface temperature.

    Its closure predicts two things that can be checked against observations
    other than the surface temperature: the lapse rate (5.5 K/km, shallower
    than the global-mean 6.5, because a global mean is not all convecting) and
    the free-tropospheric humidity (0.74, wetter than the observed 0.4 to 0.5,
    because nothing in one column dries air by making it sink).
    """
    kwargs.setdefault("run_id", f"REFERENCE-COLUMN-{kappa:g}")
    kwargs.setdefault(
        "purpose",
        "Reference column: full water cycle with moist convection, the "
        "configuration the headline surface-temperature error is quoted from.",
    )
    kwargs.setdefault("radiative_vapour", True)
    kwargs.setdefault("vapour_transport", True)
    kwargs.setdefault("surface_based_convection", True)
    kwargs.setdefault("convective_target", "entraining_plume")
    return moist_radiative_convective(kappa, **kwargs)


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
    layer: int | None = None,
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
    n_layers = int(kwargs.get("n_layers", CALIBRATION_N_LAYERS))
    # Mid-level at any resolution: the middle layer, which is layer 4 of 8. A
    # fixed index would move the cloud toward the ground as the grid refines.
    if layer is None:
        layer = n_layers // 2
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
