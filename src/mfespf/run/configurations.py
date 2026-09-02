"""Reference configurations for the validation ladder.

Charter section 19 stages 0-1 and requirement KER-012.

Each builder returns a fully assembled ``Kernel``: modules, registries, coupling
schedule, dependency graph and resource budget.  Configurations are functions
rather than data files so that the declarations and the objects they describe
cannot drift apart; a YAML front end can call these and record what it passed.
"""

from __future__ import annotations

import numpy as np

from mfespf.contracts.cost import ResourceBudget, ResourceClass
from mfespf.contracts.module import ModuleStatus
from mfespf.contracts.spectral import SpectralGrid
from mfespf.coupling.coupler import CouplingSchedule, ModuleSchedule
from mfespf.coupling.graph import DependencyGraph, Edge, EdgeType, Node, NodeType
from mfespf.physics.boundary import SPACE_ENERGY_STATE, SpaceBoundaryModule, space_state_entry
from mfespf.physics.equations import build_equation_registry
from mfespf.physics.orbit import OrbitalGeometryModule
from mfespf.physics.solar import ConstantSpectrumSolarSource, banded_grid, broadband_grid
from mfespf.physics.surface import (
    SURFACE_TEMPERATURE_STATE,
    SlabSurfaceModule,
    surface_state_entry,
)
from mfespf.registries.states import StateRegistry
from mfespf.run.driver import Kernel

__all__ = [
    "airless_radiating_body",
    "build_dependency_graph",
    "build_state_registry",
    "null_control",
]

# A 50 m ocean mixed layer: rho * c_p * depth, with both properties taken from
# the constants registry rather than folded into one tuned number.
_DEFAULT_AREAL_HEAT_CAPACITY = 1025.0 * 3985.0 * 50.0


def build_state_registry(
    *, surface_owner: str = "surface", space_owner: str = "space"
) -> StateRegistry:
    """Register the Stage 1 prognostic states with their authorisations."""
    states = StateRegistry()
    states.register(surface_state_entry(owner=surface_owner, readers=frozenset({space_owner})))
    states.register(space_state_entry(owner=space_owner))
    return states


def build_dependency_graph(
    *,
    configuration_id: str,
    solar_id: str = "solar",
    orbit_id: str = "orbit",
    surface_id: str = "surface",
    space_id: str = "space",
    coupling_frequency: str = "every_master_step",
    solar_frequency: str = "every_master_step",
) -> DependencyGraph:
    """The implemented dependency graph for the airless-body configuration (CAU-001)."""
    graph = DependencyGraph(configuration_id)
    for module_id, label in (
        (solar_id, "solar source"),
        (orbit_id, "orbital geometry and TOA boundary"),
        (surface_id, "slab surface"),
        (space_id, "system boundary sink"),
    ):
        graph.add_node(Node(module_id, NodeType.MODULE, label))
    graph.add_node(Node(SURFACE_TEMPERATURE_STATE, NodeType.STATE, "slab temperature"))
    graph.add_node(Node(SPACE_ENERGY_STATE, NodeType.STATE, "cumulative outgoing energy"))

    graph.add_edge(
        Edge(
            source=solar_id,
            destination=orbit_id,
            edge_type=EdgeType.BOUNDARY_DEPENDENCY,
            exchanged="solar_total_irradiance_1au",
            units="W m^-2",
            interface_version="1.0.0",
            time_semantics="interval_average",
            spatial_semantics="global_scalar",
            coupling_frequency=solar_frequency,
            equation_ids=("EQ-PLANCK-SPECTRAL-RADIANCE",),
        )
    )
    graph.add_edge(
        Edge(
            source=solar_id,
            destination=orbit_id,
            edge_type=EdgeType.BOUNDARY_DEPENDENCY,
            exchanged="solar_spectral_irradiance_1au",
            units="W m^-2 m^-1",
            interface_version="1.0.0",
            time_semantics="interval_average",
            spectral_semantics="carried_with_its_own_grid",
            coupling_frequency=solar_frequency,
        )
    )
    graph.add_edge(
        Edge(
            source=orbit_id,
            destination=surface_id,
            edge_type=EdgeType.PHYSICAL_TRANSFER,
            exchanged="toa_incident_irradiance",
            units="W m^-2",
            interface_version="1.0.0",
            time_semantics="interval_average",
            spatial_semantics="global_mean_over_planetary_surface",
            coupling_frequency=coupling_frequency,
            equation_ids=("EQ-SPHERICAL-MEAN-INSOLATION", "EQ-INVERSE-SQUARE-IRRADIANCE"),
        )
    )
    graph.add_edge(
        Edge(
            source=orbit_id,
            destination=surface_id,
            edge_type=EdgeType.BOUNDARY_DEPENDENCY,
            exchanged="toa_spectral_irradiance",
            units="W m^-2 m^-1",
            interface_version="1.0.0",
            spectral_semantics="carried_with_its_own_grid",
            coupling_frequency=coupling_frequency,
        )
    )
    for exchanged, equation in (
        ("surface_reflected_shortwave", "EQ-SPECTRAL-SURFACE-ABSORPTION"),
        ("surface_emitted_longwave", "EQ-GREYBODY-EMISSION"),
    ):
        graph.add_edge(
            Edge(
                source=surface_id,
                destination=space_id,
                edge_type=EdgeType.PHYSICAL_TRANSFER,
                exchanged=exchanged,
                units="W m^-2",
                interface_version="1.0.0",
                time_semantics="interval_average",
                spatial_semantics="global_scalar",
                coupling_frequency=coupling_frequency,
                equation_ids=(equation,),
            )
        )
    graph.add_edge(
        Edge(
            source=SURFACE_TEMPERATURE_STATE,
            destination=surface_id,
            edge_type=EdgeType.STATE_DEPENDENCY,
            exchanged="surface_temperature",
            units="K",
            equation_ids=("EQ-GREYBODY-EMISSION",),
            uncertainty_rule="initial condition propagates through the emission term",
        )
    )
    graph.add_edge(
        Edge(
            source=surface_id,
            destination=SURFACE_TEMPERATURE_STATE,
            edge_type=EdgeType.NUMERICAL_DEPENDENCY,
            exchanged="forward_euler_update",
            units="K",
            equation_ids=("EQ-NUM-FORWARD-EULER",),
        )
    )
    graph.add_edge(
        Edge(
            source=space_id,
            destination=SPACE_ENERGY_STATE,
            edge_type=EdgeType.NUMERICAL_DEPENDENCY,
            exchanged="accumulation_update",
            units="J m^-2",
            equation_ids=("EQ-BOUNDARY-ACCUMULATION",),
        )
    )

    # The surface temperature feeds its own emission and is updated from it: a
    # real, single-module feedback loop, declared rather than left implicit.
    from mfespf.coupling.graph import Cycle

    graph.declare_cycle(
        Cycle(
            nodes=(SURFACE_TEMPERATURE_STATE, "surface"),
            coupling_scheme="explicit_forward_euler_within_one_master_step",
            iteration_order=("surface",),
            convergence_criterion=(
                "single explicit pass; stability requires dt < 2 C / (4 eps sigma T^3)"
            ),
        )
    )
    return graph


def airless_radiating_body(
    *,
    run_id: str = "STAGE1-AIRLESS-001",
    spectral_grid: SpectralGrid | None = None,
    albedo: float | np.ndarray = 0.3,
    emissivity: float = 1.0,
    areal_heat_capacity: float = _DEFAULT_AREAL_HEAT_CAPACITY,
    initial_temperature: float = 200.0,
    total_irradiance: float | None = None,
    eccentricity: float = 0.0,
    master_dt: float = 3600.0,
    solar_step_multiple: int = 24,
    conservation_tolerance: float = 1e-12,
    budget: ResourceBudget | None = None,
) -> Kernel:
    """KER-012 configuration 1: an airless sphere in radiative balance.

    Multirate by default (KER-011): the surface advances hourly while the solar
    source and orbital geometry are exchanged daily, with the intervening steps
    holding the last exchanged values.  The hold is declared in the schedule and
    reported in the run warnings, because at non-zero eccentricity it is a real
    approximation and TIM-004 requires it be visible.
    """
    grid = spectral_grid if spectral_grid is not None else broadband_grid()

    solar = ConstantSpectrumSolarSource(
        grid=grid, total_irradiance=total_irradiance, status=ModuleStatus.PRESCRIBED
    )
    orbit = OrbitalGeometryModule(eccentricity=eccentricity)
    surface = SlabSurfaceModule(
        areal_heat_capacity=areal_heat_capacity,
        emissivity=emissivity,
        albedo=albedo,
        initial_temperature=initial_temperature,
    )
    space = SpaceBoundaryModule()

    schedule = CouplingSchedule(
        master_dt=master_dt,
        schedules=[
            ModuleSchedule(solar, step_multiple=solar_step_multiple),
            ModuleSchedule(orbit, step_multiple=solar_step_multiple),
            ModuleSchedule(surface, step_multiple=1),
            ModuleSchedule(space, step_multiple=1),
        ],
        splitting_order=("solar", "orbit", "surface", "space"),
        convergence_criterion="explicit_single_pass",
    )

    solar_frequency = (
        "every_master_step"
        if solar_step_multiple == 1
        else f"every_{solar_step_multiple}_master_steps"
    )
    graph = build_dependency_graph(configuration_id=run_id, solar_frequency=solar_frequency)

    return Kernel(
        run_id=run_id,
        purpose=(
            "Stage 1 airless radiating body: analytic radiative equilibrium, transient "
            "relaxation, and closure of the coupled energy budget."
        ),
        schedule=schedule,
        states=build_state_registry(),
        equations=build_equation_registry(),
        graph=graph,
        budget=budget
        or ResourceBudget(
            resource_class=ResourceClass.ENGINEERING_WORKSTATION,
            max_memory_bytes=1 * 1024**3,
            max_work_units=1.0e8,
        ),
        initial_state={
            SURFACE_TEMPERATURE_STATE: initial_temperature,
            SPACE_ENERGY_STATE: 0.0,
        },
        conservation_tolerance=conservation_tolerance,
        spectral_bins=grid.n_bins,
    )


def null_control(**kwargs: object) -> Kernel:
    """KER-012 null control: the same configuration with the star switched off.

    VAL-004 and MOD-014 require a disabled-coupling control.  Here the incident
    irradiance is zero, so the surface can only cool, the outgoing energy is
    exactly what the surface emits, and any temperature rise is a bug rather than
    a result.
    """
    kwargs.setdefault("run_id", "STAGE1-NULL-CONTROL-001")
    kwargs.setdefault("total_irradiance", 0.0)
    kernel = airless_radiating_body(**kwargs)  # type: ignore[arg-type]
    kernel.purpose = "Null control: zero solar input, surface may only cool."
    return kernel


def banded_variant(n_bands: int = 16, **kwargs: object) -> Kernel:
    """The same physics at a different spectral resolution (KER-009).

    Interchangeable with ``airless_radiating_body`` in every respect except the
    spectral grid, which is what makes the comparison a test of the spectral
    machinery rather than of two different models.
    """
    kwargs.setdefault("run_id", f"STAGE1-AIRLESS-BANDED-{n_bands}")
    kwargs["spectral_grid"] = banded_grid(n_bands)
    return airless_radiating_body(**kwargs)  # type: ignore[arg-type]
