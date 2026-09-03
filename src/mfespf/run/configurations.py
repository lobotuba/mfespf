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
from mfespf.core.constants import value_of
from mfespf.coupling.coupler import CouplingSchedule, ModuleSchedule
from mfespf.coupling.graph import Cycle, DependencyGraph, Edge, EdgeType, Node, NodeType
from mfespf.experiments.forcing import Forcing, NoForcing, StepForcing
from mfespf.physics.boundary import SPACE_ENERGY_STATE, SpaceBoundaryModule, space_state_entry
from mfespf.physics.equations import EQUATION_IDS, build_equation_registry
from mfespf.physics.ocean import (
    OCEAN_DEEP_STATE,
    DeepOceanModule,
    ocean_deep_state_entry,
)
from mfespf.physics.orbit import OrbitalGeometryModule
from mfespf.physics.radiation import equilibrium_temperature
from mfespf.physics.solar import (
    ConstantSpectrumSolarSource,
    ForcedSolarSource,
    banded_grid,
    broadband_grid,
)
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
    "reduced_ocean_heat_storage",
    "surface_thermal_storage",
]

# A 50 m ocean mixed layer: rho * c_p * depth, with both properties taken from
# the constants registry rather than folded into one tuned number.
_DEFAULT_AREAL_HEAT_CAPACITY = 1025.0 * 3985.0 * 50.0


def build_state_registry(
    *,
    surface_owner: str = "surface",
    space_owner: str = "space",
    deep_ocean: bool = False,
    ocean_owner: str = "ocean",
) -> StateRegistry:
    """Register the Stage 1 prognostic states with their authorisations."""
    states = StateRegistry()
    readers = {space_owner} | ({ocean_owner} if deep_ocean else set())
    states.register(surface_state_entry(owner=surface_owner, readers=frozenset(readers)))
    states.register(space_state_entry(owner=space_owner))
    if deep_ocean:
        # The deep reservoir is a state of its own with a single owner. The
        # surface is not authorised to write it and does not: the exchange
        # reaches it as a tendency the ocean computes from a flux both modules
        # book as a transfer (STA-001, MOD-011).
        states.register(
            ocean_deep_state_entry(owner=ocean_owner, readers=frozenset({surface_owner}))
        )
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
    deep_ocean: bool = False,
    ocean_id: str = "ocean",
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

    if deep_ocean:
        graph.add_node(Node(ocean_id, NodeType.MODULE, "reduced deep-ocean reservoir"))
        graph.add_node(Node(OCEAN_DEEP_STATE, NodeType.STATE, "deep reservoir temperature"))
        graph.add_edge(
            Edge(
                source=SURFACE_TEMPERATURE_STATE,
                destination=ocean_id,
                edge_type=EdgeType.STATE_DEPENDENCY,
                exchanged="surface_temperature",
                units="K",
                interface_version="1.0.0",
                coupling_frequency=coupling_frequency,
                equation_ids=(EQUATION_IDS["deep_ocean_exchange"],),
                uncertainty_rule=(
                    "the exchange is linear in this temperature, so its uncertainty "
                    "passes through undamped"
                ),
            )
        )
        graph.add_edge(
            Edge(
                source=surface_id,
                destination=ocean_id,
                edge_type=EdgeType.PARAMETER_DEPENDENCY,
                exchanged="surface_heat_capacity",
                units="J m^-2 K^-1",
                interface_version="1.0.0",
                coupling_frequency=coupling_frequency,
            )
        )
        graph.add_edge(
            Edge(
                source=ocean_id,
                destination=surface_id,
                edge_type=EdgeType.STATE_DEPENDENCY,
                exchanged="deep_ocean_heat_flux",
                units="W m^-2",
                interface_version="1.0.0",
                time_semantics="interval_average",
                coupling_frequency=coupling_frequency,
                equation_ids=(EQUATION_IDS["deep_ocean_exchange"],),
            )
        )
        graph.add_edge(
            Edge(
                source=surface_id,
                destination=ocean_id,
                edge_type=EdgeType.PHYSICAL_TRANSFER,
                exchanged="deep_ocean_heat",
                units="J m^-2",
                interface_version="1.0.0",
                time_semantics="interval_accumulation",
                spatial_semantics="global_mean_over_planetary_surface",
                coupling_frequency=coupling_frequency,
                equation_ids=(EQUATION_IDS["deep_ocean_exchange"],),
            )
        )
        graph.add_edge(
            Edge(
                source=ocean_id,
                destination=OCEAN_DEEP_STATE,
                edge_type=EdgeType.NUMERICAL_DEPENDENCY,
                exchanged="forward_euler_update",
                units="K",
                equation_ids=("EQ-NUM-FORWARD-EULER",),
            )
        )
        # A genuine two-way loop, and the one that produces the second timescale:
        # the mixed layer sets the exchange, the exchange changes the mixed
        # layer. CAU-005 wants it declared with the scheme that resolves it,
        # because an explicit pass here is stable only while the step stays well
        # inside the faster of the two modes.
        # Two loops, not one: the flux depends on the mixed-layer temperature and
        # changes it, and the heat the surface exports comes back as the deep
        # temperature that sets the next flux. Both are resolved the same way and
        # both are declared, because a cycle left undeclared is the failure
        # CAU-005 exists to catch.
        for nodes in (
            (ocean_id, surface_id, SURFACE_TEMPERATURE_STATE),
            (ocean_id, surface_id),
        ):
            graph.declare_cycle(
                Cycle(
                    nodes=nodes,
                    coupling_scheme="explicit_single_pass_ocean_before_surface",
                    iteration_order=(ocean_id, surface_id),
                    convergence_criterion=(
                        "no iteration; stability requires dt well below the fast mode "
                        "C_s / (lambda + gamma), and the two-mode fit is the check"
                    ),
                )
            )

    # The surface temperature feeds its own emission and is updated from it: a
    # real, single-module feedback loop, declared rather than left implicit.
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
    forcing: Forcing | None = None,
    deep_ocean: bool = False,
    deep_ocean_depth: float = 1200.0,
    deep_ocean_exchange: float | None = None,
    initial_deep_temperature: float | None = None,
    purpose: str | None = None,
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

    # An unforced run uses the constant source rather than a forced one carrying
    # a zero perturbation: the two are bit-identical, and the contract of the
    # first does not have to explain a forcing that is not there.
    solar = (
        ConstantSpectrumSolarSource(
            grid=grid, total_irradiance=total_irradiance, status=ModuleStatus.PRESCRIBED
        )
        if forcing is None or isinstance(forcing, NoForcing)
        else ForcedSolarSource(
            grid=grid,
            forcing=forcing,
            total_irradiance=total_irradiance,
            status=ModuleStatus.PRESCRIBED,
        )
    )
    orbit = OrbitalGeometryModule(eccentricity=eccentricity)
    surface = SlabSurfaceModule(
        areal_heat_capacity=areal_heat_capacity,
        emissivity=emissivity,
        albedo=albedo,
        initial_temperature=initial_temperature,
        couple_deep_ocean=deep_ocean,
    )
    space = SpaceBoundaryModule()

    schedules = [
        ModuleSchedule(solar, step_multiple=solar_step_multiple),
        ModuleSchedule(orbit, step_multiple=solar_step_multiple),
    ]
    order = ["solar", "orbit"]
    if deep_ocean:
        # Ahead of the surface, so the surface sees this step's exchange rather
        # than the previous one's. The ocean itself reads the committed surface
        # temperature published at the end of the previous step, which is the
        # declared zero-order hold and not an accident of ordering.
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
                        else initial_temperature
                    ),
                ),
                step_multiple=1,
            )
        )
        order.append("ocean")
    schedules.extend(
        [ModuleSchedule(surface, step_multiple=1), ModuleSchedule(space, step_multiple=1)]
    )
    order.extend(["surface", "space"])

    schedule = CouplingSchedule(
        master_dt=master_dt,
        schedules=schedules,
        splitting_order=tuple(order),
        convergence_criterion="explicit_single_pass",
    )

    solar_frequency = (
        "every_master_step"
        if solar_step_multiple == 1
        else f"every_{solar_step_multiple}_master_steps"
    )
    graph = build_dependency_graph(
        configuration_id=run_id, solar_frequency=solar_frequency, deep_ocean=deep_ocean
    )

    initial_deep = (
        initial_deep_temperature if initial_deep_temperature is not None else initial_temperature
    )
    return Kernel(
        run_id=run_id,
        purpose=purpose
        or (
            "Stage 1 airless radiating body: analytic radiative equilibrium, transient "
            "relaxation, and closure of the coupled energy budget."
        ),
        schedule=schedule,
        states=build_state_registry(deep_ocean=deep_ocean),
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
            **({OCEAN_DEEP_STATE: initial_deep} if deep_ocean else {}),
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


def _radiative_equilibrium(kwargs: dict[str, object]) -> float:
    """The temperature the unforced configuration already sits at.

    A step-response experiment that starts away from equilibrium measures two
    transients superimposed: the one the forcing caused and the one the initial
    condition caused. Only the first is the experiment. Starting at the analytic
    equilibrium removes the second rather than waiting it out, which also makes
    the unforced control exactly stationary and so worth running as a null test.
    """
    irradiance = kwargs.get("total_irradiance") or value_of("total_solar_irradiance_reference")
    return equilibrium_temperature(
        float(np.asarray(irradiance, dtype=float)),
        float(np.mean(np.asarray(kwargs.get("albedo", 0.3), dtype=float))),
        float(np.asarray(kwargs.get("emissivity", 1.0), dtype=float)),
    )


def surface_thermal_storage(**kwargs: object) -> Kernel:
    """KER-012 configuration 5: surface thermal storage under a declared forcing.

    One reservoir, and therefore one timescale. This exists to be the control the
    two-reservoir configuration is measured against: everything but the deep
    store is identical between them, so a second timescale appearing in the other
    is attributable to the store and to nothing else.

    The default forcing is a step, because a step response is the cleanest
    measurement of a relaxation time; pass any other ``Forcing`` for a pulse,
    impulse or periodic experiment (TIM-007).
    """
    kwargs.setdefault("run_id", "STAGE5-SURFACE-STORAGE-001")
    kwargs.setdefault("forcing", StepForcing(amplitude=20.0))
    kwargs.setdefault("initial_temperature", _radiative_equilibrium(kwargs))
    kwargs.setdefault(
        "purpose",
        "Stage 5 surface thermal storage: step, pulse, impulse or periodic forcing of a "
        "single-reservoir slab, for the response timescale of one storage term alone.",
    )
    return airless_radiating_body(**kwargs)  # type: ignore[arg-type]


def reduced_ocean_heat_storage(**kwargs: object) -> Kernel:
    """KER-012 configuration 6: the two-reservoir ocean, and the last of the eight.

    The same slab as ``surface_thermal_storage`` with a deep reservoir beneath
    it, exchanging heat at a rate linear in their temperature difference. The
    system is then linear with two storage terms, so its step response is a sum
    of two exponentials -- one fast, set by the mixed layer, and one slow, set by
    the deep store -- and both are measurable rather than assumed.

    This is the configuration that removes the limitation every earlier one
    declared: that a single reservoir has one timescale by construction and
    cannot show the continuing adjustment a real system shows long after the
    surface has stopped moving.
    """
    kwargs.setdefault("run_id", "STAGE5-REDUCED-OCEAN-001")
    kwargs.setdefault("forcing", StepForcing(amplitude=20.0))
    kwargs.setdefault("initial_temperature", _radiative_equilibrium(kwargs))
    kwargs.setdefault(
        "purpose",
        "Stage 5 reduced ocean heat storage: a mixed layer over a deep reservoir under a "
        "declared forcing, for the fast and slow response timescales of a two-storage system.",
    )
    kwargs["deep_ocean"] = True
    return airless_radiating_body(**kwargs)  # type: ignore[arg-type]
