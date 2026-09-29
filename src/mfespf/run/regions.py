"""Several columns in one run: regions.

Why this exists
---------------

A single horizontally uniform column cannot hold air that sinks, because its
mass-weighted vertical motion is zero. Earth's free troposphere is dry because
most of it is sinking air last saturated somewhere colder, and every failure of
the observed-cloud test ran back to a free troposphere too wet for that
(``docs/protocols/observed-cloud-field.md``). Representing it needs a second
column, and a second column needs the framework to hold two atmospheres.

How
---

Regions belong to the coupler, not to the modules. ``regionalize`` takes an
ordinary kernel and files it under a region: every module id becomes
``region/id`` -- including the ids modules hold for one another, so transfers
and budgets stay paired within the region -- and every state becomes
``region:state`` in the registry, the graph and the initial state. The coupler
then routes a regional module's fields and tendencies under its region, so each
module still publishes and reads plain names and needs no change at all.
``combine_regions`` merges regional kernels, and any global modules, into one
run.

What makes this checkable is that two unconnected regions built from one
configuration must each reproduce that configuration run alone, to the bit.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import replace
from typing import Any

from mfespf.core.errors import ConfigurationError
from mfespf.coupling.coupler import CouplingSchedule, ModuleSchedule, qualify
from mfespf.coupling.graph import Cycle, DependencyGraph, NodeType
from mfespf.registries.states import StateRegistry
from mfespf.run.driver import Kernel

__all__ = ["combine_regions", "regionalize"]

REGION_MODULE_SEPARATOR = "/"


def _renamed(old: str, mapping: Mapping[str, str], region: str) -> str:
    return mapping.get(old, qualify(region, old))


def regionalize(kernel: Kernel, region: str) -> Kernel:
    """File a single-column kernel under ``region``, modules and states alike."""
    if not region or ":" in region or REGION_MODULE_SEPARATOR in region:
        raise ConfigurationError(f"region name {region!r} must be non-empty and plain")
    if kernel.coupler is not None:
        raise ConfigurationError("regionalize a kernel before it runs, not after")
    mapping = {
        s.module_id: f"{region}{REGION_MODULE_SEPARATOR}{s.module_id}"
        for s in kernel.schedule.schedules
    }

    schedules = []
    for s in kernel.schedule.schedules:
        module = s.module
        new_id = mapping[s.module_id]
        cost = module.contract.cost
        module.contract = replace(
            module.contract,
            module_id=new_id,
            cost=replace(cost, module_id=new_id) if cost is not None else None,
        )
        # The ids modules keep for one another: a surface's atmosphere, a space
        # boundary's sender. Renamed in step, or a transfer would pair a
        # region's atmosphere with a module that no longer exists.
        for name, value in list(vars(module).items()):
            if name.endswith("_module_id") and isinstance(value, str) and value in mapping:
                setattr(module, name, mapping[value])
        schedules.append(replace(s, region=region))

    states = StateRegistry()
    for entry in kernel.states.current():
        states.register(
            replace(
                entry,
                record_id=qualify(region, entry.record_id),
                owner_module_id=mapping.get(entry.owner_module_id, entry.owner_module_id),
                authorized_readers=frozenset(mapping.get(m, m) for m in entry.authorized_readers),
                authorized_tendency_providers=frozenset(
                    mapping.get(m, m) for m in entry.authorized_tendency_providers
                ),
                authorized_boundary_providers=frozenset(
                    mapping.get(m, m) for m in entry.authorized_boundary_providers
                ),
            )
        )

    graph = DependencyGraph(f"{region}:{kernel.graph.configuration_id}")
    for node in kernel.graph.nodes:
        new_id = (
            mapping.get(node.node_id, node.node_id)
            if node.node_type is NodeType.MODULE
            else _renamed(node.node_id, mapping, region)
        )
        graph.add_node(replace(node, node_id=new_id))
    for edge in kernel.graph.edges:
        graph.add_edge(
            replace(
                edge,
                source=_renamed(edge.source, mapping, region),
                destination=_renamed(edge.destination, mapping, region),
            )
        )
    for cycle in kernel.graph.declared_cycles:
        graph.declare_cycle(
            Cycle(
                nodes=tuple(_renamed(n, mapping, region) for n in cycle.nodes),
                coupling_scheme=cycle.coupling_scheme,
                iteration_order=tuple(_renamed(n, mapping, region) for n in cycle.iteration_order),
                convergence_criterion=cycle.convergence_criterion,
                residual=cycle.residual,
            )
        )

    return replace(
        kernel,
        run_id=f"{kernel.run_id}@{region}",
        schedule=replace(
            kernel.schedule,
            schedules=schedules,
            splitting_order=tuple(mapping.get(m, m) for m in kernel.schedule.splitting_order),
        ),
        states=states,
        graph=graph,
        initial_state={qualify(region, k): v for k, v in kernel.initial_state.items()},
        boundary_state={qualify(region, k): v for k, v in kernel.boundary_state.items()},
        manifest=None,
    )


def combine_regions(
    regions: Sequence[Kernel],
    *,
    run_id: str,
    purpose: str,
    global_modules: Sequence[ModuleSchedule] = (),
    middle_modules: Sequence[ModuleSchedule] = (),
    split_after: str = "radiation",
    budget: Any = None,
) -> Kernel:
    """One run from regional kernels, plus any modules that span them.

    Each region keeps its own splitting order. Middle modules run after every
    region's ``split_after`` module and before the rest of any region: a
    circulation needs this step's radiative cooling from both columns, and the
    columns' water and convection must compose with its transport. Trailing
    global modules run last and read what every region published.
    """
    if not regions:
        raise ConfigurationError("combine at least one region")
    first = regions[0]
    if any(k.schedule.master_dt != first.schedule.master_dt for k in regions):
        raise ConfigurationError("regions must share one master timestep")
    if any(k.operating_mode is not first.operating_mode for k in regions):
        raise ConfigurationError("regions must share one operating mode")

    def split(k: Kernel) -> tuple[list[ModuleSchedule], list[ModuleSchedule]]:
        ordered = list(k.schedule.schedules)
        ids = [s.module_id.split(REGION_MODULE_SEPARATOR, 1)[-1] for s in ordered]
        if not middle_modules or split_after not in ids:
            return ordered, []
        cut = ids.index(split_after) + 1
        return ordered[:cut], ordered[cut:]

    parts = [split(k) for k in regions]
    schedules = (
        [s for head, _ in parts for s in head]
        + list(middle_modules)
        + [s for _, tail in parts for s in tail]
        + list(global_modules)
    )
    states = StateRegistry()
    graph = DependencyGraph(run_id)
    for k in regions:
        for entry in k.states.current():
            states.register(entry)
        for node in k.graph.nodes:
            graph.add_node(node)
        for edge in k.graph.edges:
            graph.add_edge(edge)
        for cycle in k.graph.declared_cycles:
            graph.declare_cycle(cycle)

    return replace(
        first,
        run_id=run_id,
        purpose=purpose,
        schedule=CouplingSchedule(
            master_dt=first.schedule.master_dt,
            schedules=schedules,
            splitting_order=tuple(s.module_id for s in schedules),
            master_time_units=first.schedule.master_time_units,
        ),
        states=states,
        graph=graph,
        # Several columns hold several columns' trajectories: the declared
        # memory is the regions' declarations summed, not one column's reused.
        budget=budget
        if budget is not None
        else replace(
            first.budget,
            max_memory_bytes=sum(k.budget.max_memory_bytes for k in regions),
        ),
        initial_state={key: v for k in regions for key, v in k.initial_state.items()},
        boundary_state={key: v for k in regions for key, v in k.boundary_state.items()},
        conservation_tolerance=max(k.conservation_tolerance for k in regions),
        manifest=None,
        coupler=None,
    )
