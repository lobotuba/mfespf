"""Implemented dependency graph.

Amendment A2 and requirements CAU-001 through CAU-008, plus KER-004.

The graph records the information and physical-influence paths that this
configuration actually implements.  CAU-008 is not a footnote: this is a map of
the software, and it is reported as such.  It cannot establish real-world
causation, and the ``describe`` output says so wherever the graph is printed.

Its practical value is intervention (A2.5, CAU-007): a declared edge can be
disabled, clamped or replaced for a controlled experiment, and because every
edge is explicit, disabling one cannot silently change another.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from enum import StrEnum

from mfespf.core.errors import ConfigurationError, CouplingError

__all__ = ["Cycle", "DependencyGraph", "Edge", "EdgeType", "Node", "NodeType"]


class NodeType(StrEnum):
    """CAU-002."""

    MODULE = "module"
    STATE = "state"
    PARAMETER = "parameter"
    DATASET = "dataset"
    DIAGNOSTIC = "diagnostic"
    BOUNDARY = "boundary"


class EdgeType(StrEnum):
    """A2.2 / CAU-003."""

    PHYSICAL_TRANSFER = "physical_transfer"
    STATE_DEPENDENCY = "state_dependency"
    PARAMETER_DEPENDENCY = "parameter_dependency"
    BOUNDARY_DEPENDENCY = "boundary_dependency"
    INITIAL_CONDITION_DEPENDENCY = "initial_condition_dependency"
    DIAGNOSTIC_DEPENDENCY = "diagnostic_dependency"
    CALIBRATION_DEPENDENCY = "calibration_dependency"
    NUMERICAL_DEPENDENCY = "numerical_dependency"
    UNCERTAINTY_DEPENDENCY = "uncertainty_dependency"

    @property
    def is_physical(self) -> bool:
        return self is EdgeType.PHYSICAL_TRANSFER


@dataclass(frozen=True, slots=True)
class Node:
    node_id: str
    node_type: NodeType
    label: str = ""


@dataclass(frozen=True)
class Edge:
    """A2.3 / CAU-004.  Every edge carries the semantics of what crosses it."""

    source: str
    destination: str
    edge_type: EdgeType
    exchanged: str
    interface_version: str = ""
    units: str = ""
    time_semantics: str = ""
    spatial_semantics: str = ""
    spectral_semantics: str = ""
    coupling_frequency: str = ""
    equation_ids: tuple[str, ...] = ()
    uncertainty_rule: str = "unresolved"
    enabled: bool = True
    intervention: str = ""

    @property
    def key(self) -> tuple[str, str, str, str]:
        return (self.source, self.destination, self.edge_type.value, self.exchanged)

    def describe(self) -> str:
        state = "" if self.enabled else f" [DISABLED: {self.intervention or 'intervention'}]"
        units = f" [{self.units}]" if self.units else ""
        return (
            f"{self.source} --{self.edge_type.value}:{self.exchanged}{units}--> "
            f"{self.destination}{state}"
        )


@dataclass
class Cycle:
    """A feedback loop, with the resolution scheme A2.4 and CAU-005 require."""

    nodes: tuple[str, ...]
    coupling_scheme: str = "undeclared"
    iteration_order: tuple[str, ...] = ()
    convergence_criterion: str = "undeclared"
    residual: float | None = None

    @property
    def is_declared(self) -> bool:
        return self.coupling_scheme != "undeclared" and self.convergence_criterion != "undeclared"

    def describe(self) -> str:
        loop = " -> ".join([*self.nodes, self.nodes[0]])
        detail = f"scheme={self.coupling_scheme}, convergence={self.convergence_criterion}"
        residual = f", residual={self.residual:.3e}" if self.residual is not None else ""
        return f"{loop} ({detail}{residual})"


class DependencyGraph:
    """A directed, typed graph of implemented dependencies for one configuration."""

    def __init__(self, configuration_id: str = "") -> None:
        self.configuration_id = configuration_id
        self._nodes: dict[str, Node] = {}
        self._edges: dict[tuple[str, str, str, str], Edge] = {}
        self._declared_cycles: list[Cycle] = []

    # --- construction -----------------------------------------------------
    def add_node(self, node: Node) -> Node:
        existing = self._nodes.get(node.node_id)
        if existing is not None and existing.node_type is not node.node_type:
            raise ConfigurationError(
                f"node {node.node_id!r} is already registered as {existing.node_type.value}, "
                f"cannot re-add as {node.node_type.value}"
            )
        self._nodes[node.node_id] = node
        return node

    def add_edge(self, edge: Edge) -> Edge:
        for endpoint in (edge.source, edge.destination):
            if endpoint not in self._nodes:
                raise ConfigurationError(
                    f"edge {edge.describe()} references unknown node {endpoint!r}; every "
                    f"dependency must name declared nodes (CAU-002, CAU-006)"
                )
        self._edges[edge.key] = edge
        return edge

    def declare_cycle(self, cycle: Cycle) -> Cycle:
        self._declared_cycles.append(cycle)
        return cycle

    # --- intervention (CAU-007, A2.5) -------------------------------------
    def disable_edge(self, source: str, destination: str, exchanged: str, *, reason: str) -> Edge:
        """Disable exactly one declared edge for a controlled experiment.

        Matching is on the triple, so an intervention names the specific
        influence it removes.  Disabling an edge that does not exist is an error
        rather than a no-op: a control experiment that silently did nothing is
        worse than one that failed to start.
        """
        matches = [
            e
            for e in self._edges.values()
            if e.source == source and e.destination == destination and e.exchanged == exchanged
        ]
        if not matches:
            raise CouplingError(
                f"no declared edge {source} -> {destination} carrying {exchanged!r} to disable"
            )
        if len(matches) > 1:
            raise CouplingError(
                f"{len(matches)} edges match {source} -> {destination} carrying {exchanged!r}; "
                f"disambiguate by edge type"
            )
        edge = matches[0]
        disabled = Edge(**{**edge.__dict__, "enabled": False, "intervention": reason})
        self._edges[edge.key] = disabled
        return disabled

    def is_enabled(self, source: str, destination: str, exchanged: str) -> bool:
        return all(
            e.enabled
            for e in self._edges.values()
            if e.source == source and e.destination == destination and e.exchanged == exchanged
        )

    # --- queries ----------------------------------------------------------
    @property
    def nodes(self) -> list[Node]:
        return [self._nodes[k] for k in sorted(self._nodes)]

    @property
    def edges(self) -> list[Edge]:
        return [self._edges[k] for k in sorted(self._edges)]

    def active_edges(self) -> list[Edge]:
        return [e for e in self.edges if e.enabled]

    def successors(self, node_id: str) -> list[str]:
        return sorted({e.destination for e in self.active_edges() if e.source == node_id})

    def predecessors(self, node_id: str) -> list[str]:
        return sorted({e.source for e in self.active_edges() if e.destination == node_id})

    def find_cycles(self) -> list[tuple[str, ...]]:
        """Every elementary cycle among active edges.

        Depth-first with an explicit stack; the graphs here are small (tens of
        nodes), so exhaustive enumeration is affordable and exact.
        """
        adjacency: dict[str, list[str]] = defaultdict(list)
        for e in self.active_edges():
            adjacency[e.source].append(e.destination)

        cycles: set[tuple[str, ...]] = set()

        def walk(start: str, node: str, path: list[str], seen: set[str]) -> None:
            for nxt in adjacency.get(node, []):
                if nxt == start:
                    canonical = _canonical_cycle(path)
                    cycles.add(canonical)
                elif nxt not in seen and nxt > start:
                    # Only explore nodes ordered after the start to visit each
                    # elementary cycle from a single representative.
                    walk(start, nxt, [*path, nxt], seen | {nxt})

        for node in sorted(adjacency):
            walk(node, node, [node], {node})
        return sorted(cycles)

    def undeclared_cycles(self) -> list[tuple[str, ...]]:
        """Cycles present in the graph that no ``Cycle`` record explains (CAU-005)."""
        declared = {_canonical_cycle(list(c.nodes)) for c in self._declared_cycles if c.is_declared}
        return [c for c in self.find_cycles() if c not in declared]

    def validate(self) -> None:
        """Fail on a feedback loop with no declared coupling scheme (A2.4, CAU-005)."""
        undeclared = self.undeclared_cycles()
        if undeclared:
            listed = "\n  - ".join(" -> ".join(c) for c in undeclared)
            raise CouplingError(
                "the configuration contains feedback cycles with no declared coupling scheme, "
                "iteration order or convergence criterion (CAU-005):\n  - " + listed
            )

    def to_dict(self) -> dict[str, object]:
        """Serialisable form for the run manifest."""
        return {
            "configuration_id": self.configuration_id,
            "disclaimer": (
                "Implemented dependency graph. It documents information and physical influence "
                "paths in this configuration and is not evidence of real-world causation "
                "(CAU-008)."
            ),
            "nodes": [
                {"id": n.node_id, "type": n.node_type.value, "label": n.label} for n in self.nodes
            ],
            "edges": [
                {
                    "source": e.source,
                    "destination": e.destination,
                    "type": e.edge_type.value,
                    "exchanged": e.exchanged,
                    "units": e.units,
                    "interface_version": e.interface_version,
                    "time_semantics": e.time_semantics,
                    "spatial_semantics": e.spatial_semantics,
                    "spectral_semantics": e.spectral_semantics,
                    "coupling_frequency": e.coupling_frequency,
                    "equation_ids": list(e.equation_ids),
                    "uncertainty_rule": e.uncertainty_rule,
                    "enabled": e.enabled,
                    "intervention": e.intervention,
                }
                for e in self.edges
            ],
            "cycles": [
                {
                    "nodes": list(c.nodes),
                    "coupling_scheme": c.coupling_scheme,
                    "iteration_order": list(c.iteration_order),
                    "convergence_criterion": c.convergence_criterion,
                    "residual": c.residual,
                }
                for c in self._declared_cycles
            ],
        }

    def describe(self) -> str:
        lines = [
            f"implemented dependency graph ({self.configuration_id or 'unnamed'})",
            "  NOTE: implemented dependencies only; not evidence of causation (CAU-008)",
            f"  {len(self._nodes)} nodes, {len(self._edges)} edges "
            f"({len(self.active_edges())} enabled)",
        ]
        lines.extend(f"  {e.describe()}" for e in self.edges)
        for c in self._declared_cycles:
            lines.append(f"  cycle: {c.describe()}")
        return "\n".join(lines)


def _canonical_cycle(path: list[str]) -> tuple[str, ...]:
    """Rotate a cycle so it starts at its lexicographically smallest node."""
    if not path:
        return ()
    i = min(range(len(path)), key=lambda k: path[k])
    return tuple(path[i:] + path[:i])
