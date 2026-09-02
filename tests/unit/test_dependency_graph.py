"""Stage 0: the implemented dependency graph and controlled interventions."""

from __future__ import annotations

import pytest

from mfespf.core.errors import ConfigurationError, CouplingError
from mfespf.coupling.graph import (
    Cycle,
    DependencyGraph,
    Edge,
    EdgeType,
    Node,
    NodeType,
)

pytestmark = pytest.mark.stage(0)


def _graph() -> DependencyGraph:
    graph = DependencyGraph("test")
    for name in ("a", "b", "c"):
        graph.add_node(Node(name, NodeType.MODULE))
    return graph


def _edge(source: str, destination: str, exchanged: str = "energy") -> Edge:
    return Edge(
        source=source,
        destination=destination,
        edge_type=EdgeType.PHYSICAL_TRANSFER,
        exchanged=exchanged,
        units="W m^-2",
        interface_version="1.0.0",
    )


@pytest.mark.req("CAU-001", "CAU-002", "KER-004")
def test_graph_records_typed_nodes_and_edges():
    graph = _graph()
    graph.add_edge(_edge("a", "b"))
    assert len(graph.nodes) == 3
    assert graph.successors("a") == ["b"]
    assert graph.predecessors("b") == ["a"]


@pytest.mark.req("CAU-004")
def test_every_edge_carries_its_quantity_and_interface_version():
    edge = _edge("a", "b")
    assert edge.units == "W m^-2"
    assert edge.interface_version == "1.0.0"
    assert "energy" in edge.describe()


@pytest.mark.req("CAU-006")
def test_an_edge_to_an_undeclared_node_is_refused():
    graph = _graph()
    with pytest.raises(ConfigurationError, match="unknown node"):
        graph.add_edge(_edge("a", "nowhere"))


@pytest.mark.req("CAU-005")
def test_an_undeclared_feedback_cycle_fails_validation():
    graph = _graph()
    graph.add_edge(_edge("a", "b"))
    graph.add_edge(_edge("b", "a", "back_reaction"))
    assert graph.find_cycles() == [("a", "b")]
    with pytest.raises(CouplingError, match="no declared coupling scheme"):
        graph.validate()


@pytest.mark.req("CAU-005")
def test_a_declared_cycle_passes_validation():
    graph = _graph()
    graph.add_edge(_edge("a", "b"))
    graph.add_edge(_edge("b", "a", "back_reaction"))
    graph.declare_cycle(
        Cycle(
            nodes=("a", "b"),
            coupling_scheme="explicit_sequential",
            iteration_order=("a", "b"),
            convergence_criterion="single pass; residual reported",
        )
    )
    graph.validate()
    assert graph.undeclared_cycles() == []


@pytest.mark.req("CAU-005")
def test_longer_cycles_are_found():
    graph = _graph()
    graph.add_edge(_edge("a", "b"))
    graph.add_edge(_edge("b", "c"))
    graph.add_edge(_edge("c", "a"))
    assert graph.find_cycles() == [("a", "b", "c")]


@pytest.mark.req("CAU-007")
def test_disabling_one_edge_leaves_the_others_untouched():
    graph = _graph()
    graph.add_edge(_edge("a", "b", "energy"))
    graph.add_edge(_edge("a", "b", "moisture"))
    graph.add_edge(_edge("b", "c", "energy"))

    graph.disable_edge("a", "b", "moisture", reason="isolate the thermal pathway")

    assert not graph.is_enabled("a", "b", "moisture")
    assert graph.is_enabled("a", "b", "energy")
    assert graph.is_enabled("b", "c", "energy")
    assert len(graph.active_edges()) == 2


@pytest.mark.req("CAU-007")
def test_disabling_an_edge_that_does_not_exist_is_an_error_not_a_no_op():
    """A control experiment that silently did nothing is worse than one that fails."""
    graph = _graph()
    graph.add_edge(_edge("a", "b"))
    with pytest.raises(CouplingError, match="no declared edge"):
        graph.disable_edge("a", "c", "energy", reason="typo in the intervention")


@pytest.mark.req("CAU-008")
def test_the_graph_disclaims_causal_interpretation_wherever_it_is_reported():
    graph = _graph()
    graph.add_edge(_edge("a", "b"))
    assert "not evidence of real-world causation" in graph.to_dict()["disclaimer"]
    assert "not evidence of causation" in graph.describe()
