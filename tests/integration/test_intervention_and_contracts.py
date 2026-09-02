"""Module contract completeness (MOD-*) and controlled intervention (CAU-007)."""

from __future__ import annotations

import pytest

from mfespf.contracts.module import ModuleContract, ModuleStatus
from mfespf.core.errors import ConfigurationError
from mfespf.core.versioning import Version
from mfespf.run.configurations import airless_radiating_body

pytestmark = pytest.mark.stage(0)

C_SMALL = 1.0e6


@pytest.mark.req("MOD-001", "MOD-002", "MOD-003", "MOD-004", "SYS-003")
def test_every_kernel_module_publishes_a_complete_contract():
    kernel = airless_radiating_body()
    for schedule in kernel.schedule.ordered():
        contract = schedule.module.contract
        contract.validate()
        assert contract.module_id and contract.module_name
        assert contract.governing_equation_ids
        assert contract.physical_basis and contract.validity_domain
        assert contract.correctness_tests
        assert contract.cost is not None
        assert schedule.module.status in set(ModuleStatus)


@pytest.mark.req("MOD-002", "EQR-005")
def test_every_cited_equation_exists_in_the_registry():
    """A module cannot claim to use an equation nobody registered."""
    kernel = airless_radiating_body()
    registered = set(kernel.equations.ids())
    for schedule in kernel.schedule.ordered():
        for equation_id in schedule.module.contract.governing_equation_ids:
            assert equation_id in registered, (
                f"{schedule.module_id} cites unregistered {equation_id}"
            )


@pytest.mark.req("MOD-013", "CMP-001")
def test_a_contract_without_a_cost_declaration_is_refused():
    with pytest.raises(ConfigurationError, match="MOD-013"):
        ModuleContract(
            module_id="m",
            module_name="m",
            module_version=Version(1, 0, 0),
            interface_version=Version(1, 0, 0),
            fidelity_level=__import__(
                "mfespf.contracts.module", fromlist=["FidelityLevel"]
            ).FidelityLevel.REDUCED,
            physical_basis="basis",
            governing_equation_ids=("EQ-X",),
            assumptions=(),
            validity_domain="domain",
            correctness_tests=("test_x",),
        ).validate()


@pytest.mark.req("MOD-014", "VAL-001")
def test_a_contract_without_a_null_or_limiting_case_test_is_refused():
    from mfespf.contracts.cost import CostDeclaration
    from mfespf.contracts.module import FidelityLevel

    with pytest.raises(ConfigurationError, match="MOD-014"):
        ModuleContract(
            module_id="m",
            module_name="m",
            module_version=Version(1, 0, 0),
            interface_version=Version(1, 0, 0),
            fidelity_level=FidelityLevel.REDUCED,
            physical_basis="basis",
            governing_equation_ids=("EQ-X",),
            assumptions=(),
            validity_domain="domain",
            correctness_tests=(),
            cost=CostDeclaration(module_id="m", time_complexity="O(1)", memory_complexity="O(1)"),
        ).validate()


@pytest.mark.req("MOD-009", "UNC-001", "UNC-004")
def test_every_module_reports_its_unresolved_uncertainty():
    """UNC-004: an applicable contribution that is not quantified must be marked."""
    kernel = airless_radiating_body(areal_heat_capacity=C_SMALL, solar_step_multiple=1)
    coupler, _ = kernel.prepare(n_steps=2)
    coupler.step()

    for module_id, report in coupler.uncertainty.items():
        assert report.module_id == module_id
        assert report.unresolved or report.quantified, module_id
    surface = coupler.uncertainty["surface"]
    assert any("structural" in u for u in surface.unresolved)
    assert "UNRESOLVED" in surface.describe()


@pytest.mark.req("CAU-007", "VAL-004")
def test_disabling_a_declared_edge_withholds_the_field_rather_than_zeroing_it():
    """An intervention removes an influence; it does not fabricate a zero one.

    DAT-007 forbids substituting zero for a missing physical input, so a module
    that genuinely requires the withheld field must fail loudly. Otherwise the
    control experiment would silently be a different experiment.
    """
    kernel = airless_radiating_body(areal_heat_capacity=C_SMALL, solar_step_multiple=1)
    kernel.graph.disable_edge(
        "orbit", "surface", "toa_incident_irradiance", reason="isolate the emission term"
    )
    coupler, _ = kernel.prepare(n_steps=2)

    with pytest.raises(ConfigurationError, match="did not receive required inputs"):
        coupler.step()


@pytest.mark.req("CAU-001", "AUD-002")
def test_the_run_manifest_carries_the_dependency_graph():
    kernel = airless_radiating_body()
    manifest = kernel.build_manifest()
    graph = manifest.dependency_graph
    assert graph["nodes"] and graph["edges"]
    assert "not evidence of real-world causation" in graph["disclaimer"]
    edge_pairs = {(e["source"], e["destination"]) for e in graph["edges"]}
    assert ("orbit", "surface") in edge_pairs
    assert ("surface", "space") in edge_pairs


@pytest.mark.req("SYS-003", "AUD-002")
def test_the_manifest_declares_a_status_for_every_known_module():
    kernel = airless_radiating_body()
    manifest = kernel.build_manifest()
    assert set(manifest.module_status) == {"solar", "orbit", "surface", "space"}
    assert manifest.module_status["solar"] == "prescribed"
    assert manifest.module_status["surface"] == "active"


@pytest.mark.req("AUD-002")
def test_module_parameters_reach_the_configuration_hash():
    """A parameter that does not enter the hash makes two runs indistinguishable."""
    baseline = airless_radiating_body(run_id="A").build_manifest()
    changed = airless_radiating_body(run_id="B", emissivity=0.9).build_manifest()
    assert baseline.configuration_hash() != changed.configuration_hash()
    assert baseline.parameters["surface"]["emissivity"] == 1.0
    assert changed.parameters["surface"]["emissivity"] == 0.9
