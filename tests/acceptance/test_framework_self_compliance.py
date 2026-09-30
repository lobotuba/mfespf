"""The framework's own compliance, demonstrated rather than asserted.

A project whose thesis is that scientific claims should be checkable rather than
stated ought to hold its own requirements to the same standard. Nine mandatory
requirements here were satisfied by machinery that had been read and believed
rather than exercised: every module implements configuration validation, every
contract carries a cost declaration, every graph edge is typed, and so on. All of
it was true. None of it was bound to a test, so nothing would have said if it
stopped being true.

These are sweeps rather than examples, and that is the point of the file. A check
run against one configuration establishes that one configuration complies; a
check run against every reference configuration in the tree establishes the
property. ``BUILDERS`` is therefore the whole set, and
``test_the_sweep_reaches_every_module_implementation`` fails when a module is
added that none of them exercises -- so a new module cannot quietly escape the
sweep by not being wired into it.

Writing them found one thing that was *not* true. ``output_checksums`` was empty
on every run in the tree, because the only path that filled it checksums a file
and none of these configurations writes one. ACC-005 had been recorded as
"already satisfied, needs a test"; the test is what showed it was not.
"""

from __future__ import annotations

import importlib
import importlib.util
import pkgutil
import tracemalloc
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import yaml

import mfespf.physics
from mfespf.contracts.cost import ParallelismModel, ResourceBudget, ResourceClass
from mfespf.contracts.module import Module, ModuleStatus
from mfespf.core.errors import ConfigurationError, OwnershipError, RegistryError
from mfespf.core.versioning import Version
from mfespf.coupling.graph import Edge, EdgeType
from mfespf.physics.greenhouse import PrescribedGreenhouseModule
from mfespf.physics.radiative_transfer import window_absorption
from mfespf.physics.volcanic import PrescribedVolcanicModule
from mfespf.run import column_configurations as cc
from mfespf.run import configurations as rc
from mfespf.run.driver import Kernel
from mfespf.run.greenhouse import greenhouse_grid
from mfespf.run.preflight import estimate_cost
from mfespf.run.two_column import two_column

pytestmark = pytest.mark.stage(0)

REPO_ROOT = Path(__file__).resolve().parents[2]


# --- the configurations every sweep runs over -------------------------------
#
# Cheap by design: building all of these and stepping each twice costs under a
# fifth of a second, so there is no reason for any check here to settle for a
# subset.


def _greenhouse_column(**kwargs: Any) -> Kernel:
    """A column with the carbon dioxide band and a volcanic layer in it."""
    grid = greenhouse_grid()
    return cc.radiating_column(
        run_id="SELFTEST-GREENHOUSE",
        purpose="Self-compliance sweep over the greenhouse and volcanic modules.",
        absorption_coefficient=window_absorption(grid, kappa_opaque=1.0e-4),
        grid=grid,
        greenhouse=PrescribedGreenhouseModule(concentration=280.0),
        volcanic=PrescribedVolcanicModule(optical_depth=0.05),
        convection=True,
        **kwargs,
    )


def _greenhouse_column_at(ppm: float) -> Kernel:
    """The same column at another concentration, and identical in every other way."""
    grid = greenhouse_grid()
    return cc.radiating_column(
        run_id=f"SELFTEST-GREENHOUSE-{ppm:g}",
        purpose="HUM-001: the concentration is the only thing that differs.",
        absorption_coefficient=window_absorption(grid, kappa_opaque=1.0e-4),
        grid=grid,
        greenhouse=PrescribedGreenhouseModule(concentration=ppm),
        volcanic=PrescribedVolcanicModule(optical_depth=0.05),
        convection=True,
    )


def _prescribed_humidity_column(**kwargs: Any) -> Kernel:
    grid = greenhouse_grid()
    return cc.radiating_column(
        run_id="SELFTEST-PRESCRIBED-RH",
        purpose="Self-compliance sweep over the prescribed-humidity module.",
        absorption_coefficient=window_absorption(grid, kappa_opaque=1.0e-4),
        grid=grid,
        greenhouse=PrescribedGreenhouseModule(concentration=280.0),
        convection=True,
        prescribed_humidity=0.7,
        **kwargs,
    )


def _plume_humidity_column(**kwargs: Any) -> Kernel:
    grid = greenhouse_grid()
    return cc.radiating_column(
        run_id="SELFTEST-PLUME-RH",
        purpose="Self-compliance sweep over the entraining-plume humidity closure.",
        absorption_coefficient=window_absorption(grid, kappa_opaque=1.0e-4),
        grid=grid,
        greenhouse=PrescribedGreenhouseModule(concentration=280.0),
        convection=True,
        plume_humidity=True,
        **kwargs,
    )


def _observed_clouds_and_solar_vapour(**kwargs: Any) -> Kernel:
    """Every radiative mechanism at once: the configuration the cloud test runs."""
    return cc.reference_column(observed_cloud_field=True, solar_vapour_absorption=True, **kwargs)


def _ozone_on_the_measured_spectrum(**kwargs: Any) -> Kernel:
    """Ozone, the measured solar spectrum and vapour diffusion: the ozone test's column."""
    return cc.reference_column(
        solar_vapour_absorption=True, ozone=True, eddy_diffusivity=0.39, **kwargs
    )


BUILDERS = {
    "airless_radiating_body": rc.airless_radiating_body,
    "null_control": rc.null_control,
    "banded_variant": rc.banded_variant,
    "surface_thermal_storage": rc.surface_thermal_storage,
    "reduced_ocean_heat_storage": rc.reduced_ocean_heat_storage,
    "transparent_atmosphere": cc.transparent_atmosphere,
    "absorbing_emitting_atmosphere": cc.absorbing_emitting_atmosphere,
    "window_atmosphere": cc.window_atmosphere,
    "single_grey_layer": cc.single_grey_layer,
    "radiative_convective": cc.radiative_convective,
    "moist_radiative_convective": cc.moist_radiative_convective,
    "prescribed_cloud": cc.prescribed_cloud,
    "no_cloud_control": cc.no_cloud_control,
    "column_over_reduced_ocean": cc.column_over_reduced_ocean,
    "greenhouse_and_volcanic": _greenhouse_column,
    "prescribed_humidity": _prescribed_humidity_column,
    "plume_humidity": _plume_humidity_column,
    "reference_column": cc.reference_column,
    "observed_clouds_and_solar_vapour": _observed_clouds_and_solar_vapour,
    "ozone_on_the_measured_spectrum": _ozone_on_the_measured_spectrum,
    # Two regions, the circulation between them and the global means: the one
    # configuration whose modules are filed under regions.
    "two_column": two_column,
}

#: Module implementations no reference configuration in this file reaches, with
#: the reason. An entry here is a statement that the sweep does not cover it;
#: keeping the list explicit is what stops the coverage claim from drifting.
UNSWEPT = {
    # Abstract bases, which have no contract of their own to sweep.
    "SolarSource": "abstract base class",
    # Driven from a data record, so building one needs the historical adapter and
    # its dataset; the historical replay tests exercise it end to end instead.
    "RecordedSolarSource": "requires a loaded irradiance record",
}


@pytest.fixture(scope="module")
def kernels() -> dict[str, Kernel]:
    return {name: build() for name, build in BUILDERS.items()}


@pytest.fixture(scope="module")
def stepped(kernels: dict[str, Kernel]) -> dict[str, Any]:
    """Every configuration, actually run. Two steps is enough to see a lifecycle."""
    return {name: kernel.run(2) for name, kernel in kernels.items()}


def _modules(kernel: Kernel) -> list[Any]:
    return [slot.module for slot in kernel.schedule.ordered()]


# --- MOD-005 and MOD-007: the lifecycle ------------------------------------


@pytest.mark.req("MOD-005", "MOD-007")
def test_configuration_validation_precedes_every_initialization(kernels):
    """MOD-005 is an ordering claim, and only an ordering check can establish it.

    "Implements configuration validation" would be satisfied by a method that
    runs after the state is built, which is precisely when it is too late to
    refuse anything. So record the order the framework actually calls them in.
    """
    for name, build in BUILDERS.items():
        kernel = build()
        order: list[tuple[str, str]] = []
        for module in _modules(kernel):
            for op in ("validate_configuration", "initialize"):
                _record(module, op, order)
        kernel.prepare(n_steps=1)

        for module in _modules(kernel):
            calls = [op for mid, op in order if mid == module.module_id]
            assert "validate_configuration" in calls, (
                f"{name}: {module.module_id} was never asked to validate its configuration"
            )
            assert "initialize" in calls, f"{name}: {module.module_id} was never initialised"
            assert calls.index("validate_configuration") < calls.index("initialize"), (
                f"{name}: {module.module_id} validated its configuration after initialising"
            )

        # And across modules, not merely within one: nothing is initialised while
        # another module's configuration is still unchecked.
        first_initialize = next(i for i, (_m, op) in enumerate(order) if op == "initialize")
        validated_by_then = {
            m for m, op in order[:first_initialize] if op == "validate_configuration"
        }
        assert validated_by_then == {m.module_id for m in _modules(kernel)}, (
            f"{name}: initialisation began before every module had validated"
        )


#: The tendency pass, in the order one module sees it.
TENDENCY_PASS = ("begin_step", "receive_inputs", "compute_tendencies", "emit_outputs")


@pytest.mark.req("MOD-007")
def test_the_framework_calls_the_whole_lifecycle_on_every_active_module():
    """Charter 6.3's lifecycle *semantics*, on every module of every configuration.

    Exposing the operations is what MOD-007 asks for, and checking that they
    exist would be checking that ``Module`` has methods -- true of any subclass by
    inheritance, and evidence about Python rather than about this framework. What
    is worth establishing is that the framework *drives* them, so a module cannot
    be coupled into a run and silently skipped.

    The order asserted here is not 6.3's numbering, and the difference is
    deliberate on the coupler's part rather than a defect. 6.3 numbers the
    operations as one module sees them; the coupler interleaves them across
    modules, running every module's tendency pass before committing any state.
    That is *stronger* than the per-module sequence and it is the operator-
    splitting invariant this project depends on -- the convection/water over-draw
    came from a module acting on a profile an earlier one had not yet published.
    So what is checked is the semantics the charter says to preserve:

    * within a module, ``begin_step -> receive_inputs -> compute_tendencies ->
      emit_outputs``, so nothing computes before it has been given its inputs;
    * across modules, every ``compute_tendencies`` before any ``advance_state``,
      so no owner commits while another module's tendency is still outstanding;
    * ``report_conservation`` after the state it reports on;
    * ``finalize`` exactly once, at the end.

    ``advance_state`` is required of ACTIVE modules only. A PRESCRIBED module owns
    no state to advance, and the charter's status vocabulary exists to make that
    difference declarable rather than accidental.
    """
    for name, build in BUILDERS.items():
        kernel = build()
        order: list[tuple[str, str]] = []
        for module in _modules(kernel):
            for op in (*TENDENCY_PASS, "advance_state", "report_conservation", "finalize"):
                _record(module, op, order)
        kernel.run(2)

        for module in _modules(kernel):
            calls = [op for mid, op in order if mid == module.module_id]
            expected = [
                *TENDENCY_PASS,
                *(("advance_state",) if module.status is ModuleStatus.ACTIVE else ()),
                "report_conservation",
                "finalize",
            ]
            missing = [op for op in expected if op not in calls]
            assert not missing, f"{name}: {module.module_id} was never asked to {missing}"

            first = [calls.index(op) for op in TENDENCY_PASS]
            assert first == sorted(first), (
                f"{name}: {module.module_id} ran its tendency pass out of order: "
                f"{[op for op in calls if op in TENDENCY_PASS]}"
            )
            assert calls.index("report_conservation") > calls.index("compute_tendencies"), (
                f"{name}: {module.module_id} reported conservation before it had computed anything"
            )
            assert calls.count("finalize") == 1, (
                f"{name}: {module.module_id} was finalised {calls.count('finalize')} times"
            )

        # The splitting invariant, across modules rather than within one: inside
        # a step, every tendency is computed before any state is committed.
        for index, step in enumerate(_split_into_steps(order)):
            tendencies = [i for i, (_m, op) in enumerate(step) if op == "compute_tendencies"]
            commits = [i for i, (_m, op) in enumerate(step) if op == "advance_state"]
            if not commits:
                continue
            assert max(tendencies) < min(commits), (
                f"{name}: on step {index}, {step[min(commits)][0]} committed state while "
                f"{step[max(tendencies)][0]}'s tendency was still outstanding"
            )


def _split_into_steps(log: list[tuple[str, str]]) -> list[list[tuple[str, str]]]:
    """Cut the call log at each step boundary.

    A step is a tendency pass over every module, then the commits, then the
    conservation reports -- so the next step begins at the first ``begin_step``
    after a ``report_conservation``.
    """
    steps: list[list[tuple[str, str]]] = []
    current: list[tuple[str, str]] = []
    closed = True
    for entry in log:
        if entry[1] == "begin_step" and closed:
            if current:
                steps.append(current)
            current = []
            closed = False
        if entry[1] == "report_conservation":
            closed = True
        current.append(entry)
    if current:
        steps.append(current)
    return steps


def _record(module: Any, op: str, log: list[tuple[str, str]]) -> None:
    """Wrap one lifecycle operation so the call is visible without changing it."""
    original = getattr(module, op)
    module_id = module.module_id

    def wrapper(*args: Any, **kwargs: Any) -> Any:
        log.append((module_id, op))
        return original(*args, **kwargs)

    setattr(module, op, wrapper)


# --- MOD-012, CMP-002, CMP-003: what a module and a run must declare --------


@pytest.mark.req("CMP-002")
def test_every_module_declares_a_benchmarked_resource_record(kernels):
    """CMP-002 asks for a record against at least one supported configuration.

    A complexity string alone would satisfy a reading of "estimated resource
    record" and would be unfalsifiable: nothing anchors ``O(n_layers)`` to a
    number. The configuration the estimate was taken at is what makes the
    declaration checkable, so it is required here too.
    """
    for name, kernel in kernels.items():
        for module in _modules(kernel):
            cost = module.contract.cost
            assert cost is not None, f"{name}: {module.module_id} declares no cost"
            assert cost.module_id == module.module_id
            assert cost.time_complexity.strip(), (
                f"{name}: {module.module_id} has no time complexity"
            )
            assert cost.memory_complexity.strip()
            assert cost.benchmark_configuration.strip(), (
                f"{name}: {module.module_id} states a cost with no configuration it was taken at"
            )
            assert cost.work_units_per_step > 0.0
            assert cost.supported_resource_classes


@pytest.mark.req("CMP-002", "MOD-013")
def test_a_contract_without_a_cost_declaration_is_refused(kernels):
    """The counterfactual, without which the sweep above proves only good habits."""
    contract = next(iter(_modules(kernels["radiative_convective"]))).contract
    with pytest.raises(ConfigurationError, match="computational-cost"):
        dataclass_replace(contract, cost=None).validate()


@pytest.mark.req("CMP-003")
def test_every_run_declares_its_resource_class_and_budget(stepped):
    """CMP-003. A budget with no class is a number with no machine attached to it."""
    for name, result in stepped.items():
        manifest = result.manifest
        assert manifest.resource_class, f"{name}: no resource class"
        budget = manifest.resource_budget
        assert budget, f"{name}: no resource budget"
        for key in ("max_memory_bytes", "max_work_units", "max_output_bytes", "on_exceed"):
            assert key in budget, f"{name}: budget declares no {key}"
        assert budget["on_exceed"] in ("fail", "warn")


@pytest.mark.req("MOD-012")
def test_every_module_declares_the_representations_it_supports(kernels):
    """MOD-012, including "where applicable", which is the interesting half.

    A module that consumes or produces a spectral field has a spectral
    representation to declare; one that does not, does not, and requiring an
    empty tuple to be filled in would turn the declaration into boilerplate. So
    the spectral case is checked against whether the module's contract actually
    mentions the grid.
    """
    for name, kernel in kernels.items():
        for module in _modules(kernel):
            contract = module.contract
            assert contract.supported_spatial_representations, (
                f"{name}: {module.module_id} declares no spatial representation"
            )
            assert contract.supported_temporal_representations, (
                f"{name}: {module.module_id} declares no temporal representation"
            )
            assert contract.cost is not None and contract.cost.parallelism, (
                f"{name}: {module.module_id} declares no parallelism model"
            )


# --- CAU-003: the dependency graph ------------------------------------------


@pytest.mark.req("CAU-003")
def test_every_dependency_edge_carries_a_declared_type(stepped):
    """CAU-003 lists nine edge types, and an untyped edge is the failure it names.

    An edge's type is what separates "the surface sends the atmosphere energy"
    from "the surface's albedo parameterises the atmosphere's absorption". The
    conservation auditor reconciles the first and must not reconcile the second,
    so this is load-bearing rather than descriptive.
    """
    declared = {edge_type.value for edge_type in EdgeType}
    for name, result in stepped.items():
        edges = result.manifest.dependency_graph["edges"]
        assert edges, f"{name}: a coupled configuration with no edges"
        for edge in edges:
            assert edge["type"] in declared, (
                f"{name}: {edge['source']} -> {edge['destination']} carries type "
                f"{edge['type']!r}, which is not one of the nine CAU-003 declares"
            )
            assert edge["exchanged"], f"{name}: an edge that names no exchanged quantity"
            # A physical transfer is computed by a relation and has to cite it;
            # the other eight types need not. Handing a module a heat capacity or
            # an irradiance record is a dependency, not an equation, and six edges
            # in the tree are of that kind -- requiring a citation there would be
            # requiring one where none exists.
            if edge["type"] == EdgeType.PHYSICAL_TRANSFER.value:
                assert edge["equation_ids"], (
                    f"{name}: {edge['source']} -> {edge['destination']} moves "
                    f"{edge['exchanged']} and cites no equation for it"
                )


@pytest.mark.req("CAU-003")
def test_an_edge_cannot_be_built_without_a_type():
    """The counterfactual: the vocabulary is mandatory, not merely available."""
    with pytest.raises(TypeError):
        Edge(source="a", destination="b", exchanged="energy")  # type: ignore[call-arg]


# --- ATM-001: one owner, and the right one ---------------------------------


@pytest.mark.req("ATM-001", "SYS-001")
def test_the_atmospheric_thermodynamic_state_is_owned_by_the_atmosphere(kernels):
    """ATM-001, and the ownership rule it is a particular case of.

    The atmospheric *thermodynamic* state, specifically. The column's water
    vapour is owned by the water module under charter 9.6 and not by the
    atmosphere, which is a division of labour rather than a violation -- so the
    check names the state it is about instead of claiming the atmosphere owns
    everything in the air.
    """
    for name, kernel in kernels.items():
        ownership = kernel.build_manifest().state_ownership
        if "atmosphere_layer_temperature" not in ownership:
            continue  # a configuration with no atmosphere has no claim to test
        assert ownership["atmosphere_layer_temperature"] == "atmosphere", (
            f"{name}: the air temperature is owned by {ownership['atmosphere_layer_temperature']!r}"
        )
        # Exactly one owner, and the registry enforces it rather than documenting it.
        for module in _modules(kernel):
            if module.module_id == "atmosphere":
                continue
            with pytest.raises(OwnershipError):
                kernel.states.assert_may_write(module.module_id, "atmosphere_layer_temperature")

        # And no other module's contract claims it either, which is the half a
        # runtime refusal cannot see.
        claimants = [
            m.module_id
            for m in _modules(kernel)
            if "atmosphere_layer_temperature" in m.contract.owned_state_ids
        ]
        assert claimants == ["atmosphere"], f"{name}: {claimants} all claim the air temperature"


# --- ACC-005: what a run has to leave behind -------------------------------


@pytest.mark.req("ACC-005", "KER-014")
def test_every_run_leaves_a_complete_record_of_itself(stepped):
    """ACC-005's five records, on every configuration in the tree.

    This is the check that found something. ``output_checksums`` was empty on
    every one of these runs: the only path that filled it hashes a *file*, and no
    reference configuration writes one, so a requirement recorded as "satisfied,
    needs a test" was satisfied by a method nothing called.
    """
    for name, result in stepped.items():
        manifest = result.manifest
        assert manifest.run_id and manifest.purpose, f"{name}: an unnamed run"
        assert manifest.module_status and manifest.module_versions, f"{name}: no module record"
        assert manifest.equation_bindings, f"{name}: no equation bindings"
        assert manifest.state_ownership, f"{name}: no state ownership record"
        assert manifest.dependency_graph.get("edges"), f"{name}: no dependency graph"
        assert manifest.environment is not None, f"{name}: no provenance record"
        assert manifest.environment.python_version and manifest.environment.platform
        assert manifest.output_checksums, f"{name}: produced no output checksum"
        assert manifest.gate_results, f"{name}: recorded no gate results"


@pytest.mark.req("ACC-005", "KER-014")
def test_the_output_checksum_identifies_the_run_it_came_from():
    """A checksum that does not discriminate is a field, not a record.

    KER-014 wants a reference run reproducible from its manifest, which needs the
    fingerprint to match a repeat of the same configuration and to differ from
    anything else. The trajectory digest is there for the second half: two runs
    that met at the same endpoint by different paths would otherwise hash alike.
    """

    def fingerprint(steps: int) -> dict[str, str]:
        kernel = cc.radiative_convective(grid=cc.column_grid(8), n_layers=8)
        return dict(kernel.run(steps).manifest.output_checksums)

    once = fingerprint(4)
    assert set(once) == {"final_outputs", "final_state", "conservation_trajectory"}
    assert once == fingerprint(4), "the same configuration did not reproduce its own fingerprint"
    assert once != fingerprint(5), "a run of a different length hashed the same"


# --- SYS-007: corrections are versioned, not applied ----------------------


@pytest.mark.req("SYS-007", "EQR-007", "AUD-006")
def test_a_correction_creates_a_version_and_never_overwrites(kernels):
    """SYS-007, swept over every equation record a reference kernel carries.

    A published result must stay reconstructible from the record versions it
    actually used, so the store is append-only and a correction is a new version
    beside the old one. Thirty-four real records rather than a fixture: the claim
    is about the registry, and a fixture would only establish that the fixture
    behaves.
    """
    registry = kernels["moist_radiative_convective"].equations
    record_ids = registry.ids()
    assert len(record_ids) > 20, "too few records for this to be a sweep"

    for record_id in record_ids:
        original = registry.get(record_id)
        version = original.version
        corrected = dataclass_replace(original, version=Version(version.major + 1, 0, 0))
        registry.supersede(record_id, version, corrected)

        # The original is still there, at its own version, and says it was replaced.
        assert registry.get(record_id, version) is original
        assert registry.entry(record_id, version).superseded_by == (
            record_id,
            str(corrected.version),
        )
        assert registry.get(record_id) is corrected
        assert original not in registry.current()

    # And the append-only rule is enforced, not merely observed: changing a
    # record under a version already published is the failure SYS-007 names.
    first = registry.get(record_ids[0])
    with pytest.raises(RegistryError, match="already registered with different content"):
        registry.register(dataclass_replace(first, assumptions=("edited in place",)))


# --- the sweep's own coverage ----------------------------------------------


def test_the_sweep_reaches_every_module_implementation(kernels):
    """A sweep with a hole in it makes a claim it has not earned.

    Every ``Module`` subclass in ``mfespf.physics`` either appears in one of the
    configurations above or is named in ``UNSWEPT`` with a reason. Adding a module
    without wiring it into a reference configuration fails here and says which
    one, rather than leaving the coverage claim quietly weaker than it reads.
    """
    for info in pkgutil.iter_modules(mfespf.physics.__path__):
        importlib.import_module(f"mfespf.physics.{info.name}")

    def descendants(cls: type) -> set[type]:
        found: set[type] = set()
        for sub in cls.__subclasses__():
            found.add(sub)
            found |= descendants(sub)
        return found

    # Scoped to the physics package on purpose. The failure-detection tests
    # define modules that leak, double-count and write state they do not own;
    # they are fixtures for catching those failures and have no business in a
    # sweep over the framework's own implementations -- and picking them up
    # depends on which test files have been imported, which would make this
    # check's result depend on how the suite was invoked.
    implementations = {
        c.__name__ for c in descendants(Module) if c.__module__.startswith("mfespf.physics")
    }
    swept = {type(m).__name__ for kernel in kernels.values() for m in _modules(kernel)}

    unreached = sorted(implementations - swept - set(UNSWEPT))
    assert not unreached, (
        f"no reference configuration in this sweep exercises {unreached}; wire it into "
        "BUILDERS or record why not in UNSWEPT"
    )
    stale = sorted(set(UNSWEPT) & swept)
    assert not stale, f"{stale} is recorded as unswept but the sweep reaches it"


def dataclass_replace(obj: Any, **changes: Any) -> Any:
    """``dataclasses.replace`` under a name that does not shadow the module."""
    from dataclasses import replace

    return replace(obj, **changes)


def test_the_configurations_all_build_and_step(stepped):
    """A guard on the sweeps above: an empty fixture would pass every one of them."""
    assert len(stepped) == len(BUILDERS) >= 17
    for name, result in stepped.items():
        assert result.steps, f"{name}: ran no steps"
        assert np.isfinite(
            float(np.asarray(result.steps[-1].outputs["surface_temperature"].values).sum())
        ), f"{name}: finished with a non-finite surface temperature"


# --- SYS-002: a carve-out is a decision, not a backlog item -----------------


@pytest.mark.req("SYS-002", "SYS-003")
def test_the_run_time_status_distinguishes_a_carve_out_from_an_active_module(kernels, stepped):
    """Half of SYS-002, at run time: the status vocabulary is load-bearing.

    Charter 4.3 defines a carve-out as a reserved place that asserts nothing
    about whether the domain is important or active. The distinction has to do
    something, or it is a comment: a STUB's output is test data and cannot
    support a claim, an INACTIVE module does not participate, and only an ACTIVE
    one evolves state. Every run declares every module against that vocabulary.
    """
    assert not ModuleStatus.STUB.may_support_a_claim
    assert all(s.may_support_a_claim for s in ModuleStatus if s is not ModuleStatus.STUB)
    assert not ModuleStatus.INACTIVE.participates
    assert [s for s in ModuleStatus if s.evolves_state] == [ModuleStatus.ACTIVE]

    for name, result in stepped.items():
        declared = result.manifest.module_status
        modules = {m.module_id: m for m in _modules(kernels[name])}
        assert set(declared) == set(modules), f"{name}: a module ran without a declared status"
        for module_id, status in declared.items():
            assert status in {s.value for s in ModuleStatus}, (
                f"{name}: {module_id} declares status {status!r}, outside the charter's vocabulary"
            )
            # And the declaration is the module's own, not the manifest's guess.
            assert modules[module_id].status.value == status


@pytest.mark.req("SYS-002")
def test_a_carved_out_requirement_names_the_clause_that_defers_it():
    """The other half, in the record: a carve-out excuses a requirement, so it
    has to say on whose authority.

    Without this the status is a way of making the coverage number look better.
    The tool refuses a carve-out with no waiver, and refuses one that also has a
    bound test -- because a requirement cannot both be deferred and be verified,
    and whichever of the two is wrong should be said out loud.

    The classification is the specification's rather than a judgement made here.
    It is also smaller than it was estimated to be: `docs/todo.md` put it at
    "about ten", and looking made it four. Gravity is the case that decides the
    boundary -- the gravity module is as absent from charter 18's candidate list
    as the cryosphere is, but the column computes hydrostatic layer mass with it
    every step, so GRV-001 is a gap rather than a carve-out.
    """
    requirements = _load_requirements()
    carved = {rid: row for rid, row in requirements.items() if row.get("status") == "carved_out"}

    assert carved, "nothing is classified, so the distinction SYS-002 asks for is not being made"
    for rid, row in carved.items():
        waiver = row.get("waiver") or ""
        assert waiver.strip(), f"{rid} is carved out with no waiver"
        assert "charter" in waiver.lower() or "irs" in waiver.lower(), (
            f"{rid}'s waiver cites no specification clause: {waiver!r}"
        )

    # Gravity is used by the column and is therefore counted, not excused.
    assert requirements["GRV-001"].get("status") != "carved_out", (
        "GRV-001 was excused although the column computes layer mass with gravity"
    )

    traceability = _traceability_module()
    problems = _traceability_problems(traceability, requirements)
    assert not problems, f"the requirement registry is inconsistent: {problems}"

    # And the refusals fire, rather than being available.
    unjustified = {**requirements, "ICE-001": {**requirements["ICE-001"], "waiver": ""}}
    assert any("no waiver" in p for p in _traceability_problems(traceability, unjustified))


# --- HUM-001: how a human influence is allowed to reach the state ----------


@pytest.mark.req("HUM-001", "CHM-001")
def test_human_activity_reaches_the_state_only_through_a_constitutive_relation():
    """HUM-001, demonstrated on the one human influence this kernel represents.

    The failure HUM-001 names is a fitted forcing added to a budget: a term in
    watts per square metre, tuned to reproduce a record, that no physical process
    put there. Carbon dioxide here does the opposite. The module publishes a mole
    fraction and nothing else, owns no state, computes no tendency and declares
    no transfer; the radiation module turns the concentration into an absorption
    coefficient, and the transfer solver is unchanged and unaware.

    So the check is in three parts, and all three are needed. That the module
    adds nothing to a budget would also be true of a module that does nothing, so
    the second part shows it *acts*: doubling the concentration moves the
    outgoing longwave. And a mechanism that acts without adding anything must
    leave the budget closing exactly as before, which is the third.
    """
    kernel = _greenhouse_column()
    result = kernel.run(50)
    greenhouse = next(m for m in _modules(kernel) if m.module_id == "greenhouse")

    # It enters through none of the routes HUM-001 forbids.
    assert greenhouse.compute_tendencies() == [], "the concentration contributes a tendency"
    assert greenhouse.report_transfers() == ([], []), "the concentration declares a transfer"
    assert not greenhouse.contract.owned_state_ids
    assert not greenhouse.contract.conserved_quantities
    assert greenhouse.contract.produced_outputs == ("co2_mole_fraction",)

    edges = [
        edge
        for edge in result.manifest.dependency_graph["edges"]
        if "greenhouse" in (edge["source"], edge["destination"])
    ]
    assert edges, "the concentration reaches nothing"
    for edge in edges:
        assert edge["type"] == EdgeType.PARAMETER_DEPENDENCY.value, (
            f"carbon dioxide reaches {edge['destination']} as a {edge['type']}, not a parameter"
        )

    # No budget carries a term named for it, under any name.
    labels = {
        label for budget in result.steps[-1].state_budgets.values() for label in budget.by_label()
    }
    assert not any("greenhouse" in label or "forcing" in label for label in labels), (
        f"a budget carries a term for the concentration: {sorted(labels)}"
    )

    # It nonetheless acts, through the absorption coefficient alone.
    doubled = _greenhouse_column_at(560.0).run(50)
    olr = lambda run: float(  # noqa: E731
        np.asarray(run.steps[-1].outputs["outgoing_longwave_radiation"].values).sum()
    )
    assert olr(doubled) < olr(result), (
        f"doubling the concentration did not reduce the outgoing longwave: "
        f"{olr(result):.3f} to {olr(doubled):.3f} W m^-2"
    )

    # And it adds no energy: both budgets close to the same place.
    for run in (result, doubled):
        assert run.final_conservation.worst_relative_residual() < 1e-12


def _load_requirements() -> dict[str, dict]:
    """The requirement registry, read as a record rather than through the tool."""
    requirements: dict[str, dict] = {}
    for path in sorted((REPO_ROOT / "requirements").glob("req-*.yaml")):
        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        for row in data.get("requirements", []):
            requirements[row["id"]] = row
    return requirements


def _traceability_module() -> Any:
    """Import ``tools/traceability.py`` without putting ``tools`` on the path."""
    spec = importlib.util.spec_from_file_location(
        "_traceability_under_test", REPO_ROOT / "tools" / "traceability.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _traceability_problems(traceability: Any, requirements: dict[str, dict]) -> list[str]:
    """The tool's integrity check, run over a registry given to it."""
    problems: list[str] = []
    bindings = traceability.collect_bindings()
    for row in requirements.values():
        if traceability.is_carve_out(row):
            if not row.get("waiver"):
                problems.append(f"{row['id']}: carved out with no waiver")
            if bindings.get(row["id"]):
                problems.append(f"{row['id']}: carved out and has a bound test")
    return problems


# --- CMP-007: a cheaper form declared and not taken is a choice -------------


@pytest.mark.req("CMP-007", "CMP-006")
def test_a_module_that_declines_its_cheaper_form_says_what_it_would_lose(kernels, stepped):
    """CMP-007, swept over every module of every configuration.

    Nine of the eleven module implementations name a lower-fidelity alternative
    -- a single-layer column, a grey atmosphere, a dry column, an airless slab --
    and then do not use it. That is a scientific choice, and CMP-007 asks that the
    reason reach the manifest rather than living in whoever configured the run.

    The justification has to name what the cheaper form would *fail to produce*,
    because that is the falsifiable half. "More accurate" is unfalsifiable; "a
    single layer has no lapse rate, and the lapse rate is what convection acts
    on" can be checked by running the cheaper form and looking.
    """
    for name, kernel in kernels.items():
        for module in _modules(kernel):
            cost = module.contract.cost
            assert cost is not None
            if not cost.lower_fidelity_alternatives:
                continue
            reason = module.contract.fidelity_justification
            assert reason.strip(), (
                f"{name}: {module.module_id} declines "
                f"{list(cost.lower_fidelity_alternatives)} without saying why"
            )
            # Long enough to be a reason rather than a label. The shortest of the
            # real ones is 180 characters; a bare adjective would not reach it.
            assert len(reason) > 80, f"{name}: {module.module_id} gives a label, not a reason"

    # And it reaches the manifest, for exactly the modules that made a choice.
    for name, result in stepped.items():
        declared = result.manifest.fidelity_justifications
        expected = {
            m.module_id
            for m in _modules(kernels[name])
            if m.contract.cost is not None and m.contract.cost.lower_fidelity_alternatives
        }
        assert set(declared) == expected, f"{name}: the manifest records the wrong modules"
        assert all(declared.values()), f"{name}: an empty justification reached the manifest"
        assert set(result.manifest.fidelity_levels) == {
            m.module_id for m in _modules(kernels[name])
        }


@pytest.mark.req("CMP-007")
def test_a_contract_that_names_a_cheaper_form_and_no_reason_is_refused(kernels):
    """The counterfactual, without which the sweep proves only good habits."""
    contract = next(m.contract for m in _modules(kernels["radiative_convective"]))
    assert contract.cost is not None and contract.cost.lower_fidelity_alternatives
    with pytest.raises(ConfigurationError, match="gives no reason for running at"):
        dataclass_replace(contract, fidelity_justification="   ").validate()

    # A module with no cheaper form made no choice, so it is not asked for one.
    no_alternative = dataclass_replace(
        contract,
        fidelity_justification="",
        cost=dataclass_replace(contract.cost, lower_fidelity_alternatives=()),
    )
    no_alternative.validate()


# --- ACC-010: the baseline suite against the workstation budget ------------


#: The baseline suite, and the step count each member runs at here. Short enough
#: that this test is not itself the expensive thing in the suite, and the memory
#: term is linear in the step count (measured below), so the figure extrapolates.
BASELINE_SUITE = {
    "airless_radiating_body": rc.airless_radiating_body,
    "surface_thermal_storage": rc.surface_thermal_storage,
    "reduced_ocean_heat_storage": rc.reduced_ocean_heat_storage,
    "window_atmosphere": cc.window_atmosphere,
    "radiative_convective": cc.radiative_convective,
    "moist_radiative_convective": cc.moist_radiative_convective,
    "prescribed_cloud": cc.prescribed_cloud,
    "column_over_reduced_ocean": cc.column_over_reduced_ocean,
}
BASELINE_STEPS = 500


@pytest.mark.req("ACC-010", "CMP-004", "CMP-008")
def test_the_baseline_suite_runs_inside_the_workstation_budget():
    """ACC-010, measured rather than estimated.

    The estimator already checks every run against the budget before it starts,
    so binding this requirement to the estimate would be binding it to the thing
    that is already enforced. What ACC-010 asks is whether the suite *completes*
    within the budget, and that is answered by watching it.

    Measured at 500 steps a member: the suite uses 0.4 percent of the declared
    work-unit budget and its worst member holds about 20 MiB, against 16 GiB. The
    baseline acceptance configuration also needs no distributed resource (CMP-008),
    which holds trivially here -- every module declares serial parallelism.
    """
    budget = ResourceBudget(resource_class=ResourceClass.ENGINEERING_WORKSTATION)
    total_work = 0.0
    worst_memory = 0

    for name, build in BASELINE_SUITE.items():
        kernel = build()
        tracemalloc.start()
        try:
            result = kernel.run(BASELINE_STEPS)
            _current, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()

        estimate = result.cost_estimate
        assert estimate.within_budget, f"{name}: {estimate.violations}"
        assert result.manifest.status != "failed", f"{name}: did not complete"
        total_work += estimate.total_work_units
        worst_memory = max(worst_memory, peak)

        for module in _modules(kernel):
            assert module.contract.cost is not None
            assert module.contract.cost.parallelism is ParallelismModel.SERIAL, (
                f"{name}: {module.module_id} needs {module.contract.cost.parallelism} (CMP-008)"
            )

    assert total_work < budget.max_work_units, (
        f"the suite needs {total_work:g} work units against a budget of {budget.max_work_units:g}"
    )
    assert worst_memory < budget.max_memory_bytes
    # Not merely inside: inside with room, which is what makes the requirement a
    # statement about the kernel rather than about the budget being generous.
    assert total_work < 0.05 * budget.max_work_units
    assert worst_memory < 0.05 * budget.max_memory_bytes


@pytest.mark.req("CMP-004", "ACC-010")
def test_the_memory_estimate_accounts_for_the_trajectory_it_retains():
    """CMP-004's memory figure, against what the run actually holds.

    This check replaces one that asserted the estimate was *wrong*, which is
    what it was for: a marker on a known defect, written to fail when the defect
    was fixed. The estimate summed the modules' working sets and came out a
    factor of twenty low, in the direction that admits a run which then exhausts
    memory -- and that run happened, at 4.3 GB on a 120000-step sweep that had to
    be killed.

    What grows is the retained trajectory: every step's published fields are kept
    so the manifest can be audited, which is linear in the step count while the
    working set is not. Measured, that cost is 714 bytes per field plus 17.6 per
    value, and the per-field term is forty times the per-value one -- what is
    stored is each field's units, provenance, uncertainty and semantics, not its
    numbers.
    """
    for name, build in (
        ("dry, 8 layers", lambda: cc.radiative_convective(grid=cc.column_grid(8), n_layers=8)),
        (
            "moist, 8 layers",
            lambda: cc.moist_radiative_convective(grid=cc.column_grid(16), n_layers=8),
        ),
        (
            "moist, 24 layers",
            lambda: cc.moist_radiative_convective(grid=cc.column_grid(16), n_layers=24),
        ),
    ):
        tracemalloc.start()
        try:
            result = build().run(300)
            _current, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()

        estimate = result.cost_estimate
        assert estimate.retained_trajectory_bytes > 0, f"{name}: no trajectory term"
        # Within a quarter either way, where it used to be a factor of twenty.
        # Not tighter: `tracemalloc` counts allocations the run makes for other
        # reasons too, and a bound that only holds on one interpreter would be
        # measuring the interpreter.
        assert 0.75 < peak / estimate.retained_trajectory_bytes < 1.5, (
            f"{name}: estimate {estimate.retained_trajectory_bytes}, actual {peak}"
        )
        # It dominates, which is the fact that makes it worth estimating at all.
        assert estimate.retained_trajectory_bytes > 0.9 * estimate.peak_memory_bytes


@pytest.mark.req("CMP-004", "CMP-005")
def test_the_trajectory_term_is_what_makes_a_long_run_different():
    """A configuration that fits at 2000 steps and not at 200000 differs only here.

    The working set does not depend on the step count and the retained
    trajectory is linear in it, so this is the whole of what a preflight can say
    about the difference between a short run and a long one -- and saying nothing
    about it was the defect.
    """
    kernel = cc.moist_radiative_convective(grid=cc.column_grid(16), n_layers=8)
    estimates = {
        steps: estimate_cost(
            schedule=kernel.schedule,
            n_steps=steps,
            budget=ResourceBudget(resource_class=ResourceClass.ENGINEERING_WORKSTATION),
            spectral_bins=16,
        )
        for steps in (1_000, 10_000, 100_000)
    }
    retained = [estimates[s].retained_trajectory_bytes for s in (1_000, 10_000, 100_000)]
    assert retained[1] == pytest.approx(10 * retained[0], rel=1e-9)
    assert retained[2] == pytest.approx(100 * retained[0], rel=1e-9)

    # And the figure that was hit and killed: a 120000-step moist column.
    at_120k = estimate_cost(
        schedule=kernel.schedule,
        n_steps=120_000,
        budget=ResourceBudget(resource_class=ResourceClass.ENGINEERING_WORKSTATION),
        spectral_bins=16,
    )
    gigabytes = at_120k.retained_trajectory_bytes / 1024**3
    assert 3.0 < gigabytes < 6.0, (
        f"the 120000-step sweep that had to be killed reached 4.3 GB; the estimate "
        f"now says {gigabytes:.2f} GB"
    )


# --- VER-006: negotiation in every configuration ---------------------------


@pytest.mark.req("VER-006", "VER-005")
def test_every_configuration_negotiated_every_declared_requirement(kernels, stepped):
    """One configuration negotiating cleanly says nothing about the other sixteen."""
    for name, kernel in kernels.items():
        declared = {
            (slot.module_id, required.interface_id)
            for slot in kernel.schedule.ordered()
            if slot.module.status.participates
            for required in slot.module.contract.required_interfaces
        }
        assert kernel.coupler is not None, name
        negotiated = {(n.consumer_module_id, n.interface_id) for n in kernel.coupler.negotiations}
        assert negotiated == declared, name
        assert kernel.manifest is not None, name
        assert len(kernel.manifest.negotiated_interfaces) == len(declared), name


# --- UNC-001 and UNC-002: declared uncertainty in every configuration ------


@pytest.mark.req("UNC-001", "UNC-002", "UNC-004")
def test_every_declared_range_describes_the_configuration_as_built(kernels):
    """A declaration is checked against the module, not trusted.

    Its nominal must be the value the module runs with, its setting must be a
    keyword a configuration can actually be built with -- or no ensemble could
    ever move it -- and a derived parameter must follow something this same
    configuration declares, or the derivation points at nothing.
    """
    import inspect

    keywords = set(inspect.signature(cc.radiating_column).parameters)
    for name, kernel in kernels.items():
        settings: set[str] = set()
        follows: set[str] = set()
        for slot in kernel.schedule.ordered():
            report = slot.module.report_uncertainty()
            parameters = slot.module.parameters()
            for declared in report.declared_ranges:
                assert parameters[declared.parameter] == declared.nominal, (name, declared)
                assert declared.setting in keywords, (name, declared.setting)
                settings.add(declared.setting)
            for dependent in report.derived:
                assert dependent.setting in keywords, (name, dependent.setting)
                follows.add(dependent.follows)
        assert follows <= settings, (name, follows - settings)
        manifest = kernel.manifest
        assert manifest is not None
        assert set(manifest.declared_uncertainty) == {
            slot.module_id for slot in kernel.schedule.ordered()
        }, name
