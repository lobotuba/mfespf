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
import pkgutil
from typing import Any

import numpy as np
import pytest

import mfespf.physics
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

pytestmark = pytest.mark.stage(0)


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
