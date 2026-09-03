# MFESPF

**Modular Falsifiable Earth-System Physics Framework** — a computational
laboratory for deriving, testing, and falsifying physical explanations of
Earth-system behaviour.

> Every subsystem is optional.
> Every approximation is visible.
> Every claim is measurable.
> Every prediction is falsifiable.
> Every failure is reportable.

The framework exists to evaluate hypotheses through conservation laws, physical
mechanisms, numerical validation, and comparison with observations. Observations
falsify or support claims; they do not dictate governing equations.

The normative documents live in [`specifications/`](specifications/): the
project charter (v0.2 + Amendment 1) and the v1.0 Interface Requirements
Specification. **Where this README and a specification disagree, the
specification governs.**

## Status

Version 0.3.0 — the charter's validation ladder (§19) through **Stage 3**:
the accounting kernel, the airless radiating body, the radiating atmospheric
column, and atmospheric thermodynamics with convection and water — on top of the
full Version 1.0 interface spine.

| Ladder stage | State |
| --- | --- |
| 0 — units, coordinates, transfers, conservation accounting, restart | implemented |
| 1 — airless radiating body: geometry, absorbed radiation, thermal emission, analytic equilibrium, transient | implemented |
| 2 — atmospheric radiation: transparent, absorbing and emitting, multilayer, variable spectrum | implemented |
| 3 — atmospheric thermodynamics: hydrostatic balance, dry convection, water vapour, phase-change energy | implemented |
| 4 — clouds | not started |
| 5 — thermal storage and measured response timescales | partially: single-reservoir response measured |
| 6 — historical replay | machinery present (withheld-data path), no datasets |
| 7 — hindcast and registered prediction | machinery present (freeze records), not exercised |

189 mandatory requirements are tracked; 129 currently have a bound acceptance
test. See [`docs/traceability.md`](docs/traceability.md) — the uncovered ones are
listed there rather than left implicit.

## Install and run

```bash
python -m pip install -e ".[dev]"     # Python 3.11+
python -m pytest                       # 185 tests, ~90 s on a workstation
mfespf run airless --steps 4800        # Stage 1: 200 simulated days
mfespf run window --steps 6000 --heat-capacity 1e6            # Stage 2 column
mfespf run moist  --steps 4000 --heat-capacity 1e6            # Stage 3 moist column
```

```
$ mfespf run airless --steps 240
run STAGE1-AIRLESS-001 [complete]
  gates:
    PASS  configuration_validation
    PASS  conservation_reconciliation
    PASS  resource_preflight
  worst relative conservation residual: 8.938e-16
```

Configurations: `airless`, `airless-banded`, `null-control` (Stage 1);
`transparent`, `absorbing`, `window`, `grey-layer` (Stage 2);
`radiative-convective`, `moist` (Stage 3).
Add `--manifest run.json` to write the run record.

## What the kernel actually does

**Stage 1** — a zero-dimensional airless sphere in radiative balance:

```
solar ──irradiance──▶ orbit ──energy──▶ surface ──reflected──▶ space
(PRESCRIBED)          (ACTIVE)          (ACTIVE)  ──emitted───▶ (ACTIVE)
                      owns the TOA      owns T_s               owns E_out
                      boundary
```

**Stage 2** — a hydrostatic radiating column above it:

```
solar ─▶ orbit ─▶ radiation ─▶ atmosphere ⇄ surface        atmosphere ─▶ space
                  solves the   owns T_a[]   owns T_s        owns the two
                  whole column                              column interfaces
```

**Stage 3** — thermodynamics, convection and water on the same column:

```
solar ─▶ orbit ─▶ radiation ─▶ water ─▶ convection ─▶ atmosphere ─▶ surface ─▶ space
                               owns     sensible       commits      owns T_s and
                               vapour   flux and       four         its water store
                               and its  dry adiabatic  labeled
                               latent   adjustment     tendencies
                               energy
```

Three modules now contribute labeled tendencies to one state, and the framework
moves a conserved quantity that is not energy. Water leaves the surface as
vapour carrying latent energy, condenses aloft releasing it as sensible heat,
and returns as precipitation carrying none — with **both** budgets required to
close at the same time, which is what CON-007 asks of a phase change.

The physics is deliberately small; the point is that it is **accounted**. The
coupled energy budget closes to ~1e-16 relative at every coupling interval.
Stage 1 reaches the closed-form equilibrium to 1e-9 relative and the analytic
relaxation timescale to 2%. Stage 2's single grey layer reproduces the
closed-form greenhouse temperatures — 295.42 K against 295.62 K analytic, with
the 0.07% gap attributable to a measured spectral-overlap term rather than to
tolerance.

There is no shortwave/longwave split in the Stage 2 physics. There is radiation
in band *b*, with two sources: sunlight at the top and thermal emission from the
layers and surface. A greenhouse is then just an absorption coefficient that is
small where the solar spectrum lives and large where the terrestrial one does —
which is what the real one is. The SW/LW *diagnostics* the charter requires are
recovered by partitioning bands after the solve.

## Architecture

The interface spine is the substantial part. See
[`docs/architecture.md`](docs/architecture.md) for the design and the reasoning.

| Concern | Module | Enforces |
| --- | --- | --- |
| Dimensions and units | [`core/units.py`](src/mfespf/core/units.py) | incompatible quantities fail; conversions are logged, never silent |
| Field metadata | [`core/quantities.py`](src/mfespf/core/quantities.py) | units, sign convention, time semantics and support travel with every value |
| Interface versions | [`core/versioning.py`](src/mfespf/core/versioning.py) | a changed sign convention cannot hide inside a compatible patch |
| State ownership | [`registries/states.py`](src/mfespf/registries/states.py) | one owner per state; others submit *labeled* tendencies |
| Equations | [`registries/equations.py`](src/mfespf/registries/equations.py) | a calibrated correction cannot be declared a fundamental law |
| Data provenance | [`registries/datasets.py`](src/mfespf/registries/datasets.py) | calibration data cannot be re-reported as independent validation |
| Claims | [`registries/claims.py`](src/mfespf/registries/claims.py) | no claim without a failure criterion; thresholds freeze before evaluation |
| Transfers | [`coupling/transfers.py`](src/mfespf/coupling/transfers.py) | both parties record each exchange; overlapping records are double counts |
| Conservation | [`coupling/conservation.py`](src/mfespf/coupling/conservation.py) | numerical sources tracked apart from physical ones |
| Dependency graph | [`coupling/graph.py`](src/mfespf/coupling/graph.py) | undeclared feedback cycles fail; edges can be disabled for controls |
| Multirate coupling | [`coupling/coupler.py`](src/mfespf/coupling/coupler.py) | modules on different clocks; the hold is declared and reported |
| Spectral interface | [`contracts/spectral.py`](src/mfespf/contracts/spectral.py) | no fixed spectrum; regridding conserves energy and logs its residual |
| Cost | [`contracts/cost.py`](src/mfespf/contracts/cost.py) | preflight refuses a run that will not fit its declared budget |
| Withheld data | [`evaluation/withheld.py`](src/mfespf/evaluation/withheld.py) | unreadable until the claim is frozen |
| Convergence | [`evaluation/convergence.py`](src/mfespf/evaluation/convergence.py) | an answer that moves with the time step is a solver property |
| Run audit | [`run/manifest.py`](src/mfespf/run/manifest.py) | configuration hashing, parent/child diffs, immutable evidence |
| Radiative transfer | [`physics/radiative_transfer.py`](src/mfespf/physics/radiative_transfer.py) | flux divergence closes by construction; scattering absent, not approximated |
| Atmospheric column | [`physics/atmosphere.py`](src/mfespf/physics/atmosphere.py) | owns layer temperatures; radiation contributes, the column commits |
| Thermodynamics | [`physics/thermodynamics.py`](src/mfespf/physics/thermodynamics.py) | convective adjustment conserves column enthalpy by construction |
| Water | [`physics/water.py`](src/mfespf/physics/water.py) | owns the latent energy, so water and energy close together |

### The design decisions worth knowing

**Reconciliation happens at coupling-interval boundaries, not every step.**
Under multirate coupling a daily module writes one transfer record where an
hourly one writes twenty-four. Comparing them per step reports a mismatch that
is an artifact of the schedule. Transfers accumulate and reconcile at the
synchronization point (charter §8.4), which is the least common multiple of the
module step multiples.

**Double counting is detected by interval overlap, not by count.** Many records
tiling an interval end to end are one exchange at finer granularity; two records
whose intervals *overlap* are the same energy entering the budget twice.

**The residual is scaled against the largest term in the budget, not the net
change.** A budget whose imports and exports nearly cancel is still wrong if it
is wrong by a percent of those flows.

**Space is a module.** Energy leaving the system accumulates in a boundary sink
that records its own side of every transfer, computed from the fluxes the
surface published rather than from the surface's internals. That independence is
what makes the sender/receiver reconciliation a real check.

**A disabled dependency edge withholds the field entirely.** It does not
substitute zero — a module that requires the input fails loudly, so a control
experiment cannot silently become a different experiment.

**Radiation runs before the state owners, and owners republish state after
advancing.** `Module.emit_state()` publishes committed state at initialisation
and after every `advance_state`, so a solver placed earlier in the splitting
order reads the state as of the start of *its* step rather than one step stale.
Fluxes still come from `emit_outputs()` mid-step, because a flux needs an
averaging interval and state does not.

**Emissivity follows from albedo, band by band.** That is Kirchhoff's law, and
it caught a real inconsistency: Stage 1's pairing of a 0.3 albedo with unit
emissivity is only self-consistent if the reflectance is confined to short
wavelengths. Read as a genuinely grey surface it gives 278 K, not 254.6 K. Both
are legitimate configurations; the framework now records which one a run used.

**The water module owns the latent energy, not the atmosphere.** Its energy
store is `L_v` times the vapour it holds. Without that store the budget could
not close between evaporation and condensation — energy would leave the surface
at one event and appear in the atmosphere at the other with nothing holding it
in between.

**Convective adjustment conserves enthalpy by construction, not to a
tolerance.** The mixed potential temperature of each unstable group is *defined*
as the value leaving `sum(C_i T_i)` unchanged. Convection crosses no module
boundary, so its contribution to the column total must be exactly zero — and
that is asserted per step, not assumed.

## Requirement traceability

Requirements are extracted from the IRS rather than retyped:

```bash
python tools/extract_requirements.py   # specifications/ -> requirements/*.yaml
python tools/traceability.py --check   # tests -> docs/traceability.md
```

Tests bind to requirements with markers, and the matrix is generated from them:

```python
@pytest.mark.req("CON-004", "KER-005")
def test_matched_transfers_reconcile_to_zero(): ...
```

A test citing a requirement ID that does not exist fails the build. Coverage
below 100% is reported, not enforced — a gate that failed for unimplemented
requirements would only encourage marking them done.

## Development method

Per charter §20, for each new module, equation or coupling: define the physical
claim, derive the governing equation, declare its validity domain, implement the
simplest valid form, verify units and limiting cases, verify conservation,
verify numerical convergence, compare against an independent process-level
observation, couple it, re-run the gates, run a null control, quantify
uncertainty, record the claim status — and preserve negative results.

## Licence

MIT. See [`LICENSE`](LICENSE).
