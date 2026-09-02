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

Version 0.1.0 — the **Stage 0 accounting kernel and Stage 1 airless radiating
body** of the charter's validation ladder (charter §19), on top of the full
Version 1.0 interface spine.

| Ladder stage | State |
| --- | --- |
| 0 — units, coordinates, transfers, conservation accounting, restart | implemented |
| 1 — airless radiating body: geometry, absorbed radiation, thermal emission, analytic equilibrium, transient | implemented |
| 2 — atmospheric radiation (transparent → absorbing → emitting, multilayer) | not started |
| 3 — atmospheric thermodynamics | not started |
| 4 — clouds | not started |
| 5 — thermal storage and measured response timescales | partially: single-reservoir response measured |
| 6 — historical replay | machinery present (withheld-data path), no datasets |
| 7 — hindcast and registered prediction | machinery present (freeze records), not exercised |

189 mandatory requirements are tracked; 122 currently have a bound acceptance
test. See [`docs/traceability.md`](docs/traceability.md) — the uncovered ones are
listed there rather than left implicit.

## Install and run

```bash
python -m pip install -e ".[dev]"     # Python 3.11+
python -m pytest                       # 136 tests, ~18 s on a workstation
mfespf run airless --steps 4800        # 200 simulated days
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

Configurations: `airless`, `airless-banded`, `null-control`.
Add `--manifest run.json` to write the run record.

## What the kernel actually does

A zero-dimensional airless sphere in radiative balance, coupled through the full
interface machinery:

```
solar ──irradiance──▶ orbit ──energy──▶ surface ──reflected──▶ space
(PRESCRIBED)          (ACTIVE)          (ACTIVE)  ──emitted───▶ (ACTIVE)
                      owns the TOA      owns T_s               owns E_out
                      boundary
```

The physics is deliberately small; the point is that it is **accounted**. The
coupled energy budget closes to ~1e-16 relative at every coupling interval, and
the run reaches the closed-form equilibrium temperature to 1e-9 relative and the
analytic relaxation timescale to 2%.

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
