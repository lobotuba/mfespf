# Architecture

This document explains how the code realises the charter, and why particular
design choices were made. It is subordinate to
[`specifications/`](../specifications/); where they disagree, the specification
governs.

## The problem the architecture solves

The charter's hard requirement is not that the physics be good. It is that the
physics be **auditable**: that a reader can determine, from the artifacts a run
leaves behind, which mechanisms were active, which equations they used, where
every parameter came from, what the model was told before it was tested, and
whether the answer survives refinement.

Most of that cannot be achieved by discipline. An author who intends to keep
calibration and validation separate will occasionally fail to, and the failure
is invisible in the physics — the model simply fits better. So the architecture
makes these conditions **mechanically checkable**, and fails the run rather than
reporting a result it cannot vouch for.

## Layers

```
                        ┌─────────────────────────┐
  configurations ──────▶│  Kernel (run/driver)    │──▶ RunResult
                        │  registries + schedule  │    manifest, gates,
                        │  + graph + budget       │    conservation reports
                        └───────────┬─────────────┘
                                    │
                     ┌──────────────▼──────────────┐
                     │  Coupler (coupling/coupler) │
                     │  multirate schedule,        │
                     │  tendency routing, ledger   │
                     └──────────────┬──────────────┘
                                    │
        ┌───────────────┬───────────┴────────┬──────────────────┐
        ▼               ▼                    ▼                  ▼
   Module contracts  Transfers +        Dependency graph   State registry
   (contracts/)      conservation       (coupling/graph)   (registries/states)
                     (coupling/)
        │                                                        │
        └──────────────────────┬─────────────────────────────────┘
                               ▼
                    Core data model (core/)
              dimensions, units, fields, versions, hashing
```

Nothing in the upper layers can bypass the lower ones. A module has no method
that writes another module's state; the only path is a labeled tendency the
owner commits.

## Stage 2: the radiating column

```
solar ─▶ orbit ─▶ radiation ─▶ atmosphere ─▶ surface ─▶ space
```

| Module | Status | Owns | Role |
| --- | --- | --- | --- |
| `radiation` | ACTIVE | — | solves the whole column, hands out labeled heating |
| `atmosphere` | ACTIVE | `atmosphere_layer_temperature` | hydrostatic layers, commits what radiation contributes |
| `surface` | ACTIVE | `surface_temperature` | publishes optical properties, consumes its fluxes |

### One spectrum, not two

The solver has no shortwave/longwave split. There is radiation in band *b* with
two sources — sunlight arriving at the top, and thermal emission from the layers
and the surface — and a greenhouse configuration is an absorption coefficient
that is small where the solar spectrum lives and large where the terrestrial one
does. The SW/LW diagnostics RAD-002 and CLD-001 require are recovered by
partitioning bands at a declared 4 µm cutoff *after* the solve.

That choice has a consequence the tests measure rather than hide: the two
spectra overlap. About 0.9% of a 5772 K spectrum lies above 4 µm and about 0.07%
of a 288 K spectrum lies below it. Any comparison against a textbook identity
that assumes a clean split is therefore accurate only to that overlap, and the
Stage 2 tests quantify the term instead of widening a tolerance until they pass.

### Why radiation runs before the state owners

A radiation solver needs current temperatures. If it were ordered after the
owners it would see the previous step's values, and the lag would be a silent
first-order error in the coupling.

`Module.emit_state()` solves this: it publishes committed state at
initialisation and again after every `advance_state`, so the bus always carries
state as of the start of the next step. Fluxes stay in `emit_outputs()`,
published mid-step in splitting order, because a flux must declare an averaging
interval (DAT-006) and a state need not.

### Closure by construction

Layer heating is computed as the divergence of the interface fluxes. Whatever
leaves one interface and does not arrive at the next is exactly what heated the
layer between them, so the column's radiative budget closes identically rather
than approximately. This is deliberate: it means a conservation failure in a
coupled run points at the coupling or at another module, never at arithmetic
inside the solver. `ColumnRadiationSolution.energy_residual` measures it anyway,
because "by construction" is a claim about the code and the test is the evidence.

### Kirchhoff caught a real inconsistency

The column derives surface emissivity from albedo band by band. Running the
transparent configuration with a genuinely grey albedo of 0.3 gives 278 K, not
the textbook 254.6 K — because a grey surface reflecting 30% has an infrared
emissivity of 0.7, and Stage 1's pairing of albedo 0.3 with unit emissivity is
only self-consistent if the reflectance is confined to short wavelengths.

Both readings are legitimate configurations. What the framework will not do is
let the difference stay implicit: `shortwave_albedo` names the shortwave-only
choice, a plain float gives the grey one, and the run manifest records which was
used band by band.

### The analytic gate, and where it lives

For one layer transparent to sunlight and grey in the infrared, the steady state
is `T_s⁴ = 2 T_e⁴ / (2 − ε)` and `T_a⁴ = T_s⁴ / 2`. The rigorous test evaluates
the solver *at* those temperatures and requires the residuals to vanish — no
coupler, no time stepping, no convergence question, so a failure can only be the
solver. The coupled configuration is checked separately and more loosely, because
it violates two of the identity's assumptions through the spectral overlap above.

### Spectral convergence has two regimes

A grey atmosphere is **exactly** band-count invariant: every band shares one
emissivity and emission is normalised to σT⁴ whatever the grid, so the transfer
is linear in band-summed quantities. Agreement to the last digit across 8 to 256
bands is the strongest available statement that the spectral machinery adds no
bias of its own.

A spectrally structured atmosphere is not, and must not be. Its coefficient is
discontinuous at the window edges, so a coarse grid genuinely cannot resolve it;
refinement converges to within 0.1 W m⁻² between 256 and 512 bands, while 16
bands is wrong by more than 10 W m⁻². `window_absorption` uses each band's
*overlap fraction* with the window rather than sampling at the band centre —
centre sampling makes the coefficient jump as edges move and produces a sawtooth
that looks like non-convergence of the solver but is an artifact of how the
coefficient was built.

## The Stage 1 configuration

```
solar ──▶ orbit ──▶ surface ──▶ space
```

| Module | Status | Owns | Role |
| --- | --- | --- | --- |
| `solar` | PRESCRIBED | — | irradiance at 1 au on a declared spectral grid |
| `orbit` | ACTIVE | — | Kepler distance, spherical-mean geometry, **the top-of-atmosphere boundary** |
| `surface` | ACTIVE | `surface_temperature` | spectral absorption, grey emission, slab energy balance |
| `space` | ACTIVE | `space_cumulative_energy_received` | accounting sink at the system boundary |

### Why the orbital module owns the boundary

Charter §9.1 says the solar module "provides external energy entering the
Earth-system boundary. It does not itself add absorbed energy to Earth modules."
Something must still debit an account when energy crosses into the system, or
the top-level budget cannot close. The orbital module — which is where the
distance and geometry are applied, and therefore where the incident flux is
finally determined — declares that inflow as an external physical source and
exports exactly the same amount to the surface. Its own budget balances to zero,
and the system-level identity becomes checkable:

```
ΔE_surface + ΔE_space  =  solar input across the boundary
```

That is `test_system_energy_is_conserved_across_the_boundary`, and it holds to
1e-12 relative.

### Why space is a module rather than a number

The surface could simply subtract what it emits. Then the surface would be the
only witness to its own losses, and a sign error would be undetectable. Instead
the boundary sink recomputes its side of each transfer from the flux fields the
surface *published*, and the auditor compares the two records. A mismatch is
reported as a coupling residual (CON-004/CON-005) rather than absorbed.

### Where the spectral consistency check lands

The surface imports the broadband incident energy the orbital module sent, but
computes absorption and reflection from the **spectral** field on its own grid.
The two must agree, so the surface's conservation residual is

```
residual = (spectral_incident − broadband_incident) × dt
```

which is zero to round-off when the spectral and broadband representations are
consistent, and non-zero the moment they are not. The consistency of the two
representations is therefore enforced by the energy budget rather than asserted.

## Multirate coupling

Modules declare a `step_multiple`: how many master steps pass between their
invocations. The surface runs every hour while the solar source and orbital
geometry are exchanged daily. Between exchanges the bus holds the last published
values — a zero-order hold, declared in the schedule and reported once per
coupling interval with its maximum staleness (TIM-004).

**Reconciliation happens at synchronization points**, which fall at the least
common multiple of the step multiples. This was not the first design. The first
reconciled every master step, and every multirate run failed: the daily module
had written one transfer covering 24 hours while the hourly module had written
one covering the first hour, and the auditor correctly reported a 23-hour
mismatch that was entirely an artifact of the schedule. Charter §8.4 says
reconciliation belongs at the end of the coupling interval, and it is right.

Two consequences follow:

- **Transfer pairing excludes the time bounds.** What must match between a
  sender's and a receiver's record is the transaction's identity — parties,
  quantity, species, channel — not its granularity. The interval is checked
  separately, against the reconciliation window.
- **Double counting is detected by overlap.** Twenty-four records tiling a day
  end to end are one exchange reported finely. Two records whose intervals
  overlap are the same energy counted twice. Counting records would flag the
  first case and miss nothing about the second.

A run that stops mid-interval leaves transfers unaudited, and says so in its
warnings rather than reporting a budget it did not close.

### The `channel` field

Reflected shortwave and emitted longwave both flow from surface to space, as
energy, over the same interval. Without a discriminator they look like one
transaction recorded twice. `channel` names the physical pathway and
participates in the pairing key.

## Conservation accounting

For each module and conserved quantity, over each coupling interval:

```
residual = ΔQ_stored − ( imports − exports
                       + physical_sources − physical_sinks
                       + numerical_sources − numerical_sinks )
```

Numerical sources and sinks — clipping, limiters, filters, stabilisation — are
tracked in their own lines (CON-008). A budget that closes only because a
limiter injected energy is not a closed budget, and the separation is what makes
that visible.

The residual is judged **relative to the largest term in the budget**, not to
the net change. A module whose imports and exports nearly cancel has a near-zero
stored change and can still be wrong by a percent of the flows passing through
it.

CON-010 is implemented as `ConservationReport.can_support_effect`: a run whose
unexplained residual is comparable to the effect being claimed cannot decide
that claim, whatever the trend looks like.

## Information provenance

The dataset registry tracks evidentiary flow separately from physical flow. Each
use of a dataset declares a role for a specific claim and, optionally, specific
variables. Contamination is then a query rather than a judgement call:

- **calibration → independent validation** is refused outright (VAL-006);
- **prescribed boundary → evaluated against** is flagged as partly definitional;
- **disjoint variables of the same dataset** are not contamination — using a
  reanalysis' winds as a boundary condition does not spend its radiation record.

Withheld data goes further. A `WithheldDataset` holds a loader, not values, and
raises on any read until a `FreezeRecord` exists pinning the source revision,
configuration hash, parameter hash and criteria hash. The evolution path
therefore *cannot* read it, whatever the author intended — which is the only way
INF-008 means anything, since the person most likely to peek wrote the model.

Claim criteria hash over exactly the fields HYP-009 forbids changing after
freeze. Moving a threshold afterwards raises rather than silently widening the
success condition.

## Numerical validity

Two gates, catching different failures:

- **Impossible states.** An unphysical value (negative temperature, emissivity
  outside [0,1]) raises where it is produced. The run fails and reports why.
- **Time-step convergence** (`evaluation/convergence.py`). The harder case is a
  configuration just past its stability limit: every individual number is
  physically plausible, the trajectory oscillates instead of relaxing, and
  nothing about a single run reveals it. Refining the step does. The gate judges
  *the run's own step* — whether its answer agrees with the refined limit — and
  reports an observed order of convergence and a Richardson estimate alongside
  the verdict.

Both are exercised in `tests/acceptance/test_failure_detection.py`, together
with a deliberately leaky surface (ACC-006) and a deliberately double-counting
one.

## Run audit

`RunManifest` records what the charter's §16 asks for and hashes the subset that
determines the numbers — modules, implementations, equation bindings, state
ownership, **module parameters**, initial and boundary conditions, resolutions,
coupling schedule, solver settings and seeds — while excluding wall-clock time,
host and output checksums.

Module parameters reaching the hash matters more than it sounds. In an early
version modules published no parameters, so two runs differing only in surface
albedo hashed identically: the lineage record would have shown two runs as the
same configuration with different answers. `Module.parameters()` now exists for
exactly this, and `test_module_parameters_reach_the_configuration_hash` guards
it.

Writing a manifest over a *different* manifest at the same path raises: retained
evidence is immutable (AUD-005, A8.4), and a correction creates a child run
linked to the superseded record.

## Extending the framework

To add a module:

1. Register its equations (`registries/equations.py`) with category, validity
   domain, limiting cases and at least one correctness test. A calibrated
   parameter inside a `FUNDAMENTAL_LAW` is refused.
2. Register the states it will own, with units, authorised readers and
   authorised tendency providers.
3. Subclass `Module` with a complete `ModuleContract`, including a cost
   declaration and at least one null or limiting-case test.
4. Implement `compute_tendencies` (labeled contributions), `advance_state`
   (owned state only, retaining a `StateBudget`), `emit_outputs`,
   `report_transfers` (your own side only) and `report_conservation`.
5. Declare the new dependency edges, and any feedback cycle they close.
6. Override `parameters()` with everything that determines your numbers.

The configuration validator will refuse the module if any of that is missing.
That is the intent: an incomplete declaration is caught before initialisation
rather than discovered in a result.

## Known limitations of the current kernel

- Zero-dimensional global mean only; no diurnal, seasonal or latitudinal
  structure, and the geometry module says so in its validity domain.
- The slab has one temperature and therefore one response timescale *by
  construction*. The charter's timescale philosophy expects a real system to
  have several; the Stage 5 work is where that assumption gets tested rather
  than assumed.
- Uncertainty is declared and marked UNRESOLVED, not propagated. UNC-005's
  ensemble and Monte Carlo methods have interfaces but no implementations.
- The column has no dynamics, convection or water, so it reaches *pure radiative
  equilibrium* — which is convectively unstable in the real lower troposphere.
  That is a known wrong answer for the lapse rate, registered as a known failure
  mode of `EQ-LAYER-ENERGY-BALANCE` rather than presented as a result.
- No scattering anywhere. A scattering atmosphere is outside scope rather than
  approximated badly, and the radiation contract says so.
- Absorption coefficients are free parameters of the reduced model, not measured
  properties of any gas. Nothing in the framework yet connects them to a
  molecular spectroscopic database.
- One diffusivity factor serves both the solar beam and the diffuse thermal
  field, which is a parameterization with a declared few-percent accuracy.
- Nothing in Stage 3 onward (moist thermodynamics, clouds, ocean transport) is
  exercised.
