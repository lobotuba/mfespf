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

Version 0.6.0 — the charter's validation ladder (§19) through **Stage 6**:
the accounting kernel, the airless radiating body, the radiating atmospheric
column, atmospheric thermodynamics with convection and water, clouds, thermal
storage with measured response timescales, and historical replay of a declared
record — on top of the full Version 1.0 interface spine.

| Ladder stage | State |
| --- | --- |
| 0 — units, coordinates, transfers, conservation accounting, restart | implemented |
| 1 — airless radiating body: geometry, absorbed radiation, thermal emission, analytic equilibrium, transient | implemented |
| 2 — atmospheric radiation: transparent, absorbing and emitting, multilayer, variable spectrum | implemented |
| 3 — atmospheric thermodynamics: hydrostatic balance, dry convection, water vapour, phase-change energy | implemented |
| 4 — clouds: no-cloud control, prescribed and time-varying cloud state, separate SW and LW effects | implemented |
| 5 — thermal storage: surface heat capacity, two-reservoir ocean, step/pulse/impulse/periodic forcing, measured timescales | implemented |
| 6 — historical replay: declared inputs, frozen equations, withheld records, whole-run budget | implemented, on a **synthetic** record — no observational dataset is registered |
| 7 — hindcast and registered prediction | machinery present (freeze records), not exercised |

189 mandatory requirements are tracked; 147 currently have a bound acceptance
test. All eight KER-012 reference configurations exist and run. See
[`docs/traceability.md`](docs/traceability.md) — the uncovered requirements are
listed there rather than left implicit.

**Nothing in this framework has yet been compared with an observation.** The
Stage 6 replay runs against a record the framework generates itself, and the
provenance machinery refuses to let that support any empirical claim. That
refusal is the deliverable; the physics is Stage 1's.

## Install and run

```bash
python -m pip install -e ".[dev]"     # Python 3.11+
python -m pytest                       # 305 tests, ~2 min on a workstation
mfespf run airless --steps 4800        # Stage 1: 200 simulated days
mfespf run window --steps 6000 --heat-capacity 1e6            # Stage 2 column
mfespf run moist  --steps 4000 --heat-capacity 1e6            # Stage 3 moist column
mfespf run cloud  --steps 3000 --heat-capacity 1e6            # Stage 4 cloudy column
mfespf run ocean  --steps 16000 --master-dt 864000 --steps-per-solar-exchange 1 \
                  --forcing step                              # Stage 5 step response
mfespf run replay --steps 500                                 # Stage 6 historical replay
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
`radiative-convective`, `moist` (Stage 3); `no-cloud`, `cloud` (Stage 4);
`surface-storage`, `ocean`, `column-ocean` (Stage 5); `replay` (Stage 6).
Add `--manifest run.json` to write the run record. The Stage 5 configurations
take `--forcing step|pulse|impulse|periodic|none` and report the response
timescales they measure; `--reconstruct` adds the whole-run energy budget to any
configuration.

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

**Stage 4** adds clouds, and with them scattering. A cloud's shortwave effect
*is* reflection, so the absorption-only solver had to become a two-stream one
first. The reference run reports the two effects separately, as the charter
requires, before any net:

```
$ mfespf run cloud --steps 3000 --heat-capacity 1e6
cloud radiative effect  shortwave: -55.8632 W m^-2   longwave: +48.1483 W m^-2   net: -7.7149 W m^-2
outgoing longwave: 184.6212 W m^-2   reflected shortwave: 156.9432 W m^-2
```

Against the no-cloud control the cloud raises reflected sunlight from 101 to
157 W m⁻² and lowers outgoing longwave from 240 to 185. The signs are the
declared ones; the magnitudes are not tuned to anything.

**Stage 5** adds a second place to put heat, and with it a second timescale.
Every configuration up to here stored energy in exactly one reservoir and
carried the same declared limitation: one temperature and one relaxation time by
construction. The two-reservoir ocean removes it, and the removal is measured
rather than asserted:

```
$ mfespf run surface-storage --steps 16000 --master-dt 864000 --forcing step
  1-mode relative residual: 0.0004138          # one reservoir, one exponential
  radiative_relaxation_of_the_mixed_layer: tau = 1.702 yr

$ mfespf run ocean --steps 16000 --master-dt 864000 --forcing step
  1-mode relative residual: 0.5436             # the same fit now fails
  2-mode relative residual: 0.0002401          # by a factor of 2260
  radiative_relaxation_of_the_mixed_layer: tau = 1.433 yr
  deep_reservoir_heat_uptake:              tau = 262.9 yr
```

The two measured timescales match the analytic eigenvalues of the linearised
two-reservoir system — 1.455 yr and 263.7 yr — to 1.5% and 0.3%, the rest being
forward-Euler error and the fact that the emission is quartic rather than linear
in temperature. The two configurations are the same builder with one flag
changed, so the second timescale is attributable to the reservoir and to nothing
else.

**Stage 6** replays a declared record. The physics is Stage 1's — an airless
body under a prescribed irradiance — and everything new is around it: where the
boundary condition came from, what was done to it on the way in, whether the
equations were still the frozen ones when the evaluation ran, and what the run
is allowed to be called afterwards.

```
$ mfespf run replay --steps 500
the record supports 479 master steps at 2.592e+06 s; running 479 instead of 500

dataset SYNTHETIC-TSI-001@1.0.0 (synthetic) from GENERATED by mfespf.data.synthetic.
  total_solar_irradiance [W m^-2] interval_mean, 480 samples over [0, 1.24157e+09] s
  transformation: 480 -> 480 samples by interval_mean/piecewise_constant;
                  variance retained 0.9762

energy reconstruction over [0, 1.24157e+09] s (479 coupling intervals)
  residual    -0.00012207031 (relative 2.890e-16)
  as a mean imbalance: -9.832e-14 W m^-2

  operating mode: prescribed_history; inputs taken from the compared period
```

The record is **synthetic**, which the framework treats as a fact with
consequences rather than a caveat in prose: the evaluator caps what this stage
may conclude at `EMPIRICALLY_UNRESOLVED` however well the comparison goes, and
the run refuses to be described as a hindcast or a prediction. Registering an
observational record in its place changes the provenance class and nothing else.

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
| Clouds | [`physics/clouds.py`](src/mfespf/physics/clouds.py) | a separate variable set, never a global albedo parameter |
| Ocean | [`physics/ocean.py`](src/mfespf/physics/ocean.py) | a second reservoir in its own module, so the exchange is an auditable transfer |
| Forcing experiments | [`experiments/forcing.py`](src/mfespf/experiments/forcing.py) | a forcing delivers its declared energy at any time step |
| Declared experiments | [`experiments/protocol.py`](src/mfespf/experiments/protocol.py) | the threshold is hashed before the run; editing it invalidates the evaluation |
| Timescales | [`evaluation/timescales.py`](src/mfespf/evaluation/timescales.py) | measured per mechanism; one universal climate index is refused |
| Historical records | [`data/records.py`](src/mfespf/data/records.py) | resampling reports the variability it removed, and replays from its own record |
| Data adapter | [`data/adapter.py`](src/mfespf/data/adapter.py) | a run cannot read past a record's coverage, or a variable the record does not declare |
| Operating modes | [`run/modes.py`](src/mfespf/run/modes.py) | a reconstruction cannot be reported as a prediction |
| Whole-run budget | [`evaluation/reconstruction.py`](src/mfespf/evaluation/reconstruction.py) | a drift too small to see per interval is visible over the run |

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

**Clouds forced the solver to scatter.** A cloud cools by reflecting sunlight;
represented as absorption it would warm the column instead — a wrong response
sign, which the charter treats as refutation rather than inaccuracy. The
two-stream layer coefficients reduce *exactly* to `t = exp(−Dτ)`, `r = 0` when
nothing scatters, so every earlier result is unchanged, and give `r + t = 1`,
`e = 0` for conservative scattering.

**A prescribed cloud's water comes from nowhere, and the budget says so.** Its
condensate is not drawn from the water cycle, so the module books the change
against a source named `prescribed_cloud_water_source`. The fiat appears in the
conservation report rather than being left for a reader to notice.

**The deep ocean is a module, not a second temperature inside the surface.** It
would have been fewer lines to give the surface two reservoirs. Then the heat
moving between them would have been an internal number no audit could see, the
reservoir could not have been disabled as a single edge, and the same module
could not have been dropped under both the airless slab and the radiating column
— which it now is, unchanged. As written, the exchange is a transfer both parties
record and the ledger reconciles.

**A forcing declares energy, not an instantaneous value.** Each forcing
integrates itself exactly over the step about to be taken, so a pulse whose
edges fall inside a step still delivers exactly what it declared. Sampling at
the step start instead would make the delivered energy a function of the time
step, and a convergence study run against it would be measuring the forcing.
An impulse takes this furthest: in discrete time the narrowest perturbation is
one step wide, so `ImpulseForcing` fixes the *energy* and lets the flux follow
from the step. The impulse response is then identical across time steps to
nine digits, and the part that is not resolved is named in the declaration.

**A timescale longer than the record is marked as extrapolated.** A two-mode fit
will happily return a 124-year mode from a 44-year run, and the number looks
exactly like a measured one. It is an estimate under an assumed model, and
`ExponentialFit.constrained` says which modes the record was long enough to
have measured — the same separation between operating modes that TIM-006 asks
for between runs, applied inside a single analysis.

**A second timescale is claimed from a residual, not from a second mode.**
Granting any fit an extra mode always lowers its residual. The question is by
how much: on the single-reservoir control a two-mode fit buys a factor of five
and returns two near-degenerate timescales, and on the two-reservoir run it buys
a factor of two thousand and returns two separated ones. The framework reports
two timescales only when the residual justifies it.

**A record says what it is, not only what it was used for.** A dataset's *role*
says the irradiance was supplied to the run; its *provenance class* says the
numbers were generated. Neither implies the other, and a reader given only the
first has no way to tell that agreement with the record means nothing about the
Sun. The registry refuses an undeclared provenance in any evidence role, and the
evaluator caps a claim tested against a generated record at
`EMPIRICALLY_UNRESOLVED` — the metric is still reported, because the numerical
check is real; what is withheld is the empirical conclusion.

**A resampling reports the variability it removed.** Charter §11.4 asks for
aggregation, interpolation, uncertainty *and lost frequency content*. The fourth
is the one that decides whether a comparison means anything: a record averaged
onto a coarse model clock has had variance taken out of it, and a model that
reproduces the smoothed version has reproduced less than it appears to. The
report says which band of periods is gone and what fraction of the variance went
with it.

**Frozen equations means the whole set, not the ones a claim names.** A claim
frozen against the three equations it cites says nothing about a fourth that was
edited before the evaluation ran. `FreezeRecord` pins one hash over the entire
registry, and the evaluator refuses when it has moved.

**Running past the end of a record fails.** Holding the last value would produce
a longer run that looks the same and is partly fabricated. The steps taken before
the failure are kept — the manifest marks the run unusable, and a failed run
whose trajectory was discarded cannot be diagnosed.

**The framework refuses to state a single climate timescale.** `TimescaleReport`
carries observable, subsystem, mechanism and forcing regime with every
measurement, and `single_index()` raises when the measurements span more than a
declared factor. An average of a 1.4-year adjustment and a 263-year one
describes no mechanism in the system, and a reader given it has no way to know
that.

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
