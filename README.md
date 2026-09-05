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

Version 0.11.0 — the charter's validation ladder (§19) through **Stage 6**:
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
| 6 — historical replay: declared inputs, frozen equations, withheld records, whole-run budget | implemented, and running on real irradiance, temperature and ocean-heat records |
| 7 — hindcast and registered prediction | machinery present (freeze records), not exercised |

189 mandatory requirements are tracked; 147 currently have a bound acceptance
test. All eight KER-012 reference configurations exist and run. See
[`docs/traceability.md`](docs/traceability.md) — the uncovered requirements are
listed there rather than left implicit.

Nine real global-mean records are declared in
[`data/sources.yaml`](data/sources.yaml) and load through the framework's
provenance path. **Four claims have been registered against them and all four have
been refuted** — see below.

## Install and run

```bash
python -m pip install -e ".[dev]"     # Python 3.11+
python -m pytest                       # 382 tests, ~4 min on a workstation
python tools/fetch_data.py             # optional: fetch the declared datasets
mfespf run airless --steps 4800        # Stage 1: 200 simulated days
mfespf run window --steps 6000 --heat-capacity 1e6            # Stage 2 column
mfespf run moist  --steps 4000 --heat-capacity 1e6            # Stage 3 moist column
mfespf run cloud  --steps 3000 --heat-capacity 1e6            # Stage 4 cloudy column
mfespf run ocean  --steps 16000 --master-dt 864000 --steps-per-solar-exchange 1 \
                  --forcing step                              # Stage 5 step response
mfespf run replay --steps 500                                 # Stage 6, synthetic record
mfespf run replay-observed --steps 2000                       # Stage 6, real NRLTSI2 record
mfespf experiment solar-only                                  # the first registered claim
mfespf experiment solar-and-co2                               # with a greenhouse pathway
mfespf experiment solar-co2-ocean                             # and with deep-ocean heat uptake
mfespf experiment solar-co2-volcanic                          # and with volcanic aerosol
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
`surface-storage`, `ocean`, `column-ocean` (Stage 5); `replay`, `replay-observed`
(Stage 6).
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

That record is **synthetic**, which the framework treats as a fact with
consequences rather than a caveat in prose: the evaluator caps what a comparison
against it may conclude at `EMPIRICALLY_UNRESOLVED` however well the metric
passes, and the run refuses to be described as a hindcast or a prediction.

`mfespf run replay-observed` is the same configuration driven by the real
NRLTSI2 irradiance reconstruction — 1727 monthly steps from 1882, budget closing
to 4.7e-16 over the whole replay. Only the record changed, and because its
provenance class is `PROXY` rather than `SYNTHETIC`, the evaluator's cap lifts by
itself. That was the design claim Stage 6 made with a generated record, and it is
now tested rather than asserted.

The operating-mode rule does *not* lift with it. A replay is still a run handed
the irradiance of the period it covers, so it still cannot be reported as a
hindcast. Two guards against two different mistakes, and a real record relaxes
only one.

## Four registered claims, and what separated them

Carbon dioxide enters through the **absorption coefficient**, so the radiative
transfer computes its effect rather than being handed one. Four runs differ in
exactly one driver each, and the control with both drivers held is flat to
machine precision:

| configuration | range K | rms K | variance explained | status |
| --- | ---: | ---: | ---: | --- |
| solar only | 0.065 | 0.400 | +0.034 | refuted |
| solar + CO₂ | 0.830 | **0.222** | **+0.709** | refuted |
| solar + CO₂ + deep ocean | 0.717 | 0.243 | +0.646 | refuted |
| *observed* | *2.397* | | | |

All fail the same threshold, 0.1137 K, taken from the record's own
uncertainty. What changes is how far, and what the residual points at.

Adding volcanic aerosol is the fourth, and it makes the fit **worse**. Over its
own period (1882–2012, where the aerosol record ends) the variance explained
falls from 0.602 with the aerosol held to 0.534 with it varying:

```
$ mfespf experiment solar-co2-volcanic
run                    range K     rms K  var. expl.
everything varying      1.2717    0.2399      0.5342
volcanic held           0.5694    0.2118      0.6017
solar only              0.0583    0.3372      0.0402
CO2 only                1.2551    0.2464      0.5179
neither (control)       0.0000    0.3459      0.0000
observed                1.9379
```

That is the most useful thing the framework has produced, because a short pulse
of independently known amplitude measures the *response* on its own. A
stratospheric optical depth of 0.15 gives −2.87 W m⁻² here against an accepted
≈−3, with nothing calibrated to it — the record supplies the optical depth and
the two-stream solver does the rest. The cooling that followed was 0.68 K
against an observed 0.3–0.5 K, and chasing that gap found a defect in the CO₂
band. Recalibrated, the modelled Pinatubo cooling is **0.489 K**, inside the
observed range — and the volcanic configuration *still* fits worse than without
it, now because the underlying warming is too weak rather than because the
cooling is too strong.

**The band was calibrated at one opacity and run at another.** The wing decay
was fitted against a published forcing per doubling using
`co2_band_absorption`'s default background of 1e-4 m² kg⁻¹, while every
configuration ran at `GREENHOUSE_BACKGROUND = 3e-5`. A thinner background masks
less of the band, so the forcing per doubling in the runs was **6.28 W m⁻²** —
1.75× the calibrated figure and 70% above the accepted 3.7. Narrowing the band
to fix that then exposed a second defect: the spectral grid no longer resolved
the wings, and gave 3.700 W m⁻² per doubling where the converged answer was
3.601. A calibration performed there would have fitted the grid.

Both are now single constants shared by the calibration and the runs —
`INFRARED_BACKGROUND_ABSORPTION`, and the band region and resolution that are
`greenhouse_grid`'s defaults. The recalibrated band gives **3.698 W m⁻² per
doubling**, and everything below it is a prediction: 283.6 K, a 28.6 K
greenhouse effect, and a feedback parameter of **2.45 W m⁻² K⁻¹** against a
no-feedback Planck response of ≈3.2. Still over-sensitive, by about 30% rather
than by a factor of two.

**Every number in the table above changed, and one of them changed sign.** With
the inflated forcing, adding the deep ocean *improved* the fit — 0.609 to 0.759
variance explained — and this file said so, and said I had predicted the
opposite. That benefit was an artefact. The old column over-amplified, so
damping its amplitude helped; the recalibrated one already undershoots the
observed range by a factor of three, and damping it further now costs more than
the better timing gains: 0.709 falls to 0.646. The mechanism did not change. The
error it was compensating for did.

That is the sharpest thing the framework has done to its own results. A slow CO₂
ramp could hide a forcing error and a sensitivity error cancelling; a deep
reservoir and a volcanic pulse could not. The version 1.0.0 claims are
superseded rather than deleted, because charter §3.5 requires an abandoned
formulation to stay in the record when it affects interpretation — and this one
affects every conclusion drawn before it.

Adding CO₂ still takes the model from explaining **3%** of the observed variance
to **71%**, and it remains refuted at 0.222 K against a 0.1137 K threshold. The
protocol said in advance that the configuration could not settle
which way it would go, which is the only reason the answer means anything.

**The logarithm is a result, not an input.** The forcing per doubling is nearly
constant — within 3.2% across a factor of eight in concentration — because the
band core is saturated, so more gas widens the opaque region instead of
deepening it. Nothing asserts that; the transfer solver produces it and the
tests measure it. One parameter is calibrated, the wing decay scale, and it is
calibrated against a **published forcing value** and never against a temperature
record, which is what keeps HadCRUT5 available as independent validation.
`fetch_data.py` and the dataset registry keep the two roles apart.

**The base state is a prediction, and it is 4 K low.** Nothing but the forcing
per doubling was fitted, so the column's 283.6 K against an observed 288 — and
its 28.6 K greenhouse effect against 33 — are results. The shortfall is about
the size of the water-vapour greenhouse effect, which this configuration does
not have: the vapour is tracked and radiatively inert. Closing the gap by
raising the background opacity would mean calibrating against the temperature
record that serves as the validation, so it stays open and declared, in the
hypothesis record's validity domain rather than in a footnote.

**It used to be 27 K too warm, and that was a bug, not a limitation.** The CO₂
band was calibrated at one background opacity and run at another, so every run
used a forcing 1.75× the calibrated one; and once the band was narrowed to fix
that, the spectral grid turned out not to resolve it. Both are now single
constants shared by the calibration and the runs. The version 1.0.0 claims are
superseded rather than deleted — charter §3.5 keeps failed formulations in the
record when they affect interpretation.

**What the residual points at.** The configuration has no aerosols, no volcanic
forcing, no water-vapour feedback and no ocean heat uptake below a 50 m slab.
Two of those push the model warm and one pushes it cool, which is why the
outcome was not predictable from the configuration alone — the first time that
has been true here.

## The first registered claim, and its refutation

```
$ mfespf experiment solar-only
threshold:    0.113716 K   (the record's own mean 95% half-width)
measured rms: 0.400214 K
claim status: refuted

observed range   2.397 K        modelled range  0.0638 K
variance explained  +0.033      control (constant Sun)  +0.000
```

The hypothesis: *variation in total solar irradiance is sufficient to account
for the observed change in global-mean surface temperature since 1882*. It is
registered, frozen, and refuted against HadCRUT5 by a factor of 3.5 on a
threshold taken from **the record's own uncertainty**, not from anything the
model did.

The kernel has no carbon dioxide, no aerosols and no volcanic forcing, so this
outcome was expected. Recording it properly is the point:

- **The threshold came first.** `ExperimentProtocol` hashes the declaration
  before the run and refuses to evaluate against one that moved. The number is
  the record's mean 95% half-width over the compared period — and the compared
  period is the overlap of the two records' coverage, a property of the data.
- **The control produces exactly nothing.** A constant-irradiance run gives zero
  modelled range and zero variance explained. That is what caught a real defect:
  the run had been starting away from equilibrium, and its first decade of
  spin-up was showing up as a modelled range eight times larger than the actual
  solar response.
- **A refutation from a reconstruction is decisive where a pass would not be.**
  The run was handed the irradiance of the period it is judged on. Succeeding
  under those conditions earns `QUALIFIED` at best; failing under them earns
  `REFUTED`. `OperatingMode.status_for_result` enforces the asymmetry rather
  than leaving it to the reader.
- **The result carries what the experiment could not settle.** It observes one
  global-mean scalar, so it cannot separate the mechanism from greenhouse,
  aerosol or volcanic forcing, or from internal variability. That list is a
  required field of the hypothesis record, inside the frozen hash, and it is
  printed with the verdict.

What is refuted is *this configuration's sufficiency*, not the physics of solar
forcing. The claim record says so in its own limitations field.

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
| Hypotheses | [`registries/hypotheses.py`](src/mfespf/registries/hypotheses.py) | a mechanism must name its control and what its experiment cannot separate it from |
| Model/record comparison | [`evaluation/comparison.py`](src/mfespf/evaluation/comparison.py) | anomalies differenced against the record's own baseline; the model moves to the data |
| Greenhouse gases | [`physics/greenhouse.py`](src/mfespf/physics/greenhouse.py) | CO₂ enters as an absorption coefficient, never as a fitted forcing |
| Volcanic aerosol | [`physics/volcanic.py`](src/mfespf/physics/volcanic.py) | a scatterer; the record *is* the optical property, so little is constructed |
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
| Record readers | [`data/loaders.py`](src/mfespf/data/loaders.py) | four formats onto one declared epoch; the checksum is verified before the file is parsed |
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

**A mechanism enters through a constitutive relation, not as a forcing.** CO₂
changes an absorption coefficient and the transfer solver computes the
consequence — it is unchanged and unaware. Adding a fitted forcing to an energy
budget instead would make the greenhouse effect an input to the model rather
than a result of it, and the question this framework asks is whether the physics
produces the observed behaviour.

**A right forcing does not make a right response.** The volcanic forcing was
correct to within a few percent of a published value — it needed no calibration,
because the record *is* the optical depth — and the modelled cooling was still
half again too large. A short pulse of known amplitude separates forcing from
response where a slow ramp cannot, because a ramp lets a forcing error and a
sensitivity error cancel. That is how a CO₂ band calibrated at one opacity and
run at another was found, four claims after it started affecting results — so a
mechanism that made the fit worse was worth adding.

**A calibration target must not be the thing being predicted.** The CO₂ band's
one fitted parameter is calibrated against a published forcing per doubling.
Calibrating it against the temperature record would have spent that record, and
`check_contamination` would have said so at the next evaluation.

**A hypothesis must say what its experiment cannot separate it from.** HYP-008
is the field most likely to be left vague and the one that decides the
inference: a mechanism supported by an experiment that could not distinguish it
from four alternatives is supported jointly with those four. It is required, it
is inside the frozen hash — quietly shortening the list afterwards would
strengthen the claim as surely as loosening the threshold — and it is printed
with the verdict.

**A failure is decisive where a pass would only be qualified.** A run handed the
boundary conditions of the period it is judged on has been told part of the
answer, so reproducing the observations is weak evidence. Failing to reproduce
them is not weak at all. `OperatingMode.status_for_result` returns `REFUTED` for
a failure in any comparison mode, `QUALIFIED` for a pass that used target-period
inputs, and `SUPPORTED_WITHIN_SCOPE` only for a prospective one.

**An uncertainty column says what kind of uncertainty it is.** Real records
forced this. HadCRUT5 publishes a 95 percent interval, NOAA's ocean heat content
a standard error, NRLTSI2 a standard uncertainty — and combining them as though
they were the same quantity would understate some by a factor of two. Nothing in
the numbers distinguishes them, so `UncertaintyKind` travels with the array and a
series carrying uncertainties without declaring their meaning is refused.

**A temperature anomaly is in kelvin, not degrees Celsius.** HadCRUT5's column
heading says "deg C" and the value is a *difference*; a difference of one degree
Celsius is one kelvin exactly. Declaring it `degC` would make the unit registry
apply the 273.15 offset at the first conversion and be wrong by 273.15 K.

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
