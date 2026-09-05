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

## The greenhouse pathway

```
CO2 record ─▶ greenhouse ─▶ radiation ─▶ atmosphere ─▶ surface
              publishes a   rebuilds     the transfer solver is
              concentration its kappa    unchanged and unaware
```

Carbon dioxide enters through a constitutive relation on the absorption
coefficient. That is the charter's answer to where a mechanism may act (HYP-005),
and it is the difference between a model that produces the greenhouse effect and
one that is handed it: a fitted forcing added to an energy budget would make the
effect an input, and the question this framework asks is whether the physics
produces the observed behaviour.

### Where the logarithm comes from

The near-constant forcing per doubling is not written in. Absorption is strongest
at the band centre and falls exponentially into the wings *in wavenumber*, which
is how a vibration-rotation band's line strengths fall away. The optical depth at
a displacement `dv` from the centre is `tau_0 exp(-dv / D)`, so the displacement
at which the atmosphere stops being opaque is `D ln(tau_0)`, and the opaque width
grows by `2 D ln 2` per doubling whatever the starting concentration.

The transfer solver produces the consequence and the tests measure it: four
doublings spanning a factor of eight in concentration give forcings within 35
percent of each other, and the spread grows when the grid is coarsened until the
wings are no longer resolved. `greenhouse_grid` refuses fewer than eight bands
across the band for that reason -- with the wings unresolved the model gives
linear growth or complete saturation, and neither is what carbon dioxide does.

The forcing is converged in spectral resolution: 3.5786 W m^-2 per doubling at 67
grid bands against 3.5780 at 109.

### What is calibrated, and against what

One parameter: the wing decay scale, 36 cm^-1. It is calibrated against a
*published radiative-forcing value* of 3.7 W m^-2 per doubling and achieves 3.58
-- a three percent shortfall reported rather than tuned away by adding a second
knob.

The calibration target matters more than the number. Fitting it to a temperature
record would have spent that record: the dataset registry enforces INF-003 one
way and at variable granularity, so HadCRUT5 would have stopped being available
as independent validation for the very claim it was needed for. Calibrating
against a forcing value keeps the two roles apart, and the registry's
`assert_clean` is what confirms it.

The band-centre absorption is *not* fitted. It is set to make the core deeply
saturated -- an optical depth of order ten thousand through the column -- which
is a qualitative fact about the real 15 micron band. Once the core is saturated
its exact value stops mattering: the forcing depends on where the wings cross
unit optical depth, which depends on the logarithm of that number.

### The base state is wrong, and the reason is structural

This column equilibrates near 315 K rather than 288 K. The cause is not the
concentration but the band's *width*. A reduced column with eight layers and no
stratosphere has a weaker contrast between the surface and the level it radiates
from than the real atmosphere, so reproducing the accepted forcing per doubling
requires a band opaque over a wider range of wavelengths than the real one. The
forcing is right, the width is not, and the surplus width closes part of the
atmospheric window.

Two consequences, and they pull in opposite directions. The comparison is of
anomalies differenced against the model's own baseline, so a constant offset does
not enter it. But the Planck response at 315 K is about a third stronger than at
288 K, so the sensitivity is suppressed relative to what this configuration would
give at the right base state. That is in the hypothesis record's validity domain.

### Four runs, differing in one thing each

Both drivers varying, solar held, CO2 held, both held. The last produces a range
of exactly zero and explains exactly zero of the observed variance, which is what
makes a difference between any pair attributable to the driver that changed.

Getting there took two corrections, both found by the control rather than by
inspection. The replay was starting from a guessed initial state and spending its
first decades relaxing -- seventeen kelvin of it, appearing as a trend. And the
controls were holding their drivers at the record *mean* while the spin-up had
equilibrated to the record's *first* value, which left a transient belonging to
neither run in the difference between them. The control is now flat because both
were fixed.

### Putting the ocean under it

The two-reservoir ocean drops in unchanged -- it was written against whatever
owns the mixed layer, and the greenhouse column owns one. It changes no
equilibrium, because at equilibrium the reservoirs are at one temperature and
the exchange is exactly zero, so the spin-up is unaffected and the control stays
flat to machine precision. What it changes is the path.

The result was the opposite of what I expected. The deep reservoir reduces the
modelled amplitude, as it must. But the amplitude was not the larger error: the
*timing* was. Without it the column tracks the concentration record too closely
and warms too early; with it the response is delayed and back-loaded like the
observations, and the variance explained rises from 61 percent to 76 while the
range falls from 2.11 K to 1.83 against an observed 2.40.

That the protocol declared this undecidable in advance is what makes the answer
worth anything. Had the direction been written down after the run, it would have
been a description.

**One hypothesis, two versions.** The mechanism is identical, so it is one
hypothesis. What differs is the answer to HYP-008: a column with no deep
reservoir cannot distinguish its result from one in which ocean heat uptake
mattered, and a column with one can, so that entry leaves the list. The list is
inside the frozen criteria hash, so the difference has to be a new version --
which is what HYP-009 is for, and the first time the registry's versioning has
been used for its actual purpose rather than tested.

## The volcanic pathway, and the result it exposed

Volcanic aerosol was the cheapest mechanism left to add and it was expected to
help: the observed record has decadal dips the model had none of. It made the
fit worse, and that is the most useful thing the framework has produced.

### Almost nothing had to be constructed

Unlike carbon dioxide, whose absorption coefficient had to be built and
calibrated, the volcanic record *is* the optical property -- 550 nm
stratospheric aerosol optical depth. What is assumed is only how to spread it
across wavelength and where to put it vertically, and both are declared.

The two-stream solver needed no change at all. It has handled scattering since
Stage 4, and the aerosol contributes an optical depth and an asymmetry in the
same layers a cloud would. What did change is the solver's *argument names*:
they had been `cloud_*`, and once two different scatterers occupy one layer that
naming is a lie. They now describe a layer, and `combine_scatterers` adds the
contributions -- optical depths summing, asymmetry weighted by how much each one
scatters, because a thin aerosol beneath a thick cloud must not drag the layer's
asymmetry towards its own value.

One correction came from a magnitude check. An uncapped Angstrom power law puts
two or three times the reference optical depth into the blue and ultraviolet,
where there is a great deal of solar energy. Post-eruption droplets are large
compared with visible light, so their extinction efficiency has already reached
the geometric-optics limit at 550 nm and stops rising as the wavelength falls.
Capping it there is physics, not tuning, and it removed about a third of a
too-large cooling.

### The forcing is right and the response is not

A stratospheric optical depth of 0.15 -- roughly Pinatubo -- gives -2.87 W m^-2
here, against an accepted value near -3. Nothing was calibrated to that.

The modelled cooling after Pinatubo is 0.68 K, against an observed 0.3 to 0.5.
The cause is measurable rather than speculative: this column's feedback
parameter is 3.578 W m^-2 per doubling divided by a 3.945 K response, which is
0.907 W m^-2 K^-1 against an accepted figure near 1.2. The column is about
thirty percent over-sensitive to *any* forcing, and that traces back to the same
band-width compromise that put the base state 27 K too warm: a wider band raises
the emission level, which lowers the feedback parameter.

### Why a mechanism that made the fit worse was worth adding

A slow carbon dioxide ramp and an over-sensitive column are partly degenerate: a
forcing that is too small and a response that is too large produce a curve that
looks about right. Nothing in the carbon dioxide comparison could separate them,
and its 76 percent of variance explained was partly that cancellation.

A volcanic pulse breaks the degeneracy. Its amplitude is known independently,
its duration is short compared with the mixed layer's response time, and the
model's answer to it therefore measures the sensitivity almost directly. Adding
it lowered the variance explained from 0.674 to 0.571 -- and told us why.

That is what the charter means by preserving a negative result. The right
response is not to remove the volcanic pathway; it is to record that the
configuration is over-sensitive, which the earlier claims could not have shown.

## The first registered claim

```
hypothesis ─▶ claim ─▶ protocol ─▶ run ─▶ comparison ─▶ verdict
   frozen      frozen    frozen           anomaly,      status set by
                                          alignment     the operating mode
```

The framework's first experiment against a real observational record is
`mfespf experiment solar-only`, and it was registered in order to fail. The
kernel has no carbon dioxide, no aerosols and no volcanic forcing, so an
irradiance record alone cannot account for the observed warming. What the
experiment demonstrates is the machinery for saying so.

### The threshold has to come from somewhere that is not the model

An experiment whose threshold is chosen after the numbers arrive tests nothing.
Here it is the observational record's own mean 95 percent half-width over the
compared period -- 0.1137 K -- and the compared period is the overlap of the two
records' coverage, which is a property of the data and knowable before anything
runs. `ExperimentProtocol` hashes the whole declaration at `freeze()` and
refuses to evaluate against one that has moved since.

### A failure is decisive where a pass would only be qualified

This is the asymmetry the whole comparison path exists to encode, and it is in
`OperatingMode.status_for_result` rather than in a reader's judgement.

The run is handed the irradiance of the period it is compared against. Under
those conditions, reproducing the observations is weak evidence -- the model was
told part of the answer -- so a pass returns `QUALIFIED`. Failing under those
conditions is not weak at all: the mechanism was given every advantage the
experiment could offer and still could not account for what happened. A failure
returns `REFUTED` in any comparison mode. Only a prospective comparison can earn
`SUPPORTED_WITHIN_SCOPE`.

### What the control is for

HYP-007 requires a null configuration, and the reason is not bookkeeping. The
constant-irradiance control gives a modelled range of exactly zero and explains
exactly zero of the observed variance. Anything else would mean something other
than the solar variation was moving the modelled temperature -- and every number
in the forced run would carry the same contamination.

It caught one. The replay had been starting from an arbitrary initial
temperature, and its first decade of relaxation towards equilibrium was
appearing as a modelled range of 0.55 K -- eight times the actual solar response
of 0.064 K, and almost all of it initial condition. The forced run alone looked
plausible. The control is what made it obvious.

### Two transformations stand between a model and a record

Neither is performed silently, because both change what the comparison means.

*Differencing to a baseline.* An observational anomaly is a temperature minus
that record's mean over a stated window; a model producing absolute temperature
must be differenced against its own mean over the same window. Skipping it
compares an absolute temperature with a difference and is out by 250 K, which is
obvious. Using a convenient window instead of the record's produces a constant
offset that looks like a physical bias, which is not.

*Alignment.* The model moves onto the record's times and never the reverse.
Resampling an observational record onto a model's clock would let the model
choose how much of the observations survive.

Both emit a `TransformationStep`, so the chain from model output to verdict is
in the record.

### What the experiment could not settle is part of the result

HYP-008 asks a hypothesis to name the competing explanations its experiment can
and cannot distinguish it from. It is the field most likely to be left vague and
the one that decides the inference: a mechanism supported by an experiment that
could not separate it from four alternatives is supported jointly with those
four.

This experiment observes one global-mean scalar. It separates the solar
hypothesis from a null control and from a response of the opposite sign, and it
separates it from nothing else -- not from greenhouse, aerosol or volcanic
forcing, and not from internal variability, all of which are absent from the
configuration entirely. That list is inside the frozen criteria hash, because
shortening it afterwards would strengthen the claim as surely as loosening the
threshold, and it is printed with the verdict.

What was refuted is this configuration's *sufficiency*. Not the physics of solar
forcing, and the claim record's limitations field says so.

## Stage 6: historical replay

```
record ─▶ adapter ─▶ solar ─▶ orbit ─▶ surface ─▶ space
declared  resamples  plays     the physics is Stage 1's
input     and        the
          reports    record
```

The physics of this stage is the physics of Stage 1. Everything that is new sits
around it: where the boundary condition came from, what was done to it on the way
in, whether the equations were still the frozen ones when the evaluation ran,
whether the budget closes over the whole replay rather than interval by interval,
and what the run is allowed to be called afterwards.

### The record the framework ships is synthetic, and that is enforced

There is no observational dataset in this repository. Inventing numbers that look
like a solar composite would be worse than having none: a plausible record with
an instrument's name on it is the one thing a provenance system cannot protect a
reader from.

So the shipped record is generated, says so in its source field, its note and its
quality flags, and is classified `SYNTHETIC`. Three things follow, and none of
them is a comment:

- the dataset registry refuses a record with no declared provenance class in any
  evidence role;
- the withheld evaluator caps a claim tested against a generated record at
  `EMPIRICALLY_UNRESOLVED` however well the metric passes, with the metric still
  reported -- the numerical check is real, the empirical conclusion is not;
- the run's operating mode refuses to be described as a hindcast or a prediction.

What that buys is a Stage 6 that is honestly what it is: the replay pathway run
end to end and shown to work, with no claim about the Earth in it. Registering an
observational record in its place changes the provenance class and nothing else,
and the cap lifts on its own.

### Provenance class and role answer different questions

INF-002 asks the framework to distinguish dataset roles: calibration, boundary
condition, validation, withheld evaluation. That is what a dataset was *used
for*. It does not say what the dataset *is*.

A synthetic record can be correctly declared as withheld evaluation data. The
role is right, the declaration is honest, and a reader who sees only the role has
no way to tell that agreement with it means nothing about the world. So both
facts travel: a `DatasetRole` and a `ProvenanceClass`, neither implying the
other.

### The resampling reports what it removed

Charter 11.4 asks a resampling to declare aggregation, interpolation, uncertainty
and lost frequency content. The first three are bookkeeping. The fourth decides
whether a later comparison means anything.

A record averaged onto a coarser model clock has had variance taken out of it. A
model that then reproduces the resampled series has reproduced a smoother signal
than the instrument recorded, and the difference is invisible in the output. The
report therefore states which band of periods no longer exists and what fraction
of the record's variance went with it -- and converts to a `TransformationStep`
on the dataset record, so the registry alone says what was done and can be
replayed to the same numbers (DAT-008).

The sampling convention is carried for the same reason. Averaging instantaneous
samples and interpolating interval means are both wrong, both produce a plausible
array, and neither leaves a trace in it. A method the convention does not permit
is refused rather than accepted.

### Frozen equations means the set, not the citations

A claim frozen against the three equations it names says nothing about a fourth
that was added, edited or superseded before the evaluation ran. `FreezeRecord`
therefore pins one hash over the whole equation registry, and the evaluator
refuses when it has moved between the freeze and the evaluation: the model being
tested is not the model that was frozen.

### Reconstructing the budget over the whole run

The coupler checks conservation at every synchronization point and those checks
pass to round-off. This asks a question they cannot: over the entire replay, does
the change in stored energy equal the net energy that crossed the boundary?

The distinction matters because a residual of 1e-16 per interval accumulated over
fifty thousand intervals is 1e-16 *per interval*, and whether it matters against
the run depends on whether the errors are biased. A drift that leans one way
grows linearly; a random one grows as a square root. Neither is visible in a
per-interval report, and a slow drift is exactly what would first be mistaken for
a trend.

Nothing is recomputed. The coupling residuals say what crossed the boundary and
the module budgets say what is stored, so a disagreement is between the run and
itself. Which stores count as *inside* is a modelling decision and is named
rather than inferred -- counting the space boundary inside the system would make
every run conserve trivially and check nothing.

The result is reported as a mean imbalance in W m^-2 as well as a residual,
because that is the form a reader can weigh: a drift of 0.01 W m^-2 would matter,
and the replay's 1e-13 does not.

### Running past the end of a record fails

A framework that held the last value would produce a longer run that looks the
same and is partly fabricated. The recorded source requires the *whole* step
interval to be covered, not merely its start, because a step is an integration
over `[t, t + dt]` and running one whose second half the record does not reach
would extend the record by half a step under the zero-order hold.

The steps taken before the failure are kept. They are not evidence -- the
manifest marks the run unusable and that is what a consumer checks -- but a
failed run whose trajectory has been discarded cannot be diagnosed, and charter
3.5 asks for negative results to be preserved.

### Reading four formats onto one clock

The four declared records agree on nothing. HadCRUT5 writes `1850-01`, Berkeley
Earth writes percent-commented fixed columns with `NaN` where it has no value,
NOAA writes decimal years, LISIRD writes days since 1610. One parser each, and
each converts to seconds since a declared epoch of 1850-01-01. Two records
loaded through this module are then on the same clock and can be compared
without anyone doing calendar arithmetic by hand, which is the arithmetic that
silently puts a series a year out.

Interval means are timestamped at interval *midpoints*, because `resample`
builds its cells centred on the sample times. A series timestamped at interval
starts resamples half an interval early -- a small error, and one that survives
every check that does not look at a calendar.

Three things real data forced that the synthetic record had not:

- **An uncertainty column has to say what it is.** HadCRUT5 publishes a 95
  percent interval, NOAA's heat content a standard error, NRLTSI2 a standard
  uncertainty. Combining them as though they were one quantity would understate
  some by a factor of two, and nothing in the numbers distinguishes them. So
  `UncertaintyKind` travels with the array, and a series carrying uncertainties
  without one is refused.
- **An anomaly is in kelvin.** HadCRUT5's heading says "deg C" and the value is a
  difference; a difference of one degree Celsius is one kelvin exactly.
  Declaring it `degC` would make the unit registry apply the 273.15 offset at
  the first conversion.
- **The clock origin is a boundary condition.** A record covering 1882 onward is
  not at zero on a coordinate whose epoch is 1850, so a run replaying it starts
  its master clock where the record does -- and that origin is in the manifest
  and the configuration hash, because two replays over different periods of one
  record are different runs.

The checksum is verified before the file is parsed. Reading first and noticing
second would mean the values had already reached the caller by the time anyone
knew they were the wrong ones.

### What a run may be called

TIM-005 lists seven operating modes; TIM-006 adds the rule that makes the list
worth keeping, that a report may not describe one mode in the evidentiary
language of another; VAL-007 names the case that matters most.

The modes are ordered by how much information about the compared period was
available when the model ran. A run may be described as a mode of equal or lower
standing -- claiming less than you are entitled to is not a misrepresentation --
and never a higher one. Prescribed history is deliberately the lowest of the
comparison modes: a run handed the boundary conditions of the period it
reproduces has been told part of the answer, and its agreement is correspondingly
cheap.

The ordering does not rank quality. A controlled experiment is a better thing to
run than a bad hindcast. What it ranks is the strength of the claim the mode
licenses.

## Stage 5: thermal storage and response timescales

```
solar ─▶ orbit ─▶ ocean ─▶ surface ─▶ space
 forced by a     computes  applies the exchange,
 declared        the       exports the transfer
 experiment      exchange
```

### The reservoir is a module, and that is the whole design decision

Charter 9.10's fidelity path reaches "mixed-layer and deep-ocean reservoirs" at
its third rung. There were two ways to build it.

The short way is to give the surface module a second temperature. It works, it
is fewer lines, and it is wrong for this framework: the heat moving between the
two reservoirs would be an internal number, invisible to the conservation
auditor, impossible to disable as a single dependency edge under CAU-007, and
welded to whichever surface implementation grew it.

The way taken puts the deep reservoir in its own module owning one state. The
mixed layer stays with the surface. The heat between them is an ordinary
transfer, recorded independently by both parties and reconciled by the ledger
like every other exchange (charter 8.4), and the exchange can be disabled to
recover the single-reservoir configuration exactly.

The visible consequence is that the same `DeepOceanModule` runs unchanged under
`SlabSurfaceModule` and under `ColumnSurfaceModule`, which share nothing but an
interface. Neither surface knows anything about the ocean beyond one input field
and one exported transfer.

### Two reservoirs, two timescales, one prediction to check against

Linearised about a radiative equilibrium the coupled system is

```
C_s dT_s/dt = F - lambda T_s - gamma (T_s - T_d)
C_d dT_d/dt =                  gamma (T_s - T_d)
```

which is linear, so its step response is a sum of two exponentials whose
timescales are minus the reciprocal eigenvalues, and its frequency response is a
transfer function in closed form. Both are in `physics/ocean.py` to be compared
against -- never used by the modules, which integrate the nonlinear system.

That comparison is what distinguishes a coupling that is right from one that
merely produces a curve with two bends in it. The measured timescales agree with
the eigenvalues to about one percent, which is the forward-Euler error plus the
fact that the emission is quartic rather than linear in temperature.

One result is worth stating because it looks like a bug. A two-reservoir system
lags a periodic forcing *less* than a slab does, not more. The deep branch
contributes a zero to the transfer function as well as extra damping, and the
zero reduces the phase lag. A measured lag smaller than the control's is
therefore the prediction, and the test asserts it against the analytic value
rather than against intuition.

### Forcing: the energy is the declaration

Every forcing implements `mean_amplitude(t0, t1)` as an exact integral rather
than a sample, so the energy delivered over a run is what was declared at any
time step. A pulse whose edges fall inside a step delivers exactly its declared
energy; a periodic forcing integrates to exactly zero over whole periods and so
cannot warm the system by accident.

An impulse takes this to its conclusion. In discrete time the narrowest possible
perturbation is one step wide, so `ImpulseForcing` fixes the energy and lets the
flux follow from the step. The impulse response is then reproducible across time
steps -- identical to nine digits at 2, 4 and 8 day steps -- while the shape
below one step is unresolved and declared as such.

The forcing perturbs total solar irradiance rather than being injected inside
the system, so its energy enters through the top-of-atmosphere boundary and is
accounted by transfers that already existed. No fiat source term was needed and
none was added.

### Measuring a timescale, and refusing to invent one

`evaluation/timescales.py` offers two things that are not interchangeable.
`efolding_time` measures: it reports when the response first covered `1 - 1/e`
of its change, assuming nothing about its shape. `fit_exponential_modes` fits:
it assumes a sum of exponentials and returns the modes with the residual that
says how well the assumption held.

The residual is the part that matters. Granting a fit an extra mode always
lowers it, so a second timescale is reported only when the improvement is large:
on the single-reservoir control a two-mode fit buys a factor of five and returns
two near-degenerate timescales; on the two-reservoir run it buys a factor of two
thousand and returns two separated ones.

A fit will also return a timescale several times longer than the record it saw.
`ExponentialFit.constrained` marks those, and a `TimescaleMeasurement` carrying
one prints EXTRAPOLATED. A 44-year run of the reference ocean returns a
124-year slow mode; the same configuration run for 438 years returns 263 years,
which is the right answer. Without the marking the first would have been
reported in the same words as the second.

`TimescaleReport` carries observable, subsystem, mechanism and forcing regime
with every measurement, because TIM-008 is a prohibition as much as a
capability. `single_index()` exists only to refuse: it raises when the
measurements span more than a declared factor, rather than averaging a 1.4-year
adjustment and a 263-year one into a number that describes no mechanism.

### The experiment was declared before it was run

EXP-001 lists six things an experiment must declare, and adds the phrase that
gives the requirement teeth: *before final evaluation*. `ExperimentProtocol`
hashes the whole declaration at `freeze()`, and `evaluate()` refuses against a
protocol that was never frozen or whose declaration has moved since -- an
`InvalidForScientificInferenceError`, not a warning, because a threshold chosen
after seeing the result is not a weaker test but a different kind of thing.
Deviations (EXP-002) are appended to the outcome, never folded back into the
declaration, so what was promised and what happened stay separately readable.

## Stage 4: clouds

```
solar ─▶ orbit ─▶ cloud ─▶ radiation ─▶ ... ─▶ space
                  publishes  solves clear-sky and
                  a state    overcast, then mixes
```

### Scattering was a prerequisite, not an option

A cloud cools the planet by reflecting sunlight. The Stage 2 solver absorbed and
emitted but did not scatter, and a cloud represented as an absorber would have
*warmed* the column — a wrong response sign, which the charter's 15.4 taxonomy
calls a refutation rather than an inaccuracy. So the solver was replaced before
the cloud module was written.

The replacement is a two-stream layer solved by the adding method. Two limits
made it safe to put underneath already-validated work:

- **no scattering** gives `t = exp(−Dτ)`, `r = 0`, `e = 1 − t` — algebraically
  identical to the sweeps it replaced, so Stages 1–3 are unchanged;
- **conservative scattering** gives `r + t = 1` and `e = 0` to one machine
  epsilon, so a purely scattering layer neither absorbs nor emits.

Both are asserted directly, and a third test recomputes an absorption-only
column by an independent route and requires agreement to 1e-12.

Delta scaling comes first, because a two-stream that treated a droplet's forward
lobe as isotropic would make every cloud far too bright.

### Two solves, and why the cost is worth declaring

The clear-sky column is always solved; it is the reference the cloud radiative
effect is defined against. When a cloud is present an overcast column is solved
too, and the all-sky answer is their linear mix weighted by the cloud fraction —
the independent-column approximation under maximum overlap.

That doubles the radiation cost, which the module's limitations state. What it
buys is that the cloud radiative effect is a *difference of two computed
columns* rather than an estimate, and that horizontal photon transport between
the clear and cloudy parts is the only thing being neglected.

Linearity is what keeps the mix legitimate: the transfer is linear in its source
terms, so the mixed fluxes are still a consistent solution of the same budget
and the flux divergence closes identically after mixing.

### Reporting the components before the net

`CRE_sw` and `CRE_lw` have opposite signs and comparable magnitudes, so the net
is a small difference of two large numbers. A model with both components badly
wrong can produce a plausible net. SPC-010 and CLD-001 therefore require the
components to be available first, and the module publishes them as separate
fields — along with the clear-sky fluxes, so a reader can recompute the
difference rather than take it on trust.

### The honesty problem with a prescribed cloud

A prescribed cloud's condensate is not drawn from the modeled water cycle. Left
alone, that is silent non-conservation of water.

The module instead reports a water budget in which the change is booked against
a source called `prescribed_cloud_water_source`. The budget closes, the run
stays admissible, and the conservation report shows a named line telling the
reader the water was supplied by fiat. That is the same pattern CON-008 uses for
numerical sources: the way to handle a term that is not physics is to name it,
not to hide it or to refuse to account for it.

## Stage 3: thermodynamics, convection and water

```
solar ─▶ orbit ─▶ radiation ─▶ water ─▶ convection ─▶ atmosphere ─▶ surface ─▶ space
```

| Module | Owns | Role |
| --- | --- | --- |
| `water` | `atmosphere_water_vapour` | vapour, evaporation, condensation, **and the latent energy** |
| `convection` | — | surface sensible heat flux and dry convective adjustment |
| `surface` | `surface_temperature`, `surface_water_store` | as before, plus a water reservoir |

### Splitting order is physics, not bookkeeping

Water runs before convection so that convection stabilises the profile
condensation has already produced, rather than one it is about to destabilise.
Both run after radiation, so the profile each acts on is the one its predecessor
left. The order is declared in the schedule and appears in the run manifest.

### The latent energy has to live somewhere

Evaporation removes energy from the surface; condensation delivers it to the
atmosphere. Between the two events that energy is in the vapour, and if nothing
owns it the budget cannot close at either moment.

So the water module's energy store is `L_v` times the vapour it holds.
Evaporation becomes an ordinary transfer of *two* quantities from the surface —
water mass and its latent energy — condensation transfers the energy onward to
the atmosphere as sensible heat, and precipitation returns the mass carrying
none. That last part is a declared approximation: the sensible heat content of
the falling water is neglected, and the equation record says so.

This is what CON-007 means by a phase change reconciling matter *and* latent
energy: the two budgets close simultaneously, between the same modules, over the
same interval.

### Convection was inert until the surface flux existed

Convective adjustment was written first, on the expectation that pure radiative
equilibrium would be convectively unstable — the classical Manabe–Strickler
result. Measured in the Stage 2 configurations it was not. Adding the adjustment
changed the answer by 0.00 K, to the last digit, at every opacity tested.

The reason is instructive. The instability the textbooks describe lives in the
27 K jump between the surface and the air immediately above it, and a column
whose only surface-to-air coupling is radiative has no mechanism that can act on
that jump. The equal-pressure layers are also far too thick to resolve the
shallow super-adiabatic sub-layer.

Adding a bulk sensible heat flux supplied the missing coupling. With it, the
adjustment mixes three layers, the surface cools by about 5 K and the lowest
layer warms by nearly 20 K — the expected radiative-convective behaviour. The
flux is a closure with a free coefficient and no wind or stability dependence,
and the equation record carries all of that.

The general lesson is the one the charter is built around: a module that
produces no measurable change is not evidence that the process is unimportant.
It is evidence that something else is missing.

### Three tendencies, one owner

The atmospheric layer temperature now receives labeled contributions from three
modules at once — radiation, convection (twice: sensible heating and mixing) and
water. The column owns the state and commits their sum, and every label survives
into the per-state budget, so a reader can see how much of a layer's warming was
radiative, how much turbulent and how much latent. This is A1.3 exercised
properly for the first time.

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

Kept current rather than accumulated: a limitation that has since been removed
is deleted from this list, and where a stage removed one the entry says so.

- Zero-dimensional global mean only; no diurnal, seasonal or latitudinal
  structure, and the geometry module says so in its validity domain.
- Uncertainty is declared and marked UNRESOLVED, not propagated. UNC-005's
  ensemble and Monte Carlo methods have interfaces but no implementations.
- Absorption coefficients are free parameters of the reduced model, not measured
  properties of any gas. Nothing in the framework yet connects them to a
  molecular spectroscopic database.
- One diffusivity factor serves both the solar beam and the diffuse thermal
  field, which is a parameterization with a declared few-percent accuracy.
- Clouds are prescribed, so they respond to nothing. There is no cloud feedback
  of any kind, and the condensate the water module produces still precipitates
  immediately rather than becoming the cloud.
- Cloud overlap is maximum only, with one column fraction; partial overlap
  between layers is not represented.
- Gas scattering is absent, so Rayleigh scattering and its contribution to
  planetary albedo are missing. Cloud scattering *is* present as of Stage 4.
- One effective radius per phase, no size distribution, no mixed phase. Ice
  differs from liquid only through a larger radius, which is a crude analogy.
- Water vapour is radiatively inert: the absorption coefficient does not respond
  to it, so no configuration has a water-vapour feedback and no claim about one
  is in scope.
- The ocean has two reservoirs and no vertical structure, so nothing about the
  profile of heat uptake is in scope. Its exchange coefficient is calibrated
  rather than measured, and constant, so the model cannot represent a change in
  heat-uptake efficiency under sustained warming.
- The ocean has no salinity, circulation, meridional transport, sea ice,
  freshwater, momentum or chemical exchange, all of which charter 9.10 lists for
  the full-fidelity module.
- The forcing perturbation is spectrally flat, so no experiment about the
  wavelength dependence of solar forcing can be run with it.
- **Four claims have been registered against observations, and all four were
  refuted.** Solar alone explains 3 percent of the observed variance, adding
  carbon dioxide takes that to 61, and adding deep-ocean heat uptake to 76.
  Adding volcanic aerosol *lowers* it, because the column is about thirty
  percent over-sensitive and a short pulse exposes that where a slow ramp does
  not. Nothing has
  been *supported* by observation, and nothing can be until a prospective
  comparison exists -- a reconstruction pass earns QUALIFIED at most.
- The greenhouse configuration's base state is about 27 K too warm, which
  suppresses its sensitivity; the cause is the band width required to match the
  forcing in a column with no stratosphere.
- No aerosol, volcanic or land-use forcing exists in any configuration, and
  there is no water-vapour or cloud feedback in the greenhouse one. Two of those
  omissions push the model warm and one pushes it cool.
- Aerosol and volcanic forcing remain absent, and there is still no
  water-vapour or cloud feedback. Those are what the remaining residual is
  attributed to, and none of them is implemented.
- The recorded solar source is broadband, so the spectral shape of the imposed
  irradiance is assumed rather than recorded, and no claim about spectral solar
  forcing is admissible from a replay.
- Stage 7 is unexercised. The freeze records and the information-cutoff field
  exist and are tested, but no prediction has been registered before its
  evaluation period.

Removed by later stages, recorded because the reasoning still matters:

- *"The slab has one temperature and therefore one response timescale by
  construction."* Stage 5's two-reservoir ocean removes this, and the removal is
  measured rather than asserted -- see the residual comparison above.
- *"No scattering anywhere."* Stage 4 replaced the absorption-only solver with a
  two-stream one that reduces to it exactly when nothing scatters.
- *"The column reaches pure radiative equilibrium, which is convectively
  unstable."* Stage 3's convective adjustment removes this; the wrong answer
  stays registered as a known failure mode of `EQ-LAYER-ENERGY-BALANCE`, because
  a configuration can still be built without convection.
