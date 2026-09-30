# Open work

Deferred deliberately, with the reason and the measurement that motivates each.
An item earns a place here only if something has been *measured* that says it
matters; a wish without a number belongs in a discussion, not a backlog.

Items are removed when done, not ticked — git history is the record of what
happened. Negative results that closed an item stay in
[`architecture.md`](architecture.md) under charter §3.5.

Ordered by what the evidence says is worth doing first, not by area.

## Missing physics, in the order measured

These are mechanisms, each with an observable other than the surface
temperature that it has to reproduce. None is a correction.

### The headline is not unique without something carrying water aloft

The reference column's equilibrium depended on the vapour its initial state put
above the convective top: 291.94 K dry, 297.34 K saturated, with no diffusion.
Vertical eddy diffusion of vapour at Massie and Hunten's tracer minimum
(`eddy_diffusivity=0.39`, off by default) removes that. Dry and saturated
starts land 0.001 K apart, at **297.93 K, +9.93 K**
(`docs/protocols/stratospheric-water.md`). The +3.9 K quoted as the headline is
the dry-start member of the family, with an upper troposphere at exactly zero
relative humidity. Whether diffusion becomes the reference column's default is
a decision still open. Two things weigh on it:

- its stratospheric-water prediction fails, 17.1 ppmv against 4.0 (next item);
- the uniqueness verdict is fragile: 0.137 K apart at K = 3.9.

### Ozone heats the stratosphere, but three layers cannot hold an inversion

With vapour diffusing upward, the temperature fell all the way to 27 hPa
(190 K). So every layer above 130 hPa sat saturated, and the stratospheric
water was 17.1 ppmv against SWOOSH's 4.0.

**What has been added.** Ozone's solar absorption (`ozone=True`: Chou and
Suarez 1999 coefficients, a SWOOSH profile of 287 DU) on the measured solar
spectrum (`measured_solar_spectrum=True`: TSIS-1 HSRS v2). Both are off by
default.

**What it does.** Ozone takes 11 W m^-2 and warms every stratospheric layer by 7
to 9 K. But the temperature still falls to the top, so the frozen test fails:

- 60 hPa passes (203.5 K against 209.0);
- 27 hPa fails (198.6 K against 218.0);
- the cold point is the top, not between 60 and 160 hPa;
- stratospheric water is 60 ppmv
  (`docs/protocols/ozone-solar-heating.md`, attempt 2).

**Resolution is not the reason** (`docs/protocols/stratospheric-resolution.md`).
A log-pressure cap aloft (`max_log_step`, off by default) resolves the
stratosphere with up to 23 layers above 100 hPa.

- **The column converges.** Between caps 0.2 and 0.1, the temperatures differ
  by at most 0.26 K and the water by 0.5 percent. The refits move the wing
  decays by under 0.5 percent.
- **Still no inversion.** The minimum sits near 22 hPa at 194.6 K, against an
  observed 100 hPa at 205.7 K, and 27 hPa is 195 K against 218 K.
- **Refinement makes it colder, not warmer.**
- **No cap was adopted.** On the fine grids, the thin top layer carries the
  ozone above 10 hPa and reached 412 K. The validity gate rightly invalidated
  those runs.

**The model top is now above the ozone, and the deficit is the water**
(`docs/protocols/model-top.md`).

- **What the new top did.** At 1 hPa (`top_pressure=100.0`, with the cap), the
  ozone's heating lands where it happens and makes a stratopause of 268-270 K
  against about 270 observed.
- **Why nothing was adopted.** Those runs are invalid: the top layers overshoot
  to 342 K during spin-up, past water's declared 340 K.
- **What is left.** The middle stratosphere is 21 to 32 K too cold. Under fixed
  dynamical heating, replacing its water (34 to 44 ppmv) with the observed 4
  brings 60, 27 and 10 hPa within 8 K of NCEP. The whole deficit is the water.

**In the order to look at them:**

1. **Where the stratosphere's water comes from.** A global-mean column is
   saturated at its own cold point: 205.7 K at 100 hPa gives 37.6 ppmv, and
   the observed 4 needs the tropical cold point (195 K gives 7.3, and 190 K
   gives 3.5). This is regional: a tropical column whose tropopause sets the
   water entering a global stratosphere. The region framework is built.
2. **The upper stratosphere's cooling.** With realistic water it runs about 50 K
   hot at 1.5 hPa, so its CO2 cooling is too weak or its ozone heating too
   strong. The CO2 band core was never calibrated for this.
3. **Water's saturation relation above 340 K**, measured against an independent
   reference before its domain is moved, or a realistic initial profile. Both
   stop spin-up overshoots invalidating runs.
4. **Ozone's 9.6-micron band,** still declared absent.

The first attempt at this test ran with water vapour's solar absorption and
reached 319 K. That is the surface-reflector confound
`solar-vapour-absorption.md` had already measured, and it says nothing about
ozone.

### Clouds reflect from altitude: shortwave right, longwave waits on the air

The observed cloud field (`observed_cloud_field`, ISCCP D-series, Rossow and
Schiffer 1999) with water vapour's solar absorption puts the planetary albedo at
0.283 against CERES 0.291 and the shortwave cloud effect at -52 against -46,
nothing fitted. It fails the longwave cloud effect (+11.6 against +27.9) and the
surface-energy test, so it was not adopted
(`docs/protocols/observed-cloud-field.md`). The model's air at ISCCP's cloud-top
pressures is 13 to 28 K warmer than ISCCP's cloud tops, so the clouds trap little;
the cause is the over-moist, over-warm free troposphere, which is the next item.
Both mechanisms stay in the code, off by default, with their tests.

### Sinking air: the two columns need something that tells them apart

The free troposphere carries 0.74 relative humidity against an observed 0.4 to
0.5, and subsiding regions observe below 0.25 with a C-shaped profile
(Couhert, Schneider et al. 2010). The framework now holds several columns
(`run/regions.py`), and `two_column` links two reference columns by a relaxed,
spectral WTG circulation. It fails its frozen test
(`docs/protocols/spectral-wtg-circulation.md`). Two identical columns never
develop a temperature difference, so nothing rises or sinks: omega500 = 0 and
RH 0.74 at every relaxation distance. The strict balance before it diverged
near the cold top of convection (`two-column-circulation.md`).

The relaxed-WTG literature always imposes a surface contrast. Sunlight cannot
supply one globally: on the NCEP omega500 masks, rising regions receive 2.75
W m^-2 *less* than sinking ones, because the global partition mixes the tropical
overturning with the midlatitude storm tracks. The next attempt is the tropics
alone (30S-30N, sigma = 0.446), with whatever sets their sea-surface-temperature
contrast -- ocean heat transport -- as a measured input or a mechanism.

### The global-mean lapse rate is not all convective

The moist-convective closure gives 5.5 K/km; the global mean is 6.5, because
the extratropics are set by baroclinic eddies rather than convection (Stone and
Carlson 1979). Needs meridional structure; after the two columns.

## The number being quoted

### The plume column does not settle, and warms with every refinement

At the measured evaporating fraction, after 24000 hourly steps, the
self-consistent plume column is at +7.4 K (8 layers, settled), and still
warming at 12, 20, 30, 40 and 60 layers -- past +11.3, +15.8, +19.0, +20.9 and
+23.1 K, carrying 0.1 to 1.6 W m^-2 of imbalance that grows with resolution. It
is refuted at every one of them, so no conclusion is waiting on this. What is
open is whether it has a resolution limit at all: a column whose equilibration
time grows with the grid and whose temperature has not stopped climbing by 60
layers may be a closure that amplifies grid-scale structure, which the plume's
mass-flux derivative has done before.

### The convective top snaps to layer edges; two fixes have been tried

The full column's error scatters by 0.84 K across 25 to 45 layers and follows
the time-mean top of the convective mixing, because the adjustment's groups end
at layer edges. Both attempts so far were judged by a rule frozen before they
ran (`docs/protocols/`), with S <= 0.3 K to be adopted:

| attempt | S | adopted |
|---|---:|---|
| relaxation over 1 to 16 h | 1.1 to 2.2 K (worse) | no |
| fractional top at the level of neutral buoyancy | 0.47 K | no |

The fractional top got most of the way: five of six layer counts within 0.16 K,
the sixth (36 layers) 0.3 K above them. Its one known remaining discontinuity is
at onset -- a newly unstable layer mixes part-way into the one above rather than
from zero. That is the next thing to look at, with its own frozen test, and the
36-layer case is not to be dropped from S to get under the line.

### Closure choice is the largest uncertainty measured

Within the published relaxation range the error runs from +10 to +18 K, and the
fractional top moves it to +11.6. That is model-form uncertainty larger than any
parameter range declared so far (the
sensible heat coefficient's is 1.2 K) and larger than the resolution scatter.
It belongs in the headline as a statement about the closure, not as a
parameter to pick.

### Tests that still pin 8 layers

The Stage 3, 4 and 5 integration tests pass `n_layers=8`
explicitly: their assertions are about mechanisms (signs, orderings, closures)
and were measured there. The historical replay runs at 8 by declaration
(`REPLAY_N_LAYERS`), and whether its *response* is layer-converged has not been
checked -- its forcing is (3.698 to 3.701 W m^-2 per doubling from 10 to 60
layers). The plume column runs at 8 because it does not settle above it.

## Framework integrity

These are the items where the framework fails to do the thing it exists to do.
They come before the physics, because the physics is what they are supposed to
be protecting.

### A configuration refusal is reported as a failed conservation gate

A module raising `ConfigurationError` mid-run -- measured with the vapour eddy
diffusion refusing an unstable step at K = 11 m^2 s^-1 -- leaves a run with zero
steps whose manifest says "mandatory gate 'conservation_reconciliation'
failed". Nothing failed to conserve. The manifest should name the refusal and
its module, or a reader chasing a conservation bug will find none.

### `initial_relative_humidity` is validated and hashed but never applied

`WaterModule` accepts it, checks it lies in [0, 1] and puts it in
`parameters()`, and its `initialize` promises to apply it on the first step. No
code does, so the vapour starts at zero whatever it is set to, and two runs that
differ only in it hash differently while computing the same numbers. Starts
that need vapour set it through `initial_state` today, which works.

### Convection is an all-or-nothing closure and the column limit-cycles

At a physically consistent exchange coefficient the column stops settling. The
time-mean is a proper equilibrium -- mean imbalance zero to three decimals, and
timestep-converged -- but it oscillates about it with a surface span of 1.75 K
and a top-of-atmosphere span of 7.4 W m^-2.

The timestep study says this is not numerical. Halving the step twice moves the
mean surface temperature 0.15 K while the amplitude converges to a finite value
and the imbalance span slightly *grows*. It is a relaxation oscillation of the
closure: moisture accumulates, destabilises the group, the adjustment homogenises
all of it at once, the lifted vapour condenses in a burst that restabilises the
column, and it rebuilds.

The cure is a convective closure with a timescale rather than instantaneous
complete mixing -- a different scheme, not a tuning of this one. Until then,
anything read off this column must be a time-mean over at least a full cycle, and
a single step's numbers mean very little.

### Three modules declare ranges; the rest of the parameters do not

Modules now report through the ensemble machinery (`docs/architecture.md`, "A
module declares its ranges, and an ensemble is held to them"), and a run knows
when it is a member of one. What is declared is small, on purpose: the
evaporating fraction (measured), the bulk-plume parameter (asserted) and the
sensible heat coefficient (asserted), with the evaporation coefficient declared
*derived* from the last.

Still open, each for a stated reason rather than an oversight:

- **The ocean exchange coefficient has no range.** Its 0.7 W m^-2 K^-1 sits
  beside the 0.73 multi-model mean of two-layer fits to CMIP5 (Geoffroy et al.
  2013), but the spread across those models has not been looked up, and an
  asserted range there would repeat the evaporating-fraction mistake. It stays
  UNRESOLVED until someone reads the table.
- **Four of the five hand sweeps are still hand sweeps.** Only the exchange
  coefficient has been rerun through `across`. The evaporating fraction and the
  bulk-plume parameter now have declarations to sweep across, and their
  equilibrated sweeps on the plume column have not been rerun. The timestep and
  gravity are a numerical sweep and a model-form comparison respectively, and
  belong in an undeclared-setting sweep and `model_forms` rather than in a
  module's ranges.
- **Nothing declares a distribution.** Every declared range is "a range, not a
  distribution", which is honest and means no interval in this project can yet
  state a coverage. A Monte Carlo needs a parameter whose distribution was
  measured.

### Six mandatory requirements are still unbound

179 of 185 have a bound acceptance test, with 4 architectural carve-outs counted
separately rather than folded into the denominator. What is left is two kinds of
thing, and neither is a backlog item that code alone closes:

| | requirements | what it needs |
|---|---|---|
| a spatial decision | SPA-001, SPA-003, SPA-004, SPA-005, SPA-006 | see below |
| true only when the rest are | ACC-001 | everything above |

VER-008 closed when negotiation was wired in -- see `docs/architecture.md`,
"Interface negotiation existed and nothing called it".

**The spatial five need a decision before they need code.** The kernel has two
spatial representations, global scalar and vertical column, and it changes the
column's resolution in every convergence study — by rebuilding the configuration
rather than by regridding. So SPA-003 and SPA-004 (regridding preserves
integrated quantity; intensive quantities use a documented method) have nothing
to regrid, SPA-005 (resolution changes appear in run lineage) is arguably already
true through the manifest's `spatial_resolution`, and SPA-001 (the architecture
permits seven representations) is satisfied by a contract field that nothing
constrains, which is not a demonstration. They were loosely counted as carve-outs
and they are not: the kernel uses vertical structure. Deciding which of the five
are honestly satisfiable now, and by what, is the work.

## The ladder

### Stage 7 is unexercised

`FreezeRecord`, `WithheldDataset` and `WithheldEvaluator` are implemented and
tested. Nothing on this list blocks a registered prediction, and none of the
physics work does either.

The uncomfortable part is the argument for doing it *now*. Four claims are
registered and all four are refuted; the column has no working water-vapour
feedback and no cloud feedback. A prediction registered today would be one the
model is expected to fail. That is not a reason to wait — it is the closest
thing this project has to its own point, and waiting until the model is good
enough to pass is how prospective prediction quietly turns into fitting. The
ladder's terminus has never been climbed, and climbing it badly on purpose is
worth more than not climbing it.

## Convection and the vapour pathway

### The global-mean column itself

**It is a ceiling, not a bug.** The subsiding-branch mechanism works and the
parameters it needs are not obtainable: the two parameter-free limits are 0.015
and 0.725 against an observed 0.4 to 0.5, so they bracket nothing, and what sets
the real value is horizontal transport that a single column does not contain.
See [`architecture.md`](architecture.md) for the three honest responses, none of
which is a closure.

**Now the only large item.** One lever moves the base state and the sensitivity
together -- 305.18 K and 8.76 K per doubling at RH 0.725, against 288.84 K and
3.81 K at RH 0.40 -- so the warm bias and the two-to-three-fold excess
sensitivity are one error, not two and not many. The table is in
[`architecture.md`](architecture.md). Everything else in this file is small by
comparison, and the vertical-structure path that looked like an alternative is
closed.


Promoted from a caveat to the main obstacle by measurement. With every parameter
taken from the literature rather than assumed, the bulk-plume closure runs the
column 13 to 26 K too warm, and no value inside the measured range escapes it —
the entrainment floor requires an evaporating fraction *above* where the column
has already failed, so the two literature constraints do not overlap.

The closure is not wrong. At those parameters it predicts a free-tropospheric
humidity of 0.73 to 0.84, and cloud-resolving simulations of
radiative-convective equilibrium sit near 0.80. It describes a convecting column
correctly. This framework's column is a global mean, whose free troposphere is
far drier because it is dominated by subtropical subsidence and by where
descending air was last saturated — neither of which a bulk plume contains.

So this is the deepest structural assumption in the model, and it now has a
price tag: about seventeen kelvin at the central parameter value. The options
are a second column for the subsiding branch, an explicit last-saturation
treatment, or an honest restriction of the configuration's declared validity
domain to convecting conditions. The third is much the cheapest and is arguably
what the evidence already supports. Everything else in this section is
downstream of the choice.

### A radiative criterion for the convective top

Convection currently stops where the profile becomes stable with respect to the
target. That is defensible and it lands in a sensible place — 186 hPa and
12.1 km at 20 layers, which is a reasonable tropopause — but it is a *local*
criterion, and the physically correct statement in a radiative-convective column
is that convection extends to where radiation alone can sustain the profile.

The two coincide in the current column. They will not in general, and in
particular they will diverge in any configuration whose stratosphere is doing
something, so this is worth having before a column acquires ozone or any other
stratospheric heating.

**Why not the alternative.** An undilute-parcel level of neutral buoyancy was
considered and rejected on measurement rather than taste: in a column that is
neutral to an *entraining* plume, an undilute parcel is buoyant everywhere the
entraining profile holds, because its temperature excess is the accumulated
entrainment deficit and grows monotonically with height (Singh and O'Gorman
2013, eq. 4). Measured here it runs from +0.10 K at the surface to +9.04 K at
10.7 km and only turns negative at layer 17 — one layer above where the
adjustment already stops. The parcel's level of neutral buoyancy is a
*consequence* of the convective top rather than a criterion for it, which is why
Singh and O'Gorman took theirs from simulations. Implementing it would have
changed nothing.

**Cost.** The criterion needs the radiative equilibrium temperature, which
convection does not currently receive. That is a new module dependency and a
graph edge, not just a function.

### The humidity's vertical structure is closed, not pending

The general form is implemented, tested and wired behind
`plume_vertical_structure`, and it stays **off**. It does not converge -- that part is this
implementation's own -- and closing it properly needs an independent entrainment
law that the literature does not have to give. Romps takes the rates from
cloud-resolving output, and his direct measurement finds the bulk-plume rate too
small by about a factor of two. The published rates also rise with model
resolution, but that is a diagnostic artefact rather than physics: the core and
environment populations are separated by a definition that moves with the grid.
The mixing rate itself is well posed.

Supplying those profiles would be permitted -- no temperature record is involved
-- and would still be the wrong trade. This closure is here because one number
predicts the humidity and the lapse rate together, so it can fail; profiles read
off simulations can reproduce anything. Kept in the tree under charter 3.5 as a
negative result. See [`architecture.md`](architecture.md).

## Feedbacks the model does not have

### Clouds respond to nothing

Clouds are prescribed, so there is no cloud feedback of any kind, and the
condensate the water module produces precipitates immediately rather than
becoming the cloud.

What has been ruled out is clouds as the *remedy for the warm bias*. The column
carries Earth's planetary albedo, 0.30, on its surface, so the clouds' shortwave
cooling is already present under another name; the same column at a true surface
albedo of 0.12 runs at 338 K. Making the cloud explicit at matched absorbed
sunlight adds its longwave warming and leaves the column 0.9 K *warmer*, and a
feedback cannot move a base state in any case. See
[`architecture.md`](architecture.md).

What remains genuinely open is the cloud module's fidelity. Its single liquid
deck produces a shortwave effect of -65 W m^-2 against an observed -47, and a
longwave effect of +5 against an observed +27, because Earth's longwave cloud
effect comes mostly from high cold ice cloud. A configuration that reproduced
both components would be worth having on its own terms, and it is a prerequisite
for any later cloud response -- but it should be built to match the observed
components, not to move the surface temperature.

The honest order is still to settle what the column *is* before deciding what
should respond inside it, and the global-mean question above is that decision.

### Surface albedo is carrying the clouds

Filed separately because it outlives the cloud question. `shortwave_albedo`
defaults to 0.30, which is a planetary albedo standing in a surface's place. It
is a defensible idealisation and it is nowhere declared as one: the parameter is
named for the surface, documented as the surface, and worth about 33 K.

Any configuration that acquires real clouds has to drop it to something like
0.12 in the same change, or count their shortwave effect twice. That coupling
between two parameters in different modules is exactly what a validity domain
would record, which is the item at the top of this file.
