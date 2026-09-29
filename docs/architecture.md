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

One parameter: the wing decay scale, 14 cm^-1. It is calibrated against a
*published radiative-forcing value* of 3.7 W m^-2 per doubling and achieves 3.698
-- at the background opacity and on the spectral grid the configurations use,
because both of those are `physics/greenhouse.py` constants that the calibration
and the runs share.

**They were once not shared, and that was the largest defect this project has
found in itself.** ``co2_band_absorption`` defaulted to a background of
1e-4 m^2 kg^-1 and every greenhouse configuration passed 3e-5. A thinner
background masks less of the band, so the same wing decay gave 6.28 W m^-2 per
doubling in the runs -- 1.75 times the calibrated figure, 70 percent above the
accepted value, and the number this document once quoted was the one no run had
ever used. It was found by asking why a volcanic pulse of the right forcing
produced twice the right cooling.

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

### The base state is a prediction, and it is four kelvin low

This column equilibrates at 283.6 K against an observed 288, with a greenhouse
effect of 28.6 K against 33. Neither is a calibration target. The only fitted
number in the band is the wing decay, and it is fitted to a published forcing per
doubling, so the base state is what the physics gives once that is fixed.

Four to five kelvin is about the greenhouse warming water vapour supplies. The
vapour is tracked -- the water module has owned it since Stage 3 -- but it is
radiatively inert, so the shortfall is exactly where the missing mechanism is.
Closing it by raising the background opacity would be calibrating against the
temperature record, which is the one thing INF-003 forbids here, because that
record is the validation.

**An earlier version of this section said the opposite, and the reason is worth
keeping.** The column used to equilibrate near 315 K, twenty-seven kelvin *too
warm*, and this document explained it as a band too wide. The width was not the
problem. The band was calibrated with `co2_band_absorption`'s default background
opacity and every configuration passed a different one, so the runs used a
forcing 1.75 times the calibrated figure. Two things fixed it: one constant used
in both places, and a spectral grid on which the answer has stopped moving. The
band that results is 2.6 times narrower, and the base state crossed from 27 K
high to 4 K low.

### A calibration is only worth anything on the grid it was performed on

Narrowing the band made the old grid inadequate without changing anything about
it. Forty-eight bands over 6.5 to 50 microns gave 3.700 W m^-2 per doubling where
the converged answer was 3.601 -- so a calibration performed there would have
fitted a property of the grid, and the fit would have silently drifted with any
later change of resolution.

So `CO2_BAND_LOWER_M`, `CO2_BAND_UPPER_M` and `CO2_BAND_RESOLUTION` are constants
next to the wing decay, and they are `greenhouse_grid`'s defaults, so a
configuration that takes the defaults gets the calibrated forcing and one that
does not is visibly choosing otherwise. The region is 8 to 30 microns because
that is where the band is: the opaque core spans about 12 to 20 microns at the
highest concentration run here, and the wings have fallen by nineteen orders of
magnitude by 8 microns. Sixty-four bands across it agree with a 179-band
reference to 0.05 percent.

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

The protocol declared the direction undecidable in advance -- whether the
reservoir helps depends on whether the amplitude or the timing is the larger
error, and the configuration alone does not say which. That declaration is what
makes the answer worth anything, and it is also what let the answer *change*
without any special pleading.

It has been measured twice, with opposite results. Under the mis-calibrated band
the reservoir raised the variance explained from 61 percent to 76: the column
over-amplified, so damping it helped, and the delayed back-loaded response
gained more than the lost amplitude cost. Recalibrated, the same reservoir
*lowers* it from 71 percent to 65, with the range falling from 0.83 K to 0.72
against an observed 2.40. The mechanism did not change and neither did the
reasoning; the column no longer over-amplifies, so there is nothing left for the
damping to correct and only amplitude to lose.

The first answer is superseded rather than deleted, because a reader who takes
"the timing mattered more than the amplitude" from this project's history should
also learn that it was an artefact of a defect elsewhere.

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

### The pulse that found the calibration defect

A stratospheric optical depth of 0.15 -- roughly Pinatubo -- gives -2.87 W m^-2
here, against an accepted value near -3. Nothing was calibrated to that: the
record *is* the optical depth.

The response was 0.68 K of cooling against an observed 0.3 to 0.5, and that gap
is what exposed the CO2 band's calibration defect. A slow concentration ramp
could not have: an inflated forcing and an inflated sensitivity partly cancel
against the mechanisms the configuration omits, and the result still looks
plausible. A short pulse of independently known amplitude cannot be compensated
that way, so it measured the response on its own.

Recalibrated, the modelled Pinatubo cooling is **0.489 K**, at the upper edge of
the observed range, and the feedback parameter is 2.45 W m^-2 K^-1 against a
no-feedback Planck response near 3.2. The column is still over-sensitive, by
about thirty percent rather than by a factor of two.

Two earlier versions of this section got this wrong in opposite directions. The
first quoted 0.907 W m^-2 K^-1 against an accepted 1.2 and called the column
thirty percent over-sensitive; both halves were wrong, because the forcing was
the calibrated figure rather than the one the runs used, and 1.2 is the *net*
parameter including feedbacks this configuration does not have. The second
quoted 1.59 against 3.2, which was right about the reference and right about the
defect, and is now superseded by the fix.

### Adding it still makes the fit worse, for a different reason than before

Over 1882 to 2012 the variance explained falls from 0.602 with the aerosol held
to 0.534 with it varying. That was true before the recalibration too, and the
drop is nearly the same size -- but the cause is not. It used to be an
over-sensitive column turning a right forcing into too much cooling. Now the
cooling is about right and the *background trend* is too weak: with the aerosol
held the modelled range is 0.5694 K against an observed 1.9379, so adding
correctly-sized cooling spikes to an under-warming record moves it further from
the observations rather than nearer.

That is a better failure than the one it replaces, because it points somewhere
specific: at the missing water-vapour feedback rather than at the aerosol.

### Why a mechanism that made the fit worse was worth adding

A slow carbon dioxide ramp and an over-sensitive column are partly degenerate: a
forcing that is too small and a response that is too large produce a curve that
looks about right. Nothing in the carbon dioxide comparison could separate them,
and its 76 percent of variance explained was partly that cancellation.

A volcanic pulse breaks the degeneracy. Its amplitude is known independently,
its duration is short compared with the mixed layer's response time, and the
model's answer to it therefore measures the sensitivity almost directly. Adding
it lowered the variance explained from 0.674 to 0.571 -- and told us why, which
is the part that mattered: the trail from that gap led to a band calibrated at
one background opacity and run at another.

That is what the charter means by preserving a negative result. The right
response was not to remove the volcanic pathway. It was to record that the
configuration was over-sensitive, which the earlier claims could not have shown,
and then to find out why. After the fix the aerosol still lowers the variance
explained -- 0.602 to 0.534 -- but the modelled Pinatubo cooling is 0.489 K
inside an observed 0.3 to 0.5, so what the pathway now exposes is the missing
warming rather than an excess of cooling.

## The water-vapour pathway, and the two walls it hit

Water vapour is the largest greenhouse gas and the largest feedback, and until
now this framework had neither. The pathway was built in five steps, each of
which produced a measurement, and two of the five produced a wall.

### The relation is calibrated against radiation, never against temperature

`EQ-VAPOUR-BAND-ABSORPTION` gives a mass absorption coefficient across the
spectral grid: a window continuum, a rotational band below 200 cm^-1 and the
bending band at 1595 cm^-1, each falling off as an exponential wing. Two numbers in it
are fitted, and both are fitted to *published radiation observations* -- the
CERES all-sky outgoing longwave near 240 W m^-2, and a window transmittance of
0.65 at 10 micron -- never to a temperature record. That is INF-003, and it is
what makes the feedback the relation produces a prediction rather than a
restatement.

Nothing in the relation knows about temperature. The feedback is emergent, and
in isolation it measures **+2.29 W m^-2 K^-1**.

### The first wall: a feedback the column could not use

Coupled to a real column the same relation produced an amplification of
**0.99x** and a surface 12 K *colder*. The diagnosis was structural rather than
numerical. The water module puts all of the vapour in the lowest layer, because
nothing in the framework lifted it; a layer at the surface's own temperature
re-emits almost everything it absorbs, so a large opacity there traps nearly
nothing. Radiatively active but confined vapour gives a 63 W m^-2 greenhouse
effect against 125 W m^-2 for the same water spread through the column.

### The control: what the feedback would be worth if the water were where it belongs

`PrescribedHumidityModule` imposes a uniform relative humidity of 0.7 and lets
Clausius-Clapeyron decide the amount -- the Manabe-Wetherald idealisation. It
owns no state, conserves nothing and may not carry a claim; it exists to
measure the target. Against it, dT_2x is 4.394 K with the humidity free and
1.768 K with it held: amplification **2.49x**, feedback **+1.251 W m^-2 K^-1**.
That is below the radiative +2.29 because convective adjustment supplies
something like a lapse-rate feedback, and published water-vapour-plus-lapse-rate
estimates are about +1.2, which is the right neighbourhood for the wrong-looking
reason to be the right one.

### The second wall: the transport that could not start

The obvious next step was to let convection carry the water, reusing the closure
that already mixes heat so that no new parameter enters. It moved nothing, and
the reason is a loop rather than a bug. Vapour with nothing to lift it collects
in the lowest layer; radiatively active, that layer becomes an efficient
radiator; it cools far below the layer above; and the inversion decouples it
from the neutral region higher up, so the adjustment is inactive exactly where
the water is. **The missing transport sustained the stratification that blocked
the transport.** At 20 layers the potential-temperature step across the lowest
interface measured 50.23 K with the vapour radiatively active and 0.00 K with it
inert -- which locates the cause in the trapped vapour and not in the closure.

### Surface-based convection removes the blockage

The layer-only adjustment compares each layer with its neighbour, so it cannot
see the instability that actually drives convection: the one between the ground
and the air on top of it. Admitting the surface as element zero of the profile,
with its own heat capacity and its own Exner factor, fixes that without adding a
parameter -- the mixed value is still defined by enthalpy conservation, now over
the surface and the air together.

| | layer-only | surface-based | surface-based, carrying water |
|---|---:|---:|---:|
| surface temperature | 273.02 K | 272.63 K | 278.57 K |
| lowest layer | 240.69 K | 270.68 K | 276.59 K |
| potential-temperature step | 50.23 K | 0.00 K | 0.00 K |
| layers mixed | 0 | 5 | 5 |
| vapour above the lowest layer | 0.0% | 0.0% | 78.7% |

It also cost something the framework insisted on. Mixing the surface into the
column moves energy across a module boundary, and the first version balanced its
own arithmetic perfectly while failing conservation reconciliation at 3e-4. The
exchange is now published as a flux and declared as a transfer by both the
surface and the atmosphere.

### The third wall is the one still standing

With the water finally able to move, the column produces an emergent humidity
profile for the first time -- and the feedback it gives is not safely converged
in the layer grid.

| layers | 10 | 20 | 30 | 40 | 60 |
|---|---:|---:|---:|---:|---:|
| surface temperature (K) | 275.82 | 278.57 | 279.56 | 280.75 | 280.94 |
| precipitable water (kg m^-2) | 4.52 | 4.89 | 5.20 | 5.28 | 5.32 |
| amplification | 1.10x | 1.22x | **1.88x** | 1.28x | 1.30x |
| feedback (W m^-2 K^-1) | +0.214 | +0.402 | **+1.014** | +0.466 | +0.495 |

The spin-up is not the explanation: at 20 layers, 2000, 4000 and 8000 steps give
the same surface temperature and the same dT_2x to three decimals. Nor is it
simple divergence -- read only to 30 layers the sequence looks like one, which
is how it was first written down here, and 40 and 60 layers then land within two
percent of each other in the amplification and six in the feedback. What is left is a 45 percent excursion at one resolution
in an otherwise smooth sequence, which is an unexplained switch rather than a
limit, and VAL-002 says a number in that state is recorded and not quoted.

The physics behind it is visible in the humidity profile. Convective mixing is
complete rather than entraining, so it drives the top of the mixed group past
saturation; condensation clips it back; and which layer that happens in is set
by the grid. Above the group, relative humidity is exactly zero, because
evaporation is the only source of water in this framework and convection its
only transport. So the emergent feedback comes to about +0.5 W m^-2 K^-1 against
+1.251 for the prescribed-humidity control -- roughly the fraction of the column
the mixing reaches.

Two numbers say the same thing more bluntly. The column settles at 280.9 K and
5.3 kg m^-2 of precipitable water; the relation in it was calibrated on a column
carrying 23.2. **The relation is being evaluated far outside the conditions it
was fitted at** -- the same class of defect this project has already found twice,
in the carbon dioxide band's background opacity and then in this relation's own
clear-sky against all-sky pairing, appearing a third time. Closing that gap means
entrainment and detrainment, and the cheapest form of those carries a fitted
eddy diffusivity -- which is exactly the choice deferred when this pathway
started, arriving now with a measurement attached to it.

### Entrainment, and a humidity that is predicted rather than chosen

The wall above was an entrainment problem wearing a numerical disguise, so the
fix had to come from the entrainment literature rather than from the grid. The
cheapest form is a bulk plume in the limit of zero cloud buoyancy, and it was
chosen over the alternatives for one property: **the same dimensionless
parameter fixes both the humidity and the lapse rate**. A closure that predicts
two things from one number can be caught being wrong. A prescribed relative
humidity cannot.

The environment's water budget in radiative-convective equilibrium balances
subsidence drying against detrainment of saturated plume air and evaporating
condensate, and closes for the relative humidity:

```
RH = (delta + alpha gamma - alpha epsilon) / (delta + gamma - alpha epsilon)
```

with `gamma = -d ln q* / dz` following from Clausius-Clapeyron and hydrostatic
balance alone. Holding the mass flux constant forces entrainment to equal
detrainment and collapses this to `RH = (a + alpha) / (1 + a)`, while the same
`a` fixes the stratification an entraining plume leaves neutral -- which reduces
exactly to the moist adiabat at `a = 0`, checked against this framework's own
Stage 3 implementation rather than against a second copy of the algebra.

**Nothing is fitted.** `a = 0.2` is the published value for Earth's tropics,
derived from convective available potential energy rather than from any humidity
or temperature measurement, and `alpha` is a declared assumption carrying the
range 0.1 to 0.5 rather than a fitted number.

Sources: Romps (2014), *J. Climate* **27**, 7432, for the humidity; Singh and
O'Gorman (2013), *GRL* **40**, 4398, for the zero-buoyancy plume; Romps (2016),
*JAS* **73**, 3719, for the constant-mass-flux closure.

#### The grid dependence goes away, which confirms the earlier diagnosis

| layers | 10 | 20 | 30 | 40 | 60 |
|---|---:|---:|---:|---:|---:|
| surface temperature (K) | 286.94 | 288.11 | 288.51 | 288.72 | 288.93 |
| amplification, transport scheme | 1.10x | 1.22x | **1.88x** | 1.28x | 1.30x |
| amplification, plume closure | 1.86x | 1.80x | 1.92x | 1.93x | 1.86x |
| feedback (W m^-2 K^-1) | +1.038 | +0.886 | +0.984 | +0.977 | +0.894 |

A factor-of-1.7 range with a 45 percent excursion becomes scatter of about three
percent with no trend. That is the point worth keeping: the transport scheme's
feedback moved with resolution because the humidity aloft depended on which
layer a mixed group's top fell in, and a humidity that comes from a water budget
rather than from a mixing operator does not inherit that. The diagnosis was
right, and replacing the mechanism it blamed fixed the symptom it predicted.

#### Both predictions land, and the model disagrees with one of them

The over-determination is the reason this closure was worth building, so it has
to be tested rather than asserted. At `a = 0.2` it predicts a relative humidity
of **0.417**, which is inside the observed free-tropospheric band, and a lapse
rate of **6.33 K/km**, against roughly 6.5 in the standard atmosphere. Two
observables, one unfitted parameter, and both land.

The column they were measured in does not.

| | K/km |
|---|---:|
| the column's diagnosed lapse rate | 9.68 |
| the plume prediction at a = 0.2 | 6.33 |
| moist adiabat | 5.98 |
| dry adiabat | 9.77 |

Convective adjustment still targets the dry adiabat, so the column sits within
one percent of it while carrying a humidity derived from `a = 0.2`. **The model
is internally inconsistent**: its stratification corresponds to an entrainment
rate far larger than the one its own humidity closure assumes. Wiring
`EQ-ENTRAINING-LAPSE-RATE` into the adjustment is the obvious next step, and the
measurement above is what makes it obvious.

#### Making the column consistent, and what that revealed

Wiring `EQ-ENTRAINING-LAPSE-RATE` into the adjustment removes the
contradiction: the convective closure and the humidity closure now share one
parameter, and both are used.

| target | surface (K) | precipitable water | diagnosed | predicted |
|---|---:|---:|---:|---:|
| dry adiabat | 288.11 | 8.02 | 9.68 K/km | 6.33 K/km |
| entraining plume | 289.35 | 14.28 | 5.63 K/km | 5.48 K/km |

Three things came out of it, and only the first was expected.

**The lapse-rate prediction is weaker than it first appeared.** Measured against
the dry-adiabat column it looked like a hit -- 6.33 K/km against a standard
atmosphere's 6.5 -- but that was the closure evaluated on a profile it
disagreed with. Run self-consistently it settles at 5.5, about fifteen percent
shallower than observed. That is a real disagreement, and it only became visible
once the closure was made to live with its own answer.

**The two errors that were cancelling were not independent.** Fixing the
stratification fixed most of the dryness by itself: precipitable water went from
8.0 to 14.3 kg m^-2, from a factor of 2.9 below observation to 1.75, because a
shallower lapse rate warms the upper column and raises the saturation humidity
the closure fills. The surface moved 1.2 K through all of that. **The surface
temperature was the least sensitive quantity in the column**, which is exactly
why it looked right while the interior was wrong in two ways at once.

**The feedback decomposition broke, and its numbers are not quotable.** With the
dry target, freezing the vapour changed one thing. It no longer does: the
entraining lapse rate depends on the saturation humidity, so a frozen-vapour run
still gets the full lapse-rate response, and the control is no longer null. It
shows up as an absurdity rather than a subtlety -- the held run gives
dT_2x = 0.299 K, a feedback parameter of 12.4 W m^-2 K^-1 against a Planck
response of 3.2 -- so the amplification of 13x it implies measures nothing. A
null control that is not exactly null is a bug, and this one is recorded rather
than reported.

The first explanation offered for it was wrong, and the correction is worth
keeping. It was attributed to convection reaching too high -- the adjustment
mixes sixteen of twenty layers, which was read as the column having no real
tropopause. Layer sixteen is at 186 hPa and 12.1 km, which is a perfectly
ordinary tropopause; the layer count had been read without checking what
pressure it corresponded to. Convective depth had nothing to do with it, and the
broken control was the whole of the problem.

#### The control repaired, and the number it now yields

Holding the vapour has to hold everything the vapour does. It acts on this
column twice -- radiatively through the absorption, and dynamically because the
entraining lapse rate depends on the saturation humidity -- so the control now
freezes both through one switch, since holding the vapour is one physical
intervention even if it has two consequences. What stays free is the geometry:
the *rate* is frozen and not the offsets it integrates to, because a warming
column expands and that is hypsometry rather than the vapour's doing.

| | broken control | repaired |
|---|---:|---:|
| dT_2x held (K) | 0.299 | 1.663 |
| lambda held (W m^-2 K^-1) | 12.43 | 2.22 |
| amplification | 13.23x | 2.38x |
| feedback (W m^-2 K^-1) | +11.433 | +1.289 |

The repaired number is grid-stable -- +1.289 at 20 layers and +1.241 at 40 --
and it agrees to within four percent with the +1.251 the prescribed-humidity
control measures and with published water-vapour-plus-lapse-rate estimates near
+1.2. From a closure with nothing fitted, that is the best agreement this
project has produced.

It is still not a result, for the reason the closure's own uncertainty makes
plain.

| alpha | 0.1 | 0.3 | 0.5 |
|---|---:|---:|---:|
| surface temperature (K) | 284.47 | 289.35 | 295.54 |
| dT_2x free (K) | 2.843 | 3.957 | 7.010 |
| dT_2x held (K) | 1.834 | 1.663 | 1.516 |
| feedback (W m^-2 K^-1) | +0.716 | +1.289 | +1.911 |

Eleven kelvin of surface temperature and a factor of 2.7 in the feedback, across
a range nobody has narrowed. The held column is steady at 1.5 to 1.8 K
throughout, which is what a working control looks like -- the spread is the
feedback responding to the assumption, not the control wobbling.

There is a temptation here that has to be named and refused. The upper end of
the range gives dT_2x of 7 K, which is outside any defensible sensitivity, so
the model's own temperature response would narrow `alpha` considerably. Using it
that way would be calibrating against the thing being predicted, which is
exactly what INF-003 forbids and what this project has already caught itself
doing twice. `alpha` has to come from precipitation-efficiency measurements or
stay wide.

#### Looking the parameter up refuted the configuration

The evaporating condensate fraction had been *assumed* -- 0.3, asserted over 0.1
to 0.5, on the reasoning that a declared range is more honest than a fitted
number. The reasoning was sound and the range was wrong. It placed the
precipitation efficiency at 0.5 to 0.9; the literature places it at 0.2 to 0.4.
Langhans, Yeo and Romps (2015) measure 33 and 38 percent in large-eddy
simulations of cumulus congestus, using the definition this closure needs, and
report Braham (1952) at 19 percent and other numerical work near 20. So `alpha`
is near two thirds, not near a third, and the range asserted for it excluded the
measured value rather than bracketing it.

**Two checks agreed the new value was right and the old one wrong**, which is
the part worth keeping. The identity `a = epsilon PE / gamma` ties the
bulk-plume parameter to the entrainment rate and the precipitation efficiency,
so a value chosen on humidity grounds can be tested against entrainment rates it
was not derived from. At a precipitation efficiency of 0.33 it returns
`a = 0.204` against the 0.2 Romps quotes, and an entrainment rate of 0.198 per
kilometre against the 0.2 to 0.8 that cloud-resolving simulations diagnose. At
the assumed value it returned 0.09 per kilometre, a factor of two below the
bottom of that range — and that disagreement sat in the equation registry as a
declared known failure mode for two commits before anybody looked the efficiency
up. The over-determination was working the whole time; nobody was reading it.

**At the measured value the column is refuted.**

| alpha | 0.62 | 0.67 | 0.81 |
|---|---:|---:|---:|
| predicted relative humidity | 0.683 | 0.725 | 0.842 |
| surface temperature (K) | 301.58 | 305.18 | 313.95 |
| dT_2x (K) | 9.50 | 8.76 | 5.25 |
| feedback (W m^-2 K^-1) | +2.259 | +2.312 | +2.221 |

Thirteen to twenty-six kelvin too warm, across the whole measured range, with
sensitivities no configuration of this model should produce. And it cannot be
escaped by choosing within the range: the entrainment floor of 0.2 per kilometre
requires `alpha` of at least 0.70, which is *above* where the column has already
failed. **The two literature constraints do not overlap.**

That also disposes of the agreement reported one commit earlier. A feedback of
+1.29 W m^-2 K^-1, matching the prescribed-humidity control and published
estimates to four percent, was obtained at `alpha = 0.3` -- a value a factor of
two below the entrainment floor. It was a wrong parameter producing a
right-looking answer, which is the third time that pattern has appeared in this
pathway and the third time it was caught by checking something other than the
headline number.

#### What the refutation is actually about

The closure is not broken. At the measured parameters it predicts a
free-tropospheric relative humidity of 0.73 to 0.84, and cloud-resolving
simulations of radiative-convective equilibrium sit near 0.80 -- Singh and
O'Gorman took exactly that as a constant. It is describing a convecting column
correctly.

This framework's column is a *global mean*, and the global-mean free troposphere
is far drier than a convecting one, because it is dominated by subtropical
subsidence and by where descending air was last saturated. No bulk plume
contains that. Applying a radiative-convective closure to a global mean is a
category error, and the model has now measured its size: about seventeen kelvin
at the central value.

The caveat was on the record before the work began. When the three candidate
mechanisms were weighed, a fourth thread was noted -- Sherwood (1996),
Pierrehumbert and Roca (1998), Galewsky, Sobel and Held (2005) -- holding that
free-tropospheric humidity is set by last saturation and large-scale transport,
so that a single column cannot get it right in principle. It was recorded as
bounding what an entrainment closure could buy. It turns out to decide the
outcome rather than bound it.

Charter section 3.5 governs what happens next: `alpha` stays at the measured
value and the configuration stays refuted. Tuning it back down to 0.3 would
restore the agreement and would be calibrating against the temperature record,
which is the thing this project exists to refuse.

#### The base state looks right and that means nothing yet

The column reaches 288.1 K and a greenhouse effect of 33.1 K, against observed
values of 288 and 33. Neither is evidence, for two independent reasons.

The first is cancellation. The column carries 8.0 kg m^-2 of precipitable water
against an observed 25, a factor of 2.9 outside the conditions the vapour band
was calibrated at, *and* its lapse rate is 3.2 K/km too steep. Too little water
cools the surface and too steep a profile warms it. Two large errors of opposite
sign are agreeing to a tenth of a kelvin, which is a coincidence and not a
result.

The second is the assumption. Sweeping the evaporating condensate fraction
across the range it is asserted over moves everything:

| alpha | 0.1 | 0.3 | 0.5 |
|---|---:|---:|---:|
| predicted relative humidity | 0.250 | 0.417 | 0.583 |
| surface temperature (K) | 284.08 | 288.11 | 292.26 |
| precipitable water (kg m^-2) | 3.64 | 8.02 | 14.89 |
| feedback (W m^-2 K^-1) | +0.470 | +0.886 | +1.375 |

Eight kelvin of surface temperature and a factor of three in feedback, from one
number nobody has measured tightly. The midpoint landing on the observed
temperature is what a midpoint does. Choosing `alpha` because it produced 288 K
would have been calibrating against the temperature record, which INF-003
forbids and which this project has already caught itself doing in other forms;
0.3 was written down as the midpoint of the declared range before any column ran
and has not been touched since.

So the honest state of the vapour pathway is: the closure's two predictions both
match observations with nothing fitted, the grid dependence that disqualified
the previous mechanism is gone, no number from the coupled column is quotable
until the evaporating fraction is constrained, and the column contradicts its
own closure until the lapse rate is wired in.

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

### The cloudy column is the clear column plus the cloud

That phrasing is load-bearing, and the code did not obey it for as long as the
cloud existed. The cloudy solve was assembled by *overwriting* the base optical
depths with the cloud's rather than adding to them, which deleted the water
vapour from the cloudy fraction of the column. The cloudy fraction then radiated
to space through an atmosphere that was not there.

The symptom was a longwave cloud effect with the wrong sign: a cloud that
*increased* the outgoing longwave, which no cloud in a profile cooling with
height can do. Measured across cloud height in a moist column it ran from
-34.7 W m^-2 for a cloud at the surface to +11.4 W m^-2 for one near the
tropopause, growing more negative the more vapour lay above the cloud -- the
signature of the missing absorber rather than of any cloud physics. After the
fix the same sweep runs -0.3 to +16.5 and rises monotonically with height, which
is the expected ordering: a higher, colder cloud traps more.

Two things kept it alive. Carbon dioxide enters through
``absorption_coefficient``, which *both* solves receive, so only the vapour was
lost -- and with vapour active the band's flat background is deliberately
switched off, because the infrared background is then the vapour itself. The
cloudy fraction of a moist column was therefore very nearly bare. And no test
and no configuration had ever put a cloud and an active vapour continuum in the
same column: the cloud tests ran dry, the vapour work ran clear, and each passed
in the half of the space where the defect is invisible. No registered claim used
a cloudy configuration, so nothing published depended on it.

The comment above the code stated the correct rule the whole time. It said that
water vapour and volcanic aerosol are in the clear sky as well as under the
cloud, and the volcanic term was passed to the cloudy solve; the vapour was not.
A comment describing an invariant is not the same as a test asserting it.

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
| `convection` | — | surface sensible heat flux and dry convective adjustment; optionally with the surface as the column's bottom element, and optionally carrying the vapour it mixes |
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

## Saturation over the phase that is actually there

The first violation the domain gate found, fixed rather than argued away.

The saturation vapour pressure was integrated Clausius-Clapeyron at constant
latent heat, over liquid water, at every temperature. A column whose top is at
172 K spends most of its depth below freezing, where the stable condensed phase
is ice and the relevant latent heat is that of sublimation. Two corrections,
both first-principles:

* **Kirchhoff's relation** makes the latent heat vary, `dL/dT = c_pv - c_cond`.
  Integrating Clausius-Clapeyron with that linear form is exact rather than
  numerical and setting the heat-capacity difference to zero recovers the old
  expression, which is the limiting case the tests use.
* **Below the triple point the vapour is in equilibrium with ice**, so the
  latent heat is `L_v + L_f` and the heat capacity is the ice one. Both branches
  are anchored at the triple point, so the result is continuous with no
  blending.

The old form was **five times too high at 172 K**, 2.7 times at 200 K, and 1.45
times at 240 K. Against Murphy and Koop (2005) the new one is within 0.6 percent
from 240 to 330 K, 2 percent to 200 K, and 6.3 percent at 172 K, which is what
the declared domain of 170 to 340 K now says.

**Two calibrations had to be redone, and the framework's own rule said so.** A
calibration is only valid at the conditions it was performed at, and the vapour
those constants were fitted against had just changed. The CO2 forcing per
doubling moved from 3.70 to 3.81 W m^-2 and the outgoing longwave from 238 to
247. Refitting both against their original published targets -- 3.7 W m^-2 per
doubling and the CERES all-sky 238 W m^-2, neither of them a temperature record
-- gives a vapour wing decay of 3970 m^-1 in place of 3617, and a CO2 wing decay
of 1272.5 in place of 1273. The CO2 number returning to where it was is worth
noticing: only the vapour band was actually mis-set, and the joint fit shows it.

**What it cost the column**, none of it aimed at:

| | before | after |
|---|---|---|
| surface temperature | 305.18 K | 303.85 K |
| precipitable water | 98.31 kg m^-2 | 85.25 kg m^-2 |
| dT_2x | 8.761 K | 8.584 K |

Thirteen percent of the column's water was vapour the old formula put into cold
layers that ice would never have held. The configuration passes the validity
gate and may support a claim again.

## Ocean heat capacity does not move the equilibrium

Recorded because the opposite is easy to measure by accident, and the numbers
look convincing.

The question was whether more ocean compartments would hold the column's
temperature down -- there is a great deal of cold water not in contact with the
sun for half the year. At eleven simulated years the answer looks like yes: a
slab-only column sits at 297.3 K and one with a 1200 m deep layer at 294.5 K,
nearly three kelvin cooler.

It is not cooler. It is unfinished. The deep-ocean runs carry a
top-of-atmosphere imbalance of about 3 W m^-2 -- still absorbing -- and the
temperature climbs as the run lengthens: 294.51 K at eleven years, 294.73 at
forty-one, 295.16 at a hundred and ten, while the imbalance decays from 3.00 to
2.29. A 1200 m layer has a time constant of centuries, so a century-long run is
a third of the way there.

Shrink the deep layer until it can finish and the point is unambiguous:

| | surface (K) | TOA imbalance (W m^-2) |
|---|---|---|
| slab only | 297.335 | -0.0000 |
| 50 m deep layer, 27 years | 296.813 | 0.5323 |
| 50 m deep layer, 110 years | 297.332 | 0.0031 |

Three thousandths of a kelvin apart. Heat capacity sets *when* a column
arrives, never *where*: at equilibrium the reservoir is full and drops out of
the budget entirely.

**What is real in the idea, and where it belongs.** Ocean heat uptake genuinely
does hold transient temperature below the equilibrium value -- that is the whole
gap between transient response and sensitivity -- so compartment structure
matters a great deal for a run with a moving forcing, and would matter for a
Stage 7 hindcast. It cannot touch an equilibrium bias.

The mechanism named, cold water shielded from the sun for half the year, is
also not in this configuration at all: it runs at constant insolation with no
seasonal cycle. A seasonal cycle can rectify into the annual mean, because
outgoing radiation is convex in temperature while saturation humidity is
exponential in it, and those push opposite ways; heat capacity modulates that by
damping the cycle. It is a tenths-of-a-kelvin effect, not a fifteen-kelvin one.

And the deep ocean is cold because dense polar water sinks and fills it, which
is transport rather than capacity -- the same horizontal transport a single
column cannot have, and the same gap that stops it setting its own
free-tropospheric humidity.

## Three declared optical-depth ceilings were guesses, and a real column passed them

Recorded because the measurement reverses the intuition the records were written
from, and because the run that exposed it is the one the carbon dioxide forcing
claims rest on.

Publishing the optical depth made four checks possible at once and immediately
found the greenhouse configuration running at **tau = 1110 in a layer**, against
declared ceilings of 10 (diffusivity), 100 (Beer-Lambert) and 1000 (two-stream).
Every one was exceeded, so the historical replay was reporting itself unable to
support a claim.

The measurement says the ceilings were wrong, not the run.

**The diffusivity approximation.** Its record said accuracy "degrades toward
both the optically thin and optically thick limits", and its domain was floored
at 0.01 and capped at 10 to match. That reads the *relative* error as though it
were the absolute one. Against the exact hemispheric integral, computed by
quadrature rather than from a library:

| | max absolute transmittance error |
|---|---|
| tau in [0, 0.01], below the old floor | 0.0030 |
| tau in [0.01, 10], the old declared range | 0.0316 |
| tau = 10 | 7e-6 |
| tau = 30 and beyond | below double precision |

An order of magnitude *more* accurate below the floor than inside the range the
floor was protecting, with the worst error at tau = 1.28, in the middle. The
shape is two humps rather than one, because the signed error changes sign near
tau = 0.4: the approximation over-transmits a thin layer and under-transmits a
thick one. A first attempt at the test asserted a single monotonic rise to the
peak, which looked right and is false; the test now pins the structure.

What survives of the thin-limit concern is real but is not a bound: the
absorptance of a thin layer is short by the fixed ratio D/2 = 0.83, a 17 percent
underestimate however thin the layer is. It is 17 percent of a vanishing
quantity, the flux error goes to zero with tau, and the radiation module already
declares `diffusivity_approximation` as an uncertainty output. That is where a
known bias belongs. A validity domain answers whether a relation may be
evaluated; an uncertainty answers how wrong it is there. Conflating them
declared an ordinary column invalid.

**The two-stream coefficients and Beer-Lambert.** Measured rather than argued:
reflectance, transmittance and emissivity sum to one to the last bit at every
optical depth from zero to 1e8, at every single-scattering albedo up to and
including conservative scattering, with and without forward asymmetry. The
transmittance underflows smoothly to zero -- the correct value, not a failure --
and the layer saturates at its right limits, emissivity 1 for an opaque absorber
and 0 for a conservative scatterer. The solver carries an explicit large-tau
branch; it was written for this and the ceiling denied it. Beer-Lambert was
contradicting itself outright: the record declares the limiting case "tau ->
infinity gives zero transmittance" and then refused to be evaluated past 100.

All three now carry no upper bound on tau, with the measured error bound stated
as a condition. Negative optical depth is still excluded, because that is a
defect rather than a thick layer.

**And a caution that came with it.** The approximate layer transmittance is
exactly multiplicative, so splitting a path into layers leaves the approximate
answer unchanged while moving every per-layer optical depth. Sixteen layers of
tau = 0.08 each -- comfortably inside any range -- carry the path optical depth
1.28 at which the approximation is at its worst. The accuracy is a property of
the path, not of a layer, so a per-layer range is necessary and not sufficient.

## The cost of a diagnostic is part of the diagnostic

The optical depth was first published as the full `(sky, layer, band)` array
every step. That is 80 MB across the 25,916 steps of the historical replay,
which was enough to exhaust memory on a run already near the limit -- in order
to carry information that is consumed as a minimum and a maximum.

A declared range is a statement about the *extent* a relation was evaluated
over, which is two numbers whatever the resolution behind them. The modules
publish the extent; the full arrays stay on each `ColumnRadiationSolution` for
the current step, where anyone who wants to know which layer is thick can still
ask. Coverage is identical either way, because the check only ever took the
minimum and the maximum.

The general point: a diagnostic added to make something checkable is still a
diagnostic that runs every step of every configuration, and publishing at a
resolution nothing consumes is a cost with no claim attached to it.

## What the column stores is its metadata, not its physics

The preflight memory estimate was wrong by a factor of twenty, and fixing it
turned up a better fact than the fix.

The estimate summed the modules' declared working sets, which is a real quantity
and not the one that matters. What grows is the **retained trajectory**: the
driver keeps every step's published fields so the manifest can be audited, and
that is linear in the step count while the working set is not. A configuration
that fits at 2000 steps and not at 200000 differs only there, and the estimate
said nothing about it -- which is how a 120000-step sweep reached 4.3 GB and had
to be killed.

Measured across six configurations from 4 to 24 layers and 8 to 32 bands, the
retained cost fits

    bytes per step = 714 * n_fields + 17.6 * n_values

to within 8.8 percent. **The per-field term is forty times the per-value term.**
The moist column's 49 fields carry 1.7 kB of floating-point numbers per step out
of 37.7 kB retained; the other 36 kB is each field's units, quantity kind, time
semantics, spatial support, provenance, uncertainty, sign convention and
missing-value policy.

That is worth stating plainly, because it is a fact about what this framework
*is*. The metadata is not overhead around the model -- it is most of what the
model stores, by a factor of twenty, and it is the part that makes a run
auditable. The cost of being able to say what a number means is about 700 bytes
per field per step.

**The field count needs no declaration.** It is ``len(produced_outputs)`` summed
over the participating modules, and it matches the retained count exactly on
every configuration tested -- 35 for the dry column, 49 for the moist one. The
value count does not follow from the contracts, so the preflight omits that term
and says so in a warning; the estimate then comes out 9 to 19 percent low instead
of a factor of twenty low, and the driver records the exact figure once the
fields exist. Measured after the fix, estimate against actual: 1.08 at 8 layers
and 1.003 at 24.

The 120000-step sweep now estimates 4.2 GB before it starts, against the 4.3 GB
it reached.

## The error in kelvin, reported by the suite instead of by hand

The number the project is judged on had no test. It was recomputed in a scratch
script whenever anyone wanted it, so it moved between commits without anyone
noticing and whatever figure appeared in this file was whatever the last person
to run a script had seen. It is now measured by
`tests/acceptance/test_surface_temperature_error.py`, and every test that
equilibrates a column reports its own through the `report_bias` fixture, with
the whole set printed together at the end of the run.

**Where the column stands**, against an observed global-mean surface temperature
of 288.0 K, each as a time-mean over the last fifth of a 16000-step run:

| configuration | error | evaporation | PW | Bowen |
|---|---|---|---|---|
| airless body, no atmosphere | **-33.42 K** | — | — | — |
| dry radiative-convective | **-4.49 K** | — | — | — |
| vapour present, radiatively inert | **-5.22 K** | 0.99 | 8.6 | 0.810 |
| vapour radiatively active | **-2.55 K** | 2.46 | 8.4 | 0.754 |
| vapour active and transported | **+6.73 K** | 4.05 | 28.3 | 0.339 |
| observed | 0 | 2.7 | 25 | 0.26 |

The last row of the model is the answer: **+6.7 K**. The sequence above it is
what the answer is made of, and the shape matters more than any single entry.
The dry column is *cold*; every step toward a complete water cycle warms it,
because the vapour is a greenhouse gas and its absorption beats its own
evaporative cooling. So the bias changes sign along the sequence, and a
configuration that happens to land near zero in the middle of that is passing
through rather than right -- which is exactly what "vapour radiatively active"
at -2.55 K is doing while carrying a third of the observed water.

The -33.42 K of the airless body is the useful end of the scale: it is the
greenhouse effect, measured, and it is the quantity the other rows are eating
into.

### Three ways this number has been quoted wrongly, and what stops each

**From a run that had not equilibrated.** The dry radiative-convective column at
a large heat capacity reads -16.3 K at 20000 steps. It is not 16 K too cold; it
is unfinished, carrying +43.8 W m^-2 and still warming. A deep ocean under the
same column looked 3 K cooler at eleven years for the same reason. So
`BiasReport.error_k` is **`None`** unless the mean top-of-atmosphere imbalance is
below 0.1 W m^-2, and the report says what the imbalance was instead. `None`
rather than a wide error bar, because a bias from an unfinished run is not an
uncertain result -- it is a different quantity, and returning a number would let
it be compared against a real one.

The threshold was set against what the imbalance *is* rather than what would be
convenient: an equilibrated column here sits at 0.001 W m^-2 and an unfinished
one at 43.8, so anything between the two separates them and 0.1 leaves three
orders of magnitude of margin on the side that matters.

**From an endpoint rather than a mean.** The column limit-cycles with a surface
span of 2.08 K that does not shrink with the timestep, so a single step is a
point on a cycle. Every figure here is a tail mean.

**Without its process quantities.** This column has sat 3 K from the observed
temperature while its evaporation was twelve times too small and its Bowen ratio
thirty-seven times too large -- two large errors of opposite sign cancelling at
the one number anyone looks at. So the report carries the water cycle beside the
temperature and the table above has four columns rather than one.

One more guard, in the other direction: the observational records disagree among
themselves by about 0.5 K, so a modelled bias smaller than that is not a
distinguishable result and `BiasReport` marks it as inside the spread rather
than quoting three decimals of it.

## Two requirements that were about declaring rather than computing

**CMP-007**: higher fidelity shall record its scientific justification in the
experiment manifest. Nine of the eleven module implementations name a
lower-fidelity alternative -- a single-layer column, a grey atmosphere, a dry
column, an airless slab -- and then do not use it, which is a scientific choice
that was living in whoever configured the run. Each now states what the cheaper
form would *fail to produce*, because that is the falsifiable half: "a single
layer has no lapse rate, and the lapse rate is what convection acts on" can be
checked by running the cheaper form and looking, where "more accurate" cannot.
`ModuleContract.validate` refuses a contract that names an alternative and gives
no reason, so the requirement is structural rather than a habit.

**ACC-010**: the reference kernel shall complete its baseline suite within the
ENGINEERING_WORKSTATION budget. Bound to *measured* resources rather than to the
estimate, since the estimator already gates every run and binding the requirement
to it would bind it to the thing already enforced. Measured: the suite uses 0.42
percent of the work-unit budget and its worst member holds 80 MiB against 16 GiB.

Writing it found something. **The preflight memory estimate is wrong by two
orders of magnitude**, in the direction that admits a run which then exhausts
memory -- which has already happened here once, at 4.3 GB on a 120000-step sweep
that had to be killed. Every module declares
`persistent_output_bytes_per_step = 0` while the driver retains every step's
outputs for the manifest, so the term that actually grows is invisible to the
estimate. Measured on the moist column the retained trajectory costs 37.7 kB per
step and scales linearly to four decimal places; 120000 steps is 4.5 GB, which is
the figure that was hit. Recorded as open work, with a test that asserts the
discrepancy so it fails when the estimator is fixed.

## The six uncertainty methods are one mechanism, and the refusal is the point

UNC-005 asks for bounds, sensitivity derivatives, parameter sweeps, ensembles,
Monte Carlo propagation and model-form comparison. There were interfaces and no
implementations, and the evidence that it mattered was accidental: the bulk-plume
parameter sweep that established a factor of 2.7 in the feedback and eleven
kelvin of surface temperature was run by hand in scratch scripts. Four more have
been run by hand since -- the evaporating fraction, the timestep, the exchange
coefficient, and gravity in the section above.

**All six are the same mechanism.** Each runs a configuration several times with
something changed and collects an observable; they differ only in *how the
members are chosen*. So there is one `Ensemble` and six constructors:

| | members are |
|---|---|
| `bounds` | the two ends of a declared range |
| `sweep` | points along one parameter |
| `ensemble` | joint assignments, chosen |
| `monte_carlo` | joint assignments, sampled from declared distributions |
| `model_forms` | different configurations |
| sensitivity derivative | a secant over a numeric `sweep` |

Writing six mechanisms would have made the difference between an ensemble and a
Monte Carlo sample a matter of which function was called. It is not: it is a
matter of how the members were chosen, and that is the only thing carrying any
information about what the spread means.

**Which is UNC-007, and is why this module refuses more than it computes.** A
three-point sweep across an asserted range tells you what the model does across
that range. It says nothing about how likely any part of it is, so a number
between 0 and 1 printed beside it will be read as a probability it has not
earned -- the same false precision as quoting an assumed parameter's midpoint,
which this project has warned about since the evaporating condensate fraction.
So a `PredictiveInterval` carries its meaning and its nominal coverage, and a
method that cannot support a coverage must state `None` **and say why**. An
absence with no reason reads as an oversight.

The refusal runs both ways, which is what makes it a distinction rather than a
prohibition. A Monte Carlo sample from declared distributions *can* state a
coverage, and refusing it one would be as wrong as granting one to a sweep.

`Ensemble` also requires a `basis`: where the values came from. A range asserted
from taste and a range read out of the literature produce identical numbers and
mean different things, and the evaporating condensate fraction is the case --
swept over 0.1 to 0.5 when published precipitation efficiencies put it near two
thirds, so the sweep excluded the answer and nothing in its output said where the
bounds had come from.

**UNC-006 is the factor of two that leaves no trace.** Root-sum-square is valid
only for independent standard uncertainties. A 95 percent half-width is about two
sigma, so passing one in unconverted halves the combined figure and nothing about
the result looks wrong. `UncertaintyKind` already existed because real records
disagree about which they publish -- HadCRUT5 a 95 percent interval, NOAA's heat
content a standard error -- and this is the place that disagreement would do its
damage. Both conditions are checked, and independence has to be *justified*
rather than asserted.

### The test corrected me, which is the reason to write it that way

The sensitivity check runs against a closed form: an airless body sits at
`T = (S(1-a)/4 e sigma)^(1/4)`, so `dT/dS = T/4S` exactly. Landing on that would
be worth something and could be a coincidence. Missing it by a *predictable*
amount cannot be, so the test asserts the finite-difference error too -- and I
wrote it with the sign backwards.

The reasoning that fails: `T` is concave in `S`, so the secant should fall below
the tangent. The measurement says otherwise, and it is right. A secant measures
the *mean of the derivative* over the interval, `f'` is convex here (`f'' < 0`
but `f''' > 0` for a fourth root), and by Jensen the mean of a convex function
exceeds its midpoint value. The excess is `h^2 f'''(S_0) / 6`, which is
8.837e-6 predicted against 8.839e-6 measured -- four significant figures, and a
far stronger statement than the inequality I set out to write.

Traceability 173 to 176 of 185.

## The gravity constant was an approximation nobody had to declare

GRV-001 asks that gravity be provided through a replaceable interface beginning
with a declared approximation. It was provided by
`value_of("earth_standard_gravity")`, read straight from the constants table at
four call sites, and the constant's own note said that the gravity module
supersedes it.

The requirement turns out not to be about gravity. **An approximation nobody has
to declare is one nobody has to measure**, and the hydrostatic equation's record
proved the point by asserting a bound that was wrong for the column it was
attached to. It said constant gravity is "about 1.5 percent at 100 km and
negligible through the troposphere". Both halves are true and neither is about
this column, whose top at 1000 Pa is **19 to 27 km up** depending on the
resolution.

Held against the inverse square law over the column's own heights:

| | constant g | Newtonian g(z) | |
|---|---|---|---|
| column mass | 10230.30 | 10252.15 kg m^-2 | **+0.21%** |
| topmost layer | — | — | **+0.61%** at 8 layers, +0.85% at 32 |
| surface temperature | 272.444 | 272.481 K | +0.037 K |
| outgoing longwave | 198.449 | 198.365 W m^-2 | -0.085 |

The surface numbers are the time-mean over the last tenth of a 20000-step run and
repeat to 0.001 K at 8 and 16 layers, so they are the approximation's cost and
not the limit cycle's noise. The error sits in the worst place available: the
topmost layer is the one that sets the outgoing longwave, and it is the layer the
constant misses by the most.

**Two of charter 9.15's four rungs are implemented.** `ConstantGravity` is what
the framework has always done, now saying so. `NewtonianAltitudeGravity` is
`g(z) = g_0 (R/(R+z))^2`, checked against its own series expansion term by term
rather than against the closed form -- agreeing with the closed form would
establish only that the code computes the formula it was handed. The spatial
field and the relativistic correction are carve-outs in the charter's sense.

**Latitude is deliberately not implemented**, although charter 9.15 pairs it with
altitude on the second rung, and not because it is small: equator to pole is
about 0.5 percent, the same size as the altitude effect measured above. It is
omitted because every configuration here is a global mean or a single column,
neither of which has a latitude to evaluate it at, and because the conventional
value is already a latitude-weighted compromise. Supplying one would be inventing
a coordinate the physics does not have.

**The default stays constant, and that is a decision rather than an omission.**
GRV-001 says "beginning with a declared approximation", and charter 9.15's first
rung is constant g. The switch is a keyword away and the downstream calibration
was checked before recommending against it: the carbon dioxide forcing per
doubling moves from 3.6977 to 3.6943 W m^-2, which is 0.0034 of a +/-0.05 window,
so a physics change here would *not* invalidate the fit the way the saturation
correction did. What argues against promoting it is the trade, not the risk. The
column carries a +6.6 K warm bias; spending a renumber of every recorded result
in this file on 0.037 K is not where the accuracy is.

The interface refuses rather than falls back. A height-dependent field asked for
layer masses without layer heights raises, instead of quietly using the surface
value -- which would put the approximation back to being applied without anyone
declaring it, which is the whole thing this was written to stop.

## A carve-out is a decision, and counting it as a backlog item hid four of them

`docs/todo.md` complained that the carve-out requirements were recorded as
`not_started`, so the headline figure read as twenty-nine things left undone when
"about ten of them are architectural decisions". The complaint was right and the
estimate was not. Looking made it **four**.

The classification is the specification's rather than a judgement made here.
Charter 4.3 defines an architectural carve-out as a reserved and documented place
for a domain, guaranteeing an implementation may be added and asserting nothing
about whether the domain is important or active. Charter 18 lists the fifteen
candidate active modules of the first executable kernel and says that exclusion
from it "does not remove the corresponding architectural interface". IRS 23
repeats it for nine named extensions.

**The rule that decides the hard cases: a carve-out is a domain the kernel does
not use.** A domain it uses without the declared interface is a gap.

Gravity is what forces the distinction and is the reason the number is four
rather than five. The gravity module is as absent from charter 18's candidate
list as the cryosphere is, so the mechanical reading would excuse GRV-001 with
the rest. But the column computes hydrostatic layer mass with gravity on every
step, out of a bare constant whose own note says the module supersedes it. That
is not a reserved place for an unused domain; it is a used domain reaching the
physics without the replaceable interface the requirement asks for. GRV-001 stays
counted.

So: BIO-001, ICE-001, LND-001 and PAR-001 are carved out, each naming the clause
that defers it. The nine spatial and uncertainty requirements that were loosely
lumped in with them are **not** — the kernel has a vertical column whose
resolution changes in every convergence study, and it reports uncertainty from
every module, so both are domains it uses.

**The status has to be load-bearing or it is a way of making the number look
better.** `tools/traceability.py` refuses a carve-out with no waiver, and refuses
one that also has a bound test, because a requirement cannot both be deferred and
be verified and whichever of the two is wrong should be said out loud. The matrix
reports carve-outs against their own line rather than folding them into the
denominator: **172 of 185**, with 4 carved out, where the same tree read 160 of
189 before any of this week's work.

## Carbon dioxide is the only human influence here, and it enters as an opacity

HUM-001 asks that human activity alter model state only through explicit physical
inputs, sources, sinks, transfers or boundary changes. The failure it names is a
fitted forcing added to a budget: a term in watts per square metre, tuned to
reproduce a record, that no physical process put there.

The carbon dioxide pathway is the opposite, and it is now demonstrated rather
than described. `PrescribedGreenhouseModule` publishes a mole fraction and
nothing else: it owns no state, computes no tendency, declares no transfer and
conserves nothing, and its single graph edge is a `parameter_dependency` to the
radiation module. No state budget carries a label naming it.

That much would also be true of a module that does nothing, so the check has two
more parts. It *acts*: doubling the concentration lowers the outgoing longwave,
through the absorption coefficient alone, with the transfer solver unchanged and
unaware. And it adds no energy: both runs close their budget to one part in
10^16, which is what a mechanism entering through a constitutive relation must
do and what a forcing added to a budget could not.

## Nine requirements were satisfied by machinery nobody had run

`docs/todo.md` carried an item saying that about ten unbound mandatory
requirements needed a test rather than an implementation, because the machinery
they describe is enforced elsewhere. It listed them in a table. That table was
written from reading the code, and reading is not measuring.

Nine of them were true. Every module implements configuration validation, every
contract carries a cost declaration, every dependency edge is typed, the state
registry enforces single ownership, the registries supersede rather than
overwrite. All of it held on every reference configuration in the tree, which is
what `tests/acceptance/test_framework_self_compliance.py` now establishes.

**One was not.** ACC-005 requires every mandatory run to produce manifests,
dependency graphs, equation bindings, provenance records and *output checksums*.
The first four were there. `output_checksums` was **empty on every run in the
tree**: the only path that fills it, `RunManifest.record_output`, hashes a file,
and no reference configuration writes one. The requirement had been recorded as
already satisfied on the strength of a method nothing called.

So the run now fingerprints itself in memory, with three digests:

| digest | what it identifies |
|---|---|
| `final_outputs` | every field the last step published |
| `final_state` | every module's committed state at the end |
| `conservation_trajectory` | the worst relative residual, per step |

The third is the one worth arguing for. An endpoint digest cannot distinguish
two trajectories that met at the end, and KER-014 asks for a reference run to be
*reproducible* from its manifest, which is a claim about the path. A scalar per
step costs nothing -- 2000 steps of the dry column fingerprint in the noise of
the run -- and it catches a divergence that recovered.

**On writing sweeps rather than examples.** Each of these checks runs over all
seventeen reference configurations, and
`test_the_sweep_reaches_every_module_implementation` fails when a module exists
in `mfespf.physics` that none of them exercises, naming it. A sweep with a hole
in it makes a claim it has not earned, and the hole is exactly what nobody
notices.

Two of the checks were written over-strict and had to be corrected by
measurement rather than by taste, which is worth recording because both
corrections went the same way -- the framework was right and the test was wrong:

* **The lifecycle is not the charter's numbering, deliberately.** Charter 6.3
  numbers `advance_state` before `emit_outputs`, and the coupler calls them the
  other way round, publishing every module's diagnostics in the tendency pass and
  committing state afterwards. That is *stronger* than the per-module sequence:
  it is the operator-splitting invariant, and the convection/water over-draw
  earlier in this project came from violating it. The test asserts the semantics
  6.3 says to preserve -- inputs before computation, every tendency before any
  commit -- rather than the order it lists them in.
* **Six edges cite no equation, correctly.** Handing a module a heat capacity, an
  albedo or an irradiance record is a dependency, not a relation. Only a
  `physical_transfer` is computed by an equation and has to name it.

Traceability moves from 160 to 170 of 189.

## The domain check was reading the last step, and a run is not its last step

The validity check had a limitation that had been written down and left open:
it judged a run by its **final state**. That is a different and much weaker
question than the one a declared domain asks, and it gets the answer wrong in
both directions. A closure that ran far outside its conditions while the column
span up and settled back inside them is reported clean. A closure that is merely
*idle* on the last step leaves no footprint to be judged on at all.

The convective adjustment is the case that forced it. It mixes on 79 percent of
steps on the dry column and on none of the first four hundred, the column being
dry-stable while it spins up, so the identical configuration came out resolvable
at 2000 steps and unresolvable at 400 -- a fact about where the run stopped
rather than about the physics. The repair at the time was to widen the module's
NaN mask to the *union over the run*, which left something behind at the end for
the check to find. That worked and it was the wrong shape: the union answers
"which layers does this relation govern" while the question is "was it in range
when it ran", so a layer mixed at 250 K during spin-up was judged at whatever it
held when the run stopped.

**The fix is to accumulate rather than reconstruct.** `EvaluatedExtents` folds
each step's published fields into a running minimum and maximum, and the ranges
are judged against that. The per-step rules are unchanged and now apply per
step: a NaN says the relation was not evaluated there and is skipped, an
infinity is a defect that disqualifies the field for the whole run rather than
being masked away, and a field that is NaN at every step never enters -- which
is the correct report that its relation never ran, and is what the dry column
gives at 200 and 400 steps.

Resolution still happens once, against the final state. Which field a relation's
`T` refers to is a property of the declared contracts, not of any one step; only
the values vary. Separating the two is what makes the change small.

**What it was hiding.** No reference configuration acquires a violation, which
is worth saying since accumulating is strictly wider and a new violation would
have meant the run had been failing all along unseen. But the extent being
judged changes substantially:

| | at the final step | over the run |
|---|---|---|
| moist column, `EQ-DRY-CONVECTIVE-ADJUSTMENT` `T` | [276.1, 287.5] K | [235.3, 287.5] K |
| dry column at 6000 steps, same | [234.7, 256.1] K | [228.9, 256.1] K |

The moist column's range is **4.6 times wider** over the run. It stays inside
the declared [100, 400], so nothing is refused -- but the check had been
examining a fifth of the interval the relation was actually evaluated over, and
that it passed was not something the old check had established.

**And the mask goes back to meaning what it says.** Convection publishes this
step's mixed groups again; `_ever_mixed` is gone. A repair that removes the
machinery an earlier repair needed is the sign the second one was at the right
level.

The cost is 7.3 percent of a 20000-step run, after taking the two-pass route: a
plain minimum and maximum propagate a NaN, and that propagation is the signal to
pay for building a mask rather than a reason to build one on every field on
every step. The first version did it unconditionally and cost 17.5 percent.

## The evaporation coefficient was never free, and that is the Bowen ratio

The plan here was to calibrate the evaporation coefficient against the observed
evaporation rate, which would at least have been INF-003-compliant: set it from
something other than the temperature and report the temperature as a prediction.
Checking the number first made the calibration unnecessary.

The bulk formulae are

    H = rho c_p C |U| (T_s - T_1)   =  k_H (T_s - T_1)
    E = rho     C |U| (q_s - q_1)   =  k_E (q_s - q_1)

and the ``C`` in them is the *same* ``C``: the eddies that carry heat across the
surface layer are the eddies that carry moisture, which is the Reynolds analogy
and why most texts write one coefficient. So

    k_E = k_H / c_p

exactly, with no new number. The configured pair was ``k_H = 10 W m^-2 K^-1``
against ``k_E = 2.0e-4 kg m^-2 s^-1``, and ``10 / 1004.6 = 9.95e-3``. **Fifty
times apart.** The column was moving heat fifty times more efficiently than
moisture using eddies it had declared to be the same ones.

That is the whole of the Bowen-ratio error, which had been recorded as a
separate symptom -- 9.6 against an observed 0.26, "37x too large" -- for as long
as the hydrological cycle had been measured. It was one inconsistency, not two
errors.

**What it predicts.** Tying them makes the partition a consequence rather than a
setting, so it can be wrong. Measured, timestep-converged:

| | configured, 2.0e-4 | derived, 9.95e-3 | observed |
|---|---|---|---|
| Bowen ratio | 9.6 | **0.349** | 0.26 |
| evaporation | 0.22 mm/day | 4.00 mm/day | 2.7 |
| latent flux | 6.5 W m^-2 | 118 W m^-2 | ~82 |
| sensible flux | 62.0 W m^-2 | 37.5 W m^-2 | ~21 |
| precipitable water | 3.1 kg m^-2 | 29 kg m^-2 | 25 |
| surface minus air | 6.2 K | 3.0 K | ~1.5 |
| surface temperature | 291.1 K | 294.6 K | 288 |

The partition error collapses from a factor of 37 to a factor of 1.34, with
nothing fitted to it. The precipitable water lands at 29 against 25, having been
3.1. And the surface bias goes from +3.1 K to +6.6 K, which is the compensating
pair separating again and is the right direction for the same reason as before:
the +3.1 was two large errors nearly cancelling.

**The remaining evaporation excess is not an independent error.** 4.0 mm/day
against 2.7 is a factor of 1.48, and the surface is 6.6 K too warm. Clausius-
Clapeyron over 6.6 K is a factor of about 1.6 in the saturation humidity that
drives the bulk formula, so the excess is what a column this warm must evaporate
at this coefficient. Fix the warm bias and the evaporation follows; there is
nothing separate left to explain here.

## The column stops settling and starts cycling, and the timestep says it is real

At the derived coefficient the column no longer converges to a steady state. Over
the last fifth of a 60000-step run the surface spans 3.4 K, the top-of-atmosphere
imbalance spans -4.8 to +11.1 W m^-2 and the latent flux spans 100 to 237, while
the *mean* imbalance is 0.3 W m^-2. An endpoint read off that is not an
equilibrium whatever its imbalance happened to be on the last step.

Convection names its own suspect -- "instantaneous mixing means the result is
independent of the timestep, which is convenient and physically wrong at short
timesteps" -- so the test is the timestep. The same physical span at three:

| dt (s) | steps | T_s mean | T_s span | imbalance mean | imbalance span | mm/day | Bowen |
|---|---|---|---|---|---|---|---|
| 3600 | 20000 | 294.735 | 2.053 | -0.001 | 6.854 | 4.051 | 0.339 |
| 1800 | 40000 | 294.628 | 1.842 | 0.004 | 7.228 | 4.012 | 0.346 |
| 900 | 80000 | 294.581 | 1.756 | -0.000 | 7.443 | 3.999 | 0.349 |

Two findings, and they point opposite ways.

**The mean state is converged.** 294.735, 294.628, 294.581 K; evaporation 4.051,
4.012, 3.999 mm/day; Bowen 0.339, 0.346, 0.349. Halving the step twice moves the
surface by 0.15 K in total and the mean imbalance sits at zero to three decimals.
So the numbers in the table above are a resolved time-mean equilibrium and may be
quoted as one.

**The oscillation is not a timestep artefact.** The surface span converges to
about 1.75 K rather than to zero, and the imbalance span *grows* slightly as the
step shrinks. This is a limit cycle of the coupled closures, not numerical noise:
moisture accumulates in the boundary layer, destabilises it, the adjustment fires
and mixes the whole group at once, the lifted vapour condenses in a burst that
restabilises the column, and it rebuilds. All-or-nothing instantaneous mixing
against a continuously recharging moisture supply is a relaxation oscillator.

It is a closure artefact rather than a numerical one, which is a distinction the
timestep study is what establishes. Real convection is continuous and spread
across a domain; a single column with an adjustment that homogenises completely
or not at all has nowhere to put the intermediate state. Recorded as open work
rather than repaired, because the cure is a convective closure with a timescale
-- a different scheme, not a tuning of this one.

## The saturation was over the stable phase and the energy was not

Recorded twice before as a percentage and left open both times because "fixing
it moves the climate and needs a recalibration". What had never been identified
was the obstacle, and once identified it turned out to have a clean answer.

`saturation_vapour_pressure` has saturated over the *stable* phase since the ice
branch was added, so the model already decides that a layer below the triple
point makes ice. The energy bookkeeping then released the latent heat of
vaporisation for it, understating the release by the whole heat of fusion. In the
configured moist column **half the water sits below the triple point** -- the
mixed group spans 284.9, 273.6, 261.0 and 246.7 K -- so this was not a fringe
case.

**The declared domain hid it instead of catching it.** `T` over [230, 340] with
no stated condition. Measured, the constant-L_v error is not a gradual
degradation with a floor at 230; it is a *step* at the triple point and then a
plateau:

| T (K) | L(T) | constant L_v | error |
|---|---|---|---|
| 288.0 | 2466378 | 2501000 | +1.40% |
| 273.2 (liquid) | 2500907 | 2501000 | +0.00% |
| 273.1 (ice) | 2834713 | 2501000 | **-11.77%** |
| 240.0 | 2842028 | 2501000 | -12.00% |
| 230.0 (the floor) | 2844238 | 2501000 | -12.07% |
| 170.0 | 2857498 | 2501000 | -12.48% |

The floor sat in the middle of a region where the error is flat at about 12
percent, so crossing it meant nothing: the check fired at 230 K and passed at
260 K on an identical error. The constant-L_v form was within 1 percent only
over [273.2, 283.8] K -- a ten-kelvin window above the triple point.

**The obstacle, which is why this was not a swap.** Making the latent heat vary
with temperature makes the module's energy *store* vary with temperature, and
then a layer that simply warms changes the store with no water going anywhere
and no transfer to balance it. That is why the constant-L_v approximation is
self-consistent rather than merely convenient, and it is the reason two previous
attempts stopped here.

**The answer is to follow the water.** The store keeps a constant reference
latent heat, so it never moves with temperature. The phase enters as the *excess*
released at deposition, and that excess is traced to where it comes from:

    liquid -> vapour at the surface      absorbs L_v
    vapour -> ice aloft                  releases L_v + L_f
    ice -> liquid at the surface         absorbs L_f

Over the cycle the surface pays the fusion energy and the atmosphere receives
it, which is exactly what freezing aloft and melting at the ground do. So it
enters as a declared transfer, `deposition_fusion_energy`, surface to water,
with both sides declaring it, and the surface carries a labeled
`melting_of_frozen_precipitation` cooling. Nothing is created and the
reconciliation closes to the same tolerance as before.

**What it moved.** At the configured evaporation coefficient, -0.02 K: the phase
error multiplies a mass flux twelve times too small, which is what was predicted
of it and is now measured rather than estimated. At a coefficient calibrated to
the observed evaporation it matters -- the column ends up 3 kg m^-2 wetter and
the run needs re-equilibrating -- which is why this had to be fixed *before* the
calibration rather than after.

**And the declared domain was inverted.** With the phase term the error against
the Kirchhoff-varying latent heat is below 1 percent from 170 to 284 K, below 2
percent to 294 K, and 6.7 percent at 340 K. The old [230, 340] excluded the
region where the relation is now most accurate and admitted the region where it
is worst. The new range is [170, 340] with the accuracy stated as a condition,
and the floor is the saturation relation's own, because below it no vapour is
placed and there is nothing for this relation to act on.

Two existing integration tests failed on this change, correctly: both asserted
`L_v` per kilogram for every layer, which was the bug rather than the
specification. They now assert the phase-weighted sum, layer by layer, and one
of them additionally requires that some layer condensed below the triple point,
so the phase term cannot silently stop being exercised. A third had to stop
taking the *first* surface-to-water energy residual, because there are two
energy channels now and picking one was unambiguous only while there was one.

## What the warm bias actually is, assembled in one place

The figure this project has carried -- seventeen kelvin -- had never been set
beside the other columns' biases, and it could not be, because they are
different radiative configurations. Assembled, with the top-of-atmosphere
imbalance beside every temperature:

**Window-absorption family**, 20000 hourly steps, equilibrated:

| configuration | T_s | bias | PW (25) | mm/day (2.7) |
|---|---|---|---|---|
| dry radiative-convective | 283.50 | -4.50 | 0 | 0 |
| moist, vapour radiatively inert | 283.47 | -4.53 | 7.8 | 0.04 |
| moist, vapour radiative | 288.04 | **+0.04** | 7.4 | 0.09 |
| moist + transport (configured) | 291.12 | +3.12 | 3.1 | 0.22 |
| moist + transport, C_E calibrated | 297.85 | **+9.85** | 23.0 | 2.53 |

**Greenhouse-grid family**, 20 layers, 3000 daily steps:

| configuration | T_s | bias | PW | imbalance |
|---|---|---|---|---|
| prescribed RH 0.40 | 288.31 | +0.31 | 7.5 | -0.38 |
| prescribed RH 0.45 (near observed) | 289.51 | +1.51 | 9.2 | -0.30 |
| prescribed RH 0.50 | 290.70 | +2.70 | 11.1 | -0.20 |
| prescribed RH 0.725 | 296.15 | +8.15 | 23.3 | +0.65 |
| bulk plume, alpha 0.3 (assumed) | 289.79 | +1.79 | 14.3 | -0.28 |
| bulk plume, alpha 0.67 (measured) | 299.84 | **+11.84** | 60.4 | +2.55 |

**The answer is about ten kelvin, and the two routes to it agree.** The only
window-family column with Earth-like process quantities is +9.85 K. The plume
column at its measured parameter is +11.84 K and still carrying 2.55 W m^-2 of
imbalance, so it is not converged -- it is on its way to the 305 K recorded
elsewhere in this file, which is the seventeen.

Three things this makes visible that no single number did.

**The best-looking column is the least trustworthy.** "Moist, vapour radiative"
is +0.04 K -- dead on observation -- while holding 7.4 kg m^-2 of water against
25 and evaporating 0.09 mm/day against 2.7. It is the compensating-error pattern
in its purest form and it is the configuration a reader would pick out of this
table if the process columns were not next to the temperature.

**The bias is not a sum of small terms.** Across the prescribed-humidity rows one
lever moves the answer from +0.31 to +8.15 K, and the plume rows move from +1.79
to +11.84 K on the value of a single parameter. A factor of 1.6 in relative
humidity produces a seven-fold error in precipitable water and of order ten
kelvin, because the vapour feedback compounds it. A small cause with a large
effect is the signature of one mechanism under strong feedback.

**The dry column is cold, not warm.** -4.50 K with no water at all. So the bias
is not a background offset that water sits on top of; the entire warm bias and
more is the water, which is why every route to it goes through the humidity
closure.

## Refining the grid walked the column out of its own equation's domain

Found while assembling the table above, and it invalidated every row of the
greenhouse-grid family: `EQ-CLAUSIUS-CLAPEYRON-SATURATION: T spans [159.1,
294.0] K, outside the declared [170, 340]`, so `may_support_a_claim` was False
for all of them.

Unlike the three optical-depth ceilings, this floor is not a guess. The record
measures it: against Murphy and Koop (2005) the constant-heat-capacity form is
within 2 percent to 200 K, 6.3 percent at 172 K and 13.6 percent by 150 K. The
floor is right and the column really was evaluating the relation where it is
wrong.

What puts it there is the grid. A radiative-convective top layer gets colder as
it gets thinner:

| layers | T_s | coldest layer | water below 170 K | as a fraction of the column | refused |
|---|---|---|---|---|---|
| 8 | 293.15 | 174.40 | 0 | 0 | no |
| 12 | 294.63 | 166.99 | 2.6e-5 | 1.3e-6 | yes |
| 20 | 295.81 | 159.07 | 3.7e-6 | 1.6e-7 | yes |
| 30 | 296.39 | 153.70 | 1.4e-5 | 5.7e-7 | yes |
| 40 | 296.69 | 150.30 | 5.0e-6 | 2.1e-7 | yes |

**VAL-002 asks for the refinement and the domain check refuses the result, and
both are right.** That is worth stating plainly, because it is the first time
two of this framework's own gates have pulled against each other. The refinement
is necessary; the relation really is 13.6 percent wrong at 150 K; and the water
under dispute is one part in ten million of the column.

The repair is neither to widen the floor nor to coarsen the grid. It is for the
closure to stop evaluating a relation outside its domain:
`saturation_humidity_within_domain` returns the humidity where the relation
holds and exactly zero below it, together with the mask of layers it was
evaluated at, and the humidity modules publish
`saturated_layer_temperature` -- NaN where they declined -- and bind
`EQ-CLAUSIUS-CLAPEYRON-SATURATION` to that rather than to the whole profile.

Declining is a statement about the closure rather than a repair of the relation,
and its price is measured: the surface temperature moves by 0.01 K at every
resolution from 12 to 40 layers, and the precipitable water by about a tenth of
a percent. The real upper atmosphere is not in equilibrium with ice anyway --
its few parts per million of water are set by transport and by methane
oxidation, neither of which a single column contains -- so a saturation closure
has nothing to say up there and now says nothing.

In the prognostic water module the same change has a consequence worth naming:
saturation of exactly zero means whatever vapour a too-cold layer holds
condenses out and precipitates. That is a transfer by the ordinary route with
the budget closing, not a clip -- the water leaves as water.

All six greenhouse-grid configurations now report 0 violated and can support a
claim.

## The vapour splitting was additive where the temperature splitting was sequential

The over-draw recorded below -- two modules jointly taking 3.79 kg m^-2 from a
layer holding 3.17 -- was not a hard coupling problem. It was an asymmetry
inside convection, and finding it took reading the module rather than the
symptom.

Convection has always adjusted the *provisional* temperature: the profile after
radiation and after the water module's latent release. Its own docstring says
why -- adjusting the pre-heating profile would stabilise a column that
condensation is about to destabilise. The vapour was not treated the same way.
`_vapour_transport` mixed the committed, start-of-step vapour, so convection's
tendency and the water module's condensation were *additive* rather than
sequential, and both charged the same layer for water only one of them could
have. Each was affordable alone. Together they were not.

The repair is composition, not a limiter:

    w_provisional = w_committed + vapour_phase_change_rate * dt
    contribution  = mix(w_provisional) - w_provisional

so the committed vapour is exactly `mix(w_provisional)`. That is non-negative
because the water module's provisional profile is (it condenses at most the
excess over saturation) and because mixing to a uniform specific humidity cannot
take a non-negative group negative. The group total is untouched. Nothing is
clipped and no demand is scaled back: the over-draw becomes *structurally
impossible* rather than detected.

The water module publishes `vapour_phase_change_rate` for this, as a rate,
exactly as it already published `latent_heating_rate` for the temperature half
of the same composition. And convection now declares the vapour, the phase-change
rate and the latent heating as required inputs when vapour transport is on. It
was reading all three through `.get`, which is the same contract defect the
radiation module had when it derived cloud optics from three fields it declared
nowhere: optional in general is not optional in the configuration that switches
it on.

The rate is converted over the interval the water module published it for rather
than over convection's own step, so the owner's commit reproduces the mixed
profile exactly rather than only where the two modules happen to share a
schedule.

**What it unlocked.** Every evaporation coefficient strong enough to approach the
observed evaporation rate used to over-draw a layer and halt. They now run:

| C_E | T_s | latent | sensible | Bowen | mm/day | precipitable water | imbalance |
|---|---|---|---|---|---|---|---|
| 2.0e-4 (configured) | 291.12 | 6.5 | 62.0 | 9.54 | 0.224 | 3.1 | 0.000 |
| 1.0e-3 | 294.72 | 36.1 | 60.5 | 1.68 | 1.247 | 11.0 | 0.571 |
| 2.0e-3 | 297.85 | 73.1 | 53.9 | 0.74 | 2.526 | 23.0 | -0.786 |
| observed | 288 | ~82 | ~21 | ~0.26 | 2.7 | 25 | ~0.7 |

20000 hourly steps each, with the top-of-atmosphere imbalance in the last
column: the first row is equilibrated and the last still carries most of a watt,
so read 297.85 as a close waypoint rather than a converged equilibrium.

At 2.0e-3 the column carries an Earth-like hydrological cycle for the first
time: evaporation 2.53 mm/day against 2.7, precipitable water 23.0 against 25,
latent flux 73 against 82. And the surface temperature goes the wrong way, from
291.1 to 297.8 -- a bias of +9.9 K where the too-dry column had +3.1.

That is the compensating pair, measured from the other side. Repairing the dry
half alone does exactly what was predicted of it: the process quantities become
right and the temperature becomes worse. The remaining error is now
concentrated in the half that is still wrong -- a surface sitting far too far
above the air beneath it -- rather than split between two errors of opposite
sign. That is a better place to be even though the headline number is worse.

## A relation is judged where it ran, which is not everywhere it was evaluated

Convection computes its stability criterion for every layer and mixes only
inside the unstable groups, so judging its declared temperature range against
the whole profile let a stratospheric layer flag a relation whose output never
reaches it. It now publishes `convectively_mixed_layer_temperature`, NaN outside
the mixed layers, and binds `EQ-DRY-CONVECTIVE-ADJUSTMENT` and
`EQ-ENTRAINING-LAPSE-RATE` to it. The potential temperature stays bound to the
full state, because that one *is* diagnosed for every layer.

The mask is the union over the run rather than the current step's groups, and
the reason is a limitation of the check worth stating: **it judges a run by its
final state, so a closure that is inactive at the end has no footprint to be
judged on.** Measured on the dry radiative-convective column, the adjustment
mixes on 79 percent of steps and on *none* of the first four hundred -- the
column is dry-stable while it spins up. Masking per step made the identical
configuration resolvable at 2000 steps and unresolvable at 400, which is a fact
about where the run stopped and not about the physics. The union answers the
question actually being asked -- which layers does this relation govern in this
configuration -- and it errs wide, which is the right direction for a check
whose job is to catch a relation evaluated out of range.

Worth recording as a negative result: the first version of this mask *lowered*
coverage, from 23 inside to 22, and the drop was correct. The check it removed
had been passing on a relation that had never run.

## The hydrological cycle is an order of magnitude too weak, and the surface temperature hides it

The most consequential measurement in this file, because the quantity that looks
closest to right is the one being produced by two large errors cancelling.

Set beside the observed global mean, on an equilibrated column with interactive
water, radiative vapour and convective vapour transport:

| | model | Earth | ratio |
|---|---|---|---|
| evaporation | 0.225 mm/day | 2.7 mm/day | 12x too small |
| surface latent heat flux | 6.5 W m^-2 | ~82 W m^-2 | 13x too small |
| surface sensible heat flux | 62.2 W m^-2 | ~21 W m^-2 | 3x too large |
| Bowen ratio | 9.6 | ~0.26 | 37x too large |
| precipitable water | 3.1 kg m^-2 | ~25 kg m^-2 | 8x too small |
| surface minus lowest-air temperature | 6.2 K | ~1-2 K | 4x too large |
| **surface temperature** | **291.2 K** | **288 K** | **+3.2 K** |

Every process quantity is wrong by roughly an order of magnitude and the
temperature is out by three kelvin. That is not a model that is nearly right.

**Why evaporation stops.** The water module has no vertical transport of vapour,
which it declares. What was not measured until now is the consequence: the
evaporated water piles into the lowest layer until that layer saturates --
relative humidity exactly 1.0 -- and evaporation then shuts itself off, because
the bulk formula is driven by the saturation deficit. The cycle throttles to
whatever trickle condensation removes. Switching on convective vapour transport
helps and does not cure: the Bowen ratio falls from 41 to 12, and the free
troposphere above the convective top stays at relative humidity exactly zero.

**Why the surface then runs hot.** The surface can only shed its net radiation
radiatively or as sensible heat, and the bulk sensible flux is a coefficient
times the surface-to-air temperature difference. The coefficient is not the
suspect: Earth's `rho c_p C_H |U|` is about 7 W m^-2 K^-1 and the configured
value is 10. So pushing 62 W m^-2 through that term *requires* a six-kelvin gap,
and in the prescribed-humidity configuration, which moves no water at all and
so carries a latent flux of exactly zero, pushing 90 W m^-2 requires nine.

**What closing the gap is worth.** Stiffening the surface-air coupling until the
gap nearly vanishes, in the prescribed-humidity column that sits 13.7 K warm:

| coefficient (W m^-2 K^-1) | surface | lowest air | gap |
|---|---|---|---|
| 10 (configured) | 301.66 | 292.61 | 9.05 |
| 20 | 299.93 | 294.27 | 5.66 |
| 50 | 298.22 | 295.56 | 2.67 |
| 200 | 297.04 | 296.31 | 0.73 |

**4.6 K of the 13.7**, and not the 9 K the gap itself suggests, because the air
warms as the surface cools and the two converge near 296.7 rather than the
surface falling to meet a fixed air temperature. Arithmetic on the gap would
have overstated it by a factor of two.

**And the hypothesis this refutes.** It is natural to conclude that restoring
the missing evaporation would cool the surface by the amount the missing latent
flux represents. Measured, it does the opposite. Equilibrated, with the
imbalance below 0.07 W m^-2 in every row:

| evaporation coefficient | surface | latent | precipitable water | evaporation |
|---|---|---|---|---|
| 2.0e-4 (configured) | 291.15 | 6.5 | 3.12 | 0.225 mm/day |
| 5.0e-4 | 291.25 | 16.2 | 3.58 | 0.560 mm/day |
| 1.0e-3 | 294.68 | 36.1 | 10.64 | 1.248 mm/day |

Five times the evaporation makes the surface **3.5 K warmer**, because the
vapour that evaporation puts into the column is a greenhouse gas and its
absorption beats its own evaporative cooling. So the weak hydrological cycle is
not a cooling term this model is missing. It is a pair of large opposing errors
-- a column eight times too dry, which should be far colder, and a surface four
times too decoupled from the air above it, which is far too warm -- that happen
to leave the surface temperature three kelvin out.

That is the shape of a compensating error, and it is worth saying plainly that a
surface temperature within a few kelvin of observation is not evidence for this
column. The evidence is in the process quantities, and they are out by ten.

## Clipping a state variable creates the quantity it was protecting

`WaterModule.advance_state` ended in `np.maximum(initial + total * dt, 0.0)`.
Clipping a negative vapour to zero *creates water*, and it bought nothing: the
conservation reconciliation catches the invented mass anyway, but hundreds of
steps later and as a percentage residual with no indication of where it came
from.

Driving the column toward a realistic evaporation rate produced exactly one
clip, of 0.62 kg m^-2, which surfaced 1475 steps on as a 5.4 percent water
residual and a halted run. The clip fires nowhere in the test suite, so its
entire contribution was to obscure the one case it fired in.

The underlying defect is operator splitting. Convection computes how much vapour
to mix out of a layer from the state at the start of the step, the water module
computes how much to condense out of the same layer from the same state, and
neither sees the other. Each is affordable alone; together they draw 3.79 from a
layer holding 3.17. The state is refused now, at the step that creates it,
naming the layer, the over-draw and the labels responsible:

    the committed vapour tendencies would take a layer below zero, so more water
    was removed than the layer holds. Clipping this to zero would create water
    rather than conserve it. layer 3 holds 3.16686 kg m^-2 and the step draws
    3.78658. Labels this step: convective_vapour_transport,
    evaporation_minus_condensation

That is a complete diagnosis in one line, in place of a percentage residual a
thousand steps downstream. The splitting itself is not fixed and is recorded as
open: the honest repair is for the vapour tendencies to be reconciled against
what the layer holds before they are committed, not for either module to guess.

## Validity domains, enforced

Thirty-four equation records carried a `ValidityDomain` and
`ValidityDomain.contains` was called from the registry's tests and from no
running code. The conditions each relation was established for were declared and
never checked, and five relations had by then been evaluated outside them --
each caught by hand, late, and by accident.

`coupling/validity.py` closes that. Every module's cited equations are judged
against the run's final state, and the driver records the result as a
**mandatory** gate, so a violation invalidates the run for inference rather than
annotating it. That is the harsher reading and the right one; making it advisory
would rebuild the situation it was written to end.

**How a variable is resolved.** Equations declare ranges against physical
symbols and published fields carry symbols too, so the vocabularies already meet
and no lookup table has to be invented and kept in step. Resolution is scoped to
the citing module's own declared inputs and outputs, then falls back to a field
of the same symbol root and identical units -- an equation writes `T` where a
column publishes `T_a`.

Three refusals matter more than the matching:

* **Never across modules.** `alpha` is the evaporating condensate fraction in
  the plume equations and the surface albedo in the radiation fields. A global
  match by symbol would check one against the other's range and report
  confidently on nothing.
* **Never across units.** A range in degC against a field in kelvin differs by
  273.15 and neither says so, so the pair is refused rather than converted.
* **Never when ambiguous.** A module seeing both `T_a` and `T_s` in kelvin gets
  an unresolved report naming both, because nothing here knows which the
  equation meant.

**Unchecked is a third outcome, not a pass.** A check that quietly succeeded
when it could not find its variable would licence exactly the confidence it
fails to earn. Coverage on a cloudy column is 31 of 32 declared ranges
inside, 0 violated, 1 unresolvable, and that one is a stated category error
rather than a gap -- so a clean report is no longer the much weaker statement it
was, though the summary still reports the three outcomes separately.

**Where inference is not enough, the module declares the binding.**
`ModuleContract.equation_variables` maps an equation's symbol to one of the
module's fields or parameters, because only the module knows that the `T` in
Schwarzschild is its layer temperature and not its surface temperature, and
because the bulk-plume parameter is never a field at all. The binding states its
own units rather than inheriting the equation's: a parameter carries none, so
binding `concentration_ppm` to a range declared in `1` would compare parts per
million against a fraction and report confidently -- the same silent mismatch
the checker refuses everywhere else. Where the binding names a field, both units
are compared, so a module that mis-states its own is caught too.

**Three ways to be unresolvable, and they do not mean the same thing.** A
variable can be one the checker could not find, one whose units do not match the
declaration, or one the module *integrates over* rather than evaluating at any
value. The surface's thermal emission is Planck integrated over every
wavelength, so its `lambda` has no value a range could be checked against;
binding it to the surface's grid edges would state a wavelength support for a
quantity that has none left. A module says so with `INTEGRATED_OVER`, which
keeps the check honest without inviting anyone to "close" it.

**A unique inferred match is not a correct one.** The ambiguity refusal only
fires when a second candidate exists to notice. The deep ocean's slab energy
balance was being checked against the *mixed layer's* heat capacity, which was
the only field carrying the symbol `C` in the right units anywhere in that
module's scope -- so inference resolved it and reported confidently on the wrong
reservoir. Nothing in the report looked wrong. Only the module can say which
quantity plays a variable, and where it has not said, a clean check is evidence
about the naming rather than about the physics.

**A known coarseness, and its cure.** The check judges the extent a quantity
*spans*, not the region the equation was actually *applied* to. Convection computes its
lapse-rate target for every layer and only adjusts below the convective top, so
a stratospheric temperature can flag a relation whose output is never used up
there. The check is conservative in that direction, which is the right way for
it to be wrong, but a violation still has to be read rather than obeyed.

A module that wants the narrower statement can now make it: **a NaN means the
relation was not evaluated there.** The latent energy of vapour forced this.
`E = L_v sum_i w_i` was flagged for reaching 225 K against a declared floor of
230, and the column holds *exactly zero* vapour below 230 K -- every temperature
that broke the range multiplied a weight of zero. Publishing the layer
temperature masked to the layers that carry vapour turns the extent into the
region, and the field declares `MissingValuePolicy.MASKED` so the mask is a
statement rather than NaNs smuggled past a default that exists to catch them. An
infinity is still refused, because masking is something a module does on purpose
and an infinity is a defect -- if non-finite simply meant "skip me", a field that
had blown up could pass by being broken.

**It found a sixth violation on its first run.** The bulk-plume column evaluates
the Clausius-Clapeyron saturation formula down to 172 K against a declared floor
of 230 K, and feeds the resulting vapour to the radiative transfer. The
configuration is now `invalid_for_scientific_inference` in consequence. The
stratospheric vapour involved is small, but "small enough not to matter" was an
assumption nobody had checked, which is the whole complaint. The remedies are to
extend the formula to the ice phase, to widen the declared domain if 230 K is
too conservative, or to restrict the configuration -- and until one is chosen
the column cannot support a claim.

The remaining 20 unresolvable ranges are mostly internal quantities never
published as fields -- optical depths, orbital eccentricity, the bulk-plume
parameter. Closing those needs modules to declare which of their quantities
plays each equation variable, which is real work and is in
[`todo.md`](todo.md).

## Known limitations of the current kernel

### The C-shaped humidity closure does not converge

Implemented, measured, and left switched off. Charter 3.5: it is in the tree
because the negative result is worth more than the absence of it.

The constant-mass-flux closure gives one relative humidity for the whole column
because it assumes `epsilon = delta`. Dropping that needs a second relation, and
radiative-convective equilibrium supplies one: the environment subsides at
whatever rate balances its own radiative cooling, so `M = rho |Q| / (Gamma_d -
Gamma)` is determined rather than chosen, and `d_z M = (epsilon - delta) M` then
fixes the detrainment. No new parameter, and Romps makes the same argument.

Diagnosed on a fixed profile it works, and produces the observed shape: 0.96 in
the boundary layer, a minimum of 0.44 near 2 km, climbing to 0.81 by 13 km.

Coupled and run to equilibrium it fails, in two separate ways.

**It does not converge.** The humidity sets the radiative cooling, the cooling
sets the mass flux, and the closure differentiates that mass flux -- so a
layer-scale kink feeds back on itself. Refining the grid makes it worse:

| layers | 10 | 20 | 30 | 40 | 60 |
|---|---|---|---|---|---|
| surface (K), per-layer gradient | 299.9 | 301.5 | 302.4 | 303.4 | 305.9 |
| surface (K), gradient over 2 km | 298.5 | 302.1 | 303.2 | 303.6 | 304.9 |

Taking the gradient over a fixed physical depth rather than layer by layer is
the right thing to do on its own terms -- detrainment belongs to a cloud layer,
not to a model layer -- and it halves the drift without removing it. The
equilibrium profile is a sawtooth rather than a C, alternating between the
zero-detrainment floor of 0.444 and values near saturation, and the fraction of
layers pinned on that floor grows with resolution. VAL-002 is not satisfied and
the configuration is therefore not admissible evidence for anything.

**It buys much less than the humidity sensitivity implies.** Holding a uniform
humidity and varying it, the column moves 13 K between 0.725 and 0.5. The C
shape moves it 3.7 K, because it *raises* the humidity where the water actually
is: 42 percent of the column's vapour sits in the lowest three layers, and the
closure puts those near saturation. It is a free-tropospheric relation being
asked about a boundary layer, where cumulus mass flux is not what sets the
humidity -- a validity-domain violation of exactly the kind this framework
declares and does not check.

**Neither failure is fixable by patching, and the reason closes the path.** The
zero-detrainment floor is not a regime needing a clip: in a consistent bulk
plume `d_z M = e - d` with `d >= 0` gives `d_z ln M <= epsilon`, so `delta >= 0`
is automatic. Negative detrainment is proof that the mass flux and the
entrainment rate are mutually inconsistent, and they are, for an identifiable
reason -- `epsilon = a gamma / (1 - alpha)` is *derived under* `epsilon = delta`,
so using it here imports the constant-mass-flux assumption into the model that
drops it. It shows up as an entrainment profile that is qualitatively inverted:
`epsilon` proportional to `gamma` gives 0.06 per km near the surface rising to
0.51 aloft, where plume similarity and every cloud-resolving diagnosis put
entrainment *largest* low down.

Replacing it needs an independent entrainment law, and there is not one to
have. Romps' general model takes `epsilon` and `delta` from cloud-resolving
output rather than deriving them, and his own direct measurement finds the
bulk-plume `epsilon` underestimates the true entrainment by about a factor of
two.

What that does *not* mean is that entrainment is an ill-posed quantity. The
diagnosed rates rise with model resolution, and it is tempting to read that as
the physics being resolution-dependent; it is not. The bulk-plume diagnostic is
`epsilon = -d_z phi_c / (phi_c - phi_e)`. The numerator is the core's physical
dilution. The denominator is a contrast between two populations whose boundary
is a *definition*, and refining the grid moves it: the cloud's outer,
partially-mixed shell is reclassified from core to environment, `phi_e` moves
towards `phi_c`, and `epsilon` inflates while the mixing is unchanged. Zhao et
al. (2024) frame the same thing from the other side, as the scale-awareness
requirement that subgrid convective transport fall to zero as the flow becomes
explicitly resolved.

So the underlying mixing rate is a well-defined, resolution-independent
quantity. It is only harder to measure than the published spread suggests, and
the grid divergence above is entirely this implementation's own.

It also raises a question this closure has not asked: the relative humidity the
bulk plume predicts is the humidity of the environment *as the bulk plume
partitions it*, which -- by the same argument -- is not obviously the clear-sky
humidity that the radiative transfer wants to consume. Unexamined here.

Taking `epsilon(z)` from simulations would not violate INF-003 -- no temperature
record is involved -- but it costs the thing that justified this closure. The
bulk plume earns its place here by being over-determined: one number fixes the
humidity *and* the lapse rate, so it can be caught being wrong, and it was. A
height-varying `epsilon(z)` and `delta(z)` read off simulations has enough
freedom to reproduce any humidity profile, and a closure that cannot fail is
worth less here than one that is wrong in a legible way.

So the vertical structure is not the seventeen kelvin, and it is not the next
thing to fix either. What the measurement establishes is that a uniform 0.725 is
wrong in both directions at once, and that the error is in the column rather
than in the profile.

### The evaporating fraction has no slack in it

Recorded under charter 3.5 because it closes an escape route rather than opening
one, and because the hypothesis it refutes was this project's own.

The bulk-plume column fails across the whole measured range of the evaporating
condensate fraction, and the entrainment floor sits above where it fails, so the
two literature constraints do not overlap. The proposed escape was that the two
constraints might not be measuring the same quantity: Romps' `alpha` is the
fraction of condensate evaporating *at the height it formed*, whereas a measured
precipitation efficiency counts everything reaching the surface, which also
excludes detrained anvil ice and rain evaporating below cloud base. If so,
`1 - PE` would be an upper bound on `alpha` rather than `alpha`, and a smaller
true value might reopen the window.

It does not. Romps states the identity directly (Les Houches lectures, section
1.11): "alpha the fraction of condensates formed at height z that evaporate into
the environment. For analytic solubility, we will assume that this fraction of
condensates formed at z also evaporates into the environment at height z. Note
that the precipitation efficiency PE is equal to 1 - alpha." The
same-height assumption governs *where* the moisture is deposited, for analytic
solubility. It does not restrict *how much*.

The measurement agrees from the other side. Langhans, Yeo and Romps (2015)
define `PE = FE x SE`, the formation efficiency times the sedimentation
efficiency, and report `0.54 x 0.62 = 0.33` and `0.69 x 0.55 = 0.38`. The
sedimentation efficiency is exactly the rain-evaporating-before-the-surface term
that the hypothesis assumed was excluded: 38 to 45 percent of hydrometeors do
not make it down, and that loss is inside `1 - PE`. Both sides count the same
thing -- condensate returning to vapour -- so `alpha = 1 - PE`, and the
refutation of the global-mean configuration stands on firmer ground than before.

Two things worth keeping from the exercise. Romps' theory gives `alpha < RH` and
therefore `PE > 1 - RH`, a self-consistency check the closure could make and
does not: at `alpha = 0.67` it returns `RH = 0.725`, which passes, but narrowly.
And the same-height assumption *is* wrong in a way that matters elsewhere --
roughly 40 percent of the loss happens during sedimentation, below the formation
height, so the closure moistens too high. That bears on the missing C shape in
the humidity profile rather than on the magnitude of `alpha`.

The numbers are also from simulations of cumulus congestus, not of deep tropical
convection, which is a validity-domain caveat on the measurement itself.

### A single column cannot set its own free-tropospheric humidity

The strongest conclusion in this file, reached twice by different routes, and it
puts a ceiling on the current configuration rather than describing a defect in
it.

The subsiding-branch mechanism is real: air that has not been saturated since it
left the convective outflow keeps the vapour it had then, and drying by
adiabatic descent is what makes the global-mean free troposphere drier than a
convecting one. Diagnosed against an imposed export level it lands on the
observed 0.4 to 0.5. But the export level and the ascending area fraction are
both free, and no published value for the area fraction was locatable.

Removing them does not help, because the two parameter-free limits do not
constrain anything:

| | free-tropospheric RH |
|---|---|
| full local remoistening (Romps' equation, what the model does now) | 0.725 |
| no remoistening at all (pure last saturation) | 0.015 |
| observed | 0.40 to 0.50 |

The dry limit is that dry for a real reason rather than a numerical one:
detrainment concentrates near 13 to 16 km where the saturation humidity is of
order 1e-5, so air subsiding from there arrives essentially bone dry. It is the
textbook result that a troposphere ventilated only from the tropopause would
hold no water. The observation sits between the limits and a bracket spanning
0.015 to 0.725 predicts nothing.

What sets the real value is how much a descending parcel is remoistened on the
way down, which is the ratio of its descent time to the time it takes to travel
laterally away from convection. That is horizontal transport. A single column
does not contain it, and there is no first-principles closure for it -- the same
finding as for the entrainment rate, arrived at independently.

So the seventeen kelvin is not a defect to be fixed inside this architecture.
The free-tropospheric humidity of a global-mean single column is not predictable
from within it, and since that humidity sets both the base state and the
sensitivity, neither is predictable either. Three honest responses, and none of
them is a closure:

* **Prescribe the humidity and drop the feedback claim.** This is what
  `PrescribedHumidityModule` already says of itself -- no feedback claim may be
  registered against it, because the feedback follows from the assumption.
* **Restrict the configuration's validity domain to convecting conditions**,
  where the bulk plume is right and the comparison to a global mean is the
  error. Cheap, but it does not fully rescue the numbers: the column still
  carries about twice the observed tropical precipitable water.
* **Build resolved horizontal transport.** That is a different model, not a
  closure swapped into this one.

### The warm bias and the excess sensitivity are one error

Worth stating separately because it decides how the remaining work is scoped. It
is tempting to read a seventeen-kelvin base-state error and a climate
sensitivity two to three times too large as two problems, and from there to
conclude that no single mechanism accounts for either, so a long tail of small
corrections must. Measured, that is wrong. One lever moves both, and they arrive
together:

| RH | 0.725 | 0.60 | 0.50 | 0.45 | 0.40 | observed |
|---|---|---|---|---|---|---|
| surface (K) | 305.18 | 296.36 | 292.13 | 290.42 | 288.84 | 288 |
| dT_2x (K) | 8.761 | 7.577 | 5.002 | 4.294 | 3.811 | 2.5-4.0 |
| amplification | 6.48x | 5.04x | 3.15x | 2.63x | 2.27x | -- |

The column is too warm *and* too sensitive for the same reason: it holds far too
much water vapour. What makes the error look like it needs many contributions is
the amplification. The humidity error itself is modest -- 0.725 against a
global-mean free troposphere nearer 0.45, a factor of 1.6 -- and it produces a
7.5-fold error in precipitable water and seventeen kelvin, because the vapour
feedback compounds it. A small cause with a large effect is the signature of one
mechanism under strong feedback, not of many summing.

The single mechanism is available and is not a closure that can be swapped in: a
global mean is not a convecting column. The bulk plume is right about what it
describes. The global-mean free troposphere is drier because it is dominated by
subsidence and by where descending air was last saturated, which is horizontal
transport and which no single column contains.

Two qualifications. There is a second, smaller discrepancy: at the humidity that
gives the right temperature the column carries about 14 kg m^-2 of water against
an observed 25, so it reaches 288 K with roughly half the observed vapour and
something else is supplying the difference -- the band calibration, the absent
cloud longwave effect, or the lapse rate. And none of this licenses prescribing
0.45. That is the same move as tuning the evaporating fraction back to 0.3, and
it is calibrating against the temperature record through a different door.

The methodological point is the one that would cost the project. If seventeen
kelvin genuinely required tuning many small principles, the model would stop
being falsifiable, because enough free parameters reproduce any target. INF-003,
the over-determination preference, and the rule that a mechanism enters as a
constitutive relation rather than a fitted forcing all exist to detect exactly
that failure, and the measurement above is what keeps this out of it.

### Clouds cannot repair the global-mean column's warm bias

Recorded under charter 3.5 because it closes a proposal rather than opening one.
The bulk-plume column runs about 17 K too warm, and the obvious remedy is that
clouds cool. Measured, they do not help here, for a reason that is a property of
the configuration rather than of clouds.

The column's surface reflectance defaults to 0.30. That is Earth's *planetary*
albedo, which already contains the clouds; Earth's actual surface reflects about
0.12. So the clouds' shortwave cooling is already in the column, applied at the
ground. The size of the substitution is measurable: the same column at a true
surface albedo of 0.12 with no cloud equilibrates at 338 K rather than 305 K.

Replacing the shortcut honestly -- surface at 0.12, an explicit cloud reflecting
the difference -- leaves the absorbed sunlight where it was and *adds* the
cloud's longwave warming, which the shortcut never contained. At matched
reflection the explicit-cloud column is 306.1 K against the shortcut's 305.2 K:
warmer by 0.9 K, not cooler. The 17 K survives clouds intact.

A cloud *feedback* cannot help either, for a separate reason. A feedback is zero
at the base state by construction; it changes the response to a forcing, not the
equilibrium the column sits at. And the assessed sign is positive -- Sherwood et
al. (2020) give +0.45 W m^-2 K^-1 -- so making clouds respond would amplify the
warming rather than damp it.

None of which says the cloud module is adequate. Its single liquid deck gives a
shortwave effect of -65 W m^-2 against Earth's -47 and a longwave effect of
+5 W m^-2 against Earth's +27, because most of Earth's longwave cloud effect
comes from high cold ice cloud that one mid-level liquid layer cannot represent.
That is a real gap. It is just not the gap that the 17 K is in.


Kept current rather than accumulated: a limitation that has since been removed
is deleted from this list, and where a stage removed one the entry says so.

- Zero-dimensional global mean only; no diurnal, seasonal or latitudinal
  structure, and the geometry module says so in its validity domain.
- Uncertainty is declared and marked UNRESOLVED, not propagated. UNC-005's
  ensemble and Monte Carlo methods have interfaces but no implementations. The
  parameter sweep that decided the bulk-plume closure's fate -- a factor of 2.7
  in the feedback across one assumed parameter -- was run by hand in scratch
  scripts, which is the clearest demonstration of the gap the project has
  produced and was entirely accidental.
- **Every equation record declares a validity domain and nothing checks it.**
  `ValidityDomain.contains` is tested and called from no running code, so a
  configuration that evaluates a relation outside the conditions it was
  calibrated at produces a plausible number and no complaint. That has happened
  three times: the carbon dioxide background opacity behind four claims, the
  vapour band evaluated at a third of its calibration humidity, and the
  evaporating fraction outside the measured range. All three were caught by
  hand, and the second was already outside a domain this framework had written
  down.
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
- Water vapour is radiatively inert **in the claim-bearing configurations**. The
  spectral relation and the coupling exist and are exercised, but every
  configuration that carries a claim still runs with the fixed background
  opacity it was calibrated against, so no registered claim depends on a
  water-vapour feedback and none is in scope.
- The vapour feedback the coupled column produces is **not safely converged in
  the layer grid**. Amplification runs 1.10x, 1.22x, 1.88x, 1.28x and 1.30x at
  10, 20, 30, 40 and 60 layers, with the spin-up converged to three decimals at
  every one of them. The two finest grids agree to within two percent, so there
  is probably a limit near 1.3x -- but a 45 percent excursion at 30 layers, in a
  sequence that is otherwise smooth, is an unexplained resolution-sensitive
  switch and not something to quote around. Under VAL-002 the number stays
  recorded and unquoted until the excursion is understood.
- Above the convecting layer the column is bone dry, and that is where the
  feedback goes missing: the emergent scheme delivers about +0.5 W m^-2 K^-1
  against +1.251 for the fixed-relative-humidity control, which is roughly the
  fraction of the column the mixing reaches. Complete mixing without entrainment
  drives the top of the group past saturation, condensation clips it there, and
  nothing moistens anything above it.
- Relative humidity above the convecting layer is exactly zero, because
  evaporation is the only source of water and convection the only transport. A
  real atmosphere moistens its middle troposphere by detrainment from deep
  convection and by large-scale ascent, neither of which exists here.
- The zero-buoyancy plume predicts a stratification but not the height
  convection reaches, so the convective top is set by where the profile becomes
  stable with respect to the target rather than by a criterion of its own. It
  lands at 186 hPa and 12.1 km, which is reasonable, and an undilute-parcel
  criterion would change nothing -- see [`todo.md`](todo.md).
- **The bulk-plume humidity configuration is refuted at its measured
  parameters.** With the evaporating condensate fraction taken from published
  precipitation efficiencies rather than assumed, the column runs 13 to 26 K too
  warm across the whole measured range, and the entrainment floor requires a
  value above where it has already failed, so the two literature constraints do
  not overlap. The closure describes a convecting column correctly -- its
  predicted humidity matches radiative-convective simulations -- and a
  global-mean column is not one. The parameter is not being tuned back down.
- The +1.29 W m^-2 K^-1 feedback previously reported for this pathway was
  obtained at an evaporating fraction a factor of two below the entrainment
  floor, and is superseded rather than deleted (charter 3.5).
- The bulk-plume humidity closure once ran in a column that relaxed to the dry
  adiabat at 9.68 K/km while its humidity assumed 6.33, so the model
  contradicted its own closure. That is fixed -- `convective_target=
  "entraining_plume"` makes both use the same parameter -- and the fix cost the
  lapse-rate prediction its apparent accuracy, since 6.33 was the closure
  evaluated on a profile it disagreed with and 5.5 is what it gives when made to
  live with its own answer, against an observed 6.5.
- The evaporating condensate fraction is assumed rather than measured, over the
  range 0.1 to 0.5, and that range spans 8 K of surface temperature and a factor
  of three in the vapour feedback. No result from a plume-humidity
  configuration may be quoted without it.
- The plume column reaches 288.1 K and a 33.1 K greenhouse effect while carrying
  a third of the observed water and a lapse rate 3.2 K/km too steep. **The
  agreement is cancellation between two large errors of opposite sign**, and is
  recorded here so that it is never mistaken for validation.
- **The vapour relation is being evaluated outside the conditions it was
  calibrated at.** The fit was performed on a column carrying 23.2 kg m^-2 of
  precipitable water; the coupled column settles at 5.3, and its surface at
  280.9 K against an observed 288. This is the same defect pattern that the
  carbon dioxide band had twice, in a third place, and it is why no claim has
  been registered against the vapour pathway.
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
  refuted.** Solar alone explains 3 percent of the observed variance and adding
  carbon dioxide takes that to 71. Adding deep-ocean heat uptake *lowers* it to
  65, and adding volcanic aerosol lowers it again, because the column
  undershoots the observed amplitude by a factor of three and both mechanisms
  reduce amplitude further. Nothing has been *supported* by observation, and
  nothing can be until a prospective comparison exists -- a reconstruction pass
  earns QUALIFIED at most.
- **Those numbers replace a set that was wrong, and the deep-ocean one changed
  sign.** Under the mis-calibrated band the deep reservoir raised the variance
  explained from 61 to 76 percent, and that was read as the timing mattering
  more than the amplitude. It was an artefact of an over-sensitive column. The
  version 1.0.0 claims are superseded rather than deleted (charter 3.5).
- The greenhouse configuration's base state is about 4 K too cool, which is
  roughly the size of the water-vapour greenhouse effect it does not have. It
  was previously 27 K too warm, and that was a calibration defect rather than
  physics; see "The base state is a prediction".
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
  evaluation period. Nothing blocks one; the obstacle is that a prediction
  registered now is one the model is expected to fail, and waiting for it to be
  good enough to pass is how prospective prediction turns into fitting.

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

## Interface negotiation existed and nothing called it

VER-006 asks the framework to negotiate interface compatibility before
initialisation, and `core/versioning.negotiate` did exactly that -- for a
provider and a consumer handed to it by its own unit tests. No running code
called it. In the moist column eight provided interfaces carried forty-four
declared semantic keys, seven consumer declarations named what they required,
and nothing compared the two. It was the validity-domain situation again: the
conditions were written down, and a run was never checked against them.

Two things were missing, and the second was the more important.

**The call.** `Coupler.validate_configuration` now resolves every participating
consumer's `required_interfaces` against the participating providers, after
each module validates its own configuration and before any module initialises.
A missing provider, a version outside the declared range or a disagreement in
meaning fails the run there, with the configuration-validation gate recorded as
failed. Every resolved pair is written to `RunManifest.negotiated_interfaces`.
Across all seventeen reference configurations every declared pair resolves
directly -- airless 2, dry 5, moist 7, cloud 5, deep ocean 8.

**Something to compare against.** A consumer that states an interface id and a
version range can only be checked for identity and version, and neither changes
when a provider turns its column upside down inside a compatible minor release
-- the case VER-007 exists for. So a `RequiredInterface` also states what the
consumer *assumes* the numbers mean: `layer_ordering=index_0_is_the_lowest_layer`,
`flux_sign_convention=positive_upward_away_from_surface`,
`toa_irradiance_basis=global_mean_over_planetary_surface`. Only load-bearing keys
belong there. Copying a provider's whole semantics block would make every
provider change a consumer change and say nothing about which assumptions the
consumer's code is built on.

Writing those nine declarations caught a guess of mine: I first gave the
radiation module an `insolation_units` expectation on `orbital_geometry`, and the
provider publishes `toa_irradiance_units`. The negotiation refused the run on
its first execution, which is the check doing its job on the person who wrote
it.

**One declared way across a gap (VER-008).** Refusing is the right default and
cannot be the only option, because the way round a refusal that cannot be
satisfied is to loosen the refusal. A `MigrationAdapter`, declared on the kernel
rather than on a module because it reconciles two modules and belongs to
neither, may bridge a version or semantics gap, and must state the approximation
it introduces, the information it loses and the conservation residual it leaves.
Each is required; "none" is a valid declaration that has to be *made*, and an
empty string is refused -- a silent omission and a declared absence must not
look alike. A nonzero `residual_magnitude` described as "none" is refused too.

An adapter also names the semantic keys it `translates`. One written for a
pressure-units change cannot carry a reversed layer ordering that arrived in the
same release: the negotiation names the untranslated key and refuses. No
shipped configuration uses an adapter. The requirement binds to the mechanism
now, so that the first time someone needs one, the way through is a declaration
rather than a change to `negotiate`.

Traceability moves from 178 to 179 of 185.

## A module declares its ranges, and an ensemble is held to them

`evaluation/uncertainty.py` implemented UNC-005's six methods and nothing
reported through them. Every module's `report_uncertainty` had two ways to talk
about a parameter it knew was uncertain -- a number, or the word UNRESOLVED --
and neither fitted what the project actually knew. The evaporating fraction has
a measured range, and it lived in a comment beside a constant. The evaporation
coefficient is *not* independently uncertain, because it is `k_H / c_p` exactly,
and that lived in a docstring. Five sweeps had been run from scratch scripts,
and none of their runs knew it was part of anything.

**Two new states between a number and UNRESOLVED** (`contracts/uncertainty.py`):

| state | meaning | declared now |
|---|---|---|
| `DeclaredRange` | bounded, with its basis: measured, published or asserted | evaporating fraction 0.62-0.81 (measured); bulk-plume parameter 0.1-0.3 (asserted); sensible heat coefficient 6-14.5 W m^-2 K^-1 (asserted) |
| `DerivedParameter` | fixed by another parameter, so not to be perturbed alone | evaporation coefficient follows the sensible heat coefficient |

A report that declares a range must say what distribution it assumes across it,
and one that declares a derivation must state the correlation; both are refused
as "unresolved". "A range, not a distribution" is a valid answer and is the one
every module gives, which is also why no interval here states a coverage.

**An ensemble is checked against the declarations before it runs.**
`run_ensemble` reads them from the configuration and refuses a member that
perturbs a derived parameter or leaves a declared range. The two refusals are
the project's own mistakes, replayed as tests: the 0.1-0.5 sweep of the
evaporating fraction that was actually run, which excluded the measured value,
and a sweep of the evaporation coefficient alone, which is how the fifty-fold
Bowen-ratio inconsistency would come back.

The first version of that check read each *member's* configuration, and the
derived-parameter refusal did not fire. Passing the evaporation coefficient
explicitly is precisely what stops a configuration deriving it, so the
perturbed configuration no longer declared the rule its own override broke --
asking it whether the perturbation was allowed got the answer yes. The
declarations are read from the unperturbed configuration now, which is why an
ensemble's `build()` with no overrides must produce its baseline.

**Every member knows it is one.** A run built by `run_ensemble` records the
design's hash, the method, the basis, its label and the member count in
`RunManifest.ensemble_membership`, and the `EnsembleResult` lists each member's
configuration hash. Membership is kept out of the configuration hash: a member
computes exactly what the same configuration computes alone, and the test
checks that they hash alike. A standalone run configured outside a declared
range is allowed, and the manifest says so in a warning.

**The first result through it, in kelvin.** The sensible heat coefficient fixes
both surface fluxes, so it is the first suspect for a warm surface. Across its
declared range, equilibrated at 20000 steps:

| k_H (W m^-2 K^-1) | error | evaporation (mm/day) | PW (kg m^-2) | Bowen |
|---|---:|---:|---:|---:|
| 6.0 | +8.42 K | 3.65 | 27.0 | 0.343 |
| 10.0 (as configured) | +6.74 K | 4.05 | 28.3 | 0.339 |
| 14.5 | +6.92 K | 4.47 | 32.8 | 0.308 |

The error never gets closer than +6.7 K, and the weakly coupled end is the
warmest. Whatever the warm bias is, no value of this parameter anyone has
declared plausible removes it. The 0.19 K between the nominal and the strong
end is inside the limit cycle's time-mean convergence of about 0.15 K, so it is
not evidence of a minimum inside the range, and nothing here claims one. The
interval is stated as the full range with no coverage, because the members were
chosen rather than drawn.

The ocean exchange coefficient was left undeclared. Its value sits beside the
published multi-model mean, but the spread across models has not been looked
up, and writing one down from memory is the mistake this section began with.

## The headline error was an eight-layer number

Rerunning the plume refutation through `across` gave +5.9 to +11.5 K across the
measured evaporating fraction, against a recorded +13 to +26 K. It read as a
refutation half the size of the recorded one, and nearly went into the record
as that. The recorded table had been made at 20 layers; the reference builder
defaults to 8. Before explaining the difference by the physics changes since, it
was checked -- and the grid explained it, and more.

**The full water-cycle column**, 16 bands, hourly steps, equilibrated:

| layers | 8 | 12 | 20 | 30 | 40 | 60 |
|---|---:|---:|---:|---:|---:|---:|
| error (K) | +6.73 | +10.46 | +11.12 | +12.07 | +11.98 | +12.46 |
| evaporation (mm/day, obs 2.7) | 4.05 | 4.50 | 4.65 | 4.86 | 4.86 | 4.92 |
| precipitable water (kg m^-2, obs 25) | 28.3 | 35.2 | 35.2 | 37.4 | 36.7 | 37.5 |
| Bowen ratio (obs 0.26) | 0.339 | 0.253 | 0.196 | 0.155 | 0.142 | 0.125 |

**The error is about +12 K, not +6.7.** Thirty layers and finer agree to half a
kelvin, and the suite now reports the 40-layer figure on every run. The process
quantities move with it: the column is wetter and evaporates more, and the
Bowen ratio that sat reassuringly near the observed 0.26 at 8 and 12 layers
passes straight through it and keeps falling. An agreement that holds at one
resolution and not at the next is a coincidence of the grid.

No layer study had been run on this column. One existed -- the vapour
feedback's, 10 to 60 layers -- but on an older configuration carrying 5 kg m^-2
of water, and it was never repeated after the transport and the derived
exchange coefficient changed what the column is. So the quoted number rested on
a convergence result about a different column.

**The plume column** at the measured fraction does not settle at any
resolution above 8 layers in 24000 steps:

| layers | 8 | 12 | 20 | 30 | 40 | 60 |
|---|---:|---:|---:|---:|---:|---:|
| error (K) | +7.40 | > +11.3 | > +15.8 | > +19.0 | > +20.9 | > +23.1 |
| imbalance (W m^-2) | 0.03 | 0.10 | 0.50 | 0.90 | 1.21 | 1.63 |

It is refuted at every resolution, so the verdict stands and the recorded +13
to +26 K is not contradicted. Whether it has a limit is open (`docs/todo.md`).

What this changes in the sweeps just added: nothing about their conclusions,
both of which strengthen with resolution. The exchange coefficient cannot
remove a +12 K bias any more than a +6.7 K one, and the plume column is further
from observation, not closer. Both still run at 8 layers because that is what
they can afford, and both now say so in every label they print.

## Thirty layers, and the calibration moved with it

The default column resolution is now 30 layers, and it is a constant beside the
fitted wing decays, `CALIBRATION_N_LAYERS`, rather than a default of its own.
The vapour bands had been fitted on a 10-layer prescribed profile while every
configuration ran 8, which was the same defect the background opacity once had:
two conditions that should be one, separated by nothing but the fact that they
lived in different files.

**What the fits depend on.** The same calibration profile, at increasing layer
count, with the old constants:

| layers | CO2 forcing, fixed background | outgoing longwave (target 238.2) | CO2 forcing, vapour background |
|---|---:|---:|---:|
| 10 | 3.698 | 238.20 | 3.700 |
| 30 | 3.700 | 235.69 | 3.640 |
| 60 | 3.701 | 235.12 | 3.626 |

The fixed-background CO2 fit does not depend on it and was left alone. The two
vapour-dependent fits do, because the vapour a prescribed humidity places in a
layer depends on where the midpoints fall. Refit jointly at 30 layers against
the original published targets -- the 238.2 W m^-2 the column absorbs, and 3.70
W m^-2 per doubling -- the vapour wing decay moves from 3970 to 3827 m^-1 and the
masked CO2 one from 1272.5 to 1274.9. The method was checked first by rerunning
it at 10 layers, where it returned 3971 and 1272.7.

**What it did to the headline**, the full water-cycle column:

| layers | 8 | 20 | 30 (default) | 40 | 60 |
|---|---:|---:|---:|---:|---:|
| error (K) | +5.58 | +9.70 | **+9.81** | +10.62 | +10.63 |
| Bowen ratio (obs 0.26) | 0.363 | 0.211 | 0.177 | 0.158 | 0.142 |

**+9.8 K at the default, +10.6 K at the limit.** The narrower vapour wing lets
the column radiate more, and it is 2.3 K cooler at 30 layers than the same
resolution was under the 10-layer fit. The shape is not a smooth convergence:
20 and 30 agree, 40 and 60 agree, and there is a 0.8 K step between. The test
holds the step's size rather than pretending it is not there.

Everything else that was measured moved the way this suggests and changed no
conclusion. The dry column, which has no vapour, reads -3.50 K at 30 layers
(-4.50 at 8, -3.42 at 40). Across the declared sensible heat coefficient range
the full column is +9.8 to +11.0 K, the weakly coupled end still the warmest.
The plume column, pinned at 8 layers, is +4.1 to +9.6 K across the measured
evaporating fraction and so still refuted at every point, and at 8 layers only
because above that it does not settle.

**The one exception is declared.** The historical replay takes two-day steps
with forward Euler, and at 30 layers its thinnest layers relax faster than that:
the run goes unstable, a layer reaches -63 K, and the conservation gate refuses
it -- which is the gate working. It stays at 8 layers as `REPLAY_N_LAYERS`,
because the calibration it depends on is the resolution-independent one. Its
response has not been checked for layer convergence.

The suite now takes about twenty minutes, most of it in the equilibrated runs
that report the error.

### The step was a scatter, and the scatter is the convective top

The 0.8 K "step" between 30 and 40 layers came from three points. Six more
dissolved it:

| layers | 20 | 25 | 30 | 32 | 34 | 36 | 40 | 45 | 60 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| error (K) | +9.70 | +9.91 | +9.81 | +10.51 | +10.35 | +10.21 | +10.62 | +10.17 | +10.63 |

No plateaus and no trend: from 25 layers up the error scatters over 0.8 K
around +10.3. The time-means are stable to 0.03 K against the tail length, so
the scatter is in the column, not in the averaging. It follows one thing:

| layers | mean top of mixing | error |
|---|---:|---:|
| 25 | 462 hPa | +9.89 K |
| 30 | 446 hPa | +9.81 K |
| 45 | 415 hPa | +10.16 K |
| 36 | 413 hPa | +10.23 K |
| 32 | 400 hPa | +10.53 K |
| 40 | 391 hPa | +10.64 K |

Deeper mixing, warmer surface, about 0.015 K per hPa, with the outgoing
longwave at 240.0 W m^-2 throughout -- so it is the column's vertical structure
moving under a fixed energy balance, not the balance moving. The adjustment
mixes whole layers or none, and the top of what it mixes is a layer edge; where
the edges fall relative to the true neutral level differs with every grid.

So the headline is **+10.3 +/- 0.4 K**, the scatter being the resolution
uncertainty, and the 30-layer default's +9.8 sits at the low edge of it by
chance of where its edges fall. The test asserts the band and the mechanism --
of 30 and 40 layers, the one that mixes deeper is the warmer -- rather than a
limit that does not exist.

## A relaxation timescale, chosen by convergence, and not chosen

The layer-edge scatter looked like a consequence of instantaneous mixing, and
the suggestion was to let convergence choose a relaxation timescale: relax the
adjustment over `tau`, and take the `tau` at which the answer stops depending on
the grid. That is a numerical criterion that never consults the observed
temperature, so it was a legitimate way to set a parameter -- provided the rule
was fixed before the numbers were seen. It was: the protocol was committed
first, and the result appended after.

The relaxed adjustment applies `1 - exp(-dt / tau)` of the instantaneous one
each step (`EQ-RELAXED-CONVECTIVE-ADJUSTMENT`). One scalar per step keeps every
group's adjustment summing to zero, so energy and water close per step exactly
as before; `tau -> 0` is the old scheme, checked against it run separately; and
the applied fraction matches the closed form to 1e-9 at the first mixing step.

The sweep, errors in K:

| tau | 25 | 30 | 32 | 36 | 40 | 45 | S(tau) | span at 30 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 (instantaneous) | +9.89 | +9.81 | +10.53 | +10.23 | +10.64 | +10.16 | 0.84 | 3.55 |
| 1 h | +12.07 | +12.62 | +12.23 | +11.69 | +11.47 | +12.10 | 1.15 | 4.89 |
| 2 h | +12.13 | +13.04 | +12.84 | +13.16 | +13.59 | +12.89 | 1.46 | 3.93 |
| 4 h | +13.73 | +13.76 | +14.62 | +14.33 | +14.10 | +14.85 | 1.12 | 3.52 |
| 8 h | +14.22 | +14.74 | +14.71 | +15.48 | +16.22 | +16.42 | 2.21 | 1.80 |
| 16 h | +15.96 | +16.48 | +16.48 | +16.55 | +17.47 | +17.88 | 1.92 | 0.79 |

**Nothing met the rule, so nothing was adopted.** The published 2 h gives a
scatter of 1.46 K against the 0.3 K threshold, and every `tau` in range is worse
than the instantaneous scheme's 0.84. The negative result has three parts:

* The scatter is spatial. The unstable groups are whole layers at every `tau`,
  and smoothing the mixing in time does nothing to where its top lies -- which
  the equation record said before the sweep, as an approximation.
* The limit cycle *is* temporal: it falls from 3.55 K to 0.79 K as `tau` grows.
  The two symptoms of the all-or-nothing closure had looked like one problem,
  and they are two.
* The closure choice is worth 2 to 8 K of surface temperature inside the
  published range, most plausibly because a partly adjusted column carries a standing
  instability (not yet measured). That is now the largest model-form uncertainty
  measured here.

Charter 3.5 keeps it in the record: the relaxation stays in the code as a
declared option at zero, with its tests, and this section and the protocol's
result section are the evidence for not using it.

