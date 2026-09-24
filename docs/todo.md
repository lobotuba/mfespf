# Open work

Deferred deliberately, with the reason and the measurement that motivates each.
An item earns a place here only if something has been *measured* that says it
matters; a wish without a number belongs in a discussion, not a backlog.

Items are removed when done, not ticked — git history is the record of what
happened. Negative results that closed an item stay in
[`architecture.md`](architecture.md) under charter §3.5.

Ordered by what the evidence says is worth doing first, not by area.

## Framework integrity

These are the items where the framework fails to do the thing it exists to do.
They come before the physics, because the physics is what they are supposed to
be protecting.

### Validity-domain coverage is closed, and the latent heat is not

Coverage is 31 of 32 declared ranges on a cloudy column and 0 violated on every
reference configuration. The one remaining unresolvable is not a gap: the
surface's thermal emission is Planck integrated over every wavelength, so its
`lambda` has no value a range could be checked against, and the module now says
so through `INTEGRATED_OVER` instead of leaving the checker to report that it
could not find the variable. Two unresolvables that mean different things no
longer read alike.

What is left here is the coarseness rather than the coverage. The check still
judges the extent a quantity spans, not the region an equation was applied to,
except where a module masks the field. Convection is the remaining case: it
computes its lapse-rate target for every layer and adjusts only below the
convective top, so a stratospheric temperature can flag a relation whose output
is never used up there. The cure now exists -- publish the target masked to the
adjusted layers -- and has not been applied.

### The saturation is over ice and the latent heat is not

Found while measuring the latent-energy domain, and left open because fixing it
moves the climate and needs a recalibration.

`saturation_vapour_pressure` saturates over the phase that is actually present,
so below the triple point it uses the latent heat of sublimation. The latent
*energy* bookkeeping in `WaterModule` still uses one constant latent heat of
vaporisation at 273.16 K. The column's vapour sits at 259.6 K, below the triple
point, so its condensate should be depositing as ice and releasing the heat of
sublimation.

Measured on the reduced-ocean column: respecting the phase raises the column's
latent energy by **13.5 percent**, against **1.3 percent** for merely letting
L_v vary with temperature by Kirchhoff. So the phase term is ten times the
temperature term, and the module's declared
`latent_heat_constant_approximation` uncertainty is carrying both.

This is an energy the column releases on condensation, so it is not a
bookkeeping detail -- it warms the atmosphere where the condensation happens.
The fix is to make the latent heat phase-aware the way the saturation already
is, and then to recalibrate whatever it moves.

### Uncertainty is declared everywhere and propagated nowhere

UNC-005 asks for bounds, sensitivity derivatives, parameter sweeps, ensembles,
Monte Carlo propagation and model-form comparison. There are interfaces and no
implementations, and every module's `report_uncertainty` returns `unresolved`
for both distributional assumptions and correlations.

The evidence that this matters is accidental and recent. Establishing what the
bulk-plume closure's assumed parameter was worth needed a three-point sweep at
two resolutions, run by hand in scratch scripts, and the answer — a factor of
2.7 in the feedback and eleven kelvin of surface temperature — was the single
most decision-relevant number produced that day. That is UNC-005, done manually
because the framework cannot do it.

### Ten mandatory requirements are unbound but already satisfied

160 of 189 mandatory requirements have a bound acceptance test. About a third of
the remainder need a test rather than an implementation, because the machinery
they describe is already enforced elsewhere:

| requirement | already true because |
|---|---|
| MOD-005 | every module implements `validate_configuration` |
| MOD-007 | the `Module` lifecycle is the contract every ACTIVE module runs under |
| CMP-002 | `ModuleContract.validate` refuses a module with no cost declaration |
| CMP-003 | runs declare a resource class and budget |
| ATM-001 | the state registry enforces single ownership |
| CAU-003 | dependency-graph edges carry a type |
| SYS-002 | `ModuleStatus` distinguishes carve-outs from active modules |
| ACC-005 | manifests, graphs, equation bindings and checksums are produced per run |
| SYS-007 | the claims registry supersedes rather than overwrites |

A framework whose thesis is that compliance should be demonstrated rather than
asserted ought to demonstrate its own. These are sweep tests over the
registries, not new physics.

### The carve-out requirements are counted as undone

BIO-001, ICE-001, LND-001, GRV-001, PAR-001 and SPA-001 through SPA-006 are
requirements for modules that deliberately do not exist. They are recorded as
`not_started`, so the headline figure reads as twenty-nine things left undone
when about ten of them are architectural decisions. SYS-002 requires
distinguishing architectural carve-outs from active modules; the traceability
record does not yet do it, which makes the number mean less than it should.

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
