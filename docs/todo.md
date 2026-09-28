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

### Uncertainty has methods now and nothing uses them

`evaluation/uncertainty.py` implements UNC-005's six methods, UNC-006's
quadrature refusal and UNC-007's declared intervals, and the sensitivity
derivative is checked against a closed form. What is still open is that **no
module and no configuration reports through it**. Every module's
`report_uncertainty` still returns `unresolved` for distributional assumptions
and correlations (UNC-002), and no manifest records that a run was a member of an
ensemble.

That is the half with the evidence behind it. Five sweeps have been run by hand
in scratch scripts -- the bulk-plume parameter, the evaporating fraction, the
timestep, the exchange coefficient, gravity -- and each produced a decision-
relevant number that lives only in a commit message. The machinery exists to stop
the sixth from doing the same, and using it means a run knowing it was a member
of something.

### Nine mandatory requirements are still unbound

176 of 185 have a bound acceptance test, with 4 architectural carve-outs counted
separately rather than folded into the denominator. What is left is three kinds
of thing, and only the first is a backlog:

| | requirements | what it needs |
|---|---|---|
| unimplemented | VER-008, CMP-007, ACC-010 | the machinery they describe |
| a spatial decision | SPA-001, SPA-003, SPA-004, SPA-005, SPA-006 | see below |
| true only when the rest are | ACC-001 | everything above |

Of the three unimplemented, ACC-010 is the odd one: the reference kernel shall
complete its baseline suite within the configured ENGINEERING_WORKSTATION budget.
The preflight estimator and the budget both exist, so this is a matter of
measuring the suite against them rather than building anything.

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
