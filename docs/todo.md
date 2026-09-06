# Open work

Deferred deliberately, with the reason and the measurement that motivates each.
An item earns a place here only if something has been *measured* that says it
matters; a wish without a number belongs in a discussion, not a backlog.

Items are removed when done, not ticked — git history is the record of what
happened. Negative results that closed an item stay in
[`architecture.md`](architecture.md) under charter §3.5.

## Convection and the vapour pathway

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

### The humidity has no vertical structure

The constant-mass-flux closure gives one relative humidity for the whole column,
so the observed C shape — moist near the surface and near the tropopause,
driest in the middle — is absent. The general form
(`bulk_plume_relative_humidity`) already supports height-varying entrainment and
detrainment and is tested; what is missing is a defensible profile for those
rates. The vapour aloft is where the feedback is most sensitive, so this is the
least trustworthy part of the current answer.

### The global-mean column itself

Promoted from a caveat to the main obstacle by measurement. With every parameter
taken from the literature, the bulk-plume closure runs the column 13 to 26 K too
warm, and no value inside the measured range escapes it. The closure is right
about a convecting column and this framework's column is a global mean.

This is the deepest structural assumption in the model and it now has a price
tag: about seventeen kelvin, at the central parameter value. Anything that
addresses it — a second column for the subsiding branch, an explicit
last-saturation treatment, or an honest restriction of the configuration's
validity domain to convecting conditions — is a larger piece of work than
anything else on this list, and everything else in the vapour pathway is
downstream of it.

### A convecting-column closure in a global-mean column

The bulk plume describes a radiative-convective column, where cloud-resolving
simulations put the free-tropospheric humidity near eighty percent. This
framework's column is a global mean, and the global-mean free troposphere is far
drier — dominated by subtropical subsidence and by where air was last saturated,
neither of which a bulk plume contains.

So the closure should be expected to run *moist* here. That is the opposite of
the direction assumed while the evaporating fraction was wrong, and it is a bias
of known sign rather than unknown magnitude, which is the most that can be said
from inside a single column. Resolving it properly needs either a second column
or an explicit representation of the subsiding branch; neither is small.

## Ladder

### Stage 7 is unexercised

The freeze records and the information-cutoff field exist and are tested, but no
prediction has been registered before its evaluation period. Nothing else on
this list blocks it.
