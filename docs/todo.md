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

### Constrain the evaporating condensate fraction

`EVAPORATING_CONDENSATE_FRACTION` is asserted over 0.1 to 0.5. It is the single
largest source of uncertainty in the plume pathway and the reason no number from
it is quotable.

Measured across that range in the self-consistent column, the feedback runs
+0.72 to +1.91 W m^-2 K^-1 and the surface temperature spans eleven kelvin. The
midpoint, +1.29, agrees to within four percent with the prescribed-humidity
control and with published water-vapour-plus-lapse-rate estimates -- which is
why narrowing the range is now the highest-value item here rather than a tidy-up.

It must not be fitted to a temperature record (INF-003), and the temptation is
concrete: the upper end gives dT_2x of 7 K, so the model's own sensitivity would
narrow it considerably. That is calibrating against the thing being predicted.
Precipitation efficiency is directly diagnosable from cloud-resolving
simulations and is reported in the literature, so this is a reading task before
it is a coding one.

### The humidity has no vertical structure

The constant-mass-flux closure gives one relative humidity for the whole column,
so the observed C shape — moist near the surface and near the tropopause,
driest in the middle — is absent. The general form
(`bulk_plume_relative_humidity`) already supports height-varying entrainment and
detrainment and is tested; what is missing is a defensible profile for those
rates. The vapour aloft is where the feedback is most sensitive, so this is the
least trustworthy part of the current answer.

### The implied entrainment rate disagrees with the literature

At `a = 0.2` and `alpha = 0.3` the closure implies 0.06 to 0.16 per kilometre,
several times below the 0.2 to 0.8 that cloud-resolving simulations diagnose for
bulk plumes. Either the parameter pair or the definitions disagree with that
literature. Recorded as a known failure mode on `EQ-BULK-PLUME-HUMIDITY` rather
than adjusted away, and unresolved.

### A tropical closure in a global-mean column

The bulk plume is derived for the tropics, where deep convection sets the
humidity. The extratropical free troposphere is moistened by large-scale ascent
and eddy transport that no bulk plume represents, so the column should be
expected to run dry — and does, at 14.3 kg m^-2 against an observed 25. The
framework cannot estimate the shortfall from within, which is the honest
statement and also the uncomfortable one.

## Ladder

### Stage 7 is unexercised

The freeze records and the information-cutoff field exist and are tested, but no
prediction has been registered before its evaluation period. Nothing else on
this list blocks it.
