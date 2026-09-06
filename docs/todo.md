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

### Validity domains are declared and never enforced

Thirty-four equation records carry a `ValidityDomain`. `ValidityDomain.contains`
is exercised in `tests/unit/test_registries.py` and called from `src/` **zero
times**. The framework declares the conditions each relation is good for and
then never checks whether a run is inside them.

This is not hypothetical. *Evaluating a relation outside its calibration
conditions, and getting a plausible-looking answer,* has now happened three
times:

- the carbon dioxide band's background opacity, which diverged between the fit
  and every configuration, and stood behind four registered claims;
- the vapour band, fitted on a column carrying 23.2 kg m^-2 and evaluated on one
  carrying 8.0, right at the edge of its own declared domain;
- the evaporating condensate fraction, asserted over a range that excluded the
  measured value.

Each was caught by hand, late, and by accident. The second was inside a domain
this framework had already written down. A run-time check — every equation a
running module cites, evaluated against its declared ranges, with violations
reported the way conservation residuals are — would have caught it
automatically. The machinery exists and is wired to nothing.

**Cost.** Small. `contains` already returns `(ok, violations)`; what is missing
is the pass that collects the state each equation's variables refer to, and the
decision about whether a violation halts a run or annotates its manifest.

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

### The evaporating fraction in Romps' own definition

The cheapest genuine research item here, and it may decide whether the closure
above is recoverable.

Romps' `alpha` is the fraction of condensate that evaporates **at the height it
formed**. The measured precipitation efficiencies used to set it — Langhans, Yeo
and Romps (2015) — count the fraction of condensate that reaches the surface,
which also excludes detrained anvil ice and rain evaporating below cloud base.
So `1 - PE` is an *upper bound* on `alpha`, not `alpha`, and the size of the gap
is unknown.

It matters because the column fails at the upper bound and the entrainment
identity fails below it. If Romps' `alpha` is materially smaller than `1 - PE`,
the window between the two constraints may reopen. If it is not, the refutation
stands on firmer ground than it currently does. Either outcome is worth having,
and neither needs code.

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

## Feedbacks the model does not have

### Clouds respond to nothing

Clouds are prescribed, so there is no cloud feedback of any kind, and the
condensate the water module produces precipitates immediately rather than
becoming the cloud. Together with the global-mean column above, this is why no
sensitivity this model produces should be trusted: two of the three feedbacks
that set the real number are absent or refuted, and the third — water vapour —
has no working implementation in a global-mean column.

This is recorded rather than scheduled. Making clouds respond is a larger piece
of work than anything else in this file, and the honest order is to settle what
the column *is* before deciding what should respond inside it.
