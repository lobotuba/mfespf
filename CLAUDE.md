# Working in this repository

MFESPF is a falsifiability framework first and a climate model second. The
specifications in [`specifications/`](specifications/) are normative and govern
every disagreement — charter v0.2 plus Amendment 1, and the v1.0 Interface
Requirements Specification (197 numbered requirements, 189 mandatory).

Read [`docs/architecture.md`](docs/architecture.md) before changing anything in
`coupling/` or `registries/`. [`docs/todo.md`](docs/todo.md) is the open-work
list: an item goes there only when something measured says it matters, and it
carries that measurement with it.

## Non-negotiables

These are enforced in code, and the enforcement is the point of the project. Do
not weaken a check to make a test pass — if a check fires, either the code under
it is wrong or the check's reasoning needs to change, and the second case needs
saying out loud.

- **One owner per prognostic state.** Other modules submit labeled tendencies;
  the owner commits. There is no method on `Module` that writes foreign state.
- **Conservation closes or the run is invalid.** Never widen a tolerance to get
  a green budget. `assert_conservation=False` is for tests that are *about* a
  failure, never for making one go away.
- **Numerical sources stay separate from physical ones.** Anything a limiter,
  filter or clip does to a conserved quantity is a numerical source.
- **No fixed spectra.** Spectral fields travel with their grid. Regridding
  conserves the integral and logs its residual. There is no shortwave/longwave
  split in the radiative transfer -- only radiation per band, with solar and
  thermal sources. Keep it that way; the SW/LW split is a diagnostic formed
  after the solve.
- **Emissivity follows from albedo band by band** (Kirchhoff). Do not add a
  free-standing emissivity parameter to a spectral surface; if a configuration
  wants a reflective-but-black surface, it says so spectrally.
- **A process that moves matter must move its energy in the same interval.**
  The water module owns the latent energy for exactly this reason. Do not add a
  phase change whose energy is applied separately from its mass.
- **A closure that redistributes within a module must sum to zero.** Convective
  adjustment conserves column enthalpy by construction, and the test asserts it
  per step. Do not replace it with a scheme that only conserves on average.
- **Cloud shortwave and longwave effects are published separately, before any
  net.** They are of opposite sign and comparable size; a net alone cannot tell
  a right model from a doubly wrong one (SPC-010, CLD-001).
- **A term that is not physics gets a name, not a hiding place.** A prescribed
  cloud's water enters the budget as `prescribed_cloud_water_source`, the same
  way a limiter's energy enters as a numerical source.
- **A second storage reservoir is a second module.** Heat crossing between the
  mixed layer and the deep ocean is a transfer both parties record. Do not fold
  a second temperature into a surface module — it would hide the exchange from
  the auditor and weld the reservoir to one surface implementation.
- **A forcing delivers its declared energy at any time step.** Forcings
  integrate themselves exactly over the interval. Do not add one that samples at
  the step start, and do not give `ImpulseForcing` a fixed flux — the energy is
  the invariant, the sub-step shape is not (TIM-007).
- **A second timescale is claimed from a residual, not from a second mode.**
  Granting any fit an extra mode always lowers its residual; report two only
  when the improvement is large and the modes are separated.
- **Timescales carry observable, subsystem, mechanism and forcing regime, and
  one longer than the record is marked extrapolated.** `single_index()` raises
  rather than averaging measurements that disagree (TIM-008).
- **An experiment's threshold is frozen before the run.** `ExperimentProtocol`
  hashes the declaration; editing it invalidates the evaluation. Do not add a
  path that evaluates an unfrozen protocol (EXP-001).
- **A dataset offered as evidence declares its provenance class.** Role and
  provenance answer different questions. A generated record cannot support an
  empirical claim, and the evaluator caps the status rather than trusting the
  reader to notice (INF-002).
- **Resampling reports its lost frequency content.** Do not add a path that
  resamples without producing a report, and do not average instantaneous samples
  or interpolate interval means — the convention is carried for a reason
  (charter §11.4).
- **A run cannot read past a record's coverage.** No held last value, no
  extrapolation. The failure is the feature (DAT-007).
- **An uncertainty column declares its kind.** A standard error and a 95 percent
  half-width differ by a factor of two and the numbers do not say which. A
  series carrying uncertainties without an `UncertaintyKind` is refused
  (UNC-002).
- **Loaders put every record on one declared epoch** (`EPOCH`, 1850-01-01) and
  timestamp interval means at interval *midpoints*. A series timestamped at
  interval starts resamples half an interval early, which survives every check
  that does not look at a calendar.
- **Anomalies are kelvin, not degC.** A difference carries no 273.15 offset, and
  declaring one `degC` makes the unit registry add it.
- **Verify a file's checksum before parsing it, not after.** Reading first means
  the values reached the caller before anyone knew they were the wrong ones.
- **The operating mode governs what a run may be called.** A prescribed-history
  or reconstruction run cannot be reported as a hindcast or a prediction
  (TIM-006, VAL-007).
- **A failure is decisive; a pass is worth what the mode licenses.** Do not add
  a path that returns `SUPPORTED_WITHIN_SCOPE` for a run that used inputs from
  the period it was judged on — `OperatingMode.status_for_result` is where that
  decision belongs.
- **A hypothesis declares what its experiment cannot separate it from.**
  Required, and inside the frozen hash: shortening the list afterwards
  strengthens the claim as surely as loosening the threshold (HYP-008, HYP-009).
- **An experiment's threshold comes from the data, never from the model's
  error.** The solar-only threshold is the record's own mean 95% half-width over
  the compared period, and the compared period is fixed by the records' overlap
  before anything runs.
- **A mechanism enters through a constitutive relation, a tendency, a boundary
  or a source — never as a fitted forcing added to a budget.** CO₂ changes an
  absorption coefficient; the transfer solver is unchanged and unaware.
- **Calibrate against something other than the thing being predicted.** The CO₂
  band's wing decay is fitted to a published forcing per doubling, never to a
  temperature record (INF-003, EQR-004).
- **A grid must resolve what the physics depends on.** The logarithmic forcing
  comes from the band wings; a grid that cannot resolve them gives linear growth
  or complete saturation instead, and `greenhouse_grid` refuses too few bands.
- **A calibration is only valid at the conditions it was performed at, so the
  calibration and the runs share one constant.** The background opacity and the
  band region and resolution live beside `CO2_WING_DECAY_PER_METRE` and are
  `greenhouse_grid`'s defaults. They were once separate — the fit used one
  background and every configuration passed another — and the runs used a
  forcing 1.75x the calibrated one for four claims before a volcanic pulse
  exposed it. Do not add a call site that passes its own value for either.
- **Check convergence again after changing the physics, not only after changing
  the grid.** Narrowing the band made a previously converged grid inadequate
  without anything about the grid changing.
- **Convergence applies to what the model produces, not only to what it is
  given.** The spectral forcing was converged; the vapour feedback the coupled
  column produced was not, and the spin-up being converged to three decimals at
  every resolution said nothing about it. Take the study far enough to see the
  shape: stopping at 10, 20 and 30 layers showed amplification climbing 1.10x,
  1.22x, 1.88x and read as divergence, and 40 and 60 layers then landed at 1.28x
  and 1.30x, making the 30-layer point a spike rather than a trend. An emergent
  equilibrium quantity gets its own resolution study, with enough points to tell
  a limit from an excursion, before it is quoted (VAL-002).
- **Prefer a closure that predicts two observables from one parameter.** A
  prescribed relative humidity cannot be wrong; the bulk-plume closure that
  replaced it fixes the humidity and the lapse rate from the same number, so
  both can be checked and one of them promptly caught the convective adjustment
  still targeting a dry adiabat. Over-determination is what makes a closure
  falsifiable rather than merely tunable.
- **Saturation is over the phase that is present.** Below the triple point the
  vapour is in equilibrium with ice, so the latent heat is that of sublimation
  and the heat capacity is the ice one; the latent heat also varies with
  temperature by Kirchhoff's relation. Using the liquid branch at constant L all
  the way down was five times too high at 172 K and put thirteen percent more
  water in the column than ice would hold. Both branches anchor at the triple
  point, so nothing has to be blended.
- **A physics change invalidates every calibration downstream of it, and the
  refit is part of the change.** Correcting the saturation formula moved the CO2
  forcing per doubling to 3.81 W m^-2 and the outgoing longwave to 247, both
  outside their fitted targets. Refit against the original published targets
  rather than widening a tolerance; here it moved the vapour wing decay by ten
  percent and left the CO2 one where it was, which is itself the evidence that
  only one of them was wrong.
- **Validity domains are enforced, as a mandatory gate.**
  `coupling/validity.py` judges every module's cited equations against the run's
  final state and the driver invalidates a run that steps outside one. It was
  written after five relations had been evaluated outside their conditions and
  caught by hand; it found a sixth within minutes. Do not downgrade it to
  advisory -- that rebuilds the situation it ended. Resolution is scoped to the
  citing module's own fields and refuses ambiguity, cross-unit comparison and
  cross-module symbol matches rather than guessing, because `alpha` is an
  evaporating fraction in one module and an albedo in another. Coverage is 9 of
  30 ranges, so a pass is a weaker statement than a failure, and the summary
  reports unresolvable separately for exactly that reason.
- **A declared range is only as honest as its bounds — asserting one is not a
  substitute for looking the number up.** The evaporating condensate fraction was
  assumed at 0.3 over 0.1 to 0.5, on the reasoning that a stated range beats a
  fitted number. The reasoning was right and the bounds were wrong: published
  precipitation efficiencies put it near two thirds, so the range excluded the
  answer instead of bracketing it, and every sensitivity study run across it was
  measuring the wrong interval.
- **Read the declared failure modes you wrote down.** The implied entrainment
  rate disagreeing with cloud-resolving simulations by a factor of two sat in the
  equation registry as a known failure mode for two commits. It was the symptom
  of the wrong parameter the whole time. An over-determined closure only helps if
  somebody looks at the second prediction.
- **An assumed parameter's range is part of the result, not a footnote to it.**
  The evaporating condensate fraction is asserted over 0.1 to 0.5, and across
  that range the column moves 8 K and its feedback changes by a factor of three.
  Quoting the midpoint would be false precision, and *choosing* the midpoint
  because it reproduced the observed temperature would be calibrating against
  the temperature record by the back door. Sweep the range and report it.
- **Evaluate a closure on a profile it agrees with, or not at all.** The
  entraining lapse rate looked accurate at 6.33 K/km against an observed 6.5 --
  measured in a column that was relaxing to a dry adiabat, which the closure
  disagreed with. Made self-consistent it gives 5.5, and the apparent agreement
  was an artefact of the profile it was evaluated on.
- **Adding a temperature dependence can break a null control that used to
  work.** Freezing the vapour isolated the feedback while the lapse rate did not
  depend on humidity. The entraining target depends on saturation humidity, so
  the frozen run keeps the full lapse-rate response and the control stops being
  null -- visibly, at 12.4 W m^-2 K^-1 against a Planck response of 3.2. Re-check
  every control after adding a coupling, not only the runs it was controlling.
- **A base state that looks right is not evidence until you know why.** The
  plume column lands at 288.1 K with a 33.1 K greenhouse effect while carrying a
  third of the observed water and a lapse rate 3.2 K/km too steep -- two large
  errors of opposite sign cancelling. Check the pieces, not the total.
- **A mechanism that moves a conserved quantity between two modules has to be
  declared by both of them.** Mixing the surface into the convective adjustment
  balanced its own arithmetic exactly and still failed reconciliation at 3e-4,
  because an internal redistribution had quietly become a boundary crossing.
  Publish the flux and declare the transfer on each side; the reconciler exists
  to catch exactly this and it did.
- **Check whether two quantities differ before building on the difference.** A
  research item was filed on the premise that Romps' evaporating fraction and a
  measured precipitation efficiency count different things, and that the gap
  might rescue a refuted configuration. Romps defines `PE = 1 - alpha` in one
  sentence, and the measurement's sedimentation efficiency is precisely the term
  the premise assumed it excluded. The hypothesis cost nothing to check and
  would have cost a great deal to build on.
- **A quantity that varies with resolution may be a moving definition rather
  than moving physics.** Diagnosed entrainment rates climb with model
  resolution, which reads as the physics being scale-dependent. It is not: the
  bulk-plume rate divides the core's dilution by a core-to-environment contrast,
  and refining the grid reclassifies the cloud's partially-mixed shell as
  environment, shrinking the denominator. Before inheriting someone else's
  resolution dependence as an excuse for your own, check whether their
  diagnostic's boundary moved.
- **Bracket a missing closure with its parameter-free limits before building
  it.** The global-mean free-tropospheric humidity lies between full local
  remoistening (0.725) and none at all (0.015). Both limits are free of
  parameters and the bracket is useless, which establishes in one calculation
  that the quantity is not obtainable from a single column -- and it costs an
  afternoon rather than an architecture.
- **Do not import a relation derived under the assumption you are dropping.**
  `a = epsilon (1 - alpha) / gamma` is defined where `epsilon = delta`. Using it
  to set `epsilon` in a variable-mass-flux column carried the constant-mass-flux
  assumption straight back in, and showed up as an entrainment profile inverted
  against every cloud-resolving diagnosis. The zero-detrainment clipping that
  followed was not a regime to model -- `delta >= 0` is automatic in a
  consistent plume -- but a signal of that inconsistency.
- **Generalising a closure can cost the falsifiability that justified it.** The
  bulk plume is here because one number fixes the humidity and the lapse rate
  together. Height-varying rates read off simulations would fit any profile, so
  the general form is more realistic and worth less. Prefer the version that can
  be caught being wrong.
- **A closure that differentiates a field it also controls will amplify grid
  noise.** The profile humidity takes the gradient of a mass flux derived from
  the radiative cooling that the humidity itself sets, so a layer-scale kink
  feeds back on itself and the surface temperature drifts upward with resolution
  instead of settling. Differentiating over a fixed physical depth rather than
  per layer is both physically right and half the cure; the rest is the
  clipping. Run the resolution study on any closure with a derivative in it
  before believing a single-resolution result (VAL-002).
- **A relation extrapolated past its own domain can move a result the wrong
  way.** The bulk plume describes a free troposphere and was applied to the
  ground, where it returns 0.96 against an observed 0.8 -- and since 42 percent
  of the column's water is in the lowest three layers, giving the column its
  C shape made it wetter overall. The error was not in the closure.
- **The cloudy column is the clear column plus the cloud.** Build the cloudy
  solve by adding the cloud to the base optical depths, never by passing the
  cloud's alone -- that overwrites them and deletes the water vapour from under
  the cloud, which reverses the sign of the longwave cloud effect. It survived
  because carbon dioxide travels by a different route that both solves receive,
  and because no test had put a cloud and an active vapour continuum in the same
  column: the cloud tests ran dry and the vapour work ran clear. Where two
  couplings have never been switched on together, neither one's tests cover the
  pair.
- **A planetary albedo on a surface is a cloud you have already counted.** The
  column's 0.30 is worth about 33 K and is doing the clouds' job at the ground,
  so adding real clouds without dropping it counts their cooling twice.
- **Optical properties describe a layer, not a contributor.** Cloud droplets and
  volcanic aerosol scatter in the same layers; combine them with
  `combine_scatterers`, which weights asymmetry by how much each one scatters.
- **A right forcing does not make a right response.** Check both separately: the
  volcanic forcing matched a published value while the modelled cooling was half
  again too large, and only measuring the feedback parameter found it — a
  calibration defect four claims had already been reported against. A short
  pulse of known amplitude measures a response where a slow ramp cannot, because
  a ramp lets a forcing error and a sensitivity error cancel.
- **A null control that is not exactly null is a bug.** The constant-irradiance
  control gives zero modelled range. Anything else means something other than
  the mechanism is moving the model, and it contaminates the forced run too.
- **Interval-mean cells are centred on their sample times.** `cell_edges` puts
  the boundaries halfway between samples, so resampling a series onto its own
  clock is the identity. A zero-order hold is a different reconstruction and is
  not interchangeable with it.
- **Calibration data is never independent validation.** The dataset registry
  enforces it per claim, at variable granularity.
- **Withheld data is unreadable before the claim freezes.** Do not add a
  convenience accessor.
- **A failure criterion is mandatory on every claim.** A claim that cannot fail
  cannot be registered.

## Conventions

- SI internally. Non-SI units convert only at declared boundaries, and the
  conversion is logged. Unknown unit symbols raise — register them in
  `core/units.py` with an explicit dimension and scale.
- Every new equation gets a registry record *before* a module cites it:
  category, mathematical form, assumptions, approximations, validity domain,
  limiting cases, and at least one correctness test.
- Every module contract needs a cost declaration and at least one null,
  limiting-case or disabled-coupling test. The validator refuses it otherwise.
- Override `Module.parameters()` with everything that determines the module's
  numbers. A parameter missing from it makes two different runs hash alike.
- Physics notation wins over PEP 8 naming: `T`, `C`, `S`, `E` are fine (ruff's
  N803/N806 are disabled for this reason).

## Tests

Every test binds to requirement IDs:

```python
@pytest.mark.req("CON-004", "KER-005")
@pytest.mark.stage(1)
def test_matched_transfers_reconcile_to_zero(): ...
```

Citing an ID that does not exist in `requirements/` fails
`tools/traceability.py --check`. Layout:

- `tests/unit/` — one component, no coupled run; analytic identities belong
  here, evaluated directly against the solver rather than through a coupled run
- `tests/integration/` — a coupled configuration end to end
- `tests/acceptance/` — the IRS §22 acceptance criteria, including the failure
  detection ones

Prefer a test that checks physics against an **independently derived** result —
a closed form, a limiting case, an analytic identity — over one that checks the
code against itself.

## Before committing

```bash
python -m ruff check src tests tools && python -m ruff format --check src tests tools
python -m mypy
python -m pytest -q
python tools/traceability.py            # regenerate docs/traceability.md
```

If you amended a specification, regenerate the requirement registry first:
`python tools/extract_requirements.py` (verification method, status and evidence
fields are preserved across regeneration — they are the human record).

## Adding a validation-ladder stage

Stage 7 is unimplemented. Follow charter §20 rather than writing the
physics first: define the claim, derive the equation, declare its validity
domain, implement the simplest valid form, verify units and limiting cases,
verify conservation, verify convergence, compare against an independent
process-level observation, couple it, re-run the gates, run a null control,
quantify uncertainty, record the claim status. Preserve negative results —
charter §3.5 requires failed mechanisms and abandoned formulations to stay in
the record when they affect interpretation.
