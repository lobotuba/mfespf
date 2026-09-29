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
- **Check whether a parameter is free before calibrating it.** The evaporation
  coefficient was about to be fitted to the observed evaporation rate. It is not
  free: the same eddies carry heat and moisture, so `k_E = k_H / c_p` exactly,
  and the configured pair was fifty times apart. That single inconsistency was
  the entire Bowen-ratio error, which had been recorded as a separate symptom
  for as long as the water cycle had been measured. A derivation beats a
  calibration, and costs less.
- **A time-mean is not a steady state, and the timestep tells you which you
  have.** At a consistent exchange coefficient the column limit-cycles: the mean
  is converged to 0.15 K across a fourfold step refinement, and the amplitude
  converges to 1.75 K rather than to zero. That makes the oscillation a property
  of the closure, not of the integration. Quote time-means over a full cycle and
  never an endpoint.
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
- **Never clip a conserved state variable.** `np.maximum(vapour, 0.0)` creates
  water. The conservation reconciliation catches the invented mass anyway, so a
  clip buys nothing and costs the diagnosis: one clip of 0.62 kg m^-2 surfaced
  1475 steps later as a 5.4 percent residual with no cause attached. Refuse the
  state at the step that creates it and name what over-drew it.
- **Check the process quantities before believing the temperature.** This column
  sits 3 K from the observed surface temperature while its evaporation is 12x too
  small, its Bowen ratio 37x too large and its precipitable water 8x too small.
  Two large errors of opposite sign -- far too dry, and a surface far too
  decoupled from the air -- nearly cancel at the one number anyone looks at.
- **Strengthening the water cycle warms this column, it does not cool it.** The
  vapour evaporation adds is a greenhouse gas and its absorption beats its own
  evaporative cooling: 5x the evaporation is +3.5 K. Repairing either half of the
  compensating pair alone makes the surface temperature worse first. Expect that
  and do not read it as a regression.
- **Compose split operators, do not cap them.** Convection took the provisional
  temperature and the *stale* vapour, so its mixing and the water module's
  condensation were additive rather than sequential and jointly over-drew a
  layer. The cure is for the later module to act on the profile the earlier one
  published -- `vapour_phase_change_rate` beside `latent_heating_rate` -- which
  makes the over-draw structurally impossible instead of detected. If two
  modules contribute to one state, check that the second sees the first's
  tendency for *every* quantity, not just the one someone thought about.
- **A module that reads a field through `.get` still has to declare it.** Optional
  in general is not optional in the configuration that switches it on. Convection
  read the vapour and the latent heating undeclared, exactly as radiation once
  derived cloud optics from three fields it declared nowhere.
- **A run is not its last step, and a check that reads one is answering a
  different question.** The domain check judged the final state, so a closure
  that ran out of range during spin-up and recovered came out clean and one
  merely idle at the end had no footprint at all -- the convective adjustment
  mixes on 79 percent of steps and on none of the first four hundred, and the
  same configuration was resolvable at 2000 steps and unresolvable at 400. It
  accumulates each relation's extent step by step now. Where a repair needs a
  field widened to the union over the run to leave evidence behind, suspect the
  check rather than the field.
- **Two gates can pull against each other, and neither is the one to weaken.**
  VAL-002 asks for grid refinement; refining the grid makes the top layer colder
  until it leaves Clausius-Clapeyron's measured 170 K floor, and the domain check
  refuses the run. The floor was right, the refinement was right, and the water
  in dispute was one part in ten million. The repair was to stop *evaluating* the
  relation outside its domain, which cost 0.01 K -- not to widen the floor and
  not to coarsen the grid.
- **Declining to evaluate is not the same as clipping.** A humidity closure that
  places no vapour below the relation's floor is making a statement about its own
  scope. In the prognostic module the same rule means a too-cold layer's vapour
  condenses and precipitates -- it leaves as water, by the ordinary route, with
  the budget closing.
- **A declared domain is a guess until someone measures it.** Three optical-depth
  ceilings -- 10, 100 and 1000 -- were all exceeded by the greenhouse
  configuration at tau = 1110, and all three were wrong rather than the run. The
  diffusivity approximation is *most* accurate in the thin limit and in the
  thick one, worst in the middle of the range its floor and ceiling were
  protecting. Measure against an independent reference before moving a bound,
  and before trusting one: a record that declares a limiting case its own domain
  excludes is contradicting itself, and that has now happened three times.
- **A domain says whether a relation may be evaluated; an uncertainty says how
  wrong it is there.** The diffusivity approximation really does understate a
  thin layer's absorptance by 17 percent, and that belongs in the module's
  declared uncertainty, not in a floor that declares an ordinary column invalid.
- **A unique inferred match is not a correct one.** The ambiguity refusal only
  fires when a second candidate exists. The deep ocean's slab balance was
  checked confidently against the mixed layer's heat capacity because it was the
  only field of that symbol and units in scope. Where a module has not bound a
  variable, a clean check is evidence about the naming, not about the physics.
- **A diagnostic added to make something checkable still runs every step.**
  Publishing the optical depth per layer per band cost 80 MB on the historical
  replay and exhausted memory, to carry what is consumed as a min and a max.
  Publish at the resolution the claim needs, not the resolution you have.
- **Heat capacity sets when a column arrives, not where.** A deep ocean makes
  the column look three kelvin cooler at eleven years, and it is not cooler, it
  is unfinished -- the run still carries a 3 W m^-2 top-of-atmosphere imbalance
  and goes on warming for centuries. Shrink the layer until it can equilibrate
  and it lands within 0.003 K of the slab. Check the imbalance before reading
  any temperature difference between configurations with different reservoirs;
  ocean structure belongs to transient response, never to an equilibrium bias.
- **A conserved store must not depend on temperature.** Making the latent heat
  vary with temperature makes the water module's energy store vary with it, and
  then a layer that merely warms changes the store with no water moving and no
  transfer to balance it. That is why constant-L_v was self-consistent and why
  two attempts at this stopped. The way through is to keep a constant reference
  in the store and let the phase enter as a declared transfer -- the surface pays
  the fusion energy when the ice it made aloft lands and melts.
- **A bound that fires on one side of a plateau and passes on the other is
  measuring nothing.** The latent-energy domain ran from 230 K while the
  constant-L_v error is a step at the triple point and then flat at 12 percent
  all the way down. It refused at 230 K and passed at 260 K on an identical
  error, and with the phase term added the true bound is inverted: most accurate
  below the triple point, worst at the hot end the old range admitted.
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
- **Report the error in kelvin from every test that equilibrates a column.**
  `report_bias` collects it and the suite prints them together at the end. The
  number the project is judged on used to be recomputed by hand in a scratch
  script, so it moved between commits without anyone noticing. Currently, at the
  default 30 layers: -33.42 K airless, -3.50 K dry, **+9.8 K** with the full water
  cycle, against 288.0 K. Across 25 to 60 layers the full column scatters from
  +9.8 to +10.6 K around **+10.3 K**, with no trend. Quote it as +10.3 +/- 0.4 K,
  and label any error computed at another resolution with its layer count.
- **The calibration resolution and the run resolution are one constant.**
  `CALIBRATION_N_LAYERS = 30` sits beside the fitted wing decays and is the
  default of every column configuration. The vapour bands were fitted on a
  10-layer profile while every run used 8, and neither was converged; refitting
  at 30 moved the vapour wing decay 3.6 percent and the headline from +12.1 to
  +9.8 K. The historical replay is the declared exception (`REPLAY_N_LAYERS =
  8`): its two-day explicit step is unstable at 30 layers, and the calibration it
  uses is resolution-independent. A new exception needs the same two reasons.
- **An error that scatters with resolution does not have a limit to quote.** The
  full column's error follows the top of the convective mixing -- 0.015 K per
  hPa, at identical outgoing longwave -- and that top snaps to a layer edge,
  because the adjustment mixes whole layers or none. Three resolutions read as
  "a step between two plateaus"; nine showed a scatter, and a diagnostic showed
  its cause. Take a resolution study far enough to tell a limit from a scatter
  before naming either, and when it is a scatter, find what it follows.
- **The headline was quoted for weeks at a resolution nobody had tested.** With
  the old 10-layer calibration the full water-cycle column read +6.73 K at 8
  layers, +11.1 at 20, +12.0 at 30 and 40, +12.5 at 60. The layer study had been run on an older, drier column
  and never repeated after the water cycle changed it. The finding surfaced only
  because a sweep disagreed with a recorded table and the difference had to be
  explained -- the same lesson as the feedback's resolution study, one level up:
  the quantity being quoted gets its own resolution study, and gets it again
  after the physics under it changes (VAL-002).
- **A bias from a run that has not equilibrated is a different quantity, not an
  uncertain one.** `BiasReport.error_k` returns `None` above 0.1 W m^-2 of
  imbalance and reports the imbalance instead. The dry column reads -16.3 K at
  20000 steps while carrying +43.8 W m^-2, and it is unfinished rather than cold.
  Returning a number with a wide error bar would let it be compared against a
  real one.
- **A cheaper form declared and not used is a choice, so say what it would fail
  to produce.** Nine module contracts name a lower-fidelity alternative;
  `ModuleContract.validate` now refuses one that gives no reason (CMP-007). "More
  accurate" is unfalsifiable -- "a single layer has no lapse rate, and the lapse
  rate is what convection acts on" can be checked by running it.
- **What a run stores is its metadata, not its physics.** The retained
  trajectory costs 714 bytes per field per step plus 17.6 per value -- the
  per-field term is forty times the per-value one, because each `Field` carries
  its units, provenance, uncertainty and semantics. The moist column holds 1.7 kB
  of numbers per step out of 37.7 kB. That is the price of being able to say what
  a number means, and it is the memory profile.
- **A preflight that only sums working sets says nothing about a long run.** The
  working set is fixed and the retained trajectory is linear in the step count,
  so the second is the whole of the difference between 2000 steps and 200000.
  Leaving it out made the estimate a factor of twenty low, in the direction that
  admits a run which then exhausts memory -- and that run happened, at 4.3 GB.
- **A spread means nothing until you say how its members were chosen.** The six
  UNC-005 methods are one mechanism and differ only in that, so a sweep across an
  asserted range gets a `FULL_RANGE` with **no** nominal coverage and a note
  saying why, while a Monte Carlo from declared distributions gets a real one.
  Granting a coverage to a sweep and refusing one to a sample are the same
  mistake. Every ensemble also declares its `basis`: 0.1 to 0.5 asserted and 0.1
  to 0.5 measured produce identical numbers.
- **Check a finite difference against its error term, not just against the
  answer.** Landing on the closed form could be luck; missing it by exactly the
  term that was dropped could not. Doing this caught me asserting the wrong sign:
  a secant measures the *mean of the derivative*, so for a fourth root -- concave,
  but with a convex derivative -- it lies *above* the midpoint tangent by
  `h^2 f'''/6`, verified to four significant figures.
- **Quadrature has two preconditions and both fail silently downward.** Root-sum-
  square needs independent contributions that are each one sigma. A 95 percent
  half-width passed in unconverted halves the answer and leaves no trace, which
  is why `UncertaintyKind` travels with the value. Independence must be
  justified, not asserted.
- **An approximation nobody has to declare is one nobody has to measure.**
  Gravity was a constant read from a table at four call sites, and the
  hydrostatic record's bound for it -- "negligible through the troposphere" --
  was true of a troposphere and wrong about this column, whose top is 19 to 27 km
  up. Measured, the constant puts 0.21 percent too little mass in the column and
  0.6 to 0.85 percent too little in the topmost layer, which is the one that sets
  the outgoing longwave. Put a constant behind an interface that states its
  approximation, and the measurement follows.
- **A height-dependent field with no heights refuses; it does not use the
  surface value.** A silent fallback would put the approximation back to being
  applied without anyone declaring it, which is what the interface exists to
  stop. Same rule as declining to evaluate Clausius-Clapeyron below its floor.
- **Check the downstream calibration before recommending against a change, not
  only before making one.** Promoting `NewtonianAltitudeGravity` to the default
  moves the CO2 forcing per doubling by 0.0034 of a +/-0.05 window, so the fit is
  safe. The reason to leave it off is the trade -- 0.037 K against a +6.6 K bias,
  for a renumber of every recorded result -- and that is a different and better
  argument than "it might break something".
- **A carve-out is a domain the kernel does not use; one it uses without the
  declared interface is a gap.** That rule is what keeps the classification from
  being an excuse. The gravity module is as absent from charter 18's candidate
  list as the cryosphere, but the column computes hydrostatic layer mass with
  gravity every step out of a bare constant, so GRV-001 is counted and ICE-001 is
  not. The estimate before anyone applied the rule was "about ten"; the answer is
  four. A carved-out requirement must name the clause that defers it, and
  `tools/traceability.py` refuses one that does not, or that also has a bound
  test.
- **A human influence enters as a constitutive relation, and the test has three
  parts.** That the carbon dioxide module adds nothing to a budget would also be
  true of a module that does nothing, so show it acts -- doubling the
  concentration lowers the outgoing longwave -- and show the budget still closes
  to 1e-16, which a forcing added to a budget could not do.
- **"Already satisfied, just needs a test" is a reading, not a measurement.**
  Ten mandatory requirements were listed in `docs/todo.md` as enforced by
  machinery elsewhere. Nine were. ACC-005's output checksums were empty on every
  run in the tree, because the only path that fills them hashes a file and no
  reference configuration writes one -- a requirement counted as met on the
  strength of a method nothing called. Write the test before believing the table.
- **A sweep with a hole in it makes a claim it has not earned.** The
  self-compliance checks run over every reference configuration, and
  `test_the_sweep_reaches_every_module_implementation` fails and names any module
  in `mfespf.physics` that none of them exercises. A check run against one
  configuration establishes that one configuration complies.
- **The coupler's lifecycle order is not charter 6.3's numbering, on purpose.**
  6.3 numbers `advance_state` before `emit_outputs`; the coupler publishes every
  module's diagnostics in the tendency pass and commits state afterwards, so no
  owner commits while another module's tendency is outstanding. That is the
  operator-splitting invariant, and it is what the convection/water over-draw
  violated. Preserve the semantics, not the numbering.
- **Interface negotiation runs before initialisation, and a consumer states
  what it assumes, not just what version it wants.** `negotiate` had unit tests
  and no caller, so forty-four declared semantic keys in the moist column were
  compared against nothing. Put only load-bearing keys in
  `RequiredInterface.expects`. The one way across a gap is a `MigrationAdapter`
  declaring its approximation, information loss and conservation residual, and
  translating only the keys it names (VER-006, VER-008) -- never a looser
  `negotiate`.
- **A module declares a parameter's range with its basis, or marks it derived,
  or leaves it UNRESOLVED -- and an ensemble is held to the declaration.**
  `run_ensemble` refuses a member that leaves a declared range or perturbs a
  derived parameter, reading the declarations from the *unperturbed*
  configuration: overriding the evaporation coefficient is exactly what stops a
  configuration deriving it, so a member cannot be asked whether its own
  override is allowed. Do not declare a range you have not looked up; the ocean
  exchange coefficient stays UNRESOLVED for that reason.
- **A parameter chosen by a criterion has the criterion frozen first, and "none"
  is an allowed outcome.** The convective relaxation timescale was to be chosen
  by resolution convergence; the rule was committed before the sweep, and the
  sweep found that no timescale in the published range meets it -- relaxation
  worsens the layer-edge scatter and warms the column by 2 to 8 K. The option
  stays at zero, with the result appended to the protocol, not the rule edited.
- **Validity domains are enforced, as a mandatory gate.**
  `coupling/validity.py` judges every module's cited equations against the extent
  each one was evaluated over as the run proceeded, and the driver invalidates a
  run that steps outside one. It was written after five relations had been
  evaluated outside their conditions and caught by hand; it found a sixth within
  minutes. Do not downgrade it to
  advisory -- that rebuilds the situation it ended. Resolution is scoped to the
  citing module's own fields and refuses ambiguity, cross-unit comparison and
  cross-module symbol matches rather than guessing, because `alpha` is an
  evaporating fraction in one module and an albedo in another. Where inference cannot decide, the module binds the symbol to
  one of its own fields or parameters through
  `ModuleContract.equation_variables`, and the binding states its own units -- a
  parameter has none, so inheriting the equation's would compare ppm against a
  fraction. Coverage is 31 of 32 ranges on a cloudy column and 0 violated
  everywhere, so the summary's three outcomes still matter but a pass now means
  something. The check judges the extent a quantity spans, not where the
  equation was applied, unless the module masks the field with NaN to say where
  the relation actually ran -- so read a violation before obeying it. The mask is
  per step, because the accumulation is.
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
