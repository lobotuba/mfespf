# Protocol: diagnostic branch, the top-layer drift and the trigger

**Frozen before any run of this branch.** Results are appended below a line
at the end; nothing above it changes afterwards.

**This is a diagnostic deviation authorised by the user (2026-10-13), not a
physics run.** No production conclusion is drawn from it, and nothing in it
changes a module or a default.

## Disposition of the mass-flux protocol, recorded as the user worded it

> Mass-flux protocol implementation completed. Conservation and isolation
> gates A1–A5 passed. A6 failed its frozen endpoint convergence criterion, with
> nonmonotonic RH differences despite convergence of precipitation and mass
> flux. A7 failed because reduced evaporation increased free-tropospheric RH.
> P1–P3 and the sensitivity ensemble remain blocked. All current trajectories
> are additionally invalidated by upper-layer departure from the declared
> RRTMG domain. A bounded diagnostic deviation is authorized to identify the
> upper-layer tendency source and test whether the CAPE trigger controls
> convective duty cycle and free-tropospheric humidity. No production physics
> conclusions will be drawn from the diagnostic branch.

**A6, classified as the user directs.** It failed its frozen endpoint
criterion, and that stands. Its failure mode is to be established here:

- **event-phase sensitivity:** non-monotonic, as intermittent convection
  shifts its events with the step;
- **or secular non-convergence:** the attractor itself changes with the step.

Nothing here makes A6 pass.

**A7 is a real failure.** Reducing evaporation made the free troposphere
wetter. It does not yet show that the CAPE threshold is the cause. Reduced
evaporation also changes the surface temperature, the boundary layer's moist
static energy, the convective frequency and the radiation.

## Part R: the top-layer drift

The 1.1 hPa layer drifts to 36 K below RRTMG's reference in every mass-flux
run. It starts at -27.7 K in R4's state. B0 (`free-troposphere-moisture-budget.md`)
used R4's own convection from the same state for 10000 steps and was not
flagged by the gate.

**R0a, without mass-flux convection.** From R4's exact state, 240 steps:

- **R0a-i:** R4's configuration as it is, the entraining-plume adjustment,
  mass flux off.
- **R0a-ii:** R0a-i with the convection module removed: no adjustment, no
  sensible flux. Surface turbulent exchange is then zero, not held at a
  baseline. That is declared: the column has no other route to hold it, and
  the question is the top layer's radiation, which does not see the
  surface's turbulent flux.

**R0b, with it.** A6's 3600 s run repeated, recording at every step, for the
top three layers:

- every labelled temperature tendency: radiation, mass-flux transport, dry
  adjustment, latent heating, sensible, and any other;
- the layer's vapour, and its labelled vapour tendencies.

**R0c, the top-boundary audit,** from R0a-i and R0b, at every step:

- the top interface's upward and downward longwave;
- the net shortwave at the top;
- the radiative heating of the top three layers;
- the convective enthalpy flux through their interfaces;
- the plume's termination pressure and the mass it detrains there;
- the stratospheric water, 100-1 hPa.

There is no upper-boundary relaxation in this model. That is recorded, not
added.

**The R verdict:** for the top layer, over each run:

- the time-integrated departure from RRTMG's reference;
- each labelled tendency's integral.

The operator whose integral carries the drift past -30 K is named.

| finding | means |
|---|---|
| R0a-i drifts past -30 K too | not the mass-flux scheme: R4's own state drifts |
| R0a-i stays inside and R0b does not, and the top layer's labelled tendencies differ mainly in radiation | mass-flux convection changes what the top layer sees radiatively, for example through stratospheric water. The difference in water is reported |
| R0b's top layer receives a non-radiative tendency | that operator is named, and it reaches the top layer |

**No clamp, relaxation or bound is added to the top layer as a remedy.**

## Part T: the trigger ladder, under one radiation

**The common radiation.** The radiation module of R4's configuration solves
once, on R4's saved tail-mean state. Its solution, every flux and the layer
heating, is then frozen and replayed unchanged on every step of every T case.

- **This is a diagnostic configuration, not a climate run.** It is built in
  the run script by replacing the module's solve with the frozen solution.
  No module code changes.
- **It holds the radiation identical across cases,** so different humidity
  paths cannot produce different radiative paths, or different departures
  from RRTMG's tables.
- **The surface still emits by its own temperature.** The surface module
  owns that.
- **The check:** the replayed heating and surface fluxes are identical,
  value for value, across T0-T3 at every step.

**The cases,** mass-flux convection with everything else as frozen in
`mass-flux-moisture-transport.md`:

- ε = 1 km^-1, δ = ε, the closure, τ_c = 3600 s, α = 0.67;
- the surface-flux law, the initial state (R4's) and dt = 3600 s;
- the replayed radiation.

| case | trigger C_min | purpose |
|---|---:|---|
| T0 | 70 J kg^-1, frozen | the controlled baseline |
| T1 | 35 | threshold sensitivity |
| T2 | 10 | near-continuous activation |
| T3 | 0, convecting whenever a plume solution exists | the most opportunity for circulation drying |

**Length.** Each runs 480 steps (20 days). Means are over days 10-20.

**No case is tuned toward any RH.**

### What each case records, every step

- **The trigger:** dilute CAPE before triggering, and the trigger state.
- **The convection:**
  - cloud-base mass flux;
  - the pressure-integrated environmental subsidence over 300-700 hPa,
    Σ M_e Δp across those interfaces.
- **The environment over 300-700 hPa:** q_v, T, q_sat(T, p) and RH.
- **Tendencies by label:** the radiative temperature tendency, and the
  subsidence (mass-flux) temperature and vapour tendencies.

### Derived quantities

- **Duty cycle:** the share of steps convecting.
- **Events:** the number of contiguous convecting spells, and their timing.
- **Cumulatives:** the convective mass transport Σ M_b dt, the subsidence
  and the precipitation.
- **RH over 300-700:** its mean and its variance over days 10-20.
- **The RH decomposition,** d ln RH = d ln q_v - d ln q_sat, accumulated over
  days 10-20 in the 300-700 hPa layers:
  - d ln q_v by vapour label: mass-flux transport, condensation, diffusion;
  - d ln q_sat by temperature label, through the model's own ∂ln q_sat/∂T:
    radiation, mass-flux transport, latent heat and the adjustment.

  This separates genuine drying (dq_v < 0), subsidence warming
  (dq_sat > 0) and radiative cooling (dq_sat < 0).

## Part A6s: the supporting quantities for A6

**The runs:** A6's three runs (3600, 1800 and 900 s, interactive radiation
as frozen) repeated, recording:

- RH over 300-700: its time-mean over days 0-10 and over days 5-10, and its
  variance;
- the number and timing of trigger events;
- the cumulative mass transport, subsidence and precipitation.

**What they decide:** whether A6 failed by phase shift or by a changed
attractor.

- **A changed attractor:** the day 5-10 mean RH differences behave as the
  endpoint ones did.
- **A phase shift:** they shrink, while the event counts and timings
  differ.

**What they do not do:** make A6 pass.

## The decision rule for Part T, frozen

**The trigger hypothesis is confirmed** if lowering the trigger (T0 → T3)
gives *all* of:

1. a duty cycle that increases monotonically;
2. cumulative subsidence that increases;
3. subsidence warming that offsets the radiative cooling increasingly, in
   the d ln q_sat terms;
4. mean RH over 300-700 that decreases materially, by at least 0.03 from T0
   to T3;
5. a decrease the recorded d ln q_v and d ln q_sat terms account for;
6. water and energy gates that stay closed.

Reaching Romps' 0.81 is not required.

**It is rejected** if convection becomes nearly continuous (duty cycle at
least 0.8 in T3) and cumulative subsidence increases substantially (at
least doubles from T0 to T3), but mean RH over 300-700 stays within 0.02 of
T0's.

**It is implicated but not uniquely identified** if RH falls but the d ln
q_sat terms show it carried mainly by surface or boundary-layer temperature
changes, not by the duty cycle and the subsidence.

**Anything else** is reported as found, without a classification.

## What follows

**If the trigger is confirmed,** it is an emergent humidity closure, and it
gets an explicit redesign under its own protocol, not a tuned threshold.

**If it is rejected,** the single-column mass-flux route has failed its
central drying test, and that is recorded.

**Either way, the user decides.**

---

## Result (appended; nothing above this line changed)

**Corrected while running, before any result was read:** two exactness
defects in the water module's composition with convection, both found by the
gates.

1. **The arithmetic (own + c) - c** left a layer emptied below the
   saturation floor at about -1e-21 kg m^-2. Fixed by keeping the module's
   own condensation apart (commit 9289034).
2. **Committing from summed rates** leaves the same residue on any exactly
   emptied layer. The module now commits such a layer as exactly empty, only
   where the residue is within 1e-12 of that layer's own gross flux. A real
   over-draw is still refused.

Both are option-off neutral. The T runs were rerun after them.

### Part R: the top-layer drift is the convective adjustment, not radiation or the mass flux

The top layer is 1.07 hPa. Integrals over 240 steps (10 days):

| run | radiative tendency | dry convective adjustment | mass-flux transport | dT_ref, start → end | valid |
|---|---:|---:|---:|---|---|
| R0a-i, R4's scheme | -70.9 K | **+70.9 K** | | -27.66 → -27.67 | yes |
| R0b, mass flux | -19.7 K | +11.3 K | 0.0 | -27.66 → **-36.07** | no (RRTMG domain) |

**R0a-ii could not be built as specified:** `moist_radiative_convective`
forces the convection module on. It ran with convection present, and only
its water-composition switches differed: -70.88 and +70.85 K, inside the
domain. It does not decide anything, and that is recorded.

**Verdict as frozen,** third row: R0b's top layer is governed by a
non-radiative tendency, and that operator is named.

- **The operator:** `convection:dry_convective_adjustment`.
- **In R4's configuration** (entraining-plume target), it mixes the layers
  near 1 hPa, where no convection exists, and heats the top layer by about
  7 K/day, exactly balancing its radiative cooling. That is what held R4's
  top layer inside RRTMG's tables.
- **In the mass-flux configuration** the target is the dry adiabat. The
  adjustment supplies only +11.3 K, and the top layer falls toward its own
  radiative equilibrium, 36 K below RRTMG's reference.
- **The drift is not caused by the mass flux.** The plume tops out at
  103-154 hPa, and its transport at the top layer is zero.
- **The stratospheric water** changes little: 26.4 → 25.3 ppmv.

**What it implies, not acted on:**

- R4's upper stratosphere is held by a convective closure acting where it
  has no business.
- The model's own radiative equilibrium at 1 hPa lies outside RRTMG's
  tables.
- `tropical-cloud-field.md` found the tropical stratosphere 12-20 K too cold
  above 70 hPa, which is consistent.

### Part T: the trigger ladder under one radiation

**The common radiation held.** The replayed heating was identical value for
value at every step of every case (checksum -2646.2636 in all four). The
water and energy residuals were 4.9-5.8e-16.

**All four runs are invalid on the mass-flux record's domain.** With the
radiation frozen, cold detrained air drove the plume down to 143 K. That is
a property of the diagnostic configuration, and the cases are compared
within it. Days 10-20:

| case | C_min | duty cycle | events | mass transport, kg m^-2 | subsidence, kg m^-2 Pa | precipitation, kg m^-2 | mean CAPE_d | RH 300-700 (variance) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| T0 | 70 | 0.908 | 23 | 7207 | 2.88e8 | 43.3 | 382 | 0.9372 (7.5e-5) |
| T1 | 35 | 0.883 | 29 | 7204 | 2.88e8 | 43.3 | 372 | 0.9372 (6.6e-5) |
| T2 | 10 | 0.908 | 22 | 7187 | 2.88e8 | 43.2 | 379 | 0.9373 (6.9e-5) |
| T3 | 0 | 1.000 | 1 | 12922 | 4.81e8 | 43.8 | 348 | 0.9351 (8.6e-5) |

**The RH decomposition,** days 10-20 over 300-700 hPa, mass-weighted, the
same in all four cases to 0.01. T0's:

- **d ln q_v:**
  - convective transport **+1.53**;
  - condensation **-1.59**;
  - diffusion +0.05.
- **d ln q_sat:**
  - radiation -1.27;
  - mass-flux transport (its dry part) -0.22;
  - latent heat **+1.45**.
- **Reconstruction:** d ln RH 0.037, against an actual 0.038.

**Verdict as frozen: not confirmed, and not rejected.**

- **Not confirmed:** RH falls only 0.002 from T0 to T3, not the 0.03
  required. The duty cycle does not rise monotonically (0.908, 0.883,
  0.908, 1.000).
- **Not rejected:** T3 is continuous (1.000 ≥ 0.8), but its subsidence rose
  by a factor of 1.67, not the 2 the rule requires.
- **So it is reported as found, without classification.**

**What was found.**

- **Under this radiation, convection was nearly continuous at every
  trigger.** The duty cycle was 0.91 at the frozen 70 J/kg, with mean CAPE
  near 380, so the ladder had little activation left to vary.
- **Continuous convection with 67 percent more subsidence left RH at 0.935.**
- **In every case convective transport moistens the free troposphere:**
  - detrainment of buoyant, saturated plume air;
  - that air is warmer than the environment, so it carries more water than
    the environment can hold.
- **Condensation removes it:** the water module's saturation adjustment, at
  saturation.
- **The plume-subsidence circulation therefore does not dry this column's
  free troposphere.** Its net effect there is to moisten, and saturation
  holds the result.

### Part A6s: the supporting quantities

Interactive radiation, from R4's state, 10 days:

| dt | RH, days 0-10 | RH, days 5-10 | variance, days 5-10 | events | mass transport | subsidence | precipitation |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3600 s | 0.9191 | 0.8941 | 1.04e-5 | 71 | 5484 | 2.151e8 | 36.52 |
| 1800 s | 0.9170 | 0.8945 | 1.20e-5 | 86 | 5735 | 2.252e8 | 36.78 |
| 900 s | 0.9204 | 0.8979 | 1.51e-5 | 139 | 5733 | 2.249e8 | 36.68 |

**Classification as frozen: a changed attractor.**

- **RH:** the day 5-10 mean differences do not shrink, 0.0004 then 0.0034,
  just as the endpoint ones grew.
- **Mass transport, subsidence and precipitation converge** from 1800 s on.
- **The event count rises with the resolution of the step,** 71, 86 and 139,
  so events are shorter and more frequent at short steps.
- **The size:** the RH dependence on the step is 0.003, against the
  0.5-0.8 gap between this column and the observed.

**A6 still fails.**

### What the branch says, and what it does not

- **The top-layer blocker is not the mass flux.** It is the removal of a
  convective adjustment that R4 applied at 1 hPa.
- **The trigger is not the mechanism holding the free troposphere moist.**
  With convection running continuously, the plume's detrainment of warm,
  saturated air moistens the free troposphere faster than its compensating
  subsidence dries it. Condensation then holds the result near saturation.
- **This is the branch's "rejected" outcome in substance,** though not by
  its frozen numbers, since subsidence rose ×1.67, not ×2.
- **The likely physics,** not measured here:
  - the detrainment law δ = ε at 1 km^-1 detrains a great deal of saturated
    air throughout the free troposphere;
  - Romps' (2014) RH = (a + α)/(1 + a) assumes the detrained air carries the
    environment's saturation humidity at the environment's temperature, not
    the warmer plume's.
- **No production conclusion is drawn.** The user decides.

### Addendum (appended after the user's review, 2026-10-13)

**The user's disposition of Part T, recorded as worded:**

> Formally not classified under the preregistered threshold; scientifically
> inconsistent with the trigger hypothesis and supportive of
> detrainment-dominated moistening.

The frozen verdict above is not altered.

**Interpretive erratum, both readings kept.**

- **Reading 1 (above):**
  - "convective transport moistens the free troposphere (+1.53) and
    condensation removes it (-1.59)";
  - so "the plume's detrainment of warm, saturated air moistens faster than
    its compensating subsidence dries it".
- **Reading 2, after re-reading the operator's bookkeeping:**
  - the transport label `convective_vapour_transport` is
    (q - q0)·m + formed;
  - it therefore includes the plume's own condensate formed in each layer;
  - the label `evaporation_minus_condensation` subtracts the same condensate
    (less what evaporates);
  - so the pair contains the plume's in-layer rain-out on both sides;
  - and it does not by itself measure detrainment against subsidence.
  - Also, the operator carries the plume's dry static energy D·s_u with its
    vapour. The "more water than the environment can hold" phrasing holds
    only to second order in the plume's temperature excess.
- **What stands:** the RH numbers, the ladder's frozen verdict and the
  reconstruction (0.037 against 0.038).
- **What is open:** the attribution to detrainment.
- **Where it is resolved:** `detrainment-thermodynamic-state.md` splits the
  operator into subsidence and net detrainment and measures it.
