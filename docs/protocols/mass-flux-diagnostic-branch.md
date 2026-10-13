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
