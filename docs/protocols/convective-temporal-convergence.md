# Protocol: convective temporal and plume-top convergence (F4)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The design is the
user's (2026-10-13), recorded here with its implementation fixed.

**It follows from** `fall-path-closure-candidate.md`. F4 failed there: the RH
profile was a changed equilibrium at 254-288 hPa. Duty cycle and CAPE fell
steeply with the host step, while precipitation and RH 300-700 converged.

**Constraints:**

- F0 stays not qualified until this protocol resolves F4.
- Interactive-radiation F0 stays blocked.
- F4's criterion is not loosened.
- Plume physics and timestep handling are never changed in the same
  experiment.

## Question

Why does reducing the host timestep change the equilibrium RH specifically
around the plume's termination, while bulk precipitation and 300-700 hPa RH
are nearly converged?

## Hypotheses, tested in sequence

**H1, temporal coupling.** Trigger, CAPE removal, plume integration and
thermodynamic response occur at the host timestep. A common internal
convective timestep should bring the host timesteps to the same state.

**Read before this was frozen:**

- `convect` limits its substep only by the half-mass rule. In F0 that
  usually lets one substep span the whole host step.
- So the closure M_b = CAPE_d/(τ_c Γ_M) is evaluated once, from the
  start-of-step CAPE, and applied for Δt.
- With τ_c = 3600 s that is an explicit step at Δt/τ_c = 1, 0.5 and 0.25
  for 3600, 1800 and 900 s.
- **The fraction of CAPE an event removes is predicted to depend on Δt.**
  That is the user's f_C.

**H2, discrete plume-top placement.** The physical plume is similar, but
small timestep-dependent changes move the termination between adjacent
layers, and with it where detrainment, water and heat are deposited.

## Common setup (identical to F0)

- **F0's configuration:**
  - carried cloud water, c0 = 0.045;
  - the rain law at ke = 1e-6, "rate_law_only";
  - `detrained_state="plume"`;
  - recycling diagnosed;
  - the provenance tracer on.
- **R4's state and the common frozen radiation.** It is replayed unchanged.
  Its forcing is a rate, so the energy it delivers over each host step is the
  same in every case.
- **20 days, read over days 10-20.**
- **Radiation is never subcycled.**

## T0: reproduce F0, with event diagnostics

**No physics changes.** Host steps 3600, 1800 and 900 s. The diagnostics
added are outputs only, and must leave the state bit-identical (G2).

**Per step:**

- CAPE_d before convection: the published `dilute_cape`;
- **CAPE_d after convection:** the plume relaunched on the state convection
  leaves. It is computed and published, never fed back;
- M_b requested, the closure's first substep, and M_b applied, the step mean;
- the terminal layer and the plume-top pressure;
- **the continuous plume top.** The termination variable is the buoyancy
  work integrated from the level of free convection, X:
  - with X_k > 0 at the last layer before termination and X_(k+1) ≤ 0 at
    the terminal layer, λ = X_k/(X_k − X_(k+1));
  - p_top = p_k + λ(p_(k+1) − p_k);
  - published as a diagnostic, deposition unchanged;
- **detrained mass by layer,** ∫M_b·D_k dt, and **the plume's mass flux by
  interface,** ∫M_b·M_k dt, which is the subsidence mass;
- what already exists, unchanged:
  - detrained vapour, condensate and dry static energy by layer;
  - subsidence of vapour and of dry static energy;
  - precipitation formed and evaporated by layer;
  - every labelled tendency, from which the 254-288 hPa layers' temperature
    and vapour changes are taken by operator.

**Per event** (a maximal run of steps with M_b > 0):

- its start and duration;
- CAPE_before (at its first step) and CAPE_after (at its last);
- ΔCAPE_event, and f_C = ΔCAPE_event / CAPE_before;
- 𝓜 = ∫M_b dt;
- 𝓓_k = ∫M_b D_k dt;
- the plume tops it reached.

**Reported per host step:**

- the distribution of f_C, its mean and median;
- the distribution of terminal layers and of continuous p_top;
- daily 𝓜 and 𝓓(p);
- F0's quantities.

## T1: a common convective substep (decisive for H1)

**Δt_conv = 300 s for all three host steps:** 12, 6 and 3 substeps.

**The implementation, frozen:**

- **`convective_substep` (s)** on `mass_flux.convect`, and the convection and
  configuration option of that name.
- **Within a host step,** convection integrates in windows of Δt_conv.
- **Every substep ends at the next window boundary at the latest,** and
  sooner only if the existing half-mass rule requires it.
- **Every substep, as now:**
  1. relaunch the plume into the current state;
  2. re-evaluate CAPE_d and the trigger;
  3. close M_b;
  4. integrate the plume;
  5. apply entrainment, detrainment and compensating subsidence;
  6. apply the fall-path treatment, so convection's own condensation and
     evaporation are inside its operator;
  7. update the temperature and vapour.
- **This ordering is identical in every case.**
- **When CAPE_d falls below the trigger,** convection ends for the rest of the
  host step, as now. Nothing else changes the state within convection's call,
  so a later window would find the same CAPE.
- **Everything else stays at the host step:** the water module's end-of-step
  adjustment, the dry adjustment, diffusion, surface exchange and the frozen
  radiation.
- **Off, the option is today's code, bit for bit.**

**H1 supported if** G3 passes. F4 was then host-timestep coupling, not a
failure of the fall-path closure.

## T2: internal-step convergence (if T1 passes G3)

**Host step 1800 s, with Δt_conv = 600, 300 and 150 s.**

**H1 confirmed if** G4 passes: the 300 s result is not just another attractor.

## T3: plume-top quantization (only if T1 fails G3)

**From T1's runs,** the distributions of continuous p_top are compared with
those of the terminal layer.

**H2 is implicated if:**

- the three host steps' continuous p_top agree by F4's rule, applied to its
  late mean and median;
- and their terminal-layer distributions differ.

## T4: fractional terminal deposition

**It is frozen in its own protocol, and only if T3 implicates H2.** The
user's design is recorded so that protocol starts from it:

- **The weights.** Terminal detrainment is split between the two bracketing
  layers, (1 − λ)·D and λ·D.
- **The quantities split.** Mass, total water, enthalpy, the provenance
  classes and condensate, all with the same weights.
- **What is never interpolated.** RH and humidity are not interpolated.

## Gates, frozen

**G1, conservation:**

- the conservation gate passes;
- the worst relative residual is ≤ 1e-12;
- **provenance:** the successor F6a, F6b and F6c frozen in
  `fall-path-closure-candidate.md`.

**G2, option-off identity:**

- with the diagnostics and the substep option off, the code reproduces the
  previous code bit for bit;
- T0's runs reproduce F0's results to every printed digit.

**G3, host-step independence (T1).**

- **The requirement:** the three host-step runs at Δt_conv = 300 s pass F4
  exactly as frozen. A6's rule on the late means, then the classification;
  F4 passes if no quantity is a changed equilibrium.
- **The quantities:** F4's list, including the RH profile from 100 to
  850 hPa.

**G4, internal convergence (T2),** by F4's rule and classification across
Δt_conv = 600, 300 and 150 s:

- the RH profile;
- CAPE_d;
- daily 𝓜;
- 𝓓(p), as its largest layer difference;
- precipitation;
- FE, SE and PE;
- daily subsidence mass;
- terminal detrainment.

**A changed equilibrium in any of them fails G4.**

**G5, plume-top convergence,** by the same rule:

- the late mean and median of the continuous p_top;
- the vertical detrainment distribution 𝓓(p).

It is applied across the host steps in T1, and across the internal steps in
T2. RH can converge by coincidence while terminal transport does not.

**G6, F0's invariants:**

- R_recycle ≤ 0.0042;
- evaporation never exceeds its wet-bulb limit (by construction, checked);
- F6c's molecule-once identities;
- no negative states.

**Not a gate:**

- **RH = 0.71.**
- **The duty cycle.** If event frequency changes while 𝓜, 𝓓 and the state
  converge, that is the numerical representation of intermittency.

**Secondary, reported, not gated:** FE, SE and PE by the pressure of first
condensation, for T0, T1 and T2. If these move around 254-288 hPa, the
timestep is changing the water's path, not only its timing.

## Decision tree, frozen

| outcome | what follows |
|---|---|
| **T1 passes G3** | Run T2. If G4 and G5 pass, the internal step is chosen from the converged regime, F0's qualification is rerun, and interactive radiation follows when the stratospheric configuration permits. |
| **T1 fails, T3 implicates H2** | Freeze T4 in its own protocol. |
| **T1 fails, and the continuous tops also differ** | No repair of the deposition. Find which of CAPE → M_b → B(z) → p_top stays step-dependent. |
| **T1 and T2 converge the state but not the event timing** | Acceptable, if the physical-state gates pass. |
