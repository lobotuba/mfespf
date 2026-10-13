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

---

## Result (appended; nothing above this line changed)

### Erratum to the framing, both readings kept

- **Reading 1 (above, and in the F0 disposition):** F4's failure is in "the
  254-288 hPa plume-top region".
- **Reading 2, measured in T0** (= F0 to every digit):
  - the plume terminates at 57-126 hPa (terminal-layer shares 0.02, 0.18-0.28,
    0.38-0.42, 0.26-0.36 and 0.04-0.06 at 57, 69, 85, 103 and 126 hPa);
  - the continuous top is at 96-103 hPa (median).
  - So 254-288 hPa is the upper troposphere, about 150 hPa below the
    termination, where the plume detrains continuously (δ = ε).
- **The label was mine, from F0's write-up.** It is corrected there too.
- **H2 as written, terminal-layer quantization, is not a candidate cause** of
  a failure 150 hPa below the terminal layers. T3's reading below is recorded
  as frozen anyway.

### G2: option-off identity, and T0 reproduces F0

**Code identity:** with the diagnostics and the substep option off, the code is
bit-for-bit the previous code. The digest matched the ozone-state code.

**Run identity:** T0's runs reproduce F0's late results to every digit at all
three host steps:

| dt | RH 300-700 | residual |
|---:|---|---:|
| 3600 s | 0.7110118708707003 | 5.9e-16 |
| 1800 s | 0.7105124822995943 | 5.7e-16 |
| 900 s | 0.7135748151498488 | 6.4e-16 |

Precipitation, FE, SE, CAPE and M_b reproduce as well. **G2 passes.**

### T0: the event diagnostics

The failure reproduces: the RH profile is a changed equilibrium, 0.029 then
0.042.

**Per event, late period:**

| dt | events | f_C mean | f_C median | 𝓜 per event, kg m^-2 | steps per event |
|---:|---:|---:|---:|---:|---:|
| 3600 s | 95 | 0.80 | 0.985 | 44.1 | 1.09 |
| 1800 s | 99 | 0.81 | 0.987 | 42.7 | 1.09 |
| 900 s | 135 | 0.66 | 0.620 | 31.6 | 1.18 |

- **At 3600 and 1800 s a typical event removes almost all its CAPE in one
  host step.** This is the explicit step at Δt/τ_c = 1 and 0.5 that was read
  from the code before freezing.
- **At 900 s it removes about two-thirds.**
- **Daily 𝓜:** 419, 423 and 426 kg m^-2, phase-sensitive.
- **𝓓(p) and the subsidence-mass profile** converge.

### T1: a common 300 s convective substep

All three runs completed, with residuals ≤ 6.4e-16. They are invalid on the
validity-domain gate only, as F0's were, under the frozen radiation.

| quantity | 3600 s | 1800 s | 900 s | class |
|---|---:|---:|---:|---|
| **RH profile** (largest layer difference) | | 0.025 | 0.010 | **converged** |
| RH 300-700 | 0.7123 | 0.7107 | 0.7114 | converged |
| precipitation, mm/day | 4.171 | 4.155 | 4.152 | converged |
| FE | 0.9953 | 0.9954 | 0.9953 | phase-sensitive |
| SE | 0.9664 | 0.9656 | 0.9655 | converged |
| α_total,unique | 0.0381 | 0.0389 | 0.0391 | converged |
| cloud-base mass flux | 0.0049 | 0.0049 | 0.0049 | phase-sensitive |
| subsidence per day | 1.704e7 | 1.699e7 | 1.704e7 | phase-sensitive |
| R_recycle | 2.8e-7 | 2.8e-7 | 1.9e-7 | phase-sensitive |
| **CAPE_d**, J/kg | 156.9 | 102.0 | 73.8 | **changed equilibrium** |
| **duty cycle** | 0.817 | 0.621 | 0.364 | **changed equilibrium** |

**Integrated transport and the plume top:**

| quantity | 3600 s | 1800 s | 900 s | class |
|---|---:|---:|---:|---|
| daily 𝓜, kg m^-2 | 425.9 | 424.8 | 425.9 | phase-sensitive |
| 𝓓(p) | | | | converged (13.6, then 6.6) |
| subsidence-mass profile | | | | converged |
| continuous p_top, daily mean, hPa | 99.9 | 101.0 | 101.4 | converged |

**Per event:**

| dt | events | f_C mean | f_C median | 𝓜 per event | steps per event |
|---:|---:|---:|---:|---:|---:|
| 3600 s | 45 | 0.98 | 0.993 | 94.6 | 4.4 |
| 1800 s | 95 | 0.74 | 0.987 | 44.7 | 3.1 |
| 900 s | 253 | 0.45 | 0.386 | 16.8 | 1.4 |

### Gates, as frozen

**G1 passes:**

- conservation residuals ≤ 6.4e-16, and the gate passed;
- **F6a:** no failing cells (worst 0.042 of its tolerance);
- **F6b:** no failing steps (worst 0.0013);
- **F6c:** no ordering violations, |PE − FE·SE| ≤ 1.1e-16, inventories
  ≤ 7e-16 of throughput;
- ledger quality ≤ 5.7e-16.

**G3 fails as frozen.**

- **CAPE_d is a changed equilibrium.** Every pairwise difference exceeds twice
  its combined standard error, and the 900-vs-1800 s difference, 28.2, is
  more than half the 1800-vs-3600 s one, 54.9.
- **The duty cycle is a changed equilibrium too.** The frozen text conflicts
  on whether it counts: G3 takes F4's list, which includes it, but the
  protocol also says it is not a gate. Both readings are recorded; CAPE fails
  G3 under either.
- **Everything else in F4's list converges or is phase-sensitive,**
  including the RH profile that failed F4.

**G5 passes:** the continuous p_top and 𝓓(p) both converge across the host
steps.

**G6 passes:**

- R_recycle ≤ 2.8e-7;
- no negative class;
- the molecule-once identities hold;
- the fall-path evaporation stays limited by construction.

### The decision tree, applied as frozen

**T1 fails G3, so T2 is not run.** It was conditional on G3 passing.

**T3, read from T1's runs:**

- **The continuous p_top agrees:** converged, 99.9-101.4 hPa daily mean.
- **The terminal-layer shares differ.** The 103 hPa layer takes 0.37, 0.49
  and 0.53 of the convecting steps; the 69 hPa layer 0.13, 0.09 and 0.08.
- **By the letter of T3, H2 is implicated,** which points to T4.

**What stops that from being acted on:**

- **H2's symptom is gone.** The RH profile converges under T1.
- **The terminal layers are not where F4 failed.**
- **The terminal-layer shares are per convecting step,** and the number of
  convecting steps is itself step-dependent.

T4 is therefore not frozen here. Whether to pursue it is the user's decision.

### What the result says, and what it does not

**1. H1 is supported for the failure F4 actually recorded.**

- A common convective step removes the RH profile's host-step dependence.
- Precipitation, FE, SE, α, 𝓜, 𝓓(p) and the continuous plume top converge
  with it.
- **The RH failure was host-step coupling of the closure,** not the fall-path
  treatment.

**2. What remains step-dependent is CAPE as sampled, and the duty cycle.**

- **The data:** CAPE_d is published from the start of each host step, before
  convection. Under T1, convection reduces CAPE below the trigger within the
  first window or two (per-event f_C 0.98 at 3600 s), and is then idle for
  the rest of the host step. Meanwhile the frozen radiation and the surface
  keep destabilising the column.
- **The explanation, not measured here:** the start-of-step CAPE samples a
  build-up whose length is the host step. That would also make the duty
  cycle, counted per host step, step-dependent while the transport converges.
- **Two things would test it:**
  - CAPE time-averaged over the substeps;
  - convection allowed to re-trigger within a host step as the column
    destabilises. That is not possible while everything but convection acts
    once per host step.
- **That is a further temporal-coupling question:** destabilisation at the
  host step against stabilisation at the substep. It is not resolved by T1's
  design.

**Not concluded:**

- F0 is still not qualified: G3 fails;
- no default changed;
- no substep length chosen;
- T2 not run.
