# Protocol: convective-substep convergence (T2, successor gates)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The design is the
user's (2026-10-13), recorded here with its implementation fixed.

**It follows from:**

- `convective-temporal-convergence.md` (F4, T0, T1);
- `convective-cape-sampling.md`, whose disposition authorises T2 as a
  documented deviation from F4's tree. Frozen G3 on T1 remains failed.

**Constraints:**

- F0 stays not qualified until the new qualification run that follows this
  protocol.
- Interactive-radiation F0 stays blocked, and P1-P3 stay blocked.
- No plume physics changes.
- Radiation is never subcycled.

## Question

Is the 300 s convective substep internally converged, or is it one point on
another timestep-dependent sequence?

## Setup

T1's setup, except for the substep:

- F0's configuration, with R4's state and the common frozen radiation rate;
- **host step 1800 s;**
- **Δt_conv = 600, 300 and 150 s** (3, 6 and 12 windows);
- 20 days, read over days 10-20;
- forcing, ordering, initialisation and diagnostics otherwise unchanged.

## Code, output only

Labelled CAPE outputs, published by convection with the mass flux. None is
fed back.

| output | definition |
|---|---|
| `cape_host_start` | CAPE at the start of the host step's first window, on the state convection receives. It is the same value as `dilute_cape`: the disclosed operator-splitting diagnostic. |
| `cape_time_mean` | (1/Δt)·Σ_n CAPE_n·h_n over convection's substeps, CAPE_n at each substep's start, plus the idle remainder × CAPE at the end. Renamed from `dilute_cape_time_mean`, the same definition. |
| `cape_pre_convection` | the mean over the host step's windows of CAPE at each window's start, before that window's convection |
| `cape_post_convection` | the mean over the host step's windows of CAPE at each window's end, after that window's convection |
| `cape_host_end` | CAPE on the state convection leaves at the end of the host step. Renamed from `dilute_cape_after_convection`, the same definition. |
| `cape_window_start` | each window's starting CAPE, one value per window. The per-run CAPE quantiles are taken from these. |
| `convective_active_time_fraction` | the time convection integrated with M_b > 0, over Δt |

**Nothing else acts within convection's call.** So:

- a window's end CAPE is the next window's start CAPE;
- a window after convection stops has CAPE_host_end at both ends.
- Without a substep, the step is one window.

**The analysis script adds the late-mean temperature profile** over the RH
profile's layers. That is output only.

## Gates, frozen

**The rule:** A6's rule and F4's classification across Δt_conv = 600, 300
and 150 s.

- A6's rule on the late means, |X_150 − X_300| ≤ 0.5·|X_300 − X_600|. It
  passes as "converged".
- Otherwise, "phase-sensitive" if every pairwise difference is within 2
  combined daily SEs.
- Otherwise, "changed equilibrium".
- Profiles are judged on their largest layer difference, and their SE test
  layer by layer.

**U1, conservation and tracer identities, as G1 and G6 of F4:**

- the conservation gate passes;
- residual ≤ 1e-12;
- the successor F6a, F6b and F6c pass;
- R_recycle ≤ 0.0042;
- no negative class.

**U2, identity:** the 300 s run reproduces T1's 1800 s run (`c_1800`,
`convective-cape-sampling.md`) to every printed digit, for every value the
two share.

**U3, the gated list.** No quantity may be a changed equilibrium:

- `cape_time_mean`;
- `cape_post_convection`;
- the RH profile and the temperature profile, 100-850 hPa;
- daily 𝓜 = ∫M_b dt;
- the detrainment profile 𝓓(p);
- the continuous plume top, as its daily-mean late mean;
- precipitation;
- FE, SE, PE and α_total,unique;
- daily subsidence;
- D_t, the late mean of `convective_active_time_fraction`.

**Reported, not gated:**

- `cape_host_start` and `cape_pre_convection`;
- CAPE's 10th, 50th and 90th percentiles from `cape_window_start`;
- the host-step duty cycle and the event counts;
- f_C;
- FE/SE/PE by the pressure of first condensation.

**U4, the tolerance between 300 and 150 s.** For every quantity in U3:

- |X_300 − X_150| ≤ 2·√(SE_300² + SE_150²), with daily SEs over the late
  period;
- every layer, for the profiles.

**Monotonic.** A quantity changes non-monotonically if (X_300 − X_600) and
(X_150 − X_300) have opposite signs and each exceeds 2 combined daily SEs.
For profiles, any layer.

## Outcomes, frozen

| outcome | what follows |
|---|---|
| **U1, U2, U3 and U4 pass, no non-monotonic quantity** | 300 s is selected as F0's convective substep for the new qualification run. |
| **U1, U2 and U3 pass, no non-monotonic quantity, U4 fails** | 300 s is not selected. 150 s is selected for the qualification run. Its own convergence is untested, and a 75 s run would test it; the user decides whether that runs first. |
| **A non-monotonic quantity, or U3 fails** | F0 stays not qualified. The next target is the closure's integration at the convective substep. |
| **U1 or U2 fails** | The runs are invalid; nothing is read from them. |

**Selection is for F0's configuration.** The library default of
`convective_substep` (none) is not changed by this protocol.

## What this cannot separate

- **Convergence of the closure at fixed host step.** It does not test
  convergence with the host step, which T1 and C1 addressed.
- **A dependence on Δt_conv below 150 s** that the three points do not
  resolve.
- **The frozen radiation.** Interactive radiation could couple differently.
  That is the interactive-radiation F0's question.

---

## Result (appended; nothing above this line changed)

**All three runs completed** 20 days at host 1800 s. They are invalid on the
validity-domain gate only, as F0's, T1's and C's were, under the frozen
radiation.

### U1: conservation and tracer identities, passes

| Δt_conv | residual | F6a worst / tolerance | F6b worst / tolerance | ledger | R_recycle |
|---:|---:|---:|---:|---:|---:|
| 600 s | 5.7e-16 | 0.025 | 0.00089 | 6.7e-16 | 2.1e-7 |
| 300 s | 7.5e-16 | 0.042 | 0.00085 | 5.7e-16 | 2.8e-7 |
| 150 s | 5.5e-16 | 0.030 | 0.00094 | 7.0e-16 | 3.9e-7 |

**F6c holds:**

- no ordering violations;
- no negative class;
- |PE − FE·SE| ≤ 1.1e-16;
- the inventories close to ≤ 7e-16 of throughput.

### U2: identity, passes

The 300 s run reproduces `c_1800` in all 709 shared printed values. The new
labelled outputs leave the state unchanged.

### U3: the gated list, fails

| quantity | 600 s | 300 s | 150 s | class | U4 (300 vs 150) |
|---|---:|---:|---:|---|---|
| `cape_time_mean`, J/kg | 47.3 | 51.8 | 57.1 | phase-sensitive | within |
| `cape_post_convection`, J/kg | 36.3 | 41.6 | 48.0 | phase-sensitive | within |
| RH profile (largest layer difference) | | 0.018 | 0.012 | phase-sensitive | within |
| **temperature profile** (largest layer difference), K | | 0.057 | 0.088 | **changed equilibrium** | within |
| daily 𝓜, kg m⁻² | 424.0 | 424.8 | 427.1 | phase-sensitive | within |
| 𝓓(p) | | 12.0 | 3.7 | converged | within |
| continuous p_top, hPa | 101.2 | 101.0 | 100.1 | phase-sensitive | within |
| precipitation, mm/day | 4.145 | 4.155 | 4.172 | phase-sensitive | within |
| FE | 0.99531 | 0.99536 | 0.99538 | converged | within |
| SE | 0.9664 | 0.9656 | 0.9654 | converged | within |
| PE | 0.9618 | 0.9611 | 0.9609 | converged | within |
| α_total,unique | 0.0382 | 0.0389 | 0.0391 | converged | within |
| daily subsidence | 1.696e7 | 1.699e7 | 1.708e7 | phase-sensitive | within |
| D_t | 0.125 | 0.143 | 0.161 | phase-sensitive | within |

**The temperature profile is the failing quantity.**

- **A6's rule fails on its largest layer difference,** at the upper
  troposphere: 0.088 K at 126 hPa against 0.5 × 0.057.
- **One layer's 600-vs-150 s difference lies outside 2 combined SEs.** At
  557.7 hPa, −0.078 K against 0.063: that layer's daily SE is 0.02 K,
  against 0.1-0.5 K elsewhere.
- **The magnitudes are small.** Every layer is within 0.09 K, against daily
  SEs of 0.02-1.2 K. They are recorded as found, and the criterion is not
  reopened.

**No quantity is non-monotonic.** U4 passes for every gated quantity.

**Most differences grow as the substep shrinks,** without leaving the noise:

- time-mean CAPE, 4.4 then 5.3 J/kg;
- post-convection CAPE, 5.3 then 6.4;
- 𝓜, D_t, subsidence, precipitation and the continuous top, likewise.

**The quantities that converge** are FE, SE, PE, α and 𝓓(p).

### Reported, not gated

| quantity | 600 s | 300 s | 150 s |
|---|---:|---:|---:|
| `cape_host_start`, J/kg | 86.0 | 102.0 | 154.4 |
| `cape_pre_convection`, J/kg | 52.9 | 51.9 | 57.1 |
| `cape_host_end`, J/kg | 36.2 | 40.2 | 46.0 |
| CAPE p10 / p50 / p90 (window starts), J/kg | 2.5 / 37.2 / 74.3 | 2.5 / 41.2 / 74.3 | 2.9 / 49.5 / 77.7 |
| host-step duty cycle | 0.379 | 0.621 | 0.825 |
| events (late) | 155 | 95 | 85 |
| steps per event | 1.17 | 3.14 | 4.66 |
| f_C mean | 0.58 | 0.74 | 0.97 |

### Outcome, as frozen

**U3 fails.** So:

- **F0 stays not qualified.**
- **No substep is selected.**
- **The next target is the closure's integration at the convective
  substep.**

### What else T2 shows: host-start CAPE follows the windows per host step

Across all eight runs with these outputs, `cape_host_start` and the
host-step duty cycle follow N = Δt_host/Δt_conv, not the host step.

| run | host step | Δt_conv | N | host-start CAPE, J/kg | duty cycle | time-mean CAPE, J/kg |
|---|---:|---:|---:|---:|---:|---:|
| C1 | 300 s | 300 s | 1 | 52.0 ± 4.9 | 0.148 | 51.9 |
| C1 | 600 s | 300 s | 2 | 62.6 ± 2.9 | 0.254 | 51.9 |
| T1/C0 | 900 s | 300 s | 3 | 73.8 ± 3.4 | 0.364 | 52.6 |
| **T2** | **1800 s** | **600 s** | **3** | **86.0 ± 5.7** | **0.379** | 47.3 |
| C1 | 1200 s | 300 s | 4 | 82.1 ± 3.4 | 0.458 | 51.1 |
| T1/C0 = T2 | 1800 s | 300 s | 6 | 102.0 ± 8.5 | 0.621 | 51.8 |
| T1/C0 | 3600 s | 300 s | 12 | 156.9 ± 24.4 | 0.817 | 48.6 |
| **T2** | **1800 s** | **150 s** | **12** | **154.4 ± 25.1** | **0.825** | 57.1 |

**At equal N and different host steps,** host-start CAPE and the duty cycle
agree:

- N = 12: 156.9 against 154.4 J/kg, and 0.817 against 0.825;
- N = 3: 73.8 against 86.0 J/kg, within 2 combined SEs of 13.3, and 0.364
  against 0.379.

**At an equal host step of 1800 s and different N,** they do not agree:
86.0, 102.0 and 154.4 J/kg.

**Time-mean CAPE follows the substep length, if anything:**

- 47.3 at 600 s;
- 48.6-52.6 at 300 s;
- 57.1 at 150 s.

**Consequence:** `convective-cape-sampling.md`'s interpretation (HS) was
confounded. Its erratum is appended there.

**A hypothesis for the next protocol, not measured here.** The trigger is
re-tested only at window starts, and convection stops for the rest of the
host step once CAPE is below it. The number of re-tests per host step then
shapes how events start and stop. That fits the event table: longer events
and higher f_C as N rises.
