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
