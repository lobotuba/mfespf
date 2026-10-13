# Protocol: CAPE sampling and host-step destabilisation (follows F4/T1)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. Opened on the user's
decision (2026-10-13), option 2 after `convective-temporal-convergence.md`.

**It follows from** T1 in `convective-temporal-convergence.md`:

- **What converged.** A common 300 s convective substep made the RH profile,
  precipitation, FE, SE, α, 𝓜, 𝓓(p) and the continuous plume top
  host-step independent.
- **What did not.** CAPE_d (156.9, 102.0 and 73.8 J/kg at 3600, 1800 and
  900 s) and the duty cycle are a changed equilibrium, so G3 failed as frozen.

**Constraints:**

- **G3's verdict on T1 stands.** Nothing here re-grades it.
- F0 stays not qualified, interactive-radiation F0 stays blocked, and P1-P3
  stay blocked.
- No plume physics changes. No default changes. Radiation is never subcycled.

## Question

Is T1's step-dependent CAPE_d the sample of a destabilisation that arrives
once per host step? If so, it is a first-order splitting error with a limit,
and the column's state is the same. Or is it a different equilibrium of the
column?

## Read before this was frozen

**1. What convection receives.** F0's configuration sets
`condense_after_convection=True`, so:

- **Temperature:** convection acts on the committed temperature plus
  Δt·(frozen radiation + surface sensible heat), the dry adjustment, then
  the mass flux.
- **Moisture:** the water module applies surface evaporation, eddy diffusion
  and large-scale condensation after convection, at the host step. They reach
  convection at the next host step's start.
- **So all destabilisation arrives as one impulse per host step,** of size
  proportional to Δt. Convection's 300 s windows then remove it and idle for
  the rest of the host step.

**2. What CAPE_d is.** The published `dilute_cape` is the first window's CAPE,
on that post-impulse state: the peak of each host step.

**3. T1's three CAPE_d are a line in Δt.** A fit weighted by the daily
standard errors (24.4, 8.5 and 3.4 J/kg) gives CAPE_d = a + g·Δt:

- a = 46.0 ± 8.0 J/kg;
- g = 0.0310 ± 0.0069 J kg⁻¹ s⁻¹;
- χ² = 0.0015 on one degree of freedom.

**A6's rule and first order.** The rule's ratio of 0.5 is exactly what an
error proportional to Δt gives on a halving ladder. T1's ratio, 28.2/54.9 =
0.514, is that line plus noise. Recorded here as a reading, **not used to
re-grade G3.**

**4. The fit uses data already seen.** So the test below is out of sample:
new host steps, predicted before they run.

## Hypotheses

**HS, impulse sampling.** CAPE_d(Δt) = a + g·Δt.

- g is the destabilisation rate the host step accumulates before the first
  window samples it.
- a is the limit as the host step shrinks.
- The state is otherwise the same.

**HE, changed equilibrium.** The column's CAPE genuinely depends on the host
step. That would show as departure from the line at new host steps.

**What this experiment cannot separate HS from:**

- **Any other process that is linear in Δt over 300-3600 s** and is sampled at
  the same point.
- **Whether a = 46 J/kg is the continuous-time CAPE.** At a 300 s host step
  the impulse is still 300 s long, and the limit is extrapolated from there.

## Setup (T1's, unchanged)

- F0's configuration, with R4's state and the common frozen radiation rate.
- `convective_substep = 300`.
- 20 days, read over days 10-20.
- The same scripts as T1 (`t.py`, `t_gates.py`), with only the host step
  changed.

## C0: diagnostics on T1's runs (output only)

T1's three runs (3600, 1800 and 900 s) are rerun with one added output.

**`dilute_cape_time_mean`,** the host step's time-mean CAPE:

- the sum over convection's substeps of CAPE at the substep's start × its
  length;
- plus the remaining idle time × CAPE_after;
- all divided by Δt.
- It is a left-endpoint rule, published and never fed back.

**Reported per run:** late means of CAPE_d, CAPE_after (already published by
T0's code) and CAPE_mean, with daily standard errors.

## C1: a new host-step ladder, 300 s substep

**Host steps 1200, 600 and 300 s:** 4, 2 and 1 windows.

**Predictions, frozen from point 3,** with the fit's standard error:

| host step | CAPE_d predicted, J/kg |
|---:|---:|
| 1200 s | 83.2 ± 3.2 |
| 600 s | 64.6 ± 4.5 |
| 300 s | 55.3 ± 6.1 |

## Gates, frozen

**K1, conservation (as F4's G1):**

- the conservation gate passes;
- the worst relative residual is ≤ 1e-12;
- the successor F6a, F6b and F6c pass.

**K2, identity:** C0 reproduces T1's results to every printed digit.

**K3, the prediction.** At each of C1's host steps:

- |CAPE_d measured − predicted| ≤ 2·√(fit SE² + measured daily SE²).
- **HS is supported if all three pass** and refuted if any fails.

**K4, the state on C1's ladder:**

- **The rule:** A6's rule and F4's classification, applied to F4's list
  exactly as G3 applies them.
- **The pass condition:** no quantity is a changed equilibrium, apart from
  CAPE_d and the duty cycle.
- **Why CAPE_d is excluded here:**
  - under HS, CAPE_d's difference ratio on a halving ladder is 0.5 in
    expectation, exactly the rule's threshold, so its pass or fail there
    carries no information about HS;
  - K3 judges CAPE_d instead, with a sharper test.
- **The duty cycle** stays not a gate, as in F4.
- **Both are reported with their classification.**

**K5, agreement across the ladders** (reported, gated only through K4):

- C1's 1200 s values against T1's 900 s values, for F4's list;
- CAPE_after and CAPE_mean classified on both ladders.

**Secondary, reported, not gated:** per-event f_C, steps per event, the
continuous plume top, and FE/SE/PE by the pressure of first condensation.

## Decision, frozen

| outcome | what follows |
|---|---|
| **K1, K2, K3 and K4 pass** | HS supported. T1's G3 failure is recorded as host-step splitting of the CAPE sample, with limit a. F4's decision-tree row "converge the state but not the event timing" is offered to the user. Whether to run T2 and rerun F0's qualification is the user's decision. |
| **K3 fails** | HS refuted. F4's row 3 follows: find which of CAPE → M_b → B(z) → p_top stays step-dependent. |
| **K4 fails on a quantity other than CAPE_d or the duty cycle** | The state is not converged at short host steps. Recorded; no repair in this protocol. |
| **K1 or K2 fails** | The run is invalid; nothing is read from it. |

## Not in this protocol

**Re-triggering convection within a host step** is not tested here.

- **What it needs.** Convection would have to see the destabilisation
  arrive across its windows.
- **Why that is a separate change.** In this configuration the moisture part
  (surface evaporation, diffusion and large-scale condensation) is applied
  by the water module after convection. Presenting it within convection's
  windows means reordering or lagging another module's tendencies: a
  module-coupling change, and its own protocol.
- **What C1 offers instead.** Its 300 s host step is the reference where
  destabilisation and convection share one step, and so tests the same
  physics without the reorder.

**Literature on the class of error:**

- Williamson (2002), Mon. Wea. Rev. 130, 2024-2041: time-split against
  process-split coupling.
- Beljaars et al. (2004), ECMWF Seminar on numerical methods, 113-134: the
  numerics of physical parametrization.
- Wan et al. (2015), JAMES 7, 215-225, doi:10.1002/2014MS000368. Their
  process-coupling time-step convergence in CAM5 was of order 0.4, not 1.

---
