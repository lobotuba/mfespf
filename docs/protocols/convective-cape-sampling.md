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

## Result (appended; nothing above this line changed)

**All six runs completed** the planned 20 days. Their conservation residuals
are ≤ 7.5e-16. They are invalid on the validity-domain gate only, as F0's and
T1's were, under the frozen radiation.

### K2: identity, passes

C0's three runs reproduce T1's to every digit: all 701 shared printed values
are identical at 3600, 1800 and 900 s. The only additions are `cape_after`
and `cape_mean`.

### K1: conservation, passes

| host step | residual | F6a worst / tolerance | F6b worst / tolerance | ledger |
|---:|---:|---:|---:|---:|
| 3600 s | 5.4e-16 | 0.010 | 0.0013 | 2.1e-16 |
| 1800 s | 7.5e-16 | 0.042 | 0.00085 | 5.7e-16 |
| 900 s | 6.6e-16 | 0.0059 | 0.00058 | 3.1e-16 |
| 1200 s | 6.6e-16 | 0.033 | 0.00062 | 1.4e-15 |
| 600 s | 5.9e-16 | 0.011 | 0.00046 | 1.1e-15 |
| 300 s | 6.2e-16 | 0.0072 | 0.00018 | 5.5e-16 |

**F6c holds in all six:**

- no failing cells or steps;
- no ordering violations;
- |PE − FE·SE| ≤ 1.1e-16;
- the inventories close to ≤ 1.3e-15 of throughput;
- no negative class;
- R_recycle ≤ 4.9e-7.

### K3: the prediction, passes at every host step

| host step | predicted | measured | measured − predicted | tolerance, 2·√(SE²+SE²) | verdict |
|---:|---:|---:|---:|---:|---|
| 1200 s | 83.2 ± 3.2 | 82.1 ± 3.4 | −1.1 | 9.4 | pass |
| 600 s | 64.6 ± 4.5 | 62.6 ± 2.9 | −2.0 | 10.7 | pass |
| 300 s | 55.3 ± 6.1 | 52.0 ± 4.9 | −3.3 | 15.7 | pass |

**HS is supported.**

**What the test could discriminate.** If CAPE had stopped falling below
900 s, it would have stayed at T1's 73.8 J/kg. That lies 21.8 J/kg from the
300 s measurement against a 2-SE width of 11.9, so the test could have
failed and did not.

**A fit to all six points,** reported, not gated:

- a = 43.1 ± 4.0 J/kg;
- g = 0.0328 ± 0.0042 J kg⁻¹ s⁻¹;
- χ² = 0.22 on 4 degrees of freedom.

### K4: the state on C1's ladder (1200/600/300 s), passes

| quantity | 1200 s | 600 s | 300 s | class |
|---|---:|---:|---:|---|
| RH profile (largest layer difference) | | 0.021 | 0.017 | phase-sensitive |
| RH 300-700 | 0.7107 | 0.7114 | 0.7108 | phase-sensitive |
| precipitation, mm/day | 4.153 | 4.172 | 4.193 | phase-sensitive |
| FE | 0.99532 | 0.99531 | 0.99538 | phase-sensitive |
| SE | 0.9658 | 0.9653 | 0.9654 | converged |
| α_total,unique | 0.0387 | 0.0393 | 0.0391 | converged |
| cloud-base mass flux | 0.00493 | 0.00495 | 0.00497 | phase-sensitive |
| subsidence per day | 1.703e7 | 1.712e7 | 1.716e7 | phase-sensitive |
| R_recycle | 1.6e-7 | 1.7e-7 | 4.9e-7 | phase-sensitive |
| *CAPE_d (excluded, judged by K3)* | 82.1 | 62.6 | 52.0 | changed equilibrium |
| *duty cycle (not a gate)* | 0.458 | 0.254 | 0.148 | changed equilibrium |
| **CAPE_after** | 40.5 | 41.6 | 41.6 | **converged** |
| **CAPE_mean** | 51.1 | 51.9 | 51.9 | **converged** |

**G5's quantities, by the same rule:**

- continuous p_top 101.0/100.6/101.2 hPa, phase-sensitive;
- 𝓓(p) phase-sensitive;
- subsidence-mass profile converged;
- daily 𝓜 425.8/427.9/429.1 kg m⁻², phase-sensitive.

**No quantity other than CAPE_d and the duty cycle is a changed
equilibrium.** CAPE_d's difference ratio on this ladder is 10.55/19.5 =
0.541. As stated before the run, a ratio near 0.5 is what HS expects, and it
carries no information here.

### K5: across the ladders, reported

**C1's 1200 s against T1's 900 s:** every quantity in F4's list, CAPE_d,
CAPE_after, CAPE_mean and the duty cycle lie within 2 combined daily SEs. The
RH profile has none of its 23 layers outside.

**CAPE_after and CAPE_mean on C0's ladder** (3600/1800/900 s):

- **CAPE_after:** 37.4, 40.2 and 41.8 J/kg, phase-sensitive.
- **CAPE_mean:** 48.6, 51.8 and 52.6 J/kg, converged.
- **Across all six runs from 1800 s down,** CAPE_mean lies within 51.1-52.6
  J/kg and CAPE_after within 40.2-41.8 J/kg.

### Secondary, reported

**Events on C1:**

- f_C means 0.48, 0.45 and 0.46, medians 0.39-0.44;
- steps per event 1.7, 1.2 and 1.4.

**Origin-conditioned efficiencies:**

- FE by the pressure of first condensation varies by ≤ 8e-4 across all six
  runs;
- PE by ≤ 0.005, largest at 154 hPa;
- PE at 288 hPa is 0.9456-0.9476.

### What the result says, and what it does not

**1. T1's G3 failure is the CAPE sample, not the column.**

- T1's CAPE_d dependence is first-order in the host step. It predicted CAPE_d
  at three unseen host steps within their uncertainty.
- The host-step time-mean CAPE converges, and so do the state convection
  leaves and every state quantity in F4's list.
- What CAPE_d measures is the destabilisation that accumulates over a host
  step before convection's first window samples it. Its slope, about
  0.03 J kg⁻¹ s⁻¹ (≈ 110 J/kg per hour), is that rate.
- **G3's verdict on T1 stands as frozen.** This is the reading of it.

**2. The intercept a ≈ 43-46 J/kg is an extrapolation, not a state the model
reaches.**

- At a 300 s host step the host step is one convective window. There the
  sampled CAPE, 52.0, equals the time mean, 51.9.
- Below 300 s the line has not been tested, and the convective substep would
  have to shrink with it.

**3. The duty cycle** falls with the host step on both ladders, while 𝓜,
𝓓(p) and the state converge. That is F4's "numerical representation of
intermittency" case: events are counted per host step, and a shorter host
step splits one event's timing more finely.

**4. Re-triggering within a host step** was not tested and is not needed to
explain the result. It stays a separate, module-coupling question.

**Decision, as frozen:** K1, K2, K3 and K4 pass.

- **What is offered to the user:** F4's row "T1 and T2 converge the state but
  not the event timing".
- **The user's decision:**
  - whether to run T2 (host 1800 s; Δt_conv = 600, 300 and 150 s) and then
    rerun F0's qualification;
  - and whether a time-mean CAPE should be F4's CAPE quantity going forward.
    That would be a successor criterion applied prospectively, as with F6.
- **Not done:** no default changed, no substep length chosen, F0 not
  qualified.
