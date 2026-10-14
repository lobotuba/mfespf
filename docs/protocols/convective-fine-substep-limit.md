# Protocol: the explicit closure's fine-substep limit (L1)

**Frozen before any run of this protocol.** Results are appended below a line
at the end; nothing above it changes afterwards. The user authorised it
(2026-10-13): "1 then 2". 2 is `convective-cape-discontinuity.md`, run after
this one.

**It follows from** `convective-trigger-overshoot.md`. O1 found no
demonstrated limit: time-mean CAPE was 75.1, 93.5 and 127.2 J/kg at 75, 37.5
and 18.75 s, with the differences growing. The fractional runs held about
212 J/kg independent of the substep, but failed O-G3.

**Constraints:**

- No code changes.
- Every default unchanged.
- F0 stays not qualified.
- Nothing is adopted.

## Question

Does the explicit closure approach a limit as the substep shrinks further?
And if it does, is that limit near the fractional runs' ~212 J/kg?

## Runs

O1's setup, unchanged:

- F0's configuration, with R4's state and the common frozen radiation rate;
- host step 1800 s;
- "explicit" crossing, "stop" continuation;
- 20 days, read over days 10-20;
- the same script (`o.py`).

**Δt_conv = 9.375 and 4.6875 s** (N = 192 and 384), added to O1's 75, 37.5
and 18.75 s.

**The substep cap.** The integrator refuses a host step needing more than 500
substeps. At 4.6875 s a fully active host step takes 384 window-bounded
substeps, so the cap can bind only if the half-mass rule adds more than 116.
**A refusal is a result, not a reason to raise the cap.**

## The quantities

**O1's primary list:**

- time-mean CAPE;
- post-convection CAPE;
- D_t;
- the temperature and RH profiles;
- 𝓜;
- subsidence;
- precipitation;
- 𝓓(p);
- FE, SE, PE and α_total,unique.

**Also reported:**

- the undershoot statistics;
- active windows per host step;
- host steps with a stop.

## Verdicts, frozen

**L1-limit, a demonstrated limit:**

- A6's rule and F4's classification on 18.75, 9.375 and 4.6875 s give no
  changed equilibrium in any primary quantity, profiles judged layer by
  layer;
- and no quantity is non-monotonic over 37.5-4.6875 s, by T2's definition.

**Otherwise: no demonstrated limit.**

**L1-212, agreement with the fractional runs.** Reported whichever L1-limit
says.

- For time-mean CAPE, post-convection CAPE, the temperature and RH profiles,
  𝓜 and precipitation: does the 4.6875 s explicit run agree with the 150 s
  fractional run (`of_150`) within 2 combined daily SEs, every layer for the
  profiles?
- "Agrees" only if every listed quantity does.

**Reported, not judged:** if L1-limit holds, the first-order Richardson
estimate 2·X_4.6875 − X_9.375 for time-mean CAPE.

## What this cannot separate

- **Any limit found is the explicit closure's, with its CAPE definition.**
  `convective-cape-discontinuity.md` asks whether that definition is
  discontinuous, and why.
- **Agreement with the fractional runs** would not make them right. They
  failed O-G3.

---
