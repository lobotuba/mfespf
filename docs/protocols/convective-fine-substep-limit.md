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

## Result (appended; nothing above this line changed)

**Both runs completed** 20 days, with no substep refusal.

- Residuals are 7.7e-16 and 1.1e-15.
- F6a and F6b have no failing cells or steps.
- Both are invalid on the validity-domain gate only, as every
  frozen-radiation run is.

### L1-limit: no demonstrated limit, on one RH layer only

| quantity | 18.75 s | 9.375 s | 4.6875 s | class |
|---|---:|---:|---:|---|
| time-mean CAPE, J/kg | 127.2 | 189.6 | 214.8 | converged (62.3, then 25.3) |
| post-convection CAPE, J/kg | 126.0 | 189.3 | 214.8 | converged |
| D_t | 0.447 | 0.639 | 0.703 | converged |
| temperature profile | | 0.061 | 0.040 | phase-sensitive |
| **RH profile** | | 0.046 | 0.028 | **changed equilibrium** |
| daily 𝓜, kg m⁻² | 430.7 | 433.0 | 433.7 | converged |
| daily subsidence | 1.723e7 | 1.732e7 | 1.735e7 | converged |
| precipitation, mm/day | 4.160 | 4.166 | 4.165 | converged |
| 𝓓(p) | | 6.7 | 2.4 | converged |
| FE | 0.99532 | 0.99533 | 0.99534 | converged |
| SE | 0.9565 | 0.9514 | 0.9496 | converged |
| PE, α_total,unique | 0.9521, 0.0479 | 0.9470, 0.0530 | 0.9451, 0.0549 | converged |

**No quantity is non-monotonic over 37.5-4.6875 s.**

**The RH profile fails on one layer, 321.5 hPa:** 0.761, 0.716 and 0.743.

- The largest-difference ratio is 0.0275/0.0455 = 0.60, above 0.5.
- The 18.75-vs-9.375 s difference is outside 2 combined SEs: 0.0455 against
  0.037.
- The sign reverses, but only the first difference exceeds 2 SEs, so it is
  not non-monotonic by T2's definition.

**Verdict, as frozen: no demonstrated limit.** It rests on that one layer.
Every other quantity converges or is phase-sensitive, including CAPE, whose
differences now shrink.

### L1-212: agrees

The 4.6875 s explicit run agrees with the 150 s fractional run (`of_150`)
within 2 combined daily SEs, in every listed quantity and every layer:

| quantity | explicit 4.6875 s | fractional 150 s |
|---|---:|---:|
| time-mean CAPE | 214.8 | 213.1 J/kg |
| post-convection CAPE | 214.8 | 212.4 J/kg |
| 𝓜 | 433.7 | 424.4 kg m⁻² |
| precipitation | 4.165 | 4.153 mm/day |
| temperature and RH profiles | every layer within 2 combined SEs | |

**Richardson estimate:** not reported. The frozen condition, L1-limit, does
not hold.

### Context

| Δt_conv | active windows per host step | host steps with a stop (of 480) | mean undershoot | active windows ending below 70 |
|---:|---:|---:|---:|---:|
| 37.5 s | 15.2 of 48 | 414 | 1.5 J/kg | 0.057 |
| 18.75 s | 42.9 of 96 | 335 | 0.39 | 0.016 |
| 9.375 s | 122.7 of 192 | 237 | 0.067 | 0.0040 |
| 4.6875 s | 269.9 of 384 | 169 | 5.2e-5 | 0.0013 |

### What it says

**1. The explicit closure does turn over.**

- CAPE's differences shrink: 62.3, then 25.3 J/kg.
- The state it approaches is the one the fractional runs held at every
  substep, about 213-215 J/kg, and agrees with it everywhere listed.

**2. The fine-step behaviour is convection active for most of each host
step.** D_t is 0.70, and CAPE is held near 215 J/kg, about three times the
trigger.

**3. Not a limit by the frozen rule.** It fails on one RH layer, so it is not
a fine-step reference under the frozen criteria.

**4. Comparisons that read "converged" across F4, T1, T2 and C were made in
the overshoot regime.** Those runs used substeps of 150 s and more, where the
step routinely jumps the CAPE discontinuity
(`convective-cape-discontinuity.md`). That is now visible: the regime
changes below about 75 s.
