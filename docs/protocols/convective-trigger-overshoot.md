# Protocol: trigger overshoot in the explicit closure (O1, O2)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The design is the
user's (2026-10-13), recorded here with its implementation fixed.

**It follows from** `convective-window-integration.md`:

- The window stop rule is eliminated: I1 was bit-identical.
- Time-mean CAPE, post-convection CAPE, the temperature profile and event
  duration group by the substep length.
- Time-mean CAPE rose from 47-50 J/kg at 600 s to 79 J/kg at 75 s, where
  post-convection CAPE was 70.5 J/kg, against a trigger of 70.

**Constraints:**

- F0 stays not qualified, and interactive-radiation F0 and P1-P3 stay
  blocked.
- Every default is unchanged.
- T2's verdict is not reopened, and the layer-by-layer criterion, 557.7 hPa
  included, is not relaxed.
- The trigger stays at its frozen 70 J/kg. No parameter is fitted to RH,
  precipitation, temperature or the fine-step outcome (O-G8).
- If this protocol supports the mechanism, it still needs a new F0 successor
  qualification. No earlier run is qualified retroactively.

## Question

Does the explicit step, with M_b closed at the substep's start and applied
over the whole substep, carry CAPE past the trigger? And is that overshoot
the source of the substep dependence?

## Setup (T2's, unchanged)

- F0's configuration, with R4's state and the common frozen radiation rate;
- host step 1800 s;
- 20 days, read over days 10-20;
- the window stop rule ("stop");
- I0's instrumentation on.

## O1: is there a fine-step limit?

**Runs:** Δt_conv = 75, 37.5 and 18.75 s (N = 24, 48 and 96), with the
explicit closure. The 150 s case is T2's `u_150` and I0's 1800/150. The 75 s
run at 1800 s is new.

**Per window, from the existing outputs:**

- C_before (`cape_window_start`) and C_after (`cape_window_end`);
- the trigger undershoot U_C = max(0, 70 − C_after).

**Undershoot, reported per run over the late period:**

- the mean undershoot over active windows (status 0, 1 or 6);
- the maximum undershoot;
- the event-integrated undershoot: U_C at the end of each event's last
  active window, its mean over events;
- the fraction of active windows ending below the trigger.

**The primary quantities:**

- time-mean CAPE;
- post-convection CAPE;
- D_t;
- the temperature and RH profiles, 100-850 hPa;
- daily 𝓜;
- daily subsidence;
- precipitation;
- FE, SE, PE and α_total,unique;
- 𝓓(p).

**Each is judged by A6's rule and F4's classification** on 75, 37.5 and
18.75 s, with profiles judged layer by layer. **Non-monotonic** is as in
T2, applied over the four-step ladder 150, 75, 37.5 and 18.75 s.

**O1 verdict, frozen:**

| verdict | condition |
|---|---|
| **Finite limit supported** | No primary quantity is a changed equilibrium on 75/37.5/18.75, and none is non-monotonic over 150-18.75. The 18.75 s run is then the fine-step reference. |
| **No demonstrated limit** | Anything else. Includes time-mean or post-convection CAPE still drifting without reduction, a profile failing the rule, a new trend in transport, or a non-monotonic quantity. |

**Not gated:** where post-convection CAPE sits relative to the trigger. It is
reported as one of:

- **at the trigger:** |C_post − 70| ≤ 2 daily SEs at 18.75 s;
- **above it;**
- **below it.**

The limit is not required to equal 70 J/kg.

**Conservation and tracer identities (U1's)** must pass in every O1 run.

## O2: trigger-limited fractional integration

**The option: `convective_trigger_crossing`** on `convect`, convection and
the configuration:

- **"explicit":** today's code, and the default;
- **"fractional":** the variant below.

**Its meaning, chosen and used consistently:** the plume and its mass flux
are diagnosed at the substep's start, and convection acts for the fraction of
the substep needed to reach the trigger. The effective substep is
f·h, where h is the substep the explicit integrator chose.

**Each active substep under "fractional":**

1. **Keep the starting state** (T₀, w₀, the classes), and C₀ = the plume's
   CAPE on it. C₀ > 70, or the closure gives M_b = 0 and nothing happens, as
   now.
2. **Take the full trial step** from the starting state: `apply_mass_flux`
   over h. Then relaunch the plume on the trial state to get C_full.
3. **If C_full ≥ 70,** accept the full step: f = 1, and the update is exactly
   the explicit one.
4. **Otherwise solve C(f) = 70 for f in (0, 1).**
   - C(f) is the plume's CAPE on the state produced by `apply_mass_flux`
     over f·h from the starting state. It is the same diagnostic the trigger
     uses.
   - **Every trial starts from the starting state.** The whole coupled
     operator is recomputed:
     - plume transport, subsidence, entrainment and detrainment;
     - total water and enthalpy;
     - condensate formation, the carried cloud water, its conversion, and the
       falling precipitation's evaporation.
   - **Nothing is interpolated from the full step,** and no tendency is
     scaled on its own.
   - **The root finder is bracketed:**
     - the bracket is [0, 1], with g(f) = C(f) − 70, g(0) > 0 and g(1) < 0;
     - the first trial is the linear estimate f = (C₀ − 70)/(C₀ − C_full);
     - then Illinois false position, always keeping a sign-changing bracket;
     - it stops when |C(f) − 70| ≤ ε_C = 1e-3 J/kg, or when the bracket is
       narrower than 1e-12.
5. **Accept the fractional step.**
   - It is the operator over f·h. The provenance classes take the same step,
     `provenance_substep` at scale M_b·f·h.
   - Convection then ends for the rest of the host step, as for any stop: the
     trigger has been reached. The window is given **status 6** ("limited
     to the trigger").
   - The time bookkeeping advances by f·h, and the time-mean CAPE counts C₀
     over f·h.
6. **If the bracket narrows below 1e-12 without meeting ε_C,** C(f) is
   discontinuous at the crossing. The lower end of the bracket is accepted,
   the one with CAPE above the trigger, so there is no overshoot. The event
   is counted as an **unresolved crossing**, with its residual.

**Published per window:** `convective_window_fraction`. It is f for a
limited window, 1 for an unlimited active window, and 0 where the window
did not convect.

**Published per step:** `convective_trigger_unresolved_crossings`.

**Off ("explicit"), the code is today's,** bit for bit, apart from the new
outputs.

**O2 runs:** "fractional" at Δt_conv = 600, 300, 150 and 75 s, host 1800 s.
The 75 s run is included now, before any outcome is seen.

**The explicit comparison** is T2's 600/300/150 s runs and O1's 75 s run, all
at host 1800 s.

## Gates, frozen

**O-G1, identity.** "Explicit" reruns of 1800/150 s and 1800/600 s reproduce
I0's runs to every printed digit.

**O-G2, conservation:**

- **every O1 and O2 run** passes U1: the conservation gate, residual ≤ 1e-12,
  the successor F6a, F6b and F6c, R_recycle ≤ 0.0042, and no negative class;
- **trial updates:** a unit test checks that `apply_mass_flux` over trial
  fractions conserves energy and water. Only accepted updates enter a run.

**O-G3, trigger targeting.**

- Every limited window has |C_after − 70| ≤ ε_C = 1e-3 J/kg, with no
  unresolved crossing, in every O2 run.
- Any unresolved crossing fails O-G3. They are counted and reported.

**O-G4, a fine-step reference.** O1 supports a finite limit, or the protocol
records that it does not.

**O-G5, O-G6 and O-G7, convergence of O2's own ladder.** The O2 ladder
(600/300/150 s, A6's rule and F4's classification, profiles layer by layer)
has no changed equilibrium, and no quantity is non-monotonic over
600-75 s, in:

- **O-G5:** the temperature and RH profiles;
- **O-G6:** 𝓜, subsidence, precipitation and 𝓓(p);
- **O-G7:** FE, SE, PE and α_total,unique;
- **also:** time-mean and post-convection CAPE, and D_t.

**O-G8, no hidden tuning.** The only target is the frozen trigger, 70 J/kg.

**Toward the reference** (needs O-G4). For every quantity in O-G5 to O-G7,
time-mean CAPE and D_t:

- **R1, closer than explicit.** At each of 600, 300 and 150 s:
  - |X_frac − X_ref| ≤ |X_expl − X_ref|;
  - for a profile, compared on its largest layer difference from the
    reference.
- **R2, agreement at 150 s.** The 150 s fractional run agrees with the
  reference within 2 combined daily SEs, in every layer for the profiles.

**Reduced dependence, R3:** for time-mean CAPE and the temperature profile's
largest layer difference, |X_frac(600) − X_frac(150)| < |X_expl(600) −
X_expl(150)|.

## Verdicts, frozen

| verdict | condition |
|---|---|
| **Overshoot supported as the mechanism** | O-G1 to O-G8 pass, and R1, R2 and R3 hold for every listed quantity. |
| **Rejected** | Time-mean CAPE meets R1, but the temperature or RH profile, or one of 𝓜, subsidence and precipitation, fails R1. Then the step dependence comes from holding the plume and its mass flux fixed through the window, not only from crossing the trigger. |
| **Partial** | Anything else with O-G1 to O-G3 passing. Recorded as found, with no adoption. |
| **Diagnostic only** | O1 finds no demonstrated limit. O2 is run and reported, R1 to R3 cannot be evaluated, and nothing is adopted. |

**If supported, the recorded statement:**

> The unlimited explicit closure has a timestep-dependent trigger overshoot.
> Fractional-window integration removes that overshoot and reproduces the
> demonstrated fine-step limit at coarser timesteps.

**Then** the fractional variant is an integration correction, not a
calibration. Its adoption, and the F0 successor qualification, are the user's
decision.

## What this cannot separate

- **A limit is not a solution.** A finite limit of the explicit closure is the
  limit of this closure, with its trigger and τ_c. It is not evidence that
  the closure is right.
- **Fractional integration still holds the plume and M_b fixed through
  f·h.** If coarse steps still differ after limiting, that is the remaining
  candidate. This protocol does not test it beyond the rejection row.
- **Everything uses the frozen radiation.**

---
