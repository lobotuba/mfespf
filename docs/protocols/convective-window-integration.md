# Protocol: closure-window integration, a crossed design

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The design is the
user's (2026-10-13), recorded here with its implementation fixed.

**It follows from** `convective-substep-convergence.md` (T2):

- U3 failed on the temperature profile.
- Host-start CAPE and the duty cycle followed N = Δt_host/Δt_conv.
- Time-mean CAPE drifted with the substep: 47.3, 51.8 and 57.1 J/kg.

**Constraints:**

- T2's gates are not loosened, and no substep is selected here.
- F0 stays not qualified, and interactive-radiation F0 and P1-P3 stay
  blocked.
- No module is reordered.
- No plume physics changes.
- Radiation is never subcycled.

## Question

Does the closure's result depend mainly on N, the number of convective windows
per host step, because a failed trigger ends convective evaluation until the
next host step? The design separates Δt_host, Δt_conv and N.

## Read before this was frozen

**1. What happens after a failed window.** In `mass_flux.convect`, a substep
whose closure returns M_b ≤ 0 ends the loop (`break`), and every remaining
window of the host step is skipped. M_b ≤ 0 means one of:

- no active plume;
- CAPE_d ≤ the trigger;
- convection would not reduce CAPE_d.

**2. Why re-evaluating would change nothing.** Nothing else changes T, the
vapour or the provenance classes within convection's call. So a skipped
window would see exactly the state the failed one saw, and the closure is a
deterministic function of that state.

**3. Consequence for I1.** Evaluating every remaining window (I1) is
therefore predicted to reproduce the current code bit for bit. That
prediction is frozen as I1's test below.

**4. Where the trigger is tested.** At every substep's start, not only at
window starts. With a substep set, a substep ends at the next window boundary
at the latest, and sooner only under the half-mass rule.

## I0: instrument without changing behaviour

**Per window,** published by convection as arrays over `convective_window`,
output only:

| output | meaning |
|---|---|
| `cape_window_start` | exists: CAPE at the window's start. The window number is the array index. |
| `cape_window_end` | CAPE at the window's end: the next window's start, or CAPE_host_end for the last |
| `convective_window_evaluated` | 1 if the closure was evaluated at the window's start, 0 if skipped after an earlier stop |
| `convective_window_status` | the code below |
| `convective_window_mass` | ∫M_b dt over the window, kg m⁻² |

**Status codes:**

| code | meaning |
|---:|---|
| 0 | active through the window |
| 1 | active, then stopped within the window at a later substep |
| 2 | evaluated, with CAPE_d at or below the trigger |
| 3 | evaluated, with no active plume |
| 4 | evaluated, with no CAPE removal (Γ_M ≤ 0) |
| 5 | not evaluated, skipped after an earlier stop |

- "Attempted" means evaluated: a plume is launched at every evaluated window.
- "Valid plume" is status ≠ 3.
- "Trigger passed" is status 0 or 1.

**Per host step, derived in the analysis script** from these arrays:

- the first and last active window;
- the first below-trigger window;
- the numbers of evaluated, skipped and active windows;
- the windows remaining after the first stop;
- the mass transport by window position.

**The decisive check, answered by I0:** after the first stop, are the
remaining windows skipped (status 5)? From point 1, this is predicted to be
yes in every run.

**Identity:** every I0 rerun of an existing configuration reproduces its
earlier result to every printed digit.

## The crossed matrix

All runs use T1's setup:

- F0's configuration;
- R4's state and the common frozen radiation rate;
- 20 days, read over days 10-20;
- I0 instrumentation on.

| N | Δt_host / Δt_conv | status |
|---:|---|---|
| 12 | 3600/300 | rerun of `c_3600` |
| 12 | 1800/150 | rerun of `u_150` |
| 12 | **900/75** | new; diagnostic only, never a proposed default |
| 6 | **3600/600** | new |
| 6 | 1800/300 | rerun of `c_1800` |
| 6 | **900/150** | new |
| 3 | 900/300 | rerun of `c_900` |
| 3 | 1800/600 | rerun of `u_600` |

**What the matrix contains:**

- **Same N:** {3600/300, 1800/150, 900/75}, {3600/600, 1800/300, 900/150} and
  {900/300, 1800/600}.
- **Same Δt_conv:**
  - 300 s: {900/300, 1800/300, 3600/300};
  - 150 s: {1800/150, 900/150};
  - 600 s: {1800/600, 3600/600}.
- **Same Δt_host:**
  - 1800 s: {1800/600, 1800/300, 1800/150};
  - 3600 s: {3600/300, 3600/600};
  - 900 s: {900/300, 900/150, 900/75}.

## Quantities and agreement, frozen

**The list Q:**

- host-start CAPE;
- D_t, the active-time fraction;
- event structure on the physical clock. An event is a maximal run of
  consecutive windows with mass > 0, across host-step boundaries. Two
  quantities: events per day, and mean event duration in seconds;
- time-mean CAPE;
- post-convection CAPE (`cape_post_convection`);
- daily 𝓜;
- daily subsidence;
- precipitation;
- the temperature profile and the RH profile, 100-850 hPa.

**Two runs agree on a quantity if** |X_a − X_b| ≤ 2·√(SE_a² + SE_b²), with
daily SEs over the late period. Profiles must agree in every layer.

**A factor F groups a quantity if** every pair of runs sharing a level of F
agrees on it.

**This convention has a known cost, disclosed here.** It is the one F4 and
T2 used, and it is kept so as not to loosen T2's gate. A 23-layer profile
judged layer by layer at 2 SE will disagree more often than any one layer
would by chance. That is most likely where a layer's daily SE is very small,
as at 557.7 hPa in T2. A profile disagreement is reported with the layers
that cause it.

## Classifications, frozen

| classification | condition |
|---|---|
| **HN, window count, supported** | N groups every quantity in Q. Host-start CAPE alone is insufficient. |
| **HC, substep length, supported** | Δt_conv groups every quantity in Q, and N fails to group at least one. |
| **HH, host-step coupling, supported** | Δt_host groups every quantity in Q, and N and Δt_conv each fail to group at least one. |
| **Mixed** | none of the above. No default is selected. |

**Reported with every classification:** for each quantity, which of the
three factors group it.

## I1: independent window evaluation

**Frozen now. Run only if I0 confirms that a stop skips the remaining
windows** (status 5 present after a stop).

**`convective_window_continuation`, an option on `convect`, convection and
the configuration:**

- **"stop":** today's code, and the default.
- **"evaluate":** a stop ends only that window. At every later window
  boundary convection:
  1. relaunches the plume on the current state;
  2. recomputes CAPE_d;
  3. re-evaluates the trigger;
  4. applies the operator if M_b > 0.
- **Nothing else changes.** No module is reordered.

**I1's test:** the eight-run matrix with "evaluate". A run counts as identical
if it matches its "stop" run to every printed digit, except the status codes:
under "evaluate", skipped windows become evaluated, with status 2, 3 or 4.

**What follows from I1:**

| result | what follows |
|---|---|
| **all eight identical,** as predicted | Window termination is not the source of the N-dependence. Where it lies goes to a new protocol, which the user decides. If the evidence points to destabilisation reaching convection only after its block ends, that is an operator-order protocol. |
| **any run not identical** | The reading in point 2 was wrong, and hidden state evolution within the call was being discarded. Then T2's ladder is run with "evaluate". The pass condition is that none of the following is a changed equilibrium by A6's rule and F4's classification: T and RH profiles; time-mean and post-convection CAPE; D_t; 𝓜 and subsidence; 𝓓(p); precipitation; FE, SE, PE and α_total,unique. U1's conservation and provenance identities must also pass. Host-start CAPE and events per host step are reported, not gated. |

## What this cannot separate

**Only two of the three factors are independent,** since N = Δt_host/Δt_conv.

- A dependence on two of them is classified as mixed. It is not attributed.

**Everything here uses the frozen radiation.**

---
