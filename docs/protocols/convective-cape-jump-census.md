# Protocol: what causes the remaining CAPE jumps (a diagnostic census)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The design is the
user's (2026-10-14), recorded here with its implementation fixed.

**It follows from** `convective-cape-definition.md`:

- **What the alternatives fixed.** D1 (positive area) and D2 (maximum work)
  are continuous through the internal neutral layers that made D0
  discontinuous.
- **What remains.** Large window-to-window CAPE changes remain under all
  three definitions: hundreds per run, up to about 430 J/kg.

**Constraints:**

- Diagnostic only. No physics, closure, default or trigger changes.
- The 70 J/kg trigger is a legacy control value.

## Question

Are the remaining large CAPE changes caused by:

- **plume termination** changes;
- **host-step forcing;**
- **launch or source changes;**
- **another discrete change** in the plume solution?

## Runs

**The census runs.** D1 (`positive_area`) drives the closure: explicit,
"stop" continuation, host 1800 s, 20 days, read over days 10-20. Two runs:

- **Δt_conv = 300 s,** the primary census;
- **75 s,** for how the jumps scale with the substep.

**Recording.** The run script is `o2.py` with recording wrappers on
`entraining_plume`, `cloud_base_mass_flux`, `apply_mass_flux` and `convect`.
They record and change nothing.

**Identity.** The census runs must reproduce `s2_positive_area_300` and
`s2_positive_area_75` to every printed digit.

**D2 is computed on the same plumes as a secondary diagnostic.** D0, D1 and
D2 are fields of every plume. No D2-driven census is run.

## The launch sequence and its transitions

**Launches, in order within each host step:**

- **main launches:** the plume each substep closes on, the one passed to
  `cloud_base_mass_flux`;
- **the end launch:** the plume on the state convection leaves, which is the
  last plume `convect` launches.

**Excluded:**

- the closure's displaced plume, which is hypothetical;
- any fractional-crossing trial. There are none, under "explicit".

**A transition is one of two kinds:**

- **convective:** between consecutive main launches within a host step, or
  from the last main launch to the end launch. Only convection has changed
  the state.
- **host:** from host step n's end launch to host step n+1's first main
  launch. Only the non-convective processes have changed the state:
  - step n's water module (large-scale condensation and its latent heat,
    surface evaporation, eddy diffusion), and any other step-n tendency
    committed after convection;
  - step n+1's radiation, sensible heat and dry adjustment.

**A large jump** is |J| > 25 J/kg in D1, J = C_after − C_before across a
transition. That is the census's existing criterion. Positive and negative
jumps are reported separately throughout.

## Recorded per launch

- **CAPE:**
  - D0, D1 and D2;
  - the integrated positive and negative work over the plume's extent;
  - the maximum-work layer;
  - the LFC, the source, the first internal non-buoyant layer, and the
    number and depth of internal non-buoyant layers.
- **The termination:**
  - its layer;
  - its continuous pressure (the zero crossing of W);
  - the plume's depth.
- **The source layer:** T, q and moist static energy.
- **The closure's M_b** (main launches only).
- **Totals:** condensate formation and conversion.
- **Time:** elapsed since the host step's start, from the substep lengths
  `convect` applied; the launch's index; and its kind (main or end).

**For large-jump transitions only, both sides' profiles:**

- buoyancy;
- parcel temperature;
- vapour and liquid;
- the environment's T and q.

**Per host transition, the environment change by process.** It comes from the
run's labelled tendencies: step n's contributions after convection, and step
n+1's radiation, sensible heat and dry adjustment.

- Reported as each process's largest |ΔT| (K) and |Δq| (kg m⁻²) over the
  column.
- Also its residual against the measured environment difference.

## Labels, frozen (several may apply)

| label | condition |
|---|---|
| **J1, termination** | the termination layer changes, or the continuous termination pressure moves by ≥ 10 hPa. Index-only changes (< 10 hPa) are reported as their own subset. |
| **J2, host forcing** | the transition is a host transition |
| **J3, launch** | the source layer or the LFC changes. Which one is recorded. |
| **J4, profile** | not J1 and not J3, and some layer's |Δ(b·dz)| ≥ 5 J/kg |
| **J5, unexplained** | none of J1, J3 or J4. The change is spread over many layers, each under 5 J/kg. |

## Counterfactuals (diagnostic re-evaluations only)

**A plume path P is three things:**

- its parcel profile: parcel temperature, vapour and carried liquid by layer;
- its LFC;
- its termination.

**An environment E is** the layer temperature and vapour, with dz from E's
temperature.

**C(E, P)** is the definition's CAPE from the buoyancy g·[(T_u − T)/T +
0.61(q_u − q) − l] of P's parcel against E. It is integrated over P's
extent, the LFC to the termination inclusive:

- D1: the positive part;
- D2: the maximum running work from the LFC.

**Evaluated for every large jump:**

| counterfactual | what it is |
|---|---|
| C(E_n, P_n) and C(E_n+1, P_n+1) | the actual values. The recomputation must match them; the residual is reported. |
| C(E_n+1, P_n) | the old path, its parcel and extent, in the new environment |
| C(E_n, P_n+1) | the new path in the old environment |
| C_ext | the new plume's own buoyancy integrated over the old extent, holding only the termination and LFC fixed |

**The decomposition:**

- **the environment term:** C(E_n+1, P_n) − C(E_n, P_n);
- **the plume term:** C(E_n, P_n+1) − C(E_n, P_n);
- **the interaction:** the remainder of J.

**What cannot be held fixed, and is documented, not replaced:**

- **The parcel profile does not respond to the new environment** in the
  counterfactual. The real one would, through entrainment. So C(E_n+1, P_n)
  is "the old parcel path", not a re-solved plume.
- **The source layer and the entrainment rate** are part of P. The rate is a
  constant here.

## Verdicts, frozen (on the 300 s census; D1 jumps)

| verdict | conditions |
|---|---|
| **Termination supported** | (i) ≥ 50 percent of large jumps are J1; (ii) holding the old extent reduces the jump: median \|C_ext − C(E_n, P_n)\| ≤ 0.5 × median \|J\| over large jumps; (iii) J1 is less than half as frequent among non-jump transitions as among large jumps; (iv) (i) and (ii) hold for positive and for negative jumps separately. Otherwise reported per sign. |
| **Host forcing supported** | (i) ≥ 50 percent of large jumps are J2, and the large-jump rate per host transition is more than twice that per convective transition; (ii) the environment term carries the jump: median over large jumps of env/J ≥ 0.5, same sign. |
| **Mixed** | both supported. Or one supported, with the environment term and the plume term each carrying a median of at least 0.25 of J in the same sign. |
| **Neither identified** | none of the above. The census stops, and D1 is not extended. |

**Reported, not judged:**

- the counts per label and per sign;
- the same census on the 75 s run;
- D2's jumps, labelled the same way.

## The conditional D1 extension, frozen

**Δt_conv = 18.75 and 9.375 s,** D1-driven, run only if:

- (a) the verdict is host forcing; or
- (b) the verdict is termination or mixed, and the large-jump magnitude per
  day of physical time (Σ|J| per day over the late period) is no larger at
  75 s than at 300 s.

**Otherwise it is not run.**

**If run, the extended ladder is 150, 75, 37.5, 18.75 and 9.375 s:**

- A6's rule and F4's classification on 37.5/18.75/9.375 s, and on
  150/75/37.5 s, which exists already;
- non-monotonic over all five;
- **quantities:** time-mean and post-convection D1 CAPE, D_t, the
  temperature and RH profiles, RH 300-700, the lapse rate, 𝓜, subsidence,
  𝓓, precipitation, and FE, SE, PE and α_total,unique;
- **also reported:** the large-jump rate and magnitude per day, by
  substep.

**A finite-step limit for D1 is demonstrated** if no quantity is a changed
equilibrium on 37.5/18.75/9.375 s and none is non-monotonic. That would
satisfy only the frozen rule for D1's limit. **Selecting a definition, and
any trigger protocol, remain the user's decisions.**

---
