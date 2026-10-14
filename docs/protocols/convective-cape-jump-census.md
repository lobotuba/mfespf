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

## Result (appended; nothing above this line changed)

**Both census runs reproduce their S2 runs to every printed digit:**
`s2_positive_area_300` (890 shared values) and `s2_positive_area_75` (956).
The wrappers changed nothing. The recomputed C(E_n, P_n) and C(E_n+1, P_n+1)
match the plumes' own D1 exactly: residual 0.

### The 300 s census: two-state termination switching

**1188 transitions over days 10-20,** 480 host and 708 convective. **304 are
large (|ΔD1| > 25 J/kg),** and they split exactly in half:

| | count | labels | median \|J\| | where |
|---|---:|---|---:|---|
| **positive** | 152 | J1 + J2, every one | 236 J/kg | every one at a host transition, step n's end launch to step n+1's first launch |
| **negative** | 152 | J1, every one | 241 J/kg | every one within convection, typically at the second window (300 s into the host step) |

**No J3, J4 or J5.** The source and LFC never changed. No J1 is index-only:
every termination change is ≥ 10 hPa. In the examples the termination rises
329-335 hPa at the host transition and falls 358-361 hPa one window later.

**The host forcing behind a positive jump is small.** The environment changes
by at most 0.08-0.084 K in any layer. The largest labelled contributions are:

- convection's mass-flux transport and dry adjustment;
- radiation, about 0.15 K per step;
- the water module's latent heat.

**D2, on the same plumes:** 304 large jumps, 152 each sign, labelled
identically, median 208 and 215 J/kg.

### Counterfactuals (300 s)

| | all large | positive (host) | negative (convective) |
|---|---:|---:|---:|
| median \|J\| | 238 | 236 | 241 |
| environment term / J, median | 0.015 | 0.003 | 0.136 |
| plume term / J, median | 0.958 | 0.941 | 0.986 |
| fixed-extent change (C_ext), median | 85.5 | **2.6** | 240.9 (truncated) |
| C_ext truncated | 152 | 0 | **152** |

- **Positive jumps:** holding the old extent removes 99 percent of the jump.
- **Negative jumps:** the frozen C_ext cannot be evaluated. In all 152 the
  new plume terminates below the old extent, so its parcel does not exist
  over it. This is the limitation documented above, and the median 240.9
  is that truncation, not a measurement.

**A post-hoc mirror, not frozen, labelled as such.** The old plume's own
buoyancy over the new, lower extent: defined in all 152 negative jumps (none
truncated). It reproduces the new CAPE to within a median 7.7 J/kg of a
median 241 J/kg jump.

### The 75 s census, reported

**3187 transitions, 798 large:** 399 positive (all J1 + J2) and 399 negative
(all J1). Medians are 256 and 256 J/kg.

- **Rates:** large jumps occur at 0.83 per host transition, against 0.15 per
  convective one.
- **J1 among non-jump transitions:** 0.058.
- **Counterfactuals,** as at 300 s:
  - environment term / J 0.015;
  - plume term / J 0.988;
  - positive-jump C_ext change 2.2 J/kg;
  - negative C_ext truncated in all 399.

### Verdicts, as frozen

**Termination: (i) and (iii) hold.**

- (i) J1 is 100 percent of large jumps.
- (iii) J1 occurs in 24 percent of non-jump transitions, against 100
  percent.

**(ii) holds overall:** the fixed-extent change has a median of 85.5 J/kg,
against the 119 J/kg limit.

**(iv), per sign, does not:**

- (ii) holds for positive jumps: 2.6 against 118.
- It cannot be shown for negative ones by the frozen counterfactual: 152 of
  152 truncated.

**So the verdict is reported per sign:**

- **Positive jumps: termination supported,** at host transitions. The host
  forcing is the trigger, and the CAPE change is carried by the termination
  (plume term 0.94, environment term 0.003).
- **Negative jumps:** termination supported by (i) and (iii) and by the
  post-hoc mirror. Not by the frozen (ii), which cannot be evaluated.

**Host forcing: not supported.**

- The positive jumps all occur at host transitions, and J2 is 50 percent of
  large jumps.
- But the large-jump rate per host transition (0.32) is not more than twice
  the convective one (0.21).
- The environment term carries a median of 1.5 percent of J.

**Mixed: no.** The environment term is not material.

### The D1 extension, as frozen: not run

Condition (b) needs the large-jump magnitude per day at 75 s to be no larger
than at 300 s. It is 20806 against 7704 J/kg per day, 2.7 times larger.
Condition (a), a host-forcing verdict, does not hold. **The 18.75 and 9.375 s
D1 runs are not run.**

### What the result says, and what it does not

**1. The remaining jumps are the plume's own termination switching between
two states.**

- **The cycle:** a deep state and a shallow one about 330-360 hPa apart. A
  0.08 K host-step change flips the plume deep; one convective window flips
  it back.
- **The size:** each flip moves D1 by about 240 J/kg.
- **What D1 and D2 are doing:** they follow the plume, as designed. The
  discontinuity is now the plume's own.

**2. A reading, not measured here.** For the termination to switch on a
0.08 K change, the plume's running work must have an interior minimum within
a hair of zero. The plume either carries through it, or stops there, by the
sign of W at that layer. The scheme's binary penetrate-or-stop at a marginal
work minimum is then the discontinuity left in the closure, and finer
substeps resolve the switching more often rather than less (rates 0.32 →
0.83 per host transition at 300 → 75 s).

**3. Not concluded:**

- whether the termination rule is the right physics;
- whether a plume that partly penetrates a marginal layer (a distribution of
  tops, or a continuous top) is the right representation.

Those are physics choices for the user. The continuous top
(`plume_top_pressure_continuous`) already exists as a diagnostic. It is not
used for deposition.

## User disposition (2026-10-14): the branch is closed, appended

**Option 3:** stop this branch, and record binary plume termination as a
known limitation. **Recorded as the user worded it:**

> The CAPE-definition branch is closed without adoption. D1 and D2 eliminate
> the first-neutral-layer discontinuity in D0 and therefore confirm that the
> historical CAPE definition is inconsistent with the modeled plume's
> dynamical extent. Neither alternative qualifies when used to drive the
> closure. The remaining CAPE discontinuities arise from binary switching of
> the single bulk plume between shallow and deep termination states. The
> switching is intrinsic to the plume trajectory: plume-geometry changes
> account for approximately 96% of each diagnosed jump, while host-step
> environmental changes act only as the perturbation selecting the
> trajectory. Finer temporal resolution increases the resolved switching
> frequency rather than removing the discontinuity.
>
> No interpolation correction is adopted. Distributing tendencies between the
> shallow and deep tops would require an unmeasured population weight and
> would constitute a new closure.
>
> No plume spectrum is implemented in this branch. A plume-population
> formulation is retained as a possible successor architecture but requires
> independent constraints on population statistics, mass-flux partition,
> entrainment, termination probabilities and convective memory.
>
> The current bulk-plume scheme remains unqualified. No CAPE definition,
> convective substep or trigger threshold is selected. F0,
> interactive-radiation F0 and P1-P3 remain blocked.

**The three definitions:**

| definition | disposition |
|---|---|
| D0 `first_lnb` | retained as the bit-identical default only. It has a measured first-neutral-layer discontinuity. |
| D1 `positive_area` | retained as an experimental diagnostic. It removes D0's discontinuity but inherits the shallow/deep switching, and did not pass the driven-closure gates. |
| D2 `maximum_work` | retained as an experimental diagnostic, consistent with the plume's running work. It did not converge in driven CAPE or D_t. |
| 70 J/kg trigger | a legacy control for the historical configuration only, not calibrated or validated for D1 or D2 |

**Kept as research diagnostics, all with their defaults unchanged:**

- the three definitions;
- the census instrumentation;
- the window record;
- the fractional crossing;
- the window continuation.

**The single-column route is not called disproven.** The user lists what this
branch established:

- compensating subsidence can dry the free troposphere;
- condensate-return placement strongly controls humidity;
- fall-path evaporation avoids the forced evaporation-recondensation loop;
- molecule-once precipitation efficiencies can be measured exactly;
- D0's CAPE geometry is defective;
- a single deterministic plume can switch between shallow and deep states
  under tiny perturbations.

**But this mass-flux closure cannot supply a unique, timestep-independent
convective response.** So it is no longer treated as the path to P1-P3.
