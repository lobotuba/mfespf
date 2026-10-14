# Protocol: the closure's CAPE definition (D0, D1, D2)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The design is the
user's (2026-10-14), recorded here with its implementation fixed.

**It follows from** `convective-cape-discontinuity.md`:

- **Every one of the 301 CAPE jumps is LNB truncation.** A mid-tropospheric
  layer at neutral buoyancy ends CAPE_d. A median 328 J/kg and 9.3 km of
  buoyant ascent drops out, while the plume keeps rising to the same top.
- **The scheme holds two inconsistent notions of convective reach:**
  - CAPE stops at the first layer with b ≤ 0;
  - the plume stops where its accumulated work returns to ≤ 0.

**Constraints:**

- F0 stays not qualified, and interactive-radiation F0 and P1-P3 stay
  blocked.
- Every default is unchanged.
- The plume, its termination and every other closure component stay
  unchanged. Only the CAPE quantity the closure reads is varied.
- **No parameter is tuned.** In particular, no trigger is fitted to RH, the
  lapse rate, precipitation or a duty cycle.

## The three definitions, all from the same plume solution

**The plume's own cumulative work,** read before freezing (`entraining_plume`):

- It starts at the LFC, the first layer above the source with b > 0.
- It adds one layer at a time: W_k = Σ_{j=LFC}^{k} b_j·dz_j, the layer's
  buoyancy times its hydrostatic thickness.
- The plume terminates at the first layer where W_k ≤ 0 ("top").

**Your suggested trapezoid, (B_k + B_k+1)/2·Δz, is not used.** The
quadrature is the plume's own layer sum, so that D2 and the termination use
the same W exactly. Negative buoyancy below the LFC is excluded, as the
plume excludes it.

| | definition | purpose |
|---|---|---|
| **D0, `first_lnb`** | Σ b·dz from the LFC to the first layer with b ≤ 0, excluding it: today's CAPE_d | reproduces existing behaviour |
| **D1, `positive_area`** | Σ max(b, 0)·dz from the LFC to the top, inclusive | is removing the first-LNB cutoff enough? |
| **D2, `maximum_work`** | max over the LFC to the top of W_k | the definition consistent with the plume's own propagation. **The primary candidate.** |

**D2 is not W(top).** The net work at termination is ≤ 0 by construction,
so it cannot be a closure amplitude.

**The literature behind D2's reading.** Entraining CAPE is read as a
predictor of maximum updraft speed within an entraining-plume energy
framework:

- Peters, Chavas, Su, Morrison and Coffer (2023), "An analytic formula for
  entraining CAPE in midlatitude storm environments", *J. Atmos. Sci.*,
  doi:10.1175/JAS-D-23-0003.1;
- including a Bernoulli-like form.

## Code

**`entraining_plume` computes all three in the same pass, always.** They are
diagnostic fields on `Plume`:

- `cape_first_lnb`, `cape_positive_area` and `cape_maximum_work`;
- `first_neutral_layer` and `maximum_work_layer`;
- `internal_neutral_layers` and `internal_neutral_depth`: the count and
  thickness of layers with b ≤ 0 between the LFC and the top, excluding the
  top.

**A `cape_definition` option selects which one is `plume.cape`,** the value
the trigger, the closure (M_b and Γ_M) and every published CAPE read.

- The options are "first_lnb" (the default), "positive_area" and
  "maximum_work".
- It is threaded through `entraining_plume`, `cloud_base_mass_flux`,
  `convect`, the convection module (`convective_cape_definition`) and the
  configuration.
- Every plume `convect` launches uses it: the closure's displaced plume, the
  CAPE after convection, and the fractional crossing's trials.

**"first_lnb" is today's code, bit for bit.** Its CAPE is accumulated by the
existing statements, unchanged.

**Per window, convection also publishes two arrays,**
`cape_window_start_positive_area` and `cape_window_start_maximum_work`, each
window's starting value under D1 and D2. They sit beside the existing
`cape_window_start`, which follows the selected definition.

**The trigger is held at 70 J/kg in every case, recorded as the user worded
it:**

> The numerical threshold is held constant to isolate the CAPE-definition
> effect; it is not assumed to represent an equivalent physical trigger
> under the new metric.

## Stage S1: diagnostic only (nothing driven by D1 or D2)

**S1a, identity.** With "first_lnb", explicit runs at 1800/300 and 1800/150
reproduce `ox_300` and `ox_150` to every printed digit. The new outputs are
additions.

**S1b, the census.** That 1800/300 run, with `entraining_plume` wrapped to
record every launch over days 10-20, without changing anything. For each
launch:

- the three CAPEs;
- the first neutral layer, the maximum-work layer and the top;
- the internal neutral layers' count and depth;
- for each internal neutral layer, the positive work accumulated below it and
  the negative work across it.

**Window to window,** from the published per-window arrays, the
distribution of |ΔC| between consecutive window starts, for each definition.

**S1c, continuity at the 301 crossings.** D1's diagnostic run (`d1.py`'s
setup, the fractional run at 300 s) is repeated with the new fields
recorded. At each unresolved crossing, the three definitions are evaluated at
f_lo, at f_lo + 1e-9 and at f_lo − 1e-9.

- A definition is continuous at a crossing if |C(f_lo + 1e-9) − C(f_lo)| ≤
  1e-6 J/kg. The smooth slopes seen were about 50 J/kg per unit f, which is
  5e-8 J/kg over the 1e-9 step.
- **D0 is expected to jump** in every case, as measured.

**S1d, D2's consistency with the plume, at every census launch:**

- the maximum-work layer lies between the LFC and the layer below the top;
- W at the top is ≤ 0;
- D2 equals the maximum of the plume's own running W, exactly. The census
  recomputes W from the plume's buoyancy and dz independently.

**S1e, vertical refinement and floating-point zero, a unit test:**

- **The sounding:** a synthetic conditionally unstable sounding with one
  internal layer at near-neutral buoyancy, evaluated at n and 2n layers. The
  2n case splits every layer into two, with temperature and humidity
  interpolated in ln p.
- **The perturbation:** that layer's temperature is moved by ±ε through
  neutral, with ε = 1e-6 K.
- **Passes if at both resolutions:**
  - D0 jumps;
  - D1 and D2 change by at most 1e-3 J/kg;
  - no D1 or D2 value depends on whether a buoyancy is exactly 0.

**S1 passes if** S1a, S1d and S1e pass, and D1 and D2 are continuous at all
301 crossings (S1c). Otherwise D1 or D2, whichever fails, is not driven in
S2. That is recorded.

## Stage S2: driven by each definition

**Runs, for D1 and D2:**

- the explicit closure, "stop" continuation, host 1800 s;
- Δt_conv = 600, 300, 150, 75 and 37.5 s;
- 20 days, read over days 10-20;
- otherwise F0's frozen-radiation setup.

**D0's comparators are the existing runs:** `ox_600` to `ox_37.5`, with
`ox_18.75` to `ox_4.6875` reported.

**Quantities:**

- time-mean CAPE and post-convection CAPE, under the run's own definition;
- D_t;
- the temperature and RH profiles, 100-850 hPa, and RH 300-700;
- **the lapse rate:** (T_850 − T_400)/(z_400 − z_850). It uses the layers
  nearest 850 and 400 hPa in the late-mean profile, with heights from the
  hypsometric equation on that profile. The surface and the boundary layer
  are excluded;
- 𝓜, subsidence and 𝓓(p);
- precipitation;
- FE, SE, PE and α_total,unique.

**Reported, not gated:** the window-to-window |ΔC| distribution (99th
percentile and maximum) under each definition.

**U1's conservation and provenance identities must pass in every run.**

**Convergence, by A6's rule and F4's classification:**

- on 600/300/150 s and on 150/75/37.5 s;
- profiles judged layer by layer;
- **non-monotonic** as in T2, over 600-37.5 s.

## Decision rules, frozen

| verdict | condition |
|---|---|
| **D2 supported** | All of: D2 passes S1, with continuity at all 301 crossings and S1d and S1e passing. U1 passes in every D2 run. No quantity in the list, the lapse rate and RH 300-700 included, is a changed equilibrium on either three-point ladder, and none is non-monotonic. No limiter is used: none exists in the code. |
| **D1 supported over D2** | D1 meets the same conditions and D2 does not. Or the plume's equations are shown to represent negative buoyancy elsewhere, so D2 would count it twice. **Read before freezing: they do not.** The plume's only use of negative buoyancy is its termination test on the same W that D2 maximises. |
| **Neither supported** | Both fail. D0 is kept only as the historical behaviour. It is not restored as a preference: its discontinuity at an exactly neutral layer is a measured defect. |

**If D2 is supported,** the recorded statement is:

> D2 removes the first-neutral-layer discontinuity and gives a closure whose
> CAPE follows the plume's own energy path.

**Adoption, and any trigger protocol for the new metric,** are the user's
decisions. The 70 J/kg threshold is not validated for D1 or D2 by this
protocol.

## What this cannot remove

**Other discontinuities remain under every definition:**

- **The source layer.** It is the arg-max of moist static energy, and can
  switch.
- **The LFC.** It is the first positive layer.
- **The termination,** where W itself crosses zero.

The last is the plume's own dynamical limit. Under D2, CAPE and the plume's
transport change there together, which is consistent. None of these appeared
among the 301 crossings, and the S1b census reports how often each occurs.

---

## Result (appended; nothing above this line changed)

### S1: passes for D1 and D2

**S1a, identity.** With "first_lnb", the reruns reproduce `ox_300` (881
shared values) and `ox_150` (905) with no difference. S2's five D0 runs
reproduce `ox_600` to `ox_37.5`, 869-1015 values each, likewise.

**S1b, the census.** 1576 plume launches in the late half of the 1800/300 s
run.

| | D0 | D1 | D2 |
|---|---:|---:|---:|
| mean CAPE, J/kg | 85.4 | 310.8 | 304.1 |

**Every one of the 1576 launches has internal non-buoyant layers** between
its LFC and its top:

- 6.3 layers on average, with a median total depth of 5.5 km;
- each crossed with a median −31 J/kg of work, the least favourable −126;
- the work accumulated below them has a median of 151 J/kg.

**The plume normally rises through 5 km of non-buoyant air, on accumulated
work.** D0 ends at the first such layer.

**Changes between consecutive launches:**

| what changed | count |
|---|---:|
| source | 0 |
| LFC | 0 |
| LNB | 659 |
| top | 198 |
| maximum-work layer | 52 |

**Window to window, |ΔC| over 25 J/kg:**

| | D0 | D1 | D2 |
|---|---:|---:|---:|
| windows | 337 | 52 | 89 |

The largest |ΔC| is 429 J/kg under all three.

**S1c, continuity at the 301 crossings:**

| | continuous | largest step across 1e-9 of f |
|---|---|---:|
| D0 | 0 of 301 | 427 J/kg |
| D1 | **301 of 301** | 3.1e-8 J/kg |
| D2 | **301 of 301** | 3.3e-8 J/kg |

The crossing run reproduces `of_300` in all 848 values.

**S1d holds at every census launch:**

- D2 equals the maximum of the independently recomputed running work,
  exactly;
- the maximum-work layer lies between the LFC and the layer below the top;
- W at the top is ≤ 0.

**S1e passes at 30 and 60 layers.** It is the unit test
`test_maximum_work_and_positive_area_are_continuous_through_a_neutral_layer`:

- D0 jumps by more than 10 J/kg across a 2e-6 K nudge through neutral;
- D1 and D2 move by at most 1e-3 J/kg;
- the source, LFC and top are unchanged.

**Both D1 and D2 pass S1 and were driven in S2.**

### S2: neither supported, as frozen

**All fifteen runs completed** 20 days and pass U1. Residuals are ≤ 7.5e-16,
F6a and F6b have no failures, there are no ordering violations and no
negative class. Each is invalid on the validity-domain gate only, as every
frozen-radiation run is.

| Δt_conv | | D0 | D1 | D2 |
|---:|---|---:|---:|---:|
| 600 s | time-mean CAPE | 47.3 | 33.7 | 37.7 |
| 300 s | | 51.8 | 37.8 | 44.2 |
| 150 s | | 57.1 | 46.9 | 54.2 |
| 75 s | | 75.1 | 70.9 | 74.2 |
| 37.5 s | | 93.5 | 78.9 | 118.3 |
| 600 s | D_t | 0.125 | 0.070 | 0.101 |
| 150 s | | 0.161 | 0.115 | 0.172 |
| 37.5 s | | 0.317 | 0.240 | 0.462 |
| 150 s | RH 300-700 | 0.711 | 0.701 | 0.699 |
| 150 s | lapse rate 850-400, K/km | 7.125 | 7.102 | 7.088 |
| 150 s | 𝓜, kg m⁻²/day | 427.1 | 447.6 | 452.9 |
| 150 s | precipitation, mm/day | 4.172 | 4.259 | 4.261 |
| 150 s | PE | 0.961 | 0.973 | 0.967 |

**Failing quantities on either ladder (600/300/150, 150/75/37.5) or by
non-monotonicity:**

| | failing |
|---|---|
| **D0** | time-mean and post-convection CAPE and D_t (fine ladder); the temperature profile (coarse ladder) |
| **D1** | **D_t on the coarse ladder only:** 0.070, 0.0815, 0.115, a changed equilibrium |
| **D2** | time-mean and post-convection CAPE and D_t, on both ladders. CAPE's fine-ladder differences grow: 19.9, then 44.1 J/kg. |

- **D1 on the fine ladder:** every quantity converges or is phase-sensitive.
  That is the first fine ladder to do so under any definition. CAPE's
  differences shrink, 24.0 then 7.96 J/kg.
- **D2:** every state, transport and water-path quantity converges or is
  phase-sensitive on both ladders. That includes the temperature and RH
  profiles, RH 300-700, the lapse rate, 𝓜, subsidence, 𝓓, precipitation and
  FE/SE/PE. Its CAPE and D_t do not.
- **Window-to-window jumps above 25 J/kg remain under every definition.**
  Hundreds per run, largest about 380-430 J/kg. In the census they fell from
  337 under D0 to 52 under D1 and 89 under D2.

**Verdict, as frozen: neither supported.**

- D2 fails on its own CAPE and D_t.
- D1 fails on D_t on the coarse ladder.
- D0 is kept only as historical behaviour. It is not restored as a
  preference: its discontinuity at an exactly neutral layer is a measured
  defect, confirmed again (0 of 301 continuous).

### What the result says, and what it does not

**1. D1 and D2 remove the first-neutral-layer discontinuity,** as designed:
301 of 301 continuous.

**2. Both change the state, by about the same amount:**

- RH 300-700 0.70 against 0.71;
- the 850-400 hPa lapse rate 7.09-7.10 against 7.13 K/km;
- mass transport 5 percent larger;
- precipitation 2 percent larger.

**3. CAPE still depends on the substep under D2,** and D_t under both.

- **The large window-to-window changes that remain** are not this
  discontinuity: D1 and D2 are continuous through it.
- **The candidates, from the census (not measured window by window):**
  - the plume's termination, which changed between 198 consecutive launches
    in S1b. W crosses zero, the top moves, and D2's maximum is taken over a
    different range;
  - each host step's forcing impulse at the window that follows it.
- **This is a reading, not a result.**

**4. On the fine ladder, D1 converges in everything listed.**

- This suggests D1 may have a fine-step limit that D2 does not.
- A ladder below 37.5 s (18.75 and 9.375 s) would test it.
- The decision rule as frozen needs both ladders, so this is recorded and
  not acted on.

**5. The 70 J/kg trigger was held fixed as a control.** It is not validated
for D1 or D2. Under D1 and D2 the column's mean CAPE is 3.5 times D0's for the
same state, so the same 70 J/kg is a different, lower relative threshold.

**Not done:**

- no default changed;
- no definition adopted;
- F0 not qualified.

## User disposition (2026-10-14), appended

**Recorded as the user worded it:**

> D1 and D2 successfully remove the first-neutral-layer discontinuity,
> confirming that D0 is discontinuous and geometrically inconsistent with the
> modeled plume extent. Neither alternative qualifies under the frozen closure
> gates. The remaining CAPE jumps are not caused by internal neutral layers
> and require separate attribution before further definition or trigger
> changes.

**The trigger, recorded as the user worded it:**

> The 70 J/kg threshold is defined only for the historical first-LNB CAPE
> metric. It is retained unchanged in D1/D2 experiments solely to isolate
> CAPE-definition effects. It is neither adopted nor validated for either
> alternative definition.

**No mapping such as 3.6 × 70 is adopted.** The D1/D2-to-D0 ratio is a
distribution, not a constant.

**The current state:**

| item | state |
|---|---|
| D0 | the default, retained for identity only. Its discontinuity is a measured defect. |
| D1 | the most promising experimental definition, but not supported |
| D2 | dynamically consistent, but its own CAPE and D_t are unresolved |
| 70 J/kg | a legacy control value under D1/D2 |
| F0 | not qualified |
| interactive radiation, P1-P3 | blocked |
| next | the jump census (`convective-cape-jump-census.md`) |
| D1 extension | conditional on the census |

**The sequence:** jump census and counterfactual attribution → conditional D1
fine-ladder extension → definition selection → a separate trigger protocol.

**The later trigger protocol** sets any new threshold from a
definition-independent physical event. It is never fitted to RH, the lapse
rate, precipitation, surface temperature or D_t.
