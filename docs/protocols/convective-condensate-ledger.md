# Protocol: what α counts, and placement at a fixed amount

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards.

**Authorised by the user (2026-10-13)** after
`convective-reevaporation-placement.md`. It is bounded and diagnostic. No
production conclusion is drawn, and no default changes.

- **Blocked until it is done:** P1-P3, and any reading of RH ≈ 0.71 as
  realistic.
- **D2 is not combined with Stage A.**
- **The stratospheric protocol** stays independently necessary. It does not
  block this one, because every run here uses the same frozen common
  radiation.

## The E result, recorded narrowly as the user directs

- **What E2 and E3 support:** the evaporation treatment materially controls
  free-tropospheric humidity.
- **What they do not establish:** that vertical placement is the cause.
  Neither case held the evaporated amount or the convective state fixed. The
  jump from about 0.94 to 0.71 is not yet a placement effect.
- **What E1 showed:** E0's large-scale condensation was mostly recycling. It
  was evaporation forced past the receiving air's capacity, then condensed
  and rained out. It was not condensation from resolved cooling or ascent.
  - The sequence: forced evaporation → vapour and cooling → supersaturation
    → condensation → precipitation.
  - The condensation term was concealing an evaporation-condensation loop.

**On the candidates:**

| candidate | standing |
|---|---|
| E0 | rejected as a production candidate: it forces that loop |
| E1 | a diagnostic only |
| E2 | the rain-evaporation component |
| E3 | a limiting placement diagnostic, not a closure |

## What the measured α is, from its source, read before this was frozen

`EVAPORATING_CONDENSATE_FRACTION = 0.67` is 1 − PE from Langhans, Yeo &
Romps (2015, J. Atmos. Sci. 72, 1045).

**Their definitions** (their eq. 1 and text):

- **PE** is "the fraction of condensed particles that subsequently rains
  out". Equivalently, it is the probability that a water molecule reaches
  the surface given that it condensed at least once.
- **DR = CE × PE = CE × FE × SE:**
  - FE, the formation efficiency, is the probability that condensate is
    converted to precipitation;
  - SE, the sedimentation efficiency, is the probability that precipitation
    reaches the surface.

**Their Table 1:**

| simulation | PE | FE | SE |
|---|---:|---:|---:|
| REF | 0.33 | 0.54 | 0.62 |
| MOIST | 0.38 | 0.69 | 0.55 |

**Therefore the measured 0.67 is total condensate return,** α_total = 1 − PE.
It is not falling-rain evaporation.

**Their split, as fractions of condensate:**

| | REF | MOIST |
|---|---:|---:|
| cloud condensate never converted (1 − FE) | 0.46 | 0.31 |
| rain evaporated in fall, FE·(1 − SE) | 0.21 | 0.31 |
| rain evaporated, per unit of rain formed (1 − SE) | 0.38 | 0.45 |

**This answers the definitional question by citation.** The model's job here
is to show which of its pathways can realise which part.

**What the citation does not settle:**

- **The regime.** These are cumulus congestus LES, not deep tropical
  convection. The authors say the efficiencies "might depend on the specific
  environment studied".
- **The counting.** Their PE counts molecules that condensed at least once.
  A gross formation that counts each recondensation again has a larger
  denominator, so it gives a smaller α.
- **The source of the dominant control.** They name SE as "the dominant
  control on PE for water entraining into the upper parts of the clouds".
  That is a statement about its vertical structure. It is not a statement
  that rain evaporation is most of the total.

**E2's 0.06,** read against this:

- It is rain evaporation divided by all condensate, with all condensate
  treated as rain.
- Its comparable is FE·(1 − SE) = 0.21-0.31, not 0.67.
- The category mismatch accounts for part of the factor of ten. A factor of
  3-5 remains in the rain law, to be measured, not assumed.

## The ledger

Every unit of convective condensate is tagged at formation, and accumulated
from the transfers themselves, never by differencing bulk water.

**C_form,gross** is condensate formed in the plume, counted each time it forms.

**Its destinations:**

C_form,gross = E_cloud + E_precip + P_surface + C_store

- **E_cloud** is evaporation before conversion to precipitation:
  - in the plume, when entrained air takes carried cloud water;
  - on detrainment into the receiving layer.
- **E_precip** is the evaporation of precipitation.
- **P_surface** is convective precipitation reaching the surface.
- **C_store** is condensate aloft at the end of a step. It is zero by
  construction: the plume is steady within a substep, and condensate leaves
  it by conversion or detrainment.

**Precipitation formed:** C_precip = C_form,gross − E_cloud.

**The quantities reported:**

- α_total = (E_cloud + E_precip) / C_form;
- 1 − FE = E_cloud / C_form;
- FE = C_precip / C_form;
- SE = P_surface / C_precip;
- α_rain = 1 − SE = E_precip / C_precip;
- PE = P_surface / C_form.

**Gross and net.** The water module's saturation adjustment can re-condense
water that convection just evaporated. **C_recycle** is the condensation the
convective evaporation causes in the same step. It is measured as a
counterfactual inside the water module:

- the end-of-step adjustment is solved a second time;
- that second solve uses the same provisional profile less this step's
  convective evaporation, with the evaporation's cooling removed;
- C_recycle = (own condensation) − (own condensation in that counterfactual).

Then:

- α_gross = (E_cloud + E_precip) / C_form;
- α_net = (E_cloud + E_precip − C_recycle) / C_form.

**How E0's pathway maps onto the ledger.** E0 has no conversion: all
condensate leaves the plume where it forms. Its evaporation is mapped by a
local-first rule, fixed here:

- in each layer, the layer's own α·C_form evaporates first, up to the
  limit: that is E_cloud;
- what evaporates of condensate that fell from above is E_precip.

## Stage A: what each pathway realises

**Setup.** Identical to E0 (`convective-reevaporation-placement.md`), with
`detrained_state="plume"`.

- R4's state;
- the common radiation, replayed;
- 480 steps of 3600 s, read over days 10-20;
- ε = δ = 1 km⁻¹, the trigger, τ_c, surface exchange.

**Aα0: E0 with the ledger.**

- The physics is unchanged, bit for bit.
- It reports the ledger, the local-first split, C_recycle, α_gross and
  α_net.

**Aα1: rain law only, no target.**

- Condensate leaves the plume where it forms, as E0.
- It falls under E2's rate law with no α cap: `"rate_law_only"`.
- With all condensate as rain, C_precip = C_form, and α_rain = α_total.

**The cloud pathway (`condensate_pathway="carried"`),** for Aα2 and Aα3.
It follows CAM5's Zhang-McFarlane `cldprp`, as in its WRF port.

- **Carried water.** The plume carries cloud water l. Entrained air brings
  none.
- **Mixing to equilibrium.** Each layer's mixture of moist static energy and
  total water (vapour and l) is brought to equilibrium by the plume's own
  solve:
  - where the equilibrium l exceeds the diluted inflow, the difference is
    C_form,gross;
  - where it is less, the difference is in-plume E_cloud.
- **Conversion.** Conversion is implicit, as `cldprp`:
  l_after = l_eq / (1 + c0·Δz), and the difference is precipitation formed
  in that layer.
- **c0 = 0.045 m⁻¹** (`c0_ocn` in that port), nominal, since the column has
  an ocean surface. **c0 = 0.0059 m⁻¹** (`c0_lnd`) is the other end of the
  declared range, run alongside it.
- **Buoyancy** includes the carried water l_after as loading. That is the
  one change to the plume's buoyancy and therefore to CAPE. It is declared,
  and it follows from carrying the water.
- **Detrained water,** D·l_after, enters the receiving layer:
  - it evaporates up to that layer's wet-bulb capacity (E_cloud);
  - the rest becomes precipitation in that layer. This scheme carries no
    prognostic cloud water. In CAM it would enter the stratiform scheme;
    here that is declared, not hidden.
- **The energy and water bookkeeping is unchanged in form.** The plume's net
  condensation per layer is booked where it happens. The evaporated water
  enters the environment with its cooling.

**Aα2: cloud evaporation only.** The carried pathway, with no evaporation of
falling precipitation (`"none"`). At both values of c0.

**Aα3: combined, the production candidate.** The carried pathway, with E2's
rate law on the precipitation (`"rate_law_only"`) and no α target. At both
values of c0.

**Reported for every case:**

- the ledger, as column integrals and per layer, over 300-700 and
  100-850 hPa;
- α_total, 1 − FE, FE, SE, α_rain, PE, α_gross, α_net;
- C_recycle;
- the RH budget, with re-evaporation split into vapour and cooling;
- RH 300-700, duty cycle, mass transport, subsidence span, CAPE, surface
  precipitation and temperature, validity and residuals.

**The ledger closure:**

- C_form,gross = E_cloud + E_precip + P_surface, per step, to 1e-12
  relative;
- the water module's surface precipitation = P_surface + large-scale
  precipitation.

### Stage A verdicts, frozen

The comparison band is the source's own two simulations, widened by 0.05 on
each side:

| quantity | band |
|---|---|
| α_total | [0.57, 0.72] (around 0.62-0.67) |
| 1 − FE | [0.26, 0.51] (around 0.31-0.46) |
| FE·(1 − SE) | [0.16, 0.36] (around 0.21-0.31) |

**"Total condensate return" is realisable:**

- **The condition:** Aα3, at either c0, gives α_total in the band, with
  E_precip/C_form ≤ ½·α_total.
- **The conclusion:** the measured 0.67 can be met by the model's combined
  pathways as a total-return fraction, with the cloud pathway carrying most
  of it.

**Rain-only reproduces it:**

- **The condition:** Aα1 gives α_total ≥ 0.57.
- **The conclusion:** the measurement could be met as rain evaporation
  alone.

**Not operationally defined in this model:** neither of the above.

- 0.67 is then not used as a coefficient in this column.
- Each pathway's realised fraction is reported against its band.

**The pathways, read separately:**

- Aα2's E_cloud/C_form against the 1 − FE band;
- Aα1's and Aα3's E_precip/C_form against the FE·(1 − SE) band.

Each is classified as below, within or above its band.

**Production acceptance is not available from Stage A.** The user's
acceptance list requires:

- timestep convergence;
- radiation states that are not upper-atmosphere-invalid.

These runs use frozen radiation and are invalid on the validity domain, as
every run in this branch has been. Aα3 can at most be named the candidate,
with each acceptance item marked met, not met or not tested.

## Stage B: placement at a fixed amount (diagnostic replay)

**The control.** E0 (= Aα0) is rerun recording, at every step, the mass-flux
scheme's result: its transport, energy, condensate formed and evaporation,
per layer.

**The replay.** Each Stage B run replays that record in place of the
mass-flux scheme:

- the same convective transport;
- the same energy redistribution;
- the same condensate formed, and the same duty cycle, mass flux and
  subsidence.

**The one change is where the step's evaporation goes.** Its total, E0's
own per step, is redistributed by a profile, and the net condensation
handed to the water module follows it. The energy field is independent of
the evaporation, so the evaporation's cooling moves with it. Everything
else stays interactive:

- the dry adjustment;
- the water module;
- the surface;
- diffusion.

**The radiation** is the common replay.

**Profiles,** each normalised to E0's step total:

- **Pα0, E0's own profile.** Its B1 run must reproduce E0 exactly, which
  checks the replay.
- **Pα1, formation level:** proportional to the condensate formed in each
  layer.
- **Pα2, falling path:** proportional to E2's rate-law weights,
  (1 − q/q_s)⁺·√F_k·m_k. F is the replayed condensate formed above, and
  q and q_s are the current state.
  - A step whose weights sum to zero uses Pα1's profile instead.
  - Such steps are counted.
- **Pα3, below the free troposphere:** proportional to layer mass over the
  layers centred below 850 hPa.
  - **This is a declared deviation** from "below cloud base". On this grid
    cloud base is the second layer, so "below cloud base" would be the one
    layer that failed E3.

**Modes:**

- **B1, matched total:** the profile is applied as given, even where it
  supersaturates. The water module re-condenses the excess, and C_recycle
  measures it.
  - The amount evaporated is identical across profiles by construction:
    E_gross and L·E_gross.
  - This answers causality, and is not a closure.
- **B2, capacity-limited:**
  - each layer takes no more than its wet-bulb capacity on the
    post-transport state;
  - the excess falls to the next layer below under the same limit;
  - what no layer can take reaches the surface.
  - This answers feasibility, and the totals then differ.

**B1 and B2 are not combined in any reading.**

**A replay that over-draws a layer's vapour is refused,** as any over-draw
is. It is reported, not repaired: the replayed transport was computed on
E0's state.

### Stage B verdicts, frozen

**Placement independently supported:**

- **The condition:** a B1 profile gives a 300-700 hPa RH differing from
  Pα0's by ≥ 0.03, with identical E_gross.
- **The difference is attributed** to the vapour term, the cooling term and
  recondensation, read from the RH budget and C_recycle. Convection is
  identical by construction.

**Placement not supported:** every B1 profile is within 0.03 of Pα0.

**B2** is reported against B1's result for the same profile: RH, α_net and
the surface precipitation. It is not classified.

## Code, introduced by this protocol only

- **The condensate pathway.** `condensate_pathway` ("immediate", the
  default, bit for bit; or "carried") and `cloud_water_conversion_rate`
  (c0), on the plume, the operator, the closure and `convect`, and as
  convection and configuration options. The carried pathway accepts only
  `detrained_state="plume"` and the evaporations "none" or
  "rate_law_only", and refuses the others.
- **The evaporation options.** `precipitation_evaporation` gains "none" and
  "rate_law_only".
- **The ledger fields** on `MassFluxStep`, published as rates.
- **The recycling measure.** The water module's `diagnose_convective_recycling`
  computes the counterfactual adjustment and publishes C_recycle. It reads
  convection's published evaporation.
- **Tests:**
  - A1-A5 in the new modes;
  - "immediate" identical to the code before;
  - the plume's liquid budget telescopes;
  - the ledger closes per step;
  - with no evaporation the counterfactual gives C_recycle = 0;
  - an evaporation that supersaturates gives C_recycle > 0, bounded by it.
- **Stage B's replay** is diagnostic code in the run script. It wraps
  `convect` the way the common radiation wraps radiation, and is not a
  module option.

---

## Result (appended; nothing above this line changed)

**Checks.**

- **The default pathway is bit-for-bit the previous code.** The digest of
  `convect` and six coupled steps matched the E-state code.
- **All 14 runs that completed** closed water and energy to 4.9-6.3e-16.
  They are invalid on the validity-domain gate only, as every run in this
  branch has been. The conservation gate passed, and the radiation replay
  was identical at every step.
- **The ledger closes** at every step of every Stage A run, worst 3.9e-16
  relative.
- **Aα0 reproduces E0:** RH 0.93720.
- **The recording run reproduces Aα0,** and B1-Pα0 reproduces the recording
  run's RH series to every digit. That checks the replay.
- **Pα2 never needed its fallback.**

### Stage A: what each pathway realises (days 10-20, column)

| | α_total | E_cloud/C_form | E_precip/C_form | α_rain = 1 − SE | PE | C_recycle, kg m^-2 | α_net | RH 300-700 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Aα0, E0 | 0.670 | 0.655 | 0.015 | 0.044 | 0.330 | 3.40 | 0.642 | 0.937 |
| Aα1, rain law only | 0.062 | 0 | 0.062 | 0.062 | 0.938 | 0 | 0.062 | 0.711 |
| Aα2, cloud only, c0 0.045 | 0.014 | 0.014 | 0 | 0 | 0.986 | 0 | 0.014 | 0.705 |
| Aα2, cloud only, c0 0.0059 | 0.098 | 0.098 | 0 | 0 | 0.902 | 0 | 0.098 | 0.726 |
| Aα3, combined, c0 0.045 | 0.077 | 0.014 | 0.063 | 0.064 | 0.923 | 0 | 0.077 | 0.711 |
| Aα3, combined, c0 0.0059 | 0.152 | 0.098 | 0.054 | 0.059 | 0.848 | 0 | 0.152 | 0.735 |
| source (REF-MOIST), band ±0.05 | 0.62-0.67 | 0.31-0.46 | 0.21-0.31 | 0.38-0.45 | 0.33-0.38 | | | |

**In Aα1-Aα3, convection runs as in E2:**

- duty cycle 0.40-0.47;
- CAPE 118-142;
- surface 313.4-313.9 K.

**Verdict as frozen: not operationally defined in this model.**

- Aα3's α_total is 0.077 and 0.152, not within [0.57, 0.72].
- Aα1's is 0.062, not ≥ 0.57.
- So 0.67 is not used as a coefficient in this column.

**The pathways against their bands:**

- **Cloud,** Aα2's E_cloud/C_form: 0.014 and 0.098 against [0.26, 0.51].
  **Below.**
- **Rain,** E_precip/C_form: Aα1 0.062, Aα3 0.063 and 0.054, against
  [0.16, 0.36]. **Below.**
- **α_rain per unit of precipitation formed:** 0.06 against the source's
  0.38-0.45, which carries no band.

**Aα3 as the candidate, against the acceptance list:**

| requirement | status |
|---|---|
| conserves tagged and total water | met (6e-16; ledger 4e-16) |
| conserves energy | met |
| respects local condensate availability | met: evaporation never exceeds what falls in or is detrained |
| respects local thermodynamic capacity | met: the wet-bulb limit |
| emergent, defined α | met: 0.077 or 0.152 |
| no forced evaporation-recondensation | met: C_recycle = 0 |
| timestep convergence | **not tested** |
| radiation not upper-atmosphere-invalid | **not met** under the frozen radiation |

### Stage B: placement at a fixed amount

**B1, matched total.** E_gross is 80.287 kg m^-2 in every profile, with
identical convection: duty cycle 0.908 and mass transport 7207.

| | RH 300-700 | C_recycle | E_retained | evaporation mean p | surface P | T_s | budget sum / actual |
|---|---:|---:|---:|---:|---:|---:|---:|
| Pα0, E0's profile | 0.9372 | 3.40 | 76.88 | 550 hPa | 43.32 | 311.09 | 0.037 / 0.038 |
| Pα1, formation level | 0.9344 | 3.47 | 76.82 | 547 | 43.36 | 311.09 | 0.032 / 0.034 |
| **Pα2, falling path** | **0.7421** | 0.06 | 80.23 | 581 | 39.77 | 309.66 | 0.032 / 0.035 |
| Pα3, below 850 hPa | refused at step 87 | | | | | | |

**B2, capacity-limited:**

| | RH 300-700 | C_recycle | E_retained |
|---|---:|---:|---:|
| Pα0 | 0.9871 | 0.20 | 80.08 |
| Pα1 | 0.9871 (identical to Pα0 to every digit) | 0.20 | 80.08 |
| Pα2 | 0.7432 | 0.002 | 80.29 |
| Pα3 | refused at step 87 | | |

**Pα3 was refused, as the protocol provides.**

- **The refusal:** "ConfigurationError raised by water at step 87:
  convection's condensation drew a layer's vapour below zero."
- **The cause:** the replayed condensate formation in the upper layers was
  computed on E0's state. There it was offset in the same layers by local
  evaporation. With the evaporation moved below 850 hPa, those layers could
  no longer supply it.
- **It is reported, not repaired.** Pα3 has no result.

**Verdict as frozen: placement independently supported, by B1-Pα2.**

- **The size:** RH differs from Pα0's by 0.195, with identical E_gross,
  identical L·E_gross and identical convection.
- **Pα1 does not differ materially** (-0.003).

**Attribution.**

- **Convection** is identical by construction.
- **Recondensation** falls from 3.40 to 0.06 kg m^-2.
- **Both budgets reconstruct.** They are balances, so their terms describe
  how each equilibrium is held, not why the two differ:
  - re-evaporation's vapour +3.27 against +2.98, and its cooling +2.46
    against +2.68;
  - subsidence -7.96 against -7.03;
  - detrainment +1.04 against +0.74;
  - large-scale condensation 0.00 against -0.57.
- **What differs is where.**
  - Pα0 and Pα1 put most of the evaporation into layers that sit at
    saturation: 355-455 hPa, at RH 0.98-1.0. There 3.4 kg m^-2 of it
    re-condenses and rains out, and the upper half of the band stays
    saturated.
  - Pα2's weights are zero in saturated layers. They place the same water
    lower and in drier air, and 322-423 hPa falls to RH 0.61-0.68.

**B2 is reported, not classified.**

- **Pα0 and Pα1 coincide to every digit.** The capacity-limited cascade
  fills each layer to its wet-bulb capacity from the top down. Both profiles
  put enough into the upper layers that the outcome depends only on the
  step's total and the capacities. That is an inference from the identity,
  not measured apart.
- **Both are wetter than B1-Pα0** (0.987). Without the overshoot, there is no
  recondensation raining water out, as E1 found.
- **Pα2 under the capacity limit** (0.743) matches B1-Pα2.

### What the result says, and what it does not

**1. The measured 0.67 is total condensate return,** by its source's
definition. In this column, at CAM5 Zhang-McFarlane's own constants, the
model's physical pathways realise 0.08-0.15 of it.

- **Each pathway falls short of the source's split:**
  - cloud evaporation by a factor of 3-30. c0 converts the cloud water
    within a layer: an e-folding of 22 m at the ocean value, against layers
    hundreds of metres deep.
  - rain evaporation by a factor of 3-6, at ke = 1e-6.
- **The working hypothesis, "a category mismatch rather than a coefficient
  mismatch", holds in part.**
  - The category mismatch is real: 0.67 is total, and E2 was measuring rain
    alone.
  - But both categories are also under-realised here. That remains a
    coefficient or structure question.
- **E0 implements 0.67 as cloud return at the formation level.** Its
  local-first share is 0.655, larger than the source's 0.31-0.46. Its rain
  share is 0.015, far below the source's 0.21-0.31.

**2. Placement controls the humidity on its own, at a fixed amount.**

- Moving the same evaporated water, with the same cooling and the same
  convection, from saturated formation layers into subsaturated layers
  along the fall path lowers RH 300-700 from 0.937 to 0.742. It also ends
  the recycling.
- **This is the isolation E2 and E3 could not give.**
- **The RH of about 0.71 in E2, E3 and Aα1-Aα3 is therefore both amount
  and placement.** At the full amount, placement alone gives 0.742.

**3. Recycling is placement-driven.**

- It is 3.4 kg m^-2 where the evaporation goes into saturated layers.
- It is near zero where the evaporation follows subsaturation, at either
  limit.

**Not concluded:**

- no production adoption, no default changed, no constant tuned;
- **what amount is right:**
  - the source is cumulus congestus LES, and the counting differs;
  - CAM's c0 and ke are tuned values whose meaning against FE and SE is
    not established here;
- timestep convergence and interactive or valid radiation are untested;
- the stratospheric protocol is still unopened.

### Addendum: the user's disposition (2026-10-13)

**α, recorded as the user worded it:**

> The externally measured α = 0.67 represents total condensate return,
> 1 − FE × SE, in cumulus-congestus LES. It is not operationally realizable
> as a direct coefficient in the present column model. At the existing CAM5 ZM
> constants, the modeled cloud and rain pathways produce a combined
> condensate return of 0.077–0.152. Therefore, α = 0.67 shall not be imposed,
> tuned to, or used as an acceptance target.

**c0 and ke are not fitted.** The current constants are a baseline
hypothesis, not validated parameters. An amount protocol is deferred until
there are:

- regime-matched observations or LES targets; or
- several independent constraints beyond total condensate return.

**Placement, as the user records it:** "Vertical placement independently
controls free-tropospheric RH." RH 0.742 is not the model's accepted answer.
The placement effect is causal and large.

**Pα3, kept as a valid failed diagnostic:**

> Relocating all matched evaporation below 850 hPa was incompatible with the
> replayed upper-level water tendencies and produced negative vapor at step
> 87. The prescribed profile was therefore physically inadmissible under the
> protocol.

It is not repaired or rerun under this freeze.

**The identical capacity-limited results:**

> The capacity-limited E0 and formation-level profiles produced identical
> states. A top-down capacity-filling cascade is a candidate explanation, but
> was not instrumented in this protocol.

**The sequence that follows:**

1. F0, the fall-path candidate, under frozen radiation, with a timestep
   study.
2. The stratospheric protocol.
3. F0 with interactive radiation.
4. Reconsidering P1-P3, which stay blocked until then.
