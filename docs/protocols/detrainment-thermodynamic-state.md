# Protocol: the thermodynamic state of detrained air

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards.

**Authorised by the user (2026-10-13)** after the diagnostic branch
(`mass-flux-diagnostic-branch.md`). It is bounded. No production conclusion is
drawn from it, no default changes, and P1-P3 stay blocked.

## The disposition this protocol starts from, as the user worded it

> Formally not classified under the preregistered threshold; scientifically
> inconsistent with the trigger hypothesis and supportive of
> detrainment-dominated moistening.

- **The trigger ladder's frozen verdict stands unchanged.** This wording is
  appended to that protocol as an addendum, below its results.
- **A6 stays failed,** classified as a changed attractor. No repair of A6 is
  attempted until this protocol is complete. If detrainment pins the
  environment near saturation, small differences in event timing can select
  slightly different near-saturated equilibria, and a changed closure may
  alter that.
- **The 1 hPa equilibrium is a separate, stratospheric-radiation defect.** It
  is not handled here, and upper-stratospheric convective adjustment is not
  restored.

## Hypothesis

The free troposphere is too humid because warm, saturated plume vapour is put
directly into a cooler environment.

If detrained material is instead equilibrated conservatively at the
environmental or mixed temperature:

- the plume's excess water becomes condensate, not environmental vapour;
- compensating subsidence can then lower RH.

## What the operator already does: read before this was frozen, stated so the result can be read against it

**1. The operator already carries the plume's enthalpy.** Detrainment adds
D·s_u to the receiving layer together with D·q_u (`mass_flux._operator`). The
plume's warmth therefore enters the layer along with its vapour.

**Mixing a small fraction f of saturated plume air into an environment at RH < 1:**

- q rises by f·(q_s(T_u) − q_e);
- q_s rises by about f·(dq_s/dT)·(T_u − T_e);
- these cancel to first order in the temperature excess;
- what remains is a moistening f·(1 − RH)·q_s,e;
- plus a supersaturation of second order, ½·f·q_s''·(T_u − T_e)², from the
  curvature of q_s.

**So D1 is predicted to change little,** unless D0's post-detrainment
excess, measured below, is a material fraction of the free troposphere's
condensation. Without the warmth that travels with the vapour, the large
excess the hypothesis describes would not appear.

**2. The trigger ladder's decomposition counted the plume's own condensate on both sides.**

- Convection's labelled vapour tendency is `transport = (q − q0)·m + formed`,
  so the plume's in-layer condensate is counted as transport.
- The water module's label `evaporation_minus_condensation` includes the same
  condensate, less what evaporates, as condensation.
- So "+1.53 transport, −1.59 condensation" contains the plume's own rain-out,
  added and removed.
- It does not by itself show that detrainment outpaces subsidence.

**This protocol separates the terms (§ Diagnostics).** That
reading is recorded as an interpretive erratum in
`mass-flux-diagnostic-branch.md`, keeping both readings.

**These are predictions, not criteria.** The criteria below are the user's,
fixed numerically here.

## Setup (identical to T0 of mass-flux-diagnostic-branch.md)

- **Configuration:** R4 with `mass_flux_convection=True` and
  `convective_target="dry_adiabat"`, from R4's saved state
  (`gcl_r4_200000_state.npz`).
- **Radiation:** the common prescribed radiation, solved once on R4's state
  and replayed unchanged at every step. The check that it is identical value
  for value is repeated.
- **Length and window:** 480 steps at 3600 s; results are read over days 10-20
  (steps 240-480).
- **Held fixed, exactly as T0:**
  - cloud-base closure M_b = CAPE_d/(τ_c Γ_M), τ_c = 3600 s;
  - trigger 70 J/kg;
  - ε = δ = 1 km⁻¹;
  - the plume's mass-flux profile, source layer and termination;
  - the subsidence;
  - α = 0.67;
  - surface exchange;
  - the initial state and timestep;
  - the A1-A5 conservation gates.
- **Γ_M** is measured with the case's own operator: the closure is defined as
  "the rate the same operator removes CAPE per unit mass flux". Its formula
  does not change.
- **The T-ladder runs were invalid on the mass-flux record's domain** (the
  plume reached 143 K under frozen radiation). These cases are compared
  within that configuration, as T0-T3 were. A domain refusal of the same kind
  is reported, not treated as a failure of a case.

**What changes is the state of detrained material only.** The plume
(its ascent, buoyancy, condensate and CAPE), its detrainment amounts D, and its
entrainment amounts E are the same functions of the environment in every case.

## Cases

**D0, the existing implementation.**

- The receiving layer gets D·s_u and D·q_u.
- Any supersaturation this leaves is condensed by the water module's
  end-of-step saturation adjustment and rains out.
- This is exactly T0.

**D1, conservative equilibration in the receiving cell.**

- After each substep of the operator, every layer with D > 0 that is above
  saturation is brought to saturation at constant moist static energy,
  h = s + L·q, by the same Newton solve the plume uses (`_saturate`).
- The condensate c_det:
  - is booked explicitly as detrained condensate;
  - joins the plume's condensate formed in that layer;
  - follows the existing convective pathway: the fraction α evaporates where
    it formed, cascading downward under each layer's deficit, and the rest
    reaches the surface.
- Its latent heat stays in the layer, in s.
- No water or energy leaves except as surface precipitation.

**D2, the conservative environmental-temperature limit.**

- In each layer with D > 0, the detrained vapour is capped:
  q_det = min(q_u, q_s(T_e, p)), with T_e the layer's temperature at the
  start of the substep.
- The excess x = M_b·Δt·D·(q_u − q_det) is detrained condensate. Its
  enthalpy is retained: the detrained dry static energy is
  s_u + L·(q_u − q_det), so D·h_u enters the layer exactly.
- x joins the plume's condensate in that layer and follows the same
  α pathway as D1.
- **D2 is a diagnostic limit, not proposed physics.**

**In all three cases:**

- the column sum of the energy redistribution is zero;
- the column sum of vapour transport is zero;
- net condensation sums to the surface precipitation;
- the A1-A5 unit tests pass in every mode.

## Diagnostics

Every quantity below is per layer and accumulated over substeps.
Convection publishes them as rates.

**The operator in advective form.** It splits exactly, because
M_{k+1} − M_k − E_k = −D_k:

- **subsidence:**
  - vapour: M_b·Δt·M_{k+1}·(q_{k+1} − q_k);
  - dry static energy: M_b·Δt·M_{k+1}·(s_{k+1} − s_k);
- **net detrainment:**
  - vapour: M_b·Δt·D_k·(q_det,k − q_k);
  - dry static energy: M_b·Δt·D_k·(s_det,k − s_k).

**Detrained material:**

- **detrained vapour, gross** (F_det^v): M_b·Δt·D·q_det;
- **detrained condensate** (F_det^c): zero in D0 by construction, x in D2,
  c_det in D1;
- **post-detrainment excess** (C_post-det):
  - D0 and D2: the supersaturation each substep leaves in layers with D > 0,
    measured at constant h and not removed by convection;
  - D1: c_det, which removes it;
- **re-evaporation and formation:**
  - α evaporation (W_re-evap);
  - the plume's condensate formed (the in-plume condensate).

**Read from the water module's labels:**

- its own large-scale condensation in free-tropospheric layers, as
  −(evaporation_minus_condensation) − c;
- diffusion, as `vapour_eddy_diffusion`.

**Column integrals over days 10-20,** both for the 300-700 hPa layers and for
the whole free troposphere above the plume source:

- W_det,vapor;
- W_det,condensate;
- W_post-det;
- W_large-scale condensation;
- W_precip;
- W_re-evap.

**The RH budget,** d ln RH over 300-700 hPa, mass-weighted, days 10-20. Its
terms:

- **d ln q:**
  - subsidence;
  - net detrainment;
  - α evaporation;
  - detrained-condensate evaporation;
  - large-scale condensation;
  - diffusion.
- **minus d ln q_s:**
  - radiation;
  - subsidence heating;
  - detrainment heating;
  - latent heat of net convective condensation;
  - latent heat of large-scale condensation;
  - evaporative cooling.
- **Grouped as the user named them:**
  - ΔRH_det, the net-detrainment and detrained-condensate terms;
  - ΔRH_subs, the subsidence terms;
  - ΔRH_rad, radiation.
- **The reconstruction is reported against the actual change,** as in T0.

**Also reported:**

- mean RH 300-700;
- duty cycle;
- mass transport;
- subsidence span;
- precipitation;
- mean CAPE_d;
- validity status and residuals.

## Decision criteria, frozen

"Material" means 10 percent of D0's value over the same window and layers.
An RH change counts only at 0.03 or more in the 300-700 hPa mean, the
trigger ladder's threshold.

**Premise check, read first.**

- **The condition:** D0's post-detrainment excess over 300-700 hPa, days
  10-20, is below 10 percent of D0's large-scale condensation there.
- **If it holds:**
  - D1 cannot test the hypothesis, and is reported as "not a test", whatever
    its RH;
  - the hypothesis as stated, supersaturation created by detrainment, is
    then not the mechanism in this operator;
  - D2 is still classified below, as the limit it is.

**Supported, by a case (D1 or D2) that meets all of:**

1. gross detrained vapour over 300-700 hPa falls by ≥ 10 percent from D0;
2. detrained condensate is carried explicitly and is > 0;
3. post-detrainment excess plus large-scale condensation over 300-700 hPa
   falls by ≥ 10 percent from D0;
4. mean RH 300-700 falls by ≥ 0.03 from D0, and the RH budget attributes the
   fall to the detrainment terms (ΔRH_det less moistening than D0's);
5. water and energy close: the conservation gate passes, with residuals as
   T0's (≤ 1e-12);
6. δ, ε and the trigger are unchanged (true by construction, and checked
   from the published parameters).

**Rejected, by a case that:**

- changes post-detrainment excess plus large-scale condensation over
  300-700 hPa by ≥ 10 percent from D0;
- leaves mean RH 300-700 within 0.03 of D0.

The next suspect is then the amount and vertical distribution of detrainment
(δ = ε), and a δ(z) protocol would be justified. That protocol is not part of
this one.

**Anything else is reported as found, without classification.**

**Not required:** agreement with RH = 0.81. It is a theoretical comparison,
not a target. RCE models spread in humidity, and precipitation efficiency is a
published source of that spread.

## Code, introduced by this protocol only

- **`detrained_state`** on `mass_flux.convect`, `apply_mass_flux` and
  `cloud_base_mass_flux`, and the convection option of the same name:
  - "plume" (D0, the default and the existing behaviour, bit for bit);
  - "mixed_equilibrium" (D1);
  - "environment_saturation" (D2).
- **The decomposition fields** on `MassFluxStep`, and their published rates.
- **Tests:**
  - A1-A5 in every mode;
  - D0 identical to the code before this change;
  - D1 leaves no layer with D > 0 above saturation;
  - D2 detrains no vapour above q_s(T_e);
  - the advective split sums to the operator exactly.

---

## Result (appended; nothing above this line changed)

**Setup reproduced.**

- **D0 is T0:** RH 0.93720, the same mass transport, precipitation and
  duty cycle.
- **The replayed radiation was identical** value for value in all three
  cases (checksum -2646.2636).
- **D0's code is bit-for-bit the code before this change.** A digest of
  `convect` on three soundings and of six coupled steps matched one from a
  worktree of the previous commit.
- **All three runs closed water and energy to 4.9-6.0e-16.**
- **Validity:**
  - all three are invalid on the validity-domain gate only, as T0-T3 were;
  - the conservation gate passed;
  - δ, ε, the trigger, τ_c and α read back unchanged from each run's
    published parameters.

**The split closes.** Inside 300-700 hPa, the published terms equal the
labelled totals to 1e-15:

- for vapour, -0.006739 in both;
- for temperature, through d ln q_s, -0.043387 in both.

Globally the vapour check does not close. The convection label also carries
the dry adjustment's own vapour mixing near the surface, which no band here
contains.

### The free troposphere's budget, days 10-20, 300-700 hPa

**Water, kg m^-2 over 10 days:**

| | D0 | D1 | D2 |
|---|---:|---:|---:|
| detrained vapour, gross | 410.15 | 409.78 | 362.49 |
| net detrainment, D (q_det - q) | +58.89 | +58.66 | +44.70 |
| subsidence | -111.81 | -111.73 | -103.14 |
| re-evaporation of convective condensate | +56.19 | +56.38 | +59.59 |
| plume condensate formed in these layers | 81.34 | 81.46 | 69.72 |
| detrained condensate | 0 | 0.149 | 17.30 |
| post-detrainment excess | 0.066 | 0.149 | 0.022 |
| large-scale (water-module) condensation | 3.034 | 2.891 | 1.937 |
| diffusion | +1.05 | +1.05 | +1.09 |

**d ln RH by mechanism.** Each vapour term is set against its own heating:

| | D0 | D1 | D2 |
|---|---:|---:|---:|
| subsidence (drying and warming) | **-7.03** | -7.02 | -6.46 |
| net detrainment, with its heating (ΔRH_det) | +0.74 | +0.69 | **-0.60** |
| re-evaporation, with its cooling | **+5.66** | +5.67 | +6.27 |
| radiation (ΔRH_rad) | +1.27 | +1.27 | +1.27 |
| large-scale condensation, with its latent heat | -0.57 | -0.54 | -0.44 |
| convective latent heat, net of transport | -0.08 | -0.08 | -0.09 |
| diffusion | +0.05 | +0.05 | +0.05 |
| sum / actual | 0.037 / 0.038 | 0.040 / 0.041 | 0.006 / 0.008 |

**The other measures:**

| | D0 | D1 | D2 |
|---|---:|---:|---:|
| mean RH 300-700 | 0.9372 | 0.9385 | **0.8980** |
| duty cycle | 0.908 | 0.896 | 0.854 |
| mass transport, kg m^-2 | 7207 | 7200 | 6715 |
| subsidence span, kg m^-2 Pa | 2.88e8 | 2.88e8 | 2.69e8 |
| surface precipitation, kg m^-2 | 43.32 | 43.25 | 43.65 |
| mean CAPE_d | 382 | 379 | 315 |

### Verdicts, as frozen

**Premise check: it holds, so D1 is "not a test".**

- D0's post-detrainment excess over 300-700 hPa is 0.066 kg m^-2. That is
  **2.2 percent** of its large-scale condensation there (3.034), under the
  10 percent line.
- Over 100-850 hPa it is 2.4 percent.
- **Detrainment does not create supersaturation in this operator,** as the
  analysis frozen above predicted.
- D1 moved RH by +0.0013, and the excess plus large-scale condensation by
  1.9 percent. That is not material, and it is in any case not a test.

**D2: supported, by all six frozen conditions.**

1. Gross detrained vapour falls 11.6 percent (≥ 10).
2. Detrained condensate is carried explicitly, at 17.3 kg m^-2.
3. Excess plus large-scale condensation falls 36.8 percent, from 3.100 to
   1.959 (≥ 10).
4. RH falls by 0.039 (≥ 0.03). ΔRH_det goes from +0.74 to -0.60, so the
   fall is attributed to the detrainment terms.
5. Residual 6.0e-16, and the conservation gate passed.
6. δ, ε and the trigger are unchanged.

**Neither case meets the rejection rule.**

### What the result says, and what it does not

**1. The formal verdict and its mechanism point different ways.**

- **What D2 does:** it detrains 11.6 percent less vapour by declaring the
  detrained air no moister than the environment's saturation. It does not
  remove supersaturation, because there was almost none to remove.
- **What it shows:** humidity is sensitive to the amount of vapour detrained
  air carries.
- **What it does not show:** that the operator mishandles a supersaturated
  mixture.
- **D2 remains a diagnostic limit,** as declared. Physically, the warm, moist
  detrained air would have to give up water without being cooled to the
  environment's temperature.
- **Subsidence fell 7 percent in D2, and RH still fell.** So it was not
  subsidence that dried it.

**2. The dominant moistening term is the re-evaporation of convective condensate, not detrainment.**

- **In D0:**
  - subsidence -7.03;
  - re-evaporation with its cooling +5.66;
  - radiation +1.27;
  - net detrainment +0.74.
- **Re-evaporation's mechanics:**
  - The fraction α = 0.67 of the plume's condensate evaporates in the layer
    where it formed, or cascades to the layers below. It cools them as it
    does.
  - It puts 56 kg m^-2 into 300-700 hPa, about as much as net detrainment
    (59).
  - Its cooling doubles its effect on RH.
- **In D2,** the detrained condensate adds to re-evaporation (+6.27). That
  offsets part of D2's drying.
- **This confirms the interpretive erratum** in
  `mass-flux-diagnostic-branch.md`:
  - T0's "condensation -1.59" was mostly the plume's own condensate in the
    bookkeeping;
  - the water module's large-scale condensation is small, 3.0 kg m^-2
    against 56 re-evaporated;
  - "condensation is cleaning up supersaturation created by detrainment" is
    not what happens here.

**3. A property of the α pathway was found while testing D1, not acted on.**

- `deposit_evaporation` limits each layer's evaporation by its deficit at
  the temperature before that evaporation's own cooling.
- A layer it evaporates into can therefore end slightly above saturation
  (0.2 percent in the unit sounding).
- The water module's adjustment then condenses the overshoot.
- It is unchanged in every case, and it is conservative. It is a candidate
  for part of the 3 kg m^-2 of large-scale condensation, not measured apart.

**Not concluded:**

- No production conclusion, and no default changes.
- **The rule for the rejected outcome** pointed to δ(z). The supported
  outcome names no successor.
- **The data point to the α pathway** (where, and how much, convective
  condensate re-evaporates) as the larger term. That is for the user to
  decide.

### Addendum: the user's recording of D2 (2026-10-13)

The frozen classification stands: **supported**. All six criteria were met.
The causal explanation is recorded as the user worded it:

> D2 is supported as a vapor-partition sensitivity, not as a
> supersaturation-removal mechanism. Its reduction in free-tropospheric
> relative humidity arises primarily from a 12% reduction in vapor transferred
> from the plume to the environment. The diagnosed post-detrainment
> supersaturation was too small to explain the response.

**The user's decision on what follows:**

- the next protocol targets the re-evaporation of convective condensate
  (`convective-reevaporation-placement.md`), not δ(z);
- D2 is not combined with it.
