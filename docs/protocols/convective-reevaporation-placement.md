# Protocol: where convective condensate re-evaporates

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards.

**Authorised by the user (2026-10-13)** after
`detrainment-thermodynamic-state.md`. It is bounded. No production conclusion is
drawn from it, no default changes, P1-P3 stay blocked, and A6 stays failed.

## Why

D0 split the free troposphere's RH budget (300-700 hPa, days 10-20, d ln RH)
term by term:

| term | d ln RH |
|---|---:|
| subsidence, drying and warming | -7.03 |
| re-evaporation of convective condensate, with its cooling | +5.66 |
| radiation | +1.27 |
| net detrainment, with its heating | +0.74 |
| large-scale condensation, with its latent heat | -0.57 |

**The re-evaporation term splits almost evenly,** already measured in D0
before this protocol was frozen:

- vapour added (d ln q): +2.98;
- cooling (−d ln q_s): +2.68.

**The finding:** the circulation dries the column through subsidence.
Re-evaporation returns most of that drying, in water and in cooling.

## Hypothesis

The excessive free-tropospheric humidity comes from where convective
condensate evaporates, not from the total evaporated fraction alone.

Total evaporated fraction ≠ vertical distribution of evaporation.

## What the current pathway is, read before this was frozen

**E0's mechanics.** A fraction α = 0.67 of the condensate the plume forms in
each layer:

- is offered for evaporation in that layer;
- evaporates up to the layer's saturation deficit;
- passes the remainder to the layer below, under the same limit;
- reaches the surface as precipitation only if no layer can take it.

The other 1 − α leaves as precipitation at once.

**So E0 already limits evaporation by the deficit and passes the rest down.**
E1 as first worded would be E0. What differs is the limit itself, found while
testing D1:

- `deposit_evaporation` uses the deficit before the evaporation's own
  cooling;
- a layer can therefore end above saturation (0.2 percent in the unit
  sounding);
- the water module then condenses the overshoot.

E1 makes the limit consistent with that cooling, and is otherwise E0.

**What α measures.** `EVAPORATING_CONDENSATE_FRACTION` is
1 − PE from Langhans, Yeo & Romps (2015, LES of cumulus congestus). There,
PE is the fraction of condensate reaching the surface, and the rest
evaporates anywhere:

- much of it is cloud condensate, evaporating near where it formed;
- not only rain falling below cloud base.

So E0's placement is close to what the measurement counts. A rain-only
pathway (E2, E3) is expected to realise less than 0.67. Whether that is a
model discrepancy or a change of meaning is what this protocol reports. It is
not a criterion.

**The rate law's magnitude, estimated before the run.**

- **E2's evaporation, from the numbers below:** a 3 mm/h flux through a
  33 hPa layer at RH 0.9 evaporates about 0.03 kg m^-2 per hour,
  1e-6 × 0.1 × √(8.3e-4) × 3340 × 3600.
- **What falls through that layer:** 3 kg m^-2 per hour.
- **So E2 is expected to realise far less than 0.67.**
- **The law and its constant are frozen here,** not adjusted after.

## Setup: identical to D0 (`detrainment-thermodynamic-state.md`)

- **Configuration:** R4 with mass-flux convection, `detrained_state="plume"`,
  from R4's state.
- **Radiation:** the common radiation, replayed and checked value for value.
- **Length and window:** 480 steps at 3600 s; days 10-20.
- **Held fixed:**
  - the closure, the trigger, ε = δ = 1 km⁻¹;
  - the plume and subsidence operators;
  - surface exchange;
  - the initial state and timestep;
  - the conservation tolerances.
- **α_potential = 0.67 in every case.** Γ_M is measured with the case's own
  operator.
- **D2 is not combined with any case.**

## Cases

**The option is `precipitation_evaporation`.** In every case the plume's
condensate per layer is unchanged. Every limit is evaluated on the state the
operator leaves, layer by layer.

**E0, `"cascade"`:** the existing pathway, bit for bit.

**E1, `"wet_bulb_cascade"`:** E0, with each layer's limit replaced by its
wet-bulb capacity. That is the largest δ with
q + δ ≤ q_s(T − L·δ/c_p, p), solved by bisection to 2⁻⁶⁰ of the bracket. The
cooling L·δ/c_p is applied to the layer, as in E0.

**E2, `"falling_rate_law"`:** precipitation falls interface by interface,
following the evaporation of convective precipitation in CAM5's
Zhang-McFarlane scheme (`zm_conv.F90`, subroutine `zm_conv_evap`, a
"Sundqvist type algorithm").

- **Through layer k** falls F_k, all condensate formed above, less what
  evaporated above, over the substep (kg m^-2 s^-1).
- **Rain formed in layer k is not evaporated there.** It joins the flux to
  the layer below, as in CAM.
- **Evaporation in layer k:**
  e_k = min[ ke·(1 − q_k/q_s,k)⁺·√F_k·m_k·h, wet-bulb capacity_k, F_k·h,
  α·(condensate formed above) − (evaporation above) ].
- **ke = 1.0e-6 (kg m^-2 s^-1)^-1/2 s^-1,** the value the WRF port of CAM5's
  scheme sets. CAM itself reads it from the namelist (`zmconv_ke`).
- **CAM's clear-sky factor (1 − cloud fraction) is taken as 1:** this scheme
  carries no convective cloud fraction. That is the largest evaporation the
  law allows.
- **CAM's supersaturation limit (q_s − q)/Δt** is replaced by the wet-bulb
  capacity, CAM's comment calling its own "approx".
- **Whatever reaches the surface is precipitation.**

**E3, `"below_cloud_base"`, a diagnostic limit:**

- E0's rule, allowed only in layers below cloud base;
- cloud base is the lowest layer in which the plume has condensate;
- condensate falls through the layers above without evaporating.

**In all four cases:**

- water and energy close;
- net condensation is formed less evaporated;
- the A1-A5 unit tests pass.

## Diagnostics

**Per layer and as column integrals, days 10-20,** for 300-700 and
100-850 hPa:

- precipitation entering the layer from above;
- condensate formed;
- condensate evaporated;
- precipitation leaving the layer;
- the evaporation limit before evaporation, the step-mean of each substep's
  limit:
  - E0 and E3: the deficit;
  - E1 and E2: the wet-bulb capacity;
- the vapour added by evaporation, and the latent cooling;
- d ln RH from the added vapour, and from the cooling, separately;
- subsidence drying, net detrainment, and large-scale condensation, as in D0.

**Over the column:**

- realised α = (convective evaporation) / (convective condensate formed);
- surface precipitation;
- the mass-weighted mean pressure of the evaporation;
- the fraction of it inside 300-700 hPa.

**And as in D0:** mean RH 300-700, its profile, duty cycle, mass transport,
subsidence span, CAPE, validity, residuals and parameter readback.

## Decision rules, frozen

**What the words mean:**

- **"Material" for RH:** a fall of ≥ 0.03 in the 300-700 hPa mean, from E0.
- **"Material" for the re-evaporation contribution** (its d ln RH, vapour and
  cooling together): a fall of ≥ 10 percent from E0's.
- **"Substantially changes the vertical distribution":** the fraction of
  column evaporation inside 300-700 hPa changes by ≥ 0.10. Alternatively,
  its mass-weighted mean pressure moves by ≥ 50 hPa.

**Placement supported, by E1, E2 or E3 meeting all of:**

1. RH 300-700 falls by ≥ 0.03 from E0;
2. the re-evaporation contribution falls by ≥ 10 percent from E0;
3. water and energy close: the conservation gate passes, with residuals
   ≤ 1e-12;
4. α_potential reads back 0.67;
5. evaporation moves downward (mean pressure increases), or realised α
   falls with the difference reaching the surface;
6. subsidence drying stays active: its d ln RH is at least half E0's in
   magnitude.

**Placement rejected, by a case that:**

- substantially changes the vertical distribution;
- leaves RH within 0.03 of E0;
- leaves the re-evaporation contribution within 10 percent of E0's.

**Bulk α implicated:**

- **The condition:** E1 and E2, the physically local implementations, each
  realise α ≤ 0.50 (a quarter below 0.67) and each lower RH by ≥ 0.03.
- **The conclusion then:** the prescribed bulk evaporating fraction is
  incompatible with the modelled precipitation trajectory and environmental
  saturation deficits.
- **What it would justify:** a separate protocol on α, its observational
  meaning, and whether it represents evaporation, precipitation efficiency,
  or both.

**Reporting.** More than one outcome may hold, and each that holds is
reported. Anything else is reported as found, without classification.

**Not done:**

- no tuning of α or ke;
- no target RH (0.81 is a comparison, not a target);
- no return to the trigger.

## Code, introduced by this protocol only

- **`precipitation_evaporation`** on `mass_flux.convect`, `apply_mass_flux`
  and `cloud_base_mass_flux`, and the convection and configuration option of
  that name:
  - "cascade" (the default, bit for bit);
  - "wet_bulb_cascade";
  - "falling_rate_law";
  - "below_cloud_base".
- **`MassFluxStep` fields:** precipitation entering each layer, and the
  evaporation limit. Both are published as rates.
- **Tests:**
  - A1-A5 in every mode;
  - E0 identical to the code before;
  - E1 never leaves a layer above saturation after its own cooling;
  - E2 evaporates no rain in the layer that formed it and follows the law's
    closed form in an unlimited layer;
  - E3 evaporates nothing at or above cloud base;
  - realised α ≤ α in every mode.

---

## Result (appended; nothing above this line changed)

**Erratum, corrected while running, both readings kept.**

- **What happened:** E1 and E2 first stopped at step 35. Convection raised
  "temperatures must be positive".
- **The defect, mine:** the wet-bulb bisection bracketed on
  [0, q_s(T) − q]. Near 1 hPa, in a warm layer, q_s(T) exceeds anything the
  layer could cool to take, and the trial temperature T − L·δ/c_p went below
  zero.
- **The correction:** the bracket's top is now
  min(q_s(T) − q, c_p·T/(2L)). The excess is checked to be non-negative
  there, and the step refuses if it is not.
- **The test:** `test_the_wet_bulb_capacity_is_bracketed_near_the_model_top`.
- **What it touched:** E0 and E3 never call this function. Their runs, already
  complete, stand.
- **What was read:** E0 and E3 had been read before the fix. E1 and E2 had
  produced no result.

**Setup reproduced.**

- **E0 is D0:** RH 0.93720 and every water amount identical.
- **E0's code is bit-for-bit the code before this change.** The digest of
  `convect` on three soundings and six coupled steps matched the committed
  D-state code.
- **The radiation was identical** at every step (checksum -2646.2636).
- **All four runs closed to 4.9-6.3e-16.**
- **Validity:** all four are invalid on the validity-domain gate only, as
  D0-D2 and T0-T3 were. The conservation gate passed.
- **α_potential reads back 0.67 in all four,** and ke reads back 1e-6 in E2.

### Days 10-20

| | E0 cascade | E1 wet-bulb | E2 falling, rate law | E3 below cloud base |
|---|---:|---:|---:|---:|
| **RH 300-700** | 0.9372 | 0.9699 | **0.7114** | **0.7034** |
| re-evaporation's d ln RH | +5.66 | +5.80 | +0.15 | 0 |
| of which vapour / cooling | +2.98 / +2.68 | +2.87 / +2.93 | +0.09 / +0.06 | 0 / 0 |
| realised α (evaporated / formed) | 0.670 | 0.670 | **0.062** | 0.670 |
| 1 − surface precipitation / formed | 0.64 | 0.67 | 0.02 | 0.28 |
| evaporation's mean pressure, hPa | 550 | 556 | 677 | 996 |
| fraction of it in 300-700 | 0.70 | 0.71 | 0.49 | 0 |
| large-scale condensation, 300-700, kg m^-2 | 3.03 | 0.31 | 0.13 | 0.17 |
| subsidence's d ln RH | -7.03 | -7.27 | -3.88 | -3.72 |
| net detrainment's d ln RH | +0.74 | +0.34 | +2.98 | +3.03 |
| budget sum / actual | 0.037 / 0.038 | 0.052 / 0.053 | 0.49 / -0.05 | 0.54 / -0.07 |
| duty cycle | 0.91 | 0.93 | 0.40 | 0.39 |
| mass transport, kg m^-2 | 7207 | 7359 | 4086 | 3953 |
| mean CAPE_d, J/kg | 382 | 441 | 118 | 114 |
| surface precipitation, kg m^-2 | 43.3 | 42.9 | 41.2 | 29.4 |
| surface temperature, K | 311.1 | 311.0 | 313.8 | 313.7 |

**The budget does not reconstruct in E2 and E3.** It is reported as measured
and not leaned on.

- **The vapour terms leave a residual against the labels:** -0.97 in E2 and
  -1.19 in E3. That residual is the dry convective adjustment's own vapour
  mixing:
  - it shares convection's vapour label;
  - it was zero inside 300-700 hPa in E0, E1 and D0-D2;
  - with convection off 60 percent of the time, it reaches into the band.
- **Even with it included,** the linearised sum misses the actual change by
  about 0.6. With intermittent convection the per-step changes are too large
  for Σ Δw/w to stand in for Δ ln w.

### Verdicts, as frozen

**E2: placement supported, by all six conditions as computed.**

1. RH falls 0.226.
2. The re-evaporation contribution falls 97 percent.
3. Residual 6.3e-16, and the conservation gate passed.
4. α_potential is 0.67.
5. Evaporation's mean pressure rises 127 hPa, and realised α falls to 0.062.
6. Subsidence's d ln RH is 3.88, at least half of E0's 7.03. **This is read
   from a budget that does not reconstruct.**

**E3: placement supported, by all six conditions as computed.**

1. RH falls 0.234.
2. The re-evaporation contribution falls 100 percent.
3. Residual 5.3e-16.
4. α_potential is 0.67.
5. Evaporation's mean pressure rises 446 hPa.
6. Subsidence's d ln RH is 3.72, at least 3.51. **The same caveat applies,
   and the margin is 6 percent.**

**E1: not classified.**

- RH rose 0.033, so it is not supported.
- Its distribution moved by 0.014 and 6 hPa, so it is not rejected.

**Bulk α: not implicated by the frozen rule.** E1 realises 0.67 and is
wetter; the rule needs both E1 and E2 at or below 0.50 and drier.

### What the result says, and what it does not

**1. Where the evaporation goes sets this column's humidity.**

- E2 and E3 both take the evaporation out of the free troposphere, and RH
  falls from 0.94 to 0.70-0.71.
- This is the largest response any protocol in this branch has produced.
  The trigger moved RH by 0.002, and D2 by 0.039.

**2. Neither E2 nor E3 isolates placement from amount.**

- **E2 evaporates far less.** Realised α is 0.06 against 0.67. CAM's rate law
  at its own constant cannot realise the measured precipitation efficiency,
  as estimated before the run.
- **E3 nominally evaporates 0.67,** all into the lowest layer. But that layer
  cannot hold it at the deficit before cooling, and the water module
  re-condenses much of it:
  - surface precipitation is 29.4 against 13.5 for the convective net alone;
  - so the net fraction kept out of the surface is 0.28.
- **Each case therefore changes both where and how much.**

**3. The convective regime changes with them.**

- Under the common radiation, removing the free troposphere's evaporative
  cooling cuts:
  - the duty cycle from 0.91 to 0.40;
  - mass transport by 44 percent;
  - CAPE from 382 to about 115.
- The surface warms by 2.6 K.
- **The drier free troposphere comes with less convection and less
  subsidence, not more.**
- The column is bounded by a prescribed radiative cooling that does not
  respond, and an interactive column could settle elsewhere. That was not
  run.

**4. E1 found what the large-scale condensation was.**

- E0's 3.03 kg m^-2 in 300-700 hPa was almost all the α pathway's overshoot,
  evaporation past saturation then condensed and rained out. The wet-bulb
  limit leaves 0.31.
- **Removing it makes the column wetter by 0.033.** The overshoot was a hidden
  path that rained out water the cascade had just evaporated.

**5. Re-evaporation's effect is both vapour and cooling, in about equal parts.**

- E0: +2.98 vapour, +2.68 cooling.
- E1: +2.87 vapour, +2.93 cooling.

**Not concluded:**

- no production conclusion, no default change, no α or ke changed;
- **why α and the rate law disagree** by a factor of ten in realised
  evaporation is not established here. The candidates are:
  - the measurement's α counts cloud-condensate evaporation near where it
    forms, which a rain-only law omits;
  - or ke = 1e-6 is a tuning value that does not carry that meaning;
- that question belongs to the separate α protocol this one's rules name,
  and the user decides whether to open it.

### Addendum: the user's reading of the E result (2026-10-13)

The frozen verdicts stand. Their interpretation is qualified as the user
worded it:

> E2 and E3 support the proposition that the evaporation treatment
> materially controls free-tropospheric humidity. They do not independently
> establish that vertical placement is the cause, because neither case
> preserved the evaporated condensate amount or the convective state.

**The jump from about 0.94 to 0.71 is not yet a placement effect.**

**E1, as the user records it:**

> E0's large-scale condensation was primarily recycling evaporation that
> had been forced beyond the receiving air's capacity, not condensation
> independently generated by resolved cooling or ascent.

**The candidates:**

| candidate | standing |
|---|---|
| E0 | rejected as a production candidate: it forces an evaporation-recondensation cycle |
| E1 | a diagnostic only |
| E2 | the rain-evaporation component |
| E3 | a limiting placement diagnostic |

**What follows:** `convective-condensate-ledger.md`.
