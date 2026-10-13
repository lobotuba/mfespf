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
