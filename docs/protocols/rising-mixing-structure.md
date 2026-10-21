# Protocol: what does the rising column's surface-rooted mixing homogenise? (M0-M2)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The question, the
diagnostics and the comparisons are the user's (2026-10-21), mapped onto what
this model publishes. The user approved it to freeze with three decisions:

1. **the mixed-layer lines are two model layers and 0.25;**
2. **ISCCP-H low-cloud tops stay as reference (c),** with their dependence
   on the configuration declared;
3. **M0 runs.**

**It follows from** `rising-donor-humidity.md` (D0–D2, ed1ce85; erratum 1 and
its results, d51b132 and 5e2de33):

- **Surface-rooted convective mixing owns the donor water.** It supplies 97%
  of the entry into the rising column's layers 5–8 and 98% of their vapour,
  and 0.00639 of the 0.00648 kg kg⁻¹ return-water tag in the sinking lowest
  layer.
- **The rising column's surface group always spans layers 0–7** (to
  760 hPa), and reaches layer 8 (726 hPa) on 49% of steps.
- **Within it q is uniform,** 15.588 g kg⁻¹ in layers 1–7, so RH climbs from
  0.72 to 0.99 only because the air cools.

> Should the default convection operator homogenise specific humidity from
> the surface through a group about 240 hPa deep, reaching 760 hPa every
> step and 726 hPa on half of them?

**This is a mixing-structure diagnostic, not a new scheme.** Nothing is
changed:

- not convection;
- not the mixing groups or their depth;
- not surface coupling, the circulation or condensation.

**Shortening the group to match the NCEP humidity would tune the operator to
the outcome, and is not done.** A later frozen variant may compare physically
declared structures:

- the current full homogenisation;
- mixing limited to the diagnosed boundary layer;
- entraining dilution above cloud base;
- a conservative vertical exchange that relaxes gradients without
  homogenising.

Any such variant is selected by a structural observable, never by the
sinking column's humidity.

## M0: an instrumented rerun of R0

**Why a run is needed.** R0's saved window has the vapour before and after
mixing, but not the temperature convection adjusted. So it has no
pre-mixing T, RH, MSE or saturation deficit.

**Configuration:** R0's exactly. A1's final state, k_H = 10 in each region,
the lateral split and the mixing groups on, 20000 steps, with the late fifth
saved.

**Added to the saved window, both columns:**

- the temperature state's labelled tendencies, by label (as R0 saved the
  vapour's): `dry_convective_adjustment`, `surface_sensible_heat`, and every
  other label on that state;
- `convective_temperature_tendency`;
- `atmosphere_potential_temperature`, convection's post-adjustment θ.

**Validity, all frozen:**

| check | rule |
|---|---|
| reproduction | every array in `r0_late.npz` (SHA-256 c1a602fd…f9a0) is bit-identical in M0's window, and every shared printed value of `r0_20000.json` is identical |
| temperature closure | start T + Σ labels · Δt reproduces each step's committed T in every layer of both columns, to ≤ 1e-12 · max\|T\| |
| pre/post | T_post = θ_post · (p / p₀)^κ, as convection forms it; T_pre = T_post − (dry_convective_adjustment) · Δt |
| enthalpy | over each surface-rooted group, Σ C (T_post − T_pre), plus the surface's convective adjustment if it is in the group, is zero to ≤ 1e-12 of the gross Σ C\|ΔT\| |

## M1: what the homogenisation destroys (rising column; sinking reported)

**Pre and post.**

- **Vapour:** pre is the provisional profile convection mixes, start +
  (circulation + exchange) · Δt; post is the mixed profile. Both are as
  replayed and validated in D0.
- **Temperature:** as M0.
- **Derived at each stage:**
  - RH, from that stage's T and the model's saturation function;
  - MSE, h = c_p T + g z + L q, with z hydrostatic from that stage's T;
  - **total water = vapour.** Condensate precipitates immediately and none
    is stored;
  - the saturation deficit, (q_sat(T) − q) · m.

**Per step, late means and distributions:**

- **Δq_mix(p) = q_post − q_pre,** per layer.
- **The humidity-gradient destruction,**
  G_q = Σ over the group's internal interfaces of
  \|(∂q/∂p)_pre − (∂q/∂p)_post\| · Δp. Here ∂q/∂p at interface k is
  (q_k − q_(k−1)) / (p_k − p_(k−1)) between layer centres, and Δp is the
  interface's centre-to-centre pressure difference. The same is reported
  for T and h.
- **The group's geometry:**
  - **its top pressure,** at the top interface of its highest whole layer,
    raised by the partial-top fraction times the next layer's thickness;
  - **its pressure depth,** surface pressure minus top pressure;
  - **the fraction of steps it contains layer 7 (760 hPa) and layer 8
    (726 hPa).**
- **Vapour carried into each donor layer by mixing:** D1's surface-rooted
  gross entry and net change, per layer, in kg m⁻² s⁻¹ and over E_rise.
- **The condensation the mixing causes,** a fixed-temperature estimate:
  Σ max(w − q_sat(T_post) m, 0), post minus pre, per layer.
- **The share of the surface-rooted return water carried by modes 3 and 5.**
  It is R1's per-mode lateral mass in layers 5–8, carried at D1's Cs
  composition of the donor layer, over the Cs return source. It is exactly
  additive over modes, as in R1.
- **The MSE change in the donor layers,** h_post − h_pre, layers 5–8.

## M2: structural comparisons (evaluation, never fed back, never a target)

The question is whether **the model mixes surface composition through
layers above the observed mixed or cloud layer.** Each observable is
composited on the existing ±20° annual-mean ω500 rising mask,
cos-latitude weighted. Checksums are verified before parsing.

| observable | source | provenance | model counterpart |
|---|---|---|---|
| **(a) the vertical q gradient, 850–700 hPa** (also 1000–850) | NCEP R1 1991–2020 shum LTM (registered) | reanalysis | rising column's late-mean q, log-p interpolated to 850 and 700 hPa |
| **(b) the LCL,** a parcel from 1000 hPa lifted dry-adiabatically with q conserved, saturating by the model's saturation function | NCEP R1 air and shum LTM at 1000 hPa (registered) | reanalysis | the same calculation from the rising column's lowest layer |
| **(c) the low-cloud top,** amount-weighted mean cloud-top pressure over ISCCP-H's liquid low types | ISCCP-H HGM 1991–2020 `pc_types`, `cldamt_types` (registered) | observational, satellite | the surface group's top pressure |
| **(d) return-flow-weighted q** | model only | — | a downstream consequence, reported, never a selection target |

**Declared limitations of the references:**

- **(c) is not independent of the configuration.** ISCCP's tropical cloud
  types are this run's prescribed cloud field. It was not fitted to the
  mixing, but it is an input.
- ISCCP's low types are defined by cloud tops below the 680 hPa level, so
  their mean top cannot lie above 680 hPa.
- Over the rising regime, low clouds under high cloud are not seen by
  ISCCP.
- (b) uses 1000 hPa rather than 2 m.

## Classifications, frozen

| class | condition (late means) |
|---|---|
| **Invalid** | an M0 check fails |
| **Homogenisation above the observed mixed layer** | the group's top pressure is above the NCEP-derived LCL by more than two model layers' thickness there; and the model's \|dq/dp\| over 850–700 hPa is below 0.25 × NCEP's |
| **…and above the observed cloud layer** | additionally, the group's top is above ISCCP's low-cloud-top pressure |
| **Consistent with the observed structure** | the group's top is within two layers of the NCEP LCL, or the model's \|dq/dp\| is ≥ 0.25 × NCEP's |

**Reported always:**

- **whether the group extends above the model's own LCL** by more than two
  layers. That reading needs no external data;
- G_q;
- the mixing-caused condensation;
- the modes 3 and 5 share.

**Precedence:** invalid, then exactly one of the remaining classes. The
cloud-layer class is reported only when the mixed-layer class holds.

## What this cannot separate

- **It diagnoses structure, not cause.** A group deeper than the observed
  mixed layer may follow from the entraining-plume target's stability
  criterion, or from the radiation and circulation that set the profile.
  M0–M2 do not change either.
- **The pre-mixing profile is the model's own.** It already reflects every
  step's mixing before this one.
- **The fixed-temperature condensation estimate** ignores the latent heating
  it would change.
- **Two columns stand for a continuum of regimes,** under prescribed surface
  temperatures.
- **The references are annual-mean climatologies** on a regime mask, with
  the declared limitations above.

---

## Erratum 1: the ISCCP-H period (before M2 reads any data)

**Recorded before any M2 parse; both readings kept.**

- **The frozen text** gives reference (c) as "ISCCP-H HGM 1991–2020".
- **The registered record** runs from **1991-01 to 2017-06** (318 monthly
  files, `ISCCP-H-HGM-BASIC-YYYY-MM`). That is where the ISCCP-H HGM product
  ends, as `regime-cloud-attribution.md` records.

A run cannot read past a record's coverage. So (c) uses 1991-01 to 2017-06,
the annual mean of its monthly climatology, exactly as the regime-cloud
composite formed it. The rising mask is unchanged: NCEP R1 1991–2020
ω500, bilinearly interpolated to ISCCP's 1° grid.

**The liquid low types** are cumulus_liquid, stratocumulus_liquid and
stratus_liquid.
