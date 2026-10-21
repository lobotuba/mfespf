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

## Results (2026-10-21)

The run is M0, `scratchpad/m0.py`: R0's runner plus the temperature
budget's labels, θ_post, `convective_temperature_tendency` and
`surface_convective_heat_flux`. It ran 20000 steps in 6923 s. The
analysis is `m2_ncep.py`, `m2_isccp.py` and `m1.py`, written to `m1.json`.

**An operational note.** The first M0 launch did not save
`surface_convective_heat_flux`, which the frozen enthalpy check needs
whenever the surface is in a group. It was stopped after a few minutes and
relaunched with the field added. No result was read from the stopped run.

### M0: validity

| check | result |
|---|---|
| reproduction | **59/59 arrays of `r0_late.npz` bit-identical; 791/791 shared printed values identical** |
| temperature closure | 1.9e-16 in both columns |
| **enthalpy, surface-rooted groups** | **2.8e-12 rising, 4.0e-12 sinking, against the frozen ≤ 1e-12. Fails** |
| enthalpy, elevated groups (reported, not frozen) | 2.9e-9 rising, 4.1e-10 sinking |

**The enthalpy residual is floating-point round-off, measured.** Over a
surface-rooted group, Σ C (T_post − T_pre) closes against
F_s · Δt to at most **3.5 ulp of the group's enthalpy, Σ C T** (2.9 ulp
sinking). The frozen check divides by the adjustment, Σ C \|ΔT\|, which is
far smaller than Σ C T. The residual comes from the T → θ → T round trip,
which carries about 1e-16 of T, against adjustments of order 0.1 K.

**Reported only:** committed T − T_post − (latent + sensible) · Δt reaches
0.023 K. So some other tendencies, besides latent heat and convection's own
sensible heat, enter after convection's provisional profile. This does not
affect T_pre, which is T_post − adjustment · Δt by definition.

### M1: rising column (late means, all 4000 steps)

| layer | 0 | 1 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| p, hPa | 996 | 963 | 861 | 828 | 794 | 760 | 726 | 693 |
| q_pre, g kg⁻¹ | 16.27 | 15.63 | 15.60 | 15.60 | 15.60 | 15.60 | 15.38 | 13.86 |
| q_post, g kg⁻¹ | 15.60 | 15.60 | 15.60 | 15.60 | 15.60 | 15.60 | 15.42 | 13.90 |
| RH_pre | 0.727 | 0.726 | 0.823 | 0.869 | 0.912 | 0.965 | 1.001 | 0.951 |
| RH_post | 0.695 | 0.721 | 0.818 | 0.859 | 0.906 | 0.959 | 1.003 | 0.958 |
| saturation deficit post, kg m⁻² | 2.36 | 2.08 | 1.19 | 0.88 | 0.56 | 0.23 | −0.02 | 0.21 |

**The pre-mixing profile is already uniform** from 963 to 760 hPa, because
the previous step mixed it (a declared limitation). Each step's mixing
mainly carries the surface layer's fresh evaporation through the group:
Δq is −0.67 g kg⁻¹ in layer 0 and ≤ 0.03 g kg⁻¹ elsewhere.

**Gradient destruction, per step:**

| | mean | p10 | p50 | p90 |
|---|---:|---:|---:|---:|
| G_q, kg kg⁻¹ | 0.00123 | 0.00063 | 0.00070 | 0.00339 |
| G_T, K | 0.277 | 0.202 | 0.212 | 0.344 |
| G_h, J kg⁻¹ | 3244 | 1773 | 1941 | 8588 |

**The MSE change in the donor layers,** h_post − h_pre: +186, +117, +113
and +143 J kg⁻¹ in layers 5–8.

**The group's geometry:**

| quantity | p10 | p50 | p90 | mean |
|---|---:|---:|---:|---:|
| top pressure, hPa | 675.8 | 743.3 | 743.3 | **716.6** |
| depth, hPa | 270 | 270 | 337 | 297 |

The group contains layer 7 on 100% of steps, and layer 8 on 48.3%. There is
a surface group on every step.

**Vapour carried into the donor layers by mixing, over E_rise:**

| layer | 5 | 6 | 7 | 8 |
|---|---:|---:|---:|---:|
| gross | 14.1 | 14.1 | 14.1 | 6.41 |
| net | 0.0002 | 0.0003 | 0.0034 | 0.072 |

**The condensation the mixing causes, at fixed T_post:**

- 0.073 E in layer 8, and 0.0002 E in layer 7;
- none in layers 0–6;
- 2.06 E over the whole column, about 1.6 E of it above 625 hPa, in
  elevated groups.

**The sinking column, reported:** its surface group's top is fixed at
810.8 hPa (depth 202 hPa), with G_q = 0.00076.

**The surface-rooted return water by mode** (second half, sums to 1):

| mode | 1 | 3 | 4 | 5 |
|---|---:|---:|---:|---:|
| share | −0.625 | 0.823 | 0.090 | 0.638 |

**Modes 3 and 5 together: 1.462,** against mode 1's −0.625.

### M2: structural comparisons

| quantity | model, rising | reference |
|---|---:|---:|
| q at 850 / 700 hPa, g kg⁻¹ | 15.60 / 14.18 | 9.89 / 5.51 (NCEP R1) |
| dq/dp over 850–700, kg kg⁻¹ hPa⁻¹ | 9.47e-6 | 2.92e-5 (NCEP R1); **ratio 0.324** |
| dq/dp over 1000–850 | not formed (the lowest layer is at 996 hPa) | 5.11e-5 |
| LCL, hPa | 921.8 (own lowest layer) | **948.0** (NCEP, per cell and month; 949.1 from the composite mean) |
| surface group top, hPa | **716.6** | ISCCP-H liquid low-cloud top **809.5** (amount 0.152) |
| return-flow-weighted q, g kg⁻¹ | 15.60 | — |

**Measured differences:**

- the group's top is **231 hPa above the NCEP LCL,** 6.9 layers of
  33.7 hPa, against the two-layer line;
- it is **93 hPa above ISCCP's low-cloud top;**
- it is 205 hPa above the model's own LCL.

### Classification, as frozen

- **Frozen verdict: Invalid,** on the enthalpy check (2.8e-12 and 4.0e-12
  against 1e-12).

**Recorded for the user's disposition, and not a verdict.** Were the
enthalpy check held to round-off of the group enthalpy (it is at 3.5 ulp of
Σ C T), the frozen classes would read:

- **top above the NCEP LCL by more than two layers: yes,** 231 hPa;
- **|dq/dp| over 850–700 below 0.25 × NCEP's: no,** at 0.324;
- **so, "Consistent with the observed structure",** by the frozen "or".
  The cloud-layer class is not reached;
- **reported always:** the group extends above the model's own LCL by more
  than two layers: yes.

**Recorded as found, about the gradient window.** In the rising column, q
is uniform from 963 to 760 hPa, and the group's top is at 743 hPa (median).
So the model's whole 850–700 gradient comes from above the group's top,
between about 726 and 700 hPa. Below the top it is zero. The 850–700 window
was frozen before this was known.
