# Protocol: where the overturning removes and deposits boundary-layer air (L0-L1)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The design is the
user's (2026-10-20), mapped onto the circulation operator this model actually
has. It was approved to freeze with four corrections from the draft
discussion:

1. **The 7 percent of V2 is not the lateral export.** It is the error in
   counting the circulation at one interface as M·(q_above − q_below).
2. **A settled boundary layer closes its budget by construction.** So the
   magnitude of the lateral water export is not a free quantity, and no class
   is judged on it.
3. **V0's stored outputs cannot reconstruct the lateral term,** so an
   instrumented restart is needed.
4. **The intercolumn mass identity holds by construction.** It validates the
   instrumentation; it cannot adjudicate the hypothesis.

**It follows from** `regime-surface-exchange.md` (E1, f90aa9f), whose
evaluation found the near-surface humidity contrast reversed:

- **observed (OAFlux, 2 m):** sinking 0.01616, rising 0.01853 kg kg⁻¹;
- **model (V0, lowest layer):** sinking 0.01703, rising 0.01630.

Missing inversion ventilation, convective redistribution, eddy diffusion,
the wind contrast (as the leading cause) and regime clouds are already
excluded. E2's verdict does not control whether this protocol runs.

## The question

> Does the existing WTG circulation remove moist air from the sinking
> boundary layer at the correct vertical levels and deliver it to the rising
> column at levels that can correct the modeled near-surface humidity
> contrast?

It asks about **placement and carried humidity,** not aggregate export
magnitude.

## The operator, as implemented (`circulation.py`)

- **The mass flux:** f(k), downward positive, at each interface of the
  sinking column, from the spectral WTG balance. Interface 0 is the surface,
  where the modes vanish. The rising column's flux is −(1 − σ)/σ · f(k).
- **Flux form, upwind (`upwind_tendency`).** A vertical flux carries the
  upstream layer's value. A lateral exchange carries the upstream column's.
- **Lateral mass inflow into layer k, per unit of the column's own area:**
  I(k) = f(k) − f(k + 1). Inflow takes the other column's value; outflow
  takes the layer's own.
- **The published per-layer rate** (`{region}:vapour_advection_rate`, and the
  labelled `circulation_subsidence` / `circulation_ascent` tendencies) is the
  vertical and the lateral parts summed. The lateral part is computed inside
  the module and not published.

**Three consequences, fixed by the operator before anything runs:**

1. **Lateral exchange connects layer k only to layer k of the other column.**
   The rising column's import placement is the sinking column's export
   placement, layer for layer, scaled by (1 − σ)/σ. A rising-import placement
   defect can differ from the sinking export placement only through the two
   columns' different boundary-layer definitions.
2. **The exported air carries the sending layer's own q.** The export
   humidity is the boundary-layer profile weighted by where the outflow
   happens.
3. **Lateral exchange relaxes the two columns' q in a layer toward each
   other.** In advective form, an importing layer's q moves toward the
   donor's; an exporting layer's q is unchanged, because the vertical inflow
   replaces the exported mass. So the lateral term can moisten the sinking
   lowest layer only if the sinking column imports at layer 0, which needs
   f(1) < 0.

## Code (output-only, before L0)

`WeakTemperatureGradientCirculation` publishes, per region, in addition to
what it does now:

- `{region}:lateral_mass_inflow`, I(k), in kg m⁻² s⁻¹ per unit of the
  column's area;
- `{region}:lateral_donor_specific_humidity`, the q the upwind operator
  carried in each layer (the other column's where I > 0, its own where
  I ≤ 0);
- `{region}:vapour_lateral_rate` and `{region}:vapour_vertical_rate`, where
  vertical is computed as its own flux sum, not as total minus lateral;
- `{region}:lateral_dry_static_energy_rate` and
  `{region}:vertical_dry_static_energy_rate`, the same split for s.

**Requirements:**

- no number the module already produces changes;
- a test shows the new fields exist, the vertical and lateral parts sum to
  the published rate, and the area-weighted lateral water and energy cancel
  between the columns;
- the full gate passes before L0 runs.

## L0: the instrumented restart

**Configuration.** E2-null's exactly: C0's configuration, k_H = 10 written
explicitly in each region, restarted from `v0_final_state.npz`, 20000 steps,
read over the last fifth (steps 16000–19999).

**Validity, all frozen:**

| check | rule |
|---|---|
| reproduction | every shared printed value of `e2null_20000.json` is identical |
| online split | \|vertical + lateral − published rate\| ≤ 1e-12 · max\|published rate\|, per step, per column |
| offline rebuild | `upwind_tendency` on each step's published flux and start-of-step q reproduces the online lateral rate to ≤ 1e-12 · max\|online lateral\| |
| intercolumn identity | \|σ Σ_k lat_rise(k) + (1 − σ) Σ_k lat_sink(k)\| ≤ 1e-12 · (1 − σ) Σ_k \|lat_sink(k)\|, for water and for s |
| lowest-layer closure | \|storage − Σ terms\| ≤ 1e-9 · E of that column, per step |
| conservation and holding relation | as in V0 |

**No result is inferred from steady-state closure** where the direct flux is
available.

## The boundary-layer definition, frozen

- **The boundary layer of each column, each step:** the layers below the
  observational-style inversion of the ventilation protocol. That is the
  interface of largest RH decrease between adjacent layers with p > 600 hPa,
  so BL = layers 0 … k_inv − 1.
- **Second reading only, never in a verdict:** the same quantities with the
  model's mixed-layer top (layers 0 … top).
- **Every step's quantity is taken over that step's boundary layer,** before
  any averaging.

## L1: the diagnostics (late means of per-step values)

**Notation.** D(k) = max(−I_sink(k), 0) is the sinking column's lateral mass
outflow in layer k, and C(k) = max(I_rise(k), 0) is the rising column's
inflow. m_k is the layer mass.

**A. Where the sinking outflow leaves.**

- f_m,BL = Σ_{k∈BL} D(k) / Σ_k D(k), the fraction of all lateral outflow
  leaving below the inversion;
- the lowest layer's flushing ratio Φ₁ = [D(0) / Σ_{k∈BL} D(k)] / [m₀ / Σ_{k∈BL} m_k],
  its share of the boundary-layer outflow over its share of the
  boundary-layer mass;
- the full D(k) and I(k) profiles in both columns.

**B. What humidity the outflow carries.**

- q_export = Σ_{k∈BL} D(k) q_donor(k) / Σ_{k∈BL} D(k);
- R_q,1 = q_export / q₀,sink and R_q,BL = q_export / q_BL,sink, with q_BL the
  mass-weighted boundary-layer mean.

**C. Where the rising column receives it.**

- f_import,lowest = C(0) / Σ_k C(k);
- f_import,BL = Σ_{k∈BL,rise} C(k) / Σ_k C(k);
- the imported q per layer, as the operator selected it.

**D. The lowest-layer budgets, both columns.** Per step, in kg m⁻² s⁻¹:

- surface evaporation E (`surface_evaporation_flux`);
- condensation in layer 0, E − `evaporation_minus_condensation`(0);
- convective vapour transport (0);
- vapour eddy diffusion (0);
- circulation vertical (0) and lateral (0), flux form;
- storage, (w₀ at step end − w₀ at step start) / Δt;
- the residual.

Reported normalized by that column's E. **Also in advective form,** as the
q change per unit time:

- vertical: the inflowing interface flux × (q_upstream − q₀) / m₀;
- lateral: I(0) × (q_donor(0) − q₀) / m₀ where I(0) > 0, and zero where the
  layer exports.

**E. The direct export fraction.**

- f_export = −Σ_{k∈BL} lat_sink(k) / E_sink, flux form.
- **The pre-registered prediction:** about 1.7, inferred in the draft
  discussion from V0's steady-state boundary-layer balance (E 3.74e-5,
  subsidence import 3.30e-5, eddy diffusion out 0.78e-5, kg m⁻² s⁻¹). It is a
  check on that reconstruction, not an acceptance target.

**F. Mass flux.** f(1), the flux at the top of the lowest layer, and its
sign on each step: the operator's condition for the lateral term to moisten
the sinking lowest layer.

## Classifications, frozen

**Each defect class is reported if its condition holds; more than one may
hold.** "Correct transport structure" and "exchange structurally adequate"
apply only if no defect class holds.

| class | condition (late means) |
|---|---|
| **No attribution** | any L0 validity check fails, or a circulation term is found that the decomposition does not publish |
| **Vertically misplaced sinking outflow** | f_m,BL < 0.5, or Φ₁ < 0.5 |
| **Export-humidity defect** | f_m,BL ≥ 0.5, and R_q,BL < 0.90 or R_q,1 < 0.90 |
| **Rising-import placement defect** | f_import,BL < 0.5 |
| **Wrong-direction lowest-layer exchange** | the advective lateral term moistens the sinking layer 0 or dries the rising layer 0, by at least 0.10 × that column's E / m₀ |
| **Correct transport structure** | none of the above; f_m,BL ≥ 0.5; Φ₁ ≥ 0.5; R_q,BL ≥ 0.90 and R_q,1 ≥ 0.90; f_import,BL ≥ 0.5; and the circulation's total advective q tendency in layer 0 is negative in the sinking column and positive in the rising one |
| **Exchange structurally adequate** | correct transport structure holds while the lowest-layer contrast stays reversed (q₀,sink > q₀,rise) |

**If "exchange structurally adequate" holds,** the diagnosis moves to how
each column sets its lowest-layer humidity locally, with the intercolumn
pathway ruled out:

- subsidence into the lowest layers;
- the surface layer's depth and mass;
- phase changes;
- the coupling of the convective adjustment to the lowest layer;
- the single-layer representation of near-surface air.

**Observations are evaluation only.** The observed contrast
q_sink,obs < q_rise,obs is the external test. It does not place the
boundary layer or define "material".

## What this cannot separate

- **The mass-flux profile is the spectral WTG modes' shape.** Sines between
  the surface and the declared tropopause, so where the outflow leaves is
  set by the truncation and the temperature anomaly's projection. A
  placement defect found here is a property of that closure.
- **Two columns stand for a continuum of regimes.** Lateral exchange between
  them is one number per layer.
- **Everything is under prescribed surface temperatures,** in E2-null's
  configuration.

---

## Erratum before L0 runs (nothing above this line changed)

**L0's smoke test found a coupling defect.** It is recorded and fixed under
`advection-composition.md`. In the tropical two-column configuration, the
circulation and each column's tropical exchange both published
`{region}:vapour_advection_rate` and `{region}:advective_heating_rate`. The
coupler silently kept the exchange's, so the water and convection modules
never composed with the circulation's transport in the step it was made. The
user chose to fix it and rebaseline before L0.

**Three readings change. Each is recorded with both versions.**

1. **The online-split check.**
   - **Reading 1 (above):** "the published per-layer rate
     (`{region}:vapour_advection_rate`, and the labelled
     `circulation_subsidence` / `circulation_ascent` tendencies)".
   - **Reading 2:** in this configuration those were two different fields.
     The check compares against the circulation's own rate: its labelled
     tendency, which after the fix is also
     `{region}:circulation_vapour_advection_rate`. The same holds for s,
     against `{region}:circulation_advective_heating_rate`. In the smoke run,
     vertical + lateral equalled the labelled tendency exactly.
2. **The restart.**
   - **Reading 1:** from `v0_final_state.npz`.
   - **Reading 2:** from A1's final state, on the fixed code.
3. **The reproduction reference.**
   - **Reading 1:** `e2null_20000.json`.
   - **Reading 2:** E2-null ran on the pre-fix code, so the reference is
     **L0-null**. That is the same 20000-step restart from A1's final state,
     with k_H = 10 written in each region, without the split. Every shared
     printed value must be identical.

**Unchanged:** the question, the diagnostics, the boundary-layer definition,
the tolerances and the classifications.

**Seen in the 60-step smoke run** (spin-up only, pre-fix, not a result): the
reproduction against a matching restart was identical in all 612 shared
values; the offline rebuild was exact; the intercolumn identities held to
2e-16; and the lowest-layer budgets closed to 3e-15 of E.

## L0 and L1 result (2026-10-20), on the fixed code, from A1's final state

**Validity: every frozen check passed.**

| check | result | rule |
|---|---|---|
| reproduction against L0-null | 791 of 791 shared values identical | identical |
| online split (against the circulation's own rate) | 0.0 | ≤ 1e-12 |
| offline rebuild, q and s | 0.0 and 0.0 | ≤ 1e-12 |
| intercolumn identity, water and s | 3.1e-16 and 1.3e-16 | ≤ 1e-12 |
| lowest-layer closure, rising and sinking | 4.7e-16 and 3.2e-16 of E | ≤ 1e-9 of E |

**Late-mean lateral mass inflow I(k)** (×10⁻³ kg m⁻² s⁻¹, per unit of each
column's area; negative is outflow), humidity and RH:

| k | p, hPa | I sinking | I rising | q sinking | q rising | RH sinking | RH rising |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0 | 996.4 | −2.103 | +1.665 | 0.01700 | 0.01626 | 0.816 | 0.724 |
| 1 | 962.6 | −1.069 | +0.846 | 0.01662 | 0.01561 | 0.829 | 0.721 |
| 2 | 928.9 | −0.630 | +0.499 | 0.01661 | 0.01560 | 0.864 | 0.750 |
| 3 | 895.2 | −0.485 | +0.384 | 0.01661 | 0.01560 | 0.903 | 0.782 |
| 4 | 861.4 | −0.046 | +0.037 | 0.01661 | 0.01560 | 0.948 | 0.818 |
| 5 | 827.7 | **+0.447** | −0.354 | 0.01653 | 0.01560 | 0.993 | 0.859 |
| 6 | 793.9 | **+0.776** | −0.615 | 0.00924 | 0.01560 | 0.551 | 0.906 |
| 7 | 760.2 | **+0.526** | −0.417 | 0.00613 | 0.01560 | 0.375 | 0.959 |
| 8 | 726.4 | **+0.067** | −0.053 | 0.00384 | 0.01538 | 0.244 | 0.994 |
| 9 | 692.7 | −0.059 | +0.047 | 0.00361 | 0.01384 | 0.245 | 0.943 |

The sinking column's boundary layer is layers 0–5 on every late step
(inversion at interface 6). Its other inflow is at 524–423 hPa and 220–103
hPa, carrying 0.02 E of water.

**L1, inversion reading (mixed-top reading in brackets where it differs):**

| diagnostic | value | frozen line |
|---|---:|---|
| f_m,BL, outflow below the inversion | **0.569** | ≥ 0.5 |
| Φ₁, lowest-layer flushing ratio | **2.91** (49% of the boundary-layer outflow from 17% of its mass) | ≥ 0.5 |
| q_export, kg kg⁻¹ | 0.01680 | |
| R_q,1 | **0.988** | ≥ 0.90 |
| R_q,BL | **1.008** | ≥ 0.90 |
| f_import,lowest | 0.276 | |
| f_import,BL (rising) | **0.568** [0.574] | ≥ 0.5 |
| f(1) < 0 on any step | never | |
| **f_export, direct** | **1.748** | predicted ≈ 1.7 (inferred) |

**Lowest-layer budgets, flux form, normalized by each column's E:**

| term | sinking | rising |
|---|---:|---:|
| surface evaporation | +1.000 | +1.000 |
| condensation | 0.000 | 0.000 |
| convective vapour transport | −0.964 | −1.004 |
| eddy diffusion | −0.015 | −0.015 |
| circulation, vertical | +0.928 | −0.423 |
| circulation, lateral | −0.949 | +0.442 |
| storage | 0.000 | 0.000 |
| residual | 0.000 | 0.000 |

**In advective form (the circulation's effect on q₀):**

- **sinking:** vertical −2.34e-9 s⁻¹ and lateral 0, because the layer only
  exports. Drying, at 0.021 of E/m₀.
- **rising:** vertical 0 and lateral +3.58e-9 s⁻¹. Moistening, at 0.019 of
  E/m₀.

**Classification, as frozen.**

- **No defect class holds.**
  - Not vertically misplaced: f_m,BL 0.569 and Φ₁ 2.91.
  - No export-humidity defect: R_q,1 0.988 and R_q,BL 1.008.
  - No rising-import placement defect: f_import,BL 0.568.
  - No wrong-direction lowest-layer exchange: the sinking lowest layer always
    exports, and the rising lowest layer is moistened.
- **So "correct transport structure" holds.** Its lowest-layer conditions
  hold too: the circulation's advective tendency dries the sinking layer 0
  and moistens the rising one.
- **The lowest-layer contrast stays reversed:** q₀,sink 0.01700 against
  q₀,rise 0.01626.
- **Frozen verdict: exchange structurally adequate.** The intercolumn
  pathway is ruled out as the cause of the reversed contrast, and the
  diagnosis moves to how each column sets its lowest-layer humidity locally.

**What sets each lowest layer, recorded as found.** In both columns the
lowest layer's water is a balance of surface evaporation against convection's
removal (0.96 and 1.00 of E). The circulation's net effect on q₀ is about 2
percent of E/m₀ in each.

**A pathway outside the frozen classes, recorded as found.** The sinking
column imports rising-column air laterally at and just above its inversion,
in layers 5–8 (828–726 hPa), at the rising column's humidity (0.0156,
RH 0.86–0.99 there):

- **That is 24 percent of all the sinking column's lateral mass inflow,**
  and 0.776 E of water: 0.185 E into the boundary layer's top layer, and
  0.591 E into the three layers above the inversion.
- **It is what subsides back into the boundary layer.** The water crossing
  the inversion downward is f(6)·q₆ = 0.953 E. The inflow at 828–726 hPa
  raises the sinking RH there to 0.55 and 0.38, against 0.24 above.
- **Whole-column water:** the sinking column exports 2.07 E laterally (1.94 E
  from its boundary layer), imports 1.07 E (0.78 E of it at 828–726 hPa, 0.02
  E above 300 hPa), and its net lateral loss is 1.000 E, its evaporation, as
  a steady state must.
- **So about 40 percent of the boundary-layer water the sinking column sends
  to the rising column comes back at 830–730 hPa,** and subsides into its
  boundary layer again.

This shallow return flow was named in the user's design as a profile feature
to report ("midlevel return flow"). No frozen class judges it, and nothing
here says whether its strength or its humidity is realistic. That is the
user's to decide.

## Disposition (the user's, 2026-10-21, verbatim)

> The WTG circulation exports representative sinking-column boundary-layer
> air and imports it into the rising boundary layer with the expected sign
> and placement. It directly opposes the modeled near-surface humidity
> error, though weakly in each lowest-layer net budget. The resolved
> lateral-exchange structure is therefore adequate under the frozen criteria
> and is not the cause of the reversed humidity contrast. A shallow return
> branch imports rising-column air at 828–726 hPa amounting to 0.78 of
> sinking-column evaporation; its subsequent contribution to the sinking
> lowest layer remains unclassified.

**What follows:**

- **The lateral exchange operator is not modified.**
- **No broad protocol covering the whole redirect list is opened.**
- **Next is a narrow shallow-return-flow attribution protocol:** spectral
  attribution of the return branch, then tagging of the return water.
- **Then, by the result:** a mode/placement protocol if the return flow is
  controlling, or a local lowest-layer protocol if it is not.
