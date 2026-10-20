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
