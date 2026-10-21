# Protocol: does the shallow return flow moisten the sinking boundary layer? (R0-R2)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The design is the
user's (2026-10-21), mapped onto the circulation and the operators this model
actually has. The user approved it to freeze with two decisions:

1. **"substantial" for the supported class is the frozen substantial bar,**
   ΔD_subst;
2. **the donor-humidity reference is NCEP R1 RH on the rising mask,** with a
   +0.10 excess.

**It follows from** `boundary-layer-lateral-exchange.md` (L0, 2ea893b; the
user's disposition, 09d4f73):

- **The lateral exchange is structurally adequate.** It exports
  representative sinking boundary-layer air low down, and the rising
  boundary layer receives it.
- **The circulation's net effect on each lowest layer is about 2 percent of
  E/m₀,** in the corrective direction. In both columns the lowest layer is a
  balance of evaporation against convection's removal.
- **A shallow return branch, recorded as found:** the sinking column imports
  rising-column air at 828–726 hPa (layers 5–8), carrying 0.776 E of water
  (0.185 E into its boundary layer's top layer, 0.591 E above the
  inversion). The water crossing the inversion downward is 0.953 E.
- **The near-surface contrast stays reversed:** q₀,sink 0.01700 against
  q₀,rise 0.01626.

**The lateral exchange operator and the circulation's modes are not
altered.** Everything here is diagnostic.

## The question

> Does the shallow return-flow branch supply enough moisture, at the right
> levels, to cause or materially reinforce the sinking column's overly moist
> boundary layer?

## Code (output-only, before R0)

**Convection publishes this step's mixing groups, when asked**
(`diagnose_mixing_groups`, off by default, plumbed like
`diagnose_lateral_exchange`):

- `convective_mixing_group`: per layer, the index of the group it was mixed
  in, or −1;
- `convective_partial_top_fraction`: per layer, the fraction of it a group
  rose into (zero elsewhere).

These are what `_vapour_transport` and `_mix_classes` use. With them, the
mixing can be replayed exactly offline.

**Requirements:**

- no number the module already produces changes;
- a test shows that replaying the published groups on the published
  provisional vapour reproduces `convective_vapour_transport` exactly;
- the full gate passes before R0 runs.

## R0: reproduce L0, with the groups

**Configuration:** L0's exactly. A1's final state, k_H = 10 in each region,
the lateral split on, 20000 steps, read over the last fifth. The mixing
groups are added.

**Validity, all frozen:**

| check | rule |
|---|---|
| reproduction | every shared printed value of `l0_20000.json` is identical |
| L0's checks | the online split, offline rebuild, intercolumn identity and lowest-layer closure, at L0's tolerances |
| orientation | S1's G5 holds |
| replay closure | the offline operator chain below, applied to the total vapour, reproduces each step's committed vapour in every layer of the sinking column to ≤ 1e-12 · max\|w\| |

**Saved per step over the late window:** everything L0 saved, plus the
groups, the partial tops and the end-of-step vapour.

## The offline operator chain (the model's composition, step by step)

The sinking column's vapour w goes from the start-of-step state to the
committed state as:

1. **circulation:** `upwind_tendency` on the step's published flux, from the
   start-of-step q of both columns, × Δt;
2. **convection:** mix the provisional profile (start + 1) over the step's
   groups and partial tops, exactly as `_vapour_transport` does;
3. **eddy diffusion:** `eddy_diffusion_rate` on the start-of-step vapour,
   × Δt;
4. **surface evaporation:** E·Δt into layer 0;
5. **condensation:** the water module's per-layer removal, from its labels.

The replay closure (R0) shows this chain is the model's.

## R1: spectral attribution of the return branch (offline, from R0)

**The flux by mode.** `spectral_wtg_subsidence` is a sum over modes,
f(k) = Σₙ (aₙ / n) sin(nπζₖ). Its inputs are rebuilt from each step's
saved state:

- the temperature anomaly;
- c_p;
- the rising column's bulk stratification and relaxation time, as
  `_troposphere` forms them;
- the interface pressures.

**The rebuild is checked first:** Σₙ fₙ(k) reproduces the published flux to
≤ 1e-12 · max\|f\|, every step.

**Per mode n, per layer, late means:**

- lateral mass exchange Iₙ(k) = fₙ(k) − fₙ(k+1), in both columns (exactly
  additive);
- vapour and dry-static-energy exchange, Iₙ(k) · donor(k), with the donor
  chosen by the total flow's sign, as the upwind operator chooses it. This
  is the mode's share of the mass, carried at the value the total flow
  carries, and it sums exactly to the total;
- the donor and receiving humidities;
- the mode's ω at each interface;
- the mode's share of the sinking column's import in layers 5–8;
- its subsidence flux across the inversion, fₙ(6) · q₆.

**The question R1 answers:** whether the 828–726 hPa inflow comes from one
dominant mode, several acting together, cancellation between larger
opposing modes, or the reconstruction near the inversion.

**Frozen reading:**

- **one mode dominates** if it carries ≥ 0.5 of the layer-5–8 import mass
  and the others' sum is smaller;
- **cancellation** if Σₙ \|Iₙ\| over layers 5–8 is ≥ 2 × \|Σₙ Iₙ\|.

## R1b: the tagged return water (offline, from R0)

**The tag, t:** the water the sinking column imports laterally in layers 5–8
(a fixed index range, frozen from L0). Its source each step is that import,
max(I(k), 0) · q_donor(k) · Δt, for k = 5–8.

**It is carried through the chain above as the vapour is:**

- **circulation:** upwind on t, with the rising column's t = 0. Tagged water
  leaving the sinking column laterally is an **exit**;
- **convection:** mixed over the same groups;
- **diffusion:** `eddy_diffusion_rate` on t;
- **evaporation:** adds no tag;
- **condensation:** removes t in proportion to its share of the provisional
  vapour in that layer. That is an **exit**.

**The tag affects nothing.** It starts at zero at the window's start and is
read over the window's second half, after its storage has settled.

**The ledger:**

- **exits, mutually exclusive:** lateral export, and condensation. Import =
  lateral export + condensation + Δstorage, to ≤ 1e-9 of the import;
- **transits, reported separately, not exits:**
  - the tagged water crossing the inversion downward, by operator;
  - the gross tagged inflow into layer 0, by operator;
  - diffusion upward across the inversion.

**The key quantities:**

- **f_lowest|return:** the gross tagged inflow into layer 0 over the tagged
  import;
- **R_lowest:** the gross tagged inflow into layer 0 over E_sink;
- **q₀,tag:** the tagged specific humidity in the sinking lowest layer. It
  is the part of q₀,sink that is return-flow water, net of everything the
  outward low-level circulation and convection took away;
- **the no-feedback deficit enlargement** if the return water were absent:
  Δdeficit = q₀,tag.

## R2: counterfactual replays (offline, no feedback, from R0)

**The approximation, declared.** Each replay carries a signed perturbation
tracer δ through the same chain:

- the actual flow, groups, evaporation and condensation are unchanged;
- condensation removes δ at the layer's removal fraction, linearised;
- the result is the no-feedback change in q₀.

The flow is not re-solved. Mass and energy are unchanged because the flow
is unchanged; only the water a source carries in is changed.

**R2a, same mass import, donor humidity replaced.** In layers 5–8 the import
carries:

1. the actual rising-column donor humidity (the baseline; δ = 0);
2. the sinking column's own humidity at the receiving layer;
3. the sinking boundary-layer mean humidity;
4. **not run:** the observed rising-regime near-surface humidity. It is at
   2 m and the import is at 828–726 hPa, so the vertical comparison is not
   defensible.

Each gives δ's source as max(I, 0) · (q_alt − q_donor) · Δt, and its Δq₀ and
Δdeficit.

**R2b, the same imported water entering at a different height.** The tag's
source is redistributed by mass:

- the baseline, layers 5–8;
- entirely above the inversion, layers 6–8;
- entirely within the boundary layer, layers 0–5.

The flow is not changed: this moves where the tagged water enters, not
where the mass enters. **A placement change of the mass itself** would
change the modal flow and cannot be replayed without feedback. It is not
claimed here.

## Thresholds, frozen as formulas

**From A1's baseline,** with k_E = 10 / c_p and L = 2.501e6 J kg⁻¹:

- **the deficit enlargement that would raise LE_sink by A1's substantial
  exchange amount** (63.81 W m⁻²), at fixed state, is
  ΔD_subst = 63.81 / (L · k_E);
  That is 63.81 / (2.501e6 × 0.0099542) = **0.002563 kg kg⁻¹**, about 15
  percent of q₀,sink;
- **"material"** is 0.10 of a frozen reference;
- **"negligible"** is below 0.1 of the corresponding material amount.

## Classifications, frozen

| class | condition (late means) |
|---|---|
| **Invalid** | an R0 check, the modal rebuild, the replay closure or the tag ledger fails |
| **Return-flow moistening supported** | f_lowest\|return ≥ 0.10, R_lowest ≥ 0.10, and q₀,tag ≥ ΔD_subst |
| **Return-flow contributing, not controlling** | f_lowest\|return ≥ 0.10 and R_lowest ≥ 0.10, but q₀,tag < ΔD_subst |
| **Return-flow present but not controlling** | f_lowest\|return < 0.10 or R_lowest < 0.10 |
| **Donor-humidity problem** | R2a case 2 lowers q₀ by ≥ ΔD_subst, and the rising column's RH in layers 5–8 exceeds NCEP R1's rising-regime RH at the same pressures by ≥ 0.10 (evaluation, below) |
| **Modal-placement problem** | R1 finds the layer-5–8 import dominated by one mode or one coherent set; and R2b's above-inversion placement lowers q₀,tag to ≤ 0.5 of the baseline's |
| **No attribution** | q₀,tag below 0.1 × 0.10 × (q₀,sink − q₀,rise), that is below 7.4e-6 kg kg⁻¹ with L0's contrast of 0.00074 |

**The intermediate class "contributing, not controlling"** is added to the
user's list. It covers a material reach that is still below the substantial
deficit change, which the user's classes left unassigned.

**The contrast-sign reading, reported always:** whether
q₀,sink − q₀,tag < q₀,rise. That is, whether removing the return water alone
would undo the reversed contrast, at no feedback.

**Precedence:** invalid, then the three return-flow classes (exactly one
holds). Donor-humidity and modal-placement are reported in addition when
they hold. No attribution overrides all but invalid.

## Evaluation (never fed back, never a target)

**The rising column's humidity at 828–726 hPa against NCEP R1.**

- **The source:** NCEP R1 1991–2020 monthly long-term means of specific
  humidity and air temperature (registered: NCEP-R1-SHUM-LTM-1991-2020 and
  NCEP-R1-AIR-LTM-1991-2020).
- **The method:** RH from them at 850 and 700 hPa, composited on the
  existing ±20° ω500 rising mask, cos-latitude weighted, and interpolated
  in log-pressure to the layers' pressures. RH uses the model's saturation
  function.
- **Provenance:** reanalysis. It is labelled as such and is not an
  observation.

## What this cannot separate

- **No feedback.** R1b and R2 hold the flow, the convection groups,
  evaporation and condensation fixed. A drier lowest layer would evaporate
  more, which the replays do not represent.
- **Condensation is taken in proportion,** and linearised for signed
  perturbations.
- **The modes are the spectral WTG closure's.** A modal-placement finding is
  a property of that truncation and the declared tropopause.
- **Two columns stand for a continuum of regimes.**
- **Everything is under prescribed surface temperatures.**

---

## R0–R2 result (2026-10-21)

**Code:** a483f92, with the full gate passed (685 tests).

**Validity: every frozen check passed.**

| check | result | rule |
|---|---|---|
| reproduction against `l0_20000.json` | 791 of 791 shared values identical | identical |
| L0's checks (split, rebuild, intercolumn, lowest-layer) | 0.0, 0.0, 3.1e-16, 4.7e-16 of E | as L0 |
| orientation (G5) | held | |
| replay closure, sinking column, total vapour | 1.5e-16 · max\|w\| (circulation, convection and diffusion each exact) | ≤ 1e-12 |
| modal rebuild of the flux; the relaxation time | 0.0; 0.0 | ≤ 1e-12 |
| tag ledger | 5.6e-14 of the import | ≤ 1e-9 |
| tag reaching the exchange's layers (above 85 hPa) | at most 1.2e-13 kg m⁻² | |

### R1: where the return branch comes from

**Twenty modes.** Each mode's share of the layer 5–8 import (signed;
negative is export at those layers):

| mode | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| share | −0.706 | −0.002 | **+0.810** | +0.115 | **+0.686** | +0.009 | −0.005 | −0.047 | +0.063 | −0.006 |

**Lateral inflow into the sinking column by mode,** ×10⁻³ kg m⁻² s⁻¹:

| mode | layer 0 | 2 | 4 | 5 | 6 | 7 | 8 | 10 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | −0.476 | −0.457 | −0.414 | −0.383 | −0.348 | −0.308 | −0.263 | −0.165 |
| 3 | −0.608 | −0.398 | −0.002 | +0.209 | +0.394 | +0.533 | +0.608 | +0.538 |
| 5 | −0.572 | −0.071 | +0.515 | +0.596 | +0.482 | +0.210 | −0.131 | −0.586 |

**ω by mode at interfaces 4–9, hPa/day:**

| mode | interface 4 | 9 |
|---|---:|---:|
| 1 | 15.6 | 30.1 |
| 3 | 14.9 | 0.1 |
| 5 | 6.5 | −7.7 |

**The subsidence water flux across the inversion, by mode:** 2.43e-5
(mode 1), 1.43e-5 (mode 3), −3.2e-6 (mode 5) kg m⁻² s⁻¹.

**Frozen readings, both of which hold:**

- **One mode dominates:** mode 3 carries 0.810 of the import, and the other
  modes sum to 0.190.
- **Cancellation:** Σₙ\|shareₙ\| is 2.89 times the net. Modes 3 and 5
  (+1.50) import at 828–726 hPa against mode 1 (−0.71), the gravest mode,
  which exports at every level up to 525 hPa.

So the shallow branch is the third and fifth modes' low-level structure
showing through the deep first mode's outflow.

### R1b: the tagged return water

| quantity | value |
|---|---:|
| tagged import / E_sink | 0.779 |
| into the boundary layer directly (layer 5) | 0.239 of the import |
| **f_BL\|return** (direct + net downward across the inversion) | **0.960** |
| downward across the inversion by the circulation | 0.750 |
| by convection / upward by diffusion | 0.000 / 0.029 |
| **f_lowest\|return** (gross tagged inflow into layer 0) | **0.467** |
| **R_lowest** (that inflow / E_sink) | **0.364** |
| exits: lateral export / condensation / storage change | **0.997** / 4e-6 / 0.003 |
| layer 0, net by operator / E: circulation vertical, lateral | +0.361, −0.361 |
| layer 0, net by operator / E: convection, diffusion, condensation | +0.0002, 0.000, 0.000 |
| **q₀,tag** (second half; first half 0.00592) | **0.00648 kg kg⁻¹** |
| q₀,tag / q₀,sink | 0.381 |

**What the tag does.** Return water enters at 828–726 hPa. Nearly all of it
(0.96) ends up in the sinking boundary layer, and almost all of it (0.997)
leaves again laterally, low down, back to the rising column. In layer 0 its
gross vertical arrival is matched by its lateral export at the same
concentration (±0.361 E), so it passes through. But it is the air the
boundary layer is made of: **38 percent of the sinking lowest layer's vapour
is return-flow water.**

### R2: no-feedback counterfactuals

**R2a, the same mass import with a different humidity:**

| case | Δq₀,sink | Δdeficit |
|---|---:|---:|
| the sinking column's own humidity at the receiving layer | **−0.002364** | +0.002364 (0.922 of ΔD_subst) |
| the sinking boundary-layer mean humidity | +0.000445 | −0.000445 |

**R2b, the same tagged water entering at a different height:**

| entry | q₀,tag | of baseline |
|---|---:|---:|
| layers 5–8 (baseline) | 0.00648 | 1.00 |
| layers 6–8 (above the inversion) | 0.00597 | 0.92 |
| layers 0–5 (within the boundary layer) | 0.00674 | 1.04 |

### Evaluation: the donor humidity against NCEP R1 (reanalysis)

| layer | p, hPa | model rising RH | NCEP rising RH | excess |
|---|---:|---:|---:|---:|
| 5 | 827.7 | 0.859 | 0.650 | +0.210 |
| 6 | 793.9 | 0.906 | 0.616 | +0.290 |
| 7 | 760.2 | 0.959 | 0.582 | +0.378 |
| 8 | 726.4 | 0.994 | 0.545 | +0.449 |
| mean | | | | **+0.332** |

NCEP R1 1991–2020, rising mask (fraction 0.558 on its 2.5° grid): RH 0.671
at 850 hPa and 0.516 at 700 hPa.

### Classification, as frozen

- **Invalid:** no.
- **Return-flow moistening supported:**
  - f_lowest\|return = 0.467 ≥ 0.10 **yes**;
  - R_lowest = 0.364 ≥ 0.10 **yes**;
  - q₀,tag = 0.00648 ≥ ΔD_subst = 0.002563 **yes**.
  - **It holds.**
- **Donor-humidity problem:** the NCEP excess (+0.332) meets its condition,
  but R2a case 2 lowers q₀ by 0.002364, under ΔD_subst (0.922 of it).
  **Not met.**
- **Modal-placement problem:** R1 does find one dominant mode, but the
  above-inversion entry leaves q₀,tag at 0.92 of the baseline, not ≤ 0.5.
  **Not met.**
- **No attribution:** no.
- **Frozen verdict: return-flow moistening supported.**

**The contrast-sign reading (always reported).**

- **Without the return water** (at no feedback), q₀,sink − q₀,tag = 0.01052,
  below q₀,rise = 0.01626. Removing it alone would undo the reversed
  contrast.
- **R2a case 2 would also undo it:** 0.01464 < 0.01626, against a contrast
  of 0.00074.

**Two readings of the same result, recorded as found.**

- **The frozen criterion (q₀,tag)** asks what the lowest layer would hold if
  the return water were absent. The air that sinks into the boundary layer
  is largely return air, so that share is large (38 percent).
- **R2a case 2** asks the narrower question: what if the same air carried
  the sinking column's own humidity instead of the rising column's. It
  changes the deficit by 0.92 of the substantial bar: just short of the
  donor-humidity class, but enough to reverse the contrast's sign.
- **Placement is not what matters.** Moving the entry above the inversion
  changes q₀,tag by 8 percent.
- **The donor air is 0.21–0.45 moister in RH than NCEP's rising regime at
  the same pressures.** That is evaluation, not a result.

**Correction to R1's wording above.** "Mode 1 … exports at every level up
to 525 hPa" overstates what was printed. The mode-1 profile was reported for
layers 0–11 only, and it exports at every one of them, from the surface to
625 hPa. Higher layers were not examined.
