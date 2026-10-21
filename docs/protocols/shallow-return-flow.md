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
