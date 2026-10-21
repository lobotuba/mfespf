# Protocol: what keeps the rising column's donor layers so moist? (D0-D2)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The question and the
candidate terms are the user's (2026-10-21), mapped onto the operators this
model actually has. The user approved it to freeze with one decision: **the
decisive downstream metric is the saturation-consistent removal (D2b)**, with
the proportional composition (D2a) reported beside it.

**It follows from** `shallow-return-flow.md` (R0–R2, dcd4f25; the user's
disposition and prospective deviation, eb8c85f):

- **Return-flow moistening is supported.** The sinking column imports
  rising-column air at 828–726 hPa (layers 5–8). Of that water, 0.960 enters
  the sinking boundary layer, and it makes up q₀,tag = 0.00648 kg kg⁻¹ of the
  lowest layer: 38% of q₀,sink, 2.5 times ΔD_subst = 0.002563 kg kg⁻¹.
- **Donor humidity is strongly implicated but formally not met.** Giving the
  same import the sinking column's own humidity changes the deficit by
  0.002364 kg kg⁻¹, 0.922 of the bar.
- **Placement is weak.** Moving the tagged entry above the inversion changes
  q₀,tag by 8%.
- **The donor air is moister than the reanalysis.** Rising-column RH in
  layers 5–8 is 0.859, 0.906, 0.959 and 0.994, against NCEP R1's rising
  regime at 0.650, 0.616, 0.582 and 0.545 (evaluation only).

**The circulation, its modes 1, 3 and 5, and the lateral exchange are not
altered.** Everything here is diagnostic.

## The question

> Which model operator maintains excessive rising-column humidity at
> 828–726 hPa, and how much of that operator's water is exported through the
> shallow return loop into the sinking boundary layer?

## The donor region and the data

**The donor region, D:** the rising column's layers 5–8. It is the same
fixed index range as the return import, frozen from L0. Water the sinking
column imports at layer k comes from the rising column's layer k (upwind),
so D is exactly where the return water is drawn from.

**No new model run.** Everything is offline, from R0's saved late window
(`r0_late.npz`, 4001 rows, row 0 the step before the window; SHA-256
c1a602fd40344104f8e7600647b45477f1644e1cdd2d014b6c39e332a2d9f9a0). It already
holds, for both columns:

- the start and end-of-step state;
- the circulation's flux, its split rates and its own vapour rate;
- the convective mixing groups and partial tops;
- the water module's labels;
- the surface evaporation.

**No code changes.** The offline analysis uses the model's own
`upwind_tendency` and `eddy_diffusion_rate`, and R1b's chain for the sinking
column unchanged.

## The candidate terms, mapped onto this model

| the user's term | in this configuration |
|---|---|
| convective vapour mixing | convection mixing the provisional vapour to uniform q over the step's groups (`_vapour_transport`) |
| dry adjustment or related layer mixing | **the same operator.** The vapour is mixed over the adjustment's own groups; there is no separate layer mixing. It is split instead by group type: **surface-rooted** groups (containing layer 0) and **elevated** groups |
| vertical circulation | the vertical part of the circulation's upwind transport |
| lateral circulation | its lateral part |
| eddy diffusion | `eddy_diffusion_rate` on the start-of-step vapour |
| condensation and evaporation recycling | **condensation only.** Condensate precipitates immediately and mass-flux convection is off, so nothing re-evaporates. This is checked, not assumed (D0) |
| storage | the change in D's vapour, and the water already in D at the window's start (label I) |

**One more operator, checked to be absent:** the tropical exchange's
stratospheric upwelling, which is nonzero only near 85 hPa.

## D0: the rising column's replay closure

**The chain.** The rising column's vapour w goes from start to committed
exactly as the sinking column's does (`shallow-return-flow.md`):

1. the circulation's upwind transport on the rising flux,
   −(1 − σ)/σ × the sinking flux, from both columns' start q; plus the
   exchange's published rate;
2. convection mixes the provisional profile over the step's groups;
3. eddy diffusion on the start-of-step vapour;
4. surface evaporation into layer 0, then condensation from the labels.

**Validity, all frozen:**

| check | rule |
|---|---|
| input | the SHA-256 of `r0_late.npz` is the value above |
| R1b reproduction | re-running R1b's tag in the same script gives q₀,tag = 0.006481229671725366 exactly |
| rising closure | the chain reproduces the rising column's committed vapour in every layer to ≤ 1e-12 · max\|w\|, every step |
| pieces | the replayed circulation, convection and diffusion reproduce `circulation_vapour_advection_rate`, `label:convective_vapour_transport` and `label:vapour_eddy_diffusion` to ≤ 1e-12 relative |
| exchange absent | the exchange's vapour rate is exactly zero in rising layers 4–9 at every step |
| no re-evaporation | `label:evaporation_minus_condensation` is ≤ 0 in rising layers 1 and above, at every step. If it is not, its positive part becomes an entry operator, **re-evaporation**, and is reported as one |

## D1: the donor region's budget and entry ledger

**The labels, inside D only.** Each label is the water in D that entered it
by one operator, at its latest entry:

- **I:** the water in D at the window's start;
- **V:** entered by the vertical circulation;
- **L:** entered by the lateral circulation;
- **Cs:** entered by convective mixing in a surface-rooted group;
- **Ce:** entered by convective mixing in an elevated group;
- **K:** entered by eddy diffusion.

**Water that leaves D loses its label.** Water re-entering D is relabelled
by the operator that brings it back in. So no labels are needed outside D.

**Entry and exit, exactly by linearity.** Each operator is linear in the
vapour for the step's fixed flux, groups, masses and diffusivity. Split the
vapour into the part outside D and each label inside D:

- **entry by an operator:** that operator applied to the outside water, at
  the point in the chain where it acts, restricted to D;
- **a label's exits:** that operator applied to the label, landing outside
  D;
- **the lateral export from D:** the shallow return's source, per label;
- **condensation:** removes each label in proportion to its share of the
  layer's vapour at that point, exactly as R1b does.

**Checks:**

- Σ labels = the committed vapour in each D layer to ≤ 1e-9 relative, every
  step;
- each label's ledger closes: entries = lateral export + exits by operator
  + condensation + Δstorage, to ≤ 1e-9 of its entries.

**Reported (late means, second half of the window, as R1b):**

- **the D budget by operator,** net per layer and summed over D: vertical
  and lateral circulation, convection (Cs and Ce), diffusion, condensation,
  storage. Each over E_rise;
- **gross entry into D by operator,** and each operator's share of it;
- **each label's share of D's vapour;**
- **each label's exits:** lateral export, other operators, condensation;
- **the condensing frequency:** the fraction of steps each D layer
  condenses, and its mean condensation over E_rise.

**Frozen reading.** An operator **dominates the supply** if it delivers
≥ 0.5 of the gross entry into D.

## D2: downstream, through the measured return pathway

**The return's source, split by label.** At each step the sinking import in
layer k = 5–8 is max(I(k), 0) · q_rise(k) · Δt, as in R1b. Split it by the
donor layer's composition:

source_j(k) = max(I(k), 0) · t_j(k) / m_rise(k) · Δt.

The labels sum to the donor water, so the sources sum to R1b's source
exactly.

### D2a: composition (additive)

Each label's source is carried through R1b's sinking-column chain
unchanged. Its stored concentration in the sinking lowest layer is
**q₀,tag^(j)**, and Σ_j q₀,tag^(j) = q₀,tag. Also reported per label:

- f_BL|return^(j): direct entry into the sinking boundary layer, plus net
  downward across interface 6, by circulation, convection and diffusion;
- its exits: lateral export back to the rising column, condensation,
  storage.

**What it answers:** whose water makes up the return water in the sinking
lowest layer.

### D2b: removal, saturation-consistent (no feedback)

**What it answers:** how much drier the sinking lowest layer would be if
operator j had not supplied D, at no feedback, with D's condensation
responding as a saturation adjustment does at fixed temperature.

**The counterfactual tracer δ^(j), for each operator and for all at once.**

- Its source in D is −(entry by operator j), at the same point in the
  chain, every step.
- It is carried through the rising column's chain with the step's flux,
  groups and diffusivity unchanged.
- **Condensation responds first.** A layer that condenses c in the baseline
  ends at saturation. Removing water there first reduces its condensation:
  δ_cond = max(δ, −c), and the rest remains as a humidity change. A layer
  that does not condense keeps all of δ.
- The rising δ in layers 5–8 is exported through the measured return
  pathway: max(I(k), 0) · δ(k) / m_rise(k) · Δt into the sinking column. It
  is carried through R1b's chain, with condensation linearised as in R2.
- **It is single-pass.** δ reaching the sinking column does not return to
  the rising column. The loop's own feedback is not represented.

**The result:** ΔD^(j) = −Δq₀,sink, the no-feedback enlargement of the
sinking surface deficit.

**Replacement is not run.** There is no defensible alternative operator to
substitute. Substituting an observed profile would make the reference a
target.

**Proportional removal is the same as composition.** Removing operator j's
entry with proportional condensation gives exactly −q₀,tag^(j), by
linearity. So D2a is the proportional bound and D2b the
saturation-consistent one.

## Thresholds, frozen

| term | value |
|---|---|
| ΔD_subst | 0.002563 kg kg⁻¹, frozen in `shallow-return-flow.md` |
| material | 0.10 × ΔD_subst = 0.000256 kg kg⁻¹ |
| negligible | below 0.1 × material = 0.0000256 kg kg⁻¹ |

## Classifications, frozen

The decisive metric, per operator j, is **ΔD^(j) from D2b**.

| class | condition (late means) |
|---|---|
| **Invalid** | a D0 check or a D1 ledger fails; or label I still contributes q₀,tag^(I) ≥ material in the second half (the window is too short to attribute) |
| **Operator implicated** | ΔD^(j) ≥ ΔD_subst. Per the user's prospective definition, removal produces a substantial response in the frozen improvement direction |
| **Operator contributing** | material ≤ ΔD^(j) < ΔD_subst |
| **Operator minor / negligible** | below material / below negligible |

**Verdict:**

- **Single operator implicated:** exactly one operator is implicated;
- **Several operators implicated:** more than one is;
- **Saturation-limited:** no operator is implicated, removing all of them
  together leaves ΔD^(all) < ΔD_subst, and Σ_j q₀,tag^(j) ≥ ΔD_subst
  (composition). Then the donor water is mostly supply that D would
  condense anyway, and the humidity is held by D's saturation, not by any
  one supply;
- **Shared supply:** no operator is implicated, but ΔD^(all) ≥ ΔD_subst.

**Reported always:**

- the composition q₀,tag^(j) beside each ΔD^(j);
- whether removing each alone would undo the reversed contrast,
  q₀,sink − ΔD^(j) < q₀,rise, at no feedback;
- the supply-dominance reading from D1.

**Precedence:** invalid first, then exactly one verdict.

## Evaluation (never fed back, never a target)

**The comparison is the one already recorded:** the rising column's RH in
layers 5–8 against NCEP R1's rising regime (`shallow-return-flow.md`).
Nothing in D0–D2 uses it.

## What this cannot separate

- **No feedback.** The flux, groups and temperatures are the baseline's.
  Removing a supply would also change the circulation, the convection and
  the latent heating.
- **Entry attribution is not operator removal.** Water an operator delivers
  into D might have arrived by another operator if the first were absent.
  D2b removes the delivery, not the operator.
- **The saturation response is at fixed temperature.** It ignores the
  change in latent heating that less condensation would bring.
- **Single pass.** The return loop's recirculation back into the rising
  column is not fed back.
- **Two columns stand for a continuum of regimes,** under prescribed surface
  temperatures.
