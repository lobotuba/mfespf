# Protocol: why the descending boundary layer is weakly coupled to its surface (V0-V5)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The design is the
user's (2026-10-19), mapped onto the operators this model actually has. It was
approved to freeze with the user's four revisions:

1. V2 is a labelled budget, not tracer-equivalent;
2. V3's gross exchange is limited per step;
3. the thresholds are frozen as formulas;
4. the two interface diagnoses keep separate roles.

**It follows from** `regime-cloud-attribution.md` (frozen endpoint 9b35002):

- **Clouds are rejected as the cause.** C1's regime clouds raised ΔF from
  140.5 to 156.5 W m⁻².
- **The imbalance is not radiative.** The two columns' net surface radiation
  is equal: 193.25 against 193.31 W m⁻².
- **Where it is instead.** The imbalance is in the non-radiative surface
  exchange: latent heat (rising minus sinking, 66.5) and convection's
  surface-air exchange (71.5).

**Constraints:**

- **The frozen configuration:** C0's shared clouds, prescribed surface
  temperatures, circulation closure, radiation and every rising-column
  setting.
- **No tuning** against the redistribution or the observed RH. Those are
  evaluated, never targeted.
- **No interactive-surface experiment** is opened from here, unless a
  physically defensible, conservative intervention reduces ΔF without
  warming, drying or forcing the descending surface by construction.

## The question

Why does the descending column suppress both evaporation and
convection-mediated surface heat exchange, despite a moist, capped boundary
layer? Two possibilities, kept explicit throughout:

- its thermodynamic state suppresses the fluxes;
- or it stays in that state because the model lacks a pathway to ventilate
  it through the inversion.

## Read before freezing: the descending boundary layer's operators

| # | operator (module, label) | what it does |
|---|---|---|
| O1 | surface-based convective adjustment (convection: `dry_convective_adjustment`, `surface_convective_adjustment`; vapour: `convective_vapour_transport`) | Relaxes each unstable group toward the entraining-plume profile, with the surface as the bottom element, conserving enthalpy. Mixes vapour to a uniform specific humidity over the same groups. |
| O2 | bulk sensible heat (convection: `surface_sensible_heat`) | H = k_H (T_s − T_1) |
| O3 | bulk evaporation (water: `evaporation_minus_condensation`, lowest layer) | E = k_E (q_sat(T_s) − q_1), into the lowest layer only, with k_E = k_H / c_p |
| O4 | vapour eddy diffusion (water: `vapour_eddy_diffusion`) | K = 0.39 m² s⁻¹, vapour only |
| O5 | the WTG circulation (`circulation_subsidence`) | subsidence, plus lateral exchange where the mass flux diverges. The flux is zero at the surface. |
| O6 | the water module's saturation adjustment | condensation, rained out |
| O7 | radiation | |

**What the model does not have:**

- **no boundary-layer scheme;**
- **no explicit inversion entrainment;**
- **no wind;**
- **"dry adjustment" and "default convection" are one operator here, O1.**

## The two interface diagnoses (output only; separate roles)

**The model's mixed-layer top.** The top of O1's surface-based group each
step, `surface_mixed_layer_top`, a new diagnostic-only convection output.
When O1 forms no surface group, it is the lowest layer.

- **Used for:** attributing O1.

**The observational-style inversion.** The interface of largest RH decrease
between adjacent layers below 600 hPa, diagnosed each step.

- **Used for:** cross-inversion transport.

**Always reported, for each:**

- its pressure;
- whether the two coincide;
- their separation in layers and in pressure;
- the RH and temperature jumps at each;
- the frequency with which the observational inversion changes interface.

**Neither is substituted for the other when they disagree.** The
disagreement may itself be the result: O1 could be mixing to a top that does
not coincide with the strongest moisture cap.

**Fluxes across an interface, exact.** For an operator whose labelled
tendency redistributes, the flux across interface k is the sum of its
tendencies in the layers above k, upward positive.

- **Water:** from the vapour labels.
- **Dry static energy:** Σ C·dT/dt, plus the geopotential change the label
  implies where it moves mass. O1 moves none.
- **Moist static energy:** dry static energy plus L × the vapour flux.

**Each step's flux is taken across that step's diagnosed interface,** before
any averaging.

## V0: reproduce C0 and instrument it

**Reproduce C0 bit for bit:** C0's 60000-step configuration. The new
outputs are additions, and G0 requires every shared printed value of
`c0_60000.json` to be identical. V0 saves its final per-region state (T,
vapour, surface temperature and water store) for V1.

**Recorded per column, per step:**

- LE, H and C_conv;
- **the holding relation in full,** from the surface's labelled
  contributions (net radiation, sensible heat, the convective adjustment and
  evaporative cooling): F_hold = −R_net,sfc + H + C_conv + LE. Its closure
  against the published `fixed_surface_boundary_transfer` is verified every
  output interval, to 1e-9 relative. C_conv is a labelled contribution, not
  a residual;
- boundary-layer means (the lowest layer up to the mixed-layer top): T, q,
  RH and moist static energy;
- both interface diagnoses, with their jumps;
- the cross-inversion water and MSE fluxes by operator: O1, O4 and O5 as
  transport, O6 separately;
- O1's occurrence, its group depths, and the heat and moisture it moves;
- precipitation;
- atmospheric and surface storage.

## V1: operator-off attribution (descending column only, from V0's final state)

**Each branch is diagnostic, not proposed physics:**

| branch | change, descending column only | maps to the user's |
|---|---|---|
| V1a | O1 without the surface (`surface_based=False`). The layers still adjust among themselves; the surface couples only through O2 and O3. | dry adjustment / default convection disabled, as far as the column can run without O1 |
| V1b | O1 moves no vapour (`vapour_transport=False`). This needs `condense_after_convection=False` for that column, a second change, declared. | default convection's moisture redistribution off |
| V1c | O4 off (`eddy_diffusivity=0`) | boundary-layer mixing: the only other vertical mixing of vapour |
| — | inversion entrainment | **not applicable:** no such term exists |

**New code:**

- a per-region settings override on `two_column`;
- per-region state injection for the restart.

**Each branch runs 20000 steps,** read over the last fifth.

## V2: the labelled cross-inversion transport budget

**Not called a passive-tracer test.** Vapour-labelled tendencies are not
strictly equivalent to a passive tracer:

- vapour takes part in condensation;
- O1's thermal adjustment is not necessarily literal mass transport;
- a tracer would distinguish moving air from redistributing a constituent.

**The question V2 answers:** how much water and moist static energy each
existing operator carries across the observational-style inversion.

**What counts as transport:**

- **O1, O4 and O5 are transport contributions.**
- **O1's moisture flux and its MSE tendency are reported separately,** so
  its temperature adjustment is not read automatically as an exchange of air
  mass.
- **O6's phase change is reported separately and does not count as
  ventilation.**

**The ambiguous case.** If O1 carries appreciable energy through the
inversion without the corresponding vapour (an O1 MSE flux of at least 10
percent of the descending LE, with an O1 water flux below the material
threshold), V2's result is **mechanism unresolved**. Ventilation is not
inferred. A real tracer is then the next step, under its own protocol.

## V3: a minimal shallow exchange (only if V2 finds no ventilation)

**One explicit, conservative operator.** A gross bidirectional exchange
between the mixed layer's top layer (layer b, at the model's mixed-layer top)
and the first layer above the observational inversion (layer a).

- **How it works:** each layer gives M·Δt of its air to the other and
  receives the same, carrying dry static energy and vapour at its source
  values. The net mass transfer is zero.
- **The constituent fluxes:** F_q = M (q_b − q_a), and F_h = M (h_b − h_a).
  Each is upward or downward as the gradient sets it.
- **M is a gross bidirectional exchange mass flux, not a net upward mass
  flux.**
- **Conservation:** the column's mass, water and energy are unchanged
  exactly, and the two layers' tendencies are equal and opposite, shown for
  every strength.
- **What it never prescribes:** evaporation, surface flux, RH or cloud.

**The per-step limiter.** M·Δt ≤ 0.1 × min(m_b, m_a). If the master step
violates it, the operator subcycles within the step, conservatively.

**The predeclared strengths, none fitted:**

| M, kg m⁻² s⁻¹ | what it is |
|---|---|
| 0 | identity |
| 0.003 | about a 3 mm s⁻¹ entrainment velocity |
| 0.01 | |
| 0.03 | the upper predeclared sensitivity, about the proposed BOMEX cloud-base scale |

**The literature statement for 0.03 is conditional.** Its numerical
attribution must be confirmed from the cited figure, table or archived
BOMEX statistic (Siebesma et al. 2003, "A large-eddy simulation
intercomparison study of shallow cumulus convection", *J. Atmos. Sci.* 60,
1201-1219) before that statement is relied on. The strength itself is frozen
regardless.

## V4: mixing-depth attribution (only if V1 implicates the depth)

**Only the rule ending O1's surface group changes:**

- the existing rule, the target profile;
- dry static stability;
- virtual potential temperature.

**Bulk Richardson is not applicable:** the model has no wind.

**No entrainment is added.** The code for the alternative rules is frozen in
its own erratum before use.

## V5: the convection redistribution audit (only if V1 implicates O1)

**Layer by layer in the descending column, from V0's labels:**

- O1's vapour removed and deposited;
- its temperature tendency;
- its surface enthalpy term;
- precipitation;
- conservation.

**The question:** does O1 rearrange moisture within the capped layer, or
export it through the inversion?

## Primary responses

- descending LE;
- descending C_conv;
- X = LE + C_conv, and LE + H + C_conv;
- the cross-inversion water and MSE fluxes by operator;
- boundary-layer RH and depth;
- the jumps at both interfaces;
- ΔF;
- water and energy closure.

**The observed RH profile and the observed redistribution are evaluation
quantities only.**

## Thresholds, frozen as formulas

**X = LE + C_conv.** The reference values are V0's exact late means; the
printed thresholds are generated from them, never entered independently.

**G0 = X_rise,V0 − X_sink,V0**

| term | rule |
|---|---|
| **substantial exchange response** | X_sink,test − X_sink,V0 ≥ 0.5 · G0 |
| **substantial redistribution response** | ΔF_test ≤ 0.5 · ΔF_V0 |
| **substantial change** | either of the two above |
| **material inversion transport** | F_q,inv,up ≥ 0.10 · E_surface,sink,V0 |
| **negligible** | below 0.1 of the corresponding substantial amount |

**For material inversion transport:**

- both quantities are late-mean water mass fluxes in kg m⁻² s⁻¹;
- upward is positive;
- O6 is excluded;
- the denominator is V0's sinking-column surface evaporation, not the
  branch's own;
- the ratio is reported even below the threshold.

## Decision logic, frozen

| verdict | condition |
|---|---|
| **Operator artefact** | a V1 branch gives a substantial change without new physics |
| **Mixing-depth artefact** | V4 changes only the depth rule and gives a substantial change |
| **Missing ventilation** | V2 finds no operator with material inversion transport, and V3 gives a substantial change at some declared M |
| **Convection redistribution defect** | V5 shows O1's moisture removed and redeposited below the observational inversion, with negligible export through it |
| **Mechanism unresolved** | V2 is ambiguous, as defined above |
| **No attribution** | no isolated intervention gives a substantial change |
| **Invalid** | conservation, reproduction (G0), the holding-relation closure, or a frozen-state check fails |

## What this cannot separate

- **The bulk exchange coefficients.** k_H, and k_E = k_H / c_p, are shared
  by both columns. A wind-speed difference between the regimes is not
  represented.
- **Everything is under prescribed surface temperatures.**

---
