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

## Erratum before V1's code is used (nothing above this line changed)

**Reading 1 (above).** V1b, O1 moving no vapour, "needs
`condense_after_convection=False` for that column, a second change, declared."

**Reading 2, from building it.** The configuration refuses each step short
of three changes:

- `vapour_transport=False` alone is refused: condensing after convection
  composes with convection's vapour transport;
- adding `condense_after_convection=False` is refused: the end-of-step
  saturation needs the water module to run after convection;
- adding `end_of_step_saturation=False` as well builds and runs.

**So V1b changes three settings in the descending column,** all declared:

- `vapour_transport=False`;
- `condense_after_convection=False`;
- `end_of_step_saturation=False`.

**What this means for reading V1b.** The end-of-step saturation is the
water module's removal of supersaturation left after convection, and it is
switched off with the other two. A V1b result is therefore attributed to
"O1's vapour transport off, with the condensation ordering that requires",
not to O1's vapour transport alone.

V1a and V1c are unchanged.

## Erratum 2 before V1 runs: V1a (nothing above this line changed)

**Reading 1 (above).** V1a is "O1 without the surface
(`surface_based=False`). The layers still adjust among themselves; the
surface couples only through O2 and O3", run "as far as the column can run
without O1".

**Reading 2, from building it.** The configuration refuses
`surface_based=False` on its own. O1's target is the entraining plume, a
lapse rate measured upward from the ground, and `convection.py` refuses
that target without a surface instead of anchoring it to the lowest
layer. Turning convection off entirely is not available either: the moist
column is defined with convection. Building it used to overwrite a
caller's `convection=False` silently; it now refuses it (committed
separately, with a test).

**Decision (the user's, 2026-10-19): the dry-adiabat variant.** V1a
changes two settings in the descending column, both declared:

- `surface_based_convection=False`;
- `convective_target="dry_adiabat"`.

**What this means for reading V1a.** A V1a result is attributed to "O1's
surface coupling and its moist target together", not to the surface
coupling alone. The confound is expected to sit mainly in the boundary
layer, because the weak-temperature-gradient balance largely sets the
descending free troposphere from the ascending column. That is an
expectation, not a result; V1a reports its free-tropospheric temperatures
against V0's so the reading can be checked.

**Analysis detail.** A column that is not surface-based publishes no
mixed-layer top. The analysis takes the missing value as the lowest layer,
which is the definition above ("When O1 forms no surface group, it is the
lowest layer").

V1b and V1c are unchanged.

## V0 result (2026-10-19)

**Validity.**

- **Reproduction (G0): passed.** All 80 numeric values that `v0_60000.json`
  shares with `c0_60000.json` are identical, bit for bit. The run completed
  60000 steps; conservation residual 1.0e-15.
- **Holding-relation closure: passed.** The largest relative error over all
  60000 steps is 1.6e-15 (rising) and 9.6e-15 (sinking), against the frozen
  1e-9. The labelled LE and H agree with the published fluxes to 6e-14 W m⁻².
- **A deviation from the frozen wording, recorded.** The relation as frozen,
  F_hold = −R_net + H + C_conv + LE, omits one labelled surface term: the
  surface also pays for melting frozen precipitation
  (`melting_of_frozen_precipitation`). The closure above includes it, and it
  is not small in the rising column (18.0 W m⁻²); in the sinking column it is
  4e-4 W m⁻².
- **A second deviation.** O5's moist-static-energy flux is reported as the
  exchange M·(h_above − h_below) at the interface, because the circulation
  also moves air sideways and a one-sided M·h depends on the reference of h.
  Its water is reported both as the downward flux −M·q_above and as the net
  exchange M·(q_above − q_below).

**Late means (last fifth, 12000 steps).**

| | rising | sinking |
|---|---:|---:|
| LE, W m⁻² | 159.93 | 93.46 |
| H | 1.85 | −0.69 |
| C_conv | 74.72 | 19.17 |
| melt | 18.01 | 0.0004 |
| R_net,sfc | 193.20 | 193.31 |
| F_hold | 61.31 | −81.37 |
| **X = LE + C_conv** | **234.65** | **112.62** |
| LE + H + C_conv | 236.50 | 111.93 |
| surface evaporation, kg m⁻² s⁻¹ | 6.39e-5 | 3.74e-5 |
| precipitation, kg m⁻² s⁻¹ | 9.33e-5 | 1.2e-7 |

**The thresholds, generated from these numbers:**

| quantity | value |
|---|---:|
| G0 = X_rise − X_sink | 122.03 W m⁻² |
| substantial exchange response: X_sink,test − 112.62 ≥ | 61.01 W m⁻² |
| ΔF_V0 | 142.68 W m⁻² |
| substantial redistribution response: ΔF_test ≤ | 71.34 W m⁻² |
| material inversion transport: F_q,inv,up ≥ | 3.74e-6 kg m⁻² s⁻¹ |
| negligible | below a tenth of the corresponding amount |

**The two interface diagnoses.**

| | rising | sinking |
|---|---:|---:|
| mixed-layer top, hPa (mean layer) | 732.7 (7.8) | 827.7 (5.0) |
| steps with a surface group | 1.00 | 1.00 |
| observational inversion, hPa | 701.5 | 810.8 |
| inversion's interfaces in use | 1, 9, 10, 11 | 6 only |
| inversion changes interface, fraction of steps | 0.34 | 0.00 |
| coincide (inversion = top + 1) | 0.12 | **1.00** |
| separation, layers / hPa | 0.42 / 14.3 | 0 / 0 |
| RH jump at the inversion | 0.10 | **0.45** |
| T jump at the inversion, K | 1.79 | 0.73 |

**In the sinking column the two diagnoses never disagree.** O1's surface
group ends at layer 5 on every late step, and the strongest moisture cap is
the interface directly above it on every late step. The jumps at the mixed-
layer top are therefore the jumps at the inversion. In the rising column the
two coincide on 12 percent of steps; the jumps at its mixed-layer top were
not computed separately, a gap in the V0 instrumentation.

**Boundary-layer means (lowest layer to the mixed-layer top).**

| | rising | sinking |
|---|---:|---:|
| T, K | 294.85 | 295.64 |
| q, kg kg⁻¹ | 0.01567 | **0.01665** |
| RH | 0.830 | **0.895** |
| MSE, kJ kg⁻¹ | 349.36 | 347.81 |

**The sinking boundary layer is moister than the rising one,** in specific
and relative humidity, under a surface held 1.47 K colder.

## V2 result: the labelled cross-inversion budget (from V0)

**Sinking column, across each step's observational inversion, upward
positive, late means:**

| operator | water flux, kg m⁻² s⁻¹ | ratio to E_sink,V0 | MSE flux, W m⁻² |
|---|---:|---:|---:|
| O1, convective adjustment | 7.9e-20 | 2e-15 | 2.6e-12 |
| O4, vapour eddy diffusion | 7.78e-6 | **0.208** | 19.45 |
| O5, circulation, downward flux −M·q_above | −3.30e-5 | −0.884 | — |
| O5, net exchange M·(q_above − q_below) | −2.70e-5 | −0.722 | −57.36 (exchange) |
| O6, phase change above (not ventilation) | −1.2e-9 | | |

**Verdict under the frozen rules:**

- **O4 carries material upward inversion transport** (0.208 ≥ 0.10).
- **O1 carries none.** Its water flux is zero to rounding, and so is its MSE
  flux, so the ambiguous case (O1 MSE ≥ 0.10 · LE_sink = 9.35 W m⁻² without
  the vapour) does not arise. V2 is not "mechanism unresolved".
- **O5 carries water down, not up.** Its net exchange dries the layer below
  the inversion at 0.72 of the surface evaporation.
- **So V2 finds ventilation: V3's condition ("only if V2 finds no
  ventilation") is not met, and V3 does not run.** The verdict "missing
  ventilation" is not reachable.

**What the budget does not decompose.** The terms above are measured across
the inversion only. Of the surface evaporation (3.74e-5), O4's export (0.78e-5),
O5's net exchange (2.70e-5) and precipitation (0.01e-5) account for 3.49e-5;
the remaining 0.25e-5 (7 percent) is not attributed here. It would sit in the
circulation's lateral exchange within the boundary-layer layers and in storage,
neither of which V0 sums.

**V1 is running** (branches a, b and c, 20000 steps each, from V0's final
state). It reports the same budget per branch.

## Interim disposition after V0/V2 (the user's, 2026-10-19, verbatim)

> The sinking column's low evaporation is not caused by an unventilated
> boundary layer. Eddy diffusion carries a material upward moisture flux
> across the inversion, while subsidence produces a larger downward flux.
> The inversion remains sharp because the net transport and local
> thermodynamic response maintain it, not because all upward exchange is
> absent. The remaining V1 switches test which existing operator controls
> the moist boundary-layer state and the excessive fixed-surface
> redistribution requirement.
>
> No protocol change is warranted from V0/V2.

**A correction to the V2 wording above, recorded rather than edited.** V2
says the unattributed 7 percent of the sinking surface evaporation "would
sit in the circulation's lateral exchange within the boundary-layer layers
and in storage". That was not measured. **It is unresolved** until the
circulation's moisture exchange below the inversion is totalled.

## The surface humidity deficit, from V0 (recorded before any V1 result)

**The evaporation formula.** The model evaporates into the lowest layer
only, E = k_E · (q_sat(T_s, p₁) − q₁) (`water.py`), with k_E shared by
both columns.

**Computed from V0's 50 saved late states** (every 240 steps over the last
fifth, at the state each step's tendency was evaluated from):

| | rising | sinking |
|---|---:|---:|
| T_s (held), K | 300.26 | 298.79 |
| lowest-layer T, K | 300.08 | 298.86 |
| q_sat(T_s), kg kg⁻¹ | 0.022693 | 0.020786 |
| q₁ | 0.016303 | 0.017033 |
| lowest-layer RH | 0.726 | 0.816 |
| **deficit q_sat(T_s) − q₁** | **0.006389** | **0.003753** |
| E, kg m⁻² s⁻¹ (late mean, all steps) | 6.39e-5 | 3.74e-5 |

**Consistency check.** The deficit ratio is 0.587 and the evaporation ratio
is 0.584. They agree to the sampling, and the implied k_E is 0.0100 and
0.0100 in the two columns, as it must be.

**Where the sinking column's smaller deficit comes from** (−0.002636):

- **−0.001906 (72 percent) from the colder held surface,** through q_sat(T_s);
- **−0.000730 (28 percent) from the moister lowest layer.**

**So most of the evaporation gap comes from the prescribed SST contrast,
through Clausius–Clapeyron, and is not reachable by any boundary-layer
operator.** The lowest sinking layer is also 0.07 K warmer than its
surface, which is why H is slightly negative there.

**What a substantial exchange response through evaporation alone would
need.** This is arithmetic from V0's numbers and the frozen threshold, not
a new threshold. With C_conv unchanged:

- raising X_sink by 61.01 W m⁻² needs E up by 2.44e-5 kg m⁻² s⁻¹ (L = 2.501e6 J kg⁻¹);
- that is a deficit up by 0.00245;
- that needs q₁ down to 0.01458, a lowest-layer RH of about 0.70.

That is drier than the rising column's lowest layer (0.726). Making the
sinking surface layer exactly as humid as the rising one would raise LE by
about 18 W m⁻², under a third of the threshold. A V1 branch can still reach
"substantial" through C_conv, which X includes, or through ΔF.

## V1 result (2026-10-19)

**Validity.** All three branches completed 20000 steps from V0's final
state.

- **Conservation residual:** 8.9e-16 (a), 9.2e-16 (b), 8.4e-16 (c).
- **Holding-relation closure:** 6.7e-16, 9.6e-16 and 1.3e-15.
- **Restart check:** one step from V0's final state reproduces V0's water
  store to within one step's evaporation (0.34 kg m⁻²), so the injected
  vapour and water store are used. The only other prognostic states are
  cumulative ledgers, which start at zero.
- **Settling:** the sinking column's F_hold differs by under 0.4 W m⁻²
  between the two halves of each late window.

**Late means (last 4000 steps), against V0:**

| | V0 | V1a | V1b | V1c |
|---|---:|---:|---:|---:|
| sinking LE, W m⁻² | 93.46 | 54.51 | 61.08 | 92.75 |
| sinking C_conv | 19.17 | 0 | 13.87 | 19.52 |
| sinking H | −0.69 | 4.91 | −0.40 | −0.67 |
| sinking R_net,sfc | 193.31 | 183.10 | 180.58 | 193.17 |
| **X_sink** | **112.62** | **54.51** | **74.95** | **112.26** |
| **X_sink − X_sink,V0** (substantial ≥ +61.01) | | **−58.11** | **−37.67** | **−0.36** |
| X_rise | 234.65 | 266.38 | 255.55 | 235.21 |
| **ΔF** (substantial ≤ 71.34) | **142.68** | **218.94** | **189.60** | **143.37** |
| sinking E, kg m⁻² s⁻¹ | 3.74e-5 | 2.18e-5 | 2.44e-5 | 3.71e-5 |
| lowest-layer RH, sinking | 0.816 | 0.922 | 0.880 | 0.818 |
| deficit q_sat(T_s) − q₁, sinking | 0.003753 | 0.002185 | 0.002454 | 0.003725 |
| mixed-layer top, hPa | 827.7 | 996.4 (none) | 827.7 | 827.7 |
| observational inversion, hPa | 810.8 | 979.5 | 979.5 | 810.8 |
| the two coincide | 1.00 | 1.00 | 0.00 | 1.00 |
| RH jump at the inversion | 0.45 | 0.31 | 0.32 | 0.59 |
| boundary-layer mean RH (to the top) | 0.895 | 0.931 | 0.548 | 0.897 |

**Sinking RH, layers 0–6 (late-state samples):**

| | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|---:|---:|
| V0 | 0.816 | 0.829 | 0.863 | 0.902 | 0.945 | 0.973 | 0.538 |
| V1a | 0.922 | 0.621 | 0.517 | 0.504 | 0.539 | 0.606 | 0.529 |
| V1b | 0.880 | 0.562 | 0.461 | 0.444 | 0.456 | 0.473 | 0.424 |
| V1c | 0.818 | 0.830 | 0.865 | 0.904 | 0.947 | 0.973 | 0.412 |

**Cross-inversion budget, sinking, upward positive.** Each branch is measured
across its own inversion. In V1a and V1b that inversion is interface 1, so
these fluxes are not comparable with V0's at interface 6.

| | Fq O1 | Fq O4 (ratio to E_sink,V0) | Fq O5, −M·q | O5 net exchange | Fh O4, W m⁻² | Fh O5 exchange |
|---|---:|---:|---:|---:|---:|---:|
| V1a | 3e-25 | 9.12e-6 (0.244) | −2.57e-5 | −1.27e-5 | 22.8 | −26.9 |
| V1b | 0 | 1.05e-5 (0.281) | −2.23e-5 | −1.39e-5 | 26.2 | −31.5 |
| V1c | 7e-20 | 0 | −2.44e-5 | −3.30e-5 | 0 | −72.4 |

(The `material_ratio` the V1 script printed divides by each branch's own E.
The frozen rule divides by V0's, which is what the parenthesis gives.)

**V1a's declared confound, checked.** Its sinking free-troposphere
temperature (700–200 hPa) differs from V0's by −0.03 K on average and
1.0 K at most. The lowest three layers differ by −0.56, 0.00 and +0.48 K.
The confound sits in the boundary layer, as the erratum expected.

**Under the frozen thresholds:**

- **No branch gives a substantial change.** Both thresholds are signed. V1a
  and V1b move X_sink down and ΔF up, the opposite direction.
- **V1a's magnitude is recorded as found:** |ΔX_sink| = 58.11, 95 percent
  of the threshold amount, in the direction that widens the gap.
- **V1c is negligible on both:** |ΔX_sink| = 0.36 < 6.10, and ΔF changes
  by +0.69 < 7.13.

**What the branches show.**

- **O1 is what carries the surface's evaporation up into the boundary layer.**
  With its surface coupling off (V1a), or its vapour transport off (V1b), the
  water stays in the lowest layer, that layer moistens (RH 0.92 and 0.88), the
  surface deficit falls by 35 to 42 percent, and evaporation falls with it.
  The layers above dry to RH 0.44–0.62 and the observational inversion drops to
  the top of the lowest layer.
- **O1 crosses no inversion in any branch.** It ventilates the surface layer
  into the boundary layer, not the boundary layer into the free troposphere.
- **O4 is not what sets the sinking boundary layer's state.** Removing the
  only material upward transport across the inversion (V1c) leaves the
  lowest-layer RH at 0.818 against 0.816, and evaporation within 1 percent.
  Subsidence's net exchange takes over (−3.30e-5 against −2.70e-5), and the
  inversion's RH jump sharpens from 0.45 to 0.59.
- **The rising column responds too.** X_rise rises by 32 and 21 W m⁻² in V1a
  and V1b, through the circulation, so ΔF widens on both sides.

**Taken with the deficit decomposition recorded before V1:** within these
operators, the sinking column's evaporation is held up by O1, not
suppressed by any of them. 72 percent of its deficit shortfall comes from
the held SST contrast.

**The decision logic: what the frozen text settles and what it leaves open.**

- **Settled:** no V1 branch gives a substantial change, so "operator
  artefact" is not reached. "Missing ventilation" was already excluded by V2.
- **Left open:** V4 runs "only if V1 implicates the depth" and V5 "only if
  V1 implicates O1". The protocol does not define "implicates". V1a and V1b
  show large responses to O1, all in the direction that widens the gap.
- **If neither runs:** the frozen verdict is **no attribution**.
- **If V5 runs:** its question is largely answered already. O1's export
  through the observational inversion is zero to rounding in V0 and every
  branch, and V1b shows it redistributes the surface layer's vapour within
  the boundary layer. That would meet the literal "convection redistribution
  defect" condition. But V1b shows the redistribution raises the sinking
  evaporation, so the name would not describe the result.
- **This is left to the user.**
