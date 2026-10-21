# Protocol: bracketing convective vapour mixing between its parameter-free limits (X0-X2)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The user approved this
step (2026-10-21, "go"):

- a frozen comparison of declared structures;
- selected by a structural observable fixed beforehand, never by the sinking
  column's humidity;
- with every tolerance scaled to the round-off of the quantity it checks.

**It follows from** `rising-mixing-structure.md` (M0–M2, bf9c7d9; erratum 2 and
disposition, 2807260):

- **The rising column's surface-rooted group mixes vapour to a uniform q up to
  717 hPa on average.** That is 231 hPa above the observed cloud base (NCEP R1
  LCL, 948 hPa), 93 hPa above ISCCP-H's low-cloud tops (809 hPa), and 205 hPa
  above the model's own LCL.
- **Below 760 hPa its vertical q gradient is zero.**

**The operator already declares itself an upper bound.** `_vapour_transport`'s
docstring reads: "Mixing completely is therefore an upper bound on the
moistening."

## The two structures

Both are parameter-free limits of the convective vapour transport by the
surface-rooted group. **Nothing else differs:**

- the temperature adjustment and its groups;
- the elevated groups;
- radiation, the circulation, the surface and condensation.

| structure | the surface-rooted group's vapour is mixed over | status |
|---|---|---|
| **F, full** | every layer the temperature adjustment mixed (the current operator) | the declared upper bound |
| **S, subcloud** | only the layers whose centres lie at or below the lifting condensation level of the surface parcel | the lower bound: a subcloud mixed layer homogenises, a cloud layer above it does not |

**The surface parcel, in S:** layer 0's adjusted temperature and its
provisional (pre-mixing) specific humidity, at layer 0's pressure. It is
lifted dry-adiabatically with q conserved, to the pressure where it saturates
by the model's own saturation function.

- If that pressure lies above the group's top, S equals F for that step.
- If layer 0 is already saturated, S mixes no vapour in the surface group.
- When S truncates a group, the group's partial top takes no vapour.

**Neither structure is fitted to anything,** and no third structure is
introduced to interpolate between them. **An observed structure between the
two limits is an allowed outcome,** and would be recorded as the bracket
straddling the observations, not resolved by tuning.

## Code (before X0)

**`ConvectiveAdjustmentModule` gains `vapour_mixing_top`.** It takes
`"group"` (the default, F) or `"lifting_condensation_level"` (S), plumbed
through the column configurations as `diagnose_mixing_groups` is.

**S requires:**

- vapour transport;
- surface-based convection;
- no mass-flux convection.

Any other combination is refused.

**The vapour groups are kept apart from the heat groups.** The vapour mixing,
the provenance-class mixing and the published `convective_mixing_group` all
use the vapour groups. The temperature adjustment, its validity mask and
`surface_mixed_layer_top` keep the heat groups.

**A new equation record, before the module cites it:**
`EQ-LIFTING-CONDENSATION-LEVEL`. It is a diagnostic identity: a dry-adiabatic
lift with q conserved, saturating by `EQ-SATURATION-SPECIFIC-HUMIDITY`. The
record carries a validity domain, limiting cases and a correctness test
against an independent closed form: Bolton (1980), *Mon. Wea. Rev.* 108,
1046–1053, eq. 21.

**Tests:**

1. **The default is bit-identical:** every existing test passes unchanged.
2. **S never mixes the surface group's vapour above its LCL layer,** and
   conserves the column's water exactly.
3. **The published vapour groups replay the transport.**
4. **The LCL matches Bolton's closed form** within Bolton's own stated
   accuracy, carried into pressure.

The full gate passes before X0.

## X0: the run

**Configuration:** M0's exactly (R0's runner and saved fields), with
`vapour_mixing_top="lifting_condensation_level"`. It restarts from A1's final
state, runs 20000 steps, and the late fifth is saved. **F is M0 itself,** which
reproduces R0 bit for bit, so it is not rerun.

**Validity, all frozen.** Every tolerance is relative to the magnitude of what
is summed.

| check | rule |
|---|---|
| run | completes; the coupler's conservation and validity gates pass, as for every run |
| water | the column's water budget closes to ≤ 1e-12 of the column's vapour, every step |
| the truncation holds | the rising and sinking surface vapour groups never extend above the parcel's LCL layer, recomputed offline from the saved state |
| replay | the published vapour groups, mixed over the provisional vapour as D0 does, reproduce `label:convective_vapour_transport` to ≤ 1e-12 of the column's vapour |
| stationarity | the selection metric (below) from the late window's two halves differs by ≤ 0.05. If not, the run is continued 20000 steps from its own final state, once, and judged on that window |

## X1: the selection observable (structural, fixed before the run)

**The rising column's vertical humidity gradient, against NCEP R1's
rising regime,** in two windows:

- **925–850 hPa:** inside the observed subcloud and lower cloud layer;
- **850–700 hPa:** the cloud layer.

The NCEP values come from `m2_ncep.json`'s method, extended to the 925–850
window. Both are formed exactly as `rising-mixing-structure.md` M2(a) forms
them: the late-mean q of the committed state, log-p interpolated.

**The metric:** S_X = ½ Σ_windows \|g_model − g_NCEP\| / \|g_NCEP\|,
where g is dq/dp. A zero gradient scores 1 in that window. It is computed
for F (from M0) and for S (from X0).

**Not used for selection:**

- the LCL, because S's top is set by it;
- donor RH, because it is the restated outcome;
- the sinking column's humidity, contrast or deficit.

## Classifications, frozen

| class | condition |
|---|---|
| **Invalid** | an X0 check fails |
| **Subcloud structure preferred** | S_S ≤ 0.5 · S_F |
| **Full structure preferred** | S_F ≤ 0.5 · S_S |
| **Bracket straddles the observations** | neither is preferred, and in both windows the NCEP gradient lies between F's and S's |
| **Neither discriminated** | otherwise |

**A preference is a structural finding, not an adoption.** Changing the
default is the user's decision and a separate step.

## Reported always (consequences, never selection)

**The rising column:**

- q and RH at 850 and 700 hPa;
- the donor layers' RH, against NCEP's;
- precipitation, E and the precipitable water.

**The sinking column:**

- q₀,sink;
- q₀,sink − q₀,rise, with whether the reversed contrast is undone;
- the surface deficit, and LE_sink.

**And:**

- the circulation's ΔF and the return import at layers 5–8;
- the mean depth of the vapour group in both columns, and how often S
  truncated it.

## What this cannot separate

- **Two limits, not a closure.** Real cumulus transport lies between them.
  A straddle is the expected physical answer and would point to a
  transport with entrainment, not to either limit.
- **The surface parcel is one declared choice** (layer 0's adjusted T and
  pre-mixing q). A mixed-layer-mean parcel would place the LCL differently.
- **The heat adjustment is unchanged.** S mixes heat and vapour over
  different depths, which a real overturning would not.
- **The NCEP gradients are reanalysis,** on an annual-mean regime mask.
- **One restart of 20000 steps,** with the stationarity check above.
- **Two columns stand for a continuum of regimes,** under prescribed surface
  temperatures.

---

## Erratum 1: the subcloud parcel (the user's decision, 2026-10-21, before X0)

**Found in a 40-step smoke test of the implementation, before X0.**

**S as frozen collapses.** The frozen surface parcel is layer 0's
pre-mixing q, and layer 0 receives each step's evaporation after the mixing.
So:

1. the parcel gets moister every step;
2. its LCL falls into layer 0;
3. the cut then removes the surface group's vapour mixing altogether;
4. layer 0 traps the next step's evaporation, and the loop repeats.

**The measurement:** restarted from A1's state, the rising column's layer 0
went from RH 0.72 to 0.97 within about 30 steps, and the late smoke steps
published no surface vapour group in either column. S as frozen is
therefore not a subcloud mixed layer but **no surface vapour mixing**.

**Both readings, recorded.**

- **As frozen:** the parcel is layer 0's adjusted T and pre-mixing q, at
  layer 0's pressure.
- **As corrected:** the parcel is **the subcloud layer's own mean,** found
  as a fixed point. It is still parameter-free, and it is what "a subcloud
  mixed layer homogenises" means:
  1. start with k at the top of the heat adjustment's surface group;
  2. form the parcel from layer 0's adjusted T and pressure, with the
     mass-weighted mean of the provisional (pre-mixing) specific humidity
     over layers 0..k;
  3. if that parcel is still unsaturated at layer k's pressure, the
     subcloud layer is 0..k, and S mixes the vapour over it;
  4. otherwise k becomes the highest layer whose centre lies at or below
     the parcel's LCL, and go back to step 2;
  5. k only decreases, so the loop terminates. If k reaches 0, the surface
     group mixes no vapour.

**Everything else in the protocol is unchanged:** X0's checks, X1's
selection observable and windows, the classes and the reported
consequences. The truncation check recomputes this fixed point offline.

**Also recorded, from the same smoke test.** The validity-domain gate
judged the new relation's p_0 against every layer's pressure, down to
107 Pa, because the symbol was inferred from `atmosphere_layer_pressure`.
As the project's rule requires, the module now publishes its parcel's T, q,
p and condensation pressure (NaN on steps where the level is not
evaluated) and binds the relation to them. The declared domain is not
widened.

## Recorded before X0: S collapses under erratum 1 as well (2026-10-21)

**A 40-step smoke test of erratum 1's fixed point,** restarted from A1's
state:

- **the validity-domain gate passes,** with the relation bound to the
  published parcel;
- **the cut starts where expected:** layer 2 in the rising column, a
  subcloud layer about 100 hPa deep;
- **by the late smoke steps (33–40), the fixed point has descended to
  layer 0 in both columns.** No surface vapour group is published, the heat
  groups still reach layers 13–15 (rising) and 5 (sinking), and layer-0 RH
  is 0.97 (rising) and 0.96 (sinking).

**The mechanism, as measured.** Under S, the subcloud layer's only exits
for vapour through its top are the circulation and eddy diffusion.

1. Each step's evaporation is shared by 3 layers instead of about 15.
2. The subcloud layer moistens about five times faster.
3. Its mean parcel's LCL descends, and the cut follows it down.
4. Once layer 0 is above RH ≈ 0.9, the mean parcel over layers 0–1
   condenses below layer 1's centre (963 hPa), and the subcloud layer has
   no depth left.

**What this establishes, before any X0 run.** The lower limit, with no
vapour transport above cloud base, is not a viable state of this column. It
is the transport out through the subcloud layer's top that holds the
subcloud layer below saturation. So the two-limit bracket cannot be formed
as frozen. The comparison it was to make would set a degenerate state
against F.

X0 has not been run. How to proceed is the user's decision.
