# Protocol: regime-dependent surface exchange (E0-E2)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The design is the
user's (2026-10-20), mapped onto the exchange relation this model actually
has. The user approved it to freeze with two decisions:

1. **the wind source:** OAFlux v3 monthly 10 m wind speed;
2. **the mapping:** both mappings, with the relative one primary. Only the
   relative mapping enters the verdict and E2; the absolute one is reported
   as a second, offline reading.

**It follows from** `descending-boundary-layer-ventilation.md`, closed as no
attribution (3abbd30):

- **The sinking boundary layer is ventilated,** and none of the tested
  ventilation operators explains its low evaporation.
- **Convective mixing supports the sinking evaporation.** Removing it lowered
  X_sink by 37.67–58.11 W m⁻².
- **72 percent of the sinking column's smaller humidity deficit** comes from
  the held SST contrast. That share is structural to S1's prescribed
  boundary, not an identified model error.
- **What the closed protocol could not separate:** "k_H, and k_E = k_H / c_p,
  are shared by both columns. A wind-speed difference between the regimes is
  not represented."

## The question

> Does representing the stronger near-surface winds of the subsiding
> trade-wind regime increase sinking-column evaporation enough to materially
> reduce the fixed-surface redistribution requirement, without changing the
> atmospheric circulation or adding a humidity adjustment?

## The exchange relation this model has

**One coefficient carries both surface fluxes:**

- H = k_H (T_s − T₁), applied by the convection module;
- E = k_E (q_sat(T_s, p₁) − q₁), into the lowest layer only, applied by the
  water module;
- k_E = k_H / c_p exactly (the Reynolds analogy,
  `evaporation_coefficient_for`), so the latent and sensible fluxes move
  together;
- k_H = ρ c_p C |U| = 10 W m⁻² K⁻¹, which gives k_E = 0.00996 kg m⁻² s⁻¹;
- LE = L_v E, with L_v = 2.501e6 J kg⁻¹, constant.

**Not represented separately:** ρ, C and |U| exist only inside the one
product k_H. Its declared range, 6.0–14.5 W m⁻² K⁻¹, is the arithmetic of C
from 1.0e-3 to 1.5e-3 and |U| from 5 to 8 m s⁻¹ at ρ c_p = 1205 J m⁻³ K⁻¹
(`SENSIBLE_HEAT_COEFFICIENT_RANGE`). That range is asserted, not measured.
At C = 1.2e-3, the nominal k_H = 10 corresponds to |U| ≈ 6.9 m s⁻¹.

**So a regime wind enters as a per-column k_H,** and it moves H and LE
together, in the same proportion. Changing LE alone would break the Reynolds
analogy the model derives k_E from, and is not done.

**No new physics code is needed.** `two_column(regional_settings=...)`
already sets `sensible_heat_coefficient` per region, and k_E follows from it.

## E0: the baseline and its decomposition

**The baseline is V0** (`descending-boundary-layer-ventilation.md`), which
reproduces C0 bit for bit.

**The identity check.** The E2 code path, with the shared coefficient
written explicitly in each region (`sensible_heat_coefficient=10.0` for both),
runs C0's 60000-step configuration. Every shared printed value of
`v0_60000.json` must be identical.

**The decomposition, per column, from V0's late means and its 50 saved
late states:**

| factor | how it enters | rising | sinking |
|---|---|---|---|
| ρ, C, \|U\| | not represented separately | — | — |
| k_E | shared | 0.00996 | 0.00996 |
| q_sat(T_s, p₁) | held surface | (from V0) | (from V0) |
| q₁ | lowest layer | (from V0) | (from V0) |
| deficit | q_sat − q₁ | (from V0) | (from V0) |
| L_v | constant | 2.501e6 | 2.501e6 |

**What E0 can already state, from V0:** the coefficient contributes nothing
to the evaporation contrast, because it is shared. The whole contrast is the
deficit: −0.00191 from the colder surface and −0.00073 from the moister
lowest layer. E0's table restates these with the identity check's
numbers.

## E1: the diagnostic observed-wind replay (offline, no feedback)

**The data.** OAFlux v3 monthly 10 m wind speed, `ws` (Yu and Weller 2007,
*BAMS* 88, 527–539; WHOI OAFlux project).

- **Years:** 1991–2020, the masks' period.
- **Source:** `https://scienceweb.whoi.edu/oaflux/data_v3/monthly/turbulence/ws_oaflux_YYYY.nc.gz`.
- **Registered before use:** each file in `data/sources.yaml` with its
  SHA-256, verified before parsing. That is 30 `ws` files, plus the
  evaluation fields below (`lh`, `sh`, `qa`, `ts`), 150 files in all.
- **Provenance:** OAFlux's wind is an objective analysis that blends
  satellite winds (SSM/I, QuikSCAT, ASCAT, WindSat and others) with
  reanalysis. It is registered as such, with the blend stated in its notes.

**The masks are the existing ones.** The ±20° annual-mean NCEP/NCAR R1
1991–2020 ω500 masks:

- rising, ω500 < 0, area fraction 0.558; sinking, the rest;
- cos-latitude weighting;
- the 1° OAFlux grid mapped onto the masks by the same bilinear
  interpolation of ω500 the ISCCP composite used;
- **ocean cells only.** OAFlux has no land. S1 treats the land fractions
  (0.25 and 0.21) as ocean surface.

**The composited quantity:** the annual mean of the monthly mean scalar wind
speed in each regime, W_rise and W_sink. The ratio is R = W_sink / W_rise.

**Mapping the winds onto the model, relative (primary).** The band's
area-weighted coefficient stays at its current value, and only the regime
contrast is imposed:

- k_i = k_H · W_i / W̄, with W̄ = σ W_rise + (1 − σ) W_sink and σ = 0.558;
- equivalently r_i = k_i / k_H, with σ r_rise + (1 − σ) r_sink = 1.

**Why relative:** the band-mean level of k_H is an asserted value inside a
declared range, and every recorded result sits on it. Moving the level and
the contrast at once would mix two questions. This protocol asks only about
the contrast.

**The absolute mapping (second reading, offline only, never in the verdict
or E2).** Each column's coefficient comes from its own observed wind and a
published neutral bulk formula. Large and Yeager (2004), *Diurnal to decadal
global forcing for ocean and sea-ice models: the data sets and flux
climatologies*, NCAR/TN-460+STR, eq. (6a–b) and section 2:

- 1000 C_D = 2.70 / U + 0.142 + U / 13.09, with U the 10 m neutral wind in
  m s⁻¹;
- 1000 C_E = 34.6 √C_D;
- ρ_a = 1.22 kg m⁻³.

So k_E,i = ρ_a C_E(W_i) W_i, and k_H,i = c_p k_E,i, with the model's c_p (the
Reynolds analogy the model uses). The offline relations below then apply
with r_i = k_E,i / 0.00996. The level and the contrast both change here, so
this reading answers a different question: what the observed winds imply
with a published coefficient, not what the regime contrast alone does.

**Approximations in the absolute reading:**

- the composited monthly-mean scalar wind goes into a formula that is
  nonlinear in U;
- the neutral 10 m coefficient is used without a stability correction;
- OAFlux's wind is treated as the 10 m neutral wind.

**The offline replay is exact at fixed state,** because E and H are linear in
k at a given state. With V0's late means:

- LE_i* = r_i · LE_i;
- H_i* = r_i · H_i;
- X_sink* = r_sink · LE_sink + C_conv,sink;
- F_hold,i* = F_hold,i + (r_i − 1)(LE_i + H_i);
- ΔF* = F_hold,rise* − F_hold,sink*.

**Reported, for both mappings:** W_rise, W_sink, R, r_rise and r_sink;
ΔLE_sink, ΔLE_rise, X_sink* − X_sink,V0, and ΔF*. For the absolute mapping
also C_E(W_i), k_E,i and k_H,i.

**Never done:** scaling the winds, or choosing years, cells or a mapping, to
reach a threshold or the observed ocean-transport contrast.

**Evaluation only, composited the same way, never fed back:** OAFlux `lh`,
`sh`, `qa` and `ts` per regime. These give the observed regime LE and H, and
the observed deficit q_sat(ts) − qa. OAFlux applies a 0.98 salinity factor to
q_sat; the model uses pure water. These are reported beside the model's V0
values.

**What the result would have to be, recorded before any data is read.**
Under the relative mapping, with V0's numbers:

| outcome | needs r_sink | r_rise | R = W_sink / W_rise |
|---|---:|---:|---:|
| X_sink up by 61.01 (substantial) | 1.653 | 0.483 | 3.42 |
| ΔF down to 71.34 (substantial) | 1.323 | 0.744 | 1.78 |
| X_sink up by 6.10 (material) | 1.065 | 0.948 | 1.12 |
| ΔF down by 7.13 (material) | 1.032 | 0.974 | 1.06 |

The ΔF route is easier because it acts on both columns: the rising column's
LE falls while the sinking one's rises.

## E2: the interactive fixed-SST test (only if E1 is material)

**Runs only if E1 predicts a material response in the improvement
direction:** X_sink* − X_sink,V0 ≥ 0.1 × 61.01 = 6.10 W m⁻², or
ΔF_V0 − ΔF* ≥ 0.1 × 71.34 = 7.13 W m⁻².

**The only change:** `sensible_heat_coefficient` = k_rise and k_sink from
E1, per region.

**Held fixed:**

- surface temperatures, 300.26 and 298.79 K;
- WTG coupling, L = 734 km;
- the area fractions;
- the cloud fields;
- all atmospheric physics, including convective mixing and eddy
  diffusion;
- the initial state.

**Two runs, both from V0's final state, 20000 steps, read over the last
fifth:**

- **E2-null:** k = 10 in both regions. This is the matched reference for the
  restart.
- **E2:** k_rise and k_sink from E1.

**Reported for both, against V0 and E2-null:**

- sinking LE, H and C_conv, X_sink, and ΔF;
- lowest-layer RH, q₁ and the deficit, in both columns;
- the inversion and the mixed-layer top;
- ω500 subsidence;
- free-tropospheric RH;
- precipitable water;
- the circulation's orientation;
- both fixed-surface boundary transfers;
- boundedness and settling;
- conservation and the holding-relation closure.

## Thresholds, from V0, as frozen in the ventilation protocol

| term | rule | value |
|---|---|---:|
| substantial exchange response | X_sink,test − X_sink,V0 ≥ 0.5 · G0 | ≥ 61.01 W m⁻² |
| substantial redistribution response | ΔF_test ≤ 0.5 · ΔF_V0 | ≤ 71.34 W m⁻² |
| material (E2's trigger) | at least 0.1 of either substantial amount, in the improvement direction | 6.10 / 7.13 W m⁻² |
| negligible | below 0.1 of the corresponding substantial amount | |

E2 is judged against E2-null and against V0. If the two references give
different verdicts, both are recorded and the verdict is **unresolved**.

## Decision rules

| verdict | condition |
|---|---|
| **Surface-exchange explanation supported** | E2 raises X_sink by ≥ 61.01, or brings ΔF to ≤ 71.34, with a valid, settled circulation and exact accounting |
| **Material but insufficient** | E2 moves X_sink or ΔF in the improvement direction by a material amount, but meets neither substantial threshold |
| **Rejected** | E1 or E2 is negligible, or worsens ΔF |
| **Atmospheric compensation** | E1 predicts a substantial response; E2's response is less than half of E1's prediction, with the sinking lowest-layer deficit lower than in E2-null |
| **Invalid** | the identity check, conservation, the holding-relation closure, or settling fails |

**If E1 is not material,** E2 does not run and the verdict is **rejected**,
from E1.

## What this cannot separate

- **C is held equal in both regimes.** Stability is not represented: the
  sinking surface is slightly colder than its lowest layer (H < 0), which
  would lower C there, and the rising regime's convective gustiness would
  raise its effective |U|. Both act against the sinking column.
- **The mean scalar wind is not the flux-weighted wind.** The covariance of
  |U| and the deficit within a month and within a regime is not represented.
- **OAFlux's wind is a blended analysis,** not a direct measurement, and
  satellite winds are sparse under heavy rain, which falls mostly in the
  rising regime.
- **The land in each mask** has no OAFlux wind.
- **The band-mean level of k_H** stays at its asserted 10 W m⁻² K⁻¹. The
  relative mapping holds it fixed on purpose. The absolute reading shows
  what a published coefficient implies instead, offline only.
- **Everything is under prescribed surface temperatures.**

---

## E1 result (2026-10-20)

**Data validity.**

- **Checksums:** all 150 OAFlux files and the ω500 file match
  `data/sources.yaml`, verified before parsing.
- **Two formats:** 1991–2015 are netCDF-3 (scaled integers); 2016–2020 are
  netCDF-4 (unscaled, 1e36 missing, lh and sh named `lhflx` and `shflx`).
  The registry's variable names for those ten files were corrected
  (0efbdff).
- **Grid:** identical across all years, which the script asserts per file.
- **The masks on the 1° grid:** the rising fraction is 0.576 over all cells
  and 0.565 over the ocean cells used (the masks' value on NCEP's 2.5° grid
  is 0.558; σ = 0.558 is used in the mapping, as frozen). OAFlux covers
  0.765 of the rising mask's area and 0.801 of the sinking mask's.

**The regime winds (1991–2020 monthly climatology, annual mean,
cos-latitude weighted, ocean cells):**

| | rising | sinking |
|---|---:|---:|
| W, 10 m neutral wind, m s⁻¹ | **6.143** | **6.799** |
| mean OAFlux error estimate, m s⁻¹ | 0.36 | 0.29 |

**R = W_sink / W_rise = 1.107,** and W̄ = 6.433 m s⁻¹.

**The offline replay of V0's late means:**

| | V0 | relative (primary) | absolute (second reading) |
|---|---:|---:|---:|
| k_H rising, W m⁻² K⁻¹ | 10 | 9.550 | 8.445 |
| k_H sinking | 10 | 10.569 | 9.380 |
| r_rise / r_sink | 1 / 1 | 0.955 / 1.057 | 0.844 / 0.938 |
| LE rising, W m⁻² | 159.93 | 152.73 | 135.06 |
| LE sinking | 93.46 | 98.77 | 87.66 |
| **X_sink − X_sink,V0** | | **+5.31** | **−5.79** |
| **ΔF** | **142.68** | **130.12** | **123.27** |
| **ΔF reduction** | | **12.56** | **19.41** |

For the absolute reading, Large and Yeager's C_D and C_E are 1.05e-3 and
1.12e-3 (rising), 1.06e-3 and 1.13e-3 (sinking).

**Judged under the frozen rules (relative mapping):**

- **Substantial: no.** X_sink + 5.31 < 61.01, and ΔF 130.12 > 71.34.
- **Material: yes, through ΔF.** The reduction of 12.56 is at least 7.13.
  The X_sink route, at 5.31, is below its material line of 6.10.
- **So E2 runs,** with k_rise = 9.549530196129602 and
  k_sink = 10.568692648325069.

The observed R of 1.107 sits just under the 1.12 that the pre-registered
table needed for material through X_sink, and above the 1.06 for ΔF. It is
far below the 1.78 that ΔF needed for substantial.

**The absolute reading** puts both coefficients below the model's 10,
because the published C_E at these winds is about 1.12e-3. ΔF falls further
(by 19.41), but most of the fall comes from the rising column's lower
coefficient: X_sink drops by 5.79. It answers the level question, not the
contrast question, and does not enter the verdict.

**Evaluation only, never fed back.** OAFlux composited the same way, beside
V0. OAFlux's q is at 2 m and its SST is ocean only, so it is warmer than the
held NCEP skin temperatures, which include land.

| | observed rising | observed sinking | V0 rising | V0 sinking |
|---|---:|---:|---:|---:|
| surface temperature, K | 301.50 (SST) | 299.05 (SST) | 300.26 (held) | 298.79 (held) |
| near-surface q, kg kg⁻¹ | 0.01853 | 0.01616 | 0.01630 (q₁) | 0.01703 (q₁) |
| deficit, pure-water q_sat at 996 hPa | 0.00593 | 0.00512 | 0.00639 | 0.00375 |
| deficit with OAFlux's 0.98 factor | 0.00544 | 0.00470 | | |
| LE, W m⁻² | 116.97 | 113.42 | 159.93 | 93.46 |
| H, W m⁻² | 8.33 | 7.24 | 1.85 | −0.69 |

**What the evaluation shows:**

- **The observed latent heat flux is nearly equal in the two regimes,**
  with a ratio of 0.970. The model's ratio is 0.584.
- **The observed deficit ratio, sinking over rising, is 0.863.** The
  model's is 0.587.
- **Observed sinking near-surface air is drier than rising air,** by
  0.0024 kg kg⁻¹. The model's sinking lowest layer is moister, by 0.0007.
- **The bulk relation checks against the observed fluxes.** The observed
  deficit ratio times the observed wind ratio gives 0.863 × 1.107 = 0.955,
  against the observed LE ratio of 0.970.
- **So most of the model's evaporation contrast, against observation,
  lies in the deficit contrast,** not the wind contrast. The observed wind
  contrast is 11 percent, while the model's deficit contrast exceeds the
  observed one by a factor of 0.863 / 0.587 = 1.47.
- **This is evaluation and is not judged.** The held surface temperatures
  are not OAFlux's, and q at 2 m is not a layer mean.

**E2 and E2-null are running** (20000 steps each, from V0's final state),
together with E0's 60000-step identity run.

## E2 result (2026-10-20): computed with the composition defect

**Read this result with that in mind.** E2 and E2-null ran on the code
before `advection-composition.md`'s fix (ba42c7c). In this configuration the
water and convection modules did not compose with the circulation's
in-step transport. The comparison is like for like (V0, E2-null and E2 all
carry the defect), but the baseline itself does.

**Validity.**

- **Both runs** completed 20000 steps from V0's final state.
- **Conservation residual:** 8.2e-16 (E2-null) and 7.5e-16 (E2).
- **Holding-relation closure:** 1.9e-15 and 1.5e-15.
- **Orientation (G5) held:** sinking ω500 positive on every late day (26.44
  and 26.34 hPa/day).
- **Restart offset:** E2-null against V0 is 0.28 W m⁻² in X_sink and 0.30
  in ΔF, which shows how far a restart moves the late mean.
- **Settling:** the sinking F_hold differs between the two halves of the
  late window by 0.8 (E2-null) and 1.1 W m⁻² (E2).

**Late means:**

| | V0 | E2-null | E2 |
|---|---:|---:|---:|
| k_H rising / sinking | 10 / 10 | 10 / 10 | 9.550 / 10.569 |
| LE rising, W m⁻² | 159.93 | 159.83 | 157.67 |
| LE sinking | 93.46 | 93.67 | 96.26 |
| H rising / sinking | 1.85 / −0.69 | 1.85 / −0.69 | 1.82 / −0.79 |
| C_conv rising | 74.72 | 74.79 | 76.05 |
| C_conv sinking | 19.17 | 19.24 | 17.86 |
| X_sink | 112.62 | 112.91 | 114.13 |
| ΔF | 142.68 | 142.38 | 140.95 |
| q₁ rising / sinking, kg kg⁻¹ | 0.01630 / 0.01703 | 0.01625 / 0.01701 | 0.01600 / 0.01713 |
| deficit rising / sinking | 0.00639 / 0.00375 | 0.00644 / 0.00377 | 0.00669 / 0.00366 |
| lowest-layer RH rising / sinking | 0.726 / 0.816 | 0.724 / 0.815 | 0.713 / 0.820 |

**Against the offline prediction:**

| | E1 offline | E2 − E2-null | E2 − V0 |
|---|---:|---:|---:|
| X_sink change | +5.31 | **+1.22** | **+1.50** |
| ΔF reduction | 12.56 | **1.43** | **1.73** |

**Under the frozen rules:**

- **Negligible against both references.** The material amounts are 6.10 for
  X_sink and 7.13 for ΔF. E2 moves X_sink by 1.22–1.50 and ΔF by 1.43–1.73,
  in the improvement direction.
- **The two references agree,** so the verdict is not unresolved.
- **"Atmospheric compensation" is not reachable,** because it needs E1 to
  predict a substantial response, and E1 predicted a material one.
- **Frozen verdict: rejected.**

**What the interactive atmosphere did, recorded as found.** The coupled
response was about an eighth of the offline one. Three things erased most of
the offline gain:

- **The sinking lowest layer moistened** (q₁ 0.01701 → 0.01713, RH 0.815 →
  0.820), so its deficit fell (0.00377 → 0.00366) and evaporation rose by
  2.8 percent where the coefficient rose 5.7 percent.
- **The rising lowest layer dried** (q₁ 0.01625 → 0.01600), so its deficit
  grew and its evaporation fell by 1.4 percent where the coefficient fell
  4.5 percent.
- **Convection's surface exchange moved against the latent heat change in
  both columns:** C_conv fell by 1.37 W m⁻² in the sinking column and rose by
  1.26 in the rising one.

**The near-surface humidity contrast stayed reversed,** q₁,sink − q₁,rise =
+0.00113, against the observed −0.0024.

**Whether the composition fix changes this verdict is not known.** It would
take an E2 rerun from the corrected baseline (A1). That, and closing this
protocol, are the user's decisions.

## Closure (the user's disposition, 2026-10-20, verbatim)

**Closed as rejected, exactly as frozen. E2 is not rerun on the
composition-fixed code by default.**

> The regime-dependent surface-exchange explanation is rejected. The
> observed 11% wind-speed contrast produces a material offline response but
> a negligible coupled response. Relative to E2-null, sinking-column X_sink
> increases by only 1.22 W m⁻², and ΔF decreases by only 1.43 W m⁻², both
> below the frozen material thresholds of 6.10 and 7.13 W m⁻². The coupled
> atmosphere compensates by moistening the sinking lowest layer and drying
> the rising lowest layer, preserving the modeled near-surface humidity
> contrast with the wrong sign. The dominant error is therefore not the
> shared exchange coefficient or missing wind contrast, but the
> boundary-layer humidity partition between the columns.

**The central unresolved diagnostic:** q₁,sink − q₁,rise = +0.0011 kg kg⁻¹
in the model, against −0.0024 observed.

**"Atmospheric compensation" stays descriptive.** The formal verdict needed
E1 to predict a substantial effect, and it did not. The verdict is
**rejected, with coupled boundary-layer compensation observed**:

- the sinking k_H rose 5.7 percent, but its evaporation rose only 2.8
  percent, because the lowest layer moistened;
- the rising k_H fell 4.5 percent, but its evaporation fell only 1.4
  percent, because the lowest layer dried;
- convective surface exchange opposed the latent-flux response by about
  1.3 W m⁻² in each column.

## The E2 rerun condition (frozen before any A1 result is read)

**E2 is rerun on the fixed code only if A1** (`advection-composition.md`)
**changes, relative to V0** (the baseline E2 restarted from), **at least one
of the following.** The list is the user's; the numbers marked "reading" are
added at freezing so that the condition can be judged mechanically.

| criterion | rule |
|---|---|
| sinking lowest layer | q₁,sink or its deficit q_sat(T_s) − q₁ changes by ≥ 10 percent of V0's |
| rising lowest layer | q₁,rise or its deficit changes by ≥ 10 percent of V0's |
| sinking evaporation | \|LE_sink,A1 − LE_sink,V0\| ≥ 6.10 W m⁻² |
| the redistribution | \|ΔF_A1 − ΔF_V0\| ≥ 7.13 W m⁻² |
| orientation (reading) | S1's G5 fails in A1 |
| inversion location (reading) | the sinking observational inversion's modal interface differs from V0's (interface 6, 810.8 hPa) |
| subsidence (reading) | the late-mean sinking ω500 changes by ≥ 10 percent of V0's (26.37 hPa/day) |
| quantities the fix uses (reading) | the fix acts on the exchange response only through the lowest-layer state and the surface fluxes, so this criterion is judged by the four rows above; no separate number |

**The lowest-layer quantities are computed as in E2's result:** from the 50
saved late states, at the state each step's tendency was evaluated from.

**If none holds,** E2's frozen verdict stands for the fixed code as well,
and this protocol takes no further runs. **If any holds,** E2 and E2-null
are rerun from A1's final state under this protocol's frozen rules, and both
results are recorded.

**E0's identity run** is kept as the identity requirement for the
pre-fix branch. It no longer affects this verdict.

## E0 result (2026-10-20): the identity check passed

The E2 code path, with k_H = 10 written explicitly in each region, ran C0's
60000-step configuration on the pre-fix code. **All 1130 numeric values in
`e0_60000.json` are identical to `v0_60000.json`, bit for bit.** No value is
in one and not the other. So E2 and E2-null's per-region coefficient setting
is the baseline at k = 10, as the protocol required.
