# Protocol: resolution check and radiative attribution of the reference column's moist warm bias (R1-R2)

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards. The questions and the
attribution design are the user's (2026-10-21).

**One adjustment, made for the reason stated in R2: the attribution works in
TOA net, and splits the temperature effect into Planck and lapse-rate parts.**
Both states are equilibrated, so their totals nearly cancel.

## The status this starts from

> The reference column is +3.9 K from a dry start, where zero
> upper-tropospheric humidity is an initial-condition artifact, and +9.9 K
> once externally constrained vapor diffusion removes that memory at 30
> layers. The diffusion configuration still fails its independent
> stratospheric-water prediction, so the moist-state warm bias remains the
> primary unresolved equilibrium error.

**From `stratospheric-water.md`** (reference column, 30 layers, 3600 s,
200000 steps, last fifth):

| run | T_s |
|---|---:|
| K = 0, dry start | 291.941 K |
| K = 0.39 m² s⁻¹, dry start | 297.934 K |

K = 0.39 m² s⁻¹ is Massie and Hunten's tracer minimum, from Hu et al. (2012)
Table 3.

**This protocol changes no physics and adopts nothing.** It neither tunes
nor reopens the two-column work.

## R1: is the reference-plus-diffusion equilibrium vertically converged?

**Configuration.** `strat_run.py`'s reference column exactly:
`reference_column(grid=column_grid(16), areal_heat_capacity=1.0e6,
master_dt=3600, initial_water_store=1.0e7, eddy_diffusivity=0.39)`, from the
dry start, 200000 steps, read over the last fifth. Only `n_layers` changes.

**Resolutions, frozen now:**

- **30,** the calibration resolution;
- **40,** the established grid;
- **60,** the finest grid this project has already run, in the full water
  cycle's resolution study.

The vapour bands are calibrated at 30 layers (`CALIBRATION_N_LAYERS`). Running
at 40 and 60 is a resolution study under VAL-002, not a recalibration.

**Reported, per resolution:**

- T_s, PW, OLR and ASR;
- RH(p) and T(p);
- the mean RH over layers centred at 160–230 hPa;
- the coldest tropospheric temperature and its pressure;
- stratospheric water (ppmv, layers below 100 hPa);
- the run's worst relative energy and water residual;
- the validity-domain outcome;
- the equilibrium test (`surface_temperature_bias`: imbalance ≤ 0.1 W m⁻²).

**Classes, frozen:**

| class | condition |
|---|---|
| **Invalid** | a run is refused, fails a mandatory gate, or does not equilibrate |
| **Vertically converged (at this criterion)** | across {30, 40, 60}: T_s spans ≤ 0.5 K, PW spans ≤ 5% of its 30-layer value, and the 160–230 hPa RH spans ≤ 0.05 |
| **Resolution trend** | not converged, and T_s is monotonic in the layer count |
| **Resolution scatter** | not converged, and T_s is not monotonic |

**The 0.5 K basis.** The reference column without diffusion scatters
+3.7 to +4.3 K across 25–45 layers, about ±0.3 K. That scatter follows the
convective top's layer edges. A span above 0.5 K exceeds what that known
scatter explains.

**+9.93 K stays labelled "30 layers" unless the converged class holds.**

## R2: offline radiative attribution of 291.9 → 297.9 K

**The states, from two fresh runs.**

- **A:** K = 0, dry start, 30 layers.
- **B:** K = 0.39, dry start, 30 layers. This is R1's 30-layer run.

Each saves **50 start-of-step snapshots, evenly spaced over its last fifth:**
T(p), w(p), T_s, layer pressures, and the published OLR and ASR.

**Reproduction:** A's T_s and B's T_s match `stratospheric-water.md`'s
291.941 K and 297.934 K to ≤ 0.0005 K (its rounding).

**Nothing is evolved.** One radiation module, built exactly as the
reference column builds it, is solved on substituted states. This is
`regime-cloud-attribution.md` C1's method.

**Identity check:** the unsubstituted solve reproduces each snapshot's
published OLR and ASR to ≤ 1e-9 relative. The two runs' layer pressures must
be identical.

**The temperature state is (T(p), T_s) together. The vapour state is w(p).**

**Cases, on each paired snapshot (A_i, B_i),** with the TOA net
N = ASR − OLR evaluated in each:

| case | T state | vapour |
|---|---|---|
| AA | T_A | q_A |
| AB | T_A | q_B |
| BA | T_B | q_A |
| BB | T_B | q_B |
| A + ΔT_s (Planck) | T_A + (T_s,B − T_s,A) at every layer and the surface | q_A |
| AB_r, each region r | T_A | q_A with region r replaced by q_B |
| BA_r, reciprocal | T_B | q_B with region r replaced by q_A |

**Regions, frozen by layer-centre pressure:**

- lower troposphere, p ≥ 700 hPa;
- middle troposphere, 700 > p ≥ 400 hPa;
- upper troposphere, 400 > p ≥ 100 hPa;
- stratosphere, p < 100 hPa.

**Quantities, means over the 50 pairs, for OLR, ASR and N:**

| quantity | definition |
|---|---|
| vapour effect | Δ_q = X(AB) − X(AA) |
| temperature effect | Δ_T = X(BA) − X(AA) |
| interaction (a result, not a residual to redistribute) | I = [X(BB) − X(AA)] − Δ_q − Δ_T |
| Planck | Δ_P = X(A + ΔT_s) − X(AA) |
| lapse rate | Δ_LR = Δ_T − Δ_P |
| regional vapour effects | Δ_q,r = X(AB_r) − X(AA), and the reciprocal X(BB) − X(BA_r) |
| regional additivity | Σ_r Δ_q,r against Δ_q |

**Why the adjustment.** A and B are both equilibrated, so N(AA) ≈ N(BB) ≈ 0,
and the total change is near zero. The vapour and temperature effects
largely cancel, and ranking parts of a near-zero total would mean nothing.

In N, a positive term drives warming at fixed temperature. Planck is the
restoring response to the warmer surface and atmosphere. So the
**warming-attributable terms** are Δ_q, Δ_LR and I, all in N, and Δ_P is
reported as the term that balances them. The split of Δ_T into Planck and
lapse rate is the standard kernel decomposition, done here by direct
substitution.

**Robustness.** The whole matrix is repeated with a shifted pairing (A_i,
B_(i+25 mod 50)), and on the two runs' time-mean states. The limit cycle's
phases are not matched between runs, so the pairing is arbitrary and the
spread is reported.

### Classes, frozen (descriptive; no outcome is preferred)

**On the warming-attributable terms in N (Δ_q, Δ_LR, I), primary pairing:**

| class | condition |
|---|---|
| **Invalid** | the identity, reproduction or equal-pressure check fails |
| **Vapour dominant** | \|Δ_q\| is the largest, and ≥ 2 × the next |
| **Lapse-rate dominant** | \|Δ_LR\| is the largest, and ≥ 2 × the next |
| **Interaction dominant** | \|I\| is the largest, and ≥ 2 × the next |
| **Mixed** | otherwise |

**If vapour dominant, the region,** on Δ_q,r with A as the base:

- **upper-tropospheric vapour dominant,** if it is the largest and ≥ 2 × the
  next;
- **lower/mid-tropospheric vapour dominant,** if the lower and middle
  regions combined are the largest and ≥ 2 × the next;
- **stratospheric vapour dominant,** likewise;
- otherwise **regionally mixed.**

If the reciprocal swap (B as the base) gives a different region class, the
region class is recorded as **direction-dependent**.

**Reported in addition:**

- **nonlinear interaction material,** if \|I\| ≥ 0.10 × \|Δ_P\|;
- **pairing-sensitive,** if the shifted pairing or the time-mean states
  change the primary class.

**Never claimed from R2:** which dynamical operator made T_B differ from T_A.
The hierarchy is the vapour radiative effect, plus the temperature-profile
radiative effect, plus their interaction. Operator diagnostics may be used
afterwards, and labelled as such.

## What this cannot separate

- **Offline radiation at fixed states.** It explains the TOA balance of two
  equilibria, not the path between them.
- **Clouds** are whatever the reference column prescribes, the same in A and
  B.
- **The limit cycle.** Snapshot pairing is arbitrary; the robustness checks
  bound its effect.
- **One coefficient, K = 0.39.** Its bracket (K = 3.9) was already found
  fragile in `stratospheric-water.md`.
- **A single column with prescribed physics choices,** at a calibration
  made at 30 layers.

---

## Results (2026-10-21)

Runs: `scratchpad/r_run.py`, 200000 steps each, dry start, read over the last
fifth. Attribution: `scratchpad/r2.py` → `r2.json`. **Every run completes,
passes the mandatory gates (no invalidation reasons), equilibrates, and closes
energy and water to about 1e-15.**

### R1: resolution check (K = 0.39)

| layers | T_s | error | imbalance, W m⁻² | PW, kg m⁻² | RH 160–230 hPa | coldest tropospheric T | strat. water | run time |
|---:|---:|---:|---:|---:|---:|---|---:|---:|
| 30 | 297.934 K | +9.93 K | −0.003 | 69.26 | 0.907 | 208.4 K at 127 hPa | 17.1 ppmv | 555 s |
| 40 | 298.070 K | +10.07 K | +0.004 | 69.49 | 0.891 | 207.1 K at 123 hPa | 17.0 ppmv | 609 s |
| 60 | 298.391 K | +10.39 K | −0.001 | 70.65 | 0.876 | 203.2 K at 102 hPa | 15.2 ppmv | 745 s |

| criterion | span | line | |
|---|---:|---:|---|
| T_s | 0.457 K | 0.5 K | holds |
| PW | 2.0% | 5% | holds |
| RH 160–230 hPa | 0.032 | 0.05 | holds |

**Frozen class: vertically converged (at this criterion).** The headline
with diffusion is +9.9 to +10.4 K across 30–60 layers.

**Recorded as found:** T_s rises monotonically with the layer count (+0.14 K
from 30 to 40, +0.32 K from 40 to 60). The span is inside the frozen 0.5 K
line, but its shape is a trend rather than a scatter. The stratospheric
water still fails against SWOOSH's 4.0 ppmv at every resolution.

### R2: radiative attribution (A: K = 0 → B: K = 0.39, 30 layers)

**Checks:**

| check | result |
|---|---|
| reproduction | A: T_s 291.9414 K against 291.941; B: 297.9341 K against 297.934. **Pass** |
| equal pressures | **pass** |
| identity | the offline solve reproduces all 100 snapshots' published OLR and reflected shortwave **exactly** (relative 0) |

The radiation's other inputs, including `saturated_layer_temperature` and
`vapour_bearing_layer_temperature`, were held at the solver kernel's own
values. The exact identity across 100 different states shows they do not
affect the fluxes.

**Regions, by layer index:** lower 0–8; middle 9–17; upper 18–26
(394.6–127.0 hPa); stratosphere 27–29.

**TOA net N = ASR − OLR, W m⁻², means over 50 pairs (primary pairing):**

| term | N | OLR | ASR |
|---|---:|---:|---:|
| total, BB − AA | +0.17 | −0.18 | −0.01 |
| **vapour, Δ_q** | **+27.95** | −27.95 | 0.00 |
| temperature, Δ_T | −32.74 | +32.72 | −0.01 |
| Planck, Δ_P | −24.12 | +24.11 | −0.01 |
| **lapse rate, Δ_LR** | **−8.61** | +8.61 | 0.00 |
| **interaction, I** | **+4.96** | −4.96 | 0.00 |

**The reference column's shortwave does not see the vapour** (Δ_q in ASR is
0), so the attribution is entirely longwave.

**Regional vapour effects in N:**

| base | lower | middle | upper | stratosphere | sum | Δ_q |
|---|---:|---:|---:|---:|---:|---:|
| A (q_B swapped in) | +1.62 | +11.02 | **+17.66** | +1.64 | 31.93 | 27.95 |
| B (q_A swapped out) | +0.73 | +7.51 | **+16.07** | +2.28 | — | — |

The regional sum exceeds Δ_q by 4.0 W m⁻². That is the overlap
non-additivity, recorded and not redistributed.

**Robustness:**

| | Δ_q | Δ_LR | I |
|---|---:|---:|---:|
| shifted pairing | +27.92 | −8.59 | +4.96 |
| time-mean states | +27.98 | −8.55 | +4.83 |

All classes are unchanged in both, so the result is **not pairing-sensitive**.

### Classification, as frozen

- **Not invalid.**
- **Primary: vapour dominant.** \|Δ_q\| = 27.95 is the largest, and ≥ 2 ×
  \|Δ_LR\| = 8.61.
- **Region: regionally mixed, from both bases.**
  - With A as the base, upper (17.66) is the largest, but not ≥ 2 × lower
    plus middle (12.64).
  - With B as the base, upper (16.07) against lower plus middle (8.24)
    misses the 2× line by 0.4 W m⁻².
- **Nonlinear interaction: material.** \|I\| = 4.96 ≥ 0.10 × \|Δ_P\| = 2.41.

**What the numbers say, beside the class.**

- **The 6 K comes from the water vapour diffusion carries aloft.** At fixed
  temperature it traps 28 W m⁻².
  - The upper troposphere (395–127 hPa) supplies the largest share: 18 of
    32 W m⁻² regionally from A's base, 16 from B's.
  - The middle troposphere supplies most of the rest: 11 and 8.
- **The temperature-profile change opposes it,** by 8.6 W m⁻² beyond
  Planck.
- **Planck balances the remainder,** 24.1 W m⁻² for a 6.0 K warming.
- **Which operator made B's temperature profile is not claimed here.**

## Disposition (the user's, 2026-10-21), recorded verbatim

**The interpretation, tightened:**

> R2 identifies the equilibrium separation as vapour-dominant radiatively:
> the humidity change contributes +27.9 W m⁻² at fixed temperature, opposed
> by −8.6 W m⁻² from the nonuniform temperature-profile change, with a
> material +5.0 W m⁻² interaction.

**The headline after R1 and R2:**

> The diffusion solution is numerically reproducible and vertically
> converged, but not yet physically validated. R2 shows that its ~6 K
> warm-state separation is radiatively vapour-dominant and regionally mixed,
> with the largest contribution from the upper troposphere. The remaining
> decisive failure is excessive stratospheric water.

**What R2 does not establish:**

- why diffusion put the vapour there;
- whether the excess resides in physically realistic pathways;
- whether the temperature profile is dynamically realizable;
- that diffusion should become the default operator.

The failed stratospheric-water criterion blocks the inference "diffusion
removes bistability, converges, therefore is physically superior". The
four-way balance closes: 27.9 − 8.6 + 5.0 − 24.1 ≈ +0.2 W m⁻². The 4 W m⁻²
regional overcount is kept explicit, as the interaction term is.

---

# R3: is the diffusion solution's free troposphere too humid against an independent reference? (frozen 2026-10-21, before any R3 code or parse)

**The question is the user's,** and frozen without a mechanism:

> R3 tests whether the diffusion solution's free troposphere is
> systematically too humid relative to an independent reference, particularly
> in the pressure regions R2 identifies as carrying most of the vapour
> radiative effect.

**What R3 does not test:**

- **it does not test subsidence drying.** That causal test follows R3 and is
  not built into it;
- **it does not test stratospheric water.**

**Two failures are kept separate:**

- **S1, stratospheric validity.** Excess water above the tropopause.
  `stratospheric-water.md` points to the missing temperature inversion. It
  is important for physical validity, and contributes little to R2's warming
  (+1.6 to +2.3 W m⁻² of the vapour effect).
- **T1, the tropospheric humidity question.** The near-saturated middle and
  upper free troposphere, which R2 identifies as carrying most of the vapour
  effect. Its cause is to be demonstrated.

**The user's decisions:**

1. **Staged:** R3 now, over 1000–300 hPa against NCEP R1. **R3b**, over
   300–100 hPa, follows once a satellite humidity record (AIRS v7 Level 3 or
   Aura MLS) is registered. R3b's design is frozen then, not now.
2. **Global mean primary,** with the tropical ±20° band reported beside it.
3. **"Systematically too humid":** B exceeds NCEP by ≥ +0.10 RH in both the
   middle and the upper bands.

## Data

**The model:** the late-mean RH(p) published by A (K = 0) and B
(K = 0.39), 30 layers, from R1/R2's runs (`rA30.json`, `rB30.json`). B at 40
and 60 layers is reported beside it.

**The reference:** NCEP R1 1991–2020 monthly long-term means of shum and air
(registered, checksums verified before parsing). Its humidity levels are
1000–300 hPa only.

- RH = q / q_s(T, p), with the model's saturation function: over ice below
  the triple point, as the model's own RH is.
- Formed for each month and grid point, then averaged over the year.
- Composited cos-latitude weighted:
  - **global (primary);**
  - **±20° (reported).**
- Provenance: reanalysis.

**Matching.** Model RH is log-p interpolated to NCEP's levels, with no
extrapolation. The model's lowest layer is centred at 996.5 hPa, so the
1000 hPa level is not compared. **The compared levels are 925–300 hPa.**

## Bands (R2's, cut at the reference's coverage)

| band | NCEP levels | in R2 terms |
|---|---|---|
| lower | 925, 850, 700 | R2's lower, p ≥ 700 |
| middle | 600, 500, 400 | R2's middle, 700 > p ≥ 400 |
| upper, to the reference limit | 300 | the bottom of R2's upper, 400 > p ≥ 100 |

The band value is the mean over its levels.

**160–230 hPa is reported model-only,** marked as having no reference until
R3b.

## Classes, frozen (global composite)

| class | condition |
|---|---|
| **Invalid** | a checksum fails, or a model profile does not span the compared levels |
| **Free troposphere systematically too humid** | B − NCEP ≥ +0.10 in both the middle and the upper bands |
| **Too humid in one band only** | ≥ +0.10 in exactly one of them |
| **Not systematically too humid** | otherwise |

**Reported always:**

- A − NCEP and B − NCEP in every band, global and ±20°;
- the per-level differences;
- B at 40 and 60 layers;
- the 160–230 hPa model values, with no reference.

## What R3 cannot separate

- **Reanalysis humidity is model-dependent** where there are few
  observations. R1's upper-tropospheric humidity is known to be weakly
  constrained.
- **The global column is a single profile against an area mean of
  profiles,** and the RH of a mean is not the mean of RHs. Here RH is formed
  per point and month, then averaged.
- **Below-ground levels over high terrain.** R1 extrapolates there; they
  are not masked.
- **Monthly means:** q / q_s of monthly means is not the monthly mean of
  q / q_s.

## R3 results (2026-10-21)

`scratchpad/r3.py` → `r3.json`. The NCEP files' checksums were verified
before parsing, and every model profile spans 925–300 hPa.

**RH by level, global composite (NCEP R1) and the 30-layer runs:**

| level, hPa | NCEP global | NCEP ±20° | A (K = 0) | B (K = 0.39) |
|---:|---:|---:|---:|---:|
| 925 | 0.758 | 0.725 | 0.735 | 0.759 |
| 850 | 0.641 | 0.625 | 0.851 | 0.847 |
| 700 | 0.531 | 0.460 | 0.777 | 0.879 |
| 600 | 0.502 | 0.393 | 0.738 | 0.764 |
| 500 | 0.516 | 0.362 | 0.728 | 0.835 |
| 400 | 0.528 | 0.354 | 0.748 | 0.765 |
| 300 | 0.633 | 0.593 | 0.609 | 0.829 |

**Band means, and the model minus NCEP:**

| band | NCEP global | A | A − NCEP | B | **B − NCEP** | B − NCEP, ±20° |
|---|---:|---:|---:|---:|---:|---:|
| lower (925–700) | 0.643 | 0.788 | +0.144 | 0.828 | +0.185 | +0.225 |
| middle (600–400) | 0.515 | 0.738 | +0.223 | 0.788 | **+0.273** | +0.419 |
| upper, to 300 | 0.633 | 0.609 | −0.024 | 0.829 | **+0.196** | +0.236 |

**B at other resolutions,** global:

| layers | middle | upper, to 300 |
|---:|---:|---:|
| 40 | +0.276 | +0.181 |
| 60 | +0.295 | +0.154 |

**Model-only, with no reference until R3b:** the mean RH at 160–230 hPa is
0.000 in A, and 0.907 / 0.891 / 0.876 in B at 30 / 40 / 60 layers.

### Classification, as frozen

- **Not invalid.**
- **B − NCEP is +0.273 in the middle band and +0.196 in the upper band,**
  both ≥ +0.10.
- **Frozen verdict: free troposphere systematically too humid.** It holds at
  40 and 60 layers as well.

### Recorded as found, beside the class

- **The middle-tropospheric excess is not diffusion's.** A, with no
  diffusion, is already +0.22 too humid there. B adds +0.05.
- **Diffusion's added humidity is concentrated higher up:**
  - +0.22 at 300 hPa (A 0.61 against B 0.83);
  - 0 → 0.91 at 160–230 hPa.
  - A's 300 hPa level is near NCEP's (−0.02).
- **Two parts of T1, then:**
  - **a mid-tropospheric excess common to A and B,** about +0.2 to +0.3
    globally and +0.4 in the ±20° band;
  - **an upper-tropospheric excess that diffusion creates.**
- **NCEP's 300 hPa RH is high (0.63)** relative to its 400 hPa value. That
  is consistent with RH over ice above the freezing level, and with R1's
  weakly constrained upper-tropospheric humidity. R3b's satellite record
  is the check on it.
- **The cause of either part is not claimed here.**

## R3 disposition (the user's, 2026-10-21), recorded verbatim

> R3 finds a systematic free-tropospheric moist bias, but it is not produced
> uniformly by diffusion. The middle-tropospheric bias is already present in
> the no-diffusion state A, while diffusion primarily introduces an
> additional upper-tropospheric moistening beginning near 300 hPa and
> extending above the NCEP R1 humidity ceiling.

**The two defects stay separate until the evidence says otherwise:**

- **a base-state mid-tropospheric moist bias:** at 600–400 hPa, A is
  already +0.22 RH against NCEP;
- **a diffusion-specific upper-tropospheric moistening:** at 300 hPa, A is
  on the reference while B is +0.22 higher.

The region R3 does not measure, 160–230 hPa, is where B reaches about 0.91
and where R2's largest regional radiative contribution lies.

**The two-column subsidence experiment is held until R3b.** R3b's primary
reference is Aura MLS, defined narrowly:

> Does an independent satellite record confirm or reject the diffusion
> state's very high upper-tropospheric humidity in the region above NCEP R1's
> usable humidity coverage?

---

# R3b: does Aura MLS confirm or reject B's upper-tropospheric humidity? (frozen 2026-10-21, before any download)

> Does an independent satellite record confirm or reject the diffusion
> state's very high upper-tropospheric humidity in the region above NCEP R1's
> usable humidity coverage?

**What R3b does not test:** it makes no causal claim, about subsidence or
anything else, and tunes nothing.

## The reference

**MLS/Aura Level 3 Monthly Binned RH with respect to ice, V005
(ML3MBRHI),** from NASA GES DISC.

- **Coverage:** 82°S–82°N at 4° × 5°.
- **Recommended range:** 316 hPa upward, per the product's catalog entry.
  It points to section 3.20 and chapter 4 of the EOS MLS Level 2 V5 Quality
  Document.
- **Registration:** every file is registered in `data/sources.yaml` with its
  SHA-256 before any parse, as `provenance_class: observational`.

**Screening.** Before any value is parsed, the product documentation's
statement of what screening the Level 3 production already applies is
recorded here, verbatim.

- If it prescribes user-side screening of the Level 3 files (a status,
  count or quality field), exactly that is applied.
- If the Level 3 product cannot be shown to be built from screened Level 2
  data, R3b falls back to Level 2 ML2RHI with the quality document's
  section 3.20 rules, by erratum.
- **No value is discarded for being negative, noisy or surprising.** Only
  the product's own fill value is excluded.

**The convention.** The model's published RH is formed with its own
saturation function, which is over ice below the triple point. So above the
freezing level it is RH over ice, the convention ML3MBRHI uses. The
saturation formulas differ slightly, and that is declared, not converted.

## The comparison

- **Period:** 2005-01 to 2020-12. Full years after MLS's August 2004
  start, ending with NCEP's climatology period.
- **Averaging, per bin:**
  1. the mean over the years available for each calendar month;
  2. the mean of the 12 months.

  Then a cos-latitude-weighted area mean: **global (82°S–82°N), primary;**
  ±20° reported. A bin missing in a calendar month is left out of that
  month, with weights renormalised. The fraction left out is reported.
- **Levels: the product's own pressure levels with 100 ≤ p ≤ 316 hPa.** The
  model is log-p interpolated onto them, with no extrapolation. This
  compares on the satellite's levels, not by putting the reference onto the
  model's.
- **Model:** A and B late-mean RH at 30 layers (`rA30.json`, `rB30.json`).
  B at 40 and 60 layers is reported.

**The primary band:** the mean over the product levels within
147 ≤ p ≤ 215 hPa. It is the satellite-level counterpart of the model's
160–230 hPa diagnostic.

**A consistency check, reported, not classified:** MLS at 316 hPa against
NCEP R1 at 300 hPa (R3), both against the model.

## Classes, frozen (global, primary band)

| class | condition |
|---|---|
| **Invalid** | a file is unregistered or fails its checksum; a month of 2005–2020 is missing; or a model profile does not span the compared levels |
| **Diffusion state contradicted aloft** | B − MLS ≥ +0.10 |
| **Diffusion state consistent aloft** | \|B − MLS\| < 0.10 |
| **Diffusion state too dry aloft** | B − MLS ≤ −0.10 |

**Reported always:**

- A − MLS;
- whether the pattern "A within 0.10 of MLS, and B ≥ +0.10" holds;
- every level, global and ±20°;
- B at 40 and 60 layers.

## What R3b cannot separate

- **MLS's vertical resolution is 3–6 km.** The model's sharp profile is
  compared with a smoothed retrieval, and Level 3 carries no averaging
  kernels.
- **The two saturation formulas over ice differ slightly.**
- **A global equilibrium column against a 2005–2020 climatology** of
  82°S–82°N.
- **Retrieval sensitivity** to thick cloud in the upper troposphere, as the
  quality document describes it.

## R3b erratum 1: the frozen product cannot cover the period (2026-10-26, before any MLS value was parsed; nothing above this line changed)

**The finding.** ML3MBRHI V005 has no files for 2007, 2008, 2010, 2011 or
2014–2017.

- **Sources checked:** the retired on-premises archive's year listing, and
  NASA's Common Metadata Repository (CMR). CMR returns 14 granules for
  version 005: 2004–2006, 2009, 2012, 2013 and 2018–2025.
- **The Earthdata Cloud archive** has the same gaps.
- **Nothing was parsed to establish this.** Only directory listings and
  catalogue records were read.

**Both readings, kept:**

1. **As frozen:** a month of 2005–2020 is missing, so R3b is **Invalid**.
   No constructible observational reference.
2. **By this erratum:** R3b takes the fallback the frozen text already names,
   **Level 2 ML2RHI V005 with the quality document's section 3.20 rules.**
   This is the user's choice (2026-10-26), from four offered:
   - Level 2 V005;
   - Level 3 V004 for 2005–2020;
   - Level 3 V005 for the 8 years that exist;
   - closing R3b.

### The Level 3 screening statement the frozen text asked for, verbatim

The source is the EOS MLS Level 2/3 Version 5.0x Quality document,
`Aura-MLS_DataQuality_v5-0x-revA.pdf`, chapter 4 (document pages 160 and
163). Spaces and digits are restored from two text extractions of the PDF;
no word was changed.

> As part of the generation of these data products, all the data screening
> rules described in chapters 1 and 3 are applied (along with the ClO bias
> correction, see section 3.6), providing a significant additional
> simplification to data users.

> As discussed above, these Level 3 products are generated using only the
> Level 2 data points that meet all the screening criteria described in the
> earlier chapters of this document, providing a significant simplification
> for users. In cases where there are no good Level 2 data available to
> populate a Level 3 data point, the nvalues field is set to zero and the
> other data fields are set to the _FillValue (bad data flag) for the given
> netCDF variables.

So the Level 3 product is built from screened Level 2 data. The fallback is
taken because of the missing years, not because of screening.

### The Level 2 screening statement, verbatim

**Section 3.20.6, "Data screening" (document pages 138–139):**

> **Pressure range:** Profile from 316–0.001 hPa, larger pressures represent
> mid/upper troposphere column. Values outside this range are not
> recommended for scientific use.
>
> **Estimated precision:** Only use values for which the estimated precision
> is a positive number. Values where the a priori information has a strong
> influence are flagged with negative or zero precision, and should not be
> used in scientific analyses (see Section 1.5).
>
> **Status flag:** Only use profiles for which the Status field is an even
> number. Odd values of Status indicate that the profile should not be used
> in scientific studies. See Section 1.6 for more information on the
> interpretation of the Status field.
>
> **Clouds:** Ignore status bits 16 (high clouds) or 32 (low clouds) set
> indicating the presence of clouds. See artifacts for more details.
>
> **Quality field:** The Quality fields for both the L2GP-RHI and
> L2GP-Temperature swaths need to be considered, as described below:
> For pressures of 83 hPa and smaller: Use only profiles with RHI Quality
> greater than 0.7 and Temperature Quality greater than 0.2
> For pressures of 100 hPa and larger: Use only profiles with RHI Quality
> greater than 0.7 and Temperature Quality greater than 0.9.
>
> **Convergence field:** Only profiles with a value of the RHI Convergence
> less than 2.0 and Temperature Convergence less than 1.03 should be used in
> scientific studies.
>
> **Temperature precision:** The L2gpPrecision field in the
> L2GP-Temperature file can be used to further eliminate outliers that are
> believed to be the result of thick clouds, primarily in the tropics. If
> careful screening of the troposphere is required, levels 261–178 hPa
> should be avoided if any of the following criteria are met:
> At 316 hPa: L2gpPrecision > 1.1 K and latitude > −60°
> At 261 hPa: L2gpPrecision > 0.7 K
> At 215 hPa: L2gpPrecision > 0.825 K
>
> **Additional screening for based on H2O:** Eliminate highly oscillatory
> profiles by eliminating H2O profiles having values less than 0.101 ppmv at
> any pressure < 1 hPa (vertical levels 1–37).

In "for based on", the PDF text carries no word between "for" and "based";
it is quoted as it stands.

**Section 3.9.9 (H2O), the same rule as the parent product states it:**

> H2O profiles having values less than 0.101 ppmv (0.101×10⁻⁶ VMR) at any
> pressure at or larger than 1 hPa (vertical levels 1–37) should not be
> used.

**Section 1.5:**

> Data points for which L2gpPrecision is set negative or zero should not be
> used, as this flags that the resulting precision is worse than 50% of the
> a priori precision, indicating that instrument and/or the algorithms have
> failed to provide enough useful information for that point.

### The rules applied, frozen

R3b's levels are 316, 261, 215, 178, 147, 121 and 100 hPa. All of them fall
under the "100 hPa and larger" Quality rule.

| # | rule | applies to |
|---|---|---|
| 1 | RHI `L2gpPrecision` > 0 | each value |
| 2 | RHI `Status` even. Bits 16 and 32 alone reject nothing. | each profile |
| 3 | RHI `Quality` > 0.7 and Temperature `Quality` > 0.9 | each profile |
| 4 | RHI `Convergence` < 2.0 and Temperature `Convergence` < 1.03 | each profile |
| 5 | Reject if H2O `L2gpValue` < 0.101 × 10⁻⁶ at any level with p ≥ 1 hPa | each profile |
| 6 | Exclude the product's own fill value | each value |

**Rule 5 has two readings, and both are kept.**

- **Primary: p ≥ 1 hPa.** This is section 3.9.9's statement. The level
  numbers 1–37, which both sections quote, run from 1000 to 1 hPa on the
  12-per-decade grid.
- **Reported beside it: p < 1 hPa.** This is section 3.20.6 as written.
- The user's choice, 2026-10-26.

**The conditional temperature-precision rule:**

- It is **off for the class.** Its own condition is "if careful screening is
  required", and the frozen rule is "no value is discarded for being noisy".
- With it on, the result is **reported beside the class.** A profile meeting
  any of its three criteria loses its values at 261, 215 and 178 hPa, and
  the fraction removed is reported.
- The user's choice, 2026-10-26.

**Matching profiles across products.** Temperature and H2O are in their own
daily files, ML2T and ML2H2O V005.

- Each RHI profile is matched to the Temperature and H2O profiles in the
  same day's files with an identical `Time` value.
- An RHI profile without both matches is excluded, and the count is
  reported.
- Nothing else is used from Temperature or H2O.

**Version mixing.** The files are taken as distributed and as listed by CMR.

- Within a year they mix software versions v05-00, v05-01 and v05-02, and
  collections c01 to c03.
- Each file's version is recorded in its manifest line.
- Nothing is selected by version.

### The climatology, frozen

The method mirrors the Level 3 product's simple binning (chapter 4: "simply
reporting the mean measured value within a given spatial/temporal box").

- **Bins: 4° latitude × 5° longitude.**
  - Latitude edges run −82, −78, …, 82, giving 41 bins.
  - Longitude edges run −180, −175, …, 180, giving 72 bins.
  - Profiles with latitude beyond ±82° are not used.
- **The monthly bin value** at each level is the mean of the screened values
  in that bin in that calendar month of that year.
- **Units:** RHI's `L2gpValue` is in percent and is divided by 100. Nothing
  else is changed.
- **From there, the frozen averaging stands unchanged:**
  1. the per-bin calendar-month mean over the years with a value;
  2. the mean of the 12 months;
  3. the cos-latitude area mean, global primary and ±20° reported, with
     missing bins left out and weights renormalised. The fraction left out
     is reported.
- **The frozen completeness gate stands.** R3b is Invalid if any calendar
  month of 2005–2020 has no RHI file for any of its days.
- **Partial months are used as they are.** The number of days present in
  each month is reported.

The levels, the primary band (147–215 hPa: 215, 178 and 147 hPa), the model
runs, the interpolation and the classes are unchanged.

### Registration, frozen

Per-file entries for 17,340 daily files would put about 260,000 lines into
`data/sources.yaml`. Each product is therefore registered through a manifest.

- **The manifest format:** one line per file, `<sha256>  <filename>`, sorted
  by filename, UTF-8, with LF line endings.
- **Location:** each manifest is written to `data/raw/aura_mls/` from the
  downloaded files.
- **Registration:** each manifest is registered in `data/sources.yaml` as
  `provenance_class: observational`, with its own SHA-256.
  - Because the format is fixed, the manifest is reproducible from the files
    alone.
  - So the recorded hash pins every file.
- **The products:** ML2RHI, ML2T and ML2H2O V005, 2005-001 to 2020-366.
- **The quality document** (`Aura-MLS_DataQuality_v5-0x-revA.pdf`) is
  registered as its own entry.
- **The order of checks, before any value is read:**
  1. the analysis verifies each manifest against `data/sources.yaml`;
  2. it verifies each file against its manifest line, before opening it.
- **Invalid also covers a file that is unregistered or fails its checksum.**

### What R3b cannot separate, added (the list above is not shortened)

- **The climatology is this project's binning of Level 2, not the MLS team's
  Level 3.** The rules are the same; the code is not.
- **Section 3.20.8:** "MLS produces higher relative humidities in the moist
  regions of the tropics relative to AIRS."
- **Section 3.20.7:**
  - A possible drift in lower-stratospheric water vapour, inherited by RHI.
  - Clouds "significantly degrade the precision ... but the mean bias as
    compared to AIRS changes by less than 10%."
  - Cloud-flagged profiles are kept, as the document directs.

## R3b results (2026-10-26)

**Run validity: valid by the frozen gate.**

- **Registration:** the three manifests and the quality document are
  registered (`ec59424`). Each was verified before use, and each of the
  17,340 daily files was verified against its manifest line before it was
  opened.
- **Completeness:** no calendar month of 2005–2020 is missing. The fewest
  days in a month is 12 (2011-04 and 2018-06); the most is 31.
- **Spans:** every model profile spans the compared levels.
- **Profiles:** 19,970,110 RHI profiles. 2,167 had no Temperature or H2O
  profile at the identical `Time` and were excluded. None lay beyond ±82°.
- **The script:** `r3b.py`, in this session's scratchpad, as for R3.

**Screening, as a fraction of the 19,967,943 matched profiles:**

| rule | fails |
|---|---:|
| Status odd | 0.29% |
| Quality (RHI > 0.7, T > 0.9) | 2.16% |
| Convergence (RHI < 2.0, T < 1.03) | 0.53% |
| H2O, p ≥ 1 hPa (primary) | 2.11% |
| H2O, p < 1 hPa (the other reading) | 79.1% |
| value precision ≤ 0 or fill (fraction of values) | 0.68% |
| thick-cloud rule (reported variant), profiles flagged | 2.97% |

95.4% of the values considered enter the primary climatology.

**Bins left out:**

- **Primary and cloud-rule variants:** no bin-month at any level.
- **The p < 1 hPa reading:** 11% (global) and 13% (±20°) of bin-months,
  because that reading rejects four profiles in five.

**Levels.** The nominal names are those of the erratum. The product's levels
are 316.23, 261.02, 215.44, 177.83, 146.78, 121.15 and 100.00 hPa. The
model is interpolated onto these.

### MLS and the model, level by level (RH as a fraction)

**Global (82°S–82°N):**

| hPa | MLS primary | MLS cloud rule on | MLS H2O p<1 | A (30) | B (30) | B (40) | B (60) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 316 | 63.804 | 63.804 | 322.591 | 0.687 | 0.808 | 0.806 | 0.788 |
| 261 | 3.483 | 0.429 | 31.799 | 0.981 | 0.839 | 0.835 | 0.830 |
| 215 | 1.176 | 0.489 | 0.602 | 0.000 | 0.877 | 0.861 | 0.864 |
| 178 | 0.462 | 0.402 | 0.413 | 0.000 | 0.932 | 0.927 | 0.875 |
| 147 | 0.351 | 0.351 | 0.345 | 0.000 | 0.968 | 0.978 | 0.988 |
| 121 | 0.320 | 0.320 | 0.316 | 0.000 | 0.999 | 1.000 | 1.000 |
| 100 | 0.360 | 0.360 | 0.369 | 0.000 | 1.000 | 1.000 | 1.000 |

**±20° (bins whose centres lie within ±20°):**

| hPa | MLS primary | MLS cloud rule on | MLS H2O p<1 |
|---|---:|---:|---:|
| 316 | 142.719 | 142.719 | 402.211 |
| 261 | 7.908 | 0.395 | 1.014 |
| 215 | 2.434 | 0.670 | 0.907 |
| 178 | 0.842 | 0.693 | 0.741 |
| 147 | 0.699 | 0.699 | 0.691 |
| 121 | 0.654 | 0.654 | 0.649 |
| 100 | 0.747 | 0.747 | 0.772 |

### The primary band (215, 178 and 147 hPa)

| | MLS | A − MLS | B30 − MLS | B40 − MLS | B60 − MLS |
|---|---:|---:|---:|---:|---:|
| **global, primary** | 0.663 | −0.663 | **+0.263** | +0.260 | +0.246 |
| global, cloud rule on | 0.414 | −0.414 | +0.512 | +0.508 | +0.495 |
| global, H2O p < 1 hPa | 0.453 | −0.453 | +0.472 | +0.469 | +0.456 |
| ±20°, primary | 1.325 | −1.325 | −0.400 | −0.403 | −0.416 |
| ±20°, cloud rule on | 0.687 | −0.687 | +0.238 | +0.235 | +0.221 |
| ±20°, H2O p < 1 hPa | 0.780 | −0.780 | +0.146 | +0.142 | +0.129 |

### Classification, as frozen

- **B − MLS = +0.263 in the global primary band, which is ≥ +0.10. So the
  class is "Diffusion state contradicted aloft."**
  - The same class holds at 40 and 60 layers (+0.260, +0.246).
  - It also holds under both reported variants (+0.512, +0.472).
- **The pattern "A within 0.10 of MLS, and B ≥ +0.10" does not hold.**
  - A is −0.663: A has zero RH from 215 hPa up.
  - So the observation lies between the two states (A ≪ MLS ≪ B), not
    with A.

### Recorded as found, beside the class

1. **The frozen means are dominated by a tail that passes every mandatory
   rule.**
   - The global primary MLS mean is 6380% at 316 hPa, 348% at 261 hPa and
     118% at 215 hPa.
   - A diagnostic over every 40th day (145 days, the same verification
     path, the primary rules) shows physical medians and an unphysical tail.

   | hPa | median | 99th pct | 99.9th pct | max | mean | mean, cloud-flagged excluded |
   |---|---:|---:|---:|---:|---:|---:|
   | 316 | 38% | 336% | 954% | 6.8×10⁷% | 229% | 59% |
   | 261 | 31% | 204% | 712% | 9.6×10⁵% | 52% | 44% |
   | 215 | 29% | 185% | 673% | 5.7×10⁶% | 55% | 40% |
   | 178 | 13% | 172% | 380% | 3.2×10⁵% | 33% | 31% |
   | 147 | 6% | 157% | 230% | 2.1×10⁴% | 26% | 25% |
   | 121 | 6% | 131% | 168% | 3.2×10³% | 23% | 23% |
   | 100 | 7% | 132% | 177% | 698% | 25% | 24% |

   - Of the values above 1000%, the thick-cloud rule of section 3.20.6
     flags 70% at 316 hPa, 79% at 261 and 215 hPa, 85% at 178 hPa, and all
     at 147 and 121 hPa.
   - This is the outlier population the document attributes to thick
     clouds, "primarily in the tropics".
   - The medians and cloud-excluded means are diagnostics, not frozen
     quantities, and classify nothing.
2. **The tail can only work against the class.**
   - Every tail value is large and positive, so it raises the MLS mean and
     lowers B − MLS.
   - The class is "contradicted" with the tail in. Any treatment that
     removed it would raise B − MLS, as the cloud-rule variant does (+0.263
     → +0.512).
3. **The ±20° primary band reads −0.40 only because of 215 hPa,** where the
   tropical MLS mean is 243%.
   - At 178 and 147 hPa, ±20° B − MLS is +0.09 and +0.27 (primary).
   - With the cloud rule on, the tropical band is +0.238.
4. **At the two levels the tail barely reaches** (178 and 147 hPa, where the
   primary and cloud-rule MLS means agree to 0.06 and exactly):
   - B30 − MLS is +0.47 and +0.62 globally;
   - B at 147 hPa is 0.97 against MLS 0.35.
5. **A has zero RH at and above 215 hPa** (dry start, no vapour transport,
   as R3 found for 160–230 hPa).
   - Its 261 hPa value is 0.98.
   - A's upper troposphere is too dry by the whole observed value. B's is
     too moist by 0.26 to 0.68, depending on the level and treatment, in the
     global band and the four levels above 215 hPa.
   - Neither state is consistent with MLS aloft.
6. **The 316 hPa consistency check (reported, not classified) is unusable as
   a mean.**
   - The MLS 316 hPa mean is 63.8 (6380%) globally and 142.7 at ±20°: the
     tail, at the bottom of the document's useful range (316–0.001 hPa).
   - Against NCEP R1 at 300 hPa (0.633 global, 0.593 at ±20°), A 0.687 and
     B 0.808.
   - The sample median at 316 hPa (38%) and cloud-excluded mean (59%) are
     recorded beside it.
7. **The H2O rule's other reading (p < 1 hPa) rejects 79% of profiles.**
   Above about 1 hPa, water vapour below 0.101 ppmv is the mesosphere's
   ordinary state rather than an oscillation. That reading is recorded and
   not adopted, as frozen.

**What this does and does not say.**

- B's upper troposphere is far more humid than MLS observes, under every
  frozen and reported treatment.
- A's is far drier. The observation sits between them.
- **No causal claim is made.** In particular, nothing here says what removes
  B's excess or supplies A's deficit.

## R3b disposition (the user's, 2026-10-26), recorded verbatim

**Transcription note.** Three formulas arrived doubled by the paste, once as
rendered text and once as source. Each is given once here; no other
character is changed.

- "approaching 107%10^7\%" is written "approaching 10⁷%".
- "p<1p<1" is written "p < 1".
- "p≥1p\ge1" is written "p ≥ 1".

Consecutive one-line items that arrived as separate lines are shown as
lists, and headings in bold; the words are unchanged.

> **Decision**
>
> Accept R3b as Valid and retain the frozen classification: the diffusion
> state is Contradicted aloft.
>
> Do not redesign the test or replace the mean after seeing the MLS values.
> The outlier tail is serious, but it does not threaten this classification.
> It biases the result against contradiction by increasing the MLS mean and
> reducing the apparent positive B-minus-MLS discrepancy.
>
> The registered result should therefore remain:
>
> R3b: Valid. Diffusion state B is Contradicted aloft.
> Over the global 147–215 hPa band, the frozen primary comparison gives
> B − MLS = +0.263, with corresponding 40- and 60-layer values of +0.260
> and +0.246. State A is also inconsistent with MLS, at −0.663. The
> observations fall between the two modeled states. No causal
> interpretation is assigned.
>
> **Why the result stands**
>
> The chain of custody is strong:
>
> - The erratum was frozen before observational values were read.
> - The manifest was registered separately.
> - No month from 2005 through 2020 was missing.
> - Every checksum was verified before its file was opened.
> - The primary and reported alternatives were specified before evaluation.
> - The layer-resolution variants agree closely: +0.263, +0.260, and
>   +0.246.
> - Activating the conditional cloud rule strengthens the discrepancy from
>   +0.263 to +0.512.
> - The alternate H2O interpretation also leaves B too humid at +0.472.
>
> This is not a marginal classification produced by a single aggregation or
> resolution choice.
>
> **Treatment of the outlier tail**
>
> Keep the mean-based primary design exactly as frozen, but place a
> conspicuous limitation beside the result.
>
> The correct distinction is:
>
> - Protocol validity: Valid.
> - Primary classification: Contradicted.
> - Estimator quality: The MLS arithmetic mean is badly distorted at some
>   levels.
> - Effect on classification: Conservative with respect to the
>   contradiction of B.
>
> I would use language such as:
>
> Estimator limitation: The screened MLS RHI distribution contains a small,
> extreme positive tail, including values approaching 10⁷% RH. These
> observations pass the frozen mandatory screening rules and are therefore
> retained in the primary calculation. The tail materially inflates some
> arithmetic means, including approximately 6,380% at 316 hPa and 118% at
> 215 hPa, compared with medians of 38% and 29%. Because the anomalous
> values are positive, their inclusion raises the MLS reference and reduces
> B − MLS. The primary classification is therefore not caused by the tail;
> applying the preregistered conditional cloud rule increases B − MLS from
> +0.263 to +0.512. No post hoc observations are removed from the primary
> result.
>
> That says precisely what matters without calling the values erroneous in
> a way the protocol has not demonstrated.
>
> **Do not reinterpret A as supported**
>
> The failure of the original expected pattern matters. R3b has not shown
> that A is observationally supported.
>
> The scientific conclusion is narrower:
>
> - A is much too dry, reaching zero RH from 215 hPa upward.
> - B is too humid over the registered upper-tropospheric band.
> - MLS falls between A and B.
> - Therefore, B is contradicted, but A is not validated.
> - R3b identifies an observational bracketing failure, not a successful
>   state selection.
>
> I would add this sentence prominently:
>
> Contradiction of B does not constitute confirmation of A; neither modeled
> state agrees with the MLS upper-tropospheric reference.
>
> **Treatment of the distorted secondary checks**
>
> ±20° band
>
> Retain the frozen reported value of −0.40, but mark it as tail-dominated
> and not representative of the vertical pattern.
>
> Report the diagnostic decomposition:
>
> - 215 hPa tropical mean: 243%, driving the negative band aggregate.
> - 178 hPa: +0.09
> - 147 hPa: +0.27
> - Conditional cloud-rule result: +0.24
>
> The −0.40 result should not be deleted, reclassified, or silently
> replaced. It should be labeled as an arithmetic-mean sensitivity to the
> 215 hPa tail.
>
> 316 hPa NCEP comparison
>
> Record it as:
>
> Computed but not interpretable as a mean-based comparison because the
> screened MLS arithmetic mean is dominated by the extreme positive tail.
>
> Do not call the entire R3b test invalid because one reported supporting
> comparison is uninterpretable. Keep it outside the evidence used for the
> primary class.
>
> H2O alternate reading
>
> The 79% rejection rate under p < 1 hPa is strong internal evidence that
> the selected p ≥ 1 hPa interpretation was the coherent reading. It should
> remain a reported interpretive check and should not alter the primary
> population.
>
> **Should the mean-based design stand?**
>
> For this frozen R3b execution: yes. For future protocols: no, not by
> itself.
>
> Changing the estimator now would be outcome-responsive. Instead:
>
> - Preserve the arithmetic mean as the registered primary estimator.
> - Report median, sample count, relevant quantiles, and tail diagnostics
>   descriptively.
> - Do not let any new robust estimator change the R3b class.
> - Create a separately frozen successor test if a robust observational
>   metric is desired.
>
> A future protocol could preregister, for each pressure level and spatial
> band:
>
> - arithmetic mean for continuity;
> - median as the robust central estimate;
> - 5th, 25th, 75th, 95th, 99th, and 99.9th percentiles;
> - maximum;
> - number and fraction above physically diagnostic thresholds;
> - mandatory-screen and cloud-rule populations side by side;
> - a prespecified robust classification statistic, such as median
>   model-minus-observation bias.
>
> The important point is that such a successor would answer a new question.
> It must not become a retroactive repair of R3b.
>
> **Branch disposition**
>
> The three named commits show the right separation:
>
> - 5bf677e: frozen erratum
> - ec59424: manifest registration
> - 6818da8: results
>
> My decision is to retain all three, preserve the appended result without
> reclassification, and push the full 11-commit branch after the ordinary
> repository integrity checks. Do not amend or squash the erratum,
> manifest, and result commits together. Their ordering demonstrates that
> the method and population were fixed before the result was recorded.
>
> **Final disposition text**
>
> R3b disposition: Accepted as Valid; diffusion state B Contradicted aloft.
> The frozen primary comparison over the global 147–215 hPa band gives
> B − MLS = +0.263. The result is stable across 30-, 40-, and 60-layer
> realizations and is strengthened by the preregistered conditional
> cloud-screening variant. State A is also inconsistent with MLS, so the
> result contradicts B without validating A. An extreme positive tail
> materially degrades the interpretability of some MLS arithmetic means,
> particularly the tropical 215 hPa and global 316 hPa diagnostics, but the
> tail acts conservatively on the primary B-minus-MLS discrepancy and
> therefore does not cause or weaken the registered classification. No
> causal claim is made.
>
> That is the result I would preserve as the formal R3b record.

### The R3b record, as the disposition sets it

**R3b: Valid. Diffusion state B Contradicted aloft.** Global 147–215 hPa,
B − MLS = +0.263; at 40 and 60 layers, +0.260 and +0.246. A − MLS = −0.663.

**Contradiction of B does not constitute confirmation of A; neither modeled
state agrees with the MLS upper-tropospheric reference.**

**The labels the disposition places beside the record:**

- **Estimator limitation:** the paragraph above, as worded there.
- **The ±20° band (−0.40):** retained; labeled an arithmetic-mean
  sensitivity to the 215 hPa tail, not representative of the vertical
  pattern.
- **The 316 hPa NCEP comparison:** computed but not interpretable as a
  mean-based comparison; outside the evidence for the class.
- **The H2O p < 1 hPa reading:** a reported interpretive check; the
  primary population is unchanged.

**A robust-estimator successor is not opened here.** If one is wanted, it is
frozen separately and answers a new question. It does not change this class.
