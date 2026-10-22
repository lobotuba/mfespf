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
