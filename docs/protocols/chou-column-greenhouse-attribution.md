# Protocol: is R1's warm surface its gas optics or its vapour amount?

**Frozen before the calculation.** Results are appended below a line at the
end; nothing above it changes afterwards. This is a diagnostic. It changes no
module and no default.

## Why

`global-column-chou-longwave.md` R1 (Chou's vapour with ozone's 9.6-micron
band) failed G4. Its surface is 297.91 K against 288.0, while the fitted
control is at 289.59.

**At the top of the atmosphere R1 is close to observed:**

- its clear-sky outgoing longwave is 269.6 W m^-2 against CERES's 268.1;
- both cloud effects pass.

**Its column holds far too much water:**

- precipitable water 70.6 kg m^-2 against about 25;
- relative humidity 0.89 at 300-700 hPa;
- a clear-sky greenhouse effect of about 175 W m^-2 against an observed 130.

**Two causes are possible.**

1. **The gas optics.** Chou's scheme was 6.19 W m^-2 too opaque on AFGL
   mid-latitude summer (V1). The vapour feedback amplifies that. RRTMG
   addresses this.
2. **The vapour amount.** The column's humidity closure puts too much water
   in it. Gas optics cannot fix this.

The user chose RRTMG before R1's result. This diagnostic measures, before that
build, how much of R1's greenhouse effect each cause carries.

## The state

R1's tail-mean vapour profile was not saved; the run script kept summaries
only. So it is restored.

**The restart.**

- R1's configuration, unchanged (`global_chou.py r1`).
- Started from R1's recorded tail-mean temperature profile and its surface
  temperature.
- The initial water is the configuration's own; the water module spins it up.
- 20000 steps, with the means taken over the last fifth.
- **The output:** tail means of the layer temperature, vapour, mass and
  pressure, and the surface temperature.

**S0, the state is R1's,** judged before anything else:

| check | passes if |
|---|---|
| validity | complete and valid |
| imbalance | \|imbalance\| <= 0.1 W m^-2 |
| surface | \|T_s - 297.91\| <= 0.5 K |
| water | \|precipitable water - 70.6\| <= 3 kg m^-2 |

**If S0 fails,** the restart's result is recorded and R1 is rerun in full,
200000 steps, saving its state. The calculation then uses that state.

## The calculation

**The radiation module itself solves each case,** with the state's
temperatures, surface temperature and clouds, the global ISCCP field. Only
two things change between cases:

| case | gas optics | vapour |
|---|---|---|
| A | R1's (Chou's vapour, CO2 band, ozone band) | the state's |
| B | R0's (the fitted envelope at 3411.86 m^-1, CO2 band, no ozone band) | the state's |
| C | R1's | observed |
| D | R0's | observed |

**The observed vapour,** placed on the state's layers by interpolation in log
pressure:

- **1000-300 hPa:** NCEP R1 specific humidity, 1991-2020 long-term mean,
  global and cos-latitude weighted. This is the file
  `tropical-surface-diagnostic.md` used.
- **At and above 262 hPa:** SWOOSH v2.73 water vapour, 1991-2020, global and
  cos-latitude weighted.

**Reported for each case:**

- the outgoing longwave;
- the clear-sky outgoing longwave;
- the clear-sky greenhouse effect, G = σT_s^4 - clear-sky outgoing
  longwave;
- the surface's downward longwave, all-sky and clear-sky.

## The verdict

- **ΔG_vap = G(A) - G(C):** the excess vapour's greenhouse effect, with R1's
  optics.
- **ΔG_opt = G(A) - G(B):** R1's optics against the fitted control's, on R1's
  vapour.
- **The scale of a gas-optics error:** 6.19 W m^-2, Chou's measured
  outgoing-longwave miss on mid-latitude summer.

| finding | means |
|---|---|
| ΔG_vap > 2 x 6.19 = 12.4 W m^-2 | the vapour amount carries more greenhouse effect than any gas-optics error measured for this scheme could. RRTMG alone cannot close G4, and the free-tropospheric humidity is the larger mechanism to add |
| ΔG_vap <= 12.4 W m^-2 | not established: a gas-optics error of measured size is comparable, and RRTMG may matter as much |

**Reported, not judged:**

- ΔG_opt;
- the case-D greenhouse effect;
- the observed clear-sky greenhouse effect, 398 - 268.1 ≈ 130 W m^-2. The
  surface emission is the CERES EBAF surface value (Loeb et al. 2018),
  approximate.

**Whatever the finding,** RRTMG proceeds as the user decided. This tells
the user what it can be expected to do.

---

## S0 result (appended; nothing above this line changed)

**The restart.** R1's configuration ran 20000 steps from R1's recorded
profile and surface temperature, and took 2308 s.

| check | limit | restart | verdict |
|---|---|---:|---|
| validity | complete and valid | complete, valid, no domain violations | pass |
| imbalance | \|imbalance\| <= 0.1 W m^-2 | 0.40 | **fail** |
| surface | within 0.5 K of 297.91 | 297.33 K (-0.58) | **fail** |
| water | within 3 of 70.6 kg m^-2 | 66.9 (-3.7) | **fail** |

**S0 fails.** 20000 steps from a dry start do not return the column to R1's
state. The relative humidity at 300-700 hPa is 0.88 and the clear-sky
outgoing longwave is 269.8.

**As frozen,** R1 is rerun in full, 200000 steps, saving its tail-mean
state. The calculation uses that state.

## The state from R1's full rerun, and the result (appended; nothing above this line changed)

**The full rerun reproduces R1 to every printed digit.**

- **The run:** 200000 steps, 297.90841231119094 K, the same as R1, with
  precipitable water 70.60 and imbalance 0.006. The run is deterministic.
- **The saved state:** `gcl_r1_200000_state.npz`, which is R1's.
- **S0 passes on it:** valid, imbalance 0.006, surface error 0.00 K, water
  error 0.00.

**The calculation.** Case A's clear-sky outgoing longwave, 269.65,
reproduces R1's 269.62: the same state, solved by the same module.

All values in W m^-2:

| case | optics | vapour (kg m^-2) | OLR | clear-sky OLR | G, clear sky | surface down | surface down, clear |
|---|---|---:|---:|---:|---:|---:|---:|
| A | R1's | 70.6 | 246.17 | 269.65 | **176.97** | 422.41 | 412.18 |
| B | R0's | 70.6 | 261.10 | 287.65 | **158.97** | 416.99 | 408.22 |
| C | R1's | 26.2, observed | 273.57 | 309.35 | **137.28** | 392.38 | 357.79 |
| D | R0's | 26.2, observed | 295.86 | 340.51 | **106.11** | 382.65 | 347.16 |

- **ΔG_vap = G(A) - G(C) = 39.69 W m^-2.**
- **ΔG_opt = G(A) - G(B) = 18.00 W m^-2.**

**Verdict as frozen: ΔG_vap > 12.4, so the vapour amount dominates.** R1's
excess water carries 39.7 W m^-2 of clear-sky greenhouse effect. That is
6.4 times the largest gas-optics error measured for Chou's scheme (6.19 on
MLS). RRTMG alone cannot close G4; the free-tropospheric humidity is the
larger mechanism to add.

**Reported.**

- **The observed vapour** gives 26.2 kg m^-2 on R1's layers, against the
  observed global ~25. NCEP and SWOOSH are placed as frozen.
- **ΔG_opt = 18.0.**
  - It is large because the fitted envelope is transparent, not because
    Chou's scheme is opaque: R0's clear-sky OLR was 15.5 W m^-2 above CERES.
  - The envelope was fitted to all-sky OLR with clouds that were themselves
    missing their cloud tops' cold. So it compensated, and on a wet column
    it lets out what the real vapour would not.
  - On the observed vapour the two optics differ by 31.2 (D against C).
- **Case C** (R1's temperatures, observed vapour) has a clear-sky
  greenhouse effect of 137.3 W m^-2, against the observed 130.
  - The temperatures are R1's warm ones, so this is not a prediction.
  - Its order of magnitude shows that with the observed water the clear-sky
    greenhouse effect falls back near the observed one under Chou's optics.
- **What follows is the user's decision.** As frozen, RRTMG proceeds (R3
  running). The measured lever is the humidity. A mechanism that dries the
  free troposphere (0.89 relative humidity at 300-700 hPa against an
  observed 0.4-0.5) is the next physics. It is named in
  `observed-cloud-field.md` as subsiding air, and left unmeasured until now.
