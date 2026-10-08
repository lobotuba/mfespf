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
