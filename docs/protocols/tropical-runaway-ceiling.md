# Protocol: is the tropical runaway the saturated column's longwave ceiling?

**Frozen before the calculation.** Results are appended below a line at the
end; nothing above it changes afterwards. This is a diagnostic: it changes no
module and no default.

## Why

`ocean-heat-transport-circulation.md` found that the ±20° tropical column
runs away under R4's physics: N0 alone, and both of W1's columns together.

- **The state:** a surface near 314 K and RH about 0.95 through the free
  troposphere.
- **The failure:** the run left RRTMG's tables and was refused.

**The explanation proposed, not yet measured.** The band absorbs about
315 W m^-2 of sunlight and exports 42.51, so it must emit about 272 W m^-2.
A column whose free troposphere is held near saturation has a ceiling on its
outgoing longwave, the Simpson-Nakajima limit, about 280 W m^-2 (Goldblatt
et al. 2013, *Nature Geoscience* 6, 661). With its clouds, this column may
not reach 272 at any stable surface temperature. The observed tropics, at RH
0.43, would.

**If true,** the runaway is the column's humidity problem in its sharpest
form, and it says how dry the free troposphere must be to be stable.

## The calculation

**The radiation module itself solves each state,** with no time integration.

- **The module:** the radiation module of `ocean-heat-transport-circulation.md`'s
  N0 kernel. That is RRTMG with both source options, solar vapour
  absorption, solar ozone, the measured spectrum and the tropical ISCCP
  cloud field.
- **The band:** the ±20° sunlight (factor 0.29979).
- **How it is driven:** through `compute_tendencies` on prescribed inputs, as
  the CO2-doubling diagnostic of `global-column-chou-longwave.md` drove it.
- **The layers:** its 51 layers, unchanged.

**The states.** Surface temperature T_s from 290 to 330 K in 2.5 K steps.

- **Temperature.**
  - The surface air is at T_s.
  - A pseudo-adiabat from the surface, integrated on 1 Pa-log steps with the
    model's `moist_adiabatic_lapse_rate`, up to where it reaches 197.5 K.
    That is the observed tropical cold point, `tropical-cloud-field.md`.
  - Isothermal at 197.5 K above.
  - Each layer takes the profile's value at its pressure.
- **Humidity.** The vapour is RH times saturation (Murphy & Koop 2005, over
  ice below 273.16 K), in four cases:
  - **H100:** RH 1 below the tropopause, the Simpson-Nakajima limit;
  - **H95:** RH 0.95, the runaway's free troposphere (0.94-0.96 in N0 and
    W1);
  - **HOBS:** the observed ±20° band profile (`ocean-heat-transport-circulation.md`),
    interpolated in log pressure:
    - the 1000 hPa value, 0.779, below 1000 hPa;
    - the 300 hPa value, 0.597, from 300 hPa to the tropopause;
  - **H(r), reported:** HOBS below 850 hPa and a uniform r above, for r =
    0.3, 0.4, ... 0.9.
- **Above the tropopause:** 3.6 ppmv by volume, the observed 100-68 hPa band
  value (3.45-3.75).

**For each state:**

- the all-sky and clear-sky outgoing longwave;
- the absorbed sunlight, incident minus reflected;
- N = absorbed - OLR.

**The method's check.** The same module on N0's tail-mean state (temperature,
vapour, surface) must give N within 2 W m^-2 of N0's recorded 48.90. The
tolerance allows for the nonlinearity of a mean state.

**RRTMG's domain.** Each state's largest departure from RRTMG's reference
profile is reported. A state beyond 30 K anywhere is outside the tables, and
is excluded from the verdicts.

## The verdict

**Equilibrium.** The band is in equilibrium where N = 42.51, its export; the
ocean transport sums to zero over the band. A case has a **stable
equilibrium** where N - 42.51 changes sign from positive to negative as T_s
rises, between two states inside the tables.

| test | passes if |
|---|---|
| K0: the method | the module on N0's tail-mean state gives N within 2 W m^-2 of 48.90 |
| K1: the runaway at RH 0.95 | H95 has no stable equilibrium among the states inside the tables |
| K2: the observed humidity | HOBS has a stable equilibrium between 290 and 310 K |

**What the outcomes mean:**

- **K1 and K2 pass:** the explanation is supported. At the runaway's humidity
  the column cannot shed what it absorbs at any surface temperature, and at
  the observed humidity it can.
- **K1 fails:** it is refuted. A fixed RH of 0.95 does have an equilibrium,
  so the runaway is something other than the ceiling at fixed humidity.
- **K0 fails:** the method does not reproduce N0, and neither K1 nor K2 is
  read.

**Reported, not judged:**

- the curves;
- the clear-sky ceiling;
- the largest r for which H(r) has a stable equilibrium. That is how dry the
  free troposphere must be to stay stable.

---

## Erratum, before any verdict (nothing above this line changed)

**The prescribed stratosphere is outside RRTMG's tables in every state.**

- **The cause.** The protocol holds the air isothermal at 197.5 K from the
  tropopause to the 1 hPa top. The real stratosphere, and RRTMG's reference
  profile, warm to about 270 K at 1 hPa.
- **The first run.** Every state's largest departure from RRTMG's reference
  is the same 78 K, at the top. Every state was refused by the solver, as a
  negative optical depth in the cold upper stratosphere.
- **No verdict was read.** K0 passed (N = 49.53 against N0's 48.90), so the
  method stands.

**The two readings:**

1. **As frozen:** every state is excluded, K1 passes vacuously and K2 fails
   vacuously. That measures the prescription, not the column.
2. **Above the tropopause, the larger of 197.5 K and the AFGL tropical
   temperature at that pressure,** taken.
   - The source is Anderson et al. (1986) Table 1a, as AER's
     `input_rrtm_TROP-clr` carries it (already in `tests/data/aer_rrtmg_lw/`
     and checksummed in `data/sources.yaml`), interpolated in log pressure.
   - Its cold point is 196.35 K at 109.6 hPa, and it is 268.8 K at 1 hPa.
   - The troposphere, the humidity, the tests and their tolerances are
     unchanged.

## Result (appended; nothing above this line changed)

**Under the erratum's reading, the states up to T_s = 305 K are inside
RRTMG's tables:** the largest departure is 18.6-25.7 K. From 307.5 K the
tropopause rises far enough that the pseudo-adiabat's upper troposphere
departs by more than 30 K. Those states are reported and excluded, as
frozen. Above about 317.5 K the solver refuses: the gas optics go negative.

**K0 passes.** The module on N0's tail-mean state gives N = 49.53 W m^-2,
against N0's recorded 48.90.

All values in W m^-2:

| T_s | H95: OLR | clear | absorbed | **N** | HOBS: OLR | clear | absorbed | **N** |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 290.0 | 221.9 | 254.7 | 299.4 | 77.5 | 231.8 | 269.4 | 297.1 | 65.3 |
| 295.0 | 233.5 | 260.9 | 301.9 | 68.3 | 246.0 | 279.2 | 299.1 | 53.2 |
| 297.5 | 239.5 | 263.5 | 303.2 | 63.7 | 253.5 | 283.8 | 300.3 | 46.8 |
| 300.0 | 245.2 | 265.5 | 304.7 | 59.5 | 260.9 | 288.1 | 301.5 | **40.5** |
| 302.5 | 250.8 | 267.3 | 306.3 | 55.5 | 268.2 | 291.8 | 302.8 | 34.6 |
| 305.0 | 255.8 | 268.7 | 307.9 | **52.1** | 274.9 | 294.9 | 304.1 | 29.3 |
| 310.0 (outside the tables) | 263.4 | 270.5 | 311.5 | 48.1 | 285.3 | 298.9 | 307.0 | 21.7 |
| 317.5 (outside the tables) | 272.2 | 274.4 | 317.3 | 45.1 | refused | | | |

The band needs N = 42.51.

**The verdicts as frozen:**

- **K1 passes.** H95 has no stable equilibrium among the states inside the
  tables: N falls from 77.5 to 52.1 and stays above 42.51. Even the
  extrapolated states reach only 45.1 at 317.5 K.
- **K2 passes.** HOBS crosses 42.51 between 297.5 and 300.0 K: a **stable
  equilibrium near 298.75 K**, against the observed band's 299.61 K.

**The explanation is supported.** At the runaway's humidity the column
cannot shed what the tropics absorb at any surface temperature inside the
tables. At the observed humidity, the same radiation, clouds and sunlight
equilibrate within 0.9 K of the observed tropical surface.

**Reported:**

- **How dry the free troposphere must be.** HOBS below 850 hPa and a uniform
  r above:

  | r | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 |
  |---|---|---|---|---|---|---|---|
  | stable equilibrium near | 296.25 K | 298.75 | 298.75 | 301.25 | 303.75 | **none** | **none** |

  The stability threshold lies between 0.7 and 0.8.
- **The clear-sky ceiling.** At saturation (H100) the clear-sky OLR rises by
  only 3 W m^-2 from 300 to 305 K (263.6 to 266.6), and by 5.6 more by
  317.5 K outside the tables. The all-sky OLR is lower still: the cloud field
  traps 15-30 W m^-2 of it.
- **The CERES comparison.** The observed band OLR (±30°) is 259.3 W m^-2.
  At HOBS's equilibrium the model's is about 257.
- **The caveat.** The temperature profile is a prescribed pseudo-adiabat, not
  the column's own convection, and the humidity is prescribed. This shows
  what humidity the column needs. It does not show a mechanism that delivers
  it.

**What it says.** The runaway, the global warm bias and the cloud effect's
collapse are one finding seen three ways: the free troposphere is too humid.
With the observed humidity the tropical column is stable and right. The
drying mechanism must bring the free troposphere below about 0.75, and the
two-column overturning could not, because the moist runaway erased the
surface contrast it needed.
