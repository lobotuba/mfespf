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
