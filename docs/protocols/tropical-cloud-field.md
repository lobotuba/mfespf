# Protocol: the observed tropical cloud field on the tropical column

**Frozen before the implementation, the refit and any run.** Results are
appended below a line at the end; nothing above it changes afterwards.

## Why

`tropical-surface-diagnostic.md` found where the tropical column's energy goes
on the observed state.

- **It loses 33 W m^-2 more than it gains.**
- **The largest part is longwave, 26.7.** On the observed state the model's
  outgoing longwave is 6 W m^-2 below the observed *clear-sky* value, while
  the observed band's cloud longwave effect is 32.7 (+/-25).
- **Where the cloud effect went.** The vapour's globally fitted wing decay
  stands in for cloud longwave, and it supplies little of it in the tropics.
- **Shortwave is next, 7.4.**

The mechanism to add is the tropics' own clouds, at their observed heights.

## The column

The tropical column of `tropical-column.md` attempt 2, at **+/-30 degrees**
because that is ISCCP's tropical band:

- band insolation, the 32.95 W m^-2 export and the upwelling per unit area;
- CO2 band, ozone, measured spectrum, diffusion;
- 51 layers to 1 hPa;
- **both** `condense_after_convection` and `end_of_step_saturation`, as
  decided after `end-of-step-saturation.md`.

**Two configurations:**

- **Control:** that column as it is, vapour wing decay 4601.27 m^-1.
- **Clouds:** the same column with three changes.
  1. **The observed tropical cloud field.** Rossow and Schiffer (1999), *BAMS*
     80, 2261-2287, Table 5, daytime cloud-type amounts for the tropics
     (+/-30), ice (liquid):

     | level | ISCCP types | amount | ice fraction |
     |---|---|---:|---:|
     | high | cirrus 15.6 + cirrostratus 5.5 | 0.211 | 1 |
     | deep convective | 2.7 | 0.027 | 1 |
     | middle | altocumulus 0.9 (6.1), altostratus 0.1 (4.7), nimbostratus 0.0 (1.0) | 0.128 | 1.0 / 12.8 |
     | low | cumulus 0.0 (12.3), stratocumulus 0.0 (10.7), stratus 0.0 (0.8) | 0.238 | 0 |

     - **Top pressures and optical thicknesses.** They are the global Table 4
       values the global field already uses (267, 326, 557 and 826 hPa; tau
       2.2, 35.6, 4.8 and 4.7). Table 5 gives no tropical ones. This is a
       declared approximation: tropical high cloud is observed higher than
       the global mean.
     - **The surface reflectance** is ISCCP's 0.17, as the module sets it.
  2. **Water vapour's solar absorption** (`solar-vapour-absorption.md`). It
     passed its line-by-line benchmark. `observed-cloud-field.md` found the
     cloud field over-reflects without it (albedo 0.338 against 0.283 with it).
  3. **The vapour wing decay refitted with the cloud field present.** Until
     now the wing decay has stood in for the clouds' longwave effect.
     - **The method is unchanged.** The calibration profile is the same, the
       target is still OLR 238.175 W m^-2, and the CO2 band is in.
     - **What is new.** The *global* ISCCP cloud field (`ISCCP_CLOUD_TYPES`)
       is in the calibration as amount-weighted sub-columns, solved as the
       radiation module solves them.
     - **Why the global field.** The calibration is global, so it takes the
       global field. Fitting it to the tropical field would be fitting to the
       thing being tested.

## The tests

Against the +/-30 observations, NCEP R1 1991-2020 bias-adjusted to CERES
EBAF Ed4.0 component by component, as in `tropical-surface-diagnostic.md`.
The tolerances are `observed-cloud-field.md`'s.

| test | observed | passes if |
|---|---:|---|
| K1: planetary albedo | 0.265 | within 0.03 |
| K2: shortwave cloud radiative effect | -55.6 W m^-2 | within 12 |
| K3: longwave cloud radiative effect | +30.6 W m^-2 | within 8 |
| K4: the surface | 297.8 K (band-mean 1000 hPa air) | the clouds run closes at least half the control's deficit: T_s(clouds) - T_s(control) >= 0.5 (297.8 - T_s(control)) |

Both runs must be valid (T1 of `tropical-column.md`); a run that is not fails
everything.

**Reported, not judged:**

- the refitted wing decay;
- outgoing longwave and absorbed sunlight, against 259.3 and 293.1;
- the clear-sky outgoing longwave, against 289.9;
- `tropical-column.md`'s T2 to T4 at +/-30 for both runs: 100 hPa 197.5,
  70 hPa 200.4, 50 hPa 207.4, 30 hPa 215.8 and 10 hPa 230.6 K;
- the water at 100 to 68 hPa;
- the model's temperature at each cloud top.

**Adoption.** Nothing is adopted into the reference column by this protocol.
If K1 to K4 pass, the tropical column takes the cloud field, the solar vapour
absorption and the refit as its configuration for the stage-1 retest under
`tropical-column.md`, which follows as attempt 3. If any fails, the finding is
recorded and nothing is taken.

---

## Results (appended; nothing above this line changed)

**The refit.** With the global ISCCP cloud field in the calibration, the vapour
wing decay is 3411.86 m^-1 (it was 4601.27 without), holding OLR at
238.175 W m^-2. Reported on the calibration profile:

- the cloud longwave effect is 38.9 W m^-2, against CERES's global 27.9. The
  calibration profile is cold aloft, so its cloud tops are cold.
- the forcing per CO2 doubling is 2.38 W m^-2 at the top and 4.27 at 200 hPa
  (3.03 and 4.94 without clouds), because the clouds mask part of the band.

**The runs.** +/-30 degrees, both saturation options, 200000 steps. Both are
complete and valid, with budgets closing to 1e-15.

| test | observed | control | clouds | verdict |
|---|---:|---:|---:|---|
| K1: planetary albedo | 0.265 | (0.282, surface albedo 0.30) | **0.251** | **pass** (0.014) |
| K2: shortwave cloud effect | -55.6 W m^-2 | | **-51.2** | **pass** (4.4) |
| K3: longwave cloud effect | +30.6 W m^-2 | | **+24.2** | **pass** (6.4) |
| K4: surface | 297.8 K | 283.75 K | **299.09 K** | **pass**: needs at least 290.78, and it closes 109% of the deficit |

**Reported.**

| | observed | control | clouds |
|---|---:|---:|---:|
| absorbed sunlight | 293.1 | 286.5 | 299.0 W m^-2 |
| outgoing longwave | 259.3 | 251.1 | 264.1 W m^-2 |
| clear-sky outgoing longwave | 289.9 | | 288.3 W m^-2 |
| precipitable water | | 13.4 | 78.4 kg m^-2 |
| cold point | 197.5 K at 100 hPa | 172.2 K at 103 hPa | **187.6 K at 57 hPa** |
| T at 100 hPa | 197.5 | 172.9 | 191.5 K |
| T at 70 hPa | 200.4 | 179.6 | 187.9 K |
| T at 50 hPa | 207.4 | 186.0 | 188.3 K |
| T at 30 hPa | 215.8 | 199.4 | 196.1 K |
| T at 10 hPa | 230.6 | 231.8 | 225.1 K |
| water, 100-68 hPa | 3.45-3.75 ppmv (+/-25) | 0.10 | **3.20 ppmv** |

Cloud-top temperatures in the clouds run:

- high, 267 hPa: 242.1 K;
- deep convective, 326 hPa: 254.8 K;
- middle, 557 hPa: 278.3 K;
- low, 826 hPa: 292.3 K.

**The tropical surface is fixed by its clouds.** The cold tropical column of
`tropical-surface-diagnostic.md` was missing its cloud longwave effect and
over-reflecting. With the observed tropical cloud field, and the vapour no
longer standing in for cloud, each failure is resolved:

- the surface moves 15.3 K to within 1.3 K of the observed;
- both cloud effects are within tolerance;
- the clear-sky outgoing longwave is 1.6 W m^-2 from the observed.

**The tropopause follows the troposphere, as `tropical-heating-diagnostic.md`
found it would.**

- **The cold point warms 15 K**, to 187.6 K, against the observed 197.5.
- **The stratospheric water is 3.20 ppmv**, against the band's observed 3.45
  to 3.75. That is the first valid, equilibrated run of this column to dry
  the entering air to near the observed amount. Attempt 1's 5000-step check
  reached 3.83, but it was neither.
- **The coldest layer has moved up to 57 hPa**, where the observed is
  100 hPa.
- **Above 70 hPa the column is 12 to 20 K too cold.** The heating diagnostic
  attributed the deficit there to its remainder, which it inferred, not
  measured, as most likely ozone's 9.6-micron band.

**Wetter than observed?** The precipitable water is 78 kg m^-2. No tropical
observation was declared here, so it is not judged. It is large, and it is
flagged for the next step.

**Adoption, as frozen.** K1 to K4 pass, so the tropical column takes the
tropical cloud field, the solar vapour absorption and the refit, 3411.86
m^-1, for the stage-1 retest. That is `tropical-column.md` attempt 3, frozen
before it runs. The reference column is not changed.
