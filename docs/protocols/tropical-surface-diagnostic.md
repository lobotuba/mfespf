# Protocol: why the tropical column's surface is 16 K cold

**Frozen before the diagnostic is written and before any calculation.** Results
are appended below a line at the end; nothing above it changes afterwards.

This is a diagnostic. It changes no module and adopts nothing.

## Why

`tropical-heating-diagnostic.md` found that 77 percent of the missing heating
at the tropical tropopause is the column's cold troposphere and surface. The
+/-25 degree column of `tropical-column.md` attempt 2 equilibrates at 282.16 K.
The band-mean 1000 hPa air temperature is 298.45 K (NCEP R1 1991-2020).

A column's surface temperature is set by two things:

- **its top-of-atmosphere budget.** What it absorbs from the sun, less what it
  exports, must leave as outgoing longwave.
- **its structure.** Its humidity and its lapse rate set how warm the surface
  must be for that outgoing longwave to leave.

A 16 K deficit can come from either. This protocol asks which, at the observed
state.

## The method

The model's radiation module is evaluated offline on the observed tropical
state, as in `tropical-heating-diagnostic.md`. It uses the same column, built
and stepped 10 times, and the same attempt-2 equilibrium state
(`trop2_state_25_up.npz`). Nothing is integrated.

**The observed state**, declared here:

- **Temperature, surface, stratospheric water and ozone:** exactly as in
  `tropical-heating-diagnostic.md`.
- **Tropospheric humidity:**
  - NCEP R1 1991-2020 long-term-mean specific humidity over +/-25, cos-weighted,
    12 months equal, at its levels from 1000 to 300 hPa;
  - SWOOSH v2.73 water at 261 hPa and above, as there;
  - interpolated in log pressure between.

**The observed top-of-atmosphere fluxes.** NCEP R1 1991-2020 over +/-25, with
NCEP's global bias removed component by component, the assumption
`tropical-column.md` already made for the net. For each of the absorbed
sunlight and the outgoing longwave:

  observed band value = NCEP band value + (CERES global - NCEP global).

CERES is EBAF Ed4.0 with its net constraint (Loeb et al. 2018, *J. Climate*
31, 895-918, Table 5): incoming 340.0, reflected 99.1 and outgoing longwave
240.1 W m^-2.

**The export** is the column's own: 38.464 W m^-2.

## The quantities

Each is in W m^-2, and each is positive when it makes the model lose more
energy at the observed surface temperature than the observed atmosphere does.

- **N, the model's budget on the observed state.** ASR_m(obs) - OLR_m(obs) -
  38.464, where ASR is absorbed sunlight. N < 0 means the model cannot hold
  the observed state, and must cool. Its shortfall, -N, is split into two:
  - **L, the longwave:** OLR_m(obs) - OLR_obs;
  - **S, the shortwave:** ASR_obs - ASR_m(obs).

  L + S = -N, up to the observed state's own residual, which is reported.
- **H, the column's humidity.** The column's own equilibrium relative humidity
  profile, from `trop2_state_25_up.npz`, is placed on the observed temperatures
  below 100 hPa. H is OLR_m with it minus OLR_m(obs). A drier column lets more
  longwave out at the same temperatures.
- **G, the column's lapse rate.** The column's own equilibrium temperatures
  below 150 hPa are shifted by the difference between the observed surface and
  its own. The observed relative humidity is kept. G is OLR_m with that minus
  OLR_m(obs). A profile warmer aloft than the observed one loses more.

The total is D = L + S + H + G.

## The verdict

- **Carrying the deficit.** A candidate carries it if its share of D is at
  least a third. More than one can, and the largest is named first.
- **When nothing is missing.** If D <= 0, the column is not losing more energy
  than the observed state at 298 K. Then the cold surface is something this
  diagnostic cannot see, and that is reported as the result.

| carries it | what it means | next |
|---|---|---|
| L | The model's greenhouse is too weak on the tropical state. The vapour's wing decay was fitted globally, and also stands in for the cloud longwave; the tropics have more high cloud. | the observed cloud field (`observed-cloud-field.md`) on the tropical column, with the vapour refit it forces |
| S | The model absorbs too little sunlight in the tropics. Its 0.30 surface albedo stands in for a planetary albedo, globally. | the tropical reflection: surface and cloud shortwave |
| H | The column's water cycle keeps its tropical troposphere too dry. | the humidity closure in a tropical column |
| G | The column's lapse rate is wrong for the tropics. | the convective closure's tropical lapse rate |

**Reported, not judged:**

- the observed state's own residual;
- L against NCEP's clear-sky outgoing longwave, bias-adjusted the same way with
  CERES's clear-sky 268.1;
- the column's own equilibrium ASR, OLR and relative humidity profile, against
  the observed ones.

**Not tested here:** where the export leaves from. The surface module accepts
no lateral export, so that question needs code, not this diagnostic.

---

## Results (appended; nothing above this line changed)

**The observed fluxes.** NCEP R1 over +/-25, bias-adjusted to CERES EBAF
Ed4.0 component by component:

| | NCEP band | NCEP global | CERES global | observed band |
|---|---:|---:|---:|---:|
| absorbed sunlight | 406.19 - 124.92 = 281.27 | 341.90 - 117.04 = 224.87 | 240.9 | **297.30** |
| outgoing longwave | 255.67 | 237.74 | 240.1 | **258.04** |
| clear-sky outgoing longwave | 291.76 | 269.13 | 268.1 | **290.73** |

The observed state's own residual, ASR - OLR - 38.464, is +0.80 W m^-2.

**The observed humidity.** The band's precipitable water below 100 hPa is 38.2
kg m^-2; the column's own is 11.1.

**The model on the observed state:** ASR 289.87 and OLR 284.73 W m^-2, so
**N = -33.32**. The model cannot hold the observed state: it loses 33 W m^-2
more than it gains.

| quantity | W m^-2 | share of D |
|---|---:|---:|
| L, longwave: 284.73 - 258.04 | **+26.69** | 16.0 |
| S, shortwave: 297.30 - 289.87 | +7.43 | 4.46 |
| H, the column's own relative humidity | **-37.66** | -22.6 |
| G, the column's own lapse rate, shifted 16.29 K | +5.21 | 3.12 |
| **D** | **+1.67** | |

L + S = 34.12 against -N = 33.32, with the 0.80 residual between them.

**The verdict as frozen:** D > 0, and L, S and G each exceed a third of it. So
all three carry it, and L is the largest.

**Why the shares must not be read as fractions.**

- **The rule assumed no large term of the opposite sign.** The shares were to
  be fractions of a deficit, on that assumption. H is -37.7, and D is a small
  difference of large terms, so the shares run to 16 and -22.
- **What survives.** The ordering does: L is the largest term that makes the
  model lose energy, then S, then G.
- **The rule stands.** It is not edited. It is reported as having lost its
  resolving power on this case.

**What the terms say, reported.**

- **L: the tropics have no cloud longwave effect in this model.**
  - **The model's number.** On the observed state its OLR, 284.73, is 6.0
    W m^-2 below the observed *clear-sky* 290.73.
  - **The observed number.** Observed all-sky is 258.04, so the band's cloud
    longwave effect is 32.7 W m^-2 (CERES globally: 27.9).
  - **Where the model stands.** The vapour's globally fitted wing decay stands
    in for cloud longwave, and here it supplies about 6 of those 32.7.
- **S: the tropics reflect less than the 0.30 albedo says.** The model absorbs
  7.4 W m^-2 less sunlight than the observed band.
- **H is negative because the column's own troposphere is relatively
  humid.** RH is 0.74 to 0.98 between 800 and 390 hPa, against the observed
  0.36 to 0.54. Its low precipitable water is its cold temperatures, not low
  relative humidity.
- **Above 320 hPa, H includes an artefact.** The column's time-mean vapour
  exceeds saturation at its time-mean temperature by 1.24x at 322 hPa, 4.4x at
  254 hPa and 2.2x at 187 hPa. Saturation is convex in temperature, so a
  temperature that cycles puts the mean of the saturation above saturation at
  the mean. The cycle was not recorded, so this is not measured. H above
  320 hPa is therefore not a clean humidity comparison.
- **G is small.** The column's lapse rate is steeper than observed, 7 to
  17 K colder at 500 to 250 hPa after the shift. It costs +5.2 W m^-2 once the
  drying that comes with it is counted.

**Next, by the frozen table:** L is the largest. The next step is the observed
cloud field (`observed-cloud-field.md`) on the tropical column, with the vapour
refit it forces. S follows: the tropical reflection. Before either, the
column's upper-tropospheric supersaturation in the time mean needs its cycle
measured, because it bears on every humidity result above 320 hPa.
