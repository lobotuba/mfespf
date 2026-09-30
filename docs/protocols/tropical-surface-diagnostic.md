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
