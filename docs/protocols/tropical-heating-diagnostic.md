# Protocol: where the tropical lower stratosphere's missing heat is

**Frozen before the diagnostic is written and before any calculation.** Results
are appended below a line at the end; nothing above it changes afterwards.

This is a diagnostic. It changes no module and adopts nothing. Its job is to
say which mechanism the next protocol should build.

## Why

`tropical-column.md` attempt 2 failed:

- the +/-25 degree column's cold point is 166.7 K, against the observed 196.4;
- 70 to 30 hPa are 18 to 24 K too cold;
- without the upwelling, the cold point is still 174.4 K.

So the lower stratosphere is radiatively too cold before any dynamics act.
`co2-band.md`'s global column shows the same thing, 20 K too cold at 60 hPa.

For comparison, Dacie et al. (2019, *J. Climate* 32, 6769-6782) ran a 1D
clear-sky radiative-convective model. It used RRTMG, which includes ozone's
9.6-micron band, an equatorial ozone profile and a 297.7 K surface. With no
upwelling, its cold point is 200.1 K (relaxed adjustment) or 201.3 K (hard).
Their upwelling of 0.25 mm/s cools it by 9.3 to 10.3 K. That is the order of
this column's own upwelling, 0.23 mm/s at 70 hPa, which cools it by 11 K.

The model is missing about 26 K at its radiatively determined cold point.
Three candidates are named, and a fourth can be checked.

| candidate | the argument |
|---|---|
| A. The troposphere | The lower stratosphere is heated by absorbing infrared from below. The column's surface is 282.2 K; the tropical ocean is near 300 K. |
| B. The ozone profile | The column uses SWOOSH's global mean. The tropics have less ozone near 100 hPa and more above. |
| C. Ozone's 9.6-micron band | It is absent, as declared. Gettelman et al. (2004, *JGR* 109, D07103) put ozone's longwave absorption as the heating that the CO2 and vapour cooling partly balance around the tropical tropopause. |
| D. The CO2 band near the tropopause | It passed its tests on outgoing flux and on the stratopause cooling, but not on the lower stratosphere. Gettelman et al. (2004) find CO2 heats slightly just around the tropical tropopause. |

## The method

The model's own radiation module is evaluated offline, as in `model-top.md`'s
water diagnostic. Nothing is run to equilibrium.

**The column.** The +/-25 degree tropical column of attempt 2:

- the CO2 band;
- 51 layers to 1 hPa;
- vapour wing decay 4601.27 m^-1;
- band insolation;
- the same modules.

It is built, stepped 10 times so that every input exists, and then the state
the radiation sees is replaced.

**The observed profile.** Declared here:

- **Temperature.** The NCEP R1 1991-2020 long-term mean over +/-25 degrees,
  cos-latitude weighted, at its 17 levels from 1000 to 10 hPa, interpolated in
  log pressure. Above 10 hPa it is the attempt-2 +/-25 run's own profile, which
  matches at 10 hPa to 0.6 K.
- **Surface.** The surface temperature is the band-mean 1000 hPa air
  temperature.
- **Water.** Below 100 hPa it is the column's own assumption, 0.7 of
  saturation. At and above 100 hPa it is the SWOOSH v2.73 1991-2020 mean of
  the combined anomaly-filled field over +/-25 degrees. That is the bins
  centred on +/-5 and +/-15 degrees whole, plus the half of each +/-25 bin that
  lies inside 25 degrees, with every part cos-weighted.
- **Ozone.** The column's global SWOOSH profile, unless a diagnostic changes
  it.

**The budget.**

- **Radiative heating.** Q_rad is the radiation module's layer heating, in K
  per day.
- **Dynamical heating.** Q_dyn is the tropical exchange module's heating from
  the prescribed upwelling, evaluated on the same profile.
- **The missing heating** is M = -(Q_rad + Q_dyn), at each layer from 150 to
  30 hPa. It is what the column lacks to hold the observed profile steady.
- **The diffusion.** It is left out: vapour eddy diffusion moves water, not
  heat.

## The diagnostics

All are reported at the layers nearest 150, 100, 85, 70, 50 and 30 hPa. Their
mean over the layers centred between 100 and 70 hPa decides.

- **D0, the tool.** On the attempt-2 +/-25 run's own time-mean equilibrium,
  Q_rad + Q_dyn must be within 0.05 K/day of zero at every layer above
  100 hPa, where nothing else heats. If it is not, the tool is wrong and
  nothing below is reported as a finding.
- **D1, the budget on the observed profile:** M at each level. If M <= 0.1
  K/day at 100 to 70 hPa, the model's radiation can hold the observed
  profile, and the cold bias is not missing heat. It is then reported as
  that, and D2 to D4 are reported without a verdict.
- **D2, candidate A, the troposphere.** Q_rad again, with the observed profile
  above 150 hPa but the attempt-2 run's own temperatures below it and its own
  surface of 282.16 K. The difference from D1's Q_rad is A's share of M.
- **D3, candidate B, the ozone profile.** Q_rad on the observed profile with
  the ozone replaced by SWOOSH's +/-25 degree 1991-2020 mean, formed the same
  way as the water. The
  difference from D1 is B's share.
- **D4, candidate D, the CO2 band.** Q_rad on the observed profile with the
  CO2 band removed. The difference is the band's contribution. Reported: its
  sign at 100 hPa, against the slight heating Gettelman et al. (2004) find.
- **D5, reported only.** The same removal for the vapour above 150 hPa, and the
  shortwave part of Q_rad (the incident sunlight set to zero, subtracted).

## The verdict

A candidate **carries the deficit** if its share, averaged over 100 to 70 hPa,
is at least half of M. More than one can carry it.

**Candidate C** cannot be computed by this model. Its share is inferred as the
remainder, M minus the shares of A and B. It is reported as a remainder, never
as a measurement. If the remainder is at least half of M, the ozone band is
named as the likely next mechanism, and its own protocol must measure it.

**Candidate D** is not a share of M. The band exists. If its heating at 100 hPa
is below -0.2 K/day, it is named as inconsistent with Gettelman et al. near
the tropopause, and it comes before C.

**What happens next**, by the verdict:

| verdict | next |
|---|---|
| A carries it | Why the tropical surface is 18 K cold. The heat export is taken from the atmosphere, the vapour is calibrated globally, and no clouds are resolved. Its own diagnostic follows. |
| B carries it | A band-mean ozone profile for the tropical column: an observation, not a mechanism. |
| The remainder carries it | Ozone's 9.6-micron band, from a citable, public parameterization, benchmarked like the CO2 band. |
| D is inconsistent | The CO2 band near the tropopause, before anything else. |
| M <= 0.1 | The cold bias is not missing heating, so the next place to look is the column's water and convection. |
