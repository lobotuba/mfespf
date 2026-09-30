# Protocol: a tropical column freeze-dries the air entering the stratosphere

**Frozen before the mechanism is written and before any run.** Results are
appended below a line at the end; nothing above it changes afterwards.

Stage 1 of a multi-stage build. Stage 2 lets this column's stratospheric water
feed a global stratosphere; it gets its own protocol once this one has a
verdict.

## Why

`model-top.md` found that the global-mean column's middle stratosphere is 21
to 32 K too cold, and that the whole deficit is its water. The column's
stratosphere is saturated at its own cold point: 205.7 K at 100 hPa gives 37.6
ppmv, against the observed 4. The observed 4 ppmv is made where air actually
enters the stratosphere: through the cold tropical tropopause, rising in the
Brewer-Dobson circulation. A global-mean column has no such place. A tropical
one does.

## The column

A reference column representing the band between the turnaround latitudes of
the tropical upwelling, taken as |lat| <= 25 degrees (area fraction 0.4226).
Seviour et al. (2012) give the upwelling mass flux between the turnaround
latitudes but place those latitudes only in figures, so 25 degrees is declared
and the runs are repeated at 20 and 30. It is the column of `model-top.md` as
adopted:

- measured solar spectrum;
- ozone;
- vapour eddy diffusion at 0.39 m^2 s^-1;
- cap alpha = 0.2, with the top at 1 hPa, at its refit;
- no water vapour solar absorption.

Three things change, each from observation or geometry, and none fitted.

1. **The band's sunlight, from geometry.** The annual-mean daily insolation
   averaged over the band, as a fraction of the total irradiance: circular
   orbit, obliquity 23.44 degrees.

   | band | fraction | insolation |
   |---|---:|---:|
   | global | 0.25 exactly | |
   | +/-20 | 0.29979 | |
   | +/-25 | 0.29681 | 404.1 W m^-2 |
   | +/-30 | 0.29328 | |

   NCEP R1's is 406.2 W m^-2 over +/-25. It replaces the orbit module's
   one-quarter.
2. **The band's heat export, as a declared lateral boundary.** In equilibrium a
   band exports its top-of-atmosphere surplus.
   - **Source.** NCEP R1 1991-2020, band-mean net TOA radiation minus its
     global mean: +25.6 - (-12.9) = 38.5 W m^-2 at +/-25 (42.5 at +/-20, 33.0
     at +/-30).
   - **Why the global mean comes off.** NCEP's global net is -12.9 W m^-2 and
     reflects 117 W m^-2 against CERES's 99, so its bias is assumed uniform and
     removed.
   - **In total.** 38.5 W m^-2 is 8.3 PW, or 4.2 PW across each edge at 25
     degrees. Trenberth and Caron (2001) have the total poleward transport
     peaking near 6 PW at 35 degrees.
   - **Where it is taken from.** The atmosphere below 100 hPa, in proportion
     to layer mass. That is a declared approximation: part of the real
     transport is oceanic and leaves from the surface.
3. **The Brewer-Dobson upwelling, as a prescribed vertical mass flux.**
   - **Source.** Seviour et al. (2012): 5.9e9 kg s^-1 upward at 70 hPa
     between the turnaround latitudes (ERA-Interim 1989-2009; 6.2e9 in UKMO
     analyses, 5.9e9 across chemistry-climate models). Per unit area of a
     +/-25 degree band that is 2.737e-5 kg m^-2 s^-1, about 0.23 mm/s at 70
     hPa.
   - **Vertical shape, declared.** Zero below 100 hPa; F0 from 100 to 70 hPa;
     falling linearly in pressure to zero at the 1 hPa top above that.
   - **How it moves things.** Air enters laterally where the flux begins, with
     the column's own properties there, and leaves laterally where it
     decreases, with the local ones. Between, dry static energy c_p T + g z
     and vapour are carried upward in flux form, upwind. The adiabatic cooling
     of rising air is the dry-static-energy transport, not an added term.
   - **Where it goes.** What leaves is recorded as a transfer to a lateral
     boundary module, energy and water both, with the water's latent energy.
     The budgets still close.

The circulation's transport composes with the column's water and convection as
the two-column circulation's did: published as `vapour_advection_rate` and
`advective_heating_rate`.

## The tests

200000 steps, 3600 s, equilibrated as in the earlier protocols. Observations
are the +/-25 band means, NCEP R1 1991-2020 temperatures and SWOOSH v2.73
1991-2020 water, computed here. Temperatures are interpolated in log pressure.

| test | observed | passes if |
|---|---|---|
| T1: runs valid | | every run passes the validity-domain gate |
| T2: cold point temperature | 196.4 K at 100 hPa (band minimum, NCEP levels) | the coldest layer within 5 K of it |
| T2: cold point pressure | 100 hPa | the coldest layer centred between 70 and 150 hPa |
| T3: stratospheric water, mean of the layers between 100 and 68 hPa | 3.45-3.75 ppmv | 1.8 to 7.4 ppmv |
| T4: T at 70 hPa | 199.6 K | within 8 K |
| T4: T at 50 hPa | 206.9 K | within 8 K |
| T4: T at 30 hPa | 215.6 K | within 8 K |
| T4: T at 10 hPa | 230.7 K | within 8 K |

T3's band is a factor of two around the observed range.

**Sensitivity, reported.**

- The same column at +/-20 and +/-30 degrees, each with its own insolation,
  export and upwelling per unit area.
- The +/-25 column with the export taken half from the surface.
- The +/-25 column with no upwelling. This isolates the Brewer-Dobson
  cooling.

A verdict that changes across the latitude range is reported as fragile.

**Reported, not judged:** the surface temperature (the tropical ocean is about
300 K) and the cold point's saturation mixing ratio.

---

## Status (appended; nothing above this line changed)

### Not run: the checks before the runs found a blocker, and its cause

**The code.** The mechanism is built: `physics/tropical_exchange.py` and
`run/tropical_column.py`. Its five unit tests pass, including one that rising
air cools at g / c_p and a 200-step coupled run whose budgets close through the
lateral boundary.

**The blocker.** Short checks before the frozen runs found that no run of this
column can pass T1:

| check | what happened |
|---|---|
| 2000 steps from the default isothermal 250 K start | the top layers reach 355 K, past water's measured 350 K ceiling |
| 5000 steps from the US Standard Atmosphere 1976 | again 355 K; 344 K at 1.6 hPa time-mean |

The second start was declared after the first check, as an implementation
choice (`standard_atmosphere_temperature`). So it is not a spin-up artefact.
The tropical band gets more sunlight, and its stratopause heats well past the
observed 270 K.

**What the same short run already shows** (not a verdict): the coldest layer is
201.7 K at 126 hPa, and the water between 100 and 68 hPa is **3.83 ppmv**,
against the band's observed 3.45-3.75. The global column gave 37.6.

**The cause, measured like for like on the AFGL mid-latitude summer
atmosphere:**

| | model | line-by-line |
|---|---:|---:|
| ozone's ultraviolet absorption, sun at 60 degrees | 13.3 W m^-2 | 14.0 W m^-2 (Chou and Suarez 1999, Table 8) |
| longwave cooling of the stratopause | 5.8 C/day at 1.5 hPa | about 12 C/day stratospheric maximum, all gases (Chou et al. 2001) |

- **The heating is right.** The ultraviolet absorption is within 5 percent.
- **The cooling is half what it should be.** The reference column has no
  15-micron CO2 band. Its infrared absorber outside the vapour is
  `window_absorption(kappa_opaque=1e-4)`, "a deliberately crude caricature":
  a flat 1e-4 m^2 per kg of air everywhere outside 8-12 microns. That alone
  gives 3.9 C/day at 1.5 hPa, and vapour adds 1.9.
- **Ozone's 9.6-micron band** is still absent, as declared.
- **The calibrated CO2 band cannot simply be switched in.** It exists
  (`co2_band_absorption`), but it was fitted on a 64-band grid. On this
  16-band grid, which puts the whole band in three bins, it cools the
  stratopause at 277 C/day.

**What unblocks stage 1:** an infrared absorber that is physics in the middle
atmosphere as well as at the tropopause. That is a resolved 15-micron CO2 band,
benchmarked against line-by-line cooling rates as the solar absorption was, and
then the refits it forces. It is its own mechanism, with its own protocol,
before these runs.
