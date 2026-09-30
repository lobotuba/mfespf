# Protocol: resolving the stratosphere

**Frozen before the grid option is written and before any run.** Results are
appended below a line at the end; nothing above it changes afterwards.

## Why

Ozone takes 11 W m^-2 of sunlight and warms every stratospheric layer by 7 to 9
K, yet no inversion forms (`ozone-solar-heating.md`, attempt 2).

- **What failed.** 27 hPa sits 19 K below the observed 218, the coldest layer
  is the top, and the stratosphere holds 60 ppmv of water against 4.
- **The grid.** The column has thirty equal-mass layers, so three lie above
  100 hPa. Its layers are 0.3 km thick at the surface and 10 km at the top
  (43 to 10 hPa), and that top layer holds 177 of the 287 dobson units.
- **Why that matters.** The observed inversion spans 100 to 10 hPa, about 16 km.
  It cannot form in a layer that thick. By this project's own rule, a grid must
  resolve what the physics depends on, and this one does not resolve the
  stratosphere.

## The grid

**The rule.** Layers keep the calibration grid's pressure step, (p_s - p_top) /
30 = 33.4 hPa, wherever that is thinner than a declared maximum depth in log
pressure. Above that, each layer spans exactly that maximum:

- the step is dp_i = min(33.4 hPa, p_i (1 - exp(-alpha))), built upward from the
  surface;
- a remainder under half a step at the top is merged into the last layer;
- alpha times the scale height (about 7 km) is the thickest a layer may be.

**Alpha is a resolution, not a physical parameter.** It is chosen by
convergence below, like every resolution here (VAL-002). Omitting it gives
today's grid exactly, so nothing that does not ask for it changes.

| alpha | layers | above 100 hPa | thickest layer | first interface that moves |
|---:|---:|---:|---:|---:|
| none (today) | 30 | 3 | 10.3 km | -- |
| 0.4 | 33 | 6 | 3.1 km | 51.5 hPa |
| 0.2 | 39 | 11 | 1.9 km | 145 hPa |
| 0.1 | 55 | 23 | 1.0 km | 312 hPa |

## The calibration moves with it

The vapour wing decay (outgoing longwave 238.2 W m^-2 on a prescribed profile)
and the vapour-masked CO2 wing decay (3.70 W m^-2 per doubling) are fitted
jointly on the grid the runs use. By the rule in CLAUDE.md, the refit is part of
the change.

- **The method is unchanged.** Same targets, same bisection, same prescribed
  humidity.
- **One thing is restated.** The prescribed temperature profile was defined per
  layer index, 288 K in the lowest layer to 220 K in the top one, linear. On
  equal-mass layers that is linear in pressure. It is now stated in pressure:
  - linear from 288 K at the calibration grid's lowest layer centre (996.6 hPa)
    to 220 K at its top centre (26.7 hPa);
  - held at 220 K above that.

  On today's grid this reproduces today's fit exactly.
- **Each alpha gets its own refit**, recorded with its run, so every grid
  compared runs at its own calibration.

## The tests

**Column.** Attempt 2 of `ozone-solar-heating.md`: the reference column with
the measured solar spectrum, ozone and vapour eddy diffusion at 0.39 m^2 s^-1,
and without water vapour's solar absorption. 3600 s steps, 200000 steps,
equilibrated as there. Temperatures at 60 and 27 hPa are interpolated in log
pressure between layer centres, since the layers no longer sit there.

**R1: convergence, which chooses alpha.**

- Run today's grid and alpha = 0.4, 0.2 and 0.1, each refitted.
- Compare each alpha with the next finer one on:
  - the temperatures at 60 and 27 hPa;
  - the pressure of the coldest layer;
  - the surface temperature;
  - the stratospheric water (mass-weighted mean above 100 hPa).
- **Adopt** the coarsest alpha whose results differ from the next finer one's by
  at most 1 K in each temperature, at most one layer's depth in the cold point,
  and at most 10 percent in the water.
- If 0.2 against 0.1 does not meet that, **none is adopted** and that is the
  result.

**R2: ozone's test, at the adopted alpha.** O1, O2 and O3 of
`ozone-solar-heating.md`, same thresholds:

| test | passes if |
|---|---|
| O1: temperature at 60 hPa | 209.0 +/- 8 K |
| O1: temperature at 27 hPa | 218.0 +/- 8 K |
| O2: coldest layer | centred between 60 and 160 hPa |
| O3: stratospheric water | 2 to 8 ppmv |

**R3: the calibration and the headline, re-stated.** Reported, not judged, at
the adopted alpha:

- the two refitted wing decays;
- the default reference column's surface-temperature error, with and without
  vapour diffusion.

**Declared in advance.** The ozone above 10 hPa is still placed in the top
layer, which is now thin. So that layer takes heating that on Earth is spread
over the upper stratosphere, and may run hot. It is reported, and the 27 hPa
test is evaluated below it.

The observed surface temperature decides nothing here.
