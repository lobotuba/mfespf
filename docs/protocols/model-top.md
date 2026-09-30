# Protocol: a model top above the ozone

**Frozen before any run at the new top.** Results are appended below a line at
the end; nothing above it changes afterwards.

## Why

`stratospheric-resolution.md` resolved the stratosphere and found it
converged, but no refinement could be adopted. The model top is at 10 hPa, and
the ozone above it (50 of 287 dobson units) is placed in the top layer.

- **On today's grid** the top layer is 10 km deep, and that placement is
  harmless.
- **On a fine grid** it puts the heating of the whole upper stratosphere into
  one thin layer. That layer reached 412 K, and the validity gate invalidated
  the run.

There are two ways out:

- **Remove that ozone's sunlight before it enters, declared as absorbed above
  the domain.** That takes energy out of the column's budget which on Earth
  is re-emitted, partly downward. It is a correction.
- **Raise the model top above the ozone,** so its heating is deposited where
  it happens. That is the physics, and it is what is tested here.

## The change

**The top.** The model top moves from 10 hPa to 1 hPa, where the SWOOSH profile
ends. The ozone above 1 hPa, about one dobson unit, is still placed in the top
layer, now as declared a sliver.

**The grid.** The same log-pressure cap as `stratospheric-resolution.md`. The
reference step becomes (p_s - p_top) / 30 = 33.74 hPa, up from 33.41, so the
lower troposphere shifts by 1 percent; each grid is refitted anyway.

| alpha | layers, top 1 hPa |
|---:|---:|
| 0.4 | about 39 |
| 0.2 | about 51 |
| 0.1 | about 78 |

**The calibration.** Refitted on each grid by the unchanged method. The
prescribed profile is held at 220 K above 26.7 hPa, as that protocol stated it,
so it is isothermal to 1 hPa.

Nothing else changes. The column is attempt 2 of `ozone-solar-heating.md`:

- the reference column;
- the measured solar spectrum;
- ozone;
- vapour eddy diffusion at 0.39 m^2 s^-1;
- no water vapour solar absorption.

It runs 3600 s steps for 200000 steps, equilibrated as before. Temperatures at
given pressures are interpolated in log pressure.

## The tests

**M1: the runs are valid.** Every run at the new top must pass the mandatory
validity-domain gate. If one does not, it supports no claim, as before.

**M2: convergence chooses alpha.** Exactly the rule of
`stratospheric-resolution.md`: compare 0.4 with 0.2 and 0.2 with 0.1, and adopt
the coarsest alpha within the limits below of the next finer one. "None" is
allowed.

| quantity | limit |
|---|---|
| each temperature below | 1 K |
| the pressure of the coldest layer | one layer's depth |
| stratospheric water | 10 percent |

**M3: the stratosphere, at the adopted alpha.** O1 to O3 of
`ozone-solar-heating.md` unchanged, plus 10 hPa, which the new top puts inside
the domain:

| test | observed | passes if |
|---|---|---|
| T at 60 hPa | 209.0 K (NCEP R1 1991-2020 global mean) | within 8 K |
| T at 27 hPa | 218.0 K | within 8 K |
| T at 10 hPa | 228.8 K | within 8 K |
| coldest layer | global-mean minimum at 100 hPa | centred between 60 and 160 hPa |
| stratospheric water, mass-weighted mean between 100 and 10 hPa | about 4.0 ppmv (SWOOSH) | 2 to 8 ppmv |

The water is averaged between 100 and 10 hPa, not to the new top. SWOOSH rises
to 4.8 ppmv at 10 hPa through methane oxidation, which the model does not have,
and the comparison band is kept where it was.

Reported, not judged:

- the temperature at 1-2 hPa, against the US Standard Atmosphere's stratopause
  of about 270 K;
- the surface-temperature error;
- the refitted wing decays.

The observed surface temperature decides nothing here.
