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

---

## Results (appended; nothing above this line changed)

**Refits.** Both targets held to 0.001 W m^-2.

| alpha | layers | vapour wing (m^-1) | CO2 wing with vapour (m^-1) |
|---:|---:|---:|---:|
| 0.4 | 39 | 3655.66 | 1270.201 |
| 0.2 | 51 | 3652.99 | 1270.192 |
| 0.1 | 78 | 3652.14 | 1270.183 |

The vapour wing is 4.5 percent narrower than at the 10 hPa top, because the
calibration profile now carries an isothermal 220 K layer up to 1 hPa.

**Runs.** 200000 steps, last fifth. All are equilibrated (net TOA +0.013 to
+0.016) and closed to 1e-15.

| alpha | T 60 | T 27 | T 10 | T 1.5 | coldest layer | water 100-10 hPa | valid |
|---:|---:|---:|---:|---:|---|---:|---|
| 0.4 | 203.86 | 196.71 | 197.00 | 266.7 | 17.2 hPa (20.6-13.8) | 52.86 | no |
| 0.2 | 203.43 | 196.66 | 196.96 | 268.4 | 17.1 hPa (18.8-15.4) | 45.85 | no |
| 0.1 | 203.24 | 196.59 | 196.99 | 269.0 | 16.0 hPa (16.8-15.2) | 48.73 | no |

### M1: fails. Every run is invalid, narrowly

- **What failed.** In each run, the validity gate found water's two relations
  (Clausius-Clapeyron, and the latent energy of vapour, both declared for 170
  to 340 K) evaluated up to 341.6-342.4 K.
- **Where and when.** It is the top layers, during spin-up from the 250 K
  isothermal start. Their time-mean is 270-276 K.
- **Not changed.** A run that leaves a declared domain supports no claim. The
  bound is not moved here: moving it needs the relation measured against an
  independent reference above 340 K first, and that is its own piece of work.

### M2: the numbers converge at 0.2, but nothing is adopted

- **0.2 against 0.1:**
  - the temperatures differ by 0.19, 0.07 and 0.03 K;
  - the coldest layer at 0.1 (16.0 hPa) lies inside the 0.2 one (18.8-15.4);
  - the water differs by 6.3 percent.

  All within the rule.
- **0.4 against 0.2:** the water differs by 15 percent, outside it.
- **Result.** By the numbers the rule selects 0.2. With M1 failed, no alpha is
  adopted, and M3 is not a verdict.

### What the converged, invalid runs show

- **Ozone placed where it is makes the stratopause.** The temperature rises from
  197 K at 10 hPa to 268-270 K at 1.5 hPa, against the US Standard
  Atmosphere's stratopause of about 270 K. Nothing was fitted to it.
- **What is left is the middle stratosphere.** Against NCEP R1:

  | pressure | model | observed | difference |
  |---:|---:|---:|---:|
  | 60 hPa | 203.4 K | 209.0 K | -5.6 K |
  | 27 hPa | 196.7 K | 218.0 K | -21 K |
  | 10 hPa | 197.0 K | 228.8 K | -32 K |

  The cold point sits at 17 hPa, not 100.
- **That is where the water sits.** The layers are saturated, and the
  stratosphere holds 46 to 49 ppmv of water against 4. The water's longwave
  cooling is now the leading candidate for the deficit, and it can be
  measured directly: the radiative heating of the converged state, with the
  stratospheric water set to the observed amount.

### Diagnostic after the verdict: the deficit is the water, and so is the cause

**How it was measured.** The alpha = 0.2 run was repeated (200000 steps,
invalid for claims as above) and its last-fifth mean state kept. The model's
own radiation module was evaluated on that state twice: as it is, and with the
water above 100 hPa set to the observed 4 ppmv. The stratosphere was then
adjusted under fixed dynamical heating (Fels et al. 1980; Forster and Shine
1997): every non-radiative term is held at its equilibrium value, and each
layer's temperature relaxes until its radiative heating is what it was. The
adjustment converged to 1e-4 K/day in 323 iterations.

| pressure | model | with 4 ppmv, adjusted | observed (NCEP R1) |
|---:|---:|---:|---:|
| 60 hPa | 203.4 K | 214.8 K | 209.0 K |
| 27 hPa | 196.7 K | 220.8 K | 218.0 K |
| 10 hPa | 197.0 K | 236.3 K | 228.8 K |
| 1.5 hPa | 267.1 K | 322.7 K | about 270 K |

**The middle stratosphere's whole deficit is its water.**

- **With the observed amount, all three NCEP levels are within 8 K,** the
  thresholds of M3.
- **The coldest layer moves** from 17 hPa to about 69 hPa, inside the 60 to 160
  hPa band.
- **The upper stratosphere then runs about 50 K hot.** The excess water was
  hiding too little cooling near the stratopause, which is the next thing
  behind this one.

**Why the column is that wet.** It is saturated at its own cold point. The
model's saturation function gives these mixing ratios:

| cold point | saturation mixing ratio |
|---|---:|
| global mean at 100 hPa, 205.7 K (NCEP R1) | 37.6 ppmv |
| tropical mean at 100 hPa, 195.1 K | 7.3 ppmv |
| 190 K at 90 hPa, illustrative of the seasonal tropical minimum | 3.5 ppmv |

- **The column carries 34 to 44 ppmv**, which is the global-mean value.
- **What makes the observed 4.** Air enters the stratosphere through the coldest
  tropical tropopause, and a horizontally uniform column has no such place.
- **A single global-mean column cannot produce a 4 ppmv stratosphere,** however
  good its radiation. That is the failure mode Hu et al. (2012) declared, now
  measured in this model.

**The mechanism it points to is regional:** a tropical column whose cold point
feeds the global stratosphere. The region framework (`run/regions.py`) exists
for it.
