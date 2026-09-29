# Protocol: what sets the water above convection

**Frozen before the mechanism is written and before any run.** Results are
appended below a line at the end; nothing above it changes afterwards.

## Why

`spectral-wtg-circulation.md` found that the reference column has no single
equilibrium. Nothing in it moves water above the level convection reaches, and
nothing removes water there except condensation. So whatever vapour the
initial state puts above the convective top stays, and it sets how deep
convection goes:

| start | surface |
|---|---:|
| dry (the default) | 291.95 K |
| 1e-6 kg/kg in every layer | 293.30 K |
| moistened by a transient circulation | 296.59 K |

The +3.9 K headline is the dry-start member of this family. Until something
physical sets that water, the headline is not a unique prediction.

## The physics

In the real atmosphere, water reaches the stratosphere by slow mixing and
ascent across the tropopause, and it is freeze-dried at the cold point on the
way. One-dimensional photochemical models represent this with vertical eddy
diffusion of trace gases, K_zz, derived from observed tracer profiles (Massie
and Hunten 1981), plus condensation. In such a model, the water entering the
stratosphere "is set by the tropopause temperature" (Tsai et al. 2021, their
section 3.3).

- **Eddy diffusion of vapour.** The upward vapour flux at each interface is
  F = -rho K dq/dz, from the adjacent layers' specific humidities, the
  interface density, and the hypsometric distance between the layer centres.
  - It is flux form, so it conserves the column's water exactly.
  - It is part of the water module, which owns the vapour. It enters the
    provisional profile before the saturation adjustment, so freeze-drying at
    the cold point is the existing condensation acting on diffused air. It
    reaches convection through `vapour_phase_change_rate`.
- **K = 0.39 m^2 s^-1, uniform.** This is the minimum of Massie and Hunten's
  profile, 3.9e3 cm^2 s^-1 at 17 km, as tabulated by Hu et al. (2012,
  doi:10.1088/0004-637X/761/2/166, their Table 3). That table is for their
  benchmark planets; its N2-atmosphere column is Earth's profile "at the same
  pressure level ... scaled by ... 1.0", so it is Massie and Hunten's own
  values. It is not fitted to anything.
  - Their tropospheric value, 1.0e5 cm^2 s^-1, stands for convective mixing.
    Convection is explicit here, so it is not used. At 0.39 m^2 s^-1,
    diffusion is negligible against convection in the troposphere.
  - The tabulated profile rises above the minimum to 1.1e5 cm^2 s^-1 near
    86 km. So the runs are repeated at 3.9 and 11 m^2 s^-1, which brackets it.
- **Energy.** Diffusion changes no column total of vapour, so the water
  module's latent store is unchanged. The sensible enthalpy the vapour carries
  is neglected, as it already is for precipitation.
- **Stability.** The explicit step needs K dt / dz^2 below 1/2 at every
  interface. The module refuses a configuration that exceeds it rather than
  limiting the flux.

A known failure mode, declared in advance. Hu et al. (2012) report that
one-dimensional transport-and-condensation schemes "may over-predict the
amount of water vapor above the cold trap by a few orders of magnitude". A
global-mean tropopause is much warmer than the tropical cold point that
actually dries the stratosphere: 217 K in the US Standard Atmosphere against
about 200 K needed. They repair it by capping saturation at 20 percent. That
is a correction and is not adopted here. If this column fails for that reason,
the failure is the result, and it points at the tropics.

## The tests

Reference column, 30 layers, 3600 s steps, 200000 steps (22.8 years; the
slowest diffusion mode over the stratosphere at K = 0.39 has an e-folding time
of about 8 years).

- **Equilibrated** means both of the following over the last fifth of the run:
  - net top-of-atmosphere flux within +/-0.1 W m^-2;
  - stratospheric water drifting by less than 5 percent between the first
    and second halves of that fifth.

**Test W1: one equilibrium.** Start the column dry, and start it saturated
(initial relative humidity 1 everywhere).

- **Passes** if the two end within 0.1 K of each other in surface temperature
  and within 10 percent in stratospheric water.
- **Null control:** the same two starts without diffusion.

**Test W2: the stratospheric water it predicts.** Stratospheric water is the
mean mixing ratio of the model layers centred above 100 hPa (with this grid,
at 94, 60 and 27 hPa; the layers have equal mass).

- **Observed:** SWOOSH v2.73 combined, anomaly-filled, 10-degree zonal means,
  1991-2020, area-weighted global mean, computed here. It is 3.85 ppmv at
  100 hPa, 3.57 at 83 (the minimum), 3.91 at 56 and 4.51 at 26. Interpolated
  to the model's layers and averaged, that is about 4.0 ppmv.
- **Passes** if the model is within a factor of two: 2.0 to 8.0 ppmv.
- The model has no methane oxidation, the source of the rise above the
  minimum, so if it has any bias it should be low above 60 hPa.

**Test W3: the headline, re-stated.** The surface-temperature error of the
reference column with diffusion, at the verdict K. It is reported, not judged.

**Sensitivity.** W1 and W2 are repeated at K = 3.9 and 11 m^2 s^-1. A verdict
that changes across that range is reported as fragile.

The observed surface temperature decides nothing here; it is reported.

---

## Results (appended; nothing above this line changed)

Reference column, 30 layers, 3600 s, 200000 steps. Last fifth of each run:

| K (m^2/s) | start | surface | net TOA | stratospheric water (first half -> second half) | PW |
|---:|---|---:|---:|---:|---:|
| 0.39 | dry | 297.934 K | -0.003 | 17.1 ppmv (17.1 -> 17.1) | 69.3 |
| 0.39 | saturated | 297.933 K | +0.002 | 17.1 ppmv (17.1 -> 17.1) | 69.2 |
| 0 (null) | dry | 291.941 K | -0.004 | 0 | 37.4 |
| 0 (null) | saturated | 297.339 K | -0.008 | 12.7 ppmv | 64.9 |
| 3.9 | dry | 298.304 K | 0.000 | 17.4 ppmv (17.4 -> 17.4) | 73.1 |
| 3.9 | saturated | 298.441 K | 0.000 | 17.3 ppmv (17.3 -> 17.3) | 74.0 |
| 11 | -- | refused | | | |

All runs are equilibrated by the protocol's definition. Energy and water closed
to about 1e-15.

At K = 11 the module refused the configuration at step 0: layer 1 would
exchange 1.23 of its mass per 3600 s step. That is the declared refusal, and
the frozen step was not changed to get around it. Separately, the driver's
manifest reports this refusal as "mandatory gate 'conservation_reconciliation'
failed", which it is not; that is recorded in `docs/todo.md`.

### W1: one equilibrium. Passes at the verdict K; fragile across the range

| | surface difference | stratospheric-water difference | verdict |
|---|---:|---:|---|
| K = 0.39 | 0.001 K | 0.0 % | **pass** (0.1 K, 10 %) |
| null, K = 0 | 5.40 K | 0 against 12.7 ppmv | fails, as a null should |
| K = 3.9 | 0.137 K | 0.6 % | fails on temperature |

Diffusion removes the dependence on the initial state, from 5.4 K to 0.001 K.
At K = 3.9 the two starts differ by 0.137 K, beyond the frozen 0.1 K, while
their water agrees to 0.6 percent. Both runs are equilibrated. The
difference is of the size of this column's limit-cycle means (CLAUDE.md), but
that is not measured here, so the verdict stands: **fragile**.

### W2: stratospheric water. Fails, by a factor of 4.3

17.1 ppmv against the observed 4.0 (SWOOSH 1991-2020); the band was 2.0 to 8.0.
This is the declared failure mode, and the profile says why:

| p (hPa) | T (K) | RH | H2O (ppmv) |
|---:|---:|---:|---:|
| 160 | 219.6 | 0.95 | 155 |
| 127 | 208.4 | 1.00 | 44 |
| 94 | 202.1 | 1.00 | 23 |
| 60 | 196.7 | 1.00 | 16 |
| 27 | 190.3 | 0.99 | 12 |

- **No inversion.** Temperature falls all the way to the model top, so the top
  layer is the cold point. Every layer above 130 hPa sits at saturation at its
  own temperature, and the water decreases upward instead of staying at the
  value set at a cold point.
- **Where observation differs.** The tropical cold point is near 190 K at about
  90 hPa, and the air above it warms, so it is subsaturated and carries the
  cold-point value upward. This column's 94 hPa layer is 12 K warmer than that
  cold point.
- **Cause.** The column has no stratospheric temperature inversion, and the
  inversion is made by ozone absorbing sunlight.

That is a missing mechanism, not a transport error.

### W3: the headline, re-stated

The reference column with diffusion at K = 0.39: **+9.93 K** (297.93 K against
288.0), from either start. The +3.94 K quoted until now was the dry-start
member of an initial-condition family. With nothing to carry water above the
convective top, the upper troposphere stayed at exactly zero relative humidity.
With diffusion it carries 0.85 to 0.95 from 160 to 230 hPa, and the vapour
greenhouse there warms the column by 6 K.

The null run started saturated lands at 297.34 K without any diffusion. That
confirms the 291.94 K branch was the artefact, not the physics.
