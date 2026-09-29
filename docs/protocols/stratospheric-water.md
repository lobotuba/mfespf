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
