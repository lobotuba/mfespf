# Protocol: saturation at the temperature the step ends at

**Frozen before the measurement, the implementation and any run.** Results are
appended below a line at the end; nothing above it changes afterwards.

## Why

`convective-condensation-order.md` cut the tropical column's
upper-tropospheric supersaturation from 303 times saturation to 1.056 times.
It still failed O2: 4.74 percent of judged layer-steps ended above 1.01,
against a 1 percent threshold.

**The hypothesis**, not measured there:

- The water module takes saturation at the *start-of-step* temperature.
- A layer that cools within the step (radiation, the upwelling's adiabatic
  cooling, convection) ends it above saturation at its new temperature.
- A 5.6 percent excess is about 0.4 K at these temperatures.

There is a second omission of the same kind. The latent heat of what
condenses warms the layer, and that raises its saturation. A saturation
adjustment that ignores it condenses slightly too much. The classical moist
saturation adjustment solves for the condensate c that leaves the layer
exactly saturated after its own latent heating:

    q - c / m = q_s(T* + L_eff c / C, p)

where:

- q is the provisional specific humidity;
- T* is the temperature the step would otherwise end at;
- L_eff is the latent heat the module already releases: vaporisation, plus
  fusion where the vapour deposits as ice;
- C is the layer's heat capacity.

## Part 1, the measurement (M1)

**The run.** The +/-25 degree tropical column with `condense_after_convection`,
exactly as O2, 200000 steps.

**The quantity.** For each judged layer-step in the last fifth (layers centred
330 to 150 hPa), r_start = q_end / q_s(T_start, p). Here T_start is the
committed temperature of the step before, which is the temperature the water
module used.

**The verdict.** If r_start <= 1 + 1e-9 on every judged layer-step that ends
with r > 1.01, the excess is entirely the within-step temperature change, and
Part 2 is built. If not, something else adds vapour after the adjustment. That
is reported, and Part 2 waits.

## Part 2, the change

An option, `end_of_step_saturation`, only together with
`condense_after_convection`, and off by default. Running after radiation, the
exchange and convection, the water module forms

    T* = T_start + dt / C * (Q_rad + Q_adv + H_0) + dt * (dT/dt)_conv

from their published fields:

- `layer_radiative_heating`;
- `advective_heating_rate`, when the column has an exchange;
- `surface_sensible_heat_flux`, into the lowest layer;
- `convective_temperature_tendency`.

It then solves the equation above in each layer where the provisional vapour
exceeds q_s(T*), by Newton iteration to 1e-12 relative. Nothing else changes:

- the condensate leaves as precipitation, by the ordinary route;
- its latent heat enters as the module's existing tendency;
- the saturation is still not evaluated below its 170 K floor, so a layer
  there condenses everything, as now.

## The tests

- **E1, unit.** In one step, with a prescribed cooling and a convective
  lift:
  - every condensing layer ends exactly saturated at its end-of-step
    temperature, which includes the latent heating, to 1e-10 relative;
  - water closes to 1e-12;
  - the latent energy released equals L_eff times the condensate.
- **E2, the tropical column.** The O2 run with both options. At most 1
  percent of judged layer-steps end above r = 1.01, the same threshold, and
  the run is valid.
- **E3, defaults unchanged.** The full suite passes, and the headline is
  +3.944 K.

**Reported:** E2's maximum r, and the same columns
`convective-condensation-order.md` reported, with both options on.

**Adoption** of either option as the default is a decision for after.
