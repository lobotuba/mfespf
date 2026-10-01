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

---

## Results (appended; nothing above this line changed)

**M1: the hypothesis holds.** The O2 run was repeated exactly: 11369
layer-steps above 1.01, as before.

- Every one of them is saturated at its start-of-step temperature. The
  largest r_start is 1 + 9e-16, so none exceeds 1 + 1e-9.
- These layers cooled within the step by 0.37, 0.10 and 0.07 K at the 5th,
  50th and 95th percentiles.

So the excess is entirely the within-step temperature change, and Part 2 was
built (`3a5903a`).

| test | result | verdict |
|---|---|---|
| **E1**, unit | Every condensing layer ends saturated at its end-of-step temperature, latent heating included, to 1e-10. Water closes to 1e-12, and the latent energy is L_eff times the condensate to 1e-12. | **pass** |
| **E2**, tropical column, both options | 0 of 240000 judged layer-steps end above r = 1.01, and the maximum r is 1.000. The run is valid, and its budgets close to 1e-15. | **pass** |
| **E3**, defaults unchanged | The full suite passes, and the headline is +3.944 K. | **pass** |

**Reported.** All columns at 200000 steps except the 30-layer headline, at
24000:

| | default | condensing after convection | + end-of-step saturation |
|---|---:|---:|---:|
| **30-layer headline column** | | | |
| surface error | +3.944 K | +3.538 K | **+5.078 K** |
| layer-steps r > 1.01 (whole column) | 4.8% | 1.4% | 0.0% |
| max r | 20.3 | 1.089 | 1 + 1.3e-12 |
| precipitable water | 37.3 | 35.1 | 42.5 kg m^-2 |
| **51-layer CO2-band global column** | | | |
| surface error | -9.23 K | -10.25 K | **-8.74 K** |
| T at 60 / 27 / 10 / 1.5 hPa | 189.0 / 203.4 / 227.6 / 284.4 | 189.1 / 203.5 / 227.7 / 284.5 | 189.0 / 203.3 / 227.5 / 284.4 K |
| coldest layer | 184.8 K, 103 hPa | 184.7 K, 126 hPa | 185.0 K, 103 hPa |
| stratospheric water | 1.00 | 0.97 | 1.05 ppmv |
| **+/-25 tropical column** | | | |
| surface | 282.16 K | 279.75 K | 280.91 K |
| cold point | 166.7 K | 163.3 K | 165.0 K |
| T at 70 / 50 / 30 / 10 hPa | 175.3 / 182.7 / 197.3 / 231.3 | 173.1 / 181.1 / 196.4 / 231.0 | 174.2 / 181.9 / 196.8 / 231.1 K |
| water, 100-68 hPa | about 0 | about 0 | about 0 |

- **No layer ends a step supersaturated in the two columns where it was
  measured,** the 30-layer column and the tropical one. The 51-layer global
  run did not record r.
- **The columns are wetter and warmer than with the first option alone.** The
  condensate's own latent heat raises the saturation it condenses to, and the
  end-of-step temperature replaces the start-of-step one. The 30-layer column
  holds 42.5 kg m^-2 against 35.1 and warms 1.5 K. Precipitable water is
  observed at about 25 (`observed-cloud-field.md`), so this is further from
  it.
- **The tropical column's verdicts do not change.** Its lower stratosphere is
  still 25 to 30 K too cold, and its cold point still lies below the
  saturation relation's 170 K floor.

**Not adopted.** Adopting either option moves every recorded result, so it is
a decision for after.
