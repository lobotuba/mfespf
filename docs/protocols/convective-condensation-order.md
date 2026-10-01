# Protocol: condensing what convection lifts, in the step it lifts it

**Frozen before the implementation and before any run.** Results are appended
below a line at the end; nothing above it changes afterwards.

## Why

`tropical-supersaturation-cause.md` found that every supersaturated
layer-step in the tropical column (63166 of them, up to 300 times saturation)
is a layer convection mixed on that step.

Within a step, the water module condenses its provisional vapour first.
Convection then mixes the vapour to a uniform specific humidity over each
group. The top of a group ends the step holding the group's mean humidity.
Nothing condenses it until the next step. Meanwhile radiation, which runs
first, sees the supersaturated vapour every step.

This is a composition defect of the kind `CLAUDE.md` names: "compose split
operators, do not cap them". The upwelling already composes correctly:

- it publishes `vapour_advection_rate`;
- the water module adds it to the profile it adjusts;
- so nothing the upwelling carries survives above saturation.

## The change

A configuration option, `condense_after_convection`, off by default.

- **Convection runs before the water module.**
  - It mixes the committed vapour plus the circulation's advection.
  - It does not read the water module's `vapour_phase_change_rate` or
    `latent_heating_rate`, because those are not yet formed.
  - It publishes its mixing as `convective_vapour_transport_rate`, and still
    submits it as its own vapour contribution, as before.
- **The water module runs after convection.** It adds that rate to the
  provisional profile before its saturation adjustment, exactly as it adds
  `vapour_advection_rate`, and excludes it from its own tendency. So what
  convection lifts above saturation condenses in the same step:
  - its latent heat is released in the layer where it condenses;
  - it precipitates by the ordinary route.

**What is given up, declared.** Convection no longer sees this step's latent
heating when it forms its provisional temperature. Condensation is applied
after the mixing, as in a real plume, where condensation happens as the air
rises.

The configuration's comment "water releases latent heat before convection" is
the reason for the present order. Under this option it is replaced by the
opposite argument, which is stated here: the vapour convection lifts is what
condenses, so condensation has to follow the lifting.

The saturation is still taken at the start-of-step temperature, as now.

## The tests

- **O1, unit.** A convecting group whose mixing lifts vapour above a cold
  layer's saturation, in one step under the option:
  - that layer ends at its saturation to 1e-12 relative;
  - the excess precipitates;
  - its latent heat appears in that layer;
  - water and energy both close.
- **O2, the tropical column.** The +/-25 degree attempt-2 column with the
  option on, 200000 steps. In the last fifth, the layers centred between 330
  and 150 hPa end with r > 1.01 on at most 1 percent of layer-steps, the
  threshold `tropical-supersaturation.md` used. The run must be valid.
- **O3, defaults unchanged.** With the option off, the full test suite passes
  unchanged, and the recorded headline numbers do not move.

**Reported, not judged:**

- **O2's run:** the maximum r, the cold point, T at 70/50/30/10 hPa, the water
  at 100 to 68 hPa, and the surface temperature.
- **The global columns, with the option on:**
  - `reference_column` at its default 30 layers. Its headline error is +3.9 K
    off, and with the option on it moves.
  - the 51-layer CO2-band column of `co2-band.md`.

**Adoption** as the default is a decision for after, because it moves every
recorded result.

---

## Results (appended; nothing above this line changed)

**Built** in `9fe754c`:

- `ConvectiveAdjustmentModule(condensation_follows=True)` mixes the committed
  vapour and the advection, and publishes `convective_vapour_transport_rate`.
- `WaterModule(receives_convective_transport=True)` composes that rate before
  its saturation adjustment.
- `reference_column(condense_after_convection=True)` runs convection before the
  water module and refuses a configuration without the water module, convection
  and vapour transport.
- **A defect found on the way.** `tropical_column` spliced its exchange in by
  the order its schedules had been *appended* in, not by the column's
  splitting order, so it silently restored water-before-convection. It now
  uses the splitting order. The default order is unchanged, which is tested.

| test | result | verdict |
|---|---|---|
| **O1**, unit | The lifted layer ends at its saturation to 1e-12, the condensation equals the excess exactly, the latent heat is in the condensing layers only, and water closes to 1e-12. A 24-step coupled run with the option closes its budgets. | **pass** |
| **O2**, tropical column | 11369 of 240000 judged layer-steps (4.74%) end above r = 1.01, against the 1% threshold. Max r is 1.056; it was 303 in the default order. | **fail** |
| **O3**, defaults unchanged | The full suite passes with the headline at +3.944 K, as recorded. | **pass** |

**Why O2 fails.**

- **The excess that is left is small.** It is at most 5.6 percent above
  saturation, where the default order left up to 300 times.
- **Most likely, it is the declared approximation. This is not measured.**
  The saturation is taken at the start-of-step temperature, and a layer that
  cools within the step ends it above saturation at its new temperature. A
  5.6 percent excess is about 0.5 K of cooling at these temperatures.
- **It sits where convection mixes most often.** The 321 hPa layer, mixed on
  89 percent of steps, carries 20.7 percent of its steps above 1.01; 288 hPa
  6.4 percent, 254 hPa 1.3 percent, and those above 220 hPa none.
- **The threshold stands, and the test fails as frozen.** Closing the rest
  needs a saturation taken at the end-of-step temperature, which is a further
  composition, not a tolerance.

**Reported.**

| | default order | condensing after convection |
|---|---:|---:|
| **30-layer headline column** (`reference_column`, 24000 steps) | | |
| surface error | +3.944 K | **+3.538 K** |
| layer-steps with r > 1.01 (whole column, last fifth) | 4.8% | 1.4% |
| max r | 20.3 | 1.089 |
| precipitable water | 37.3 kg m^-2 | 35.1 kg m^-2 |
| **51-layer CO2-band global column** (`co2-band.md`, 200000 steps) | | |
| surface error | -9.23 K | -10.25 K |
| T at 60 / 27 / 10 / 1.5 hPa | 189.0 / 203.4 / 227.6 / 284.4 K | 189.1 / 203.5 / 227.7 / 284.5 K |
| coldest layer | 184.8 K, 103 hPa | 184.7 K, 126 hPa |
| stratospheric water | 1.00 ppmv | 0.97 ppmv |
| **+/-25 tropical column** (O2's run) | | |
| surface | 282.16 K | 279.75 K |
| cold point | 166.7 K, 103 hPa | 163.3 K, 103 hPa |
| T at 70 / 50 / 30 / 10 hPa | 175.3 / 182.7 / 197.3 / 231.3 K | 173.1 / 181.1 / 196.4 / 231.0 K |
| water, 100-68 hPa | about 0 | about 0 |

- **The upper troposphere is now nearly saturated, not supersaturated.**
  Removing the forbidden vapour removes some greenhouse effect, and every
  column cools: by 0.4 K (30-layer), 1.0 K (51-layer global) and 2.4 K
  (tropical).
- **The tropical results do not change their verdicts.** The tropical
  column's failures of `tropical-column.md` are unchanged in kind, and its
  cold point is 3.4 K further from the observed.
- **The tropical surface diagnostic's H term overstated the column's humidity
  above 320 hPa.** Its direction stands below 320 hPa.

**Not adopted.** O2 failed. Adoption was in any case a decision for after.
