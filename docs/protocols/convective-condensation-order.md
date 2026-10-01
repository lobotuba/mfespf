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
