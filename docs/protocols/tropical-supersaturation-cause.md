# Protocol: which composition leaves the tropical upper troposphere supersaturated

**Frozen before the run.** Results are appended below a line at the end;
nothing above it changes afterwards.

## Why

`tropical-supersaturation.md` found that the tropical column's steps end with
vapour above saturation from 322 to 154 hPa, up to 300 times saturation. The
verdict was ordering.

**The hypothesis, from the code.** Three things happen in each step:

1. The water module condenses everything above saturation in its
   *provisional* profile: the committed vapour plus upwelling, diffusion and
   evaporation.
2. Convection then mixes the vapour to a uniform specific humidity over each
   group of layers it overturns, starting from the profile the water module
   left (`convection.py`, `_mix_vapour`).
3. Nothing condenses what that mixing lifts until the next step. So a layer at
   the top of a group ends the step holding the group's mean specific humidity,
   far above its own saturation.

`convection.py`'s declared failure modes say this in words: "mixing completely
instead drives the top of the group past saturation, condensation clips it
there". The clip lands one step late.

**The alternatives.**

- The upwelling's vapour advection is composed *before* the adjustment
  (`water.py`), so it should leave nothing.
- So should eddy diffusion and evaporation.

## The run

The +/-25 degree tropical column of `tropical-column.md` attempt 2, exactly as
in `tropical-supersaturation.md`, 200000 steps.

From each step of the last fifth, read:

- the committed temperature and vapour;
- the layers convection mixed: the finite entries of
  `convectively_mixed_layer_temperature`.

The layers judged are those centred between 330 and 150 hPa. r = q / q_sat(T, p),
as before.

## The test

Pool every layer-step, over the judged layers, that ends with r > 1.01.

- **X1, attribution.** The fraction of those that convection mixed on the same
  step.
  - If X1 >= 0.95, convection's vapour mixing after the saturation adjustment
    is the composition, and the fix is designed under its own protocol.
  - If X1 < 0.95, the remainder is examined next, starting with the
    upwelling's advection.

**Reported:**

- the converse: the fraction of mixed layer-steps that end supersaturated;
- X1 per layer;
- whether any layer-step with r > 1.01 has no mixing on that step or the one
  before it.
