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

---

## Results (appended; nothing above this line changed)

The run is complete and reproduces attempt 2: coldest layer 166.70 K. In the
last fifth, 63166 layer-steps end with r > 1.01 at the judged layers.

**X1 = 1.000.** Every one of them is a layer convection mixed on that same
step, at every judged layer. None has no mixing on that step or the one
before.

| layer | X1 | steps mixed |
|---:|---:|---:|
| 321.5 hPa | 1.000 | 100% |
| 287.8 hPa | 1.000 | 99.7% |
| 254.1 hPa | 1.000 | 29.8% |
| 220.3 hPa | 1.000 | 4.9% |
| 186.6 hPa | 1.000 | 8.7% |
| 154.3 hPa | 1.000 | 8.7% |

**The converse:** 62.7 percent of mixed layer-steps end supersaturated.

**Verdict.** The supersaturation is convection's vapour mixing after the
saturation adjustment. Convection mixes the group to a uniform specific
humidity, the top of the group ends far above its own saturation, and nothing
condenses it until the next step. The upwelling, diffusion and evaporation
leave none.

Convection reaches 154 hPa on 8.7 percent of steps: deep convection to the
tropical tropopause. The fix is designed under its own protocol.
