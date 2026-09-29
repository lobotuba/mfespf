# Protocol: a convective top that is not a layer edge

**Frozen before the implementation and before any run.** Result is appended
below a line at the end; nothing above it changes after the sweep.

## The defect

The full water-cycle column's error scatters by 0.84 K across 25 to 45 layers
and follows the time-mean top of the convective mixing at about 0.015 K per hPa.
The adjustment mixes whole layers, so that top is a layer edge. A relaxation
timescale made the scatter worse at every value in the published range
(`convective-relaxation-timescale.md`): the defect is spatial.

## The change

Where a mixed group stops growing, it penetrates the next layer by a mass
fraction `f` in [0, 1]. Inside that layer the stability variable is taken as
linear in mass, with the gradient it has to the layer above, so that the layer
mean is unchanged; the mixed part is the fraction below the level where the
reconstructed profile equals the group's mixed value. `f = 1` is the whole layer
and continues into the next. The mixed part joins the group's enthalpy budget and
the unmixed remainder keeps its reconstructed mean, so enthalpy is conserved
exactly and the remainder is more stable than the group by construction. Vapour
is carried in the same mass fraction.

It has **no free parameter**, so there is nothing to choose, only a test to
pass. The instantaneous whole-layer scheme is the `f = 0` case.

## The test

The configuration and layer counts of `convective-relaxation-timescale.md`:
the full water-cycle column, 16 bands, `areal_heat_capacity = 1e6`, hourly
steps, vapour radiatively active and transported, at 25, 30, 32, 36, 40 and 45
layers, each run to under 0.1 W m^-2 of imbalance, time-means over the last
fifth.

- **Adopted as the default if S <= 0.3 K**, S being max minus min of the error
  over the six layer counts. Otherwise not adopted, and the finding recorded.
- Recorded alongside and not used to decide: the limit-cycle span at 30 layers,
  the mean convective top, and each run's error.

The observed surface temperature plays no part: S is a difference between runs.
