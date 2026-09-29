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

---

## Result (appended after the sweep; nothing above this line was changed)

Errors in K against 288 K, every run under 0.1 W m^-2 of imbalance (28000
steps; 44000 at 36 layers), budgets closing to 1e-15. Raw rows:
`convective-fractional-top.results.jsonl`.

| layers | 25 | 30 | 32 | 36 | 40 | 45 | S |
|---|---:|---:|---:|---:|---:|---:|---:|
| whole-layer (reference) | +9.89 | +9.81 | +10.53 | +10.23 | +10.64 | +10.16 | 0.84 |
| fractional top | +11.53 | +11.57 | +11.44 | +11.91 | +11.60 | +11.50 | 0.47 |

**S = 0.47 K > 0.3 K: not adopted.** The whole-layer scheme stays the default
and the fractional top stays available as `convective_fractional_top`, off.

What it did:

1. The scatter fell by 45 percent, from 0.84 to 0.47 K. Five of the six layer
   counts lie within 0.16 K of each other; 36 layers, the one run that needed
   44000 steps to settle, lies 0.3 K above them. The rule counts it, and it is
   not set aside after the fact.
2. The limit-cycle span at 30 layers fell from 3.55 to 2.31 K.
3. The error moved from about +10.3 to about +11.6 K. Why has not been
   measured; the obvious hypothesis is that a group stopping at its neutral
   level mixes less of the air above than one that must take whole layers, so
   the surface keeps more heat.
   As with the relaxation, the closure's form shifts the answer by more than
   the scatter it was built to remove.

The remaining known discontinuity is at onset, where newly unstable layers
first mix part-way into the layer above rather than from zero. Whether that is
what the 36-layer run is showing has not been checked, and would need its own
frozen test before it was used to change anything.

*Reproducibility note.* The sweep ran with the vapour mixing written as
`q * mass`; the committed code writes `water * mass / total`, the same quantity
in the order the whole-layer scheme always used, so that the default path stays
bit-identical. A limit-cycling column amplifies that round-off to a few
thousandths of a kelvin (0.002 K on the headline column), which cannot move
S = 0.47 K across 0.3 K.

