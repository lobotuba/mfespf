# Protocol: choosing the convective relaxation timescale

**Frozen before any sweep was run.** This file is committed ahead of the
implementation and the measurements it governs, so that the choice it produces
cannot drift toward whichever value happens to look best. Editing the rule after
the sweep has run invalidates the choice (the same reasoning as EXP-001's frozen
thresholds).

## Why a timescale at all

The convective adjustment mixes whole layers instantaneously or not at all. Two
measured symptoms follow from that:

1. **In time**, the full water-cycle column limit-cycles, with a surface
   temperature span of 3 to 4 K over the tail of a run.
2. **In resolution**, its error scatters from +9.8 to +10.6 K across 25 to 60
   layers with no trend, following the time-mean top of the mixing at about
   0.015 K per hPa, because that top snaps to a layer edge.

A relaxed adjustment applies the fraction `1 - exp(-dt / tau)` of the
instantaneous adjustment each step. `tau -> 0` recovers the current scheme
exactly, and conservation holds per step for any `tau`, because the adjustment
is scaled by one number within each mixed group.

## The declared range

`tau` from **1 h to 16 h**, nominal **2 h**, basis **published**: Frierson
(2007, J. Atmos. Sci. 64) relaxes over a fixed 2 h timescale, tested 1 h to 16 h,
and found the simulations remarkably insensitive across that range; the 2 h value
descends from Betts (1986, Q. J. R. Meteorol. Soc. 112).

## What is measured

For each `tau` in {1, 2, 4, 8, 16} h, and `tau = 0` as the reference:

- Configuration: the full water-cycle column exactly as
  `tests/acceptance/test_surface_temperature_error.py` builds it (16 bands,
  `areal_heat_capacity = 1e6`, hourly steps, vapour radiatively active and
  transported), at **25, 30, 32, 36, 40 and 45 layers**.
- Each run long enough that its top-of-atmosphere imbalance is under
  0.1 W m^-2; errors are time-means over the last fifth.
- **S(tau)**: the resolution scatter, max minus min of the error over the six
  layer counts.
- Recorded alongside, not used to select: the limit-cycle span at 30 layers, and
  the mean convective top at each layer count.

## The selection rule

1. The published nominal, **2 h, is adopted if S(2 h) <= 0.3 K**. 0.3 K is twice
   the 0.15 K to which a limit-cycling column's time-mean is converged, so a
   scatter below it is indistinguishable from averaging noise.
2. Otherwise, the **smallest** `tau` in the declared range with S(tau) <= 0.3 K.
3. Otherwise **none is adopted**: the instantaneous scheme stays the default,
   and the result is recorded as a finding about the closure.

## What plays no part

The observed surface temperature. S is a difference between runs of the same
model and does not involve 288 K; the error of each run is reported, as every
equilibrated run's is, and is not an input to the choice. Selecting `tau` for
the error it produces would be calibrating against the temperature record,
which INF-003 forbids.
