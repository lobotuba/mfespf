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

---

## Result (appended after the sweep; nothing above this line was changed)

Errors in K against 288 K, time-means over the last fifth, every run under
0.1 W m^-2 of imbalance (28000 steps, or 44000 where 28000 was not enough).
Raw rows: `convective-relaxation-timescale.results.jsonl`.

| tau | 25 | 30 | 32 | 36 | 40 | 45 | S(tau) | span at 30 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 (instantaneous) | +9.89 | +9.81 | +10.53 | +10.23 | +10.64 | +10.16 | 0.84 | 3.55 |
| 1 h | +12.07 | +12.62 | +12.23 | +11.69 | +11.47 | +12.10 | 1.15 | 4.89 |
| 2 h | +12.13 | +13.04 | +12.84 | +13.16 | +13.59 | +12.89 | 1.46 | 3.93 |
| 4 h | +13.73 | +13.76 | +14.62 | +14.33 | +14.10 | +14.85 | 1.12 | 3.52 |
| 8 h | +14.22 | +14.74 | +14.71 | +15.48 | +16.22 | +16.42 | 2.21 | 1.80 |
| 16 h | +15.96 | +16.48 | +16.48 | +16.55 | +17.47 | +17.88 | 1.92 | 0.79 |

**S(2 h) = 1.46 K > 0.3 K, and no tau in the declared range meets 0.3 K.
Rule 3 applies: none is adopted, and the instantaneous scheme stays the
default.**

What the sweep found, recorded as the rule requires:

1. Relaxation does not reduce the resolution scatter; it increases it, from
   0.84 K to between 1.1 and 2.2 K. The scatter is not caused by the mixing
   being instantaneous. The unstable groups are still whole layers at every
   tau, which is the limitation the equation record declared in advance.
2. Relaxation does damp the limit cycle in time, from 3.55 K to 0.79 K at
   16 h -- the first of the two symptoms responds to a timescale and the
   second does not. At 8 and 16 h the error also trends upward with layer
   count rather than scattering, which is not convergence either.
3. Relaxation warms the column monotonically with tau, from +10 K to +17 K,
   most plausibly because a column only partly adjusted each step carries a standing
   instability and removes heat from the surface less efficiently. Within the
   published range the closure choice alone moves the error by 2 to 8 K -- the
   largest model-form uncertainty measured in this project.

A 1 h run at 30 layers needed 44000 steps to settle and swung 4.89 K, more
than the instantaneous scheme; a short timescale is not a gentle version of
the instantaneous one.

