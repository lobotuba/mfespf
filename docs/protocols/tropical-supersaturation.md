# Protocol: is the tropical column's upper troposphere supersaturated, or only its average?

**Frozen before the run.** Results are appended below a line at the end;
nothing above it changes afterwards.

## Why

`tropical-surface-diagnostic.md` found the tropical column's time-mean vapour
above 320 hPa is 1.24 to 4.4 times saturation at its time-mean temperature
(`docs/todo.md`). Two explanations:

- **Averaging.** Saturation is convex in temperature. A column that cycles
  (the limit cycle in `docs/todo.md`) has a mean saturation above saturation at
  its mean temperature, while every step is at or below saturation.
- **Ordering.** The upwelling's vapour tendency and the water module's
  condensation do not compose so as to remove the excess, and steps end
  supersaturated.

## The run

The +/-25 degree tropical column of `tropical-column.md` attempt 2, rerun
exactly: 200000 steps, the same script. Every step of the last fifth is read
from the run's own retained trajectory: the committed state at the step's end.

The quantity is r = q / q_sat(T, p), per layer per step, using the water
module's own saturation function, at the layers centred between 330 and 150
hPa.

## The tests

**S1, the step-wise maximum.** A layer is supersaturated in fact if r exceeds
1.01 on more than 1 percent of the steps. The 1 percent allows for the step a
cycle's burst crosses.

**S2, the averaging.** Compare the mean of r over the steps with the ratio of
the means. Averaging explains the time-mean supersaturation if two things hold
at each layer:

- the mean of r is at most 1.01;
- the ratio of the means exceeds 1.01.

**The verdict:**

| S1 | S2 | verdict | next |
|---|---|---|---|
| no layer supersaturated | holds | **averaging** | Time-mean humidities above 320 hPa are reported as means of r from now on. Nothing in the code changes. |
| any layer supersaturated | | **ordering** | Find the composition that leaves vapour above saturation, under its own protocol, before the cloud field. |
| no layer supersaturated | fails | **neither** | Reported, and the cloud field goes ahead with the humidity result above 320 hPa marked open. |

**Reported:** the temperature span (5th to 95th percentile) at each layer, and
the fraction of steps on which each layer is saturated, r > 0.99.
