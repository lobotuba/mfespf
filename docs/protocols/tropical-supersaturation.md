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

---

## Results (appended; nothing above this line changed)

The run is complete and reproduces attempt 2: coldest layer 166.70 K, surface
282.16 K. The last fifth is 40000 steps.

| layer | steps r > 1.01 | steps r > 0.99 | max r | mean of r | ratio of means | T 5th-95th percentile |
|---:|---:|---:|---:|---:|---:|---|
| 321.5 hPa | 15.7% | 15.7% | 12.3 | 1.26 | 1.24 | 216.4-217.4 K |
| 287.8 hPa | **94.5%** | 94.8% | 27.4 | 2.81 | 2.76 | 209.7-210.7 K |
| 254.1 hPa | 29.8% | 100% | 69.7 | 4.51 | 4.42 | 202.4-203.3 K |
| 220.3 hPa | 4.4% | 88.7% | 14.7 | 1.29 | 1.31 | 194.7-195.2 K |
| 186.6 hPa | 4.9% | 15.7% | 57.7 | 2.21 | 2.24 | 185.9-186.2 K |
| 154.3 hPa | 8.7% | 100% | 303 | 9.06 | 9.19 | 176.2-176.4 K |

**S1: every layer is supersaturated in fact.** Each layer ends more than
1 percent of its steps above 1.01, and the layer at 288 hPa ends 94 percent of
them there.

**S2: averaging is excluded.** The mean of r is above 1.01 everywhere, and
equal to the ratio of the means within 3 percent. Over these steps the
temperature at each layer varies by only about 1 K, so there is no cycle for
the convexity to act on.

**Verdict: ordering.** Steps end with vapour above saturation, up to 300
times saturation. By the frozen table, the composition that leaves it there is
found under its own protocol, before the cloud field.

Two more observations:

- Layers at 254 and 154 hPa never end a step below 0.99. They are held at or
  above saturation every step, with excursions far above it.
- So the humidity above 320 hPa in every tropical result so far, including the
  H term of `tropical-surface-diagnostic.md`, carries vapour that the
  thermodynamics forbids.
