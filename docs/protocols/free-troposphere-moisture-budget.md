# Protocol: what holds the free troposphere at 0.9?

**Frozen before the runs.** Results are appended below a line at the end;
nothing above it changes afterwards. This is a diagnostic: it changes no
module and no default.

## Why

`tropical-runaway-ceiling.md` showed the column's radiation is right and its
humidity is not.

- **At the observed tropical humidity,** the tropical column equilibrates
  within 0.9 K of the observed surface.
- **At the 0.95 the column makes,** it has no equilibrium.
- **The threshold:** the free troposphere must stay below about 0.75.

R4, the global column, holds 0.89 over 300-700 hPa against an observed
0.4-0.5.

**The user chose to find the process that sets it before building anything.**

**What the code says, read before freezing.** The prognostic water module
changes the vapour in three labelled ways:

| label | module | what it is |
|---|---|---|
| `convective_vapour_transport` | convection | mixes specific humidity uniformly over each convectively mixed group of layers (`_vapour_transport`) |
| `vapour_eddy_diffusion` | water | down-gradient diffusion at 0.39 m^2 s^-1 |
| `evaporation_minus_condensation` | water | surface evaporation into the lowest layer, plus the saturation adjustment, which removes everything above saturation as precipitation |

**There is no re-evaporation of falling condensate, and no descending
motion.** Convection's mixing has no compensating subsidence. On reading,
condensation is the only process that can dry a free-tropospheric layer.

**The hypothesis H.** The free troposphere is near saturation because vapour
is carried up and removed only by condensation; nothing in the column
descends and dries. In real convecting columns, the environment's
compensating subsidence does the drying (Arakawa & Schubert 1974, *J. Atmos.
Sci.* 31, 674; Romps 2014, *J. Climate* 27, 7432).

## The runs

**The state.** R4's configuration (`longwave-source-function.md`) restarted
from its saved tail-mean state, `gcl_r4_200000_state.npz`:

- temperatures, vapour and surface temperature, through `initial_state`;
- the vapour now included, which the earlier restart lacked.

**The runs, all from that state:**

| run | change | steps |
|---|---|---:|
| B0 | none | 10000 |
| Da | `eddy_diffusivity = 0` | 5000 |
| Db | `vapour_transport = False`, which requires `condense_after_convection = False` | 5000 |

**The check S:** B0's means over its last 5000 steps stay R4's:

- surface within 0.5 K of 301.49;
- RH 300-700 within 0.03 of 0.891;
- precipitable water within 5 kg m^-2 of 98.4.

**If S fails,** the budget is still reported, labelled as the restart's own
state.

## The measurements

**M1, closure.** In every step of B0, each layer's vapour change equals the
sum of its labelled contributions times dt, to 1e-12 kg m^-2. This is the
STA-006 record, checked.

**M2, the budget, reported.** For the layers centred between 300 and
700 hPa, from B0's last 5000 steps:

- the time-mean of each label, in kg m^-2 day^-1;
- each label's share of the total moistening;
- the share of the total drying.

**MH, the hypothesis, judged.** In each layer centred between 300 and
700 hPa:

- **H holds** if the only label with a negative time-mean is
  `evaporation_minus_condensation`, and that label is never positive there.
  Then it is pure condensation, since no evaporation reaches those layers.
- **H fails** if any transport label dries any of those layers on average.

**M3, denial, reported and judged.** For Da and Db, the mean RH over
300-700 hPa across the last 2500 steps, against B0's.

- **A single process sets the humidity** if one denial brings it below 0.75.
- **The source is shared, or the result is condensation-limited** if neither
  does.

**Reported, not judged:** each run's surface temperature and precipitable
water over its last 2500 steps. These are short runs from R4's state, not
equilibria, and are read only as the direction and rate of change.

## What follows

**If H holds,** the missing physics is a drying process inside a convecting
column. The established candidate is the compensating subsidence of a
mass-flux convection scheme, with detrainment. It predicts a free-tropospheric
RH from its own entrainment and detrainment, as Romps' bulk plume does. Its
protocol would come next, as the user decides.

**If H fails,** the drying term that exists is examined first.

---

## Erratum, before Db's result (nothing above this line changed)

**Db as frozen cannot be built.**

- **What was frozen:** `vapour_transport = False` with
  `condense_after_convection = False`.
- **What the configuration does:** it refuses, because `end_of_step_saturation`
  needs the water module to run after convection, which
  `condense_after_convection` arranges.
- **What the protocol missed:** that dependency, which it did not check.

**The two readings:**

1. **Drop Db.** That loses the denial of the process M2 already shows
   supplies 96 percent of the moistening.
2. **Db also sets `end_of_step_saturation = False`,** taken.
   - That is the only configuration in which convection does not carry
     vapour.
   - The saturation adjustment still runs, at the start-of-step temperature.
   - The change of two switches is declared, and Db's response is read with
     it.

B0 and Da are unaffected and already complete.

## Result (appended; nothing above this line changed)

**S passes.** B0's last 5000 steps are R4's state: surface 301.49 K, RH
300-700 0.893, precipitable water 98.5 kg m^-2. Restoring the vapour did what
the earlier restart, which lacked it, could not.

**M1 passes.** The largest per-layer residual over all 10000 steps is
6.7e-16 kg m^-2.

**M2, the budget over the twelve layers centred between 692.7 and
321.5 hPa,** time-means in kg m^-2 day^-1, summed over those layers:

| label | sum | share of moistening | share of drying |
|---|---:|---:|---:|
| `convective_vapour_transport` | **+2.497** | **96.1%** | 0 |
| `vapour_eddy_diffusion` | +0.058 | 3.9% | 1.6% |
| `evaporation_minus_condensation` | **-2.546** | 0 | **98.4%** |

Convective mixing moistens every one of these layers, and condensation dries
every one. Diffusion is a net +0.058 that moves water upward between them:
it is negative in four layers, by at most 0.027.

**MH fails as frozen,** on two counts, both measured:

1. **Eddy diffusion is negative on average in four of the twelve layers.**
   That is redistribution: water passing upward from one layer to the next,
   1.6 percent of the drying.
2. **`evaporation_minus_condensation` is positive in 4344 of 12000
   layer-steps of a 1000-step probe.** Every one is at most 4.4e-16 kg m^-2
   per step, median 1.0e-16. That is the round-off of forming the label as a
   difference (the vapour rate less diffusion), not evaporation aloft.

The frozen criterion was per layer and exact, so both count. As written, H
fails.

**What H claimed, read on the integrated budget, does hold.** The free
troposphere gains water only by transport from below: 96 percent convective
mixing. It loses it only by condensation: 98 percent, with the remaining 1.6
percent diffusion moving it up. No process in the column descends and dries.

**M3, denial,** RH 300-700 over the last 2500 steps:

| run | change | RH 300-700 | precipitable water | surface |
|---|---|---:|---:|---:|
| B0 | none | 0.893 | 98.5 | 301.49 K |
| Da | eddy diffusion off | 0.886 | 98.0 | 301.53 K |
| Db | no convective vapour transport, and (erratum) no end-of-step saturation | **0.985** | **118.6** | 303.63 K |

**Verdict as frozen: neither denial brings the free troposphere below 0.75,
so no single process sets it.**

**Db is the informative one.** With the source supplying 96 percent of the
moistening switched off, the free troposphere got *wetter*.

- Nothing in the column removes the water already aloft except
  condensation, so the water stays at saturation.
- The upper free troposphere reaches 1.0005: the supersaturation that
  `end-of-step-saturation.md` removed returns when that option is off.
- The humidity is therefore not set by how much water convection supplies.
  It is set by the absence of anything that dries.

**What follows, as frozen ("If H holds"):** the integrated reading supports
H, though the per-layer criterion fails. The missing process is drying
inside a convecting column. The established form is the compensating
subsidence of mass-flux convection with detrainment (Arakawa & Schubert
1974). Its own entrainment and detrainment would predict the free
troposphere's RH (Romps 2014). That protocol is the user's decision.
