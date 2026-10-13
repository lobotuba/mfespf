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
