# Protocol: the upper atmosphere's radiative equilibrium and the model top

**Frozen before any code or run of this protocol.** Results are appended below
a line at the end; nothing above it changes afterwards.

**Authorised by the user (2026-10-13).** It proceeds independently of F0.

**Its origin.** `mass-flux-diagnostic-branch.md` Part R found that R4's 1 hPa
layer was held by a convective adjustment heating it about 7 K/day, where no
convection exists. Removing it exposed a radiative equilibrium 36 K below
RRTMG's reference, outside its tables.

**Prohibited repair.** The upper-stratospheric convective adjustment is not
restored. That would be a compensating error, not a closure.

**The only admissible mitigations:**

- legitimate stratospheric physics, entering through the radiation;
- a declared upper boundary.

## What R4's upper atmosphere is, read before this was frozen

R4's saved state, one radiation solve. The heating is split after the solve
into its shortwave and longwave bands, and the split reproduces the published
heating to 2.5e-12 K/day.

| layer, hPa | interfaces, hPa | T, K | dT_ref, K | Q_SW, K/day | Q_LW, K/day | Q_net, K/day |
|---:|---|---:|---:|---:|---:|---:|
| 1.07 (top) | 1.14-1.00 | 247.8 | -27.7 | 13.71 | **-20.80** | **-7.09** |
| 1.27 | 1.40-1.14 | 249.9 | -24.7 | 13.17 | -11.40 | +1.77 |
| 1.55 | 1.71-1.40 | 252.4 | -19.7 | 11.98 | -10.14 | +1.84 |
| 1.90 | 2.08-1.71 | 254.0 | -14.5 | 10.31 | -10.31 | 0.00 |
| 2.31 | 2.55-2.08 | 249.9 | -14.6 | 8.37 | -8.37 | 0.00 |
| 2.83 | 3.11-2.55 | 244.5 | -16.1 | 6.56 | -6.56 | 0.00 |
| 3.45 | 3.80-3.11 | 239.2 | -17.6 | 5.18 | -5.18 | 0.00 |
| 4.22 | 4.64-3.80 | 234.5 | -18.6 | 4.23 | -4.23 | 0.00 |
| 5.15 | 5.66-4.64 | 230.2 | -19.3 | 3.53 | -3.53 | 0.00 |

**Two defects, of different kinds:**

1. **The top layer.** It absorbs as much sunlight as the layer below (13.7
   against 13.2 K/day) but cools twice as hard in the longwave (-20.8 against
   -11.4). The column ends at 1.00 hPa, with nothing above to emit back down
   in the strongly absorbing band centres. The real atmosphere above 1 hPa
   would.
2. **The broad cold bias.** Between 1.9 and 5 hPa the column is already in
   radiative balance, but 14-19 K below RRTMG's reference profile.

**The code, read:**

- **Ozone above 1 hPa is omitted.** `ozone_layer_path` says SWOOSH ends at
  1 hPa and "what lies above 1 hPa (about a dobson unit) is left out".
  EQ-OZONE-SOLAR-ABSORPTION says it is "placed in the top layer". The two
  disagree, and the code is what runs.
- **`_ozone_mixing_ratio` uses `np.interp`.** On a column whose top is above
  1 hPa it would hold SWOOSH's 1 hPa value upward: a held last value past a
  record's coverage, which the project forbids (DAT-007).
- **The shortwave absorbers are ozone and vapour only.** There is no O2:
  neither the Schumann-Runge bands nor the Herzberg continuum.
- **Sunlight is solved as diffuse radiation** with the solver's diffusivity
  factor, a declared approximation.

## Hypotheses

**H1, the open top:**

- **The claim:** the top layer's excess longwave cooling, and the cold,
  out-of-domain equilibrium it drives, come from ending the column at 1 hPa
  with nothing above it.
- **The prediction:**
  - the excess sits in the strongly absorbing longwave bands;
  - raising the model top above 1 hPa, with the atmosphere above it
    represented, removes it at 1 hPa.

**H2, the broad cold bias,** is diagnosed and reported here, not mitigated.
Its candidates, each for a later protocol:

- the missing O2 absorption;
- the diffuse treatment of sunlight;
- the ozone amount;
- the resolution.

## Method: the upper atmosphere in radiative equilibrium, with the troposphere held

**What is held and what moves.** From R4's saved state:

- held: every layer at or below 100 hPa, the surface and all vapour;
- marched: every layer above 100 hPa, by its radiative heating alone. There
  is no convection, no diffusion of heat and no circulation.

**The integration:**

- explicit, at 6 h steps;
- the radiation is solved on the current state every step, with the same
  module and options as R4;
- it runs until equilibrium, or 1000 days, whichever comes first.

**Equilibrium:** the largest |dT/dt| over the marched layers, averaged over
the final 30 days, is ≤ 0.01 K/day. A run that reaches 1000 days without it is
reported as not converged.

**The domain:** dT_ref (RRTMG's reference departure) in every layer, at every
step.

## Runs

**S1:** the current configuration. Top at 100 Pa, log step α = 0.2.

**S2:** S1 at α = 0.1, finer layers above the uniform-dp region. This is the
resolution at the top.

**S3:** the model top raised to 10 Pa (0.1 hPa), with α = 0.2.

- **The pressure range** is chosen to stay below where local thermodynamic
  equilibrium in the CO2 15-micron band fails (above about 0.05 hPa), and
  inside RRTMG's reference range, which reaches 0.01 hPa.
- **Ozone above 1 hPa** comes from the AFGL U.S. Standard profile (Anderson
  et al. 1986, AFGL-TR-86-0110, DTIC ADA175173).
  - It is taken from the machine-readable transcription in the `joseki`
    package (its `afgl_1986` tables).
  - Before any run, it is checked against the report's printed table at 3 or
    more levels, and recorded in `data/sources.yaml` with a checksum.
  - SWOOSH is kept from 316 to 1 hPa. AFGL is used only above 1 hPa,
    unscaled, and the step at 1 hPa is reported.
- **Initial state for the new layers:**
  - below 1 hPa, R4's temperature and vapour mixing ratio, interpolated in
    ln p;
  - above 1 hPa, the U.S. Standard Atmosphere 1976 temperature, and R4's
    top-layer vapour mixing ratio.
  - The equilibrium does not depend on the initial state; the initial state
    only has to start inside the domain.

**S4:** S3 at α = 0.1. This is the resolution at the raised top.

**S5, a coupled check,** for an upper configuration that passes G1, G2, G4
and G6:

- the mass-flux configuration of R0b, `mass-flux-diagnostic-branch.md`;
- 240 steps from R4's state below 100 hPa, and from that configuration's S
  equilibrium above;
- interactive radiation.

## Diagnostics

**Per layer, at the start and at equilibrium:**

- T;
- dT_ref;
- Q_SW, Q_LW and Q_net.

**In the layer containing 1.07 hPa,** Q_LW split into band groups:

- CO2 15 µm, 550-800 cm⁻¹;
- O3 9.6 µm, 980-1100 cm⁻¹;
- the rest.

**Also:**

- days to equilibrium;
- the trajectory's dT_ref range.

**Reported against the U.S. Standard Atmosphere 1976, not as a gate:** T at
1, 3, 10 and 30 hPa. It is a mid-latitude idealisation, not a global-mean
measurement.

## Gates, frozen (the user's list)

- **G1, the initial state is inside the radiation domain:** |dT_ref| ≤ 30 K in
  every layer at step 0.
- **G2, the evolving state is inside it:** |dT_ref| ≤ 30 K in every layer at
  every step.
- **G3, no convective tendency above its diagnosed top.** In S5, convection's
  temperature tendencies in layers above 50 hPa integrate to at most 0.01 K
  over the 10 days. The plume tops were at 103-154 hPa.
- **G4, a stable equilibrium:**
  - converged, by the criterion above;
  - and no layer is colder than the one below it by more than the dry
    adiabat allows, so that no dry adjustment could act.
- **G5, unchanged water and energy accounting:** S5 passes the conservation
  gate with residual ≤ 1e-12.
- **G6, the resolution at the top:** S3 and S4's equilibria, interpolated in
  ln p, agree within 2 K at every pressure from 1 to 50 hPa.

## Verdicts, frozen

**H1 supported if all of:**

- **in S1:** the top layer fails G2, or ends colder than the layer below by
  more than 5 K;
- **in S3:** the equilibrium temperature at 1.07 hPa (interpolated) exceeds
  S1's by 5 K or more;
- **in S3:** the 1.07 hPa layer's longwave cooling is within 25 percent of
  the next layer below's.

**H1 rejected** if S3's temperature at 1.07 hPa is within 2 K of S1's.

**The raised top is a candidate upper configuration** if S3 and S4 pass G1,
G2, G4 and G6, and S5 passes G3 and G5. The user decides adoption.

**Anything else is reported as found.**

## Code, introduced by this protocol only

- **`ozone_above_1hpa`,** "none" or "afgl_us_standard".
  - **Coverage:** `_ozone_mixing_ratio` refuses pressures above SWOOSH's
    coverage, that is below 1 hPa, unless the AFGL extension is declared. No
    held value.
  - **"none"** keeps today's behaviour exactly for a top at 1 hPa or below.
    A top above 1 hPa with "none" is refused.
  - **The record:** EQ-OZONE-SOLAR-ABSORPTION's ozone-profile entry, and its
    "placed in the top layer" approximation, are corrected to what the code
    does. The correction is an erratum, recorded with both readings.
- **The AFGL data:** in `data/raw`, with its `sources.yaml` entry and
  checksum, and the check against the report.
- **Tests:**
  - option-off identity;
  - the refusal above coverage;
  - the AFGL column above 1 hPa against the report's own integral;
  - continuity of the profile function.
- **The march is diagnostic code in the run script.** It drives the radiation
  module the way the common-radiation replay does, and is not a module
  option.
