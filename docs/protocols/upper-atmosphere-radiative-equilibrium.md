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

---

## Result (appended; nothing above this line changed)

### Errata, both readings kept

**1. The ozone record.**

- **Reading 1 (above):** EQ-OZONE-SOLAR-ABSORPTION and the code "disagree".
- **Reading 2, on re-reading:** the record is incomplete, not contradicted.
  - Its "the column between the model top and 1 hPa placed in the top layer"
    is what the code does, for a top below 1 hPa.
  - What it omitted is that ozone above 1 hPa was left out.
- **The record now says so.** It also says that SWOOSH's 1 hPa value is never
  held upward.

**2. The data's provenance class.** The AFGL table was first registered as
"reference_profile", which is not a class the registry has. It is
`idealised`: a standard atmosphere built from observations, input rather than
evidence about the Earth.

### The data check

`data/raw/afgl_1986_us_standard_table_1f.csv` is joseki's transcription of
AFGL-TR-86-0110 Table 1f, sha256 d2d06777...7eda.

- **Checked against** the report's printed Table 1f ("Model 6, U.S.
  Standard"), from the archive.org OCR of DTIC_ADA175173.
- **O3 matches at seven levels:**

  | km | 42.5 | 45 | 47.5 | 50 | 55 | 65 | 70 |
  |---|---:|---:|---:|---:|---:|---:|---:|
  | ppmv | 6.20 | 5.25 | 4.10 | 3.10 | 1.80 | 0.70 | 0.30 |

- **The step where SWOOSH and AFGL meet at 1 hPa** is 3.18 against 3.82
  ppmv, AFGL 20 percent higher. It is reported, not smoothed.

### The marches

Every run converged by the frozen criterion, a final 30-day mean of
max |dT/dt| ≤ 0.01 K/day.

| run | top, layers | converged, days | dT_ref range | G1 | G2 | dry-unstable layers |
|---|---|---:|---|---|---|---|
| S1 | 1 hPa, 51 | 35.8 | -38.2 to +19.4 | pass | **fail** (top layer out of domain from day 0.5) | **the top, 1.07 hPa** |
| S2 | 1 hPa, finer, 78 | 36.0 | -38.7 to +19.4 | pass | **fail** | **the top** |
| S3 | 0.1 hPa, 62 | 48.8 | -27.7 to +19.4 | pass | pass | none |
| S4 | 0.1 hPa, finer, 101 | 49.3 | -27.6 to +19.4 | pass | pass | none |

**The top layer at equilibrium:**

- S1: 237.3 K, 15.4 K colder than the layer below, dT_ref -38.2 K.
- S2: 236.9 K, also cold and dry-unstable.
- **Finer resolution changes it by 0.5 K, so the defect is not resolution.**

**At 1.07 hPa:**

| | T | Q_SW | Q_LW | longwave in the CO2 band, 550-800 cm⁻¹ | O3 band | rest |
|---|---:|---:|---:|---:|---:|---:|
| S1, the top layer | 237.3 | 13.7 | -13.7 | -9.6 | -0.7 | -3.4 |
| S3, an interior layer | 257.4 | 9.0 | -9.0 | **-5.0** | -1.5 | -2.5 |

- **The CO2 band's cooling halves with the atmosphere above present,** even
  though the layer is 20 K warmer. That is the band centre no longer seeing
  space.
- **The sunlight falls** from 13.7 to 9.0 K/day, as the ozone above takes the
  ultraviolet first.

### Verdicts, as frozen

**H1, the open top, is supported.** Every condition holds:

1. S1's top layer fails G2, and ends 15.4 K colder than the layer below.
2. S3 at 1.07 hPa is 257.4 K against S1's 237.3, +20.1 K (≥ 5).
3. S3's 1.04 hPa layer cools at 1.02 times the layer below it (within
   25 percent), against S1's top-layer ratio at the start of 20.8/11.4 = 1.82.

**G6 passes.** S3 and S4 agree within 1.46 K at every pressure from 1 to
50 hPa; the largest difference is at 1 hPa.

**S5, the coupled check:** R0b's mass-flux configuration with interactive
radiation, at the raised top, for 240 h.

- **G3 passes.** Convection's temperature tendencies above 50 hPa integrate
  to at most 3e-13 K: the dry adjustment's round-off. The mass flux and the
  surface terms contribute exactly zero. The plume tops were at 103-154 hPa.
- **G5 passes.** Residual 5.1e-16, conservation gate passed.
- **The run is complete, not invalid.** The worst dT_ref is -26.5 K, at the
  103 hPa plume-top layer. It is the first mass-flux run in this branch to
  pass every gate.
- **The 1 hPa temperature holds:** 257.79 to 257.75 K over the 10 days.

**The raised top (0.1 hPa, AFGL ozone above 1 hPa) is a candidate upper
configuration.** S3 and S4 pass G1, G2, G4 and G6, and S5 passes G3 and G5.
The user decides adoption.

### What the result says, and what it does not

**The new top at 0.11 hPa keeps a weaker form of the boundary signature,**
answering the user's question of displacement against elimination.

| run | new top | the layer below | dT_ref |
|---|---:|---:|---:|
| S3 | 210.0 K | 7.0 K warmer | -23.9 |
| S4 | 207.9 K | 5.7 K warmer | -24.3 |

- It is stable and inside the domain. So the truncation is displaced to a
  pressure where it no longer breaks a gate, not eliminated.
- A top at 0.1 hPa is still a top.

**H2, the broad cold bias, remains.** Reported against the U.S. Standard
Atmosphere 1976, which is not a gate:

| hPa | S3 | US 1976 | difference |
|---:|---:|---:|---:|
| 1 | 257.8 | 270.7 | -12.9 |
| 3 | 241.7 | 249.5 | -7.8 |
| 10 | 217.9 | 227.7 | -9.8 |
| 30 | 205.9 | 220.5 | -14.6 |

- Raising the top does not change 3-30 hPa (S1 and S3 within 0.1-1.5 K).
- The candidates named above are untouched here:
  - O2 absorption, absent;
  - the diffuse treatment of sunlight;
  - the ozone amount;
  - the stratosphere's missing circulation and water.

**Not concluded:**

- no default changed;
- **the configuration under the raised top:**
  - the temperature of layers above 1 hPa is radiative only (no
    gravity-wave drag, no circulation);
  - their vapour is held at R4's top-layer mixing ratio, not measured.

## User disposition (2026-10-13), appended

> The raised top is not adopted. It is a useful diagnostic configuration that
> restores validity near 1 hPa but displaces a weaker upper-boundary
> temperature discontinuity to 0.1 hPa and does not resolve the broad
> stratospheric cold bias.

**The reasons:**

- The new terminal layer is about 7 K colder than the layer below. The
  protocol does not permit a sacrificial terminal layer.
- The 1-30 hPa cold bias of 8-15 K remains.

**What is kept:**

- **The ozone extension** (`ozone_above_1hpa`, AFGL Table 1f). It is
  source-traceable and option-controlled, and is required for any experiment
  above 1 hPa.
- **The default top is unchanged.** The raised vertical domain is not made the
  production default.

**Next stratospheric work, when resumed:** the upper-boundary energy treatment
and the broad radiative bias. Not another increase in the top's altitude.
