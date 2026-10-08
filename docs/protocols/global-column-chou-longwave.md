# Protocol: the global column with Chou's water vapour and the observed clouds

**Frozen before the wiring and before any run.** Results are appended below a
line at the end; nothing above it changes afterwards.

## Why

`water-vapour-longwave.md` attempt 1 replaced the fitted vapour envelope with
Chou et al. (2001)'s published k-distribution.

- **Where it failed.** V1 failed on the mid-latitude summer atmosphere: OLR
  287.80 against line-by-line 293.99, a miss of 6.19 W m^-2 against a limit of
  3.
- **Where it passed.** V2-V5 passed on both atmospheres, and every band was
  within 3 W m^-2 up and 4 down. The sub-arctic winter passed everything.

**What the integration diagnostic found.** `transmittance-integration-diagnostic.md`
attributes at most 1.6 W m^-2 of the shortfall to this model's solver. Up to
about 2.3 is shared by every implementation tried on this benchmark
atmosphere, including Chou's own method. Under its frozen rule nothing was
concluded.

**The user's decision (2026-10-07).** Accept attempt 1 despite V1.

- **What it does not change.** The frozen verdict above stands as written.
  This protocol does not re-judge it.
- **What it overrides.** The "if any fails, the option is not used" clause
  of that protocol, by the user's decision.
- **Why that is legitimate.** The 3 W m^-2 limit was set by Claude, not
  derived from the data.

**What this protocol does.** It takes the step that protocol said would
follow a pass: the global reference column with this vapour and the global
cloud field. With it, nothing in the longwave is fitted:

- the vapour is Chou's Tables 3, 4, 9 and 10;
- the CO2 band is Chou's Tables 3 and 10;
- ozone's 9.6-micron band is fitted to Chou's tables, not to any flux;
- the clouds carry their own longwave effect.

The surface temperature is then a prediction with no longwave calibration.

## The wiring

**A configuration option, `chou_water_vapour`,** off by default.

- **What it requires:** `co2_band` (band 3 is shared) and the radiative
  vapour.
- **Ozone's 9.6-micron band, if on,** is built into the same grid as band 5's
  ozone sub-bands split by the water terms, rather than by
  `ozone_infrared_grid`.

**The grid.** `chou_longwave_grid` (attempt 1, `tabulated=()`) rebuilds
everything at 4 microns and longer.

- **CO2.** The CO2 band's coefficients and kinds come from that grid, in
  place of `co2_band_grid`.
- **Below 4 microns** the grid is kept and carried as before: the solar vapour
  k-distribution, ozone's solar bands, and the albedo.

**The radiation module** takes the `ChouLongwave` spec. When it is present and
water is coupled:

- **The vapour optical depth** is `water_vapour_optical_depth`, from the
  layer vapour, mass, pressure and temperature, every step. This covers lines
  and continuum, and it replaces `vapour_absorption_optical_depth`.
- **The solar vapour** is added as before. The two are disjoint by band.
- **Pressure is declared as an input.**
- **A new equation record,** `EQ-H2O-LONGWAVE-K-DISTRIBUTION`, is registered
  before the module cites it.
- **The spec's source enters `parameters()`,** so the run hashes differently.

**A declared approximation, new here.** The bin from the grid's 3.16-micron
edge to 4 microns (2500-3160 cm^-1) carries no water absorption.

- **Before.** It lay inside a 3.16-5.62-micron bin that carried the fitted
  envelope.
- **Now.** Chou's band 9 is cut at 2500 cm^-1, as declared in
  `water-vapour-longwave.md`.
- **Its effect is reported:** the sunlight in that bin at the top and at the
  surface.

**Checks before any run,** as unit tests:

| check | passes if |
|---|---|
| C1: identity | on a fixed state, the module's clear-sky upward flux at the top equals a direct solve with `water_vapour_optical_depth` + `co2_optical_depth` + `ozone_infrared_optical_depth` on the same grid, to 1e-9 relative |
| C2: dry limit | with zero vapour, the vapour contributes no optical depth in any bin |
| C3: refusals | `chou_water_vapour` without the CO2 band, or without radiative vapour, raises `ConfigurationError` |
| C4: no envelope | with the option on, `vapour_absorption_optical_depth` is not called |

The full gate (ruff, mypy, pytest, traceability) passes before the runs.

## The runs

**The global reference column,** as configured for
`tropical-cloud-field.md`'s refit and the earlier 51-layer global runs:

- `reference_column`, `column_grid(16)`;
- 51 layers to 1 hPa (`max_log_step=0.2`, `top_pressure=100` Pa), started
  from the standard atmosphere;
- `master_dt` 3600 s, areal heat capacity 1e6 J m^-2 K^-1, initial water
  store 1e7;
- measured solar spectrum, solar ozone, CO2 band at 284.7 ppmv, eddy diffusion
  0.39;
- `condense_after_convection` and `end_of_step_saturation`;
- water vapour's solar absorption;
- the global ISCCP cloud field (`ISCCP_CLOUD_TYPES`), with the 0.17 surface
  reflectance;
- 200000 steps, with the means taken over the last fifth.

**Three configurations:**

| run | longwave vapour | ozone 9.6-micron band |
|---|---|---|
| **R0, control** | the fitted envelope at 3411.86 m^-1, its refit with this cloud field (`tropical-cloud-field.md`) | off, as that refit was made |
| **R1, Chou** | `chou_water_vapour` | on, as Chou's benchmark had it |
| **R2, Chou without ozone's band** | `chou_water_vapour` | off |

- R1 against R0 is the full change.
- R2 separates the vapour scheme from ozone's band.

## The tests

**G0: validity, for each run.** The run must be complete and valid, its
budgets must close, and its top-of-atmosphere imbalance must be at most
0.1 W m^-2.

- **A run that fails G0** fails everything below.
- **An unequilibrated run** reports its imbalance, not an error.

**G1-G4 are judged on R1** against CERES EBAF Ed4.0, global, from Loeb et al.
(2018), *J. Climate* 31, Table 5. The tolerances are `observed-cloud-field.md`'s.

| test | observed | passes if |
|---|---:|---|
| G1: planetary albedo | 0.291 | within 0.03 |
| G2: shortwave cloud radiative effect | -45.8 W m^-2 | within 12 |
| G3: longwave cloud radiative effect | +27.9 W m^-2 | within 8 |
| G4: the surface | 288.0 K | \|T_s(R1) - 288.0\| <= \|T_s(R0) - 288.0\| |

**What G4 asks.** Whether replacing every fitted longwave number with
published physics makes the headline no worse. It compares two errors. It does
not tune either.

**Reported, not judged, for every run:**

- **The surface-temperature error.** For R1 it is labelled: *51 layers, both
  saturation options, observed global cloud field, no longwave calibration*.
- **The top of the atmosphere:**
  - outgoing longwave against 240.1, and its clear-sky value against 268.1;
  - reflected sunlight against 99.1, and its clear-sky value against 53.3.
- **Water:** precipitable water, and relative humidity at 300-700 hPa.
- **The clouds:** the model's temperature at each cloud top, against ISCCP's
  227.5, 234.8, 262.8 and 281.1 K (`observed-cloud-field.md`).
- **The upper atmosphere:**
  - the cold point;
  - the temperatures at 60, 27, 10 and 1.5 hPa;
  - the stratospheric water, 100-10 hPa.
- **The surface:** net radiation, evaporation and the Bowen ratio.
- **CO2 doubling, R1 and R2.** The forcing per doubling at the top and at
  200 hPa, instantaneous, on the equilibrated mean state with its clouds.
  It is the first with no fitted number in the longwave.
- **The 3.16-4.0-micron bin:** its sunlight at the top, and at the surface in
  R1.
- **The cost:** the number of bins and the wall-clock time of each run.

## What follows

**Adoption.** Nothing here changes `reference_column`'s defaults. Whether
`chou_water_vapour` becomes the reference column's longwave is the user's
decision, made on these results.

**Then:**

- the tropical column's attempt 3 under `tropical-column.md`, with this vapour
  if adopted, frozen before it runs;
- the open decisions on the saturation options and the ozone band, which
  this protocol uses but does not decide.

---

## The user's decision, recorded before R1 finished (nothing above this line changed)

Recorded on 2026-10-07. R0 had finished; R1 and R2 were still running, and no
R1 output had been read.

- **If R1 fails any of G0-G4**, Chou's vapour is not adopted. The next step
  is RRTMG's longwave gas optics (Mlawer et al. 1997; Iacono et al. 2008),
  under its own protocol, frozen before any code.
- **If R1 passes G0-G4,** adoption goes back to the user as frozen above.
  RRTMG stays a candidate refinement.
- **The user's reasoning.** Line-by-line-derived k-distributions are likely
  part of the path in any case, and the thermal infrared is a large lever on
  the surface temperature.

## Results, R0 and R1 (appended; nothing above this line changed)

**Both runs are complete and valid,** at 200000 steps:

| run | residual | imbalance, W m^-2 | domain violations |
|---|---:|---:|---|
| R0 | 9.8e-16 | 0.004 | none |
| R1 | 8.7e-16 | 0.006 | none |

**G0 passes for both.** R2 was stopped at about three hours to keep memory
for R0 and R1, and was restarted from scratch when R0 finished. Its result
is appended when it completes.

**The verdict on R1,** against CERES EBAF Ed4.0 (Loeb et al. 2018, Table 5):

| test | observed | R0 (control) | R1 (Chou) | verdict on R1 |
|---|---:|---:|---:|---|
| G1: planetary albedo | 0.291 | 0.284 | **0.277** | **pass** (0.014) |
| G2: shortwave cloud effect | -45.8 | -52.5 | **-50.9** | **pass** (5.1) |
| G3: longwave cloud effect | +27.9 | +40.0 | **+23.5** | **pass** (4.4) |
| G4: the surface, against 288.0 K | | 289.59 K (+1.59) | **297.91 K (+9.91)** | **fails**: \|9.91\| > \|1.59\| |

**Verdict as frozen: R1 fails G4.** Under the user's decision recorded
above, Chou's vapour is not adopted, and the next step is RRTMG's longwave
gas optics under its own protocol.

**Reported:**

| | observed | R0 | R1 |
|---|---:|---:|---:|
| outgoing longwave | 240.1 | 243.6 | 246.2 W m^-2 |
| clear-sky outgoing longwave | 268.1 | 283.6 | **269.6** W m^-2 |
| reflected sunlight | 99.1 | 96.7 | 94.1 W m^-2 |
| clear-sky reflected | 53.3 | 44.2 | 43.2 W m^-2 |
| sunlight absorbed in the atmosphere | | 58.2 | 69.5 W m^-2 |
| precipitable water | ~25 | 32.4 | **70.6** kg m^-2 |
| relative humidity, 300-700 hPa | 0.4-0.5 | 0.79 | **0.89** |
| latent heat flux | | 82.9 | 114.1 W m^-2 |
| sensible heat flux | | 3.8 | 1.7 W m^-2 |
| cold point | | 187.8 K at 85 hPa | 202.8 K at 69 hPa |
| T at 60, 27, 10, 1.5 hPa | | 189.1, 201.5, 225.2, 282.8 K | 203.1, 209.2, 222.8, 269.4 K |
| stratospheric water, 100-10 hPa | | 2.3 ppmv | **26.8 ppmv** |

**The cloud tops,** against ISCCP's 227.5, 234.8, 262.8 and 281.1 K:

- R0: 222.3, 232.7, 263.7 and 281.3 K;
- R1: 238.3, 251.6, 276.3 and 290.9 K.

**CO2 doubling, R1.** Instantaneous, thermal, on the tail-mean state with its
clouds:

- 2.01 W m^-2 at the top;
- 3.87 W m^-2 at 203 hPa.

**The 3.16-4.0-micron bin, R1.** It carries 2.89 W m^-2 of sunlight at the
top and 2.46 at the surface, unabsorbed by water.

**The cost.** R0 has 67 bins and took 11505 s. R1 has 201 bins and took
26542 s; it ran alongside the other runs throughout.

**What the failure is, measured.** At the top of the atmosphere R1 is close
to observed throughout:

- the clear-sky outgoing longwave is 1.5 W m^-2 from CERES, where the fitted
  control is 15.5 off;
- both cloud effects pass.

**Where it goes wrong is the atmosphere's water.** It holds 2.8 times the
observed precipitable water and 0.89 relative humidity in the free
troposphere. The stratosphere holds 26.8 ppmv.

- **The clear-sky greenhouse effect,** σT_s^4 minus the clear-sky outgoing
  longwave, is about 175 W m^-2 against an observed 130.
- **So the surface is warm under a correct top-of-atmosphere flux** because
  the column is too opaque for its surface temperature. That is the
  signature of too much vapour.
- **The mechanism suspected but not measured,** carried from
  `observed-cloud-field.md`: subsiding air that dries the free troposphere.
- **Not measured here:** how much of R1's warmth comes from the scheme's
  extra opacity (V1's -6.2 W m^-2 on mid-latitude summer), amplified by the
  vapour feedback, and how much from the humidity closure alone.

## Results, R2 (appended; nothing above this line changed)

**R2 is complete and valid.** It has 176 bins and took 28416 s. The residual
is 8.0e-16, the imbalance 0.003 W m^-2, and there are no domain violations.
G0 passes.

| | observed | R0 | R2 (Chou, no ozone band) | R1 (Chou + ozone band) |
|---|---:|---:|---:|---:|
| surface error | 0 | +1.59 K | **+7.18 K** | +9.91 K |
| planetary albedo | 0.291 | 0.284 | 0.281 | 0.277 |
| shortwave cloud effect | -45.8 | -52.5 | -51.6 | -50.9 W m^-2 |
| longwave cloud effect | +27.9 | +40.0 | +29.3 | +23.5 W m^-2 |
| outgoing longwave | 240.1 | 243.6 | 244.8 | 246.2 W m^-2 |
| clear-sky outgoing longwave | 268.1 | 283.6 | 274.0 | 269.6 W m^-2 |
| precipitable water | ~25 | 32.4 | 54.8 | 70.6 kg m^-2 |
| relative humidity, 300-700 hPa | 0.4-0.5 | 0.79 | 0.86 | 0.89 |
| cold point | | 187.8 K, 85 hPa | 194.7 K, 69 hPa | 202.8 K, 69 hPa |
| stratospheric water | | 2.3 | 8.5 | 26.8 ppmv |
| CO2 doubling, top / 203 hPa | | | 1.60 / 3.78 | 2.01 / 3.87 W m^-2 |

R2's cloud-top temperatures are 231.4, 245.1, 272.4 and 287.7 K.

**The separation the protocol asked for:**

- **The vapour scheme (R2 - R0): +5.6 K.**
- **Ozone's 9.6-micron band (R1 - R2): +2.7 K.** Ozone's band warms the
  tropopause by 8 K and triples the water entering the stratosphere. The
  surface follows the moister column.

**Each step adds vapour.** In every run the surface warms with the column's
water, 32 to 55 to 71 kg m^-2. This is consistent with R1's finding that the
humidity, not the top-of-atmosphere flux, is what is wrong.
