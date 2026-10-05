# Protocol: ozone's 9.6-micron band

**Frozen before the tables are extracted, the fit is made, or any run.** Results
are appended below a line at the end; nothing above it changes afterwards.

## Why

Two diagnostics point here.

- **`tropical-heating-diagnostic.md`.** It inferred most of the missing
  lower-stratospheric heating above 70 hPa as a remainder, 0.9 to 1.0 K/day at
  31 to 46 hPa. It named ozone's absent 9.6-micron band as the likely cause,
  without measuring it.
- **`tropical-cloud-field.md`.** It fixed the tropical troposphere and left the
  tropical column 12 to 20 K too cold above 70 hPa, where the troposphere
  supplies little.

Gettelman et al. (2004, JGR 109, D07103) put ozone's longwave absorption of
upwelling infrared as the main radiative heating that the CO2 and vapour
cooling partly balance around the tropical tropopause.

## The source

Chou, Suarez, Liang and Yan (2001), the source of the CO2 band, treat ozone's
9.6-micron band (their band 5, 980-1100 cm^-1) by table look-up, section 4.2.

- **Their tables.** Band-mean transmittance pre-computed by line-by-line from
  HITRAN 1996. The weak ozone absorption of bands 4 and 6 is folded in (their
  eq. 4.17-4.18). A quadratic temperature fit,
  T(u, p, theta) = a + b (theta - 250) + c (theta - 250)^2, makes three tables,
  `o1`, `o2` and `o3`.
- **Their grid** (Table 8):
  - log10 u from -6 in 21 steps of 0.3, with u in cm-atm;
  - log10 p from -2 in 26 steps of 0.2, with p in hPa.
- **Where the tables are.** They are not printed in the report. They are in
  the public-domain code that implements it: WRF's `phys/module_ra_goddard.F`
  (Goddard longwave, Chou and Suarez). The file is fetched from the WRF
  repository and recorded by URL, commit and SHA-256 in `data/sources.yaml`.
- **The extraction is checked.** It must give 21 by 26 by 3 values. a must
  tend to 1 as u goes to 0, and a must not increase with u at any pressure.

## The band, derived from the tables

The model's solver takes a per-layer optical depth in each spectral bin, as
the CO2 band does. A path-transmittance table cannot be entered directly, so a
k-distribution is fitted to it: an exponential sum, the method behind Chou's
own CO2 and vapour k-distributions.

**The form.** Six terms, k_n = k_1 r^(n-1):

    T(u, p, theta) = sum_n g_n exp(-1.66 k_n u (p / 30 hPa)^m h(theta)),
    h(theta) = 1 + alpha (theta - 250) + beta (theta - 250)^2

**The free parameters:**

- k_1 and r;
- the g_n, which are non-negative and sum to 1;
- m, alpha and beta.

That is ten parameters for 1638 table values.

**The fit region.**

- p from 0.1 to 1000 hPa;
- u over the whole table;
- theta at 190, 220, 250, 280 and 310 K, from a, b and c.

The fit is least squares in transmittance, with every table point weighted
equally.

**How it enters the grid,** as the CO2 band does:

- 980-1100 cm^-1 is split into six sub-bands, each sized to its term's g_n by
  the Planck flux at 250 K;
- the vapour's own absorption is added in each;
- edges outside the band are kept.

The amount is each layer's ozone in cm-atm, from the column's existing SWOOSH
profile.

**Declared approximations.**

- The tables are band means. Splitting the band by g at 250 K is the
  k-distribution assumption.
- Correlation across layers is the k-distribution's own.
- Fitting is done on homogeneous paths, at the effective pressure and
  temperature of Chou's scaling.

## The tests

- **Z0, the fit,** at theta = 250 K over the fit region: RMS transmittance
  error <= 0.01 and maximum <= 0.03. Across the five temperatures: RMS <=
  0.02. If Z0 fails, nothing else is run and the form is not changed after the
  fact.
- **Z1, like for like with line-by-line,** on the AFGL mid-latitude summer
  atmosphere as `co2-band.md`'s C1 used it: its ozone, CO2 at 350 ppmv, and
  Chou's 4e-6 g/g of water above the tropopause.
  - **The test:** the top-of-atmosphere upward flux in 980-1100 cm^-1 must be
    within 3 W m^-2 of the line-by-line 58.62 (Chou et al. 2001, Table 14).
  - **Reported:** the same without ozone, so the band's own effect is visible.
- **Z2, the stratopause,** on the same atmosphere: the maximum longwave
  cooling above 100 hPa, CO2 and ozone together. It must lie within -7 to
  -14 C/day, the window C2 used; their all-gas line-by-line value is about
  -12. Reported: the change from the CO2-only -10.67.
- **Z3, the sign at the tropical tropopause.** On the observed +/-25 profile
  of `tropical-heating-diagnostic.md`, the band's longwave heating at the
  layers centred between 100 and 70 hPa must be positive (Gettelman et al.
  2004).
  - Reported: its value at each diagnostic level, against that diagnostic's
    remainder.

**The calibration it forces.** The vapour wing decay is refitted with this band
in its calibration, by the unchanged method:

- once without clouds, for comparison with `co2-band.md`;
- once with the global cloud field, for the tropical column.

Reported: the forcing per CO2 doubling, both ways.

**Next, if Z0 to Z3 pass:** `tropical-column.md` attempt 3, frozen before it
runs, with:

- the cloud configuration of `tropical-cloud-field.md`;
- this band;
- the refit that includes both.

If any test fails, the finding is recorded and the band is not used.

---

## Results (appended; nothing above this line changed)

**The tables are extracted and checked.**

- **Source.** WRF `phys/module_ra_goddard.F` at commit
  `7042598715a710930163485deae7742ab586af05`, SHA-256
  `0917531726cbed5a128e05c45e70ee292c209dbcd1ae2bad17db6e46b8d3a5ad`.
- **Checks.** There are 26 by 21 by 3 values. a is 0.999994 to 0.999997 at the
  smallest amount and never increases with u (the largest step is -6e-6).
- **Scale.** At 1 cm-atm, a is 0.84 at 0.01 hPa, 0.41 at 100 hPa and 0.33 at
  1000 hPa.

**Z0 fails.** The best of 60 least-squares starts gives:

| | limit | fit |
|---|---:|---:|
| RMS at 250 K | 0.010 | **0.0119** |
| maximum at 250 K | 0.03 | **0.047** |
| RMS across 190-310 K | 0.02 | 0.013 |

The fitted values are:

- k_1 = 0.0105 (cm-atm)^-1 and r = 10.16;
- g = 0.181, 0.311, 0.402, 0.081, 0.020 and 0.0045;
- m = 0.419;
- alpha = 0.00463 and beta = 4.2e-6.

**Where it fails.** By pressure, the RMS error at 250 K is:

| pressure | RMS | maximum |
|---|---:|---:|
| 0.1 hPa | 0.012 | |
| 1.6 to 2.5 hPa | 0.011 | |
| 10 to 16 hPa | 0.006-0.007 (lowest) | |
| 63 to 158 hPa | 0.014-0.015 | 0.039 |
| 1000 hPa | 0.024 | 0.048 |

A single power law in pressure cannot follow the band from Doppler-broadened
lines at the top to pressure-broadened ones below. Chou et al. (2001, section
4.2) give exactly this as their reason for treating the band by tables. The
error is largest where this band is wanted: around the tropical tropopause.

As frozen, nothing else is run, the form is not changed after the fact, and
the band is not used.

## Attempt 2, frozen before its fit

**The change, and why it is not a looser test.** Attempt 1's single power law
in pressure is the scaling Chou et al. reject for this band. Attempt 2 tabulates
the absorption coefficients in pressure, which is the correlated-k form used by
table-driven schemes. Everything else is unchanged: the source, the fit region,
the temperatures and every test with its threshold.

**The form.**

- **The weights.** Six terms with weights g_n, non-negative and summing to 1,
  shared by every pressure, because they set the spectral sub-bands.
- **The coefficients.** At each of the 21 table pressures from 0.1 to
  1000 hPa, six coefficients k_n(p) are fitted independently.
- **The temperature.** Scaling h(theta) = 1 + alpha (theta - 250) + beta
  (theta - 250)^2, shared.
- **Per layer,** log k_n is interpolated linearly in log p between table
  pressures.
- **Parameters:** 5 + 126 + 2 = 133, for 2205 table values.

**The tests.**

- **Z0**, unchanged and applied to this form: RMS <= 0.01 and maximum <= 0.03
  at 250 K, and RMS <= 0.02 across the five temperatures.
- **Z0b, new, the interpolation.** The same form is fitted to the even-indexed
  table pressures only. Then the k_n are interpolated to the odd-indexed ones,
  which were held out. Their RMS error at 250 K must be <= 0.01 and the
  maximum <= 0.03. This tests what a model layer between table pressures
  relies on. The band used in the model is the fit to all 21 pressures.

Z1, Z2 and Z3, the refits and what follows are exactly as frozen above. If Z0
or Z0b fails, nothing else is run.

## Attempt 2, results (appended; nothing above this line changed)

**Z0 and Z0b pass.** One least-squares fit from attempt 1's solution, over the
same 2205 table values:

| test | limit | result |
|---|---:|---:|
| Z0: RMS at 250 K | 0.010 | **0.0015** |
| Z0: maximum at 250 K | 0.03 | **0.0103** |
| Z0: RMS across 190-310 K | 0.02 | **0.0058** |
| Z0b: held-out pressures, RMS at 250 K | 0.010 | **0.0018** |
| Z0b: held-out pressures, maximum at 250 K | 0.03 | **0.0091** |

**The fitted band.**

- **Error by pressure.** The RMS error is 0.0005 at 10-16 hPa and 0.0028 at
  1000 hPa.
- **Weights:** 4e-6, 0.549, 0.202, 0.172, 0.058 and 0.019. The fit found five
  terms enough, and the sixth stays rather than the form changing.
- **Temperature scaling:** alpha 0.00462 and beta 4.1e-6.
- **Thin air, measured in the fit, not assumed.** Going from 500 to 0.1 hPa,
  the four terms holding 98% of the weight absorb 4 to 8 times less per
  molecule, as the wings narrow. The strongest term absorbs 3.6 times more, as
  the line peaks sharpen.
- **Built in** `physics/ozone_infrared.py`, off by default as
  `ozone_infrared`.

**Z1 fails, and not because of ozone.** On the MLS atmosphere, at the vapour
wing C1 used (4789.08 m^-1):

| | model | line-by-line |
|---|---:|---:|
| TOA upward flux, 980-1100 cm^-1, with ozone | **21.14** W m^-2 | 58.62 |
| the same, without ozone | 27.88 W m^-2 | |
| C1, 540-800 cm^-1 (unchanged, a check) | 68.31 | 67.90 |
| total OLR, with ozone | 264.4 W m^-2 | 293.99 |

- **The band's own effect** is -6.7 W m^-2 at the top in its band.
- **The failure is not the band's.** Without any ozone the model's window
  already lets out 27.9 W m^-2 where line-by-line lets out 58.6 *with* ozone.
  The model's water vapour makes the 980-1100 cm^-1 window far too opaque.
- **It is not the wing decay either.** With the cloud-calibrated wing
  (3411.86), the band flux is unchanged at 21.2 with ozone and 28.0 without,
  while the total OLR rises to 295.1. So what closes the window is not the
  wing, and the vapour's window absorption has not yet been examined.

**Z2 passes.** The maximum stratospheric cooling with CO2 and ozone is
-13.22 C/day at 1.1 hPa, against about -12 line-by-line, inside -7 to -14.
With CO2 alone it was -10.69.

**Z3 passes.** On the observed +/-25 profile, the band heats the 84.7 hPa
layer by +0.237 K/day. Against `tropical-heating-diagnostic.md`'s remainder:

| layer | band heating | remainder |
|---:|---:|---:|
| 154 hPa | +0.097 | |
| 103 hPa | +0.178 | 0.232 |
| 85 hPa | +0.237 | 0.118 |
| 69 hPa | +0.306 | 0.283 |
| 46 hPa | +0.374 | 0.911 |
| 31 hPa | +0.341 | 1.018 |

- **Near the tropopause** it accounts for about all of the remainder.
- **At 31-46 hPa** it accounts for 35 to 40 percent.
- **Probably an underestimate.** It is computed through the same over-opaque
  window, which starves ozone of the upwelling infrared it absorbs.

**As frozen, Z1 failed, so the band is not used.** The finding is that the
model's vapour closes the infrared window, and that is now the first thing to
measure. The band waits on that measurement.

## Erratum to Z1 (appended after the results; nothing above changed)

**The reference value frozen for Z1 was misread from the source.** Z1 names the
TOA upward flux in 980-1100 cm^-1 in Chou et al. (2001) Table 14. The frozen
text gives that as 58.62 W m^-2, but 58.62 is the row for **800-980 cm^-1**.

- **How the table is laid out.** Its rows run 0-340, 340-540, 540-800,
  800-980, 980-1100, 1100-1215, 1215-1380, 1380-1900 and 1900-3000 cm^-1.
- **The check on the reading.** Read that way, their upward fluxes at the top
  sum to the table's total, 293.99. The 540-800 row is C1's 67.90, as
  `co2-band.md` used it.
- **The row Z1 names.** For 980-1100 cm^-1 the line-by-line value is
  **22.36** W m^-2, with Chou's parameterization at 21.41.

The threshold and the band are unchanged; only the quoted value was wrong.
Because this is found after the result, both readings are recorded:

| Z1 | reference | model, with ozone | difference | verdict |
|---|---:|---:|---:|---|
| as written (the 800-980 row's value) | 58.62 | 21.14 | 37.5 | fail |
| as named: Table 14, 980-1100 cm^-1 | **22.36** | 21.14 | **1.22** | **pass** (within 3) |

**Two earlier statements in these results are wrong and are withdrawn:**

- that the model's vapour makes the 980-1100 cm^-1 window "far too opaque";
- that "the vapour's window absorption has not yet been examined" explains a
  failure.

Without ozone the model's band flux is 27.88 W m^-2. That is 5.5 above the
line-by-line value *with* ozone, consistent with a band effect of the size
measured, 6.7.

**On the corrected reading, Z0, Z0b, Z1, Z2 and Z3 all pass.**

- **What that would mean.** By the frozen rule, a full pass makes the band
  eligible for the refits and for `tropical-column.md` attempt 3.
- **Why it is not adopted here.** Whether it is, given that the pass rests on
  a correction made after the result, is a decision for the person directing
  this work, and is not taken here.

**What remains open.** The MLS total OLR is 264.4 W m^-2 against 293.99. The
deficit is in other bands, and finding which is the next measurement.
