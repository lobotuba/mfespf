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
