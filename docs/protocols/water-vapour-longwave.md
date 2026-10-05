# Protocol: water vapour's longwave absorption from Chou et al. (2001)

**Frozen before the implementation and before any benchmark is computed.**
Results are appended below a line at the end; nothing above it changes
afterwards.

## Why

`longwave-band-comparison.md` found the model's mid-latitude summer outgoing
longwave 29.7 W m^-2 below line-by-line. Its vapour is a two-wing envelope with
one fitted decay, evaluated at the centre of coarse bins. That closes the
partly transparent regions at 1100-1380 and 340-540 cm^-1. With the
cloud-calibrated decay the total comes out right only by cancellation.

**The envelope has served two purposes.** It stands in for the vapour and,
through its fit to all-sky OLR, for the cloud longwave effect. Replacing it
with a published k-distribution removes the last fitted number in the
longwave. The clouds' effect then has to come from clouds.

## The scheme

Chou, Suarez, Liang and Yan (2001), NASA/TM-2001-104606 vol. 19, the source of
the CO2 band and of the ozone tables. Nothing is fitted or adjusted.

**Bands (Table 1), cm^-1:**

| band | 1 | 2 | 3a | 3b | 3c | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| range | 0-340 | 340-540 | 540-620 | 620-720 | 720-800 | 800-980 | 980-1100 | 1100-1215 | 1215-1380 | 1380-1900 | 1900-3000 |

Band 10 (N2O, 540-620) is minor and stays out, as before.

**Line absorption, by k-distribution, eq. 4.10.** Six terms per band,
k_n = k_1 eta^(n-1), with Planck-weighted fractions Delta-g:

- **Bands 1, 2 and 4-9:** Table 4 gives k_1/mu, eta and Delta-g.
- **Bands 3a, 3b and 3c:** Table 10 gives them. k_1/mu = 1.33e-2, eta = 8,
  and Delta-g for each sub-band, which sum to one over band 3.
- **Scaled amount,** eq. 4.4 with Table 3:
  w = q dp/g (p / 500 hPa)^m [1 + alpha (theta - 250) + beta (theta - 250)^2],
  with m, alpha and beta for each band.

**Continuum, eq. 4.19-4.24 and Table 9.**

- **Where:** bands 3a, 3b, 3c, 4, 5, 6 and 7, with k_c/mu = 109.6, 54.8,
  27.4, 15.8, 9.40, 7.75 and 8.78 g^-1 cm^2.
- **Self-broadened:** the amount is
  w_c = q (q p / 0.622 / 1013 hPa) exp[1800 (1/theta - 1/296)] dp/g.
  That is quadratic in vapour and stronger in the cold.

**Units.** Chou's k/mu carries his diffusivity, 1.66, and the solver applies
its own. So the vertical optical depth is (k/mu)/1.66 times the amount, with
g^-1 cm^2 times g cm^-2.

**Where other gases share a band.**

- **The overlap:** transmittances multiply (Chou eq. 4.27), the
  random-overlap assumption. It is represented exactly by sub-dividing:
  - in band 3, each CO2 term's sub-band (`co2-band.md`) is split again by the
    water terms of 3a, 3b or 3c;
  - in band 5, each ozone term's sub-band (`ozone-infrared-band.md`) is split
    by band 5's water terms.
  So each bin carries one water term and one CO2 or ozone term, and the bin's
  share of the band is the product of their fractions.
- **Elsewhere:** one sub-band per water term, sized by Delta-g against the
  Planck flux at 250 K, as the CO2 band's terms were.
- **Terms with Delta-g = 0** are left out.

**Grid.**

- **The longwave:** everything at wavelengths of 4 microns or more is rebuilt
  on these sub-bands, from 2500 down to 10 cm^-1.
- **The shortwave:** everything below 4 microns is kept exactly, with the
  solar vapour and ozone, and carried as the CO2 band is.

**Declared approximations.**

- 0-10 cm^-1 is left out.
- Band 9 is cut at 2500 cm^-1, the model's shortwave boundary. Line-by-line
  puts 5.13 W m^-2 of MLS OLR in all of 1900-3000.
- The k-distribution is Chou's "LOW" option for bands 1, 2 and 8, where his
  "HIGH" option uses tables. He reports it underestimates cooling above
  10 hPa.

**It replaces the envelope entirely.** It replaces the fitted wing decay, the
flat continuum, and the window caricature that the CO2 band already removed.
It is an option, `chou_water_vapour`, off by default, and it requires the CO2
band, which shares band 3 with it.

## The tests

On two atmospheres like for like with Chou's Tables 14 and 15:

- AFGL mid-latitude summer, 1b;
- AFGL sub-arctic winter, 1e.

Each runs to 1 hPa, with CO2 at 350 ppmv, Chou's 4e-6 g/g of water above the
tropopause, ozone's 9.6-micron band on the AFGL ozone, and no clouds.

| test | passes if |
|---|---|
| V1: total OLR | within 3 W m^-2 of line-by-line (MLS 293.99, SAW 204.37) |
| V2: total downward flux at the surface | within 4 W m^-2 (MLS 339.59, SAW 162.01) |
| V3: each of the nine bands, upward at the top | within 3 W m^-2 |
| V4: each of the nine bands, downward at the surface | within 4 W m^-2 |
| V5: maximum MLS cooling above 100 hPa | between -7 and -14 C/day, as C2 |

- **How bands are split.** Model fluxes go into Chou's bands exactly, because
  the grid now has his edges.
- **Band 9** is judged against line-by-line minus nothing. Its cut above
  2500 cm^-1 is part of the error and is reported.
- **Reported:**
  - the forcing per CO2 doubling with no fitted number anywhere in the
    longwave;
  - the cost, as the number of bins against the 42 now.

**What follows,** each under its own protocol:

- **If V1-V5 pass on both atmospheres:** the global reference column with this
  vapour and the global cloud field. The clouds then carry their own
  longwave effect, and that column's headline error is a prediction with no
  longwave calibration. After it, the tropical column's attempt 3.
- **If any fails:** the failing bands are recorded, and the option is not
  used.
