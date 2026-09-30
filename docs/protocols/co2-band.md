# Protocol: a CO2 band that is physics in the middle atmosphere

**Frozen before the band is written and before any run.** Results are appended
below a line at the end; nothing above it changes afterwards.

## Why

The tropical column of `tropical-column.md` cannot produce a valid run: its
stratopause heats to 355 K. Measured like for like on the AFGL mid-latitude
summer (MLS) atmosphere:

| | model | line-by-line |
|---|---:|---:|
| ozone's ultraviolet absorption | 13.3 W m^-2 | 14.0 W m^-2 |
| longwave cooling at the stratopause | 5.8 C/day | about 12 C/day |

- **Why the cooling is half.** The reference column has no 15-micron CO2 band.
  Its only non-vapour infrared absorber is
  `window_absorption(kappa_opaque=1e-4)`, a flat caricature.
- **The calibrated band will not do instead.** `co2_band_absorption` has
  pressure-independent absorption. So it cools the MLS stratopause at 88 C/day
  on its own 86-bin grid and 277 C/day on the 16-band one.
- **What is missing.** Real lines narrow as the pressure falls, and a band that
  ignores that is far too emissive in thin air.

## The band

Chou, Suarez, Liang and Yan (2001), *A Thermal Infrared Radiation
Parameterization for Atmospheric Studies*, NASA/TM-2001-104606 vol. 19, a
public-domain US Government work, section 4.5 and Tables 3 and 10. It was
derived from line-by-line calculations with the 1996 HITRAN, and nothing is
adjusted here.

- **Band 3, 540-800 cm^-1,** split into its centre (620-720 cm^-1, sub-band 3b)
  and its wings (540-620 and 720-800 cm^-1, sub-bands 3a and 3c).
- **Six k-terms each**, k_n = k_1 x 8^(n-1). In (cm-atm)^-1, k_1 / mu is
  2.66e-3 at the centre and 2.66e-5 in the wings, with the Planck-weighted
  fractions of Table 10:

  | term | centre | wings |
  |---:|---:|---:|
  | 1 | .0766 | .1395 |
  | 2 | .1372 | .1407 |
  | 3 | .1189 | .1549 |
  | 4 | .0335 | .1357 |
  | 5 | .0169 | .0182 |
  | 6 | .0059 | .0220 |

  Together they sum to 1.
- **The diffusivity.** Chou's 1 / mu is 1.66, the solver's own, so each term's
  vertical optical depth is (k_n / mu) u / 1.66.
- **The scaled amount** (their eq. 4.1 and Table 3):
  u = u_CO2 (p / p_r)^m [1 + alpha (T - 250) + beta (T - 250)^2].

  | | p_r | m | alpha | beta |
  |---|---:|---:|---:|---:|
  | centre | 30 hPa | 0.85 | 0.0042 | 2.00e-5 |
  | wings | 300 hPa | 0.50 | 0.0179 | 1.02e-4 |

  The wings take 3a's scaling, as Chou does when combining them. u_CO2 is the
  layer's CO2 in cm-atm at STP.
- **How it enters the grid.** One sub-band per term, sized to the term's share
  of the sub-band's Planck flux at 250 K, as the solar vapour's terms were
  sized to the sun.
  - The centre's six fill 620-720 cm^-1.
  - The wings' six fill 540-620 and 720-800 cm^-1, each split by the same
    shares.
  - The vapour's own absorption is added in each sub-band: Chou's
    multiplication of transmittances, since the vapour is smooth within one.
- **A check before the runs:** the centre's fractions sum to 0.389. The share of
  band 3's Planck flux at 250 K that lies in 620-720 cm^-1 should be close to
  that.

**What it replaces.** The window caricature, entirely: with a CO2 band, the
caricature's stand-in for CO2 would be counted twice.

**Absent, as declared:**

- ozone's 9.6-micron band;
- methane, N2O, CFCs and the minor CO2 bands. Chou's Table 16 puts them at
  -4.9 W m^-2 of MLS outgoing longwave.

**The calibration left.** The CO2 band is no longer fitted, so the forcing per
doubling becomes a prediction. Only the vapour wing decay stays fitted: to
238.175 W m^-2 on the calibration profile, by the unchanged method, with this
CO2 band as its background. It is refitted on each grid used.

## The tests

On the MLS atmosphere to 0.01 hPa, CO2 at 350 ppmv, like for like with Chou et
al. (2001):

| test | line-by-line | passes if |
|---|---|---|
| C1: TOA upward flux in 540-800 cm^-1 (CO2 and water vapour) | 67.90 W m^-2 (their Table 14) | within 3 W m^-2 |
| C2: maximum longwave cooling above 100 hPa | about -12 C/day, all gases, including the ozone band absent here (their section 10) | between -7 and -14 C/day |

Chou says their k-distribution version "greatly underestimates" cooling above
10 hPa. So C2 may fail by construction, and if it does, that is the result.

**Reported, not judged:**

- the refitted vapour wing decay;
- the model's forcing per doubling, at the top and at the tropopause, as a
  prediction to set beside the 3.70 W m^-2 the old band was fitted to;
- the MLS total outgoing longwave, against their 293.99;
- the global column of `model-top.md` (alpha 0.2, 1 hPa top) with this band:
  its temperatures at 60, 27, 10 and 1.5 hPa, and whether the run is valid.

If C1 and C2 pass, the tropical column of `tropical-column.md` is run under its
own frozen protocol. Adopting the band as the reference column's default is a
decision for after that, and it moves the headline.
