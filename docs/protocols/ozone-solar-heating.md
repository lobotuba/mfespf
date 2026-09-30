# Protocol: ozone absorbs sunlight, and the stratosphere warms upward

**Frozen before the mechanism is written and before any run.** Results are
appended below a line at the end; nothing above it changes afterwards.

## Why

`stratospheric-water.md` gave the reference column one equilibrium by letting
vapour diffuse above the convective top (+9.93 K from any start). Its
stratospheric-water prediction failed: 17.1 ppmv against SWOOSH's 4.0.

- **No inversion.** The column's temperature falls all the way to its top, 190 K
  at 27 hPa. So the top is the cold point and every layer above 130 hPa sits
  saturated.
- **What is observed.** The NCEP/NCAR R1 1991-2020 global mean has its minimum
  at 100 hPa (205.7 K) and warms upward: 207.1 K at 70 hPa, 211.2 at 50, 216.8
  at 30, 228.8 at 10. The column is 27 K too cold at 27 hPa.
- **Why.** Ozone absorbs sunlight, and this column has none. The solar
  absorption record says "no ozone, oxygen or carbon dioxide".

## The physics

**The sun's spectrum, measured.** The solar source is a 5772 K black body. That
puts 12.1 percent of the sunlight between 0.175 and 0.4 microns, against 8.0
percent in the real spectrum (Chou and Suarez 1999, Table 3). Ozone's heating is
mostly there, so it would be half again too strong on a black-body sun.

- **Replacement:** the TSIS-1 Hybrid Solar Reference Spectrum, version 2, full
  spectrum, 115.5 nm to 200 microns (Coddington et al. 2023,
  doi:10.1029/2022EA002637; LASP LISIRD `tsis1_hsrs_binned_fs`). It integrates
  to 1361.49 W m^-2.
- **How it is used.** It enters as the solar source's spectral shape, scaled to
  the configured total irradiance, as the black body is. The k-distribution
  sub-bands of the vapour's solar absorption are placed on its cumulative
  spectrum instead of the black body's.
- **In the repository.** A reduced table is committed: the cumulative spectrum
  at 0.5 percent logarithmic steps. The full file is declared in
  `data/sources.yaml`. The reduction's error on this grid's bands is measured
  and stated with it.
- **Checked before this freeze, so a consistency check and not a prediction.**
  HSRS integrated over Chou and Suarez's bands reproduces their Table 3 S/S0:
  - 0.00350 against 0.00367 (bands at 225-245 and 260-280 nm);
  - 0.0596 against 0.0591 (320-400 nm);
  - 0.394 against 0.391 (400-700 nm);
  - 0.320 against 0.321 (700-1220 nm).

  The black body gives 0.121 against 0.080 below 0.4 microns.

**Ozone's absorption.** Chou and Suarez (1999), section 3.5 and Table 3, give
effective ozone absorption coefficients in (cm-atm)^-1 at STP. They derived
them from WMO (1986) cross sections over 127 sub-bands.

| band | wavelengths (microns) | coefficient |
|---|---|---:|
| 1 | 0.175-0.225 | 30.47 |
| 2 | 0.225-0.245 and 0.260-0.280 | 187.24 |
| 3 | 0.245-0.260 | 301.92 |
| 4 | 0.280-0.295 | 42.83 |
| 5 | 0.295-0.310 | 7.09 |
| 6 | 0.310-0.320 | 1.25 |
| 7 | 0.320-0.400 | 0.0345 |
| 8 | 0.400-0.700 | 0.0572 |

- Band 8 already folds in ozone's weak near-infrared absorption (their eq.
  3.14).
- No temperature or pressure scaling is applied, as in the source.
- The grid's edges below 0.4 microns are replaced by these, band 2 as two
  bins. Nothing below 0.175 microns absorbs; 0.03 percent of the sunlight is
  there.

**Ozone's amount, prescribed from observation.** SWOOSH v2.73 combined,
anomaly-filled, 10-degree zonal means, 1991-2020, area-weighted global mean,
316 to 1 hPa, computed here. It is 0.061 ppmv at 316 hPa, 0.462 at 100, 3.98
at 31.6 and 7.98 at 8.25.

- Each layer takes the mixing ratio interpolated in log pressure at its centre.
- Below 316 hPa the mixing ratio is held at the 316 hPa value.
- The ozone above the model top (10 hPa), 50 DU, is placed in the top layer.
- The total column is about 297 DU.

**Declared approximations**, each a known way this can be wrong:

- Ozone's 9.6-micron longwave band is not represented. It warms the lower
  stratosphere by absorbing upwelling longwave, so its absence biases that
  region cold.
- Rayleigh scattering is not added.
- The heating that really happens above 10 hPa is deposited in the top layer.
- The ozone amount is prescribed, not predicted.

## The tests

Reference column with water vapour's solar absorption, eddy diffusion of vapour
at 0.39 m^2 s^-1, the measured solar spectrum and ozone. 30 layers, 3600 s,
200000 steps. Equilibrated as in `stratospheric-water.md`.

| test | observed | passes if |
|---|---|---|
| O1: temperature at the 60 hPa layer | 209.0 K (NCEP R1 global mean, interpolated in log p) | within 8 K |
| O1: temperature at the 27 hPa layer | 218.0 K | within 8 K |
| O2: the coldest layer | the global mean is coldest at 100 hPa and warms above | the coldest layer is centred between 60 and 160 hPa |
| O3: stratospheric water, mean of the layers above 100 hPa | 4.0 ppmv (SWOOSH, as in `stratospheric-water.md`) | 2.0 to 8.0 ppmv |

Reported, not judged:

- the ozone's solar absorption;
- the atmosphere's total absorbed sunlight (observed about 78 W m^-2; Trenberth,
  Fasullo and Kiehl 2009);
- the surface-temperature error;
- the same column with the measured spectrum and no ozone, so that the
  spectrum's effect and ozone's are separated.

The observed surface temperature decides nothing here.
