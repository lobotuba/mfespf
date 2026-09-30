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

---

## Results (appended; nothing above this line changed)

### Attempt 1: fails on every criterion, and cannot judge ozone

200000 steps, last fifth. Both runs are equilibrated and closed to 1e-15.

| | ozone (verdict) | measured spectrum, no ozone |
|---|---:|---:|
| surface | 319.08 K (+31.1 K) | 313.45 K (+25.4 K) |
| net TOA | +0.025 | -0.001 |
| precipitable water | 342 kg m^-2 | 248 kg m^-2 |
| sunlight absorbed in the atmosphere | 106.3 W m^-2 | 88.4 W m^-2 |
| reflected sunlight | 63.1 W m^-2 | 69.9 W m^-2 |
| temperature at 60 hPa (O1: 209.0 +/- 8) | 256.7 K, **fail** | 239.7 K |
| temperature at 27 hPa (O1: 218.0 +/- 8) | 237.5 K, **fail** | 209.1 K |
| coldest layer (O2: 60-160 hPa) | the top, 27 hPa, **fail** | the top |
| stratospheric water (O3: 2-8 ppmv) | 16930 ppmv, **fail** | 4094 ppmv |

The whole column ran away to a hot, wet state: 290 K and 5 percent water by
volume at 328 hPa. Diagnostic runs after the verdict isolated why (100000 steps
each):

| configuration | surface | reflected | absorbed in atmosphere |
|---|---:|---:|---:|
| measured spectrum + diffusion, no solar vapour | 297.87 K | 100.4 | 0 |
| solar vapour, black body, no diffusion | 304.31 K | 76.0 | 71.2 |
| solar vapour, measured spectrum, no diffusion | 304.82 K | 74.7 | 75.2 |
| solar vapour + diffusion, black body | 311.56 K | 72.0 | 82.1 |

- **The measured spectrum changes almost nothing.** 297.87 K against 297.93 K
  with the black body, and 304.82 K against 304.31 K.
- **Water vapour's solar absorption is the driver.** It is the confound
  `solar-vapour-absorption.md` measured and declared: with the 0.30 albedo at
  the surface, reflected sunlight crosses the absorbing vapour twice. The
  planetary albedo falls by 25 W m^-2 worth (reflected 100 -> 72-76), and a
  column at 0.75-0.95 relative humidity amplifies that into +13 to +31 K.
- **That protocol's conclusion was explicit:** keep the vapour's solar
  absorption off until the reflection is where Earth's is.

**This protocol was frozen on that configuration anyway. That is an error in
the freeze, not a property of ozone**, and the verdict above therefore says
nothing about ozone. What the pair of runs does show is ozone's own
difference on a common base:

- +17.9 W m^-2 of sunlight absorbed in the atmosphere;
- +28 K at 27 hPa and +17 K at 60 hPa;
- +5.6 K at the surface.

The test is repeated below without the confound. That is a second freeze with
its own record, not an edit of the first.

## Attempt 2, frozen before its run

**Configuration.** Identical to attempt 1 except that water vapour's solar
absorption is off. The planetary albedo then stays at the surface value's
0.295 (as in the 297.87 K diagnostic above).

- **Ozone's heating is still tested:** ozone absorbs mostly above 100 hPa, in
  the ultraviolet the stratosphere takes before any reflection.
- **What it gives up:** the troposphere absorbs no sunlight, a known deficit
  that belongs to the clouds-at-altitude item and not to this test.

**Tests.** O1, O2 and O3 exactly as above, same thresholds, same observations,
200000 steps. The same column without ozone is the 297.87 K diagnostic, rerun
at 200000 steps as the comparison.
