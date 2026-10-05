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

---

## Results (appended; nothing above this line changed)

**Built:** `physics/water_vapour_longwave.py`, with Tables 1, 3, 4, 9 and 10 as
published.

- **The transcription is checked.** Each band's Delta-g sums to one, band 3's
  three sub-bands sum to one together, and the bands tile 0-3000 cm^-1.
- **The grid.** 168 bins, 161 of them longwave; the present grid has 42.
- **Not yet wired.** It is not connected to the radiation module.

All values in W m^-2:

| atmosphere | V1: OLR | line-by-line | V2: surface down | line-by-line | V3: worst band, up | V4: worst band, down | V5: max cooling |
|---|---:|---:|---:|---:|---:|---:|---:|
| MLS | **287.80** | 293.99 | 341.95 | 339.59 | -1.91 (340-540) | +2.74 (540-800) | -11.23 C/day |
| SAW | 202.62 | 204.37 | 161.59 | 162.01 | -0.85 (340-540) | -1.11 (540-800) | (-7.12) |
| verdict | **MLS fails** (-6.19, limit 3); SAW passes | | both pass | | both pass | both pass | passes |

MLS band by band:

| band, cm^-1 | up minus line-by-line | down minus line-by-line |
|---|---:|---:|
| 0-340 | -0.93 | -0.53 |
| 340-540 | -1.91 | -0.43 |
| 540-800 | -1.54 | +2.74 |
| 800-980 | -0.28 | +1.39 |
| 980-1100 | +0.23 | +0.32 |
| 1100-1215 | +0.28 | +0.20 |
| 1215-1380 | -0.71 | +0.15 |
| 1380-1900 | -0.42 | -1.44 |
| 1900-3000 | -0.92 | -0.03 |

**Verdict as frozen: V1 fails on MLS, so the option is not used.**

**What the failure is, reported.**

- **The improvement.** The largest MLS band error falls from -11.35 to -1.91
  W m^-2, and the total from -29.70 to -6.19. The sub-arctic winter passes
  everything. Every single band is within 3 W m^-2 up and 4 down. The total
  fails only because nearly every band is slightly low.
- **It is not the layering.** Refining MLS's 1 km layers 2, 4 and 8 times
  moves the OLR to 288.06, 288.13 and 288.15: converged, and still -5.8.
- **Refining moves V2 the wrong way.** It raises the MLS surface downward flux
  to 345.81, 347.54 and 348.29, which would fail V2 (+8.7) at a resolution
  where the frozen 1 km layering passes. Recorded.
- **Where the OLR deficit is:**
  - bands 1, 2 and 8, -3.26. There Chou's own benchmark used his
    table-driven "HIGH" option and this uses his k-distribution, the declared
    "LOW" approximation;
  - band 9, cut at 2500 cm^-1, -0.92;
  - 540-800 cm^-1, -1.54.
  
  Chou's own parameterization gives 293.03.

**What would close it, not built here.**

- **Tables for bands 1, 2 and 8.** Chou's HIGH-option water tables for these
  bands are in the same WRF source as the ozone tables (h11-h13, h21-h23 and
  h81-h83). They would be fitted as pressure-tabulated k-distributions, as the
  ozone band was.
- **Band 9 above 2500 cm^-1.** It needs the longwave to extend past the
  model's 4-micron boundary.
- **Each needs its own protocol.**

## Attempt 2, frozen before the tables are extracted or fitted

**The change.** In bands 1 (0-340 cm^-1), 2 (340-540) and 8 (1380-1900), the
water vapour's line absorption comes from Chou's "HIGH" option. That is the
option his Table 14 and 15 benchmarks used.

- **Source.** His pre-computed line-by-line transmission tables, h11-h13,
  h21-h23 and h81-h83. They are in WRF's `phys/module_ra_goddard.F` at the
  commit `ozone-infrared-band.md` recorded.
- **Grid (Table 8):**
  - log10 w from -8 in 31 steps of 0.3, with w in g cm^-2;
  - log10 p from -2 in 26 steps of 0.2;
  - a quadratic in theta - 250.
- **Fitted as `ozone-infrared-band.md` attempt 2 fitted the ozone table:**
  - six terms, with weights shared across pressures;
  - six coefficients at each of the 21 pressures from 0.1 to 1000 hPa;
  - a shared quadratic temperature scaling;
  - log k interpolated in log p, held above 0.1 hPa, and extrapolated
    log-linearly below 1000 hPa;
  - Chou's tables take the raw amount at an effective pressure and
    temperature, so no further pressure scaling enters.
- **The rest of the scheme is unchanged:**
  - bands 3-7 and 9 as attempt 1;
  - the continuum;
  - the overlap;
  - band 9's cut at 2500 cm^-1.

**The tests.**

- **W0, per band,** before anything else. The same criteria as the ozone
  band's Z0 and Z0b, at 250 K over the table's whole amount range:
  - RMS <= 0.01 and maximum <= 0.03;
  - RMS <= 0.02 across 190-310 K;
  - held-out odd pressures, interpolated, RMS <= 0.01 and maximum <= 0.03.
  
  A band that fails W0 keeps attempt 1's k-distribution, and that is
  recorded.
- **V1-V5 unchanged**, on MLS and SAW at the frozen 1 km layering.
- **Reported:** the same under 4x layer refinement.

**What follows:** as frozen above.

## Attempt 2, results (appended; nothing above this line changed)

**The tables are extracted.** h11-h13, h21-h23 and h81-h83 are 26 by 31 each.
Two blemishes in the source values are recorded:

- band 1's a rises by at most 4.5e-5 with amount;
- band 2's a reaches -0.00019 at saturation.

**W0 passes for all three bands:**

| band | RMS, 250 K | max, 250 K | RMS, 190-310 K | held-out RMS | held-out max |
|---|---:|---:|---:|---:|---:|
| 1 | 0.0014 | 0.0048 | 0.0034 | 0.0021 | 0.025 |
| 2 | 0.0013 | 0.0053 | 0.0015 | 0.0016 | 0.012 |
| 8 | 0.0017 | 0.0053 | 0.0026 | 0.0025 | 0.023 |

**V1 fails again, and further off.** All values in W m^-2:

| | attempt 1 | attempt 2 | line-by-line |
|---|---:|---:|---:|
| MLS OLR | 287.80 | **286.48** (-7.51) | 293.99 |
| MLS band 1, up | -0.93 | -1.12 | |
| MLS band 2, up | -1.91 | -2.68 | |
| MLS band 8, up | -0.42 | -0.78 | |
| MLS surface down | 341.95 | 342.33 | 339.59 |
| SAW OLR | 202.62 | 201.87 | 204.37 |
| SAW surface down | 161.59 | 162.88 | 162.01 |

V2-V5 still pass on both atmospheres. The MLS maximum cooling is
-11.60 C/day. **As frozen, the option is not used.**

**What this says, reported.**

- **The tables are not the gap.** Fitted to Chou's own line-by-line
  transmittances to 0.0013, they move the MLS fluxes further from his
  line-by-line fluxes. So the remaining shortfall is not the k-distribution
  against his tables.
- **The candidate that remains is how the transmittances are integrated.**
  - Chou computes band-mean transmittance between every pair of levels, at a
    Curtis-Godson effective pressure and temperature for the path, eq. 4.13
    and 4.14.
  - This model's solver treats each bin as monochromatic, layer by layer. It
    takes each layer's own pressure, so the k-terms are assumed correlated
    from layer to layer.
  - For strong lines along strongly inhomogeneous paths the two differ, and a
    fit made independently at each pressure need not keep the terms' ordering
    across pressures.
  - This is not measured.

**Code.** The tabulated bands are kept as an opt-in, `chou_longwave_grid(...,
tabulated=("1", "2", "8"))`. The default stays attempt 1, which reproduces
287.80 and 202.62.
