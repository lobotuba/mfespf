# Protocol: RRTMG's longwave gas optics

**Frozen before the source is ported, before any benchmark is computed and
before any run.** Results are appended below a line at the end; nothing above
it changes afterwards.

## Why

`global-column-chou-longwave.md` R1 failed G4. Under the user's decision
recorded there before the result, the next step is RRTMG's longwave gas
optics.

**What RRTMG is.** It is AER's correlated-k model, derived from LBLRTM
line-by-line spectra with the terms ordered consistently across levels
(Mlawer et al. 1997, *JGR* 102, 16663-16682; Iacono et al. 2008, *JGR* 113,
D13103).

**What it should address:**

- the 1.6 W m^-2 that `transmittance-integration-diagnostic.md` attributed to
  integrating Chou's coefficients layer by layer, since RRTMG's coefficients
  are built to be integrated that way;
- the water-CO2 overlap, which RRTMG treats with a binary species parameter
  rather than random overlap;
- the gases Chou's implementation here leaves out: methane, N2O, the CFCs,
  and the O2 and N2 continua.

**What it is not expected to fix.** It is not expected to fix the column's
excess water. `chou-column-greenhouse-attribution.md` measures that
separately.

## The source

**AER-RC/RRTMG_LW** on GitHub, commit `958edb1ad1618d76d0fcc5a9dc2e3c2f791a54dd`,
BSD 3-Clause licence:

- `data/rrtmg_lw.nc`, the absorption coefficient tables, SHA-256 recorded in
  `data/sources.yaml` before use;
- `src/rrtmg_lw_setcoef.f90`, `rrtmg_lw_taumol.f90` and `rrtmg_lw_init.f90`,
  including the 256-to-140 g-point reduction `cmbgb`;
- `modules/rrlw_*.f90`, for the reference atmosphere, band edges and
  constants;
- `run_examples_std_atm/{input,output}_rrtm_{MLS,SAW,TROP,MLW}-clr`, AER's own
  input and output for four clear-sky standard atmospheres.

**The provenance class is `model_output`.**

## The port

**`physics/rrtmg_longwave.py`.** It ports `setcoef`, all sixteen `taugbN`,
and `cmbgb`, line for line.

- **The inputs:** layer pressure, temperature, and the gas amounts in
  molecules cm^-2. These are H2O, CO2, O3, N2O, CH4 and O2, plus the CFCs if
  configured.
- **The outputs:** the optical depth of each layer at each of the 140
  g-points, and RRTMG's Planck fraction for each.
- **What is not touched.** Nothing in the tables is fitted, adjusted or
  edited.

**The grid.**

- **The longwave.** Everything at 4 microns and longer, 10-2500 cm^-1, is
  rebuilt on RRTMG's bands 1-15. Each g-point is a sub-interval of its band,
  sized by its Planck fraction.
- **Band 15 (2380-2600 cm^-1)** is cut at 2500, the model's shortwave
  boundary.
- **Band 16 (2600-3250 cm^-1)** lies below 4 microns and is left out. This
  is the same cut Chou's band 9 had.
- **Below 4 microns** the grid is kept and carried unchanged: the solar
  vapour k-distribution, ozone's solar bands and the albedo.

**The fixed-fraction approximation.** This model's solver ties a bin's
emission to its wavelength edges, whereas RRTMG's Planck fractions vary with
the layer's state. Each g-point's sub-interval is therefore sized by a fixed
fraction:

- **the fraction:** RRTMG's own fraction, averaged over the layers of the
  AFGL US Standard Atmosphere (AFGL 1986, Table 1a), weighted by layer mass;
- **why that atmosphere:** it is neither benchmark atmosphere and nothing
  in the model's runs;
- **the declared approximation:** the size of this choice is measured by P1.

**Configuration.**

- **The option:** `rrtmg_longwave`, off by default. It replaces
  `chou_water_vapour`, the CO2 band and ozone's infrared band in the
  longwave, and it requires radiative vapour and ozone.
- **The radiation module** takes the spec and evaluates the optical depth
  every step.
- **A new equation record,** `EQ-RRTMG-LW-GAS-OPTICS`, is registered before
  it is cited.

**Trace gases are configuration parameters.** They default to 1850, to match
the column's 284.7 ppmv CO2: CH4 808.2 ppb, N2O 273.0 ppb (Meinshausen et al.
2017, *GMD* 10, 2057, CMIP6 historical), CFCs zero, and O2 0.209 by volume.

## The tests

**P0: is the port faithful?** A non-scattering transfer, a line-for-line port
of RRTMG's `rtrnmc` that lives in the check rather than in the model, is
driven by the ported optics on AER's four `input_rrtm_*-clr` atmospheres.

- **Passes if,** on every atmosphere, every level's total upward and downward
  flux is within 0.1 W m^-2 of AER's `output_rrtm_*-clr`, and every band's
  top-of-atmosphere upward flux is within 0.05 W m^-2.
- **Its nature.** This is a check of transcription, not of physics. A failure
  is a porting defect: it is corrected and P0 is re-run until it passes, and
  every correction is recorded. Nothing in the tables or tolerances may
  change to pass it.

**P1: what does this model's solver cost?** On AER's four atmospheres, this
model's solver runs with the ported optics, the fixed fractions and the
emissivity of AER's input. It is compared with AER's output.

- **Passes if** the OLR is within 1.5 W m^-2 and the surface downward flux
  within 2.0 W m^-2.
- **Reported:** the same with RRTMG's per-layer fractions in the offline
  transfer. That separates the solver from the fixed fractions.

**V1-V5: the benchmark against line-by-line.** These are
`water-vapour-longwave.md`'s, on the same atmospheres:

- **The atmospheres:** AFGL mid-latitude summer and sub-arctic winter in 1 km
  layers to 1 hPa, with CO2 at 350 ppmv and 4e-6 g/g of water above the
  tropopause.
- **The gases:** ozone from AFGL. **CH4, N2O and the CFCs are zero,** because
  Chou's line-by-line (Tables 14 and 15) carried only H2O, CO2 and O3. The
  O2 and N2 continua stay in, and that difference is declared.
- **The limits:**
  - V1: OLR within 3 W m^-2;
  - V2: the surface downward flux within 4;
  - V5: the maximum MLS cooling above 100 hPa between -7 and -14 C/day.
- **V3 and V4** (band by band) are not applicable: RRTMG's band edges are not
  Chou's. Each RRTMG band's flux is reported instead.

**G0-G4: the global column, run R3.** This is `global-column-chou-longwave.md`'s
configuration exactly:

- 51 layers to 1 hPa, 200000 steps;
- both saturation options;
- the global ISCCP cloud field, solar vapour absorption, measured spectrum,
  solar ozone, CO2 at 284.7 ppmv, and eddy diffusion 0.39;
- `rrtmg_longwave` in place of `chou_water_vapour`.

The tolerances are G0-G4 as frozen there. G4 is against R0's +1.59 K.

- **Reported:** everything that protocol reported, against R0, R1 and R2.
- **The state is saved:** the tail-mean temperature, vapour, mass, pressure
  and surface temperature.

## Order and outcomes

**The order:** P0, then P1, then V1-V5, then R3. A later step runs only
after the earlier ones are recorded.

- **If P0 cannot be passed,** the port is abandoned and recorded.
- **P1 and V1-V5 failures are recorded,** and R3 still runs, as an
  experiment the user can weigh.

**Adoption** is the user's decision on these results.

---

## P0 result (appended; nothing above this line changed)

**The port is `physics/rrtmg_longwave.py`.**

- **Tables.** `tools/gen_rrtmg_tables.py` generates them from `rrtmg_lw.nc`
  into `physics/data/rrtmg_lw_reduced.npz`: 129 arrays and 146631 values.
- **The reference atmosphere and Planck functions** are the netCDF's copies.
  They equal the values in `rrtmg_lw_setcoef.f90` to the last bit, checked
  once against a parse of the Fortran.
- **A g-point set the netCDF also carries.** Set 2 is the reduced set, but
  it is empty (zeros and fill values). So the reduction is ported, as
  `cmbgb1`-`cmbgb16` do it.
- **One Fortran quirk is reproduced, not corrected.** Band 16's upper
  atmosphere multiplies its table index by `nspb(16) = 0`, so it reads the
  table's first entries at every level.

**P0 passes on all four atmospheres, at the first run of the complete port.**

| atmosphere | layers | OLR, port / AER | surface down, port / AER | worst level, up / down |
|---|---:|---:|---:|---:|
| MLS | 51 | 281.5358 / 281.5358 | 348.5429 / 348.5429 | 5.0e-5 / 5.0e-5 |
| SAW | 59 | 199.4693 / 199.4693 | 172.4107 / 172.4107 | 5.0e-5 / 4.7e-5 |
| TROP | 60 | 287.5743 / 287.5743 | 398.1302 / 398.1302 | 4.9e-5 / 4.9e-5 |
| MLW | 60 | 230.6043 / 230.6043 | 223.9858 / 223.9858 | 4.8e-5 / 5.0e-5 |

- **The levels.** All are within 5e-5 W m^-2 of AER's output, which is
  that output's rounding; the limit was 0.1.
- **The bands.** The 16 MLS bands at the top, the only atmosphere AER
  printed by band, are within the same rounding; the limit was 0.05.

**Two corrections were made while writing the check, both in its own
parsing, not in the port.**

- **The layer records are fixed-width.** Some run their fields together, so
  they are read by the Fortran format's columns.
- **The output header's band range** is read with a pattern.

**The check is now `tests/unit/test_rrtmg_longwave.py`,** with AER's eight
files in `tests/data/aer_rrtmg_lw/` and registered in `data/sources.yaml`.
The test file also checks that doubling the vapour raises every g-point's
optical depth in bands 1, 2 and 10.

**Next, as frozen: P1**, this model's solver with the fixed fractions.

## Erratum, before the fixed fractions are computed (nothing above this line changed)

**The protocol contradicts itself.** "The fraction" names "the AFGL US
Standard Atmosphere (AFGL 1986, Table 1a)". In AFGL 1986, Table 1a is the
Tropical atmosphere (surface 299.7 K); the US Standard is Table 1f (surface
288.2 K).

**Both readings:**

- **Table 1a, Tropical.** This is against the protocol's own stated reason,
  "neither benchmark atmosphere and nothing in the model's runs": P1 runs
  AER's tropical atmosphere.
- **Table 1f, US Standard.** This is the atmosphere the protocol names, and
  it satisfies its reason.

**Taken: Table 1f.**

- **The source:** the rayference/afgl1986 transcription at commit
  9d170bb3e52b, the one already used for 1b and 1e (`tables/1f.csv`, SHA-256
  d2d06777c7fbb9ed5bde5bf262a854ff61f9e83ea4d9b771006f10b87468eeda).
- **How the fractions are formed:**
  - layers between consecutive levels, at the mean of their levels'
    temperatures and mixing ratios;
  - H2O, O3, N2O, CO and CH4 as tabulated;
  - CO2 at the model's 284.7 ppmv, since the table has no CO2 column, and O2
    at 0.209;
  - the mass weight is the pressure difference.
- **Within each band** the mean fractions are renormalised to sum to one. In
  bands 12 and 15, RRTMG gives upper-atmosphere layers zero fraction along
  with zero optical depth.

## P1 result (appended; nothing above this line changed)

**The fixed fractions** are generated by `tools/gen_rrtmg_fractions.py` on
Table 1f (49 layers). Each band's mean sums to one before renormalising,
except bands 12 and 15 at 0.913: RRTMG zeroes their upper atmosphere.

**The grid.** `rrtmg_longwave_grid` puts one bin per g-point, each a
sub-interval of its band sized to its fixed fraction at 250 K, in g-point
order.

**The comparison.**

- **The model.** P1 runs on that grid built to 2600 cm^-1, so the coverage
  is bands 1-15 complete. The model's flux is summed over those bins only.
- **The reference.** AER's total minus band 16, with band 16 taken from the
  offline transfer, which P0 showed equals AER.
- **One correction to the first P1 run.** It also summed the transparent
  bins shortward of 2600 cm^-1, which the reference does not contain. That
  was a bookkeeping error, corrected before anything was recorded.

**P1 fails on three of four atmospheres.** All values in W m^-2:

| atmosphere | OLR, model / AER bands 1-15 | difference | surface down, model / AER | difference | P1 |
|---|---:|---:|---:|---:|---|
| MLS | 278.90 / 281.17 | -2.27 | 342.34 / 348.35 | **-6.00** | **fail** |
| SAW | 198.36 / 199.41 | -1.06 | 173.14 / 172.39 | +0.75 | pass |
| TROP | 283.85 / 287.12 | **-3.27** | 398.31 / 397.84 | +0.47 | **fail** |
| MLW | 228.85 / 230.47 | **-1.62** | 224.07 / 223.94 | +0.13 | **fail** |

**Reported as frozen: the fixed fractions in AER's own transfer.** Against
RRTMG's per-layer fractions they change OLR by +0.64 (MLS), +0.25, +0.77 and
+0.34, and the surface downward flux by -0.03 to -0.07.

**Reported beyond the protocol: where the gap is.** The offline transfer is
switched, one choice at a time, to what this model's solver does.

- (d) **The diffusivity:** 1.66 in every band, where RRTMG uses 1.5-1.8 by
  band and column water.
- (s) **The source:** each layer emits at its mean temperature, where RRTMG
  varies the source linearly in optical depth between the level temperatures.
- (p) **The emission:** each g-point as a fixed wavelength sub-interval
  emitting that sub-interval's Planck flux, where RRTMG uses its fraction
  times the band's Planck flux. Its share of the band then varies with
  temperature.

With all three, the offline transfer reproduces the model's solver to 0.05
(MLS -2.22/-6.00 against -2.27/-6.00). So the decomposition is complete:

| atmosphere | (d) OLR / down | (s) OLR / down | (d)+(s) OLR / down | (p) OLR / down |
|---|---:|---:|---:|---:|
| MLS | -0.60 / +0.24 | +0.53 / **-7.10** | -0.07 / -6.91 | **-2.15** / +0.90 |
| SAW | -0.16 / +0.49 | +0.04 / +0.14 | -0.12 / +0.62 | **-0.94** / +0.13 |
| TROP | -0.97 / +1.12 | +0.11 / -1.57 | -0.86 / -0.47 | **-2.40** / +0.94 |
| MLW | -0.23 / +0.07 | +0.05 / -0.45 | -0.18 / -0.39 | **-1.44** / +0.51 |

**What it says.** The gas optics are RRTMG's exactly (P0). The gap is the
solver.

- **The OLR deficit is mostly (p).** That is the representation of a
  k-distribution as wavelength sub-intervals.
  - At warm temperatures the Planck function shifts toward each band's high
    wavenumbers, where the sub-intervals of the strongly absorbing g-points
    sit, so too much emission comes from opaque g-points.
  - Chou's terms were built the same way, so this is a candidate for the
    1.6 W m^-2 that `transmittance-integration-diagnostic.md` attributed to
    the solver. That is not measured.
- **The MLS surface deficit is (s).** It is the isothermal-layer source in
  AER's thick (about 1 km), moist lower layers.
- **The diffusivity matters less,** and partly offsets the others.

**Both are solver physics, not gas optics, each with a published
formulation.**

- **(p):** a source weight per g-point and layer, RRTMG's fraction times the
  band's Planck flux.
- **(s):** a linear-in-tau source (Clough et al. 1992, *JGR* 97, 15761).

Neither is built here; each would need its own protocol. **As frozen, P1's
failure is recorded and V1-V5 follow.**

## V1-V5 result (appended; nothing above this line changed)

**The setup is `water-vapour-longwave.md`'s exactly.**

- **The atmospheres:** AFGL MLS and SAW, 1 km layers to 1 hPa.
- **The gases:** CO2 350 ppmv, 4e-6 g/g of water above 100 hPa, ozone from
  AFGL, CH4, N2O, CO and the CFCs zero, and O2 0.209.
- **The model:** this model's solver on the model's grid, the longwave
  rebuilt to 2500 cm^-1. There are 145 bins, 138 of them longwave.

All values in W m^-2:

| atmosphere | V1: OLR | Chou line-by-line | difference | V2: surface down | Chou line-by-line | difference | V5: max cooling above 100 hPa |
|---|---:|---:|---:|---:|---:|---:|---:|
| MLS | 281.96 | 293.99 | **-12.03** | 342.52 | 339.59 | +2.93 | -13.55 C/day at 1.1 hPa |
| SAW | 200.06 | 204.37 | **-4.31** | 172.27 | 162.01 | **+10.26** | (-26.1) |

**Verdict as frozen.**

- **V1 fails** on both atmospheres.
- **V2 passes on MLS and fails on SAW.**
- **V5 passes.**
- **V3 and V4** are not applicable, as frozen.

**Reported beyond the protocol: the two references disagree.** RRTMG through
AER's own transfer (the offline `rtrnmr`, equal to AER to 5e-5 by P0), on
exactly these columns, nothing of this model's solver:

| atmosphere | OLR, all 16 bands | minus Chou line-by-line | surface down, all 16 bands | minus Chou line-by-line |
|---|---:|---:|---:|---:|
| MLS | 284.51 | **-9.48** | 348.09 | **+8.50** |
| SAW | 200.86 | -3.51 | 170.54 | **+8.53** |

**What this says.**

- **RRTMG disagrees with Chou's line-by-line by more than the V limits,
  before this model is involved.** RRTMG's atmosphere is more opaque by about
  9 W m^-2 at both ends on MLS.
  - Chou's line-by-line was computed for the 2001 memorandum, with the
    spectroscopy and water-vapour continuum of that time.
  - RRTMG derives from LBLRTM with the MT_CKD continuum, validated against
    surface spectral measurements (Mlawer et al. 2012, *Phil. Trans. R.
    Soc. A* 370, 2520).
  - Which is nearer the truth is not settled by this benchmark, which uses
    one of the two as truth.
- **The V limits were set for Chou's scheme against Chou's own
  line-by-line.** Against RRTMG they mostly measure the disagreement between
  references. This model's solver adds P1's -1 to -3 to the OLR.
- **The consequence for R3.** RRTMG makes a given column about 3 W m^-2 more
  opaque than Chou's scheme did on MLS:
  - RRTMG through this solver gives 281.96;
  - Chou attempt 1 through this solver gave 287.80.

  R1's warm surface came with a correct top-of-atmosphere flux and a
  column 2.8 times too wet. A more opaque gas scheme on the same humidity
  closure is expected to warm R3 further, not less. This is stated before
  R3 runs.

**As frozen, R3 follows:** the wiring, its checks and the global run.

## R3 result (appended; nothing above this line changed)

**R3** is `global-column-chou-longwave.md`'s configuration with
`rrtmg_longwave`.

- **The run:** 178 bins, 200000 steps and 18624 s.
- **Its first attempt** was stopped at about three hours for memory, while
  R1's full rerun ran, and restarted from scratch.
- **G0 passes.** The run is complete and valid, with residual 7.9e-16,
  imbalance 0.007 W m^-2 and no domain violations.

Against CERES EBAF Ed4.0, global:

| test | observed | R0 | R1 (Chou + O3 band) | **R3 (RRTMG)** | verdict on R3 |
|---|---:|---:|---:|---:|---|
| G1: planetary albedo | 0.291 | 0.284 | 0.277 | **0.268** | pass (0.023) |
| G2: shortwave cloud effect | -45.8 | -52.5 | -50.9 | **-49.5** | pass (3.7) |
| G3: longwave cloud effect | +27.9 | +40.0 | +23.5 | **+14.5** | **fail** (13.4) |
| G4: surface, against 288 K | | +1.59 K | +9.91 K | **+14.65 K** | **fail** |

**R3 fails G3 and G4.**

**Reported:**

| | observed | R0 | R1 | R3 |
|---|---:|---:|---:|---:|
| outgoing longwave | 240.1 | 243.6 | 246.2 | 248.9 W m^-2 |
| clear-sky outgoing longwave | 268.1 | 283.6 | 269.6 | **263.4** W m^-2 |
| precipitable water | ~25 | 32.4 | 70.6 | **109.6** kg m^-2 |
| relative humidity, 300-700 hPa | 0.4-0.5 | 0.79 | 0.89 | 0.91 |
| latent heat flux | | 82.9 | 114.1 | 128.8 W m^-2 |
| sunlight absorbed in the atmosphere | | 58.2 | 69.5 | 76.9 W m^-2 |
| cold point | | 187.8 K | 202.8 K | 203.6 K at 57 hPa |
| stratospheric water | | 2.3 | 26.8 | 44.4 ppmv |

**The cloud tops** are 252.1, 262.8, 283.3 and 296.3 K, against ISCCP's
227.5, 234.8, 262.8 and 281.1 K.

**CO2 doubling.** Instantaneous and thermal, on the tail-mean state with its
clouds: 2.57 W m^-2 at the top, 4.07 at 203 hPa.

**The 3.16-4.0-micron bin** carries 2.89 W m^-2 of sunlight at the top and
2.50 at the surface.

**What it says.** The warming was predicted in the V-test record, before R3
ran, and this is it.

- **RRTMG is more opaque on the same water,** and the humidity closure
  answers with more water: 109.6 kg m^-2, 4.4 times the observed.
- **The cloud effect falls.** The surface warms 4.7 K beyond R1, and the
  cloud tops warm with the atmosphere. The longwave cloud effect falls to
  14.5, because a warm cloud top traps little, as `observed-cloud-field.md`
  found.
- **The gas optics are the best-validated in the model.** They are exact to
  AER (P0). The failure is not theirs.
  `chou-column-greenhouse-attribution.md` measured where it is: the column's
  water, which carries 39.7 W m^-2 of excess clear-sky greenhouse effect on
  R1's state.

**Adoption is the user's decision.** As frozen, nothing changes the
reference column.
