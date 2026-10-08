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
