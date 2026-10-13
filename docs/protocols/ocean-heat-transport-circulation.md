# Protocol: the tropical overturning driven by the ocean's heat transport

**Frozen before the mechanism is written and before any run.** Results are
appended below a line at the end; nothing above it changes afterwards.

## Why

**The column's warm bias is its water.**

- `chou-column-greenhouse-attribution.md` put 39.7 W m^-2 of excess clear-sky
  greenhouse effect on it.
- `longwave-source-function.md` R4, the best-validated longwave in the model,
  still holds 98 kg m^-2 of precipitable water at a relative humidity of 0.89
  over 300-700 hPa. The observed global values are about 25 and 0.4-0.5.

**Nothing in a single column dries the free troposphere.** Real dryness comes
from air sinking outside convection, a horizontal structure. CLAUDE.md
records that the free-tropospheric humidity is not obtainable from one
column: its parameter-free limits are 0.015 and 0.725.

**Why the earlier attempt failed.** `spectral-wtg-circulation.md` built the
relaxed, spectral WTG overturning between two columns. Its Test A failed
because two identical columns never make a temperature difference. The
literature it follows prescribes a contrast, usually in sea-surface
temperature (Daleu et al. 2016).

**The user's decision, recorded before this protocol:** the contrast comes
from the ocean's heat transport, so the surface temperature contrast is
predicted rather than imposed.

## What was measured before freezing

The data are NCEP/NCAR R1 1991-2020 monthly long-term means (omega, air,
shum, skt, and the four surface fluxes). The independent check is OAFlux net
heat flux 1991-2009 (Yu & Weller 2007, *BAMS* 88, 527; ISCCP-FD radiation).
Weights are cos-latitude. SHA-256s are recorded in `data/sources.yaml` before
use. Rising means omega500 < 0.

**Global masks carry no contrast.**

- The ocean delivers +0.1 / -0.1 W m^-2 (NCEP) and +4.0 / -3.8 (OAFlux) to
  rising / sinking.
- The skin temperature contrast is -0.22 K, rising *colder*.

The global partition mixes the midlatitude storm tracks into the rising
region and the warm subtropics into the sinking one, as
`spectral-wtg-circulation.md`'s diagnostics found. So the domain is the
tropics.

**Monthly masks confound transport with storage.** On monthly masks the two
datasets disagree in sign at ±20 and ±30 degrees. Ascent follows the sun into
the summer hemisphere, where the ocean is storing heat, which reads as the
ocean taking heat from the rising region. On **annual-mean** fields and masks
storage cancels, and what remains is transport.

**Within ±20 degrees, on annual-mean masks, the two datasets agree in sign.**
The ocean carries heat into the rising region. At ±30 degrees they disagree
(OAFlux +2.4 / -2.1, NCEP ocean -2.1 / +1.8), so the domain is ±20 degrees.

| quantity, ±20°, annual means | rising | sinking | band | source |
|---|---:|---:|---:|---|
| area fraction (2.5° grid) | **0.558** | 0.442 | | NCEP omega500 |
| ocean heat delivered, W m^-2 (land delivers none) | +5.68 | -7.71 | 0 | OAFlux |
| ocean heat delivered, W m^-2 | +2.03 | -2.64 | 0 | NCEP, all surface |
| omega500, hPa/day | -20.2 | **+14.3** | -4.9 | NCEP |
| RH 300-700 hPa, pressure-weighted | 0.50 | **0.35** | **0.43** | NCEP shum, air |
| precipitable water 1000-300 hPa, kg m^-2 | 42.3 | 32.7 | **38.0** | NCEP shum |
| skin temperature | 300.26 K | 298.79 K | 299.61 K | NCEP skt |
| skin temperature, rising - sinking | | | **+1.46 K** (ocean only +2.52) | NCEP skt |
| T at 500 hPa, rising - sinking | | | **+0.65 K** | NCEP air |

**The sinking region's RH profile is C-shaped:**

| level, hPa | 1000 | 925 | 850 | 700 | 600 | 500 | 400 | 300 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| RH | 0.75 | 0.69 | 0.57 | 0.39 | 0.32 | 0.28 | 0.28 | 0.50 |

**How RH is computed.** RH is q / q_sat(T, p), over liquid at and above
273.16 K (Murphy & Koop 2005) and over ice below, as the model saturates.

**How the scripts did it.** `tropics_obs.py`, `oaflux_composite.py` and
`ocean_flux_obs.py` are recorded with the result. They interpolate omega500
bilinearly onto the flux grids. The ocean heat delivered is each region's
mean net downward surface flux, with the two regions' area-weighted mean
removed: the ocean only moves heat, and the datasets' global imbalances
differ (NCEP +2.5, OAFlux about +25 W m^-2 over the ocean).

## The mechanism

**The ocean's heat transport between the two regions, `OceanHeatTransport`.**
A prescribed lateral flux enters each column's surface, the slab-ocean
"q-flux" (Hansen et al. 1984, *Climate Processes and Climate Sensitivity*,
AGU Monograph 29, 130-163; Bitz et al. 2012, *J. Climate* 25, 3053).

- **The fluxes.** The rising column's surface receives Q_r = (1 - σ) ΔQ and
  the sinking column's Q_s = -σ ΔQ, in W m^-2 of each region.
  - σ Q_r + (1 - σ) Q_s = 0 exactly: the ocean moves heat and creates none.
  - ΔQ is the observed delivery contrast, rising minus sinking:
    - **OAFlux, 13.39 W m^-2: the verdict.** It is the purpose-built,
      validated air-sea flux product.
    - **NCEP, 4.67: sensitivity.** It is a reanalysis by-product with known
      surface shortwave biases.
- **What carries it.** A module scheduled between the regions, as the
  circulation is. It adds Q / C_surface to each region's surface-temperature
  tendency.
- **The bookkeeping.** It records the heat crossing from the sinking surface
  to the rising surface as a transfer, per unit of global area: (1 - σ) |Q_s|
  dt = σ Q_r dt. Both sides are declared, as CLAUDE.md requires.
- **Not used.** Neither the observed SST contrast nor the observed subsidence
  enters anywhere. Both are predictions.

**The two columns.** Each is `tropical_column` at ±20 degrees:

- **The band's sunlight:** the geometric factor, 0.29979.
- **The band's export:** NCEP R1's 42.51 W m^-2, from the atmosphere below
  100 hPa.
- **The Brewer-Dobson upwelling,** per unit area of the ±20° band.
- **The physics, R4's exactly:**
  - 51 layers to 1 hPa;
  - both saturation options;
  - eddy diffusion 0.39;
  - measured spectrum and solar ozone;
  - solar vapour absorption;
  - `rrtmg_longwave`, `k_distribution_source` and `linear_in_tau_source`;
  - mixed layer 1e6 J m^-2 K^-1;
  - 3600 s steps;
  - US Standard Atmosphere start.
- **The cloud field:** the observed *tropical* ISCCP field
  (`ISCCP_TROPICAL_CLOUD_TYPES`, Rossow & Schiffer 1999 Table 5, ±30
  degrees), as `tropical-cloud-field.md` adopted for the tropics.

**The circulation** is `spectral-wtg-circulation.md`'s, unchanged:

- σ = 0.558;
- L = 734 km, its measurement within 20 degrees;
- the tropopause at 100 hPa.

**Declared approximations:**

- one rising and one sinking region stand for the ±20 degree band;
- both regions carry the band's cloud field, its export per unit area and its
  upwelling, so the rising region's extra high cloud is not represented;
- the cloud field is ISCCP's ±30 band, and the domain is ±20;
- the ocean transport is fixed: it does not respond to the surface
  temperature it creates;
- the land fractions, 0.25 and 0.21, are treated as ocean surface;
- the OAFlux period, 1991-2009, is shorter than the masks', 1991-2020.

## The runs

| run | what | steps |
|---|---|---:|
| N0 | **null.** One `tropical_column` at ±20° with the physics above, no ocean transport. Two identical columns with ΔQ = 0 are this column, by `spectral-wtg-circulation.md`'s Test A (bit-for-bit). A 2000-step two-column ΔQ = 0 run confirms that here | 100000 |
| W1 | **verdict.** Two columns, ΔQ = 13.39 (OAFlux), L = 734 km | 100000 |
| W2 | ΔQ = 4.67 (NCEP) | 100000 |
| W3 | W1 with L = 367 km | 100000 |
| W4 | W1 with L = 1468 km | 100000 |

**Equilibrium.** A run is equilibrated if the band-mean net top-of-atmosphere
flux over its last fifth is within ±0.1 W m^-2 of the band export. If not, it
is continued from its last state, in 100000-step increments, to at most
300000. Its means are always over the last fifth of the whole run.

## The tests

All judged on W1, over the last fifth.

**Decisive: does the overturning dry the free troposphere?**

| test | observed | passes if |
|---|---:|---|
| D1: sinking column RH 300-700 hPa | 0.35 | 0.25 to 0.45 |
| D2: band RH 300-700 hPa (area-weighted) | 0.43 | 0.33 to 0.53 |
| D3: band precipitable water 1000-300 hPa | 38.0 kg m^-2 | 28.5 to 47.5 (±25%) |

**The circulation's own predictions,** none used in its construction:

| test | observed | passes if |
|---|---:|---|
| C1: omega500, sinking column | +14.3 hPa/day | 7 to 29 |
| C2: sinking column's RH profile | C-shaped, minimum 0.28 at 400-500 hPa | the minimum lies between 350 and 650 hPa, below both the 850 hPa and 300 hPa values |
| C3: surface temperature, rising - sinking | +1.46 K | +0.73 to +2.92 K (half to double) |
| C4: T at 500 hPa, rising - sinking | +0.65 K | +0.33 to +1.30 K |

All quantities are computed on the model's layers exactly as on NCEP's:

- interpolated in log pressure to NCEP's levels;
- RH pressure-weighted over 300, 400, 500, 600 and 700 hPa;
- precipitable water integrated from the surface to 300 hPa.

**The verdict:**

- **The mechanism dries the free troposphere** if D1-D3 all pass.
- **Its physics is supported** if C1-C4 also pass.
- **It is refuted as the drying mechanism** if D1 fails, since it then cannot
  dry even the air it sinks.
- **Sensitivity.** Every test is reported on W2-W4. A verdict that changes
  across them is reported as fragile.

**Reported, not judged:**

- the band-mean surface temperature, against the observed 299.61 K, and N0's;
- each column's RH and temperature profiles;
- the cloud effects and cloud-top temperatures;
- the rising column's RH, against 0.50;
- the band's outgoing longwave and clear-sky outgoing longwave;
- the stratospheric water and the cold point.

## What is expected, stated before any run

**The rising column is not expected to dry.** Even a perfectly working
overturning dries only the sinking 44 percent of the band. If the rising
column keeps R4's 0.89, then D2's 0.53 ceiling needs the sinking column
below 0.03, which is not physical.

- **The observed rising region is 0.50,** because monthly-mean ascent at
  2.5° includes much air that is not convecting.
- **D2 is therefore expected to fail** unless the overturning's low-level
  inflow of dry sinking air also dries the rising column substantially.
- **D1 is the test of the mechanism itself.**

That is a prediction, not a criterion. The criteria above stand as written.

## Order and outcomes

**The order:**

1. the code and its unit tests;
2. the 2000-step symmetry check;
3. N0, then W1;
4. W2-W4, two at a time as memory allows.

**Outcomes:**

- **A failure of D1 on W1** is recorded as the mechanism's refutation in this
  form, and W2-W4 still run as frozen.
- **Nothing here changes a default.** Adoption is the user's decision.
- **The global column.** If the mechanism dries the tropics, a separate
  protocol is needed to carry it to the global column, whose extratropics
  this does not represent.

---
