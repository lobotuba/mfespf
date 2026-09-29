# Protocol: the overturning circulation, relaxed and spectral

**Frozen before the relaxed circulation is written and before any run.** Results
are appended below a line at the end; nothing above it changes afterwards.

## Why

`two-column-circulation.md` froze the strict WTG balance: subsidence warming
equals radiative cooling layer by layer, and the mass flux is the cooling
divided by the rising column's dry-static-energy step. Its Test A could not be
completed. Near the cold top of convection the entraining-plume profile is
almost dry adiabatic, so the step goes to zero and the flux has no bound. It
reached about 1100 hPa/day at 110 hPa and emptied a layer in one step. Published
WTG schemes all avoid this with a floor on the static stability (2 K/km). That
is a numerical guard, not physics, and it is not used here.

## The physics

The weak temperature gradient is not imposed; it is the result of gravity waves
removing horizontal temperature differences. Each column relaxes toward the
area mean at the rate gravity waves cross the distance between them. Each
vertical mode relaxes at its own rate: mode *n* travels at c_1/n, so it relaxes
n times more slowly (Herman and Raymond 2014, doi:10.1002/2014MS000359; in the
form of Wang et al. 2016 as implemented in SCAM6, Cohen et al. 2024,
doi:10.1029/2023MS003866, their section 2.1):

- **Forcing.** The sinking column's temperature departure from the area mean is
  T_d - T_bar = sigma (T_d - T_a), with sigma = 0.476 the rising fraction.
  - The pressure velocity that would remove it on the time tau_1 is
    omega = c_p (T_d - T_bar) / (tau_1 dS/dp).
  - The downward mass flux is f = omega/g.
  - c_p is each layer's own heat capacity per unit mass.
- **Uniform stratification.** The modes below are exact only if the buoyancy
  frequency is uniform. So dS/dp is the rising column's *bulk* value between its
  lowest layer and its highest layer below the tropopause, not a layer-by-layer
  step. The same bulk value gives N^2 = g (dS/dz) / (c_p T_bar_trop).
  Nothing is divided by a local step, so nothing needs a floor.
- **Modes.** Normalized pressure is zeta = (p_s - p)/(p_s - p_t), with the
  tropopause p_t = 100 hPa, the SCAM6 value.
  - A_n = 2 * integral over zeta from 0 to 1 of f(zeta) sin(n pi zeta),
    evaluated with each layer's pressure thickness below p_t as the weight.
  - The flux at each interface is F = sum over n of (A_n / n) sin(n pi zeta),
    for n from 1 to M.
  - M is the number of layers at least partly below p_t, at most 20 (SCAM6).
  - F is zero at the surface and at the tropopause, and above it. No
    boundary-layer taper is needed (Herman and Raymond 2014).
- **Relaxation time.** tau_1 = L / c_1, with c_1 = N H / pi and H the
  tropopause height of the rising column.
  - L = 750 km, measured here and fitted to nothing. It is the area-weighted
    mean great-circle distance from each descending grid point (omega500 > 0)
    within 20 to 30 degrees of the equator to the nearest ascending one.
  - Source: NCEP/NCAR R1 1991-2020 monthly long-term means, the same file that
    gives sigma. Means are 734 km (20 degrees) and 770 km (30 degrees); monthly
    values span 589 to 853.
- **Everything else is unchanged from the strict protocol.**
  - The rising column carries the mirror mass flux scaled by the area ratio.
  - Transport is flux-form and upwind, carrying dry static energy, vapour and
    latent energy.
  - The circulation runs after radiation and before water and convection,
    which compose with its tendencies.
  - Between columns, amounts are recorded per unit of global area.

The 850 hPa boundary-layer top of the strict protocol is no longer a parameter.
Two quantities are declared, and both are observed: sigma and L.

## Declared approximations

- The modes are written in pressure rather than height (a Boussinesq-like
  approximation, as in the sources).
- The stratification is uniform, by construction of the modes.
- One rising and one sinking region stand for the whole globe.
- Midlatitude eddies are not represented.
- Both columns share every other setting.

## The tests

**Test A: the circulation's own predictions.** Two `reference_column`s,
30 layers, 3600 s steps. They are run to equilibrium: net top-of-atmosphere
flux within +/-0.5 W m^-2 over the last fifth of the run. None of these
quantities is used in the construction:

| quantity | observed | passes if |
|---|---|---|
| subsidence at 500 hPa, sinking column | 10-20 hPa/day, the peak of the tropical distribution (Bony et al. 2004) | within 5-40 hPa/day |
| free-tropospheric RH 300-700 hPa, area-weighted | 0.4 to 0.5 | within 0.35-0.55 |
| sinking column's humidity profile | C-shaped, minimum near 500 hPa (Couhert et al. 2010) | minimum between 400 and 700 hPa, below both ends of that range |
| **new:** rising minus sinking temperature at 500 hPa | +0.71 K within 20 degrees, +1.76 K within 30 degrees; monthly +0.51 to +2.14 (NCEP R1 1991-2020, omega500 composite, computed here) | within +0.5 to +2.2 K |

The new row is the one prediction that depends on L. The verdict is taken at
L = 750 km. The run is repeated at 375 and 1500 km, and any verdict that changes
across that range is reported as fragile.

**Test B: the hypothesis the cloud test left.** Criteria are exactly those of
`two-column-circulation.md`, run on this circulation:

- **Supported** if the area-weighted longwave cloud radiative effect rises by at
  least 8 W m^-2 above the single column's +11.6, **and** the air at ISCCP's
  cloud-top pressures moves toward ISCCP's cloud-top temperatures at every
  level.
- **Refuted** if the free troposphere dries into the observed range and neither
  happens.
- **Adoption as the reference configuration** needs all of the following,
  unchanged:

  | quantity | required |
  |---|---|
  | planetary albedo | 0.291 +/- 0.03 |
  | SW CRE | -45.8 +/- 12 W m^-2 |
  | LW CRE | +27.9 +/- 8 W m^-2 |
  | surface net radiation | at least 40 W m^-2 below 183.5 |
  | evaporation | below 4.11 mm/day |

The observed surface temperature decides nothing here; it is reported.

---

## Results (appended; nothing above this line changed)

### Test A: fails. Two identical columns never make a temperature difference

Two 30-layer `reference_column`s, 3600 s steps, 44000 steps. Last fifth of each
run:

| run | tau_1 | net TOA | omega500 | RH 300-700 | T_a - T_d at 495 hPa | surface |
|---|---:|---:|---:|---:|---:|---:|
| L = 375 km | 1.45 h | -0.027 | 0.0 | 0.741 | 0.000 K | 291.95 K |
| **L = 750 km (verdict)** | 2.90 h | -0.027 | 0.0 | 0.741 | 0.000 K | 291.95 K |
| L = 1500 km | 5.80 h | -0.027 | 0.0 | 0.741 | 0.000 K | 291.95 K |

Energy and water closed to 9.2e-16. The two columns are bit-for-bit identical in
all three runs, and the verdict does not change across the sensitivity range.

| criterion | model | verdict |
|---|---:|---|
| subsidence at 500 hPa, 5-40 hPa/day | 0.0 | fail |
| RH 300-700 hPa, 0.35-0.55 | 0.741 | fail |
| C-shaped profile, minimum between 400 and 700 hPa | minimum 0.581 at 294 hPa | fail |
| T_a - T_d at 500 hPa, +0.5 to +2.2 K | 0.000 | fail |

**Why.** Same surface, same sunlight, same start. The relaxed balance moves air
only in response to a temperature difference, and nothing makes one. The
symmetric state is an exact solution. The strict balance circulated only
because it declared which column sinks. The relaxed-WTG literature this
follows always prescribes a contrast, usually in sea-surface temperature
(Daleu et al. 2016). This framework has none.

The tau_1 the measured L gives, 2.9 h, is close to the 3 h of the GASS
intercomparison, but nothing was fitted to it.

Test B was not run. Its premise is a free troposphere dried by sinking air, and
there is none.

### Diagnostics after the verdict (not part of it)

**Sunlight does not choose the rising region.** Insolation on the NCEP omega500
masks (daily-mean top-of-atmosphere insolation at S0 = 1361 W m^-2, averaged
monthly):

- rising regions 338.9, sinking 341.7 W m^-2, difference -2.75;
- monthly differences range from -10.7 to +15.1.

The global omega500 < 0 area includes the midlatitude storm tracks, and the
omega500 > 0 area includes the sunny subtropics. A global partition mixes the
tropical overturning with midlatitude weather, which is not a WTG regime. The
declared approximations already say so, and this measures it.

**The equilibrium depends on water the initial state puts above convection.**

1. The two columns were rerun with the sinking column's atmosphere started
   0.1 K cooler. The circulation decayed (omega500 0.14 hPa/day, T_a - T_d
   -0.001 K). While it ran, it lifted water above the level convection had
   reached.
2. Both columns then settled into a different state:
   - surface 296.59 K, 4.6 K warmer than the symmetric state;
   - precipitable water 61 against 37 kg m^-2;
   - RH about 0.96 from 160 to 94 hPa, where the symmetric state has exactly
     zero;
   - a cold point of 196.9 K at 94 hPa.
3. To locate the cause, the single `reference_column` was started four ways:

   | start | surface | moist (RH > 0.05) up to | above that |
   |---|---:|---:|---|
   | as is (dry) | 291.95 K | 261 hPa | RH 0.000 |
   | +0.1 K | 291.97 K | 261 hPa | RH 0.000 |
   | -0.1 K | 291.94 K | 261 hPa | RH 0.000 |
   | 1e-6 kg/kg vapour in every layer | 293.30 K | 227 hPa | RH 0.008-0.013 |

   All four runs are equilibrated (net within +/-0.03 W m^-2) and closed to
   about 1e-15.

**So the column has no single equilibrium: the state depends on how much
water starts above the convective top.** The column starts dry by default.
Nothing in it moves water above the level convection reaches, and nothing
there removes water except condensation. So whatever vapour is there at the
start stays, and it sets the longwave cooling that decides how deep convection
goes. A trace of vapour, 1 ppm by mass and below the real lower stratosphere's
roughly 3, moves the surface 1.35 K. The circulation's brief moistening moved
it 4.6 K.

The reference column's headline figure, +3.94 K, is the dry-start member of
this family. It is not a unique prediction until a physical process sets the
water above the convective top. In the real atmosphere that is transport through
the cold point: freeze-drying of air entering the stratosphere, and large-scale
ascent through the tropical tropopause layer.
