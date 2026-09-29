# Protocol: air that rises in one column and sinks in another

**Frozen before the circulation module is written and before any run.** Results
are appended below a line at the end; nothing above it changes afterwards.

## Why

The observed-cloud test failed on its longwave half because the model's air at
ISCCP's cloud-top pressures was 13 to 28 K warmer than ISCCP's cloud tops, in a
free troposphere at 0.775 relative humidity against an observed 0.4 to 0.5
(`observed-cloud-field.md`). The leading hypothesis was that the wetness causes
the warmth. A single horizontally uniform column cannot test it: its mass-
weighted vertical motion is zero, so it cannot hold the sinking air that dries
Earth's free troposphere. Two columns can.

## The physics

Two regions, each a full `reference_column`, linked by an overturning
circulation under the weak-temperature-gradient balance (Sobel and Bretherton
2000; the two-column framing of Pierrehumbert 1995 and Nilsson and Emanuel 1999):

- **Areas from observation.** The rising region covers 0.476 of the globe: the
  area-weighted fraction with upward monthly-mean vertical velocity at 500 hPa
  in the NCEP/NCAR Reanalysis 1 1991-2020 long-term means (computed here from
  `omega.mon.ltm.1991-2020.nc`, NOAA PSL; 0.461 to 0.490 by month). The sinking
  region covers the rest.
- **Subsidence from the WTG balance.** In the sinking column's free troposphere,
  subsidence warming balances radiative cooling: the downward mass flux entering
  each layer from above is its radiative cooling divided by the dry-static-energy
  step to the layer above. Zero where a layer is radiatively heated or not
  statically stable, and zero above the highest cooled layer. Below 850 hPa (a
  declared boundary-layer top) the flux falls linearly to zero at the surface,
  and the air leaves laterally for the rising column. The rising column carries
  the same mass flux upward, scaled by the area ratio.
- **Transport is flux-form and upwind.** Dry static energy and water vapour are
  carried by the vertical fluxes from the upstream layer and by the lateral
  exchange from the upstream column. Every joule and kilogram that leaves one
  place arrives in another, exactly; between the columns it is recorded as a
  transfer per unit global area, each budget taking its own area share. Vapour
  carries its latent energy with it.
- **No new parameter is fitted.** Two are declared: the boundary-layer top
  (850 hPa) and the area fraction (observed above).

## Declared approximations

- One rising and one sinking region stand for the whole globe; midlatitude
  eddies, which are not a WTG regime, are not represented.
- Radiative heating alone balances subsidence; latent heating in the sinking
  column's free troposphere is taken as small.
- Both columns share every other setting, including, in the second test below,
  the global-mean observed cloud field -- real clouds differ by regime.

## The tests

**Test A -- the circulation's own predictions**, on two `reference_column`s,
30 layers, equilibrated. None is used anywhere in its construction:

| quantity | observed | passes if |
|---|---|---|
| subsidence at 500 hPa in the sinking column | 10-20 hPa/day, the peak of the tropical distribution (Bony et al. 2004) | within 5-40 hPa/day |
| free-tropospheric relative humidity, 300-700 hPa, area-weighted | 0.4 to 0.5 | within 0.35-0.55 |
| sinking column's humidity profile | C-shaped, minimum near 500 hPa (Couhert et al. 2010) | minimum between 400 and 700 hPa, below both ends of that range |

**Test B -- the hypothesis the cloud test left**, the same two columns with the
observed cloud field and water vapour's solar absorption in both, against the
single-column run of `observed-cloud-field.md` and CERES EBAF Ed4.0:

- **The hypothesis is supported if** the area-weighted longwave cloud radiative
  effect rises by at least 8 W m^-2 above the single column's +11.6 **and** the
  model's air at ISCCP's cloud-top pressures moves toward ISCCP's cloud-top
  temperatures in every level. It is refuted if the free troposphere dries into
  the observed range and neither happens.
- **Adoption as the reference configuration** needs everything
  `observed-cloud-field.md` required: planetary albedo 0.291 +/- 0.03, SW CRE
  -45.8 +/- 12, LW CRE +27.9 +/- 8, surface net radiation at least 40 W m^-2 below
  183.5, evaporation below 4.11 mm/day.

The observed surface temperature decides nothing here; it is reported.
