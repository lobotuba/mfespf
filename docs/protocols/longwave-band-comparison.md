# Protocol: where the model's longwave differs from line-by-line, band by band

**Frozen before the comparison is run.** Results are appended below a line at
the end; nothing above it changes afterwards. This is a diagnostic. It changes
no module.

## Why

On the AFGL mid-latitude summer (MLS) atmosphere, the model's outgoing
longwave is 264.4 W m^-2 against line-by-line 293.99 (`co2-band.md`,
`ozone-infrared-band.md`). Two of its bands have been checked and agree:

- 540-800 cm^-1 within 0.4 W m^-2;
- 980-1100 cm^-1 within 1.2, on the corrected reading.

So the 30 W m^-2 is elsewhere. Two things depend on finding it:

- the vapour's fitted wing decay, which currently absorbs whatever is
  missing;
- the forcing and feedback that come out of it.

## The comparison

Chou et al. (2001) Table 14 gives, for each of nine bands, line-by-line values
of:

- the downward flux at the surface, F_down,sfc;
- the upward flux at the top, F_up,top.

The bands are 0-340, 340-540, 540-800, 800-980, 980-1100, 1100-1215,
1215-1380, 1380-1900 and 1900-3000 cm^-1. Its absorbers are water vapour
lines and continuum, CO2 and ozone. The minor gases are left out.

**The model**, set up exactly as for C1 and Z1:

- MLS from the repository's AFGL table, down to 1 hPa;
- CO2 at 350 ppmv;
- Chou's 4e-6 g/g of water above the tropopause;
- the 16-band grid rebuilt for the CO2 band and the ozone band;
- the vapour wing decay C1 used, 4789.08 m^-1.

Each model band's flux is summed into Chou's band edges.

- **Cut bins.** Where a model bin crosses one of Chou's edges, its flux is
  split by Planck weight at the bin's emitting temperature. That is declared;
  the model grid does not share Chou's edges.
- **Second case, reported.** The same with the cloud-calibrated wing decay,
  3411.86 m^-1.

## The finding

- **The deficit.** D_b = model minus line-by-line F_up,top for each band b,
  and D = the sum over all nine.
- **What carries it.** A band carries the deficit if D_b has D's sign and is
  at least a third of D. Several can.
- **Which piece of the vapour is responsible.** For each band that carries it,
  the model's optical depth there is decomposed into the vapour's rotational
  wing, its vibrational wing, its flat continuum, CO2 and ozone. The column
  totals of each are reported.
- **Reported:** F_down,sfc per band against line-by-line, and the total.

**What it is for.** The band or bands that carry the deficit, and the vapour
term behind them, are what the next mechanism replaces. The candidates are:

- Chou's water-vapour k-distribution for those bands (Table 4, with Table 3's
  scaling);
- and/or his self-broadened continuum (Table 9 and eq. 4.19-4.21).

Each would have its own protocol.
