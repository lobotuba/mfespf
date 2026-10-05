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

---

## Results (appended; nothing above this line changed)

MLS to 1 hPa, precipitable water 29.3 kg m^-2. The split of cut bins into
Chou's bands conserves the total to 0.1 W m^-2. Values are in W m^-2.

| band, cm^-1 | model up, top | line-by-line | D_b | model down, sfc | line-by-line |
|---|---:|---:|---:|---:|---:|
| 0-340 | 32.10 | 34.56 | -2.46 | 50.49 | 51.02 |
| 340-540 | 51.60 | 61.01 | **-9.41** | 80.91 | 80.76 |
| 540-800 | 68.31 | 67.90 | +0.41 | 101.04 | 105.85 |
| 800-980 | 55.91 | 58.62 | -2.71 | 33.06 | 26.89 |
| 980-1100 | 21.14 | 22.36 | -1.22 | 18.82 | 12.36 |
| 1100-1215 | 9.61 | 20.96 | **-11.35** | 22.35 | 9.73 |
| 1215-1380 | 8.85 | 15.89 | **-7.04** | 22.71 | 19.43 |
| 1380-1900 | 11.37 | 7.56 | +3.81 | 27.58 | 30.51 |
| 1900-3000 | 5.40 | 5.13 | +0.27 | 3.02 | 3.04 |
| **total** | **264.38** | **293.99** | **D = -29.70** | 359.98 | 339.59 |

**Verdict as frozen: 1100-1215 cm^-1 carries the deficit.** Its -11.35 is
more than a third of D (9.90). Two more are reported:

- 340-540 cm^-1, -9.41, falls just short of the threshold;
- 1215-1380 cm^-1, -7.04, is next.

**What is responsible.** The model grid has one bin across each of these
regions:

| band | the model bin that covers it | column optical depth | of which |
|---|---|---:|---|
| 1100-1215 and 1215-1380 | 1100-1778 cm^-1 | 21.3 + 0.55 | the vibrational (6.3 micron) wing, plus the flat continuum |
| 340-540 | 316-540 cm^-1 | 46.1 + 0.55 | the rotational wing, plus the continuum |

- **How the wing enters.** It is evaluated once, at each bin's centre, and
  applied across the whole bin. That makes opaque both the window edge at
  1100-1300 cm^-1 and the semi-transparent 400-540 cm^-1 region, where
  line-by-line lets 21 and 9 W m^-2 out.
- **It is the grid, not the wing's value.** The exponential wing varies by
  orders of magnitude inside one bin. `CLAUDE.md` makes the same point about
  CO2: "A grid must resolve what the physics depends on."
- **The continuum is small.** Its 0.55 is not what closes these bands.

**The second case, reported.** With the cloud-calibrated wing decay
(3411.86), the total is +1.05: OLR 295.13 against 293.99. Band by band it is
not:

| band | D_b |
|---|---:|
| 1380-1900 | +10.67 |
| 1100-1215 | -5.89 |
| 800-980 | -2.68 |
| 540-800 | +2.46 |
| 0-340 | -1.54 |
| 1215-1380 | -1.27 |

The right total comes from errors that cancel. That is what a fitted
coefficient standing in for an unresolved band does.

**What this points to.** The model's water vapour is a two-wing envelope with
one fitted decay, evaluated at bin centres on a coarse grid. The deficit is
that representation, not a missing gas. The candidate that would remove the
fit altogether is Chou et al.'s water-vapour k-distribution:

- bands 1-10, Table 4, with Table 3's pressure and temperature scaling;
- the self-broadened continuum, Table 9 and eq. 4.19-4.21, for bands 3-7;
- on Chou's own band edges.

That would replace the envelope and its wing decay. The cloud longwave effect
the wing has been standing in for would then have to come from clouds
themselves.
