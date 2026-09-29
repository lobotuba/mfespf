# Protocol: clouds that reflect from altitude

**Frozen before the implementation and before any run.** Results are appended
below a line at the end; nothing above it changes afterwards.

## Why

The column's 0.30 albedo stands on its surface, doing the clouds' job at the
ground. Adding water vapour's solar absorption exposed it: reflected sunlight
crossed the absorbing vapour twice, the planetary albedo fell to 0.22 and the
column warmed to +16.3 K (`solar-vapour-absorption.md`). On Earth most sunlight
is reflected at cloud tops above most of the vapour. The mechanism to add is the
clouds, at their observed heights, over the surface's own reflectance.

## The cloud field, entirely from one observational source

Rossow and Schiffer (1999), *Advances in Understanding Clouds from ISCCP*,
BAMS 80, 2261-2287, Table 4 (D-series global annual means) and section 3:

| type | amount | top pressure | visible tau | phase (Table 5 split) |
|---|---:|---:|---:|---|
| cirrus + cirrostratus | 0.196 | 267 hPa | 2.2 | ice |
| deep convective | 0.026 | 326 hPa | 35.6 | ice |
| middle | 0.190 | 557 hPa | 4.8 | 52 % ice, 48 % liquid |
| low | 0.275 | 826 hPa | 4.7 | 9 % ice, 91 % liquid |
| surface reflectance | | | 0.17 | |

- The amounts are as seen from space, so they partition the sky: each type is
  one sub-column, and the remaining 0.313 is clear. The column's radiation is the
  amount-weighted mean of five independent solves.
- Each cloud occupies the one model layer containing its top pressure.
- Sunlight: extinction equal to the visible tau, conservative scattering,
  asymmetry 0.85. ISCCP's tau is a radiatively weighted (essentially
  logarithmic) average, chosen by its authors to represent the area-averaged
  radiative effect, which is what this needs.
- Infrared: absorption optical depth tau / 2.56 for liquid and tau / 2.13 for
  ice (ISCCP D-series), phase-weighted where a level holds both; no infrared
  scattering, as ISCCP assumes.
- The surface reflectance 0.17 replaces the 0.30 in every solar band.

Nothing is fitted. The field carries no water: it is optical properties at
observed places, and it enters no budget.

## Declared approximations

- No cloud below a visible cloud (ISCCP cannot see it); no overlap between
  types. Multi-layer cloud is real and is left out.
- One layer per cloud type; ice and liquid share their level's tau.
- ISCCP's surface reflectance is a visible-wavelength value used for the whole
  solar spectrum.
- Sunlight is still solved as diffuse radiation.

## The test

On `reference_column` at 30 layers, equilibrated, with the observed cloud field
**and** water vapour's solar absorption, against CERES EBAF Ed4.0 (Loeb et al.
2018, J. Climate 31, Table 5, July 2005 to June 2015) -- none of which any part
of this configuration has seen:

| quantity | observed | passes if within |
|---|---:|---:|
| planetary albedo | 0.291 | +/- 0.03 |
| shortwave cloud radiative effect | -45.8 W m^-2 | +/- 12 |
| longwave cloud radiative effect | +27.9 W m^-2 | +/- 8 |

and the original question, from `solar-vapour-absorption.md`: **surface net
radiation must fall by at least 40 W m^-2 from 183.5, and evaporation must
fall**.

The tolerances allow for a single global-mean column with four cloud types,
logarithmic-mean optical thicknesses and no multi-layer cloud; they were set
here, before the run.

**Adoption:** if every test passes, the observed cloud field, the 0.17 surface
reflectance and the solar vapour absorption become `reference_column`'s
defaults together. If any fails, none is adopted and the finding is recorded.
Also reported, not decided on: the surface-temperature error, the atmospheric
and surface solar absorption, the Bowen ratio and the precipitable water; and
the same configuration without solar vapour absorption, to separate the two
mechanisms.
