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

---

## Result (appended after the runs; nothing above this line was changed)

`reference_column`, 30 layers, equilibrated at 28000 steps, budgets closing to
1e-15. Raw rows: `observed-cloud-field.results.jsonl`.

| quantity | clouds only | clouds + solar vapour | observed | test |
|---|---:|---:|---:|---|
| planetary albedo | 0.338 | **0.283** | 0.291 (CERES) | +/- 0.03: **pass** |
| shortwave cloud effect (W m^-2) | -58.1 | **-52.2** | -45.8 | +/- 12: **pass** |
| longwave cloud effect (W m^-2) | +17.9 | **+11.6** | +27.9 | +/- 8: **fail** |
| clear-sky reflected (W m^-2) | 56.8 | 44.2 | 53.3 | |
| clear-sky outgoing longwave (W m^-2) | 243.3 | 255.4 | 268.1 | |
| sunlight absorbed in the atmosphere | 0.0 | 59.1 | ~78 | |
| surface net radiation (W m^-2) | 188.4 | **159.6** (-23.9) | ~104 | fall >= 40: **fail** |
| evaporation (mm/day) | 4.58 | **4.51** (was 4.11) | 2.7 | must fall: **fail** |
| precipitable water (kg m^-2) | 50.7 | 86.6 | 25 | |
| surface-temperature error (K) | +6.80 | +12.23 | 0 | |

**Not adopted.** The shortwave half passes: with the reflection moved from the
surface to cloud tops, the planetary albedo is 0.283 against 0.291 and the
shortwave cloud effect -52 against -46, nothing fitted -- the reflector is now in
the right place, which is what this step was for. The longwave cloud effect
fails, at +11.6 against +27.9, and the surface test fails: net radiation fell 24
W m^-2, not 40, and evaporation rose.

**Why the longwave effect fails, measured rather than assumed.** At the pressures
where ISCCP observes the cloud tops, the model's atmosphere is far warmer:

| cloud | pressure | model temperature | ISCCP top temperature | difference |
|---|---:|---:|---:|---:|
| cirrus and cirrostratus | 261 hPa | 252.3 K | 227.5 K | +24.8 K |
| deep convective | 328 hPa | 262.6 K | 234.8 K | +27.8 K |
| middle | 562 hPa | 281.3 K | 262.8 K | +18.5 K |
| low | 829 hPa | 293.8 K | 281.1 K | +12.7 K |

A cloud top that warm emits nearly what the clear sky does, so it traps little.
The cloud optics are not the failure; the atmosphere the clouds sit in is, and
more so with height. Its free troposphere carries 0.775 relative humidity from
300 to 700 hPa against an observed 0.4 to 0.5. The warm cloud tops are measured;
that the over-moist free troposphere causes them is the leading hypothesis, not
yet a measurement -- a wet column's strong vapour greenhouse and feedback would
produce exactly this warm, weakly contrasting profile. The next mechanism, sinking
air that dries the free troposphere, is what tests it: if it is right, the cloud
tops cool toward ISCCP's and the longwave effect rises with nothing about the
clouds changed.

*Gate note.* Both runs were first marked invalid by the validity-domain gate: the
cloud record's optical-thickness range had been inferred onto the radiation
module's *total* layer optical depth, saturated vapour bands at 1600. That was a
declaration error in this step, fixed by a distinct symbol, `tau_c`, bound
explicitly to the cloud field's own output. The binding changes which field the
gate judges and no computed number; the rerun is valid and the table stands.

