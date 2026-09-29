# Protocol: water vapour absorbs sunlight

**Frozen before the implementation and before any run.** Results are appended
below a line at the end; nothing above it changes afterwards.

## The missing mechanism

`reference_column` nets 183.5 W m^-2 of radiation at the surface against an
observed ~104 (Trenberth, Fasullo and Kiehl 2009), and sheds the excess as
evaporation (119 W m^-2 against 80) and sensible heat (55 against 17). Earth's
atmosphere absorbs about 78 W m^-2 of sunlight; this one absorbs none, because
its vapour bands stop at 4 microns. Most of Earth's atmospheric solar
absorption is by water vapour in the near infrared.

## The physics, and where its numbers come from

Chou and Suarez (1999), *A Solar Radiation Parameterization for Atmospheric
Studies*, NASA/TM-1999-104606 vol. 15 -- public domain, fitted by its authors
to line-by-line calculations, never to an energy budget:

- Table 2: ten-term k-distributions for water vapour in 0.70-1.22, 1.22-2.27 and
  2.27-10 microns, with flux weights `h` (fractions of the total solar flux,
  summing to 0.52926) and absorption coefficients `k` at 300 hPa and 240 K.
- Table 3: `k = 0.00075 g^-1 cm^2` for water vapour in 0.4-0.7 microns.
- Scaling (their 3.2-3.5): the absorber amount is weighted by
  `(p / 300 hPa)^0.8 [1 + 0.00135 (T - 240 K)]`.

Within each of their bands the grid carries one sub-band per k-term, with edges
placed so that each sub-band receives its term's share of the band's incident
sunlight. For everything that does not vary with wavelength inside the band --
the surface albedo, and a thermal emission that is negligible below 4 microns --
that is exactly the k-distribution. Nothing is fitted here.

Ozone is not part of this step: it needs an observed ozone profile, recorded as
a dataset with its own provenance, and is the next one.

## Declared in advance

- **Correctness**, checked against the source and not against the model: a
  uniform path through the sub-bands transmits `sum h exp(-k w) / sum h`, the
  source's equation 3.4; energy closes as before.
- **Two approximations that will show:** sunlight is solved as diffuse
  radiation with the solver's diffusivity factor rather than as a beam at a mean
  zenith angle; and the 0.30 albedo sits at the surface, standing in for clouds,
  so reflected sunlight crosses the absorbing vapour twice where on Earth most is
  reflected above most of it. The second lowers the planetary albedo and warms
  the column. It is measured and reported, not compensated.

## The falsification test

On `reference_column` at 30 layers, equilibrated:

- **The mechanism is falsified as the explanation of the surface excess if the
  surface net radiation does not fall by at least 40 W m^-2** -- half the gap to
  the observed ~104 -- or if evaporation does not fall.
- Reported against the same budget, none of it used in any calibration here:
  atmospheric solar absorption (obs ~78, of which clouds, ozone and aerosol are
  not in this model), surface net radiation, latent and sensible heat, planetary
  albedo, and the surface-temperature error.

The observed surface temperature plays no part in the decision. The physics is
added because it is physics with published coefficients; the test decides
whether it is the explanation claimed.

---

## Result (appended after the runs; nothing above this line was changed)

`reference_column`, 30 layers, 16-band base grid rebuilt to 43 bands, time-means
over the last fifth, both equilibrated (28000 steps without; 56000 with, under
0.03 W m^-2 -- 100000 was refused by the memory preflight at 3.7 GB against a
2 GB budget, and the budget was not raised). Budgets close to 1e-15.

| quantity | without | with vapour solar absorption | observed |
|---|---:|---:|---:|
| sunlight absorbed in the atmosphere (W m^-2) | 0 | **71.1** | ~78 |
| sunlight absorbed at the surface (W m^-2) | 233.7 | 186.9 | ~161 |
| planetary albedo | 0.295 | **0.223** | ~0.30 |
| surface net radiation (W m^-2) | 183.6 | 172.3 | ~104 |
| latent heat (W m^-2) | 118.9 | 147.6 | ~80 |
| sensible heat (W m^-2) | 55.0 | 15.5 | ~17 |
| precipitable water (kg m^-2) | 37.4 | 124.2 | ~25 |
| outgoing longwave (W m^-2) | 240.0 | 264.3 | ~239 |
| surface-temperature error (K) | +3.94 | **+16.31** | 0 |

**Falsified as the explanation of the surface excess**: surface net radiation
fell by 11 W m^-2 against the 40 required, and evaporation rose.

What the runs show, in the order the protocol declared them:

1. **The absorption itself is right.** Chou and Suarez's coefficients, untouched,
   put 71 W m^-2 of sunlight into the atmosphere against an observed ~78 that
   also includes ozone, cloud and aerosol absorption this model does not have.
   Nothing was fitted to it.
2. **The declared confound decides the outcome.** With the 0.30 albedo on the
   surface, reflected sunlight crosses the absorbing vapour twice. The planetary
   albedo fell from 0.295 to 0.223 -- about 25 W m^-2 more sunlight absorbed --
   and the column warmed until the vapour feedback carried it to 124 kg m^-2 of
   water and +16.3 K. On Earth most sunlight is reflected by clouds above most of
   the vapour, which is exactly what a surface albedo cannot do.
3. So the surface excess is not a missing absorber alone. It is a missing
   absorber *and* a reflector in the wrong place, and the second cannot be fixed
   by any setting of the first. The absorption is kept, off by default
   (`solar_vapour_absorption`), until the reflection is where Earth's is: clouds
   at altitude, with the surface's own albedo -- the next mechanism.

### Correction to point 1 (appended later)

"The absorption itself is right" was overstated, and was questioned. 71 against
an observed ~78 is not a valid check: the 78 includes ozone, oxygen, carbon
dioxide, cloud and aerosol absorption this model does not have, so vapour alone
*should* come in below it, and the 71 came from a column carrying five times
Earth's water. It could neither confirm nor refute the absorption.

The valid check is like for like: Chou and Suarez's own Table 8, the AFGL
mid-latitude summer atmosphere (29.3 kg m^-2 of water), sun at 60 degrees,
surface albedo 0.2, no scattering, against line-by-line calculations. As a
fraction of each band's sunlight absorbed in the atmosphere:

| band | model, vapour only, 60-degree path | model, diffuse (1.66) | line-by-line, all gases |
|---|---:|---:|---:|
| 0.70-1.22 microns | 0.221 | 0.207 | 0.224 |
| 1.22-2.27 microns | 0.410 | 0.399 | 0.447 |

Within 1 to 8 percent in each band, every difference in the direction and of the
size of the oxygen and carbon dioxide the model omits (about 4.3 and 3.3
W m^-2, their Table 7). The diffuse treatment of sunlight costs a further few
percent. That is now `test_the_mls_atmosphere_absorbs_what_line_by_line_says`.
The falsification verdict above does not change: it was decided by the surface
net radiation, not by this comparison.

