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
