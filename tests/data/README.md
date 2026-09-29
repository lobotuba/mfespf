# Test data

## `afgl1986_midlatitude_summer.csv`

Table 1b (Midlatitude Summer) of Anderson, G. P., Clough, S. A., Kneizys, F. X.,
Chetwynd, J. H. and Shettle, E. P. (1986), *AFGL Atmospheric Constituent
Profiles (0-120 km)*, AFGL-TR-86-0110, Environmental Research Papers No. 954
(public report: https://apps.dtic.mil/sti/tr/pdf/ADA175173.pdf).

Electronic transcription from https://github.com/rayference/afgl1986
(`tables/1b.csv`), MIT License, Copyright (c) 2021 Yvan Nollet. Columns:
altitude `z` (km), pressure `p` (mb), temperature `t` (K), number density `n`
(cm^-3), and volume mixing ratios of H2O, O3, N2O, CO and CH4 (ppmv).

Used as the atmosphere of the Chou and Suarez (1999) Table 8 benchmark in
`tests/unit/test_solar_vapour.py`: it is the profile their line-by-line
comparison was run on, so the comparison is like for like.
