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

## `afgl1986_subarctic_winter.csv`

Table 1e (Subarctic Winter) of the same report, from the same transcription
(`tables/1e.csv`, commit 9d170bb3e52b, SHA-256
686a317f549f48cffadea61ea21057d4505c4d4b207cdd52cd73f8833d2ca400). Same columns.

Used as the second atmosphere of the Chou et al. (2001) longwave benchmark
(their Table 15). It is cold and dry where mid-latitude summer is warm and moist,
so a water-vapour scheme that is right on one only by compensation shows it on
the other.
