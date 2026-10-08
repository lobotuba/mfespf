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

## `aer_rrtmg_lw/`

AER's own clear-sky example runs of RRTMG-LW for four standard atmospheres
(mid-latitude summer and winter, sub-arctic winter, tropical):
`input_rrtm_*-clr` and `output_rrtm_*-clr` from
`AER-RC/RRTMG_LW/run_examples_std_atm/` at commit
958edb1ad1618d76d0fcc5a9dc2e3c2f791a54dd, BSD 3-Clause (the licence is
reproduced in `src/mfespf/physics/data/rrtmg_lw_reduced.LICENSE`). Checksums
are in `data/sources.yaml` (`RRTMG-LW-AER-RUN-*`).

Used by `tests/unit/test_rrtmg_longwave.py` as the reference for the port
(`rrtmg-longwave.md` P0). The output is RRTMG itself, printed to 1e-4 W m^-2,
so the check is of transcription, not of physics.

## `afgl1986_us_standard.csv`

Table 1f (U.S. Standard) of the same report, from the same transcription
(`tables/1f.csv`, commit 9d170bb3e52b, SHA-256
d2d06777c7fbb9ed5bde5bf262a854ff61f9e83ea4d9b771006f10b87468eeda). Same columns.

The atmosphere on which the fixed Planck fractions of RRTMG's g-points are
averaged (`rrtmg-longwave.md`, erratum): neither benchmark atmosphere, and none
of the model's runs.
