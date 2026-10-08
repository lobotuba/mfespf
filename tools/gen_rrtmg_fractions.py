"""Add RRTMG's fixed Planck fractions to the reduced tables (rrtmg-longwave.md, erratum).

RRTMG's own per-layer Planck fractions on the AFGL US Standard Atmosphere (AFGL
1986 Table 1f, tests/data/afgl1986_us_standard.csv), averaged over its layers by
layer mass (the pressure difference) and renormalised within each band.

Layers lie between consecutive tabulated levels, at the mean of their levels'
temperatures and mixing ratios. H2O, O3, N2O, CO and CH4 are as tabulated; CO2
is the model's 284.7 ppmv (the table has none) and O2 0.209.

Usage:
    python tools/gen_rrtmg_fractions.py tests/data/afgl1986_us_standard.csv
        src/mfespf/physics/data/rrtmg_lw_reduced.npz
"""

import csv
import sys

import numpy as np

from mfespf.physics import rrtmg_longwave as rr

rows = list(csv.DictReader(open(sys.argv[1])))
level = {
    k: np.array([float(r[k]) for r in rows]) for k in ("p", "t", "H2O", "O3", "N2O", "CO", "CH4")
}
mid = {k: 0.5 * (v[:-1] + v[1:]) for k, v in level.items()}
vmr = {k.lower(): mid[k] * 1.0e-6 for k in ("H2O", "O3", "N2O", "CO", "CH4")}
vmr.update(co2=284.7e-6, o2=0.209)
coldry = rr.dry_air_column(level["p"][:-1], level["p"][1:], vmr["h2o"])
col = rr.RRTMGColumn.from_mixing_ratios(
    pressure_hpa=mid["p"], temperature=mid["t"], coldry=coldry, vmr=vmr
)
_, fracs, laytrop = rr.rrtmg_gas_optics(col)
weight = level["p"][:-1] - level["p"][1:]
mean = np.sum(fracs * weight[:, None], axis=0) / np.sum(weight)
ngs = np.concatenate(([0], np.cumsum(rr.G_POINTS_PER_BAND)))
fixed = np.empty_like(mean)
for b in range(16):
    s = slice(ngs[b], ngs[b + 1])
    fixed[s] = mean[s] / mean[s].sum()
with np.load(sys.argv[2]) as data:
    tables = {k: data[k] for k in data.files if k != "fixed_fractions"}
np.savez_compressed(sys.argv[2], **tables, fixed_fractions=fixed)
print(f"{len(rows) - 1} layers, laytrop {laytrop}; band sums before renormalising:")
print(np.round([mean[ngs[b] : ngs[b + 1]].sum() for b in range(16)], 4))
