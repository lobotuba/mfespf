"""Build the reduced (140 g-point) RRTMG-LW tables from AER's source (rrtmg-longwave.md).

Source: AER-RC/RRTMG_LW commit 958edb1ad1618d76d0fcc5a9dc2e3c2f791a54dd.
- k-tables, Planck fractions: data/rrtmg_lw.nc, g-point set 1 (16 per band),
  read exactly as src/rrtmg_lw_read_nc.f90 reads them;
- reduction to 140 g-points: lwcmbdat and cmbgb1-16 of src/rrtmg_lw_init.f90
  (k weighted by rwgt within each combined group, fractions summed);
- reference atmosphere (pref, tref, chi_mls) and integrated Planck functions
  (totplnk, totplk16): the netCDF's copies of lwatmref and lwavplank.

Usage:
    python tools/gen_rrtmg_tables.py data/raw/rrtmg_lw.nc OUT
with OUT src/mfespf/physics/data/rrtmg_lw_reduced.npz.
Output: one .npz, every array in the Fortran index order of the reduced
module arrays, flattened where taumol uses flat (absa/absb) indexing.
"""

import sys

import numpy as np
import scipy.io as sio

nc = sio.netcdf_file(sys.argv[1], "r", mmap=False)
V = {k: np.array(v[:]) for k, v in nc.variables.items()}
names = [b"".join(r).decode().strip() for r in nc.variables["AbsorberNames"][:]]

# ---- lwcmbdat ------------------------------------------------------------------
NGC = [10, 12, 16, 14, 16, 8, 12, 8, 12, 6, 8, 8, 4, 2, 2, 2]
NGN = [1,1,2,2,2,2,2,2,1,1, 1,1,1,1,1,1,1,1,2,2,2,2, *[1]*16, *[1]*13, 3, *[1]*16, *[2]*8,
       2,2,1,1,1,1,1,1,1,1,2,2, *[2]*8, 1,1,1,1,1,1,1,1,2,2,2,2, 2,2,2,2,4,4,
       1,1,2,2,2,2,3,3, 1,1,1,1,2,2,4,4, 3,3,4,6, 8,8, 8,8, 4,12]  # fmt: skip
WT = np.array([0.1527534276, 0.1491729617, 0.1420961469, 0.1316886544, 0.1181945205,
               0.1019300893, 0.0832767040, 0.0626720116, 0.0424925000, 0.0046269894,
               0.0038279891, 0.0030260086, 0.0022199750, 0.0014140010, 0.0005330000,
               0.0000750000])  # fmt: skip
NSPA = [1, 1, 9, 9, 9, 1, 9, 1, 9, 1, 1, 9, 9, 1, 9, 9]
NSPB = [1, 1, 5, 5, 5, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0]
OFF = np.cumsum([0, *NGC])


def groups(b):
    g = NGN[OFF[b] : OFF[b + 1]]
    assert sum(g) == 16, b
    return np.cumsum([0, *g])


def reduce_k(a, b, axis):
    """cmbgb: sum_k k * rwgt over each group; rwgt = wt / (sum of wt in the group)."""
    e = groups(b)
    parts = []
    for i in range(len(e) - 1):
        w = WT[e[i] : e[i + 1]] / WT[e[i] : e[i + 1]].sum()
        sl = [slice(None)] * a.ndim
        sl[axis] = slice(e[i], e[i + 1])
        shp = [1] * a.ndim
        shp[axis] = w.size
        parts.append(np.sum(a[tuple(sl)] * w.reshape(shp), axis=axis))
    return np.stack(parts, axis=axis)


def reduce_f(a, b, axis):
    """cmbgb: Planck fractions summed, unweighted."""
    e = groups(b)
    sl = lambda i: tuple(slice(e[i], e[i + 1]) if d == axis else slice(None) for d in range(a.ndim))  # noqa: E731
    return np.stack([np.sum(a[sl(i)], axis=axis) for i in range(len(e) - 1)], axis=axis)


out = {}
for b in range(16):
    n = b + 1
    kl, ku = NSPA[b], max(NSPB[b], 1)
    # Key species, lower: nc [set, band, g, p13, T5, eta9] -> keep eta[:kl]; reduce g.
    ka = reduce_k(V["KeySpeciesAbsorptionCoefficientsLowerAtmos"][0, b][..., :kl], b, 0)
    # absa(ind, ig): Fortran ka(eta, T, p, ig) flattened column-major == C-order of [p, T, eta].
    out[f"absa{n}"] = ka.reshape(ka.shape[0], -1).T.copy()  # (flat, ng)
    if NSPB[b] > 0 or n in (16,):
        kb = reduce_k(V["KeySpeciesAbsorptionCoefficientsUpperAtmos"][0, b][..., :ku], b, 0)
        out[f"absb{n}"] = kb.reshape(kb.shape[0], -1).T.copy()
    out[f"selfref{n}"] = reduce_k(
        V["H20SelfAbsorptionCoefficients"][0, b], b, 0
    ).T.copy()  # (10, ng)
    out[f"forref{n}"] = reduce_k(
        V["H20ForeignAbsorptionCoefficients"][0, b], b, 0
    ).T.copy()  # (4, ng)
    # Planck fractions: nc [set, band, eta, g] -> fracrefa(ig, eta)
    out[f"fracrefa{n}"] = reduce_f(V["PlanckFractionLowerAtmos"][0, b][:kl], b, 1).T.copy()
    out[f"fracrefb{n}"] = reduce_f(V["PlanckFractionUpperAtmos"][0, b][:ku], b, 1).T.copy()

# Minor gases: nc [set, band, absorber, g, T19, eta9] -> Fortran (eta, T, ig) or (T, ig).
MINORS = {  # band: [(module name, absorber, lower/upper, eta count)]
    1: [("ka_mn2", "N2", "L", 1), ("kb_mn2", "N2", "U", 1)],
    3: [("ka_mn2o", "N2O", "L", 9), ("kb_mn2o", "N2O", "U", 5)],
    5: [("ka_mo3", "O3", "L", 9)],
    6: [("ka_mco2", "CO2", "L", 1)],
    7: [("ka_mco2", "CO2", "L", 9), ("kb_mco2", "CO2", "U", 1)],
    8: [("ka_mo3", "O3", "L", 1), ("ka_mco2", "CO2", "L", 1), ("kb_mco2", "CO2", "U", 1),
        ("ka_mn2o", "N2O", "L", 1), ("kb_mn2o", "N2O", "U", 1)],
    9: [("ka_mn2o", "N2O", "L", 9), ("kb_mn2o", "N2O", "U", 1)],
    11: [("ka_mo2", "O2", "L", 1), ("kb_mo2", "O2", "U", 1)],
    13: [("kb_mo3", "O3", "U", 1), ("ka_mco2", "CO2", "L", 9), ("ka_mco", "CO", "L", 9)],
    15: [("ka_mn2", "N2", "L", 9)],
}  # fmt: skip
CROSS = {5: [("ccl4", "CCL4")], 6: [("cfc11adj", "CFC11"), ("cfc12", "CFC12")],
         8: [("cfc12", "CFC12"), ("cfc22adj", "CFC22")]}  # fmt: skip
for n, items in MINORS.items():
    b = n - 1
    for mod, gas, lu, ne in items:
        src = V[
            "AbsorptionCoefficientsLowerAtmos" if lu == "L" else "AbsorptionCoefficientsUpperAtmos"
        ]
        a = src[0, b, names.index(gas)][..., :ne]  # (g, T19, eta)
        r = reduce_k(a, b, 0)  # (ng, T, eta)
        out[f"{mod}{n}"] = np.transpose(r, (2, 1, 0)).copy() if ne > 1 else r[..., 0].T.copy()
for n, items in CROSS.items():
    b = n - 1
    for mod, gas in items:
        a = V["AbsorptionCoefficientsLowerAtmos"][0, b, names.index(gas)][:, 0, 0]  # (g,)
        out[f"{mod}{n}"] = reduce_k(a, b, 0)

# ---- reference atmosphere and Planck functions ------------------------------------
# From the netCDF, where AER carries the same values lwatmref and lwavplank set in
# rrtmg_lw_setcoef.f90 (equal to the last bit, checked when this was written).
# preflog is the Fortran's literal: ln(pref(1)) = 6.96, falling by exactly 0.2.
pref = V["Pressure"].astype(float)
preflog = np.round(6.96 - 0.2 * np.arange(59), 4)
tref = V["ReferenceTemperature"].astype(float)
chi = np.array([V["AbsorberAmountMLS"][:, names.index(g)] for g in
                ("H2O", "CO2", "O3", "N2O", "CO", "CH4", "O2")])  # fmt: skip
totplnk = V["IntegratedPlanckFunction"].T.astype(float).copy()  # (181, 16)
totplk16 = V["IntegratedPlanckFunctionBand16"].astype(float)
checks = {"preflog_vs_ln_pref": float(np.max(np.abs(preflog - np.log(pref))))}
WAVENUM1 = [10, 350, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080, 2250, 2380, 2600]
WAVENUM2 = [*WAVENUM1[1:], 3250]
out.update(
    pref=pref,
    preflog=preflog,
    tref=tref,
    chi_mls=chi,
    totplnk=totplnk,
    totplk16=totplk16,
    ngc=np.array(NGC),
    nspa=np.array(NSPA),
    nspb=np.array(NSPB),
    wavenum1=np.array(WAVENUM1, dtype=float),
    wavenum2=np.array(WAVENUM2, dtype=float),
)
np.savez_compressed(sys.argv[2], **out)
print({k: f"{v:.2e}" for k, v in checks.items()})
print(len(out), "arrays;", sum(v.size for v in out.values()), "values")
