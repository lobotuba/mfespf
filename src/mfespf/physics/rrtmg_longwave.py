"""RRTMG's longwave gas optics (Mlawer et al. 1997; Iacono et al. 2008).

Why this exists
---------------

`global-column-chou-longwave.md` R1 failed G4, and under the user's decision
recorded there the next step is RRTMG's gas optics
(``docs/protocols/rrtmg-longwave.md``). Its k-distributions are derived from
LBLRTM line-by-line spectra with the g-points ordered consistently across
levels, which is what a solver that treats each bin layer by layer needs.

Where the numbers and the code come from
----------------------------------------

AER-RC/RRTMG_LW at commit 958edb1ad1618d76d0fcc5a9dc2e3c2f791a54dd, BSD 3-Clause,
Copyright Atmospheric and Environmental Research; the repository's LICENSE.md is
reproduced in ``data/rrtmg_lw_reduced.LICENSE``:

- ``data/rrtmg_lw.nc``, g-point set 1, read as ``rrtmg_lw_read_nc.f90`` reads it
  and reduced from 256 to 140 g-points exactly as ``cmbgb1``-``cmbgb16`` do;
- the reference atmosphere and the integrated Planck functions from
  ``rrtmg_lw_setcoef.f90``.

The tables are generated into ``data/rrtmg_lw_reduced.npz`` and not edited.
``setcoef`` and ``taumol`` below are line-by-line ports of
``rrtmg_lw_setcoef.f90`` and ``rrtmg_lw_taumol.f90``, vectorised over layers.
Indices keep the Fortran's 1-based arithmetic and subtract one where an array is
read, so each line can be checked against its source.

What it computes
----------------

For each layer, the gaseous optical depth at each of the 140 g-points and
RRTMG's Planck fraction for each, from the layer pressure and temperature and
the molecular amounts. The optical depths are vertical (RRTMG's transfer applies
its own diffusivity); nothing here is fitted or adjusted.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np

from mfespf.contracts.spectral import SpectralCoordinate, SpectralGrid
from mfespf.core.errors import ConfigurationError

__all__ = [
    "BAND_WAVENUMBER_CM",
    "G_POINTS_PER_BAND",
    "RRTMGColumn",
    "RRTMGLongwave",
    "dry_air_column",
    "rrtmg_fixed_fractions",
    "rrtmg_gas_optics",
    "rrtmg_longwave_grid",
    "rrtmg_optical_depth",
    "rrtmg_tables",
]

_TABLE_PATH = Path(__file__).with_name("data") / "rrtmg_lw_reduced.npz"
#: g-points per band after the 256-to-140 reduction (lwcmbdat ``ngc``).
G_POINTS_PER_BAND: tuple[int, ...] = (10, 12, 16, 14, 16, 8, 12, 8, 12, 6, 8, 8, 4, 2, 2, 2)
_NGS = np.concatenate(([0], np.cumsum(G_POINTS_PER_BAND)))
#: Band edges, cm^-1 (lwdatinit ``wavenum1``/``wavenum2``).
BAND_WAVENUMBER_CM: tuple[tuple[float, float], ...] = (
    (10.0, 350.0),
    (350.0, 500.0),
    (500.0, 630.0),
    (630.0, 700.0),
    (700.0, 820.0),
    (820.0, 980.0),
    (980.0, 1080.0),
    (1080.0, 1180.0),
    (1180.0, 1390.0),
    (1390.0, 1480.0),
    (1480.0, 1800.0),
    (1800.0, 2080.0),
    (2080.0, 2250.0),
    (2250.0, 2380.0),
    (2380.0, 2600.0),
    (2600.0, 3250.0),
)
_ONEMINUS = 1.0 - 1.0e-6
_NSPA = (1, 1, 9, 9, 9, 1, 9, 1, 9, 1, 1, 9, 9, 1, 9, 9)
_NSPB = (1, 1, 5, 5, 5, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0)


@lru_cache(maxsize=1)
def rrtmg_tables() -> dict[str, np.ndarray]:
    """The reduced tables, read once."""
    if not _TABLE_PATH.exists():
        raise ConfigurationError(f"RRTMG tables not found at {_TABLE_PATH}")
    with np.load(_TABLE_PATH) as data:
        return {k: np.asarray(data[k]) for k in data.files}


@dataclass(frozen=True)
class RRTMGColumn:
    """A column in RRTMG's units, layers ordered from the surface up.

    ``pressure_hpa`` and ``temperature`` are layer means. ``coldry`` is the dry
    air column (molecules cm^-2); ``wkl`` holds H2O, CO2, O3, N2O, CO, CH4 and O2
    in molecules cm^-2, shaped (7, layers); ``wbroad`` is the broadening-gas
    column; ``wx`` holds CCl4, CFC-11, CFC-12 and CFC-22 as RRTMG's driver
    stores them, coldry x vmr x 1e-20, shaped (4, layers).
    """

    pressure_hpa: np.ndarray
    temperature: np.ndarray
    coldry: np.ndarray
    wkl: np.ndarray
    wbroad: np.ndarray
    wx: np.ndarray

    @classmethod
    def from_mixing_ratios(
        cls,
        *,
        pressure_hpa: np.ndarray,
        temperature: np.ndarray,
        coldry: np.ndarray,
        vmr: dict[str, np.ndarray | float],
    ) -> RRTMGColumn:
        """As ``rrtmg_lw_rad``'s ``inatm``: amounts from volume mixing ratios.

        ``vmr`` keys: h2o, co2, o3, n2o, co, ch4, o2, ccl4, cfc11, cfc12, cfc22;
        any absent is zero. The broadening gas is what is left of the dry air.
        """
        n = np.asarray(pressure_hpa).size
        names = ("h2o", "co2", "o3", "n2o", "co", "ch4", "o2")
        mix = np.array([np.broadcast_to(np.asarray(vmr.get(k, 0.0), float), (n,)) for k in names])
        dry = np.asarray(coldry, dtype=float)
        wbroad = dry * (1.0 - np.sum(mix[1:], axis=0))
        xs = ("ccl4", "cfc11", "cfc12", "cfc22")
        wx = np.array(
            [dry * np.broadcast_to(np.asarray(vmr.get(k, 0.0), float), (n,)) * 1.0e-20 for k in xs]
        )
        return cls(
            pressure_hpa=np.asarray(pressure_hpa, dtype=float),
            temperature=np.asarray(temperature, dtype=float),
            coldry=dry,
            wkl=dry * mix,
            wbroad=wbroad,
            wx=wx,
        )


def _int(x: np.ndarray) -> np.ndarray:
    """Fortran ``int()``: truncation toward zero."""
    return np.trunc(x).astype(int)


def _setcoef(col: RRTMGColumn) -> dict[str, np.ndarray]:
    """``setcoef``: interpolation indices, fractions and scaled amounts per layer."""
    t = rrtmg_tables()
    preflog, tref, chi = t["preflog"], t["tref"], t["chi_mls"]
    pavel, tavel, coldry, wkl = col.pressure_hpa, col.temperature, col.coldry, col.wkl
    stpfac = 296.0 / 1013.0
    plog = np.log(pavel)
    jp = np.clip(_int(36.0 - 5 * (plog + 0.04)), 1, 58)
    jp1 = jp + 1
    fp = 5.0 * (preflog[jp - 1] - plog)
    jt = np.clip(_int(3.0 + (tavel - tref[jp - 1]) / 15.0), 1, 4)
    ft = ((tavel - tref[jp - 1]) / 15.0) - (jt - 3)
    jt1 = np.clip(_int(3.0 + (tavel - tref[jp1 - 1]) / 15.0), 1, 4)
    ft1 = ((tavel - tref[jp1 - 1]) / 15.0) - (jt1 - 3)
    water = wkl[0] / coldry
    scalefac = pavel * stpfac / tavel
    lower = plog > 4.56
    laytrop = int(np.sum(lower))
    if not np.all(lower[:laytrop]) or np.any(lower[laytrop:]):
        raise ConfigurationError("RRTMG needs layers ordered from the surface up")

    forfac = scalefac / (1.0 + water)
    factor_l = (332.0 - tavel) / 36.0
    indfor = np.where(lower, np.minimum(2, np.maximum(1, _int(factor_l))), 3)
    forfrac = np.where(lower, factor_l - indfor, (tavel - 188.0) / 36.0 - 1.0)
    selffac = water * forfac
    factor_s = (tavel - 188.0) / 7.2
    indself = np.minimum(9, np.maximum(1, _int(factor_s) - 7))
    selffrac = factor_s - (indself + 7)
    scaleminor = pavel / tavel
    scaleminorn2 = (pavel / tavel) * (col.wbroad / (coldry + wkl[0]))
    factor_m = (tavel - 180.8) / 7.2
    indminor = np.minimum(18, np.maximum(1, _int(factor_m)))
    minorfrac = factor_m - indminor

    def rat(a: int, b: int, k: np.ndarray) -> np.ndarray:
        return chi[a - 1, k - 1] / chi[b - 1, k - 1]

    c = {
        "laytrop": np.array(laytrop),
        "lower": lower,
        "jp": jp,
        "jt": jt,
        "jt1": jt1,
        "rat_h2oco2": rat(1, 2, jp),
        "rat_h2oco2_1": rat(1, 2, jp1),
        "rat_h2oo3": rat(1, 3, jp),
        "rat_h2oo3_1": rat(1, 3, jp1),
        "rat_h2on2o": rat(1, 4, jp),
        "rat_h2on2o_1": rat(1, 4, jp1),
        "rat_h2och4": rat(1, 6, jp),
        "rat_h2och4_1": rat(1, 6, jp1),
        "rat_n2oco2": rat(4, 2, jp),
        "rat_n2oco2_1": rat(4, 2, jp1),
        "rat_o3co2": rat(3, 2, jp),
        "rat_o3co2_1": rat(3, 2, jp1),
        "indfor": indfor,
        "forfrac": forfrac,
        "indself": indself,
        "selffrac": selffrac,
        "indminor": indminor,
        "minorfrac": minorfrac,
        "scaleminor": scaleminor,
        "scaleminorn2": scaleminorn2,
    }
    for i, name in enumerate(("colh2o", "colco2", "colo3", "coln2o", "colco", "colch4", "colo2")):
        amount = 1.0e-20 * wkl[i]
        if name not in ("colh2o", "colo2"):
            amount = np.where(amount == 0.0, 1.0e-32 * coldry, amount)
        c[name] = amount
    c["colbrd"] = 1.0e-20 * col.wbroad
    compfp = 1.0 - fp
    c["fac10"] = compfp * ft
    c["fac00"] = compfp * (1.0 - ft)
    c["fac11"] = fp * ft1
    c["fac01"] = fp * (1.0 - ft1)
    c["selffac"] = c["colh2o"] * selffac
    c["forfac"] = c["colh2o"] * forfac
    c["coldry"] = coldry
    c["pavel"] = pavel
    return c


# ---- helpers shared by the band routines --------------------------------------------


def _take(table: np.ndarray, index_1based: np.ndarray) -> np.ndarray:
    """Rows of a (flat, ng) table at 1-based indices, clipped for unselected branches."""
    return table[np.clip(index_1based - 1, 0, table.shape[0] - 1)]


def _self(t: dict, n: int, c: dict, L: np.ndarray) -> np.ndarray:
    s = t[f"selfref{n}"]
    i = c["indself"][L]
    return c["selffac"][L, None] * (s[i - 1] + c["selffrac"][L, None] * (s[i] - s[i - 1]))


def _for(t: dict, n: int, c: dict, L: np.ndarray) -> np.ndarray:
    f = t[f"forref{n}"]
    i = c["indfor"][L]
    return c["forfac"][L, None] * (f[i - 1] + c["forfrac"][L, None] * (f[i] - f[i - 1]))


def _minor(m: np.ndarray, c: dict, L: np.ndarray) -> np.ndarray:
    """``ka_m(indm,ig) + minorfrac*(ka_m(indm+1,ig) - ka_m(indm,ig))``; m is (19, ng)."""
    i = c["indminor"][L]
    return m[i - 1] + c["minorfrac"][L, None] * (m[i] - m[i - 1])


def _minor_eta(m: np.ndarray, j: np.ndarray, f: np.ndarray, c: dict, L: np.ndarray) -> np.ndarray:
    """The minor gas interpolated in its binary parameter too; m is (eta, 19, ng)."""
    i = c["indminor"][L]
    f = f[:, None]
    m1 = m[j - 1, i - 1] + f * (m[j, i - 1] - m[j - 1, i - 1])
    m2 = m[j - 1, i] + f * (m[j, i] - m[j - 1, i])
    return m1 + c["minorfrac"][L, None] * (m2 - m1)


def _binary(
    a: np.ndarray, b: np.ndarray, ratio: np.ndarray | float, mult: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """speccomb, specparm (capped at oneminus), j = 1 + int(mult*specparm), and its fraction."""
    speccomb = a + ratio * b
    specparm = np.minimum(a / speccomb, _ONEMINUS)
    specmult = mult * specparm
    return speccomb, specparm, 1 + _int(specmult), np.mod(specmult, 1.0)


def _fk(fs: np.ndarray, specparm: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    p = np.where(specparm < 0.125, fs - 1.0, -fs)
    p4 = p**4
    return p4, 1.0 - p - 2.0 * p4, p + p4


def _major_lower(
    absa: np.ndarray,
    ind: np.ndarray,
    speccomb: np.ndarray,
    specparm: np.ndarray,
    fs: np.ndarray,
    fa: np.ndarray,
    fb: np.ndarray,
) -> np.ndarray:
    """One of tau_major/tau_major1 in a two-key lower band (nspa = 9).

    ``fa``/``fb`` are fac00/fac10 (or fac01/fac11): the three branches of the
    Fortran, selected per layer.
    """
    fk0, fk1, fk2 = _fk(fs, specparm)
    A = lambda k: _take(absa, ind + k)  # noqa: E731
    fa_, fb_ = fa[:, None], fb[:, None]
    fk0, fk1, fk2, fs_ = fk0[:, None], fk1[:, None], fk2[:, None], fs[:, None]
    low = fk0 * fa_ * A(0) + fk1 * fa_ * A(1) + fk2 * fa_ * A(2)
    low += fk0 * fb_ * A(9) + fk1 * fb_ * A(10) + fk2 * fb_ * A(11)
    high = fk2 * fa_ * A(-1) + fk1 * fa_ * A(0) + fk0 * fa_ * A(1)
    high += fk2 * fb_ * A(8) + fk1 * fb_ * A(9) + fk0 * fb_ * A(10)
    mid = (1.0 - fs_) * fa_ * A(0) + fs_ * fa_ * A(1) + (1.0 - fs_) * fb_ * A(9) + fs_ * fb_ * A(10)
    sp = specparm[:, None]
    return speccomb[:, None] * np.where(sp < 0.125, low, np.where(sp > 0.875, high, mid))


def _major_upper(
    absb: np.ndarray, ind: np.ndarray, speccomb: np.ndarray, fs: np.ndarray,
    fa: np.ndarray, fb: np.ndarray,
) -> np.ndarray:  # fmt: skip
    """A two-key upper band (nspb = 5): linear in the binary parameter."""
    fs_, fa_, fb_ = fs[:, None], fa[:, None], fb[:, None]
    B = lambda k: _take(absb, ind + k)  # noqa: E731
    return speccomb[:, None] * (
        (1.0 - fs_) * fa_ * B(0) + fs_ * fa_ * B(1) + (1.0 - fs_) * fb_ * B(5) + fs_ * fb_ * B(6)
    )


def _single(table: np.ndarray, c: dict, L: np.ndarray, upper: bool, nsp: int = 1) -> np.ndarray:
    """``fac00*abs(ind0) + fac10*abs(ind0+1) + fac01*abs(ind1) + fac11*abs(ind1+1)``.

    ``nsp`` is the band's nspa or nspb, as the Fortran multiplies by it. It is
    1 for every single-key band but band 16's upper atmosphere, where nspb is 0
    and the Fortran therefore reads the table's first entries at every level.
    That is reproduced, not corrected.
    """
    jp, jt, jt1 = c["jp"][L], c["jt"][L], c["jt1"][L]
    if upper:
        ind0 = ((jp - 13) * 5 + (jt - 1)) * nsp + 1
        ind1 = ((jp - 12) * 5 + (jt1 - 1)) * nsp + 1
    else:
        ind0 = ((jp - 1) * 5 + (jt - 1)) * nsp + 1
        ind1 = (jp * 5 + (jt1 - 1)) * nsp + 1
    f = {k: c[k][L, None] for k in ("fac00", "fac10", "fac01", "fac11")}
    return (
        f["fac00"] * _take(table, ind0)
        + f["fac10"] * _take(table, ind0 + 1)
        + f["fac01"] * _take(table, ind1)
        + f["fac11"] * _take(table, ind1 + 1)
    )


def _two_key_lower(
    t: dict, n: int, c: dict, L: np.ndarray, a: str, b: str, rat: str
) -> tuple[np.ndarray, tuple]:
    """tau_major + tau_major1 for a two-key lower band, and its binary parameters."""
    A, B = c[a][L], c[b][L]
    sc, sp, js, fs = _binary(A, B, c[rat][L], 8.0)
    sc1, sp1, js1, fs1 = _binary(A, B, c[f"{rat}_1"][L], 8.0)
    jp, jt, jt1 = c["jp"][L], c["jt"][L], c["jt1"][L]
    ind0 = ((jp - 1) * 5 + (jt - 1)) * _NSPA[n - 1] + js
    ind1 = (jp * 5 + (jt1 - 1)) * _NSPA[n - 1] + js1
    absa = t[f"absa{n}"]
    tau = _major_lower(absa, ind0, sc, sp, fs, c["fac00"][L], c["fac10"][L])
    tau += _major_lower(absa, ind1, sc1, sp1, fs1, c["fac01"][L], c["fac11"][L])
    return tau, (A, B)


def _two_key_upper(t: dict, n: int, c: dict, L: np.ndarray, a: str, b: str, rat: str) -> np.ndarray:
    A, B = c[a][L], c[b][L]
    sc, _, js, fs = _binary(A, B, c[rat][L], 4.0)
    sc1, _, js1, fs1 = _binary(A, B, c[f"{rat}_1"][L], 4.0)
    jp, jt, jt1 = c["jp"][L], c["jt"][L], c["jt1"][L]
    ind0 = ((jp - 13) * 5 + (jt - 1)) * _NSPB[n - 1] + js
    ind1 = ((jp - 12) * 5 + (jt1 - 1)) * _NSPB[n - 1] + js1
    absb = t[f"absb{n}"]
    return _major_upper(absb, ind0, sc, fs, c["fac00"][L], c["fac10"][L]) + _major_upper(
        absb, ind1, sc1, fs1, c["fac01"][L], c["fac11"][L]
    )


def _planck_eta(
    F: np.ndarray, a: np.ndarray, b: np.ndarray, ratio: float, mult: float
) -> np.ndarray:
    """fracrefa(ig,jpl) + fpl*(fracrefa(ig,jpl+1) - fracrefa(ig,jpl)); F is (ng, eta)."""
    _, _, jpl, fpl = _binary(a, b, ratio, mult)
    return F[:, jpl - 1].T + fpl[:, None] * (F[:, jpl].T - F[:, jpl - 1].T)


def _adjusted(
    col: np.ndarray, c: dict, L: np.ndarray, species: int, threshold: float, base: float,
    power: float, reference: np.ndarray | None = None,
) -> np.ndarray:  # fmt: skip
    """The empirical column adjustment for a minor gas present in quantity."""
    ref = c["chi"][species - 1, c["jp"][L]] if reference is None else reference
    ratio = 1.0e20 * (col / c["coldry"][L]) / ref
    adj = (base + np.maximum(ratio - base, 0.0) ** power) * ref * c["coldry"][L] * 1.0e-20
    return np.where(ratio > threshold, adj, col)


# ---- taumol ------------------------------------------------------------------------


def _taumol(col: RRTMGColumn, c: dict) -> tuple[np.ndarray, np.ndarray]:
    t = rrtmg_tables()
    chi = t["chi_mls"]
    c["chi"] = chi
    nlay = col.pressure_hpa.size
    taug = np.zeros((nlay, int(_NGS[-1])))
    fracs = np.zeros_like(taug)
    lo = np.flatnonzero(c["lower"])
    up = np.flatnonzero(~c["lower"])
    wx = col.wx

    def put(n: int, L: np.ndarray, tau: np.ndarray, frac: np.ndarray | float) -> None:
        s = slice(int(_NGS[n - 1]), int(_NGS[n]))
        taug[L, s] = tau
        fracs[L, s] = np.broadcast_to(frac, tau.shape)

    def R(a: int, b: int, k: int) -> float:  # noqa: N802 -- chi_mls ratio, as the Fortran's refrat
        return float(chi[a - 1, k - 1] / chi[b - 1, k - 1])

    # Band 1: 10-350 cm-1 (low key h2o, minor n2; high key h2o, minor n2).
    n = 1
    L = lo
    pp = c["pavel"][L]
    corradj = np.where(pp < 250.0, 1.0 - 0.15 * (250.0 - pp) / 154.4, 1.0)[:, None]
    scalen2 = (c["colbrd"][L] * c["scaleminorn2"][L])[:, None]
    tau = corradj * (
        c["colh2o"][L, None] * _single(t["absa1"], c, L, False)
        + _self(t, n, c, L) + _for(t, n, c, L) + scalen2 * _minor(t["ka_mn21"], c, L)
    )  # fmt: skip
    put(n, L, tau, t["fracrefa1"][:, 0])
    L = up
    corradj = (1.0 - 0.15 * (c["pavel"][L] / 95.6))[:, None]
    scalen2 = (c["colbrd"][L] * c["scaleminorn2"][L])[:, None]
    tau = corradj * (
        c["colh2o"][L, None] * _single(t["absb1"], c, L, True)
        + _for(t, n, c, L) + scalen2 * _minor(t["kb_mn21"], c, L)
    )  # fmt: skip
    put(n, L, tau, t["fracrefb1"][:, 0])

    # Band 2: 350-500 cm-1 (h2o).
    n = 2
    L = lo
    corradj = (1.0 - 0.05 * (c["pavel"][L] - 100.0) / 900.0)[:, None]
    tau = corradj * (
        c["colh2o"][L, None] * _single(t["absa2"], c, L, False)
        + _self(t, n, c, L) + _for(t, n, c, L)
    )  # fmt: skip
    put(n, L, tau, t["fracrefa2"][:, 0])
    L = up
    tau = c["colh2o"][L, None] * _single(t["absb2"], c, L, True) + _for(t, n, c, L)
    put(n, L, tau, t["fracrefb2"][:, 0])

    # Band 3: 500-630 cm-1 (h2o,co2; minor n2o).
    n = 3
    L = lo
    tau, (A, B) = _two_key_lower(t, n, c, L, "colh2o", "colco2", "rat_h2oco2")
    _, _, jm, fm = _binary(A, B, R(1, 2, 3), 8.0)
    adjn2o = _adjusted(c["coln2o"][L], c, L, 4, 1.5, 0.5, 0.65)
    tau += _self(t, n, c, L) + _for(t, n, c, L)
    tau += adjn2o[:, None] * _minor_eta(t["ka_mn2o3"], jm, fm, c, L)
    put(n, L, tau, _planck_eta(t["fracrefa3"], A, B, R(1, 2, 9), 8.0))
    L = up
    A, B = c["colh2o"][L], c["colco2"][L]
    tau = _two_key_upper(t, n, c, L, "colh2o", "colco2", "rat_h2oco2")
    _, _, jm, fm = _binary(A, B, R(1, 2, 13), 4.0)
    adjn2o = _adjusted(c["coln2o"][L], c, L, 4, 1.5, 0.5, 0.65)
    tau += _for(t, n, c, L) + adjn2o[:, None] * _minor_eta(t["kb_mn2o3"], jm, fm, c, L)
    put(n, L, tau, _planck_eta(t["fracrefb3"], A, B, R(1, 2, 13), 4.0))

    # Band 4: 630-700 cm-1 (low h2o,co2; high o3,co2).
    n = 4
    L = lo
    tau, (A, B) = _two_key_lower(t, n, c, L, "colh2o", "colco2", "rat_h2oco2")
    tau += _self(t, n, c, L) + _for(t, n, c, L)
    put(n, L, tau, _planck_eta(t["fracrefa4"], A, B, R(1, 2, 11), 8.0))
    L = up
    tau = _two_key_upper(t, n, c, L, "colo3", "colco2", "rat_o3co2")
    # Empirical modification for stratospheric CO2 cooling (g-points 8-14).
    tau[:, 7:14] *= np.array([0.92, 0.88, 1.07, 1.1, 0.99, 0.88, 0.943])
    put(n, L, tau, _planck_eta(t["fracrefb4"], c["colo3"][L], c["colco2"][L], R(3, 2, 13), 4.0))

    # Band 5: 700-820 cm-1 (low h2o,co2, minor o3, ccl4; high o3,co2).
    n = 5
    L = lo
    tau, (A, B) = _two_key_lower(t, n, c, L, "colh2o", "colco2", "rat_h2oco2")
    _, _, jm, fm = _binary(A, B, R(1, 2, 7), 8.0)
    tau += _self(t, n, c, L) + _for(t, n, c, L)
    tau += c["colo3"][L, None] * _minor_eta(t["ka_mo35"], jm, fm, c, L)
    tau += wx[0, L, None] * t["ccl45"][None, :]
    put(n, L, tau, _planck_eta(t["fracrefa5"], A, B, R(1, 2, 5), 8.0))
    L = up
    tau = _two_key_upper(t, n, c, L, "colo3", "colco2", "rat_o3co2")
    tau += wx[0, L, None] * t["ccl45"][None, :]
    put(n, L, tau, _planck_eta(t["fracrefb5"], c["colo3"][L], c["colco2"][L], R(3, 2, 43), 4.0))

    # Band 6: 820-980 cm-1 (low h2o, minor co2; high cfc11, cfc12).
    n = 6
    L = lo
    adjco2 = _adjusted(c["colco2"][L], c, L, 2, 3.0, 2.0, 0.77)
    tau = (
        c["colh2o"][L, None] * _single(t["absa6"], c, L, False)
        + _self(t, n, c, L) + _for(t, n, c, L)
        + adjco2[:, None] * _minor(t["ka_mco26"], c, L)
        + wx[1, L, None] * t["cfc11adj6"] + wx[2, L, None] * t["cfc126"]
    )  # fmt: skip
    put(n, L, tau, t["fracrefa6"][:, 0])
    L = up
    tau = wx[1, L, None] * t["cfc11adj6"] + wx[2, L, None] * t["cfc126"]
    put(n, L, tau, t["fracrefa6"][:, 0])

    # Band 7: 980-1080 cm-1 (low h2o,o3, minor co2; high o3, minor co2).
    n = 7
    L = lo
    tau, (A, B) = _two_key_lower(t, n, c, L, "colh2o", "colo3", "rat_h2oo3")
    _, _, jm, fm = _binary(A, B, R(1, 3, 3), 8.0)
    adjco2 = _adjusted(c["colco2"][L], c, L, 2, 3.0, 3.0, 0.79)
    tau += _self(t, n, c, L) + _for(t, n, c, L)
    tau += adjco2[:, None] * _minor_eta(t["ka_mco27"], jm, fm, c, L)
    put(n, L, tau, _planck_eta(t["fracrefa7"], A, B, R(1, 3, 3), 8.0))
    L = up
    adjco2 = _adjusted(c["colco2"][L], c, L, 2, 3.0, 2.0, 0.79)
    tau = c["colo3"][L, None] * _single(t["absb7"], c, L, True)
    tau += adjco2[:, None] * _minor(t["kb_mco27"], c, L)
    # Empirical modification for stratospheric O3 cooling (g-points 6-11).
    tau[:, 5:11] *= np.array([0.92, 0.88, 1.07, 1.1, 0.99, 0.855])
    put(n, L, tau, t["fracrefb7"][:, 0])

    # Band 8: 1080-1180 cm-1 (low h2o, minor co2,o3,n2o; high o3, minor co2,n2o).
    n = 8
    L = lo
    adjco2 = _adjusted(c["colco2"][L], c, L, 2, 3.0, 2.0, 0.65)
    tau = (
        c["colh2o"][L, None] * _single(t["absa8"], c, L, False)
        + _self(t, n, c, L) + _for(t, n, c, L)
        + adjco2[:, None] * _minor(t["ka_mco28"], c, L)
        + c["colo3"][L, None] * _minor(t["ka_mo38"], c, L)
        + c["coln2o"][L, None] * _minor(t["ka_mn2o8"], c, L)
        + wx[2, L, None] * t["cfc128"] + wx[3, L, None] * t["cfc22adj8"]
    )  # fmt: skip
    put(n, L, tau, t["fracrefa8"][:, 0])
    L = up
    adjco2 = _adjusted(c["colco2"][L], c, L, 2, 3.0, 2.0, 0.65)
    tau = (
        c["colo3"][L, None] * _single(t["absb8"], c, L, True)
        + adjco2[:, None] * _minor(t["kb_mco28"], c, L)
        + c["coln2o"][L, None] * _minor(t["kb_mn2o8"], c, L)
        + wx[2, L, None] * t["cfc128"] + wx[3, L, None] * t["cfc22adj8"]
    )  # fmt: skip
    put(n, L, tau, t["fracrefb8"][:, 0])

    # Band 9: 1180-1390 cm-1 (low h2o,ch4, minor n2o; high ch4, minor n2o).
    n = 9
    L = lo
    tau, (A, B) = _two_key_lower(t, n, c, L, "colh2o", "colch4", "rat_h2och4")
    _, _, jm, fm = _binary(A, B, R(1, 6, 3), 8.0)
    adjn2o = _adjusted(c["coln2o"][L], c, L, 4, 1.5, 0.5, 0.65)
    tau += _self(t, n, c, L) + _for(t, n, c, L)
    tau += adjn2o[:, None] * _minor_eta(t["ka_mn2o9"], jm, fm, c, L)
    put(n, L, tau, _planck_eta(t["fracrefa9"], A, B, R(1, 6, 9), 8.0))
    L = up
    adjn2o = _adjusted(c["coln2o"][L], c, L, 4, 1.5, 0.5, 0.65)
    tau = c["colch4"][L, None] * _single(t["absb9"], c, L, True)
    tau += adjn2o[:, None] * _minor(t["kb_mn2o9"], c, L)
    put(n, L, tau, t["fracrefb9"][:, 0])

    # Band 10: 1390-1480 cm-1 (h2o).
    n = 10
    L = lo
    tau = c["colh2o"][L, None] * _single(t["absa10"], c, L, False)
    tau += _self(t, n, c, L) + _for(t, n, c, L)
    put(n, L, tau, t["fracrefa10"][:, 0])
    L = up
    tau = c["colh2o"][L, None] * _single(t["absb10"], c, L, True) + _for(t, n, c, L)
    put(n, L, tau, t["fracrefb10"][:, 0])

    # Band 11: 1480-1800 cm-1 (h2o, minor o2).
    n = 11
    L = lo
    scaleo2 = (c["colo2"][L] * c["scaleminor"][L])[:, None]
    tau = c["colh2o"][L, None] * _single(t["absa11"], c, L, False)
    tau += _self(t, n, c, L) + _for(t, n, c, L) + scaleo2 * _minor(t["ka_mo211"], c, L)
    put(n, L, tau, t["fracrefa11"][:, 0])
    L = up
    scaleo2 = (c["colo2"][L] * c["scaleminor"][L])[:, None]
    tau = c["colh2o"][L, None] * _single(t["absb11"], c, L, True)
    tau += _for(t, n, c, L) + scaleo2 * _minor(t["kb_mo211"], c, L)
    put(n, L, tau, t["fracrefb11"][:, 0])

    # Band 12: 1800-2080 cm-1 (low h2o,co2; nothing above).
    n = 12
    L = lo
    tau, (A, B) = _two_key_lower(t, n, c, L, "colh2o", "colco2", "rat_h2oco2")
    tau += _self(t, n, c, L) + _for(t, n, c, L)
    put(n, L, tau, _planck_eta(t["fracrefa12"], A, B, R(1, 2, 10), 8.0))
    put(n, up, np.zeros((up.size, G_POINTS_PER_BAND[n - 1])), 0.0)

    # Band 13: 2080-2250 cm-1 (low h2o,n2o, minor co2,co; high minor o3).
    n = 13
    L = lo
    tau, (A, B) = _two_key_lower(t, n, c, L, "colh2o", "coln2o", "rat_h2on2o")
    _, _, jmco2, fmco2 = _binary(A, B, R(1, 4, 1), 8.0)
    _, _, jmco, fmco = _binary(A, B, R(1, 4, 3), 8.0)
    adjco2 = _adjusted(c["colco2"][L], c, L, 2, 3.0, 2.0, 0.68, reference=np.array(3.55e-4))
    tau += _self(t, n, c, L) + _for(t, n, c, L)
    tau += adjco2[:, None] * _minor_eta(t["ka_mco213"], jmco2, fmco2, c, L)
    tau += c["colco"][L, None] * _minor_eta(t["ka_mco13"], jmco, fmco, c, L)
    put(n, L, tau, _planck_eta(t["fracrefa13"], A, B, R(1, 4, 5), 8.0))
    L = up
    put(n, L, c["colo3"][L, None] * _minor(t["kb_mo313"], c, L), t["fracrefb13"][:, 0])

    # Band 14: 2250-2380 cm-1 (co2).
    n = 14
    L = lo
    tau = c["colco2"][L, None] * _single(t["absa14"], c, L, False)
    tau += _self(t, n, c, L) + _for(t, n, c, L)
    put(n, L, tau, t["fracrefa14"][:, 0])
    L = up
    put(n, L, c["colco2"][L, None] * _single(t["absb14"], c, L, True), t["fracrefb14"][:, 0])

    # Band 15: 2380-2600 cm-1 (low n2o,co2, minor n2; nothing above).
    n = 15
    L = lo
    tau, (A, B) = _two_key_lower(t, n, c, L, "coln2o", "colco2", "rat_n2oco2")
    _, _, jm, fm = _binary(A, B, R(4, 2, 1), 8.0)
    scalen2 = (c["colbrd"][L] * c["scaleminor"][L])[:, None]
    tau += _self(t, n, c, L) + _for(t, n, c, L)
    tau += scalen2 * _minor_eta(t["ka_mn215"], jm, fm, c, L)
    put(n, L, tau, _planck_eta(t["fracrefa15"], A, B, R(4, 2, 1), 8.0))
    put(n, up, np.zeros((up.size, G_POINTS_PER_BAND[n - 1])), 0.0)

    # Band 16: 2600-3250 cm-1 (low h2o,ch4; high ch4).
    n = 16
    L = lo
    tau, (A, B) = _two_key_lower(t, n, c, L, "colh2o", "colch4", "rat_h2och4")
    tau += _self(t, n, c, L) + _for(t, n, c, L)
    put(n, L, tau, _planck_eta(t["fracrefa16"], A, B, R(1, 6, 6), 8.0))
    L = up
    put(
        n,
        L,
        c["colch4"][L, None] * _single(t["absb16"], c, L, True, _NSPB[15]),
        t["fracrefb16"][:, 0],
    )
    return taug, fracs


def rrtmg_gas_optics(col: RRTMGColumn) -> tuple[np.ndarray, np.ndarray, int]:
    """Gaseous optical depth and Planck fraction, each (layers, 140), and laytrop."""
    c = _setcoef(col)
    taug, fracs = _taumol(col, c)
    return taug, fracs, int(c["laytrop"])


# ---- the grid ---------------------------------------------------------------------

#: RRTMG's driver constants (rrtmg_lw_rad: amd, amw; rrtmg_lw_init: grav, avogad).
_AMD, _AMW, _GRAV, _AVOGADRO = 28.9660, 18.0160, 9.8066, 6.02214199e23
_H, _C, _KB = 6.62607015e-34, 2.99792458e8, 1.380649e-23
#: The temperature at which a g-point's sub-interval receives its fraction of the
#: band's Planck flux, as for Chou's terms.
_SPLIT_TEMPERATURE_K = 250.0


def dry_air_column(
    p_lower_hpa: np.ndarray, p_upper_hpa: np.ndarray, h2o_vmr: np.ndarray
) -> np.ndarray:
    """``inatm``'s coldry: dry-air molecules cm^-2 in a layer between two pressures."""
    h2o = np.asarray(h2o_vmr, dtype=float)
    amm = (1.0 - h2o) * _AMD + h2o * _AMW
    dp = np.asarray(p_lower_hpa, dtype=float) - np.asarray(p_upper_hpa, dtype=float)
    return dp * 1.0e3 * _AVOGADRO / (1.0e2 * _GRAV * amm * (1.0 + h2o))


def rrtmg_fixed_fractions() -> np.ndarray:
    """Each g-point's fixed share of its band's Planck flux (rrtmg-longwave.md, erratum).

    RRTMG's own Planck fractions, averaged over the AFGL US Standard Atmosphere by
    layer mass and renormalised within each band; generated by
    ``tools/gen_rrtmg_fractions.py`` into the tables.
    """
    t = rrtmg_tables()
    if "fixed_fractions" not in t:
        raise ConfigurationError(
            "the RRTMG tables carry no fixed fractions; run gen_rrtmg_fractions"
        )
    return t["fixed_fractions"]


def _planck_split(low_cm: float, high_cm: float, weights: np.ndarray) -> np.ndarray:
    """Edges splitting [low, high] cm^-1 into pieces holding these shares of the Planck flux."""
    nu = np.linspace(low_cm, high_cm, 4001)
    x = nu * 100.0
    density = x**3 / np.expm1(_H * _C * x / (_KB * _SPLIT_TEMPERATURE_K))
    cumulative = np.concatenate(
        ([0.0], np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(nu)))
    )
    cumulative = cumulative / cumulative[-1]
    shares = np.cumsum(weights / weights.sum())[:-1]
    return np.concatenate(([low_cm], np.interp(shares, cumulative, nu), [high_cm]))


def rrtmg_longwave_grid(
    grid: SpectralGrid, *, upper_cm: float = 2500.0
) -> tuple[SpectralGrid, np.ndarray]:
    """Rebuild a wavelength grid's longwave on RRTMG's g-points.

    Everything from 10 cm^-1 to ``upper_cm`` (the model's 4-micron shortwave
    boundary by default) becomes one bin per g-point, each a sub-interval of its
    band sized to its fixed Planck fraction at 250 K, in g-point order. A band
    reaching past ``upper_cm`` is cut there; one beginning past it is left out.
    Shortward of ``upper_cm`` the grid is kept edge for edge. Returns the grid
    and each bin's g-point index (-1 outside the longwave).
    """
    if grid.coordinate is not SpectralCoordinate.WAVELENGTH:
        raise ConfigurationError("RRTMG's longwave is built on a wavelength grid")
    edges = np.asarray(grid.edges, dtype=float)
    if edges[0] > 1.0e-2 / upper_cm or edges[-1] < 1.0e-2 / 10.0:
        raise ConfigurationError("the grid must span the longwave it is to hold")
    fractions = rrtmg_fixed_fractions()
    pieces: list[tuple[float, float, int]] = []
    for b, (low, high) in enumerate(BAND_WAVENUMBER_CM):
        if low >= upper_cm:
            continue
        g = np.arange(int(_NGS[b]), int(_NGS[b + 1]))
        bounds = _planck_split(low, min(high, upper_cm), fractions[g])
        pieces += [(float(bounds[i]), float(bounds[i + 1]), int(g[i])) for i in range(g.size)]
    pieces.sort(key=lambda piece: 1.0e-2 / piece[1])
    long_edges = [1.0e-2 / piece[1] for piece in pieces] + [1.0e-2 / pieces[-1][0]]
    boundary = 1.0e-2 / upper_cm
    shortwave = [float(e) for e in edges if e < boundary and not np.isclose(e, boundary)]
    new_grid = SpectralGrid(
        grid_id=f"{grid.grid_id}+rrtmg-lw",
        coordinate=grid.coordinate,
        coordinate_units=grid.coordinate_units,
        edges=np.asarray([*shortwave, *long_edges]),
        integration_convention=grid.integration_convention,
        density_units=grid.density_units,
        resolution_identifier=f"{grid.resolution_identifier}+rrtmg-lw",
    )
    gpoint = np.array([-1] * len(shortwave) + [piece[2] for piece in pieces])
    if gpoint.size != new_grid.n_bins:
        raise ConfigurationError("internal: RRTMG's g-points do not match the grid")
    return new_grid, gpoint


# ---- the model's state in RRTMG's terms -------------------------------------------

#: Molecules per cm^2 in one cm-atm (Loschmidt's number at 273.15 K, 1 atm).
_LOSCHMIDT_PER_CM3 = 2.6867811e19


@dataclass(frozen=True)
class RRTMGLongwave:
    """RRTMG on a grid: each bin's g-point (-1 outside the longwave) and the
    well-mixed gases RRTMG needs that the model does not carry as state.

    The defaults are 1850 (Meinshausen et al. 2017, CMIP6 historical), to match
    the column's 284.7 ppmv CO2; the CFCs and CO are zero.
    """

    gpoint: np.ndarray
    co2_ppmv: float
    ch4_ppb: float = 808.2
    n2o_ppb: float = 273.0
    o2_fraction: float = 0.209

    def carried(self, home: np.ndarray) -> RRTMGLongwave:
        """The same spec on a grid rebuilt elsewhere: each new bin takes its home bin's."""
        return RRTMGLongwave(
            gpoint=self.gpoint[home],
            co2_ppmv=self.co2_ppmv,
            ch4_ppb=self.ch4_ppb,
            n2o_ppb=self.n2o_ppb,
            o2_fraction=self.o2_fraction,
        )


def rrtmg_optical_depth(
    spec: RRTMGLongwave,
    *,
    vapour: np.ndarray,
    layer_mass: np.ndarray,
    pressure: np.ndarray,
    temperature: np.ndarray,
    ozone_path: np.ndarray,
) -> np.ndarray:
    """Vertical gaseous optical depth on the grid, shaped (layers, bins).

    ``vapour`` and ``layer_mass`` in kg m^-2 (the layer mass includes the
    vapour), ``pressure`` in Pa, ``temperature`` in K, ``ozone_path`` in cm-atm.
    The dry-air column is the dry mass, the same quantity RRTMG's driver forms
    from the pressure thickness.
    """
    w = np.asarray(vapour, dtype=float)
    m = np.asarray(layer_mass, dtype=float)
    if np.any(w < 0.0):
        raise ConfigurationError("a vapour path cannot be negative")
    dry = m - w
    if np.any(dry <= 0.0):
        raise ConfigurationError("a layer's vapour cannot exceed its mass")
    coldry = dry * 0.1 / _AMD * _AVOGADRO  # kg m^-2 -> g cm^-2 -> molecules cm^-2
    col = RRTMGColumn.from_mixing_ratios(
        pressure_hpa=np.asarray(pressure, dtype=float) / 100.0,
        temperature=np.asarray(temperature, dtype=float),
        coldry=coldry,
        vmr={
            "h2o": (w / dry) * _AMD / _AMW,
            "co2": spec.co2_ppmv * 1.0e-6,
            "o3": np.asarray(ozone_path, dtype=float) * _LOSCHMIDT_PER_CM3 / coldry,
            "ch4": spec.ch4_ppb * 1.0e-9,
            "n2o": spec.n2o_ppb * 1.0e-9,
            "o2": spec.o2_fraction,
        },
    )
    taug, _, _ = rrtmg_gas_optics(col)
    tau = np.zeros((m.size, spec.gpoint.size))
    lw = spec.gpoint >= 0
    tau[:, lw] = taug[:, spec.gpoint[lw]]
    return tau
