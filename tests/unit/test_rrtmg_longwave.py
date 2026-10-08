"""RRTMG's longwave gas optics, checked against AER's own published output.

rrtmg-longwave.md P0: the port of setcoef and taumol, driven through a port of
RRTMG's clear-sky ``rtrnmr`` that lives here, in the check, must reproduce
AER's ``output_rrtm_*-clr`` for its four standard atmospheres. This is a
transcription check: the reference is RRTMG itself, not this model's solver.
The input reader is ``readprof`` (rrtmg_lw.1col.f90); the transfer, its lookup
tables and the diffusivity are ``rrtmg_lw_rtrnmr.f90`` and ``rrtmg_lw_init.f90``
of AER-RC/RRTMG_LW 958edb1a.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

from mfespf.physics.rrtmg_longwave import (
    BAND_WAVENUMBER_CM,
    G_POINTS_PER_BAND,
    RRTMGColumn,
    rrtmg_gas_optics,
    rrtmg_tables,
)

pytestmark = pytest.mark.stage(6)

ROOT = Path(__file__).parents[1] / "data" / "aer_rrtmg_lw"
NGS = np.concatenate(([0], np.cumsum(G_POINTS_PER_BAND)))
AMD, AMW, GRAV = 28.966, 18.016, 9.8066


def read_input(path):
    lines = path.read_text().splitlines()
    i = next(k for k, s in enumerate(lines) if s.startswith("$"))
    tbound = float(lines[i + 2].split()[0])
    _iform, nlay, nmol = (int(x) for x in lines[i + 3].split()[:3])
    pavel, tavel, pz, tz, wkl, wbr = [], [], [], [], [], []
    k = i + 4
    for lay in range(nlay):
        s = lines[k].ljust(85)
        # form1/form2(iform=1): g15.7,g10.4,g10.4,a3,i2, then 1x,2(g7.2,g8.3,g7.2) or 23x,(...)
        pavel.append(float(s[0:15]))
        tavel.append(float(s[15:25]))
        if lay == 0:
            pz += [float(s[48:56]), float(s[70:78])]
            tz += [float(s[56:63]), float(s[78:85])]
        else:
            pz.append(float(s[70:78]))
            tz.append(float(s[78:85]))
        amounts = [float(x) for x in lines[k + 1].split()]
        wkl.append(amounts[:7])
        wbr.append(amounts[7])
        k += 2 if nmol <= 7 else 3
    wkl = np.array(wkl).T  # (7, nlay)
    wbr = np.array(wbr)
    imix = 1 if np.all(wkl[:nmol, 0] <= 1.0) else 0
    summol = np.sum(wkl[1:nmol], axis=0)
    if imix == 1:
        coldry = wbr / (1.0 - summol)
        wkl = coldry * wkl
    else:
        coldry = wbr + summol
    amttl = np.sum(coldry + wkl[0])
    wvttl = np.sum(wkl[0])
    pwvcm = (AMW * wvttl) / (AMD * amttl) * (1.0e3 * pz[0]) / (1.0e2 * GRAV)
    col = RRTMGColumn(
        pressure_hpa=np.array(pavel),
        temperature=np.array(tavel),
        coldry=coldry,
        wkl=wkl,
        wbroad=wbr,
        wx=np.zeros((4, nlay)),
    )
    return col, np.array(pz), np.array(tz), tbound, pwvcm


def read_output(path):
    blocks, cur = {}, None
    for s in path.read_text().splitlines():
        if "Wavenumbers" in s:
            a, b = re.findall(r"([\d.]+)\s*-\s*([\d.]+)\s*cm", s)[0]
            cur = (float(a), float(b))
            blocks[cur] = []
            continue
        tok = s.split()
        if cur and len(tok) == 6 and tok[0].isdigit():
            blocks[cur].append([float(x) for x in tok])
    return {k: np.array(sorted(v, key=lambda r: r[0])) for k, v in blocks.items()}  # by level


# ---- rrtmg_lw_init lookup tables -----------------------------------------------
NTBL, TBLINT, PADE = 10000, 10000.0, 0.278
BPADE = 1.0 / PADE
tau_tbl = np.zeros(NTBL + 1)
exp_tbl = np.zeros(NTBL + 1)
tfn_tbl = np.zeros(NTBL + 1)
tau_tbl[NTBL], exp_tbl[0], exp_tbl[NTBL], tfn_tbl[NTBL] = 1.0e10, 1.0, 1.0e-20, 1.0
itr = np.arange(1, NTBL)
tfn = itr / NTBL
tau_tbl[1:NTBL] = BPADE * tfn / (1.0 - tfn)
exp_tbl[1:NTBL] = np.maximum(np.exp(-tau_tbl[1:NTBL]), 1.0e-20)
tt, ee = tau_tbl[1:NTBL], exp_tbl[1:NTBL]
tfn_tbl[1:NTBL] = np.where(tt < 0.06, tt / 6.0, 1.0 - 2.0 * ((1.0 / tt) - (ee / (1.0 - ee))))

A0 = [
    1.66,
    1.55,
    1.58,
    1.66,
    1.54,
    1.454,
    1.89,
    1.33,
    1.668,
    1.66,
    1.66,
    1.66,
    1.66,
    1.66,
    1.66,
    1.66,
]
A1 = [0.00, 0.25, 0.22, 0.00, 0.13, 0.446, -0.10, 0.40, -0.006, 0, 0, 0, 0, 0, 0, 0]
A2 = [0.00, -12.0, -11.7, 0.00, -0.72, -0.243, 0.19, -0.062, 0.414, 0, 0, 0, 0, 0, 0, 0]


def planck(table, temp):
    ind = np.clip(np.trunc(temp - 159.0).astype(int), 1, 180)
    frac = temp - 159.0 - ind
    return table[ind - 1] + frac * (table[ind] - table[ind - 1])


def rtrnmr(col, pz, tz, tbound, pwvcm, taug, fracs, bands, band16_total=True):
    t = rrtmg_tables()
    nlay = col.pressure_hpa.size
    up, dn = np.zeros(nlay + 1), np.zeros(nlay + 1)
    for ib in bands:
        n = ib + 1
        if n == 1 or n == 4 or n >= 10:
            secdiff = 1.66
        else:
            secdiff = min(max(A0[ib] + A1[ib] * np.exp(A2[ib] * pwvcm), 1.50), 1.80)
        table = t["totplnk"][:, ib] if (n < 16 or band16_total) else t["totplk16"]
        blay_all = planck(table, col.temperature)
        blev = planck(table, tz)
        bnd = planck(table, np.array(tbound))
        urad, drad = np.zeros(nlay + 1), np.zeros(nlay + 1)
        for igc in range(NGS[ib], NGS[ib + 1]):
            atrans, bbugas = np.zeros(nlay), np.zeros(nlay)
            radld = 0.0
            for lev in range(nlay - 1, -1, -1):  # Fortran lev = nlayers..1 is lev-1 here
                plfrac = fracs[lev, igc]
                blay = blay_all[lev]
                dup = blev[lev + 1] - blay
                ddn = blev[lev] - blay
                od = max(secdiff * taug[lev, igc], 0.0)
                if od <= 0.06:
                    atrans[lev] = od - 0.5 * od * od
                    odr = 0.166667 * od
                    bbd = plfrac * (blay + ddn * odr)
                    bbugas[lev] = plfrac * (blay + dup * odr)
                else:
                    it = int(TBLINT * (od / (BPADE + od)) + 0.5)
                    atrans[lev] = 1.0 - exp_tbl[it]
                    tf = tfn_tbl[it]
                    bbd = plfrac * (blay + tf * ddn)
                    bbugas[lev] = plfrac * (blay + tf * dup)
                radld = radld + (bbd - radld) * atrans[lev]
                drad[lev] += radld
            radlu = fracs[0, igc] * bnd  # semiss = 1, no reflection
            urad[0] += radlu
            for lev in range(nlay):
                radlu = radlu + (bbugas[lev] - radlu) * atrans[lev]
                urad[lev + 1] += radlu
        delwave = BAND_WAVENUMBER_CM[ib][1] - BAND_WAVENUMBER_CM[ib][0]
        up += urad * 0.5 * delwave
        dn += drad * 0.5 * delwave
    fluxfac = np.pi * 2.0e4
    return up * fluxfac, dn * fluxfac


@pytest.mark.req("EQR-002", "VAL-001")
@pytest.mark.parametrize("atmosphere", ["MLS", "SAW", "TROP", "MLW"])
def test_the_port_reproduces_aer_s_published_fluxes(atmosphere):
    """P0: every level within 0.1 W m^-2, and every published band within 0.05."""
    col, pz, tz, tbound, pwvcm = read_input(ROOT / f"input_rrtm_{atmosphere}-clr")
    ref = read_output(ROOT / f"output_rrtm_{atmosphere}-clr")
    taug, fracs, _ = rrtmg_gas_optics(col)
    up, dn = rtrnmr(col, pz, tz, tbound, pwvcm, taug, fracs, range(16))
    total = ref[(10.0, 3250.0)]
    np.testing.assert_allclose(up, total[:, 2], atol=0.1)
    np.testing.assert_allclose(dn, total[:, 3], atol=0.1)
    for ib in range(16):
        key = BAND_WAVENUMBER_CM[ib]
        if key in ref:
            u, _ = rtrnmr(col, pz, tz, tbound, pwvcm, taug, fracs, [ib], band16_total=(ib < 15))
            assert u[-1] == pytest.approx(ref[key][-1, 2], abs=0.05), key


@pytest.mark.req("EQR-002", "VAL-001")
def test_more_water_absorbs_more_in_the_water_bands():
    """A limiting case independent of AER's output: doubling the vapour raises
    the optical depth of every g-point in the pure-water bands (1, 2, 10)."""
    col, *_ = read_input(ROOT / "input_rrtm_MLS-clr")
    wetter = RRTMGColumn(
        pressure_hpa=col.pressure_hpa,
        temperature=col.temperature,
        coldry=col.coldry,
        wkl=np.vstack([2.0 * col.wkl[:1], col.wkl[1:]]),
        wbroad=col.wbroad,
        wx=col.wx,
    )
    base, _, _ = rrtmg_gas_optics(col)
    more, _, _ = rrtmg_gas_optics(wetter)
    ngs = np.concatenate(([0], np.cumsum(G_POINTS_PER_BAND)))
    for band in (1, 2, 10):
        s = slice(ngs[band - 1], ngs[band])
        assert np.all(more[:, s] > base[:, s]), band


def _rrtmg_column(**overrides):
    from mfespf.run.column_configurations import column_grid, reference_column

    options = dict(
        grid=column_grid(16),
        n_layers=10,
        master_dt=3600.0,
        measured_solar_spectrum=True,
        ozone=True,
        solar_vapour_absorption=True,
        rrtmg_longwave=True,
    )
    options.update(overrides)
    return reference_column(**options)


def _radiation(kernel):
    return next(s.module for s in kernel.schedule.ordered() if s.module_id == "radiation")


@pytest.mark.req("SPC-001", "EQR-002")
def test_the_radiation_module_solves_with_rrtmg_s_optics():
    """The module's thermal flux is a direct solve with rrtmg_optical_depth on its inputs."""
    from mfespf.physics.radiative_transfer import solve_column_radiation
    from mfespf.physics.rrtmg_longwave import rrtmg_optical_depth

    kernel = _rrtmg_column()
    kernel.run(1)
    rad = _radiation(kernel)

    def value(name, units):
        return np.asarray(rad.input_field(name).require_units(units).values, dtype=float)

    temperature = value("atmosphere_layer_temperature", "K")
    mass = value("atmosphere_layer_mass", "kg m^-2")
    tau = rrtmg_optical_depth(
        rad.rrtmg,
        vapour=value("atmosphere_water_vapour", "kg m^-2"),
        layer_mass=mass,
        pressure=value("atmosphere_layer_pressure", "Pa"),
        temperature=temperature,
        ozone_path=rad.ozone_path,
    )
    direct = solve_column_radiation(
        grid=rad.grid,
        incident_spectral_irradiance=value("toa_spectral_irradiance", "W m^-2 m^-1"),
        layer_temperature=temperature,
        layer_mass=mass,
        absorption_coefficient=rad.absorption_coefficient,
        surface_temperature=float(value("surface_temperature", "K")),
        surface_albedo=value("surface_spectral_albedo", "1"),
        diffusivity_factor=rad.diffusivity_factor,
        absorption_optical_depth=tau,
    )
    thermal = rad.rrtmg.gpoint >= 0
    np.testing.assert_allclose(
        np.asarray(rad.clear_sky.upward[-1])[thermal],
        np.asarray(direct.upward[-1])[thermal],
        rtol=1e-9,
    )
    # The vapour's own optical depth is solar only: RRTMG carries its longwave.
    assert np.all(rad.vapour_optical_depth[:, thermal] == 0.0)


@pytest.mark.req("SPC-001", "EQR-002")
def test_rrtmg_refuses_to_join_the_schemes_it_replaces_or_to_run_without_its_inputs():
    from mfespf.core.errors import ConfigurationError

    for clash in ({"co2_band": True}, {"chou_water_vapour": True, "co2_band": True}):
        with pytest.raises(ConfigurationError, match="replaces"):
            _rrtmg_column(**clash)
    with pytest.raises(ConfigurationError, match="ozone"):
        _rrtmg_column(ozone=False)
    with pytest.raises(ConfigurationError, match="radiative_vapour"):
        _rrtmg_column(radiative_vapour=False)


@pytest.mark.req("SPC-001", "EQR-002")
def test_rrtmg_replaces_the_fitted_envelope(monkeypatch):
    from mfespf.physics import greenhouse

    def refuse(*args, **kwargs):
        raise AssertionError("the fitted vapour envelope was evaluated")

    monkeypatch.setattr(greenhouse, "vapour_absorption_optical_depth", refuse)
    kernel = _rrtmg_column()
    kernel.run(2)
    assert _radiation(kernel).rrtmg is not None
