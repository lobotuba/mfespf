"""Ozone's solar absorption, and the measured solar spectrum it needs.

Each checked against something independent of the code: the measured spectrum
against its own published integral and against Chou and Suarez's band
fractions, which came from a different spectrum; the ozone bands against the
source table; the ozone column against the integral of the profile it was
built from; and the absorption against the physics -- ultraviolet taken up
high, visible barely touched.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError
from mfespf.physics.solar_absorption import (
    AFGL_US_STANDARD_OZONE_HPA_PPMV,
    OZONE_BANDS,
    OZONE_PROFILE_HPA_PPMV,
    ozone_bands,
    ozone_layer_path,
    ozone_optical_depth,
    solar_vapour_bands,
)
from mfespf.physics.solar_spectrum import HSRS_TOTAL_W_M2, hsrs_band_irradiance
from mfespf.run.column_configurations import column_grid

pytestmark = pytest.mark.stage(4)

G = 9.80665
INTERFACES = np.linspace(101325.0, 1000.0, 31)
PRESSURE = 0.5 * (INTERFACES[:-1] + INTERFACES[1:])
MASS = (INTERFACES[:-1] - INTERFACES[1:]) / G

#: Chou and Suarez (1999) Table 3, S/S0 for bands 2 to 9, from a different
#: spectrum than HSRS. Band 1 starts at 175 nm, below HSRS's full-resolution
#: product's edge, and is checked in the same way through the full spectrum.
CHOU_SUAREZ_FRACTIONS = (
    ((0.175e-6, 0.225e-6), 0.00057),
    ((0.245e-6, 0.260e-6), 0.00083),
    ((0.280e-6, 0.295e-6), 0.00417),
    ((0.295e-6, 0.310e-6), 0.00600),
    ((0.310e-6, 0.320e-6), 0.00556),
    ((0.320e-6, 0.400e-6), 0.05913),
    ((0.400e-6, 0.700e-6), 0.39081),
    ((0.700e-6, 1.220e-6), 0.32055),
    ((1.220e-6, 2.270e-6), 0.16536),
)


@pytest.mark.req("SOL-001", "SPC-001")
def test_the_measured_spectrum_integrates_to_its_total():
    whole = hsrs_band_irradiance(np.array([100.0e-9, 250.0e-6]))
    assert float(whole.sum()) == pytest.approx(HSRS_TOTAL_W_M2, rel=1e-12)
    assert HSRS_TOTAL_W_M2 == pytest.approx(1361.49, abs=0.01)
    # Nothing is invented outside the measured support.
    assert float(hsrs_band_irradiance(np.array([50.0e-9, 110.0e-9]))[0]) == 0.0


@pytest.mark.req("SOL-001", "SPC-001", "VAL-001")
def test_the_measured_spectrum_reproduces_chou_and_suarez_band_fractions():
    """The black body misses these by half in the ultraviolet; the measurement does not."""
    for (lower, upper), expected in CHOU_SUAREZ_FRACTIONS:
        fraction = float(hsrs_band_irradiance(np.array([lower, upper]))[0]) / HSRS_TOTAL_W_M2
        # Two independent spectra, different decades and instruments: a few
        # percent in the visible and infrared, up to 15 in the narrow UV bands.
        tolerance = 0.03 if lower >= 0.32e-6 else 0.15
        assert fraction == pytest.approx(expected, rel=tolerance), (lower, upper)


@pytest.mark.req("SPC-001", "EQR-002")
def test_the_ozone_bands_carry_the_source_coefficients():
    grid, _ = solar_vapour_bands(column_grid(16), measured_spectrum=True)
    rebuilt, k = ozone_bands(grid)
    edges = np.asarray(rebuilt.edges)
    centres = 0.5 * (edges[:-1] + edges[1:])
    for lower, upper, coefficient in OZONE_BANDS:
        inside = (centres > lower) & (centres < upper)
        assert inside.sum() == 1, (lower, upper)
        assert k[inside][0] == coefficient
    # Nothing absorbs below 0.175 microns or above 0.7.
    assert np.all(k[centres < 0.175e-6] == 0.0)
    assert np.all(k[centres > 0.7e-6] == 0.0)
    # The vapour's sub-bands above 0.7 microns are kept edge for edge.
    kept = np.asarray(grid.edges)[np.asarray(grid.edges) > 0.7e-6]
    assert np.all(np.isin(kept, edges))


@pytest.mark.req("EQR-002", "VAL-001")
def test_the_ozone_column_is_the_profile_integral():
    """Built layer by layer, the column converges to integrating the profile.

    Each layer takes the mixing ratio at its centre, as the protocol froze it.
    Ozone rises steeply with height, so a thick layer's centre value is below
    its mean: on the 30-layer grid the column is 3.0 percent short of the
    integral. That is the method's discretization error, and it is stated
    here and in the results rather than hidden by a looser tolerance. On a
    fine grid the construction must agree.
    """
    fine = np.geomspace(101325.0, 100.0, 20000)
    p = np.array([row[0] for row in OZONE_PROFILE_HPA_PPMV])[::-1] * 100.0
    v = np.array([row[1] for row in OZONE_PROFILE_HPA_PPMV])[::-1] * 1.0e-6
    mixing = np.interp(np.log(fine), np.log(p), v)
    # Molecules per kilogram of air, R_d / k_B; a cm-atm is n0 x 1 cm.
    column_kg = -np.trapezoid(mixing, fine) / G
    molecules = column_kg * 287.047 / 1.380649e-23
    cm_atm = 101325.0 / (1.380649e-23 * 273.15) * 0.01
    integral = molecules / cm_atm

    many = np.linspace(101325.0, 1000.0, 1001)
    path_fine = ozone_layer_path(0.5 * (many[:-1] + many[1:]), many, (many[:-1] - many[1:]) / G)
    assert float(path_fine.sum()) == pytest.approx(integral, rel=2e-3)

    path = ozone_layer_path(PRESSURE, INTERFACES, MASS)
    shortfall = 1.0 - float(path.sum()) / integral
    assert shortfall == pytest.approx(0.030, abs=0.002)
    # About 290 dobson units, the observed global mean's size.
    assert 0.27 < float(path.sum()) < 0.31


@pytest.mark.req("EQR-002", "VAL-001")
def test_ozone_absorbs_the_ultraviolet_and_little_else():
    grid, _ = solar_vapour_bands(column_grid(16), measured_spectrum=True)
    rebuilt, k = ozone_bands(grid)
    tau = ozone_optical_depth(k, ozone_layer_path(PRESSURE, INTERFACES, MASS)).sum(axis=0)
    centres = np.asarray(rebuilt.centers)
    hartley = (centres > 0.225e-6) & (centres < 0.28e-6)
    visible = (centres > 0.4e-6) & (centres < 0.7e-6)
    # Opaque in the Hartley band, and the sunlight there is taken above 100 hPa.
    assert np.all(tau[hartley] > 50.0)
    above_100 = ozone_optical_depth(k, ozone_layer_path(PRESSURE, INTERFACES, MASS))[
        PRESSURE < 1.0e4
    ].sum(axis=0)
    assert np.all(above_100[hartley] > 30.0)
    # Nearly transparent in the Chappuis band: a couple of percent.
    assert 0.01 < float(tau[visible][0]) < 0.03


@pytest.mark.req("DAT-007", "VAL-001")
def test_ozone_above_1hpa_is_refused_unless_declared_and_never_held():
    """upper-atmosphere-radiative-equilibrium.md: SWOOSH ends at 1 hPa. Above it
    no value is held: a layer there is refused, unless the AFGL extension is
    declared, which is itself refused beyond its own table."""
    centre = np.array([5000.0, 150.0, 50.0])
    interfaces = np.array([10000.0, 200.0, 100.0, 10.0])
    mass = (interfaces[:-1] - interfaces[1:]) / G
    with pytest.raises(ConfigurationError, match="not held"):
        ozone_layer_path(centre, interfaces, mass)
    path = ozone_layer_path(centre, interfaces, mass, "afgl_us_standard")
    assert np.all(path > 0.0)
    beyond = np.array([5000.0, 150.0, 1.0])
    with pytest.raises(ConfigurationError, match="AFGL ozone table ends"):
        ozone_layer_path(beyond, np.array([10000.0, 200.0, 2.0, 0.5]), mass, "afgl_us_standard")


@pytest.mark.req("VAL-001")
def test_the_two_ozone_records_meet_at_1hpa_each_exact_on_its_own_side():
    """Below 1 hPa the profile is SWOOSH's, at 1 hPa exactly SWOOSH's last
    level; above it, AFGL's log-pressure interpolation. The step between them
    is the records' disagreement, reported rather than smoothed."""
    from mfespf.physics.solar_absorption import _ozone_mixing_ratio

    at_1hpa = _ozone_mixing_ratio(np.array([100.0]), "afgl_us_standard")[0]
    assert at_1hpa == pytest.approx(OZONE_PROFILE_HPA_PPMV[-1][1] * 1e-6, rel=1e-12)
    just_above = _ozone_mixing_ratio(np.array([100.0 * (1.0 - 1e-12)]), "afgl_us_standard")[0]
    (p1, v1), (p2, v2) = AFGL_US_STANDARD_OZONE_HPA_PPMV[1], AFGL_US_STANDARD_OZONE_HPA_PPMV[2]
    afgl_1hpa = v1 + (v2 - v1) * np.log(p1 / 1.0) / np.log(p1 / p2)
    assert just_above == pytest.approx(afgl_1hpa * 1e-6, rel=1e-9)
    assert afgl_1hpa / OZONE_PROFILE_HPA_PPMV[-1][1] == pytest.approx(1.20, abs=0.01)


@pytest.mark.req("VAL-001")
def test_the_afgl_column_above_1hpa_is_its_table_integral():
    """Built layer by layer on a fine grid from 1 to 0.024 hPa, the ozone above
    1 hPa converges to the integral of the AFGL table it comes from."""
    p = np.array([row[0] for row in AFGL_US_STANDARD_OZONE_HPA_PPMV])[::-1] * 100.0
    v = np.array([row[1] for row in AFGL_US_STANDARD_OZONE_HPA_PPMV])[::-1] * 1.0e-6
    fine = np.geomspace(100.0, p[0], 20000)
    mixing = np.interp(np.log(fine), np.log(p), v)
    column_kg = -np.trapezoid(mixing, fine) / G
    integral = column_kg * 287.047 / 1.380649e-23 / (101325.0 / (1.380649e-23 * 273.15) * 0.01)
    many = np.geomspace(100.0, p[0], 801)
    path = ozone_layer_path(
        np.sqrt(many[:-1] * many[1:]), many, (many[:-1] - many[1:]) / G, "afgl_us_standard"
    )
    assert float(path.sum()) == pytest.approx(integral, rel=2e-3)
    # About a dobson unit, as the code's own note says of what it used to omit.
    assert 0.5e-3 < integral < 2.0e-3
