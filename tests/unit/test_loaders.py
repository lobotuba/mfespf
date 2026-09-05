"""Readers for the declared external records.

Requirements DAT-004, DAT-007, DAT-009, INF-002, INF-007, TIM-001, UNC-002.

Every test here builds its own fixture rather than reading `data/raw/`, so the
parsers are checked without a network and without anyone having run
`tools/fetch_data.py`. The integration suite covers the real files and skips
when they are absent.

The fixtures are short excerpts in each provider's actual layout, which is the
point: the formats agree on nothing, and a parser that only ever saw a tidy
version of its input has not been tested.
"""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pytest

from mfespf.core.errors import ConfigurationError, ProvenanceError
from mfespf.data.loaders import (
    EPOCH,
    load_berkeley,
    load_hadcrut5,
    load_noaa_ohc,
    load_nrltsi2,
    seconds_since_epoch,
)
from mfespf.data.records import SamplingConvention, UncertaintyKind

pytestmark = pytest.mark.stage(6)

YEAR = 365.25 * 86400.0

HADCRUT5 = """Time,Anomaly (deg C),Lower confidence limit (2.5%),Upper confidence limit (97.5%)
1850-01,-0.67456436,-0.98177195,-0.3673568
1850-02,-0.333416,-0.700901,0.034069024
1850-12,-0.2,-0.4,0.0
1851-01,-0.1,-0.3,0.1
"""

BERKELEY = """% Berkeley Earth, excerpt.
% Year, Annual Anomaly, Annual Unc., Five-year Anomaly, Five-year Unc., Annual Anomaly, Annual Unc., ...

  1850      -0.436        0.170              NaN             NaN            -0.419        0.156              NaN             NaN
  1851      -0.321        0.259           -0.400           0.100           -0.319        0.225           -0.390        0.090
  1852      -0.300        0.200           -0.380           0.090           -0.290        0.180           -0.370        0.080
"""

NOAA_OHC = """    YEAR      WO    WOse      NH    NHse      SH    SHse
1955.500  -3.201   1.700  -1.439   0.937  -1.762   0.771
1956.500  -2.608   0.719  -1.843   0.385  -0.765   0.533
"""

NRLTSI2 = """time (days since 1610-01-01),irradiance (W/m^2),uncertainty (W/m^2)
181.0,1360.5418701171875,0.4221000075340271
546.0,1360.5716552734375,0.4361000061035156
"""


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


# --- the shared clock ------------------------------------------------------


@pytest.mark.req("TIM-001")
def test_every_loader_puts_its_record_on_one_declared_epoch(tmp_path):
    """Four formats, four epochs in the files, one clock afterwards.

    This is what makes two records comparable without anyone doing calendar
    arithmetic by hand, which is the arithmetic that silently puts a series a
    year out.
    """
    assert EPOCH == datetime(1850, 1, 1, tzinfo=UTC)
    assert seconds_since_epoch(EPOCH) == 0.0

    hadcrut = load_hadcrut5(_write(tmp_path, "h.csv", HADCRUT5))
    berkeley = load_berkeley(_write(tmp_path, "b.txt", BERKELEY))
    noaa = load_noaa_ohc(_write(tmp_path, "n.dat", NOAA_OHC))
    nrl = load_nrltsi2(_write(tmp_path, "t.csv", NRLTSI2))

    # January 1850 is the first half-month of the epoch; mid-1850 is half a year
    # in; mid-1955 is a century later; mid-1610 is 240 years before it.
    assert hadcrut.times[0] == pytest.approx(15.5 * 86400.0)
    assert berkeley.times[0] == pytest.approx(0.5 * YEAR, rel=1e-3)
    assert noaa.times[0] == pytest.approx(105.5 * YEAR, rel=1e-4)
    assert nrl.times[0] < 0.0  # 1610 predates the epoch, and that is allowed
    assert nrl.times[0] == pytest.approx(-239.5 * YEAR, rel=1e-3)


@pytest.mark.req("TIM-001", "DAT-003")
def test_samples_land_at_interval_midpoints_not_starts(tmp_path):
    """An interval mean sits at the middle of its interval.

    ``resample`` builds its cells centred on the sample times, so a series
    timestamped at month starts would be resampled half a month early -- a small
    error, and one that survives every check that does not look at a calendar.
    Checked against midpoints computed here rather than against the loader's own
    arithmetic.
    """
    series = load_hadcrut5(_write(tmp_path, "h.csv", HADCRUT5))
    for index, (year, month) in enumerate([(1850, 1), (1850, 2), (1850, 12), (1851, 1)]):
        start = datetime(year, month, 1, tzinfo=UTC)
        end = (
            datetime(year + 1, 1, 1, tzinfo=UTC)
            if month == 12
            else datetime(year, month + 1, 1, tzinfo=UTC)
        )
        expected = seconds_since_epoch(start) + 0.5 * (end - start).total_seconds()
        assert series.times[index] == pytest.approx(expected)
        assert series.times[index] > seconds_since_epoch(start)


# --- HadCRUT5 ---------------------------------------------------------------


@pytest.mark.req("INF-007", "UNC-002")
def test_hadcrut5_anomaly_is_kelvin_and_its_interval_is_not_one_sigma(tmp_path):
    """Two corrections the file's own headings would lead you into.

    Its column says "deg C", but the value is a *difference*, and a difference of
    one degree Celsius is one kelvin exactly -- declaring it ``degC`` would make
    the unit registry apply the 273.15 offset at the first conversion. And the
    published interval is a 95 percent range, which is about two sigma;
    combining it in quadrature as though it were one would halve whatever it
    entered.
    """
    series = load_hadcrut5(_write(tmp_path, "h.csv", HADCRUT5))

    assert series.units == "K"
    assert series.uncertainty_kind is UncertaintyKind.HALF_WIDTH_95
    assert not series.uncertainty_kind.is_one_sigma
    assert "not an absolute temperature" in series.note
    assert series.values[0] == pytest.approx(-0.67456436)
    assert series.uncertainty[0] == pytest.approx(0.5 * (-0.3673568 - -0.98177195))
    assert series.sampling is SamplingConvention.INTERVAL_MEAN


@pytest.mark.req("TIM-001")
def test_hadcrut5_handles_a_december_to_january_rollover(tmp_path):
    series = load_hadcrut5(_write(tmp_path, "h.csv", HADCRUT5))
    assert series.times.size == 4
    assert np.all(np.diff(series.times) > 0.0)
    # December 1850 and January 1851 are one month apart, not eleven.
    assert (series.times[3] - series.times[2]) == pytest.approx(31.0 * 86400.0, rel=0.05)


# --- Berkeley Earth ---------------------------------------------------------


@pytest.mark.req("DAT-007")
def test_berkeley_drops_the_rows_it_has_no_value_for(tmp_path):
    """A gap is a shorter series, never a sentinel.

    The file writes NaN in its five-year columns until it has five years. Those
    rows are dropped from a five-year series rather than carried, because a
    sentinel is a number something downstream will multiply.
    """
    annual = load_berkeley(_write(tmp_path, "b.txt", BERKELEY))
    assert annual.times.size == 3  # the annual columns are complete
    assert np.all(np.isfinite(annual.values))


@pytest.mark.req("INF-002", "MOD-002")
def test_berkeley_sea_ice_choice_is_definitional_and_recorded(tmp_path):
    """The two series are different definitions, not two estimates of one thing.

    Air temperature above sea ice and water temperature below it diverge in the
    recent record because water is bounded below by the freezing point. A
    comparison against the wrong one would look like a model error, so the
    choice is a parameter and it goes in the note.
    """
    path = _write(tmp_path, "b.txt", BERKELEY)
    air = load_berkeley(path, sea_ice="air")
    water = load_berkeley(path, sea_ice="water")

    assert air.values[0] == pytest.approx(-0.436)
    assert water.values[0] == pytest.approx(-0.419)
    assert "air temperature used over sea ice" in air.note
    assert "water temperature used over sea ice" in water.note
    with pytest.raises(ConfigurationError, match="sea_ice must be"):
        load_berkeley(path, sea_ice="either")


# --- NOAA ocean heat content ------------------------------------------------


@pytest.mark.req("DAT-004")
def test_noaa_ohc_converts_its_unit_at_the_boundary(tmp_path):
    """The file is in 10^22 J, which is ten zettajoules and not one.

    The factor is applied here, at a declared boundary, rather than left in a
    header for a reader to notice. Getting it wrong is an order of magnitude,
    and nothing downstream would flag it.
    """
    series = load_noaa_ohc(_write(tmp_path, "n.dat", NOAA_OHC))
    assert series.units == "J"
    assert series.values[0] == pytest.approx(-3.201e22)
    assert series.uncertainty[0] == pytest.approx(1.700e22)
    assert series.uncertainty_kind is UncertaintyKind.STANDARD_ERROR
    assert series.uncertainty_kind.is_one_sigma
    assert "10^22 J" in series.note


@pytest.mark.req("SPA-002", "MOD-002")
def test_noaa_ohc_selects_a_named_basin_and_refuses_an_unknown_one(tmp_path):
    path = _write(tmp_path, "n.dat", NOAA_OHC)
    assert load_noaa_ohc(path, basin="NH").values[0] == pytest.approx(-1.439e22)
    assert load_noaa_ohc(path, basin="SH").values[0] == pytest.approx(-1.762e22)
    with pytest.raises(ConfigurationError, match="no 'atlantic' column"):
        load_noaa_ohc(path, basin="atlantic")


# --- NRLTSI2 ----------------------------------------------------------------


@pytest.mark.req("DAT-009", "SYS-006")
def test_nrltsi2_refuses_a_file_that_does_not_declare_the_epoch_it_assumes(tmp_path):
    """Reading it anyway would put the series out by whatever the difference is.

    The epoch is in the column heading and nowhere else; a provider changing it
    would produce a file that parses cleanly and is silently wrong by centuries.
    """
    good = load_nrltsi2(_write(tmp_path, "t.csv", NRLTSI2))
    assert good.values[0] == pytest.approx(1360.5418701171875)
    assert good.uncertainty_kind is UncertaintyKind.STANDARD_UNCERTAINTY

    bad = NRLTSI2.replace("days since 1610-01-01", "days since 1900-01-01")
    with pytest.raises(ProvenanceError, match="does not declare the epoch"):
        load_nrltsi2(_write(tmp_path, "bad.csv", bad))


@pytest.mark.req("UNC-002")
def test_a_series_carrying_uncertainties_must_say_what_they_are(tmp_path):
    """The gap real data exposed: four providers, three different meanings.

    A standard error and a 95 percent half-width differ by a factor of two, and
    nothing in the numbers distinguishes them.
    """
    kinds = {
        load_hadcrut5(_write(tmp_path, "h.csv", HADCRUT5)).uncertainty_kind,
        load_noaa_ohc(_write(tmp_path, "n.dat", NOAA_OHC)).uncertainty_kind,
        load_nrltsi2(_write(tmp_path, "t.csv", NRLTSI2)).uncertainty_kind,
    }
    assert len(kinds) == 3
    assert UncertaintyKind.UNDECLARED not in kinds
