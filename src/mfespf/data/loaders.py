"""Readers for the declared external records.

Requirements DAT-004, DAT-007, DAT-009, INF-002, INF-007, TIM-001.

One parser per file format, because the four formats agree on nothing: HadCRUT5
writes ``1850-01``, Berkeley Earth writes percent-commented fixed columns with
``NaN`` in the ones it has no value for, NOAA writes decimal years, and LISIRD
writes days since 1610. Each parser converts to the framework's master time
coordinate and declares what it did.

**The epoch is declared, not assumed.** Every loader converts calendar time to
seconds since ``EPOCH``. Two records loaded through this module are therefore on
the same clock and can be compared without anyone doing arithmetic by hand,
which is the sort of arithmetic that silently puts a series a year out.

**A temperature anomaly is in kelvin, not degrees Celsius.** A *difference* of
one degree Celsius is a difference of one kelvin exactly; only an absolute
temperature carries the 273.15 offset. Declaring an anomaly as ``degC`` would
make the unit registry add that offset at the first conversion and be wrong by
273.15 K. The anomaly series therefore carry ``K``, and their notes say the
value is a difference from a baseline rather than a temperature.

**A gap is a shorter series.** Berkeley Earth writes ``NaN`` in its five-year
columns until it has five years; those rows are dropped rather than carried,
because ``TimeSeries`` refuses a sentinel a consumer might do arithmetic on.
"""

from __future__ import annotations

import csv
from collections.abc import Callable
from datetime import UTC, datetime, timedelta
from pathlib import Path

import numpy as np
import yaml

from mfespf.core.errors import ConfigurationError, ProvenanceError
from mfespf.core.hashing import file_digest
from mfespf.core.versioning import Version
from mfespf.data.adapter import HistoricalDataAdapter
from mfespf.data.records import (
    ProvenanceClass,
    SamplingConvention,
    TimeSeries,
    UncertaintyKind,
)
from mfespf.registries.datasets import DatasetRecord, TransformationStep

__all__ = [
    "EPOCH",
    "EPOCH_LABEL",
    "LOADERS",
    "available",
    "load_berkeley",
    "load_declared",
    "load_giss_saod",
    "load_hadcrut5",
    "load_law_dome_co2",
    "load_noaa_mlo_co2",
    "load_noaa_ohc",
    "load_nrltsi2",
    "seconds_since_epoch",
]

ROOT = Path(__file__).resolve().parents[3]
SOURCES = ROOT / "data" / "sources.yaml"
RAW = ROOT / "data" / "raw"

EPOCH = datetime(1850, 1, 1, tzinfo=UTC)
EPOCH_LABEL = "1850-01-01T00:00:00Z"
"""The framework's calendar reference for historical records.

1850 because it is the start of the instrumental records here. Times before it
are negative, which is fine: a run replaying the solar reconstruction from 1610
starts its master clock at a negative value and says so, which is more honest
than shifting the epoch per dataset and losing the ability to compare them.
"""


def seconds_since_epoch(moment: datetime) -> float:
    return (moment - EPOCH).total_seconds()


def _mid_month(year: int, month: int) -> float:
    """The midpoint of a calendar month, as seconds since the epoch.

    Midpoint rather than start: the sample is an interval mean over the month,
    and ``resample`` builds its cells centred on the sample times. Using the
    month start would put every cell half a month early.
    """
    start = datetime(year, month, 1, tzinfo=UTC)
    end = (
        datetime(year + 1, 1, 1, tzinfo=UTC)
        if month == 12
        else datetime(year, month + 1, 1, tzinfo=UTC)
    )
    return seconds_since_epoch(start + (end - start) / 2)


def _mid_year(year: int) -> float:
    start = datetime(year, 1, 1, tzinfo=UTC)
    end = datetime(year + 1, 1, 1, tzinfo=UTC)
    return seconds_since_epoch(start + (end - start) / 2)


def _decimal_year(value: float) -> float:
    """``1955.5`` -> seconds since the epoch, honouring leap years.

    Interpolating within the actual year rather than assuming 365.25 days keeps
    a mid-year sample at the actual middle of that year, which matters once two
    records with different conventions are put on one clock.
    """
    year = int(np.floor(value))
    fraction = float(value) - year
    start = datetime(year, 1, 1, tzinfo=UTC)
    end = datetime(year + 1, 1, 1, tzinfo=UTC)
    return seconds_since_epoch(start + fraction * (end - start))


def _days_since(days: float, base: datetime) -> float:
    return seconds_since_epoch(base + timedelta(days=float(days)))


# --- one parser per format --------------------------------------------------


def load_hadcrut5(path: str | Path, **identity: object) -> TimeSeries:
    """HadCRUT5 global-mean monthly anomaly with its 95 percent interval."""
    times: list[float] = []
    values: list[float] = []
    half_widths: list[float] = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            stamp = row["Time"].strip()
            year, month = (int(part) for part in stamp.split("-"))
            times.append(_mid_month(year, month))
            values.append(float(row["Anomaly (deg C)"]))
            lower = float(row["Lower confidence limit (2.5%)"])
            upper = float(row["Upper confidence limit (97.5%)"])
            half_widths.append(0.5 * (upper - lower))

    return TimeSeries(
        variable="surface_temperature_anomaly",
        times=np.array(times),
        values=np.array(values),
        units="K",
        sampling=SamplingConvention.INTERVAL_MEAN,
        uncertainty=np.array(half_widths),
        uncertainty_kind=UncertaintyKind.HALF_WIDTH_95,
        note=(
            "anomaly relative to 1961-1990, not an absolute temperature; declared in K "
            "because a difference in degrees Celsius is a difference in kelvin and the "
            "273.15 offset does not apply to it"
        ),
        **identity,  # type: ignore[arg-type]
    )


def load_berkeley(path: str | Path, *, sea_ice: str = "air", **identity: object) -> TimeSeries:
    """Berkeley Earth annual global anomaly.

    ``sea_ice`` selects between the two series the file carries, which are
    different physical definitions rather than two estimates of one thing: the
    air-temperature variant extrapolates over sea ice from surface air, the
    water variant from sea-surface water, and they diverge in the recent record
    because water temperature is bounded below by the freezing point. Recorded
    in the note, since a comparison against the wrong one would look like a
    model error.
    """
    columns = {"air": (1, 2), "water": (5, 6)}
    if sea_ice not in columns:
        raise ConfigurationError(f"sea_ice must be one of {sorted(columns)}, got {sea_ice!r}")
    value_column, uncertainty_column = columns[sea_ice]

    times: list[float] = []
    values: list[float] = []
    uncertainties: list[float] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("%"):
            continue
        fields = stripped.split()
        anomaly, uncertainty = fields[value_column], fields[uncertainty_column]
        if anomaly.lower() == "nan" or uncertainty.lower() == "nan":
            continue  # a gap is a shorter series, never a sentinel
        times.append(_mid_year(int(fields[0])))
        values.append(float(anomaly))
        uncertainties.append(float(uncertainty))

    return TimeSeries(
        variable="surface_temperature_anomaly",
        times=np.array(times),
        values=np.array(values),
        units="K",
        sampling=SamplingConvention.INTERVAL_MEAN,
        uncertainty=np.array(uncertainties),
        uncertainty_kind=UncertaintyKind.HALF_WIDTH_95,
        note=(
            f"anomaly relative to 1951-1980, not an absolute temperature; {sea_ice} "
            f"temperature used over sea ice, which is a definitional choice and not an "
            f"estimate; the published uncertainty is a 95 percent interval"
        ),
        **identity,  # type: ignore[arg-type]
    )


def load_noaa_ohc(path: str | Path, *, basin: str = "WO", **identity: object) -> TimeSeries:
    """NOAA/NCEI ocean heat content, converted from 10^22 J to joules.

    The file's values are in units of 10^22 J, which is not zettajoules and is
    a factor of ten away from it. The conversion happens here, at a declared
    boundary, rather than being left for a reader to notice (DAT-004).
    """
    scale = 1.0e22
    times: list[float] = []
    values: list[float] = []
    errors: list[float] = []
    header: list[str] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0].upper() == "YEAR":
            header = [f.upper() for f in fields]
            continue
        if not header:
            raise ProvenanceError(f"{path} has no column header; its columns cannot be named")
        if basin.upper() not in header:
            raise ConfigurationError(f"{path} has no {basin!r} column; it declares {header[1:]}")
        index = header.index(basin.upper())
        times.append(_decimal_year(float(fields[0])))
        values.append(float(fields[index]) * scale)
        errors.append(float(fields[index + 1]) * scale)

    return TimeSeries(
        variable="ocean_heat_content_anomaly",
        times=np.array(times),
        values=np.array(values),
        units="J",
        sampling=SamplingConvention.INTERVAL_MEAN,
        uncertainty=np.array(errors),
        uncertainty_kind=UncertaintyKind.STANDARD_ERROR,
        note=(
            f"{basin} basin; converted from the file's 10^22 J to J at load. A whole-ocean "
            "total, not an areal density: converting to the framework's J m^-2 basis needs a "
            "declared ocean area and is a further transformation"
        ),
        **identity,  # type: ignore[arg-type]
    )


def load_nrltsi2(path: str | Path, **identity: object) -> TimeSeries:
    """NRLTSI2 total solar irradiance, times given as days since 1610-01-01."""
    base = datetime(1610, 1, 1, tzinfo=UTC)
    times: list[float] = []
    values: list[float] = []
    uncertainties: list[float] = []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if "days since 1610-01-01" not in header[0]:
            raise ProvenanceError(
                f"{path} does not declare the epoch this loader assumes; its first column is "
                f"{header[0]!r}. Reading it anyway would put the series out by whatever the "
                f"difference is"
            )
        for row in reader:
            if not row:
                continue
            times.append(_days_since(float(row[0]), base))
            values.append(float(row[1]))
            uncertainties.append(float(row[2]))

    return TimeSeries(
        variable="total_solar_irradiance",
        times=np.array(times),
        values=np.array(values),
        units="W m^-2",
        sampling=SamplingConvention.INTERVAL_MEAN,
        uncertainty=np.array(uncertainties),
        uncertainty_kind=UncertaintyKind.STANDARD_UNCERTAINTY,
        note=(
            "NRLTSI2 model reconstruction at 1 au; times converted from days since "
            f"1610-01-01 to seconds since {EPOCH_LABEL}, so values before 1850 are negative"
        ),
        **identity,  # type: ignore[arg-type]
    )


def load_law_dome_co2(path: str | Path, **identity: object) -> TimeSeries:
    """The Law Dome ice-core CO2 spline: column 6 of a ten-column file.

    The file carries three gases side by side, each with its own year column,
    because the splines were fitted separately. Reading column 1 as the year for
    column 6 happens to work here and would break silently if the CH4 record
    were ever trimmed, so the CO2 year column is used.
    """
    times: list[float] = []
    values: list[float] = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 7:
            continue
        try:
            year, co2 = float(fields[4]), float(fields[5])
        except ValueError:
            continue
        if not 0.0 < year < 2100.0 or not 150.0 < co2 < 1000.0:
            continue
        times.append(_mid_year(int(year)))
        values.append(co2)

    return TimeSeries(
        variable="co2_mole_fraction",
        times=np.array(times),
        values=np.array(values),
        units="ppm",
        sampling=SamplingConvention.INTERVAL_MEAN,
        note=(
            "a smoothed spline through ice-core air, not a direct atmospheric measurement; "
            "firn smoothing has already removed variability faster than roughly a decade, "
            "and that filter is not recorded in this file or recoverable from it"
        ),
        **identity,  # type: ignore[arg-type]
    )


def load_noaa_mlo_co2(path: str | Path, **identity: object) -> TimeSeries:
    """Mauna Loa annual means: year, mean, uncertainty, hash-commented."""
    times: list[float] = []
    values: list[float] = []
    uncertainties: list[float] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        times.append(_mid_year(int(float(fields[0]))))
        values.append(float(fields[1]))
        uncertainties.append(float(fields[2]))

    return TimeSeries(
        variable="co2_mole_fraction",
        times=np.array(times),
        values=np.array(values),
        units="ppm",
        sampling=SamplingConvention.INTERVAL_MEAN,
        uncertainty=np.array(uncertainties),
        uncertainty_kind=UncertaintyKind.STANDARD_UNCERTAINTY,
        note=(
            "one station used as a global mean; CO2 is well mixed on the timescales used "
            "here, but the seasonal cycle and the interhemispheric gradient are absent"
        ),
        **identity,  # type: ignore[arg-type]
    )


def load_giss_saod(path: str | Path, *, column: str = "global", **identity: object) -> TimeSeries:
    """GISS stratospheric aerosol optical depth: decimal year and three columns.

    ``column`` picks the global mean or a hemisphere. The default is the global
    mean, which is what a single-column model can use; the hemispheric pair is
    there because a volcanic cloud is not symmetric, and a configuration that
    could represent that difference should be able to read it.
    """
    columns = {"global": 1, "northern": 2, "southern": 3}
    if column not in columns:
        raise ConfigurationError(f"column must be one of {sorted(columns)}, got {column!r}")
    index = columns[column]

    times: list[float] = []
    values: list[float] = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 4:
            continue
        try:
            when, tau = float(fields[0]), float(fields[index])
        except ValueError:
            continue
        if not 1000.0 < when < 2200.0:
            continue
        times.append(_decimal_year(when))
        values.append(tau)

    return TimeSeries(
        variable="stratospheric_aerosol_optical_depth",
        times=np.array(times),
        values=np.array(values),
        units="1",
        sampling=SamplingConvention.INTERVAL_MEAN,
        note=(
            f"{column} mean at 550 nm. Before the satellite era this is assembled from "
            "ground-based and indirect evidence rather than measured, so the early part "
            "inherits those assumptions"
        ),
        **identity,  # type: ignore[arg-type]
    )


LOADERS: dict[str, Callable[..., TimeSeries]] = {
    "hadcrut5": load_hadcrut5,
    "berkeley": load_berkeley,
    "noaa_ohc": load_noaa_ohc,
    "nrltsi2": load_nrltsi2,
    "law_dome_co2": load_law_dome_co2,
    "giss_saod": load_giss_saod,
    "noaa_mlo_co2": load_noaa_mlo_co2,
}


# --- the declared sources ---------------------------------------------------


def _sources() -> list[dict]:
    if not SOURCES.exists():
        raise ProvenanceError(f"{SOURCES} is missing; the declared sources are the record")
    return yaml.safe_load(SOURCES.read_text(encoding="utf-8"))["sources"]


def available() -> list[str]:
    """Dataset IDs whose bytes are present and match their recorded checksum."""
    present = []
    for source in _sources():
        path = RAW / source["filename"]
        if path.exists() and file_digest(path) == source.get("sha256"):
            present.append(source["dataset_id"])
    return present


def load_declared(dataset_id: str, **loader_arguments: object) -> HistoricalDataAdapter:
    """Load a declared source into an adapter, verifying it against the record.

    The checksum is checked before the file is parsed, not after. A record whose
    bytes have changed under the entry that describes them is a different
    dataset, and reading it first and noticing second would mean the values had
    already reached the caller.
    """
    entry = next((s for s in _sources() if s["dataset_id"] == dataset_id), None)
    if entry is None:
        raise ProvenanceError(
            f"{dataset_id!r} is not a declared source; {SOURCES.name} declares "
            f"{[s['dataset_id'] for s in _sources()]}"
        )

    path = RAW / entry["filename"]
    if not path.exists():
        raise ProvenanceError(
            f"{path} is not present. The bytes are not in the repository; fetch them with "
            f"`python tools/fetch_data.py --only {dataset_id}`"
        )
    recorded = entry.get("sha256") or ""
    actual = file_digest(path)
    if recorded and actual != recorded:
        raise ProvenanceError(
            f"{path.name} does not match the checksum recorded for {dataset_id!r} "
            f"({recorded[:16]} vs {actual[:16]}); the file changed under its record, which is "
            f"a fact to be resolved rather than a mismatch to be ignored (DAT-009)"
        )

    loader_name = entry.get("loader")
    if loader_name not in LOADERS:
        raise ConfigurationError(
            f"source {dataset_id!r} declares loader {loader_name!r}; known loaders are "
            f"{sorted(LOADERS)}"
        )
    provenance = ProvenanceClass(entry["provenance_class"])
    series = LOADERS[loader_name](
        path,
        dataset_id=dataset_id,
        dataset_version=str(entry.get("dataset_version", "1.0.0")),
        provenance_class=provenance,
        **loader_arguments,
    )

    record = DatasetRecord.from_file(
        record_id=dataset_id,
        version=Version.parse(str(entry.get("dataset_version", "1.0.0"))),
        path=path,
        title=entry["title"],
        source=f"{entry['provider']} -- {entry['url']}",
        variables=(series.variable,),
        units=(series.units,),
        temporal_coverage=(f"{series.coverage[0]:g} s", f"{series.coverage[1]:g} s"),
        spatial_coverage=entry.get("spatial_coverage", ""),
        uncertainty_record=(
            f"{series.uncertainty_kind.value}; "
            + (
                "usable directly in a quadrature combination"
                if series.uncertainty_kind.is_one_sigma
                else "NOT a standard uncertainty -- roughly two sigma, so combining it as "
                "though it were one would halve whatever it entered (UNC-006)"
            )
        ),
        quality_flags=tuple(entry.get("quality_flags", ())),
        acquisition_record=(
            f"fetched by tools/fetch_data.py on {entry.get('retrieved', 'an unrecorded date')}; "
            f"licence: {entry.get('licence', 'unrecorded')}; citation: "
            f"{entry.get('citation', 'unrecorded')}"
        ),
        transformations=(
            TransformationStep(
                operation=f"parse:{loader_name}",
                parameters=(
                    f"epoch={EPOCH_LABEL}; sampling={series.sampling.value}; "
                    f"units={series.units}"
                    + "".join(f"; {k}={v}" for k, v in sorted(loader_arguments.items()))
                ),
                software=f"mfespf.data.loaders.{LOADERS[loader_name].__name__}",
                uncertainty_effect=(
                    f"published uncertainty carried through unchanged as "
                    f"{series.uncertainty_kind.value}; not propagated into anything"
                ),
            ),
        ),
        provenance_class=provenance,
    )
    return HistoricalDataAdapter(record, {series.variable: series})
