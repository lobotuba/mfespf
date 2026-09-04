"""The historical data adapter: a declared record on the model's clock.

Charter 18 candidate module 12 and section 19 stage 6; requirements DAT-007,
DAT-008, INF-001, INF-002, INF-007, TIM-001, TIM-004, KER-010.

The adapter is the only way a recorded series reaches a module, and it exists so
that three things happen in one place rather than three: the record is registered
with its checksum and provenance, the resampling onto the master clock is
performed and reported, and the transformation is written back onto the dataset
record so the preprocessing can be replayed from the registry alone.

It refuses two things that a looser adapter would let through. A variable the
dataset record does not declare cannot be read, because a series that is not in
the record's variable list is not covered by the record's checksum. And a run
that reaches outside the record's coverage fails rather than extrapolating: a
boundary condition invented past the end of a record is not a historical input,
and DAT-007 forbids substituting a value for a missing one.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from mfespf.contracts.module import FidelityLevel, ModuleStatus
from mfespf.contracts.spectral import SpectralGrid
from mfespf.core.errors import ConfigurationError, ProvenanceError
from mfespf.data.records import (
    ProvenanceClass,
    ResamplingReport,
    SamplingConvention,
    TimeSeries,
    resample,
)
from mfespf.physics.radiation import normalized_band_weights
from mfespf.physics.solar import SolarSource
from mfespf.registries.datasets import DatasetRecord, DatasetRegistry, DatasetRole, DatasetUse

__all__ = ["HistoricalDataAdapter", "RecordedSolarSource"]


class HistoricalDataAdapter:
    """One declared dataset, its series, and everything done to them."""

    def __init__(self, record: DatasetRecord, series: Mapping[str, TimeSeries]) -> None:
        if not series:
            raise ConfigurationError(f"dataset {record.record_id!r} supplies no series")
        undeclared = sorted(set(series) - set(record.variables))
        if undeclared:
            raise ProvenanceError(
                f"dataset {record.record_id!r} supplies series {undeclared} that its record does "
                f"not declare; a variable outside the declared list is outside the checksum that "
                f"is supposed to identify it (DAT-009)"
            )
        units = dict(zip(record.variables, record.units, strict=False))
        for name, entry in series.items():
            declared = units.get(name)
            if declared is not None and declared != entry.units:
                raise ProvenanceError(
                    f"dataset {record.record_id!r} declares {name!r} in {declared!r} but the "
                    f"series carries {entry.units!r}; convert at a declared boundary (DAT-004)"
                )
            if entry.provenance_class is not record.provenance_class:
                raise ProvenanceError(
                    f"series {name!r} declares provenance {entry.provenance_class.value!r} but "
                    f"its dataset record declares {record.provenance_class.value!r}"
                )

        self.record = record
        self._series = dict(series)
        self._reports: list[ResamplingReport] = []
        self._prepared: dict[str, TimeSeries] = {}

    # --- identity ---------------------------------------------------------
    @property
    def dataset_id(self) -> str:
        return self.record.record_id

    @property
    def version(self) -> str:
        return str(self.record.version)

    @property
    def provenance_class(self) -> ProvenanceClass:
        return self.record.provenance_class

    def series(self, variable: str) -> TimeSeries:
        try:
            return self._series[variable]
        except KeyError:
            raise ProvenanceError(
                f"dataset {self.dataset_id!r} has no series {variable!r}; it declares "
                f"{sorted(self._series)}"
            ) from None

    # --- preparation ------------------------------------------------------
    def prepare(
        self,
        variable: str,
        target_times: np.ndarray,
        *,
        method: str | None = None,
        target_interval: float | None = None,
    ) -> TimeSeries:
        """Resample a series onto the model clock, keeping the report.

        Called once per variable before the run rather than per step. Doing it
        once means the lost frequency content is reported at a single declared
        point instead of being smeared across thousands of silent lookups.
        """
        prepared, report = resample(
            self.series(variable),
            target_times,
            method=method,
            target_interval=target_interval,
        )
        self._reports.append(report)
        self._prepared[variable] = prepared
        return prepared

    def prepared(self, variable: str) -> TimeSeries:
        try:
            return self._prepared[variable]
        except KeyError:
            raise ProvenanceError(
                f"series {variable!r} has not been prepared onto the model clock; call "
                f"prepare() before the run so the resampling is reported once (TIM-004)"
            ) from None

    @property
    def reports(self) -> tuple[ResamplingReport, ...]:
        return tuple(self._reports)

    def transformed_record(self) -> DatasetRecord:
        """The record with every transformation attached (DAT-008, INF-007).

        Returned rather than mutated: the original record is what the checksum
        identifies, and the transformed one is a separate object that says what
        was done to it. Editing the first in place would leave no way to tell
        the raw record from the prepared one.
        """
        from dataclasses import replace

        return replace(
            self.record,
            transformations=(
                *self.record.transformations,
                *(r.as_transformation_step() for r in self._reports),
            ),
        )

    def register(
        self,
        registry: DatasetRegistry,
        *,
        role: DatasetRole,
        claim_id: str,
        run_id: str | None = None,
        influenced: tuple[str, ...] = (),
        note: str = "",
    ) -> DatasetUse:
        """Register the transformed record and declare what it was used for."""
        record = self.transformed_record()
        if record.record_id not in registry:
            registry.register(record)
        return registry.declare_use(
            DatasetUse(
                dataset_id=record.record_id,
                dataset_version=str(record.version),
                role=role,
                claim_id=claim_id,
                variables=tuple(sorted(self._prepared or self._series)),
                run_id=run_id,
                influenced=influenced,
                note=note,
            )
        )

    def describe(self) -> str:
        lines = [
            f"dataset {self.dataset_id}@{self.version} "
            f"({self.provenance_class.value}) from {self.record.source}",
            f"  checksum: {self.record.checksum[:16]}",
        ]
        lines.extend(f"  {s.describe()}" for s in self._series.values())
        lines.extend(f"  transformation: {r.describe()}" for r in self._reports)
        return "\n".join(lines)


class RecordedSolarSource(SolarSource):
    """Solar irradiance played back from a declared record (KER-010).

    A third implementation of the solar interface beside the constant and forced
    sources, and the one that makes a prescribed-history run possible: the
    irradiance is read from a dataset rather than computed, so the top-of-
    atmosphere boundary over the replay period is whatever the record says.

    The spectral shape is still a black body, because the record is broadband.
    That is a declared limitation rather than an omission: a run using this
    cannot make any claim about the wavelength dependence of solar forcing, and
    the contract says so.
    """

    def __init__(
        self,
        *,
        grid: SpectralGrid,
        adapter: HistoricalDataAdapter,
        variable: str = "total_solar_irradiance",
        effective_temperature: float = 5772.0,
        module_id: str = "solar",
        status: ModuleStatus = ModuleStatus.PRESCRIBED,
    ) -> None:
        from dataclasses import replace

        series = adapter.prepared(variable)
        if series.units != "W m^-2":
            raise ConfigurationError(
                f"irradiance record is in {series.units!r}; the solar interface requires W m^-2"
            )
        super().__init__(
            module_id=module_id,
            grid=grid,
            fidelity=FidelityLevel.PRESCRIBED,
            implementation_name=(
                f"recorded total solar irradiance from {adapter.dataset_id}@{adapter.version} "
                f"({grid.resolution_identifier})"
            ),
            status=status,
            assumptions=(
                f"total solar irradiance follows the {adapter.dataset_id} record exactly",
                f"the record is {series.provenance_class.value}, "
                f"{series.sampling.value}, resolving periods down to "
                f"{series.resolved_period:g} s",
                f"spectral shape of a black body at {effective_temperature:g} K, since the "
                "record is broadband",
                "values are held constant between record samples",
            ),
            cost_work_units=1.0,
        )
        self.adapter = adapter
        self.variable = variable
        self.series = series
        self.effective_temperature = effective_temperature
        self._weights = normalized_band_weights(grid, effective_temperature)
        self.contract = replace(
            self.contract,
            physical_basis=(
                "Total solar irradiance read from a declared dataset over the replay "
                "period, distributed over wavelength by a declared spectral shape. The "
                "irradiance is an input to this run, not a result of it."
            ),
            validity_domain=(
                f"The record's coverage, [{series.coverage[0]:g}, {series.coverage[1]:g}] s. "
                f"Outside it the module fails rather than extrapolating. Because the "
                f"irradiance is prescribed, no claim about solar variability is in scope "
                f"for a run using this (TIM-006)."
            ),
            correctness_tests=(
                "test_a_recorded_source_reproduces_its_record_exactly",
                "test_running_past_the_record_coverage_fails",
                "test_spectral_integral_matches_total_irradiance",
            ),
            known_limitations=(
                "The record is broadband, so the spectral shape is assumed rather than "
                "recorded and no claim about spectral solar forcing is admissible.",
                f"Resampling onto the model clock removed variability at periods below "
                f"{series.resolved_period:g} s; see the adapter's transformation record.",
            ),
        )

    def total_irradiance(self, time: float) -> float:
        # The whole step must be covered, not merely its start. A step is an
        # integration over [t, t + dt], and running one whose second half the
        # record does not reach would extend the record by half a step under the
        # zero-order hold -- a small fabrication, but a fabrication (DAT-007).
        dt = self._context.dt if self._context else 0.0
        self.series.require_coverage(time, time + dt)
        return self.series.value_at(time)

    def spectral_shape(self) -> np.ndarray:
        return self._weights

    def parameters(self) -> dict[str, Any]:
        return {
            **super().parameters(),
            "dataset": f"{self.adapter.dataset_id}@{self.adapter.version}",
            "dataset_checksum": self.adapter.record.checksum,
            "dataset_variable": self.variable,
            "dataset_provenance_class": self.adapter.provenance_class.value,
            "record_sampling": SamplingConvention(self.series.sampling).value,
            "record_coverage_s": list(self.series.coverage),
            "spectral_shape_effective_temperature": self.effective_temperature,
        }
