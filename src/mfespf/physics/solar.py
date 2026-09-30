"""Solar source module.

Charter 9.1 and requirements SOL-001, SOL-002, SPC-001, SPC-002, KER-009,
KER-010.

Two interchangeable implementations satisfy one interface: a broadband source
carrying the whole irradiance in a single bin, and a banded source distributing
it over an arbitrary wavelength grid.  They are the KER-009 demonstration that
spectral resolution is configurable, and the KER-010 demonstration that a module
implementation is replaceable -- both answer the same interface contract, and a
run manifest records which one was active.

Conservation role (charter 9.1): the solar module supplies irradiance.  It does
not transfer absorbed energy into Earth modules; the orbital module owns the
top-of-atmosphere boundary and accounts for what crosses it.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np

from mfespf.contracts.cost import CostDeclaration, ResourceClass
from mfespf.contracts.module import FidelityLevel, Module, ModuleContract, ModuleStatus
from mfespf.contracts.spectral import (
    IntegrationConvention,
    SpectralCoordinate,
    SpectralGrid,
)
from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError
from mfespf.core.quantities import (
    Field,
    Provenance,
    QuantityKind,
    SignConvention,
    SpatialSupport,
    TimeSemantics,
    Uncertainty,
)
from mfespf.core.versioning import InterfaceContract, Version
from mfespf.physics.equations import EQUATION_IDS
from mfespf.physics.radiation import normalized_band_weights

__all__ = [
    "SOLAR_SOURCE_INTERFACE",
    "ConstantSpectrumSolarSource",
    "ForcedSolarSource",
    "MeasuredSpectrumSolarSource",
    "SolarSource",
    "banded_grid",
    "broadband_grid",
]

V1 = Version(1, 0, 0)

SOLAR_SOURCE_INTERFACE = InterfaceContract(
    interface_id="solar_source",
    version=V1,
    physical_semantics={
        "total_irradiance_units": "W m^-2",
        "spectral_irradiance_units": "W m^-2 m^-1",
        "reference_distance": "1 au",
        "spectral_density_basis": "per_metre_of_wavelength",
        "sign_convention": "positive_definite",
        "time_semantics": "interval_average",
    },
    optional_capabilities=frozenset({"spectral_irradiance"}),
)


def broadband_grid(lower_m: float = 1.0e-7, upper_m: float = 1.0e-4) -> SpectralGrid:
    """A single bin spanning the solar spectrum: the lowest spectral fidelity."""
    return SpectralGrid(
        grid_id="broadband_100nm_100um",
        coordinate=SpectralCoordinate.WAVELENGTH,
        coordinate_units="m",
        edges=np.array([lower_m, upper_m]),
        integration_convention=IntegrationConvention.SPECTRAL_DENSITY,
        density_units="W m^-2 m^-1",
        resolution_identifier="broadband-1band",
    )


def banded_grid(
    n_bands: int = 16, lower_m: float = 1.0e-7, upper_m: float = 1.0e-4
) -> SpectralGrid:
    """A logarithmically spaced band grid over the same support.

    Log spacing rather than linear: solar and terrestrial emission span three
    decades of wavelength, and linear bins would put almost every bin in the
    far infrared where there is little solar energy to resolve.
    """
    if n_bands < 1:
        raise ConfigurationError("n_bands must be at least 1")
    edges = np.logspace(np.log10(lower_m), np.log10(upper_m), n_bands + 1)
    return SpectralGrid(
        grid_id=f"log_bands_{n_bands}",
        coordinate=SpectralCoordinate.WAVELENGTH,
        coordinate_units="m",
        edges=edges,
        integration_convention=IntegrationConvention.SPECTRAL_DENSITY,
        density_units="W m^-2 m^-1",
        resolution_identifier=f"log-{n_bands}band",
    )


class SolarSource(Module):
    """Base for solar-source implementations.

    Subclasses supply the spectrum; the contract, outputs and uncertainty
    reporting are shared so that two implementations cannot drift apart in what
    they claim to provide.
    """

    def __init__(
        self,
        *,
        module_id: str = "solar",
        grid: SpectralGrid,
        fidelity: FidelityLevel,
        implementation_name: str,
        status: ModuleStatus = ModuleStatus.PRESCRIBED,
        assumptions: tuple[str, ...] = (),
        cost_work_units: float = 1.0,
    ) -> None:
        contract = ModuleContract(
            module_id=module_id,
            module_name=implementation_name,
            module_version=V1,
            interface_version=V1,
            fidelity_level=fidelity,
            physical_basis=(
                "Measured total solar irradiance at 1 au, distributed over wavelength by a "
                "declared spectral shape."
            ),
            governing_equation_ids=(EQUATION_IDS["planck"],),
            equation_variables=(
                (
                    EQUATION_IDS["planck"],
                    (
                        # Planck is evaluated here only to shape the spectrum,
                        # never to set the total, so its temperature is the
                        # declared effective temperature rather than any state.
                        # A subclass that prescribes a measured spectrum instead
                        # does not supply this parameter and the check says so.
                        ("T", "spectral_shape_effective_temperature", "K"),
                        # The wavelengths the relation is actually evaluated at
                        # are the band edges of this module's grid.
                        ("lambda", "spectral_edges", "m"),
                    ),
                ),
            ),
            assumptions=assumptions,
            validity_domain=(
                "Irradiance at the reference distance of 1 au, over the declared spectral "
                "support. No claim is made about solar variability beyond the supplied data."
            ),
            produced_outputs=("solar_total_irradiance_1au", "solar_spectral_irradiance_1au"),
            conserved_quantities=("energy",),
            supported_spectral_representations=(grid.resolution_identifier,),
            provided_interfaces=(SOLAR_SOURCE_INTERFACE,),
            uncertainty_outputs=("total_irradiance_measurement", "spectral_shape_structural"),
            correctness_tests=(
                "test_spectral_integral_matches_total_irradiance",
                "test_solar_implementations_agree_broadband",
            ),
            known_limitations=(
                "No solar-cycle, flare, particle or magnetic output at this fidelity.",
            ),
            fidelity_justification=(
                "A broadband constant carries no spectral shape, so the surface's "
                "wavelength-dependent albedo has nothing to act on and shortwave "
                "absorption cannot be separated from thermal emission band by band."
            ),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(n_bands)",
                memory_complexity="O(n_bands)",
                scaling_variables=("n_bands",),
                work_units_per_step=cost_work_units,
                benchmark_configuration=f"{grid.n_bins} bands, global scalar",
                peak_memory_bytes=8 * grid.n_bins + 4096,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
                lower_fidelity_alternatives=("solar_constant_broadband",),
            ),
        )
        super().__init__(contract, status=status)
        self.grid = grid
        self._total_irradiance = 0.0
        self._spectral_irradiance = np.zeros(grid.n_bins)
        self._spectral_capture_fraction = 1.0

    # --- to be provided by implementations --------------------------------
    def total_irradiance(self, time: float) -> float:
        raise NotImplementedError

    def spectral_shape(self) -> np.ndarray:
        """Normalised weights per bin; they need not sum to one (see below)."""
        raise NotImplementedError

    # --- lifecycle --------------------------------------------------------
    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        # Not zero: a run replaying a dated record starts its master clock where
        # the record does, and a source that published its value for t = 0 would
        # be publishing one from outside the record entirely (TIM-001).
        self._update(float(boundary_state.get("master_start_time", 0.0)))
        self.mark_initialized()

    def compute_tendencies(self) -> list[Any]:
        self._update(self.context.time)
        return []

    def _update(self, time: float) -> None:
        total = self.total_irradiance(time)
        weights = self.spectral_shape()
        captured = float(np.sum(weights))
        if captured <= 0.0:
            raise ConfigurationError(
                f"solar spectral shape for {self.module_id!r} carries no energy"
            )
        # Scale so the spectral integral reproduces the measured total exactly.
        # The fraction of the true spectrum lying outside the grid's support is
        # reported as structural uncertainty rather than being discarded.
        self._total_irradiance = total
        self._spectral_irradiance = (weights / captured) * total / self.grid.widths
        self._spectral_capture_fraction = captured

    def emit_outputs(self) -> dict[str, Field]:
        provenance = Provenance(
            producer=self.module_id,
            # The relation that shaped the spectrum: Planck for a black body,
            # the measured-spectrum record for a measured one.
            derivation=self.contract.governing_equation_ids[0],
            notes=self.contract.module_name,
        )
        stub = self.status is ModuleStatus.STUB
        return {
            "solar_total_irradiance_1au": Field(
                field_id="solar_total_irradiance_1au",
                physical_name="total solar irradiance at 1 au",
                symbol="S_1au",
                values=np.array(self._total_irradiance),
                units="W m^-2",
                quantity_kind=QuantityKind.BOUNDARY_VALUE,
                time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                time_interval=self.context.dt if self._context else None,
                timestamp=self.context.time if self._context else None,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                provenance=provenance,
                uncertainty=Uncertainty(
                    standard_uncertainty=0.5,
                    method="reported measurement uncertainty of the reference composite",
                    unresolved=("solar_cycle_phase", "instrument_absolute_scale"),
                ),
                stub=stub,
            ),
            "solar_spectral_irradiance_1au": Field(
                field_id="solar_spectral_irradiance_1au",
                physical_name="spectral solar irradiance at 1 au",
                symbol="S_lambda",
                values=self._spectral_irradiance,
                units="W m^-2 m^-1",
                quantity_kind=QuantityKind.BOUNDARY_VALUE,
                time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                time_interval=self.context.dt if self._context else None,
                timestamp=self.context.time if self._context else None,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                spectral_coordinate=self.grid.grid_id,
                dimensions=("band",),
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                provenance=provenance,
                metadata={"spectral_grid": self.grid},
                uncertainty=Uncertainty(
                    method="none",
                    unresolved=(
                        "spectral_shape_structural",
                        f"support_truncation_fraction={1.0 - self._spectral_capture_fraction:.4g}",
                    ),
                ),
                stub=stub,
            ),
        }

    def report_uncertainty(self) -> Any:
        from mfespf.contracts.module import UncertaintyReport

        return UncertaintyReport(
            module_id=self.module_id,
            quantified={"total_irradiance_standard_uncertainty_W_m2": 0.5},
            method="declared measurement uncertainty; no propagation performed",
            distributional_assumptions="unresolved",
            correlations="unresolved",
            unresolved=(
                "spectral_shape_structural_uncertainty",
                "solar_variability_outside_supplied_record",
            ),
        )

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "spectral_grid_id": self.grid.grid_id,
            "spectral_resolution_identifier": self.grid.resolution_identifier,
            "spectral_edges": self.grid.edges.tolist(),
            "integration_convention": self.grid.integration_convention.value,
        }

    def checkpoint(self) -> dict[str, Any]:
        return {
            "total_irradiance": self._total_irradiance,
            "spectral_irradiance": self._spectral_irradiance.copy(),
        }

    def restore(self, checkpoint: dict[str, Any]) -> None:
        self._total_irradiance = float(checkpoint["total_irradiance"])
        self._spectral_irradiance = np.asarray(checkpoint["spectral_irradiance"], dtype=float)


class ConstantSpectrumSolarSource(SolarSource):
    """Constant total irradiance with a fixed black-body spectral shape.

    The same class serves both fidelities: with a one-bin grid it is the
    broadband implementation, and with an N-band grid it is the banded one.  That
    is deliberate -- if the two were separate classes, agreement between them
    would test the two implementations rather than the spectral machinery.  The
    thing under test in KER-009 is that the *grid* is interchangeable.
    """

    def __init__(
        self,
        *,
        grid: SpectralGrid,
        total_irradiance: float | None = None,
        effective_temperature: float = 5772.0,
        module_id: str = "solar",
        status: ModuleStatus = ModuleStatus.PRESCRIBED,
    ) -> None:
        fidelity = FidelityLevel.REDUCED if grid.n_bins == 1 else FidelityLevel.STANDARD
        super().__init__(
            module_id=module_id,
            grid=grid,
            fidelity=fidelity,
            implementation_name=(
                f"constant TSI with black-body spectral shape ({grid.resolution_identifier})"
            ),
            status=status,
            assumptions=(
                "total solar irradiance constant over the simulated interval",
                f"spectral shape of a black body at {effective_temperature:g} K",
                "no solar-cycle, particle or magnetic output",
            ),
            cost_work_units=1.0,
        )
        self._total = (
            total_irradiance
            if total_irradiance is not None
            else value_of("total_solar_irradiance_reference")
        )
        self.effective_temperature = effective_temperature
        self._weights = normalized_band_weights(grid, effective_temperature)

    def total_irradiance(self, time: float) -> float:
        return self._total

    def spectral_shape(self) -> np.ndarray:
        return self._weights

    def parameters(self) -> dict[str, Any]:
        return {
            **super().parameters(),
            "total_solar_irradiance_1au": self._total,
            "spectral_shape_effective_temperature": self.effective_temperature,
        }


class MeasuredSpectrumSolarSource(SolarSource):
    """Constant total irradiance with the measured solar spectrum's shape.

    The shape is TSIS-1 HSRS version 2 integrated over each band of the grid
    (``physics/solar_spectrum.py``). A black body at 5772 K puts half again
    too much sunlight below 0.4 microns, where ozone absorbs; this one
    reproduces Chou and Suarez's (1999) band fractions to within a few percent
    (``docs/protocols/ozone-solar-heating.md``). The total is still the
    configured irradiance: the spectrum supplies the distribution, not the
    amount.
    """

    def __init__(
        self,
        *,
        grid: SpectralGrid,
        total_irradiance: float | None = None,
        module_id: str = "solar",
        status: ModuleStatus = ModuleStatus.PRESCRIBED,
    ) -> None:
        from mfespf.physics.solar_spectrum import HSRS_TOTAL_W_M2, hsrs_band_irradiance

        if grid.coordinate is not SpectralCoordinate.WAVELENGTH:
            raise ConfigurationError("the measured solar spectrum is tabulated in wavelength")
        super().__init__(
            module_id=module_id,
            grid=grid,
            fidelity=FidelityLevel.REDUCED if grid.n_bins == 1 else FidelityLevel.STANDARD,
            implementation_name=(
                f"constant TSI with the TSIS-1 HSRS v2 measured spectral shape "
                f"({grid.resolution_identifier})"
            ),
            status=status,
            assumptions=(
                "total solar irradiance constant over the simulated interval",
                "spectral shape of the TSIS-1 Hybrid Solar Reference Spectrum, version 2",
                "no solar-cycle, particle or magnetic output",
            ),
            cost_work_units=1.0,
        )
        # Not Planck: nothing here evaluates a black body, and citing it would
        # have the validity gate judge a temperature this source does not have.
        self.contract = replace(
            self.contract,
            governing_equation_ids=(EQUATION_IDS["measured_solar_spectrum"],),
            equation_variables=(),
        )
        self._total = (
            total_irradiance
            if total_irradiance is not None
            else value_of("total_solar_irradiance_reference")
        )
        self._weights = hsrs_band_irradiance(np.asarray(grid.edges, dtype=float)) / HSRS_TOTAL_W_M2

    def total_irradiance(self, time: float) -> float:
        return self._total

    def spectral_shape(self) -> np.ndarray:
        return self._weights

    def parameters(self) -> dict[str, Any]:
        return {
            **super().parameters(),
            "total_solar_irradiance_1au": self._total,
            "spectral_shape": "TSIS-1 HSRS v2 full spectrum, Coddington et al. 2023",
        }


class ForcedSolarSource(ConstantSpectrumSolarSource):
    """A solar source carrying a declared forcing experiment (TIM-007).

    The forcing perturbs total solar irradiance, so its energy enters the system
    through the top-of-atmosphere boundary and is accounted by the transfers that
    were already there. No fiat source term is needed and none is added: a
    forcing experiment must not be the one place in the framework where energy
    appears from nowhere.

    Two consequences follow from where the perturbation is applied, and both are
    published rather than left implicit. A perturbation ``dS`` at 1 au becomes
    ``dS / 4`` of insolation once spread over the rotating sphere, and only the
    absorbed fraction ``1 - alpha`` of that reaches the system as heating. A
    timescale experiment measuring a gain per watt must therefore use the
    realised anomaly, not the nominal one, and both are on the bus.

    The contract is rebuilt rather than inherited: the parent declares that total
    irradiance is constant over the simulated interval, which is exactly what
    this class makes false. A contract that keeps an assumption its module has
    abandoned is worse than no contract (MOD-001).
    """

    def __init__(
        self,
        *,
        grid: SpectralGrid,
        forcing: Any,
        total_irradiance: float | None = None,
        effective_temperature: float = 5772.0,
        module_id: str = "solar",
        status: ModuleStatus = ModuleStatus.PRESCRIBED,
    ) -> None:
        super().__init__(
            grid=grid,
            total_irradiance=total_irradiance,
            effective_temperature=effective_temperature,
            module_id=module_id,
            status=status,
        )
        self.forcing = forcing
        self._anomaly = 0.0
        self.contract = replace(
            self.contract,
            module_name=(
                f"solar source under a declared {forcing.regime} forcing "
                f"({grid.resolution_identifier})"
            ),
            physical_basis=(
                "Measured total solar irradiance at 1 au perturbed by a declared "
                "experimental forcing, distributed over wavelength by a declared "
                "spectral shape. The perturbation is a boundary condition of the "
                "experiment, not a solar physics result."
            ),
            assumptions=(
                "baseline total solar irradiance constant over the simulated interval",
                f"the {forcing.regime} perturbation is applied to the total and leaves "
                "the spectral shape unchanged",
                f"spectral shape of a black body at {effective_temperature:g} K",
                "no solar-cycle, particle or magnetic output",
            ),
            validity_domain=(
                "A controlled forcing experiment. The irradiance history is imposed, so "
                "no statement about actual solar variability is in scope, and TIM-006 "
                "forbids reporting this mode in the language of a historical run."
            ),
            produced_outputs=(
                "solar_total_irradiance_1au",
                "solar_spectral_irradiance_1au",
                "solar_irradiance_anomaly",
            ),
            correctness_tests=(
                "test_unforced_source_reproduces_the_constant_one_exactly",
                "test_forcing_delivers_its_declared_energy_at_any_time_step",
                "test_spectral_integral_matches_total_irradiance",
            ),
            known_limitations=(
                "No solar-cycle, flare, particle or magnetic output at this fidelity.",
                "The perturbation is spectrally flat, so an experiment about the "
                "wavelength dependence of solar forcing cannot be run with it.",
            ),
        )

    def total_irradiance(self, time: float) -> float:
        # The interval mean, not the value at the step start: the forcing classes
        # integrate themselves exactly so that the delivered energy does not
        # depend on the time step. Sampling here would throw that away.
        dt = self._context.dt if self._context else 0.0
        self._anomaly = self.forcing.mean_amplitude(time, time + dt)
        return self._total + self._anomaly

    def emit_outputs(self) -> dict[str, Field]:
        outputs = super().emit_outputs()
        ctx = self.context
        outputs["solar_irradiance_anomaly"] = Field(
            field_id="solar_irradiance_anomaly",
            physical_name="imposed total solar irradiance anomaly at 1 au",
            symbol="dS",
            values=np.array(self._anomaly),
            units="W m^-2",
            quantity_kind=QuantityKind.BOUNDARY_VALUE,
            time_semantics=TimeSemantics.INTERVAL_AVERAGE,
            time_interval=ctx.dt,
            timestamp=ctx.time,
            spatial_support=SpatialSupport.GLOBAL_SCALAR,
            provenance=Provenance(
                producer=self.module_id,
                notes=(
                    f"declared experimental forcing: {self.forcing.describe()}; "
                    f"the insolation anomaly is a quarter of this before albedo"
                ),
            ),
        )
        return outputs

    def parameters(self) -> dict[str, Any]:
        return {
            **super().parameters(),
            "forcing": self.forcing.parameters(),
        }
