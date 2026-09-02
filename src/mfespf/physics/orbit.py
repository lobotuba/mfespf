"""Orbital geometry and top-of-atmosphere mapping.

Charter 9.2 and requirements ORB-001, ORB-002.

This module owns the top-of-atmosphere boundary of the modeled system.  It takes
the solar module's irradiance at the reference distance, applies the Earth-Sun
distance and the illumination geometry, and transfers the resulting energy across
the system boundary to the receiving surface.  Accounting the external input here
rather than inside the solar module keeps charter 9.1's rule intact: the solar
source supplies irradiance and does not itself add absorbed energy to Earth
modules.

ORB-002 -- distributed incident energy reconciling with intercepted solar power --
is a property of ``spherical_mean_factor`` and is tested directly.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from mfespf.contracts.cost import CostDeclaration, ResourceClass
from mfespf.contracts.module import FidelityLevel, Module, ModuleContract, ModuleStatus
from mfespf.contracts.spectral import SpectralGrid
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
from mfespf.core.versioning import InterfaceContract, Version, VersionRange
from mfespf.coupling.conservation import ModuleBudget
from mfespf.coupling.transfers import ConservedQuantity, Transfer
from mfespf.physics.equations import EQUATION_IDS

__all__ = [
    "ORBITAL_GEOMETRY_INTERFACE",
    "OrbitalGeometryModule",
    "orbit_radius",
    "solve_kepler",
    "spherical_mean_factor",
]

V1 = Version(1, 0, 0)
_AU = value_of("astronomical_unit")
_SIDEREAL_YEAR_S = 365.256363004 * 86400.0

ORBITAL_GEOMETRY_INTERFACE = InterfaceContract(
    interface_id="orbital_geometry",
    version=V1,
    physical_semantics={
        "distance_units": "m",
        "toa_irradiance_units": "W m^-2",
        "toa_irradiance_basis": "global_mean_over_planetary_surface",
        "geometry_factor_definition": "projected_disc_area_over_total_surface_area",
        "sign_convention": "positive_definite",
    },
)


def solve_kepler(mean_anomaly: float, eccentricity: float, *, tolerance: float = 1e-13) -> float:
    """Solve ``M = E - e sin E`` for the eccentric anomaly (EQ-NUM-KEPLER-NEWTON).

    Newton-Raphson from ``E_0 = M``, which converges for the eccentricities the
    equation's validity domain admits.  The iteration count is capped so a
    configuration outside that domain fails rather than spinning.
    """
    if not 0.0 <= eccentricity < 0.9:
        raise ConfigurationError(
            f"eccentricity {eccentricity} is outside the declared validity domain [0, 0.9)"
        )
    E = mean_anomaly
    for _ in range(64):
        f = E - eccentricity * math.sin(E) - mean_anomaly
        if abs(f) < tolerance:
            return E
        E -= f / (1.0 - eccentricity * math.cos(E))
    raise ConfigurationError(
        f"Kepler solve failed to converge for M={mean_anomaly}, e={eccentricity}"
    )


def orbit_radius(
    time_s: float,
    *,
    semi_major_axis_m: float = _AU,
    eccentricity: float = 0.0,
    period_s: float = _SIDEREAL_YEAR_S,
    perihelion_time_s: float = 0.0,
) -> float:
    """Sun-Earth distance from Kepler's equation (EQ-KEPLER-ORBIT-RADIUS)."""
    mean_anomaly = 2.0 * math.pi * (time_s - perihelion_time_s) / period_s
    E = solve_kepler(mean_anomaly, eccentricity)
    return semi_major_axis_m * (1.0 - eccentricity * math.cos(E))


def spherical_mean_factor() -> float:
    """Projected disc over total surface for a sphere: exactly one quarter.

    ``pi R^2 / (4 pi R^2)`` -- independent of R, which is the content of ORB-002.
    """
    return 0.25


class OrbitalGeometryModule(Module):
    """Maps solar irradiance onto the planetary surface and owns the TOA boundary."""

    def __init__(
        self,
        *,
        module_id: str = "orbit",
        receiver_module_id: str = "surface",
        eccentricity: float = 0.0,
        semi_major_axis_m: float = _AU,
        period_s: float = _SIDEREAL_YEAR_S,
        perihelion_time_s: float = 0.0,
        status: ModuleStatus = ModuleStatus.ACTIVE,
    ) -> None:
        contract = ModuleContract(
            module_id=module_id,
            module_name="Keplerian geometry with spherical-mean insolation",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.REDUCED,
            physical_basis=(
                "Two-body Keplerian orbit, inverse-square irradiance scaling, and the "
                "geometric ratio of a sphere's projected disc to its total surface."
            ),
            governing_equation_ids=(
                EQUATION_IDS["kepler"],
                EQUATION_IDS["kepler_solver"],
                EQUATION_IDS["inverse_square"],
                EQUATION_IDS["insolation"],
            ),
            assumptions=(
                "spherical planet",
                "orbital elements fixed over the run",
                "instantaneous global averaging of insolation",
                "no atmosphere between the boundary and the receiving surface",
            ),
            validity_domain=(
                "Zero-dimensional global-mean configurations. Provides no diurnal, seasonal "
                "or latitudinal structure; a claim requiring any of those is outside scope."
            ),
            required_inputs=("solar_total_irradiance_1au", "solar_spectral_irradiance_1au"),
            produced_outputs=(
                "earth_sun_distance",
                "toa_incident_irradiance",
                "toa_spectral_irradiance",
                "insolation_geometry_factor",
            ),
            exported_transfers=("toa_incident_energy",),
            conserved_quantities=("energy",),
            internal_sources_and_sinks=("solar_input_across_system_boundary",),
            provided_interfaces=(ORBITAL_GEOMETRY_INTERFACE,),
            required_interfaces=(("solar_source", VersionRange.parse(">=1.0.0,<2.0.0")),),
            uncertainty_outputs=("orbital_element_uncertainty", "geometry_structural"),
            correctness_tests=(
                "test_intercepted_equals_distributed_power",
                "test_inverse_square_scaling",
                "test_kepler_circular_limit",
                "test_orbit_energy_passthrough_closes",
            ),
            known_limitations=(
                "Global-mean only; obliquity, rotation and tidal terms are not implemented.",
            ),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(n_bands) per step plus O(1) Kepler iterations",
                memory_complexity="O(n_bands)",
                scaling_variables=("n_bands",),
                work_units_per_step=2.0,
                benchmark_configuration="global scalar, any band count",
                peak_memory_bytes=8192,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
            ),
        )
        super().__init__(contract, status=status)
        self.receiver_module_id = receiver_module_id
        self.eccentricity = eccentricity
        self.semi_major_axis_m = semi_major_axis_m
        self.period_s = period_s
        self.perihelion_time_s = perihelion_time_s

        self._distance = semi_major_axis_m
        self._toa_irradiance = 0.0
        self._toa_spectral = np.zeros(1)
        self._grid: SpectralGrid | None = None
        self._exported_energy = 0.0
        self._interval = (0.0, 0.0)

    def validate_configuration(self) -> None:
        super().validate_configuration()
        if self.period_s <= 0:
            raise ConfigurationError(f"orbital period {self.period_s} must be positive")
        solve_kepler(0.0, self.eccentricity)  # domain check, cheap and eager

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        self._distance = orbit_radius(
            0.0,
            semi_major_axis_m=self.semi_major_axis_m,
            eccentricity=self.eccentricity,
            period_s=self.period_s,
            perihelion_time_s=self.perihelion_time_s,
        )
        self.mark_initialized()

    def compute_tendencies(self) -> list[Any]:
        ctx = self.context
        total = self.input_field("solar_total_irradiance_1au").require_units("W m^-2").scalar()
        spectral_field = self.input_field("solar_spectral_irradiance_1au").require_units(
            "W m^-2 m^-1"
        )
        grid = spectral_field.metadata.get("spectral_grid")
        if not isinstance(grid, SpectralGrid):
            raise ConfigurationError(
                "solar_spectral_irradiance_1au arrived without its spectral grid; a spectral "
                "field must carry the grid it is defined on (SPC-004)"
            )
        self._grid = grid

        self._distance = orbit_radius(
            ctx.time,
            semi_major_axis_m=self.semi_major_axis_m,
            eccentricity=self.eccentricity,
            period_s=self.period_s,
            perihelion_time_s=self.perihelion_time_s,
        )
        distance_factor = (_AU / self._distance) ** 2
        geometry = spherical_mean_factor()

        self._toa_irradiance = total * distance_factor * geometry
        self._toa_spectral = spectral_field.values * distance_factor * geometry
        self._exported_energy = self._toa_irradiance * ctx.dt
        self._interval = (ctx.time, ctx.time + ctx.dt)
        return []

    def emit_outputs(self) -> dict[str, Field]:
        ctx = self.context
        provenance = Provenance(
            producer=self.module_id,
            derivation=EQUATION_IDS["insolation"],
            notes="global-mean top-of-atmosphere mapping",
        )
        common = {
            "quantity_kind": QuantityKind.BOUNDARY_VALUE,
            "time_semantics": TimeSemantics.INTERVAL_AVERAGE,
            "time_interval": ctx.dt,
            "timestamp": ctx.time,
            "spatial_support": SpatialSupport.GLOBAL_SCALAR,
            "sign_convention": SignConvention.POSITIVE_DEFINITE,
            "provenance": provenance,
        }
        outputs = {
            "earth_sun_distance": Field(
                field_id="earth_sun_distance",
                physical_name="Earth-Sun distance",
                symbol="d",
                values=np.array(self._distance),
                units="m",
                quantity_kind=QuantityKind.DIAGNOSTIC,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                provenance=Provenance(producer=self.module_id, derivation=EQUATION_IDS["kepler"]),
                uncertainty=Uncertainty(
                    method="none", unresolved=("orbital_element_drift", "perturbations")
                ),
            ),
            "insolation_geometry_factor": Field(
                field_id="insolation_geometry_factor",
                physical_name="projected disc over total surface area",
                symbol="f",
                values=np.array(spherical_mean_factor()),
                units="1",
                **common,  # type: ignore[arg-type]
            ),
            "toa_incident_irradiance": Field(
                field_id="toa_incident_irradiance",
                physical_name="global-mean top-of-atmosphere incident irradiance",
                symbol="S_toa",
                values=np.array(self._toa_irradiance),
                units="W m^-2",
                **common,  # type: ignore[arg-type]
            ),
            "toa_spectral_irradiance": Field(
                field_id="toa_spectral_irradiance",
                physical_name="global-mean top-of-atmosphere spectral irradiance",
                symbol="S_toa_lambda",
                values=self._toa_spectral,
                units="W m^-2 m^-1",
                dimensions=("band",),
                spectral_coordinate=self._grid.grid_id if self._grid else None,
                metadata={"spectral_grid": self._grid},
                **common,  # type: ignore[arg-type]
            ),
        }
        return outputs

    def report_transfers(self) -> tuple[list[Transfer], list[Transfer]]:
        start, end = self._interval
        return (
            [
                Transfer(
                    transfer_id=f"toa_incident@{start:g}",
                    sender_module=self.module_id,
                    receiver_module=self.receiver_module_id,
                    quantity=ConservedQuantity.ENERGY,
                    amount=self._exported_energy,
                    units="J m^-2",
                    start_time=start,
                    end_time=end,
                    spatial_mapping="global_mean_over_planetary_surface",
                    channel="incident_shortwave",
                    equation_id=EQUATION_IDS["insolation"],
                    source_or_sink_reference="solar_input_across_system_boundary",
                )
            ],
            [],
        )

    def report_conservation(self) -> list[ModuleBudget]:
        budget = ModuleBudget(
            module_id=self.module_id,
            quantity=ConservedQuantity.ENERGY,
            units="J m^-2",
            stored_initial=0.0,
            stored_final=0.0,
        )
        # The Sun is outside the modeled boundary: its energy enters here as a
        # declared external source, exactly balancing what is exported inward.
        budget.add_source("solar_input_across_system_boundary", self._exported_energy)
        return [budget]

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "eccentricity": self.eccentricity,
            "semi_major_axis_m": self.semi_major_axis_m,
            "orbital_period_s": self.period_s,
            "perihelion_time_s": self.perihelion_time_s,
            "geometry_factor": spherical_mean_factor(),
            "receiver_module_id": self.receiver_module_id,
        }

    def checkpoint(self) -> dict[str, Any]:
        return {"distance": self._distance, "toa_irradiance": self._toa_irradiance}

    def restore(self, checkpoint: dict[str, Any]) -> None:
        self._distance = float(checkpoint["distance"])
        self._toa_irradiance = float(checkpoint["toa_irradiance"])
