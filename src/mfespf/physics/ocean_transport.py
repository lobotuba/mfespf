"""The ocean's heat transport between a rising and a sinking region.

Why this exists
---------------

Two identical columns under the relaxed weak-temperature-gradient balance
never make a temperature difference, so they never overturn
(``docs/protocols/spectral-wtg-circulation.md``, Test A). The literature that
balance follows imposes a sea-surface temperature contrast. Here the contrast
is not imposed: the ocean carries heat from under the sinking air to under the
rising air, as it is observed to within 20 degrees of the equator, and the
surface temperatures that result are the columns' own.

The physics
-----------

A prescribed lateral heat flux into each region's surface, the slab-ocean
"q-flux" (Hansen et al. 1984; Bitz et al. 2012). The rising region's surface
receives ``Q_r = (1 - sigma) dQ`` and the sinking region's ``Q_s = -sigma dQ``,
each in W m^-2 of its own region, so ``sigma Q_r + (1 - sigma) Q_s = 0`` exactly:
the ocean moves heat and makes none. ``dQ`` is the observed contrast in the heat
the ocean delivers, from annual-mean net surface fluxes composited on
annual-mean omega500 masks (``docs/protocols/ocean-heat-transport-circulation.md``).

What crosses is recorded as one transfer per step, sinking surface to rising
surface, per unit of global area, each budget taking its own area share -- the
same bookkeeping as the circulation's, because it is this module, not either
surface, that moves the heat.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from mfespf.contracts.cost import CostDeclaration, ResourceClass
from mfespf.contracts.module import (
    FidelityLevel,
    Module,
    ModuleContract,
    ModuleStatus,
    UncertaintyReport,
)
from mfespf.core.errors import ConfigurationError
from mfespf.core.quantities import (
    Field,
    Provenance,
    QuantityKind,
    SpatialSupport,
    TimeSemantics,
)
from mfespf.core.versioning import Version
from mfespf.coupling.conservation import ModuleBudget
from mfespf.coupling.transfers import ConservedQuantity, Transfer
from mfespf.physics.equations import EQUATION_IDS
from mfespf.physics.surface import SURFACE_TEMPERATURE_STATE
from mfespf.registries.states import TendencyContribution

__all__ = [
    "OBSERVED_OCEAN_HEAT_CONTRAST_W_M2",
    "OceanHeatTransport",
    "ocean_heat_delivered",
]

V1 = Version(1, 0, 0)

#: The contrast in heat the ocean delivers to the surface, rising minus sinking
#: region, W m^-2, within 20 degrees of the equator: annual-mean net downward
#: surface flux composited on annual-mean omega500 masks, each region's mean less
#: the two regions' area-weighted mean (ocean-heat-transport-circulation.md).
#: OAFlux 1991-2009 (Yu and Weller 2007), land delivering none, is the verdict
#: value; NCEP/NCAR R1 1991-2020, all surface, the sensitivity.
OBSERVED_OCEAN_HEAT_CONTRAST_W_M2: dict[str, float] = {"oaflux": 13.39, "ncep_r1": 4.67}


def ocean_heat_delivered(contrast: float, ascent_fraction: float) -> tuple[float, float]:
    """Heat into the rising and the sinking surface, W m^-2 of each region.

    They sum to zero over the area: the ocean only moves heat.
    """
    if not 0.0 < ascent_fraction < 1.0:
        raise ConfigurationError("the rising region's area fraction must lie in (0, 1)")
    return (1.0 - ascent_fraction) * contrast, -ascent_fraction * contrast


class OceanHeatTransport(Module):
    """Carries heat from the sinking region's surface to the rising region's
    (EQ-OCEAN-HEAT-TRANSPORT)."""

    def __init__(
        self,
        *,
        ascent_fraction: float,
        heat_contrast: float,
        module_id: str = "ocean_transport",
        ascent_region: str = "ascent",
        descent_region: str = "descent",
        surface_module: str = "surface",
        status: ModuleStatus = ModuleStatus.ACTIVE,
    ) -> None:
        regions = (ascent_region, descent_region)
        contract = ModuleContract(
            module_id=module_id,
            module_name="prescribed ocean heat transport between two regions",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.REDUCED,
            physical_basis=(
                "Ocean currents carry heat laterally beneath the atmosphere; in a slab "
                "ocean that transport enters each region's mixed layer as a prescribed "
                "flux whose area-weighted sum is zero."
            ),
            governing_equation_ids=(EQUATION_IDS["ocean_heat_transport"],),
            assumptions=(
                "one rising and one sinking region partition the domain",
                "the transport is fixed at its observed annual mean",
            ),
            validity_domain=(
                "Equilibrium experiments in which the ocean's lateral heat transport "
                "between the regions may be taken as fixed; it does not respond to the "
                "surface temperatures it creates."
            ),
            required_inputs=tuple(f"{r}:surface_heat_capacity" for r in regions),
            produced_outputs=tuple(f"{r}:ocean_heat_convergence" for r in regions),
            conserved_quantities=("energy",),
            internal_sources_and_sinks=(),
            uncertainty_outputs=("ocean_heat_contrast_observational",),
            correctness_tests=(
                "test_the_ocean_moves_heat_and_makes_none",
                "test_each_surface_warms_at_its_flux_over_its_heat_capacity",
                "test_no_contrast_leaves_the_two_columns_identical",
            ),
            known_limitations=(
                "The transport is prescribed, not dynamic.",
                "Two regions stand for the band's distribution of surface fluxes.",
            ),
            fidelity_justification=(
                "Without it two identical columns never overturn "
                "(spectral-wtg-circulation.md Test A); a dynamic ocean would need ocean "
                "physics this framework does not have."
            ),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(1)",
                memory_complexity="O(1)",
                scaling_variables=(),
                work_units_per_step=0.01,
                benchmark_configuration="two regions",
                peak_memory_bytes=4096,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
                lower_fidelity_alternatives=("no_ocean_transport",),
            ),
        )
        super().__init__(contract, status=status)
        self.ascent_region = ascent_region
        self.descent_region = descent_region
        self.ascent_fraction = float(ascent_fraction)
        self.heat_contrast = float(heat_contrast)
        self.surface_module = surface_module
        self.into_ascent, self.into_descent = ocean_heat_delivered(
            self.heat_contrast, self.ascent_fraction
        )
        self._transfers: list[Transfer] = []

    def validate_configuration(self) -> None:
        super().validate_configuration()
        if not 0.0 < self.ascent_fraction < 1.0:
            raise ConfigurationError("the rising region's area fraction must lie in (0, 1)")
        if not np.isfinite(self.heat_contrast):
            raise ConfigurationError("the ocean heat contrast must be finite")

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        self.mark_initialized()

    def _capacity(self, region: str) -> float:
        return (
            self.input_field(f"{region}:surface_heat_capacity")
            .require_units("J m^-2 K^-1")
            .scalar()
        )

    def compute_tendencies(self) -> list[TendencyContribution]:
        ctx = self.context
        # Per unit of global area: what leaves the sinking region's surface is
        # what reaches the rising region's.
        amount = self.ascent_fraction * self.into_ascent * ctx.dt
        sender, receiver = (
            (self.descent_region, self.ascent_region)
            if amount >= 0.0
            else (self.ascent_region, self.descent_region)
        )
        fraction = {
            self.ascent_region: self.ascent_fraction,
            self.descent_region: 1.0 - self.ascent_fraction,
        }
        self._transfers = [
            Transfer(
                transfer_id=f"{self.module_id}_heat@{ctx.time:g}",
                sender_module=f"{sender}/{self.surface_module}",
                receiver_module=f"{receiver}/{self.surface_module}",
                quantity=ConservedQuantity.ENERGY,
                amount=abs(amount),
                units="J m^-2",
                start_time=ctx.time,
                end_time=ctx.time + ctx.dt,
                spatial_mapping="area_weighted_regions",
                equation_id=EQUATION_IDS["ocean_heat_transport"],
                channel="ocean_heat_transport",
                sender_area_fraction=fraction[sender],
                receiver_area_fraction=fraction[receiver],
            )
        ]
        return [
            TendencyContribution(
                state_id=f"{region}:{SURFACE_TEMPERATURE_STATE}",
                provider_module_id=self.module_id,
                label="ocean_heat_transport",
                value=flux / self._capacity(region),
                units="K s^-1",
                equation_id=EQUATION_IDS["ocean_heat_transport"],
            )
            for region, flux in (
                (self.ascent_region, self.into_ascent),
                (self.descent_region, self.into_descent),
            )
        ]

    def report_transfers(self) -> tuple[list[Transfer], list[Transfer]]:
        return list(self._transfers), list(self._transfers)

    def report_conservation(self) -> list[ModuleBudget]:
        return []

    def emit_outputs(self) -> dict[str, Field]:
        ctx = self.context
        return {
            f"{region}:ocean_heat_convergence": Field(
                field_id=f"{region}:ocean_heat_convergence",
                physical_name="heat the ocean's lateral transport delivers to the surface",
                values=np.array(flux),
                units="W m^-2",
                quantity_kind=QuantityKind.DIAGNOSTIC,
                time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                time_interval=ctx.dt,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                provenance=Provenance(
                    producer=self.module_id, derivation=EQUATION_IDS["ocean_heat_transport"]
                ),
            )
            for region, flux in (
                (self.ascent_region, self.into_ascent),
                (self.descent_region, self.into_descent),
            )
        }

    def report_uncertainty(self) -> UncertaintyReport:
        return UncertaintyReport(
            module_id=self.module_id,
            quantified={"heat_contrast_W_m2": self.heat_contrast},
            method="none propagated; the two observed values bracket it",
            unresolved=(
                "OAFlux and NCEP R1 differ by a factor of 2.9 in the contrast",
                "the transport is fixed rather than responding to the surface",
            ),
        )

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "ascent_region": self.ascent_region,
            "descent_region": self.descent_region,
            "ascent_fraction": self.ascent_fraction,
            "heat_contrast_W_m2": self.heat_contrast,
            "into_ascent_W_m2": self.into_ascent,
            "into_descent_W_m2": self.into_descent,
            "heat_contrast_source": (
                "annual-mean net surface heat flux on annual-mean omega500 masks within 20 "
                "degrees: OAFlux 1991-2009 (13.39) or NCEP/NCAR R1 1991-2020 (4.67)"
            ),
        }
