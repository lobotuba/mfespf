"""An observed global cloud field: ISCCP's clouds at ISCCP's heights.

Why this exists
---------------

The column's 0.30 albedo stood on its surface, doing the clouds' job at the
ground. That was harmless while the atmosphere absorbed no sunlight and wrong the
moment it did: reflected light crossed the absorbing vapour twice and the
planetary albedo fell to 0.22 (``docs/protocols/solar-vapour-absorption.md``).
Earth reflects most of its sunlight from cloud tops, above most of its vapour.

What it is
----------

Rossow and Schiffer (1999), *Advances in Understanding Clouds from ISCCP*,
BAMS 80, 2261-2287: Table 4's D-series global annual means for the four cloud
levels -- amount, top pressure and visible optical thickness -- with Table 5's
split of each level into ice and liquid, and section 3's infrared relations,
``tau_ir = tau / 2.56`` for liquid and ``tau / 2.13`` for ice.

ISCCP sees clouds from above, so its amounts partition the sky: each type is a
sub-column and what remains is clear. Its optical thicknesses are radiatively
weighted, essentially logarithmic, averages -- chosen by its authors to represent
the area-averaged radiative effect, which is the quantity a single column needs.
Nothing here is fitted; the test is CERES (``docs/protocols/observed-cloud-field.md``).

The field carries no water. It is optical properties at observed places, so it
enters no budget, and it responds to nothing: no cloud feedback is in scope.
"""

from __future__ import annotations

from dataclasses import dataclass
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
from mfespf.core.versioning import InterfaceContract, Version
from mfespf.physics.equations import EQUATION_IDS

__all__ = [
    "ISCCP_CLOUD_TYPES",
    "ISCCP_TROPICAL_CLOUD_TYPES",
    "CloudType",
    "ObservedCloudFieldModule",
    "cloud_type_layers",
]

V1 = Version(1, 0, 0)

#: ISCCP's infrared-to-visible optical thickness ratios (Rossow and Schiffer
#: 1999, section 3): Mie calculations for 10-micron droplets, and Minnis et al.
#: (1993a) for 30-micron ice.
TAU_IR_PER_TAU_LIQUID = 1.0 / 2.56
TAU_IR_PER_TAU_ICE = 1.0 / 2.13

#: Asymmetry parameter for scattered sunlight: strongly forward, delta-scaled by
#: the solver. One value for both phases, declared as an approximation.
CLOUD_ASYMMETRY = 0.85

#: ISCCP's D-series global annual mean surface visible reflectance (Table 4).
ISCCP_SURFACE_REFLECTANCE = 0.17


@dataclass(frozen=True)
class CloudType:
    """One ISCCP cloud level, as seen from above."""

    name: str
    amount: float
    top_pressure_pa: float
    visible_optical_thickness: float
    ice_fraction: float

    @property
    def infrared_optical_thickness(self) -> float:
        """Absorption optical depth for thermal radiation, phase-weighted."""
        per_tau = (
            self.ice_fraction * TAU_IR_PER_TAU_ICE
            + (1.0 - self.ice_fraction) * TAU_IR_PER_TAU_LIQUID
        )
        return self.visible_optical_thickness * per_tau


#: Table 4 (D2, global annual mean) for amounts, top pressures and optical
#: thicknesses. Ice fractions: high and deep convective clouds are ice (tops at
#: 227.5 and 234.8 K, below ISCCP's 260 K phase boundary); middle and low levels
#: take Table 5's global ice and liquid amounts -- middle 9.9 of 19.2 per cent
#: ice, low 2.3 of 26.5.
ISCCP_CLOUD_TYPES: tuple[CloudType, ...] = (
    CloudType("cirrus_and_cirrostratus", 0.196, 26700.0, 2.2, 1.0),
    CloudType("deep_convective", 0.026, 32600.0, 35.6, 1.0),
    CloudType("middle", 0.190, 55700.0, 4.8, 9.9 / 19.2),
    CloudType("low", 0.275, 82600.0, 4.7, 2.3 / 26.5),
)

#: The tropics (+/-30 degrees), from Table 5's daytime cloud-type amounts, ice
#: (liquid): high is cirrus 15.6 + cirrostratus 5.5; deep convective 2.7; middle
#: is altocumulus 0.9 (6.1) + altostratus 0.1 (4.7) + nimbostratus 0.0 (1.0);
#: low is cumulus 0.0 (12.3) + stratocumulus 0.0 (10.7) + stratus 0.0 (0.8). Table
#: 5 gives no tropical top pressures or optical thicknesses, so each level keeps
#: the global Table 4 values above -- declared, since tropical high cloud is
#: observed higher than the global mean (docs/protocols/tropical-cloud-field.md).
ISCCP_TROPICAL_CLOUD_TYPES: tuple[CloudType, ...] = (
    CloudType("cirrus_and_cirrostratus", 0.211, 26700.0, 2.2, 1.0),
    CloudType("deep_convective", 0.027, 32600.0, 35.6, 1.0),
    CloudType("middle", 0.128, 55700.0, 4.8, 1.0 / 12.8),
    CloudType("low", 0.238, 82600.0, 4.7, 0.0),
)

CLOUD_FIELD_INTERFACE = InterfaceContract(
    interface_id="cloud_field",
    version=V1,
    physical_semantics={
        "partition": "types_are_disjoint_sub_columns_seen_from_above",
        "optical_thickness_basis": "radiatively_weighted_mean_visible",
        "infrared_optical_thickness": "absorption_only",
        "vertical_placement": "one_layer_containing_the_top_pressure",
    },
)


def cloud_type_layers(top_pressures: np.ndarray, interface_pressures: np.ndarray) -> np.ndarray:
    """Index of the layer containing each top pressure (0 is the lowest layer).

    Interfaces run from the surface upward, so layer ``i`` spans
    ``interface[i] >= p > interface[i + 1]``. A top outside the column is
    refused rather than moved to its nearest edge.
    """
    edges = np.asarray(interface_pressures, dtype=float)
    tops = np.asarray(top_pressures, dtype=float)
    if np.any(tops > edges[0]) or np.any(tops <= edges[-1]):
        raise ConfigurationError(
            f"cloud top pressures {tops.tolist()} lie outside the column "
            f"[{edges[-1]:g}, {edges[0]:g}] Pa"
        )
    # Descending edges: count how many interfaces lie at or above each top.
    return np.array([int(np.sum(edges >= p)) - 1 for p in tops])


class ObservedCloudFieldModule(Module):
    """Publishes ISCCP's global cloud field as sub-column optical properties."""

    def __init__(
        self,
        *,
        module_id: str = "cloud_field",
        cloud_types: tuple[CloudType, ...] = ISCCP_CLOUD_TYPES,
        status: ModuleStatus = ModuleStatus.PRESCRIBED,
    ) -> None:
        contract = ModuleContract(
            module_id=module_id,
            module_name="observed global cloud field (ISCCP D-series)",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.PRESCRIBED,
            physical_basis=(
                "Clouds reflect sunlight from their tops and emit from them at the "
                "temperature there. Their observed amounts, top pressures and optical "
                "thicknesses place that reflection and emission where Earth's are."
            ),
            governing_equation_ids=(
                EQUATION_IDS["observed_cloud_field"],
                EQUATION_IDS["cloud_radiative_effect"],
            ),
            # Bound explicitly. Left to inference, the record's optical-thickness
            # range was matched by symbol to the radiation module's *total* layer
            # optical depth -- saturated vapour bands at 1600 -- and invalidated
            # a run over a cloud that never exceeded 36.
            equation_variables=(
                (
                    EQUATION_IDS["observed_cloud_field"],
                    (("tau_c", "cloud_type_visible_optical_thickness", "1"),),
                ),
            ),
            assumptions=(
                "the cloud field is prescribed from a global climatology and responds "
                "to nothing in the model",
                "ISCCP's amounts partition the sky into disjoint sub-columns",
                "each cloud type occupies the one layer containing its top pressure",
                "ice and liquid in one level share that level's optical thickness",
            ),
            validity_domain=(
                "Radiative experiments on a global-mean column in which the cloud "
                "field is an imposed boundary condition. No cloud feedback, no "
                "regional or seasonal cloud, no multi-layer cloud."
            ),
            required_inputs=(),
            produced_outputs=(
                "cloud_type_amount",
                "cloud_type_top_pressure",
                "cloud_type_visible_optical_thickness",
                "cloud_type_infrared_optical_thickness",
            ),
            conserved_quantities=(),
            internal_sources_and_sinks=(),
            provided_interfaces=(CLOUD_FIELD_INTERFACE,),
            uncertainty_outputs=(
                "logarithmic_mean_optical_thickness_structural",
                "no_multilayer_cloud_structural",
                "visible_surface_reflectance_used_for_all_solar_bands",
            ),
            correctness_tests=(
                "test_the_cloud_types_partition_the_sky",
                "test_each_cloud_sits_in_the_layer_containing_its_top",
                "test_the_infrared_optical_thickness_follows_the_phase",
                "test_an_empty_field_is_the_clear_column",
            ),
            known_limitations=(
                "Prescribed, so no cloud feedback of any kind.",
                "Global annual means: no diurnal, seasonal or regional cloud.",
                "No cloud beneath a visible cloud; ISCCP cannot see it.",
                "One asymmetry parameter for liquid and ice.",
            ),
            fidelity_justification=(
                "The surface albedo that stood in for clouds cannot place reflection "
                "above the absorbing vapour, which is the defect this module exists "
                "to remove; it is the null alternative, not a cheaper version."
            ),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(n_types)",
                memory_complexity="O(n_types)",
                scaling_variables=("n_types",),
                work_units_per_step=0.01,
                benchmark_configuration=f"{len(cloud_types)} cloud types",
                peak_memory_bytes=4096,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
                lower_fidelity_alternatives=("surface_albedo_standing_in_for_clouds",),
            ),
        )
        super().__init__(contract, status=status)
        self.cloud_types = tuple(cloud_types)

    def validate_configuration(self) -> None:
        super().validate_configuration()
        amounts = [t.amount for t in self.cloud_types]
        if any(a < 0.0 for a in amounts) or sum(amounts) > 1.0:
            raise ConfigurationError(
                f"cloud amounts {amounts} must be non-negative and sum to at most one: "
                "seen from above, they partition the sky"
            )
        for t in self.cloud_types:
            if t.visible_optical_thickness < 0.0 or not 0.0 <= t.ice_fraction <= 1.0:
                raise ConfigurationError(f"cloud type {t.name!r} is unphysical")

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        self.mark_initialized()

    def compute_tendencies(self) -> list[Any]:
        return []

    def emit_outputs(self) -> dict[str, Field]:
        ctx = self.context
        provenance = Provenance(
            producer=self.module_id, derivation=EQUATION_IDS["observed_cloud_field"]
        )
        common = {
            "quantity_kind": QuantityKind.DIAGNOSTIC,
            "time_semantics": TimeSemantics.INSTANTANEOUS,
            "timestamp": ctx.time,
            "spatial_support": SpatialSupport.GLOBAL_SCALAR,
            "dimensions": ("cloud_type",),
            "provenance": provenance,
        }
        values = {
            "cloud_type_amount": ([t.amount for t in self.cloud_types], "1"),
            "cloud_type_top_pressure": ([t.top_pressure_pa for t in self.cloud_types], "Pa"),
            "cloud_type_visible_optical_thickness": (
                [t.visible_optical_thickness for t in self.cloud_types],
                "1",
            ),
            "cloud_type_infrared_optical_thickness": (
                [t.infrared_optical_thickness for t in self.cloud_types],
                "1",
            ),
        }
        return {
            name: Field(
                field_id=name,
                physical_name=name.replace("_", " "),
                values=np.asarray(v, dtype=float),
                units=units,
                **common,  # type: ignore[arg-type]
            )
            for name, (v, units) in values.items()
        }

    def report_uncertainty(self) -> UncertaintyReport:
        return UncertaintyReport(
            module_id=self.module_id,
            quantified={"clear_sky_amount": 1.0 - sum(t.amount for t in self.cloud_types)},
            method="none propagated; the field is an observed climatology",
            unresolved=(
                "the optical thicknesses are logarithmic means, so the spread of real "
                "clouds within a type is absent",
                "multi-layer cloud is absent",
                "ISCCP's visible surface reflectance stands for the whole solar spectrum",
            ),
        )

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "source": "Rossow and Schiffer 1999, BAMS 80, Tables 4 and 5",
            "cloud_types": [
                {
                    "name": t.name,
                    "amount": t.amount,
                    "top_pressure_Pa": t.top_pressure_pa,
                    "visible_optical_thickness": t.visible_optical_thickness,
                    "ice_fraction": t.ice_fraction,
                }
                for t in self.cloud_types
            ],
            "asymmetry": CLOUD_ASYMMETRY,
        }
