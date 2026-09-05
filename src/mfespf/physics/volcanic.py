"""Stratospheric volcanic aerosol as a scattering layer.

Charter 9.8 and 9.14; requirements ATM-002, RAD-001, SPC-001, HYP-005, EQR-002.

The third mechanism to enter this framework, and the first that acts by
*scattering* rather than by absorbing. A large eruption puts sulfate aerosol into
the stratosphere, where it reflects sunlight back to space for a year or two
before settling out. The two-stream solver has handled scattering since Stage 4
and needed no change: the aerosol contributes an optical depth and an asymmetry
in the same layers a cloud would, and ``combine_scatterers`` adds the two.

**The record is an optical depth, so almost nothing has to be assumed.** Unlike
carbon dioxide, whose absorption coefficient had to be constructed and
calibrated, the volcanic record *is* the optical property -- 550 nm
stratospheric aerosol optical depth, measured or reconstructed. What is assumed
here is only how to spread it across wavelength and where to put it vertically,
and both are declared.

**Shortwave only.** Sulfate aerosol also absorbs in the longwave and warms the
stratosphere, which is why a real eruption cools the surface and warms the layer
the aerosol is in. That longwave effect is omitted: it is small in the
global-mean surface budget and this column has no stratosphere worth the name.
A claim about stratospheric temperature after an eruption is therefore not in
scope for any configuration using this.
"""

from __future__ import annotations

from collections.abc import Callable
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
from mfespf.contracts.spectral import SpectralGrid
from mfespf.core.errors import ConfigurationError
from mfespf.core.quantities import (
    Field,
    Provenance,
    QuantityKind,
    SignConvention,
    SpatialSupport,
    TimeSemantics,
)
from mfespf.core.versioning import Version
from mfespf.physics.equations import EQUATION_IDS

__all__ = [
    "REFERENCE_WAVELENGTH_M",
    "SULFATE_ANGSTROM_EXPONENT",
    "SULFATE_ASYMMETRY",
    "PrescribedVolcanicModule",
    "stratospheric_aerosol_optics",
]

V1 = Version(1, 0, 0)

REFERENCE_WAVELENGTH_M = 550.0e-9
"""The wavelength the record's optical depth is quoted at."""

SULFATE_ANGSTROM_EXPONENT = 1.0
"""How fast extinction falls with wavelength beyond the reference, as
``(lambda/lambda_0)^-alpha``, capped at the reference value below it.

Around one for stratospheric sulfate at the sizes an eruption produces, which is
smaller than the two or more typical of fine tropospheric aerosol because
volcanic droplets grow. A parameterization, not a measurement: the record gives
one number at one wavelength and the spectral shape has to come from somewhere.
"""

SULFATE_ASYMMETRY = 0.6
"""Scattering asymmetry of stratospheric sulfate in the shortwave.

Forward-peaked but far less so than a cloud droplet's 0.85, because the
particles are smaller. It matters: the fraction scattered back to space is
roughly ``(1-g)/2``, so an aerosol at 0.85 would reflect less than half as much
as one at 0.6 for the same optical depth.
"""


def stratospheric_aerosol_optics(
    grid: SpectralGrid,
    optical_depth: float,
    n_layers: int,
    *,
    stratospheric_layers: int = 1,
    angstrom_exponent: float = SULFATE_ANGSTROM_EXPONENT,
    asymmetry: float = SULFATE_ASYMMETRY,
    shortwave_cutoff_m: float = 4.0e-6,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Absorption, scattering and asymmetry arrays for a volcanic aerosol layer.

    Returns ``(absorption, scattering, asymmetry)`` shaped ``(n_layers, n_bands)``
    so they can be combined with any other scatterer in the same layers.

    The absorption array is identically zero: stratospheric sulfate is very
    nearly a pure scatterer in the shortwave, with a single-scattering albedo
    above 0.99, and representing it as partly absorbing would put energy into
    the layer that in reality leaves the planet.

    The aerosol is placed in the topmost ``stratospheric_layers`` layers. That is
    a crude vertical placement in a column whose top layer is not really a
    stratosphere, and the consequence is declared rather than tuned: the height
    matters for the longwave effect, which is omitted, and much less for the
    shortwave reflection, which is what this represents.
    """
    if optical_depth < 0.0:
        raise ConfigurationError("an optical depth cannot be negative")
    if not 0 < stratospheric_layers <= n_layers:
        raise ConfigurationError(
            f"cannot put aerosol in {stratospheric_layers} of {n_layers} layers"
        )
    if not 0.0 <= asymmetry < 1.0:
        raise ConfigurationError(f"asymmetry {asymmetry} must lie in [0, 1)")

    shape = (n_layers, grid.n_bins)
    scattering = np.zeros(shape)
    shortwave = grid.centers < shortwave_cutoff_m
    # The power law is capped at the reference value rather than continued to
    # short wavelengths. Post-eruption sulfate droplets are large compared with
    # visible light -- a size parameter around six at 550 nm -- so their
    # extinction efficiency has already reached the geometric-optics limit and
    # stops rising as the wavelength falls. An uncapped power law would put two
    # or three times the reference optical depth into the blue and ultraviolet,
    # where there is a great deal of solar energy, and the modelled cooling
    # after an eruption came out about twice what the record implies.
    spectral = np.where(
        shortwave,
        np.minimum((grid.centers / REFERENCE_WAVELENGTH_M) ** (-angstrom_exponent), 1.0),
        0.0,
    )
    # Split evenly between the stratospheric layers: the record gives a column
    # total, and nothing in it says how the aerosol is distributed within the
    # stratosphere.
    per_layer = optical_depth / stratospheric_layers
    scattering[n_layers - stratospheric_layers :, :] = per_layer * spectral
    return (
        np.zeros(shape),
        scattering,
        np.where(scattering > 0.0, asymmetry, 0.0),
    )


class PrescribedVolcanicModule(Module):
    """Stratospheric aerosol optical depth as a declared boundary condition.

    PRESCRIBED, owning no state: eruptions are an input to this framework and
    nothing in it explains when they happen. What follows from them is what is
    being tested.
    """

    def __init__(
        self,
        *,
        module_id: str = "volcanic",
        radiation_module_id: str = "radiation",
        optical_depth: float | Callable[[float], float] = 0.0,
        stratospheric_layers: int = 1,
        status: ModuleStatus = ModuleStatus.PRESCRIBED,
    ) -> None:
        contract = ModuleContract(
            module_id=module_id,
            module_name="prescribed stratospheric volcanic aerosol",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.PRESCRIBED,
            physical_basis=(
                "Sulfate aerosol injected into the stratosphere scatters sunlight back "
                "to space. The record supplies the optical depth directly, so the only "
                "constructed quantities are its spectral shape and its vertical "
                "placement."
            ),
            governing_equation_ids=(EQUATION_IDS["volcanic_aerosol_optics"],),
            assumptions=(
                "a pure scatterer in the shortwave, single-scattering albedo taken as one",
                "no longwave effect, so no stratospheric warming after an eruption",
                "spread evenly through the declared stratospheric layers",
                "globally uniform, so no hemispheric asymmetry after a one-sided eruption",
            ),
            validity_domain=(
                "Global-mean shortwave reflection by stratospheric aerosol. Nothing about "
                "stratospheric temperature, about the hemispheric asymmetry of a volcanic "
                "cloud, or about tropospheric aerosol is in scope."
            ),
            produced_outputs=("stratospheric_aerosol_optical_depth",),
            conserved_quantities=(),
            uncertainty_outputs=("prescribed_optical_depth_boundary_condition",),
            correctness_tests=(
                "test_zero_optical_depth_scatters_nothing",
                "test_the_aerosol_is_confined_to_the_shortwave_and_the_stratosphere",
                "test_a_volcanic_layer_reflects_and_cools",
            ),
            known_limitations=(
                "Prescribed, so eruptions are an input; nothing here explains their timing "
                "or size.",
                "Shortwave only: the longwave absorption that warms the real stratosphere "
                "after an eruption is absent, so no claim about stratospheric temperature "
                "is admissible.",
                "Globally uniform, so the hemispheric asymmetry of a one-sided eruption "
                "cannot be represented even though the record carries it.",
            ),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(1)",
                memory_complexity="O(1)",
                scaling_variables=(),
                work_units_per_step=0.5,
                benchmark_configuration="global scalar",
                peak_memory_bytes=4096,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
                lower_fidelity_alternatives=("no_volcanic_aerosol",),
            ),
        )
        super().__init__(contract, status=status)
        self.radiation_module_id = radiation_module_id
        self.optical_depth = optical_depth
        self.stratospheric_layers = stratospheric_layers
        self._tau = 0.0

    def validate_configuration(self) -> None:
        super().validate_configuration()
        if not callable(self.optical_depth) and float(self.optical_depth) < 0.0:
            raise ConfigurationError("an optical depth cannot be negative")
        if self.stratospheric_layers < 1:
            raise ConfigurationError("the aerosol needs at least one layer to sit in")

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        self._tau = self._at(float(boundary_state.get("master_start_time", 0.0)))
        self.mark_initialized()

    def _at(self, time: float) -> float:
        value = self.optical_depth(time) if callable(self.optical_depth) else self.optical_depth
        return float(value)

    def compute_tendencies(self) -> list[Any]:
        self._tau = self._at(self.context.time)
        return []

    def _field(self, timestamp: float) -> dict[str, Field]:
        return {
            "stratospheric_aerosol_optical_depth": Field(
                field_id="stratospheric_aerosol_optical_depth",
                physical_name="stratospheric aerosol optical depth at 550 nm",
                symbol="tau_v",
                values=np.array(self._tau),
                units="1",
                quantity_kind=QuantityKind.BOUNDARY_VALUE,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                timestamp=timestamp,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                provenance=Provenance(
                    producer=self.module_id,
                    derivation=EQUATION_IDS["volcanic_aerosol_optics"],
                    notes="prescribed boundary condition; eruptions are an input",
                ),
            )
        }

    def emit_state(self) -> dict[str, Field]:
        return self._field(self._context.time if self._context else 0.0)

    def emit_outputs(self) -> dict[str, Field]:
        return self._field(self.context.time)

    def report_uncertainty(self) -> UncertaintyReport:
        return UncertaintyReport(
            module_id=self.module_id,
            method="none; the optical depth is imposed",
            distributional_assumptions="unresolved",
            correlations="unresolved",
            unresolved=(
                "the record publishes no per-sample uncertainty, and the pre-satellite "
                "part is reconstructed rather than measured",
                "structural_uncertainty_of_a_shortwave_only_globally_uniform_aerosol",
            ),
        )

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "optical_depth": (
                "time_varying" if callable(self.optical_depth) else float(self.optical_depth)
            ),
            "stratospheric_layers": self.stratospheric_layers,
            "reference_wavelength_m": REFERENCE_WAVELENGTH_M,
            "angstrom_exponent": SULFATE_ANGSTROM_EXPONENT,
            "asymmetry": SULFATE_ASYMMETRY,
            "longwave_effect": "omitted",
        }
