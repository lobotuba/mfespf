"""Carbon dioxide as a spectral absorber.

Charter 9.5 and 9.3; requirements ATM-002, CHM-001, EQR-002, EQR-004, RAD-001,
SPC-001.

CO2 enters this framework the way the charter says a mechanism should: through a
constitutive relation that changes an absorption coefficient, so the radiative
transfer computes the consequence. It does *not* enter as a fitted forcing added
to an energy budget. The difference matters, because a fitted forcing would make
the greenhouse effect an input to the model rather than a result of it, and the
whole question this framework asks is whether the physics produces the observed
behaviour rather than whether an assumed answer can be inserted.

**Where the logarithm comes from.** The near-constant forcing per doubling is not
put in by hand and is not a fitted constant. It is a consequence of the band's
shape: absorption is strongest at the band centre and falls exponentially into
the wings, so as concentration rises the band core saturates and the *width* over
which the atmosphere is opaque grows logarithmically. Represented on a grid with
enough bands across the wings, the transfer solver produces that behaviour by
itself, and ``forcing_per_doubling`` measures it rather than asserting it.

**What is calibrated, and against what.** One number is fitted: the reference
mass absorption coefficient at the band centre, chosen so that a doubling
produces the accepted instantaneous forcing. It is calibrated against a
*published forcing value*, never against a temperature record -- which is what
keeps a temperature record available as independent validation afterwards
(INF-003). The dataset registry enforces the distinction; this module's job is
to make the calibration target something other than the thing being predicted.
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
from mfespf.contracts.spectral import (
    IntegrationConvention,
    SpectralCoordinate,
    SpectralGrid,
)
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

V1 = Version(1, 0, 0)

__all__ = [
    "CO2_BAND_CENTRE_ABSORPTION",
    "CO2_BAND_CENTRE_M",
    "CO2_MOLAR_MASS",
    "CO2_WING_DECAY_PER_METRE",
    "DRY_AIR_MOLAR_MASS",
    "PREINDUSTRIAL_CO2_PPM",
    "PrescribedGreenhouseModule",
    "co2_band_absorption",
    "co2_mass_mixing_ratio",
    "greenhouse_grid",
]

CO2_MOLAR_MASS = 44.0095e-3
DRY_AIR_MOLAR_MASS = 28.9647e-3
CO2_BAND_CENTRE_M = 15.0e-6
"""The 15 micron bending-mode band, which is where CO2 does its work.

Near the peak of terrestrial emission at terrestrial temperatures and almost
clear of the solar spectrum, which is why carbon dioxide warms rather than
shades. The 4.3 micron band is stronger per molecule and sits where there is
little of either, so it is omitted and declared.
"""

CO2_BAND_CENTRE_ABSORPTION = 2.0e3
"""Mass absorption coefficient of CO2 at the band centre, m^2 per kg of CO2.

Not fitted. It is set to make the band core *deeply saturated* at present-day
concentrations -- an optical depth of order ten thousand through the column --
which is a qualitative fact about the real 15 micron band rather than a number
tuned to an answer. Once the core is saturated its exact value stops mattering:
the forcing comes from where the wings cross unit optical depth, and that
displacement depends on the logarithm of this number, not on the number.
"""

CO2_WING_DECAY_PER_METRE = 3.6e3
"""Wing decay scale in inverse metres, equivalent to 36 cm^-1.

The one number calibrated here, and it is calibrated against a *published
radiative-forcing value* -- 3.7 W m^-2 for a doubling -- never against a
temperature record. That is what keeps a temperature record available as
independent validation afterwards, and the dataset registry enforces the
distinction (INF-003, EQR-004).
"""

PREINDUSTRIAL_CO2_PPM = 284.7
"""Law Dome ice-core value for 1850, used as the reference for forcing.

A number from a record rather than a round figure: 280 is a convention, and a
forcing quoted against a convention is not traceable to anything.
"""


def co2_mass_mixing_ratio(ppm: float | np.ndarray) -> np.ndarray:
    """Mole fraction in parts per million to mass mixing ratio, kg per kg of air.

    The radiation module multiplies its absorption coefficient by the *total*
    layer mass, so a coefficient expressed per kilogram of CO2 has to be scaled
    by this before it is handed over. Doing that conversion here rather than
    inside the radiation module keeps the radiation module ignorant of what gas
    it is transporting through, which is what lets a second absorber be added
    without touching it.
    """
    fraction = np.asarray(ppm, dtype=float) * 1.0e-6
    if np.any(fraction < 0.0):
        raise ConfigurationError("a mole fraction cannot be negative")
    return fraction * (CO2_MOLAR_MASS / DRY_AIR_MOLAR_MASS)


def greenhouse_grid(
    n_background: int = 24,
    n_band: int = 24,
    *,
    lower_m: float = 1.0e-7,
    upper_m: float = 1.0e-3,
    band_lower_m: float = 1.15e-5,
    band_upper_m: float = 1.95e-5,
) -> SpectralGrid:
    """A log-spaced grid with extra resolution across the CO2 band.

    The default column grid puts two or three bands across 12 to 20 microns,
    which cannot resolve a band whose wings are the whole point: with the wings
    unresolved the band either saturates completely or scales linearly, and
    neither is what carbon dioxide does.

    The extra bands are inside the CO2 region and nowhere else, so the cost is
    paid where it buys something. Everything outside is the same logarithmic
    spacing the other configurations use.
    """
    if n_band < 8:
        raise ConfigurationError(
            f"{n_band} bands cannot resolve the wings of an absorption band; the "
            f"logarithmic dependence on concentration comes from the wings, and a grid "
            f"that does not resolve them produces linear growth or complete saturation "
            f"instead"
        )
    if not lower_m < band_lower_m < band_upper_m < upper_m:
        raise ConfigurationError("the CO2 band must lie inside the grid's support")

    background = np.logspace(np.log10(lower_m), np.log10(upper_m), n_background + 1)
    background = background[(background < band_lower_m) | (background > band_upper_m)]
    band = np.logspace(np.log10(band_lower_m), np.log10(band_upper_m), n_band + 1)
    edges = np.unique(np.concatenate([background, band, [lower_m, upper_m]]))
    return SpectralGrid(
        grid_id=f"greenhouse_{n_background}bg_{n_band}co2",
        coordinate=SpectralCoordinate.WAVELENGTH,
        coordinate_units="m",
        edges=edges,
        integration_convention=IntegrationConvention.SPECTRAL_DENSITY,
        density_units="W m^-2 m^-1",
        resolution_identifier=f"greenhouse-{edges.size - 1}band",
    )


def co2_band_absorption(
    grid: SpectralGrid,
    ppm: float,
    *,
    centre_absorption: float = CO2_BAND_CENTRE_ABSORPTION,
    wing_decay_per_metre: float = CO2_WING_DECAY_PER_METRE,
    background: float = 1.0e-4,
    band_centre_m: float = CO2_BAND_CENTRE_M,
    shortwave_cutoff_m: float = 4.0e-6,
) -> np.ndarray:
    """Absorption coefficient per kilogram of air, with a CO2 band in it.

    ``centre_absorption`` is the mass absorption coefficient of CO2 itself at
    the band centre, per kilogram of CO2. It is scaled here by the mass mixing
    ratio, so what the radiation module receives is per kilogram of air and can
    be multiplied by the layer mass it already holds.

    The wings fall exponentially in *wavenumber*, not in wavelength, because
    that is how a vibration-rotation band's line strengths fall away from the
    band centre. The choice is what produces the logarithmic dependence of
    forcing on concentration, and it produces it as a consequence rather than by
    assertion: the optical depth at a displacement ``dv`` from the centre is
    ``tau_0 exp(-dv / D)``, so the displacement at which the atmosphere stops
    being opaque is ``D ln(tau_0)``, and the opaque width therefore grows by
    ``2 D ln 2`` per doubling of concentration whatever the starting point.

    ``background`` stands for everything else the atmosphere absorbs in the
    infrared, water vapour above all. It is a free parameter of the reduced
    model rather than a property of any gas, and it is held fixed so that a
    difference between two runs is attributable to CO2 alone.
    """
    if centre_absorption < 0.0 or background < 0.0:
        raise ConfigurationError("absorption coefficients must be non-negative")
    if wing_decay_per_metre <= 0.0:
        raise ConfigurationError("the wing decay scale must be positive")

    centres = grid.centers
    infrared = centres >= shortwave_cutoff_m
    mixing_ratio = float(co2_mass_mixing_ratio(ppm))
    displacement = np.abs(1.0 / centres - 1.0 / band_centre_m)
    # The band is confined to the infrared exactly, not merely to within an
    # exponential tail. An envelope evaluated across the whole grid leaves a
    # residue of order 1e-23 in the solar bands: numerically nothing, and still
    # a greenhouse gas absorbing sunlight, which is the wrong sign of effect and
    # not something to leave to a tolerance.
    profile = np.where(infrared, np.exp(-displacement / wing_decay_per_metre), 0.0)
    return np.where(infrared, background, 0.0) + centre_absorption * mixing_ratio * profile


class PrescribedGreenhouseModule(Module):
    """Carbon dioxide as a declared boundary condition (CHM-001, HYP-005).

    PRESCRIBED in the same sense as the cloud module: it owns no state and
    responds to nothing. The concentration is imposed on the run, so this
    configuration has no carbon cycle and cannot be asked whether the rise is
    explained -- only what follows from it.

    It publishes a concentration and nothing else. Turning that into an
    absorption coefficient is the radiation module's business, so a second gas
    could be added here without the radiation module learning about it.
    """

    def __init__(
        self,
        *,
        module_id: str = "greenhouse",
        radiation_module_id: str = "radiation",
        concentration: float | Callable[[float], float] = PREINDUSTRIAL_CO2_PPM,
        status: ModuleStatus = ModuleStatus.PRESCRIBED,
    ) -> None:
        contract = ModuleContract(
            module_id=module_id,
            module_name="prescribed carbon dioxide concentration",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.PRESCRIBED,
            physical_basis=(
                "A well-mixed greenhouse gas concentration imposed as a boundary "
                "condition. No carbon cycle: the concentration is an input, and nothing "
                "in this configuration explains why it takes the values it does."
            ),
            governing_equation_ids=(EQUATION_IDS["co2_band_absorption"],),
            assumptions=(
                "carbon dioxide is uniformly mixed through the column",
                "only the 15 micron band is represented",
                "no other greenhouse gas is represented apart from a fixed background",
            ),
            validity_domain=(
                "Concentrations within roughly a factor of four of present day, where the "
                "band model's wings stay inside the resolved spectral range."
            ),
            produced_outputs=("co2_mole_fraction",),
            conserved_quantities=(),
            supported_spectral_representations=("any_grid_resolving_the_15_micron_band",),
            uncertainty_outputs=("prescribed_concentration_boundary_condition",),
            correctness_tests=(
                "test_forcing_per_doubling_is_nearly_constant",
                "test_a_constant_concentration_control_produces_no_trend",
                "test_zero_concentration_leaves_only_the_background",
            ),
            known_limitations=(
                "Prescribed, so there is no carbon cycle and no feedback of any kind "
                "between climate and concentration.",
                "One band, one gas. Methane, nitrous oxide and the water-vapour "
                "dependence are all inside a fixed background coefficient, so this "
                "configuration has no water-vapour feedback either.",
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
                lower_fidelity_alternatives=("fixed_absorption_coefficient",),
            ),
        )
        super().__init__(contract, status=status)
        self.radiation_module_id = radiation_module_id
        self.concentration = concentration
        self._ppm = 0.0

    def validate_configuration(self) -> None:
        super().validate_configuration()
        if not callable(self.concentration) and float(self.concentration) <= 0.0:
            raise ConfigurationError("a carbon dioxide concentration must be positive")

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        self._ppm = self._at(float(boundary_state.get("master_start_time", 0.0)))
        self.mark_initialized()

    def _at(self, time: float) -> float:
        value = self.concentration(time) if callable(self.concentration) else self.concentration
        return float(value)

    def compute_tendencies(self) -> list[Any]:
        self._ppm = self._at(self.context.time)
        return []

    def _field(self, timestamp: float) -> dict[str, Field]:
        return {
            "co2_mole_fraction": Field(
                field_id="co2_mole_fraction",
                physical_name="carbon dioxide dry-air mole fraction",
                symbol="C",
                values=np.array(self._ppm),
                units="ppm",
                quantity_kind=QuantityKind.BOUNDARY_VALUE,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                timestamp=timestamp,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                provenance=Provenance(
                    producer=self.module_id,
                    derivation=EQUATION_IDS["co2_band_absorption"],
                    notes="prescribed boundary condition; no carbon cycle is modelled",
                ),
            )
        }

    def emit_state(self) -> dict[str, Field]:
        return self._field(self._context.time if self._context else 0.0)

    def emit_outputs(self) -> dict[str, Field]:
        # Republished every step. A PRESCRIBED module's emit_state runs only at
        # initialisation, so a time-varying concentration would otherwise stay
        # frozen at its first value -- the trap the cloud module fell into.
        return self._field(self.context.time)

    def report_uncertainty(self) -> UncertaintyReport:
        return UncertaintyReport(
            module_id=self.module_id,
            method="none; the concentration is imposed",
            distributional_assumptions="unresolved",
            correlations="unresolved",
            unresolved=(
                "the record's own measurement and reconstruction uncertainty is not "
                "propagated into the radiative effect",
                "structural_uncertainty_of_a_single_band_single_gas_representation",
            ),
        )

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "concentration_ppm": (
                "time_varying" if callable(self.concentration) else float(self.concentration)
            ),
            "band_centre_m": CO2_BAND_CENTRE_M,
            "centre_absorption_m2_kg": CO2_BAND_CENTRE_ABSORPTION,
            "wing_decay_per_metre": CO2_WING_DECAY_PER_METRE,
        }
