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
    "CO2_WING_DECAY_WITH_VAPOUR_PER_METRE",
    "DRY_AIR_MOLAR_MASS",
    "INFRARED_BACKGROUND_ABSORPTION",
    "PREINDUSTRIAL_CO2_PPM",
    "VAPOUR_BAND_ABSORPTION",
    "VAPOUR_WINDOW_CONTINUUM",
    "VAPOUR_WING_DECAY_PER_METRE",
    "PrescribedGreenhouseModule",
    "co2_band_absorption",
    "co2_mass_mixing_ratio",
    "co2_wing_decay_for",
    "greenhouse_grid",
    "vapour_absorption_optical_depth",
    "vapour_mass_absorption",
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

INFRARED_BACKGROUND_ABSORPTION = 3.0e-5
"""Infrared absorption of everything that is not CO2, m^2 per kg of air.

Set from the transparency of the atmospheric window rather than from any
temperature: over a column of about ten tonnes per square metre this gives an
optical depth near 0.3 in the window, which is the order the water-vapour
continuum produces.

It lives here, and is the default of ``co2_band_absorption``, because the wing
decay below is calibrated *against it*. They were once separated -- the
calibration used one value and every configuration passed another -- and the
forcing per doubling in the runs came out 1.75 times the calibrated figure. One
constant, used in both places, is what stops that recurring.
"""

CO2_BAND_LOWER_M = 8.0e-6
CO2_BAND_UPPER_M = 3.0e-5
CO2_BAND_RESOLUTION = 64
"""The band region and the resolution across it, at which the wing decay is fitted.

Three constants rather than three defaults scattered over the call sites, for
the same reason the background above is one constant: the calibration is only
worth anything on the grid it was performed on, and a configuration that resolves
the band differently gets a different forcing from the same physics.

The region is where the band actually is. At the highest concentration this
framework runs, the opaque core spans about 12 to 20 microns and the wings have
fallen by nineteen orders of magnitude by 8 microns, so nothing is lost by
resolving 8 to 30 and leaving the rest to the background spacing -- and the bands
saved there are spent where the wings are. ``CO2_BAND_RESOLUTION`` is the point
where the forcing stops moving: 64 gives 3.698 W m^-2 per doubling and a
reference grid of 179 bands gives 3.698.

The previous choice -- 48 bands over 6.5 to 50 microns -- was *not* converged
once the band narrowed, and gave 3.700 where the converged answer was 3.601. A
calibration performed there would have fitted a property of the grid.
"""

CO2_WING_DECAY_PER_METRE = 1.41e3
"""Wing decay scale in inverse metres, equivalent to about 14 cm^-1.

The one calibrated number in this module, and it is calibrated against a
*published radiative-forcing value* -- 3.7 W m^-2 for a doubling -- never against
a temperature record. That is what keeps a temperature record available as
independent validation afterwards, and the dataset registry enforces the
distinction (INF-003, EQR-004).

It is fitted at the background above and on the grid above, because it is only
meaningful with both. Nothing else is tuned: the base state and the greenhouse
effect that follow are *predictions*, and both come out about five kelvin low
against an observed 288 K and 33 K. That shortfall is where the water-vapour
greenhouse effect would be, if the vapour this framework already tracks were not
radiatively inert.
"""

# --- water vapour ------------------------------------------------------------

VAPOUR_ROTATIONAL_EDGE_PER_METRE = 2.0e4
"""Short-wavelength edge of the pure rotational band, 200 cm^-1 (50 microns).

Below this wavenumber the band is treated as flat and saturated, which is what
the real rotational band does: the far infrared is opaque in any atmosphere
carrying ordinary amounts of vapour. Above it the envelope decays.
"""

VAPOUR_VIBRATIONAL_CENTRE_PER_METRE = 1.595e5
"""Centre of the 6.3 micron bending band, 1595 cm^-1."""

VAPOUR_BAND_ABSORPTION = 100.0
"""Band-centre mass absorption coefficient, m^2 per kg of *vapour*.

Not fitted, for the same reason ``CO2_BAND_CENTRE_ABSORPTION`` is not: over a
column carrying about 23 kg m^-2 it gives an optical depth of order two
thousand, so the cores are deeply saturated and their exact value stops
mattering. What the fitted numbers below control is where the wings cross unit
optical depth.
"""

VAPOUR_WING_DECAY_PER_METRE = 3.970e3
"""Wing decay of both vapour bands, about 36 cm^-1. Calibrated.

Fitted against a *published radiation observation* -- the CERES all-sky
outgoing longwave flux of about 240 W m^-2 -- and never against a temperature
record, which is what keeps HadCRUT5 available as validation (INF-003).

**All-sky, not clear-sky, and the difference is not a detail.** It was fitted to
the clear-sky 265 first, and that was wrong: these columns carry an albedo of
0.3, which is the *all-sky* planetary albedo, so they absorb 238 W m^-2 and in
equilibrium must emit 238. Pairing a clear-sky radiative target with an all-sky
albedo gives a column the cloud's cooling without the cloud's warming, and it
settled 7.5 K below its own calibration point -- which then halved the feedback,
because Clausius-Clapeyron is exponential and a cold column is a dry one.

So the vapour opacity here stands in for the cloud longwave effect as well as
its own. That is a stand-in with a name rather than a hidden one.

**Refitted when the saturation formula gained its ice phase.** It was 3617 m^-1
while saturation was computed over liquid water all the way down; correcting the
upper column to sublimation removed vapour from the cold layers -- a factor of
five at 172 K -- and the same target then needs about ten percent more
absorption per kilogram. The number changed because the physics under it did,
which is what a calibration tied to its conditions is supposed to do.

It weakens the base state as a prediction without abolishing it. Equilibrium
requires the outgoing flux to equal the absorbed flux, and the fit sets the
outgoing flux at 288 K -- but it sets it for the *prescribed* profile the
calibration uses, and a coupled column builds its own profile by convective
adjustment. So the column does not arrive at 288 K by construction: it arrives
at 293.0, five kelvin warm, with a greenhouse effect of 38 K against an observed
33. What is constructed is the flux at one profile; what remains predicted is
everything the column does with a different one.

One decay serves both bands. Line strengths do not really fall away from the
rotational and the vibrational band at the same rate, and giving each its own
scale would have added a second fitted number for a refinement nothing here can
test.
"""

VAPOUR_WINDOW_CONTINUUM = 1.86e-2
"""Continuum absorption across the atmospheric window, m^2 per kg of vapour.

Calibrated against a published clear-sky window transmittance near 0.65 at 10
microns. It is what stops the window from being perfectly transparent, and it is
the reason the window closes as the atmosphere moistens.

Represented as flat in wavelength, which the real continuum is not -- it
strengthens towards longer wavelengths and is quadratic rather than linear in
vapour amount, because it is dominated by self-broadening. Both are declared
approximations of the reduced model and both are listed in the equation record.
"""

CO2_WING_DECAY_WITH_VAPOUR_PER_METRE = 1.2725e3
"""The CO2 wing decay for a column whose infrared background *is* the vapour.

A band's fitted width depends on what else is absorbing around it. With water
vapour active the 15 micron band is partly masked, so reproducing the same
published 3.7 W m^-2 per doubling takes a narrower band than the 1410 m^-1 that
``CO2_WING_DECAY_PER_METRE`` uses against the fixed background.

Two calibrations of one band therefore coexist, and that is not a duplication to
be tidied away: each is correct for its own configuration and neither is correct
for the other. ``co2_wing_decay_for`` selects between them so that no call site
has to remember which is which.
"""


def co2_wing_decay_for(*, active_vapour: bool) -> float:
    """The CO2 wing decay calibrated for the opacity this configuration runs at.

    The whole point of the constant this returns is that it is conditional. A
    configuration that picked the wrong one would use a forcing per doubling its
    calibration never produced -- the defect that invalidated four claims, and
    the reason this is a function rather than a default argument someone has to
    remember to override.
    """
    return CO2_WING_DECAY_WITH_VAPOUR_PER_METRE if active_vapour else CO2_WING_DECAY_PER_METRE


def vapour_mass_absorption(
    grid: SpectralGrid,
    *,
    wing_decay_per_metre: float = VAPOUR_WING_DECAY_PER_METRE,
    window_continuum: float = VAPOUR_WINDOW_CONTINUUM,
    band_absorption: float = VAPOUR_BAND_ABSORPTION,
    rotational_edge_per_metre: float = VAPOUR_ROTATIONAL_EDGE_PER_METRE,
    vibrational_centre_per_metre: float = VAPOUR_VIBRATIONAL_CENTRE_PER_METRE,
    shortwave_cutoff_m: float = 4.0e-6,
) -> np.ndarray:
    """Absorption coefficient per kilogram of water vapour, one value per band.

    Three features, and each is a spectroscopic fact rather than a fitted shape:
    a pure rotational band that makes the far infrared opaque, a vibrational
    band at 6.3 microns, and a continuum across the window between them. The
    window is what lets the surface radiate to space at all, and the continuum
    is what closes it as the atmosphere moistens.

    Confined to the infrared exactly. Water vapour does absorb in the near
    infrared, and leaving that out means this configuration cannot represent the
    shortwave absorption that warms the atmosphere directly -- declared in the
    equation record, not hidden in a tolerance.
    """
    if band_absorption < 0.0 or window_continuum < 0.0:
        raise ConfigurationError("absorption coefficients must be non-negative")
    if wing_decay_per_metre <= 0.0:
        raise ConfigurationError("the wing decay scale must be positive")

    wavenumber = 1.0 / grid.centers
    # Flat and saturated inside the rotational band, decaying above its edge.
    rotational = np.exp(
        -np.maximum(wavenumber - rotational_edge_per_metre, 0.0) / wing_decay_per_metre
    )
    vibrational = np.exp(-np.abs(wavenumber - vibrational_centre_per_metre) / wing_decay_per_metre)
    kappa = window_continuum + band_absorption * (rotational + vibrational)
    return np.where(grid.centers >= shortwave_cutoff_m, kappa, 0.0)


def vapour_absorption_optical_depth(
    grid: SpectralGrid,
    vapour_path: np.ndarray,
    **kwargs: Any,
) -> np.ndarray:
    """Absorption optical depth of the vapour in each layer, shaped (layers, bands).

    Water vapour is not uniformly mixed -- most of it is in the lowest couple of
    kilometres -- so it cannot enter as a per-band coefficient multiplied by
    layer mass the way CO2 does. It enters as an optical depth per layer, which
    is the same shape a cloud or an aerosol uses, and ``combine_scatterers``
    adds it to theirs.

    ``vapour_path`` is the vapour mass per unit area in each layer, kg m^-2 --
    exactly what the water module owns and publishes.
    """
    path = np.asarray(vapour_path, dtype=float)
    if np.any(path < 0.0):
        raise ConfigurationError("a vapour path cannot be negative")
    return path[:, None] * vapour_mass_absorption(grid, **kwargs)[None, :]


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
    n_band: int = CO2_BAND_RESOLUTION,
    *,
    lower_m: float = 1.0e-7,
    upper_m: float = 1.0e-3,
    band_lower_m: float = CO2_BAND_LOWER_M,
    band_upper_m: float = CO2_BAND_UPPER_M,
) -> SpectralGrid:
    """A log-spaced grid with extra resolution across the CO2 band.

    The default column grid puts two or three bands across 12 to 20 microns,
    which cannot resolve a band whose wings are the whole point: with the wings
    unresolved the band either saturates completely or scales linearly, and
    neither is what carbon dioxide does.

    The extra bands are inside the CO2 region and nowhere else, so the cost is
    paid where it buys something. Everything outside is the same logarithmic
    spacing the other configurations use.

    The defaults are the grid the band was calibrated on, and they are the
    constants above rather than literals here, so that a caller who takes the
    defaults gets the calibrated forcing and a caller who does not is visibly
    choosing something else.
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
    background: float = INFRARED_BACKGROUND_ABSORPTION,
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
    difference between two runs is attributable to CO2 alone -- which also means
    there is no water-vapour feedback, and the default is the value the wing
    decay was calibrated against.
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
