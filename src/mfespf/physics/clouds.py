"""Cloud state and its separate shortwave and longwave optical properties.

Charter 9.7 and section 19 stage 4; requirements CLD-001, CLD-002, RAD-002,
SPC-010, KER-012.

CLD-001 is the requirement that shapes this module: cloud state is a **separate
variable set**, never folded into a global albedo parameter. So a cloud here has
a liquid water path, an ice water path, an effective radius and a fraction, per
layer, and its radiative effect is computed from those rather than asserted.

The two effects are computed and reported separately before any net is formed
(SPC-010). They have opposite signs and comparable magnitudes, so a net cloud
radiative effect is a small difference of two large numbers; reporting only the
net would hide which of them a configuration actually got wrong.

Fidelity, on the charter's 9.7 path: this is the *prescribed* rung. The cloud
state comes from a declared function of time, not from the water cycle, so a
column using it does not conserve water through its clouds. That is stated in
the module's budget as a named source rather than left for a reader to notice.
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
from mfespf.core.versioning import InterfaceContract, Version
from mfespf.coupling.conservation import ModuleBudget
from mfespf.coupling.transfers import ConservedQuantity
from mfespf.physics.equations import EQUATION_IDS

__all__ = [
    "CLOUD_INTERFACE",
    "CloudOpticalProperties",
    "CloudProfile",
    "PrescribedCloudModule",
    "cloud_optical_properties",
    "single_layer_cloud",
]

V1 = Version(1, 0, 0)

#: Density of liquid water, used to turn a water path into an optical depth.
_WATER_DENSITY = 1000.0

#: Longwave mass absorption coefficient of cloud droplets, m^2 kg^-1. A standard
#: reduced value; it makes a 0.02 kg m^-2 cloud roughly half-emissive in the
#: infrared, which is where thin cirrus and thick stratus part company.
_LONGWAVE_MASS_ABSORPTION = 130.0

#: Asymmetry parameter for liquid droplets in the solar spectrum. Strongly
#: forward-scattering, which is why the solver delta-scales before using it.
_DROPLET_ASYMMETRY = 0.85

CLOUD_INTERFACE = InterfaceContract(
    interface_id="cloud",
    version=V1,
    physical_semantics={
        "liquid_water_path_units": "kg m^-2",
        "ice_water_path_units": "kg m^-2",
        "fraction_units": "1",
        "effective_radius_units": "m",
        "layer_ordering": "index_0_is_the_lowest_layer",
        "overlap": "maximum_with_a_single_column_fraction",
        "optical_property_basis": "per_layer_per_band_optical_depth",
    },
    optional_capabilities=frozenset({"ice_phase", "prescribed_time_variation"}),
)


class CloudProfile:
    """The cloud state a configuration hands to the radiation solver.

    Liquid and ice are carried separately because they scatter and absorb
    differently and because CLD-002 requires both to be exposed. The fraction is
    a single column value under maximum overlap: a per-layer fraction with a
    declared overlap assumption is the next rung, and pretending to it here would
    put a number in the manifest that the physics does not use.
    """

    def __init__(
        self,
        *,
        liquid_water_path: np.ndarray,
        ice_water_path: np.ndarray | None = None,
        fraction: float = 0.0,
        droplet_effective_radius: float = 1.0e-5,
        crystal_effective_radius: float = 3.0e-5,
    ) -> None:
        self.liquid_water_path = np.asarray(liquid_water_path, dtype=float)
        self.ice_water_path = (
            np.zeros_like(self.liquid_water_path)
            if ice_water_path is None
            else np.asarray(ice_water_path, dtype=float)
        )
        self.fraction = float(fraction)
        self.droplet_effective_radius = float(droplet_effective_radius)
        self.crystal_effective_radius = float(crystal_effective_radius)

        if self.ice_water_path.shape != self.liquid_water_path.shape:
            raise ConfigurationError("liquid and ice water paths must have the same shape")
        if np.any(self.liquid_water_path < 0.0) or np.any(self.ice_water_path < 0.0):
            raise ConfigurationError("cloud water paths cannot be negative")
        if not 0.0 <= self.fraction <= 1.0:
            raise ConfigurationError(f"cloud fraction {self.fraction} must lie in [0, 1]")
        if self.droplet_effective_radius <= 0.0 or self.crystal_effective_radius <= 0.0:
            raise ConfigurationError("effective radii must be positive")

    @property
    def n_layers(self) -> int:
        return int(self.liquid_water_path.size)

    @property
    def total_water_path(self) -> float:
        return float(np.sum(self.liquid_water_path) + np.sum(self.ice_water_path))

    @property
    def is_clear(self) -> bool:
        return self.fraction <= 0.0 or self.total_water_path <= 0.0


class CloudOpticalProperties:
    """Per-layer, per-band optical depths in the form the solver consumes."""

    def __init__(
        self, absorption: np.ndarray, scattering: np.ndarray, asymmetry: np.ndarray
    ) -> None:
        self.absorption = absorption
        self.scattering = scattering
        self.asymmetry = asymmetry

    @property
    def is_empty(self) -> bool:
        return not (np.any(self.absorption > 0.0) or np.any(self.scattering > 0.0))


def cloud_optical_properties(
    grid: SpectralGrid,
    profile: CloudProfile,
    *,
    shortwave_cutoff_m: float = 4.0e-6,
) -> CloudOpticalProperties:
    """Turn a cloud's water paths into optical depths per layer and band.

    Two regimes, split at the same cutoff the shortwave/longwave diagnostics use:

    * **Shortwave.** Geometric optics gives an extinction of
      ``3 W / (2 rho r_eff)``, essentially all of it scattering. This is where a
      cloud's cooling comes from, and it is why the solver had to learn to
      scatter before this module could exist.
    * **Longwave.** Droplets are comparable to the wavelength and absorb rather
      than scatter, so the optical depth is a mass absorption coefficient times
      the water path. This is where a cloud's warming comes from.

    Ice is given a larger effective radius and therefore less extinction per
    kilogram, which is the first-order difference between cirrus and stratus.
    Everything past that -- habit, size distribution, mixed phase -- is absent,
    and the module's limitations say so.
    """
    n_bands = grid.n_bins
    shape = (profile.n_layers, n_bands)
    absorption = np.zeros(shape)
    scattering = np.zeros(shape)
    asymmetry = np.zeros(shape)

    if profile.is_clear:
        return CloudOpticalProperties(absorption, scattering, asymmetry)

    shortwave = grid.centers < shortwave_cutoff_m
    water_path = profile.liquid_water_path + profile.ice_water_path

    liquid_extinction = (
        1.5 * profile.liquid_water_path / (_WATER_DENSITY * profile.droplet_effective_radius)
    )
    ice_extinction = (
        1.5 * profile.ice_water_path / (_WATER_DENSITY * profile.crystal_effective_radius)
    )
    shortwave_extinction = liquid_extinction + ice_extinction

    # Shortwave: scattering, with the small absorption of liquid water neglected.
    scattering[:, shortwave] = shortwave_extinction[:, None]
    asymmetry[:, shortwave] = _DROPLET_ASYMMETRY

    # Longwave: absorption, which by Kirchhoff is also the cloud's emission.
    absorption[:, ~shortwave] = (_LONGWAVE_MASS_ABSORPTION * water_path)[:, None]

    return CloudOpticalProperties(absorption, scattering, asymmetry)


def single_layer_cloud(
    *,
    n_layers: int,
    layer: int,
    liquid_water_path: float,
    fraction: float,
    ice: bool = False,
    droplet_effective_radius: float = 1.0e-5,
) -> CloudProfile:
    """A cloud confined to one layer -- the simplest prescribed state."""
    if not 0 <= layer < n_layers:
        raise ConfigurationError(f"layer {layer} is outside a column of {n_layers} layers")
    water = np.zeros(n_layers)
    water[layer] = liquid_water_path
    empty = np.zeros(n_layers)
    return CloudProfile(
        liquid_water_path=empty if ice else water,
        ice_water_path=water if ice else empty,
        fraction=fraction,
        droplet_effective_radius=droplet_effective_radius,
    )


class PrescribedCloudModule(Module):
    """Supplies a cloud state from a declared function of time.

    Status is PRESCRIBED (charter 6.4): the module supplies data and evolves no
    state of its own. The cloud's water is not drawn from the modeled water
    cycle, so it appears in the water budget as a named source and the run's
    conservation report shows the fiat rather than absorbing it.
    """

    def __init__(
        self,
        *,
        module_id: str = "cloud",
        grid: SpectralGrid,
        n_layers: int,
        profile: CloudProfile | Callable[[float], CloudProfile],
        shortwave_cutoff_m: float = 4.0e-6,
        status: ModuleStatus = ModuleStatus.PRESCRIBED,
    ) -> None:
        time_varying = callable(profile)
        contract = ModuleContract(
            module_id=module_id,
            module_name=(
                "prescribed time-varying cloud state"
                if time_varying
                else "prescribed static cloud state"
            ),
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.PRESCRIBED,
            physical_basis=(
                "Cloud water paths and fraction supplied from a declared function, with "
                "optical properties derived from them by geometric optics in the "
                "shortwave and a mass absorption coefficient in the longwave."
            ),
            governing_equation_ids=(
                EQUATION_IDS["cloud_optics"],
                EQUATION_IDS["cloud_radiative_effect"],
            ),
            assumptions=(
                "the cloud state is prescribed and responds to nothing in the model",
                "maximum overlap: all cloudy layers share one column fraction",
                "droplets and crystals have a single effective radius each",
                "shortwave scattering is conservative; droplet absorption is neglected",
                "no mixed-phase microphysics and no size distribution",
            ),
            validity_domain=(
                "Radiative experiments in which the cloud is an imposed boundary "
                "condition. Because the cloud does not respond to the column, no claim "
                "about cloud feedback is in scope; only the effect of a given cloud is."
            ),
            required_inputs=("atmosphere_layer_pressure",),
            produced_outputs=(
                "cloud_liquid_water_path",
                "cloud_ice_water_path",
                "cloud_fraction",
                "cloud_effective_radius",
                "cloud_shortwave_optical_depth",
                "cloud_longwave_optical_depth",
                "cloud_total_water_path",
            ),
            conserved_quantities=("water",),
            internal_sources_and_sinks=("prescribed_cloud_water_source",),
            provided_interfaces=(CLOUD_INTERFACE,),
            uncertainty_outputs=(
                "prescribed_state_structural",
                "effective_radius_parameter",
                "overlap_assumption_structural",
                "neglected_shortwave_absorption",
            ),
            correctness_tests=(
                "test_no_cloud_control_matches_the_clear_column",
                "test_cloud_shortwave_effect_cools_and_longwave_effect_warms",
                "test_cloud_optical_depth_scales_with_water_path",
                "test_a_thin_cloud_is_longwave_dominated_and_a_thick_one_shortwave",
                "test_prescribed_cloud_water_is_declared_as_a_source",
            ),
            known_limitations=(
                "The cloud is imposed, so the configuration cannot exhibit cloud "
                "feedback of any kind.",
                "Its water is not taken from the water cycle, so a run using it does "
                "not conserve water through the cloud; the budget names the source.",
                "One effective radius per phase, no size distribution, no mixed phase.",
                "Maximum overlap only; partial overlap between layers is not represented.",
            ),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(n_layers * n_bands)",
                memory_complexity="O(n_layers * n_bands)",
                scaling_variables=("n_bands", "n_layers"),
                work_units_per_step=0.5,
                benchmark_configuration=f"{n_layers} layers, {grid.n_bins} bands",
                peak_memory_bytes=8 * n_layers * grid.n_bins + 4096,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
                lower_fidelity_alternatives=("no_cloud_control",),
            ),
        )
        super().__init__(contract, status=status)
        self.grid = grid
        self.n_layers = n_layers
        self._profile_source = profile
        self.shortwave_cutoff_m = shortwave_cutoff_m
        self.profile: CloudProfile = profile(0.0) if callable(profile) else profile
        self.optics = cloud_optical_properties(
            grid, self.profile, shortwave_cutoff_m=shortwave_cutoff_m
        )
        self._water_initial = self.profile.total_water_path

    def validate_configuration(self) -> None:
        super().validate_configuration()
        if self.profile.n_layers != self.n_layers:
            raise ConfigurationError(
                f"cloud profile has {self.profile.n_layers} layers, expected {self.n_layers}"
            )

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        self._water_initial = self.profile.total_water_path
        self.mark_initialized()

    def compute_tendencies(self) -> list[Any]:
        ctx = self.context
        self._water_initial = self.profile.total_water_path
        if callable(self._profile_source):
            self.profile = self._profile_source(ctx.time)
            if self.profile.n_layers != self.n_layers:
                raise ConfigurationError(
                    f"prescribed cloud returned {self.profile.n_layers} layers at "
                    f"t={ctx.time}, expected {self.n_layers}"
                )
        self.optics = cloud_optical_properties(
            self.grid, self.profile, shortwave_cutoff_m=self.shortwave_cutoff_m
        )
        return []

    def emit_state(self) -> dict[str, Field]:
        """Publish before the first step, so radiation has a cloud to read."""
        return self._cloud_fields(self._context.time if self._context else 0.0)

    def emit_outputs(self) -> dict[str, Field]:
        """Republish within the step, in splitting order.

        A prescribed cloud may vary with time, and it sits ahead of radiation in
        the splitting order precisely so that radiation sees *this* step's cloud.
        Publishing only through ``emit_state`` would leave the state one step
        stale, because the coupler republishes that after the advance -- which a
        module evolving no state never reaches.
        """
        return self._cloud_fields(self.context.time)

    def _cloud_fields(self, timestamp: float) -> dict[str, Field]:
        provenance = Provenance(
            producer=self.module_id,
            derivation=EQUATION_IDS["cloud_optics"],
            notes="prescribed; not derived from the modeled water cycle",
        )
        profile = {
            # Charter 7.2 has no "prescribed" kind; a supplied cloud state is a
            # boundary value, and the module's PRESCRIBED status carries the rest.
            "quantity_kind": QuantityKind.BOUNDARY_VALUE,
            "time_semantics": TimeSemantics.INSTANTANEOUS,
            "timestamp": timestamp,
            "spatial_support": SpatialSupport.LAYER_AVERAGE,
            "dimensions": ("layer",),
            "coordinate_system": "pressure_layers",
            "sign_convention": SignConvention.POSITIVE_DEFINITE,
            "provenance": provenance,
        }
        shortwave = self.grid.centers < self.shortwave_cutoff_m
        return {
            # CLD-002: liquid, ice, fraction and vertical support, all exposed.
            "cloud_liquid_water_path": Field(
                field_id="cloud_liquid_water_path",
                physical_name="cloud liquid water path in each layer",
                symbol="LWP",
                values=np.array(self.profile.liquid_water_path),
                units="kg m^-2",
                phase_or_species="cloud_liquid",
                **profile,  # type: ignore[arg-type]
            ),
            "cloud_ice_water_path": Field(
                field_id="cloud_ice_water_path",
                physical_name="cloud ice water path in each layer",
                symbol="IWP",
                values=np.array(self.profile.ice_water_path),
                units="kg m^-2",
                phase_or_species="cloud_ice",
                **profile,  # type: ignore[arg-type]
            ),
            "cloud_shortwave_optical_depth": Field(
                field_id="cloud_shortwave_optical_depth",
                physical_name="cloud shortwave optical depth in each layer",
                values=np.array(
                    self.optics.scattering[:, shortwave].max(axis=1)
                    if np.any(shortwave)
                    else np.zeros(self.n_layers)
                ),
                units="1",
                **profile,  # type: ignore[arg-type]
            ),
            "cloud_longwave_optical_depth": Field(
                field_id="cloud_longwave_optical_depth",
                physical_name="cloud longwave optical depth in each layer",
                values=np.array(
                    self.optics.absorption[:, ~shortwave].max(axis=1)
                    if np.any(~shortwave)
                    else np.zeros(self.n_layers)
                ),
                units="1",
                **profile,  # type: ignore[arg-type]
            ),
            "cloud_fraction": Field(
                field_id="cloud_fraction",
                physical_name="column cloud fraction under maximum overlap",
                symbol="f_c",
                values=np.array(self.profile.fraction),
                units="1",
                quantity_kind=QuantityKind.BOUNDARY_VALUE,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                timestamp=timestamp,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                provenance=provenance,
            ),
            "cloud_total_water_path": Field(
                field_id="cloud_total_water_path",
                physical_name="total cloud condensed water path",
                values=np.array(self.profile.total_water_path),
                units="kg m^-2",
                quantity_kind=QuantityKind.EXTENSIVE_STATE,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                timestamp=timestamp,
                spatial_support=SpatialSupport.DOMAIN_INTEGRAL,
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                provenance=provenance,
            ),
            "cloud_effective_radius": Field(
                field_id="cloud_effective_radius",
                physical_name="droplet and crystal effective radii",
                values=np.array(
                    [
                        self.profile.droplet_effective_radius,
                        self.profile.crystal_effective_radius,
                    ]
                ),
                units="m",
                quantity_kind=QuantityKind.BOUNDARY_VALUE,
                time_semantics=TimeSemantics.INSTANTANEOUS,
                timestamp=timestamp,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                dimensions=("phase",),
                sign_convention=SignConvention.POSITIVE_DEFINITE,
                provenance=provenance,
            ),
        }

    def report_conservation(self) -> list[ModuleBudget]:
        """The cloud's water, and where it came from.

        A prescribed cloud's condensate is not drawn from the modeled water
        cycle. Rather than leave that outside the accounting, the change is
        booked against a named source, so the run's conservation report shows a
        line called ``prescribed_cloud_water_source`` and a reader can see the
        water was supplied by fiat rather than by physics.
        """
        budget = ModuleBudget(
            module_id=self.module_id,
            quantity=ConservedQuantity.WATER,
            units="kg m^-2",
            stored_initial=self._water_initial,
            stored_final=self.profile.total_water_path,
        )
        budget.add_source(
            "prescribed_cloud_water_source",
            self.profile.total_water_path - self._water_initial,
        )
        return [budget]

    def report_uncertainty(self) -> UncertaintyReport:
        return UncertaintyReport(
            module_id=self.module_id,
            quantified={
                "cloud_fraction": self.profile.fraction,
                "total_water_path_kg_m2": self.profile.total_water_path,
            },
            method="none propagated; the cloud state is imposed rather than predicted",
            distributional_assumptions="unresolved",
            correlations="unresolved",
            unresolved=(
                "structural_uncertainty_of_a_prescribed_cloud_state",
                "effective_radius_parameter_uncertainty",
                "structural_uncertainty_of_maximum_overlap",
                "neglected_shortwave_droplet_absorption",
            ),
        )

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "n_layers": self.n_layers,
            "time_varying": callable(self._profile_source),
            "cloud_fraction": self.profile.fraction,
            "liquid_water_path_kg_m2": self.profile.liquid_water_path.tolist(),
            "ice_water_path_kg_m2": self.profile.ice_water_path.tolist(),
            "droplet_effective_radius_m": self.profile.droplet_effective_radius,
            "crystal_effective_radius_m": self.profile.crystal_effective_radius,
            "longwave_mass_absorption_m2_kg": _LONGWAVE_MASS_ABSORPTION,
            "droplet_asymmetry": _DROPLET_ASYMMETRY,
            "overlap": "maximum",
        }
