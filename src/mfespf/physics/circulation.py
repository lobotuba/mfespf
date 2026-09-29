"""An overturning circulation between a rising and a sinking column.

Why this exists
---------------

A single horizontally uniform column has no net vertical motion, so it cannot
hold the sinking air that dries Earth's free troposphere. The observed-cloud
test failed on that: the model's free troposphere sat at 0.775 relative
humidity against an observed 0.4 to 0.5, and its air at ISCCP's cloud-top
pressures was 13 to 28 K warmer than ISCCP's cloud tops
(``docs/protocols/observed-cloud-field.md``). This module links two columns so
air can rise in one and sink in the other.

The physics
-----------

The weak temperature gradient is what gravity waves leave behind: they carry
away any horizontal temperature difference in the time they take to cross the
distance between rising and sinking air (Raymond and Zeng 2005). Each column
relaxes toward the area mean on that time, ``tau_1 = L / c_1``. Each vertical
mode relaxes on its own time: mode ``n`` travels at ``c_1 / n`` and so relaxes
``n`` times more slowly (Herman and Raymond 2014). The grid-scale structure
convection leaves in a column is then weakly coupled, as the waves make it.

The modes are sines in pressure between the surface and the tropopause. They
are exact only for a uniform buoyancy frequency, so the stratification used is
the rising column's bulk value over the troposphere, never a local
layer-to-layer step. The strict balance this replaced divided the cooling by
the local step, and near the cold top of convection that step goes to zero:
1100 hPa/day at 110 hPa, a layer emptied in one step
(``docs/protocols/two-column-circulation.md``).

Dry static energy ``s = c_p T + g z`` and water vapour move with the resulting
fluxes in flux form, upwind: vertical fluxes carry the upstream layer's value
and the lateral exchange the upstream column's. So whatever leaves one layer or
column arrives in another, exactly, and between the columns it is recorded as a
transfer per unit of global area, each budget taking its own area share.

Declared and measured in ``docs/protocols/spectral-wtg-circulation.md``.
"""

from __future__ import annotations

from collections.abc import Sequence
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
from mfespf.core.constants import value_of
from mfespf.core.errors import ConfigurationError
from mfespf.core.quantities import (
    Field,
    Provenance,
    QuantityKind,
    SpatialSupport,
    TimeSemantics,
)
from mfespf.core.versioning import InterfaceContract, Version
from mfespf.coupling.conservation import ModuleBudget
from mfespf.coupling.transfers import ConservedQuantity, Transfer
from mfespf.physics.equations import EQUATION_IDS
from mfespf.physics.thermodynamics import hydrostatic_heights
from mfespf.registries.states import TendencyContribution

__all__ = [
    "MAX_MODES",
    "OBSERVED_ASCENT_FRACTION",
    "OBSERVED_RELAXATION_DISTANCE_M",
    "TROPOPAUSE_PA",
    "AreaWeightedMeansModule",
    "WeakTemperatureGradientCirculation",
    "spectral_wtg_subsidence",
    "upwind_tendency",
]

V1 = Version(1, 0, 0)

#: Area-weighted fraction of the globe with upward monthly-mean vertical
#: velocity at 500 hPa, NCEP/NCAR Reanalysis 1, 1991-2020 long-term means
#: (NOAA PSL ``omega.mon.ltm.1991-2020.nc``); 0.461 to 0.490 by month.
OBSERVED_ASCENT_FRACTION = 0.476

#: Tropopause the modes vanish at, Pa: the nominal upper boundary of the
#: GASS WTG intercomparison and SCAM6 (Daleu et al. 2016; Cohen et al. 2024).
TROPOPAUSE_PA = 10000.0

#: Distance gravity waves cross between sinking and rising air, m: the
#: area-weighted mean great-circle distance from each descending grid point
#: (omega500 > 0) within 20-30 degrees of the equator to the nearest ascending
#: one, NCEP/NCAR Reanalysis 1 1991-2020 monthly long-term means; 734 km
#: (20 degrees), 770 km (30 degrees), monthly 589 to 853.
OBSERVED_RELAXATION_DISTANCE_M = 7.5e5

#: Vertical modes kept at most (Cohen et al. 2024 found the truncation error
#: negligible at 20); fewer when the grid resolves fewer.
MAX_MODES = 20

_TEMPERATURE = "atmosphere_layer_temperature"
_VAPOUR = "atmosphere_water_vapour"

CIRCULATION_INTERFACE = InterfaceContract(
    interface_id="overturning_circulation",
    version=V1,
    physical_semantics={
        "mass_flux_sign": "positive_downward",
        "mass_flux_basis": "per_unit_area_of_the_sinking_column",
        "balance": "relaxed_spectral_weak_temperature_gradient",
        "transport": "flux_form_upwind_dry_static_energy_and_vapour",
        "transfer_amount_basis": "per_unit_global_area",
    },
)


def spectral_wtg_subsidence(
    *,
    temperature_anomaly: np.ndarray,
    specific_heat: np.ndarray,
    bulk_stratification: float,
    relaxation_time: float,
    interface_pressure: np.ndarray,
    tropopause: float = TROPOPAUSE_PA,
    gravity: float = 9.80665,
    max_modes: int = MAX_MODES,
) -> np.ndarray:
    """Downward mass flux at each interface of a column, kg m^-2 s^-1.

    ``temperature_anomaly`` is the column's departure from the area mean, per
    layer. The flux that would remove it on ``relaxation_time`` is
    ``f = c_p T' / (g tau_1 dS/dp)``, with ``bulk_stratification`` the uniform
    ``dS/dp`` (J kg^-1 Pa^-1, negative when stable). ``f`` is projected onto
    ``sin(n pi zeta)``, ``zeta = (p_s - p) / (p_s - p_t)``, weighting each
    layer by its pressure thickness below the tropopause, and mode ``n`` is
    divided by ``n``: it travels ``n`` times more slowly. The flux vanishes at
    the surface, at the tropopause and above it.
    """
    anomaly = np.asarray(temperature_anomaly, dtype=float)
    c_p = np.asarray(specific_heat, dtype=float)
    p = np.asarray(interface_pressure, dtype=float)
    n = anomaly.size
    if c_p.size != n or p.size != n + 1:
        raise ConfigurationError(
            "anomaly, heat capacity and interfaces must describe the same layers"
        )
    if not bulk_stratification < 0.0:
        raise ConfigurationError(
            f"the troposphere's bulk stratification {bulk_stratification:g} J kg^-1 Pa^-1 "
            "is not stable, so no gravity wave relaxes it"
        )
    if not relaxation_time > 0.0:
        raise ConfigurationError("the relaxation time must be positive")
    surface = float(p[0])
    depth = surface - tropopause
    if depth <= 0.0:
        raise ConfigurationError("the tropopause must lie above the surface")
    lower = np.minimum(p[:-1], surface)
    upper = np.maximum(p[1:], tropopause)
    weight = np.clip(lower - upper, 0.0, None) / depth
    below = weight > 0.0
    centre = (surface - 0.5 * (lower + upper)) / depth
    forcing = c_p * anomaly / (gravity * relaxation_time * bulk_stratification)
    modes = min(int(np.count_nonzero(below)), int(max_modes))
    zeta = np.clip((surface - p) / depth, 0.0, 1.0)
    flux = np.zeros(n + 1)
    for mode in range(1, modes + 1):
        amplitude = 2.0 * float(
            np.sum(forcing[below] * np.sin(mode * np.pi * centre[below]) * weight[below])
        )
        flux += amplitude / mode * np.sin(mode * np.pi * zeta)
    # Exactly zero, not round-off, where the modes vanish.
    flux[0] = 0.0
    flux[p <= tropopause] = 0.0
    return flux


def upwind_tendency(
    flux: np.ndarray, own: np.ndarray, other: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Rate of change of ``phi * m`` per layer, and each layer's lateral term.

    ``flux`` is the downward mass flux at each interface of this column (per
    unit of its own area); ``own`` and ``other`` are the transported quantity per
    kilogram in this column and the other one. Vertical fluxes take the upstream
    layer's value, lateral exchange the upstream column's. Mass is conserved
    layer by layer, so the lateral mass flux is ``flux[k] - flux[k + 1]``.
    """
    f = np.asarray(flux, dtype=float)
    phi = np.asarray(own, dtype=float)
    far = np.asarray(other, dtype=float)
    n = phi.size
    rate = np.zeros(n)
    lateral = np.zeros(n)
    for k in range(n):
        top, bottom = f[k + 1], f[k]
        if top > 0.0:
            rate[k] += top * phi[k + 1]
        elif top < 0.0:
            rate[k] += top * phi[k]
        if bottom > 0.0:
            rate[k] -= bottom * phi[k]
        elif bottom < 0.0:
            rate[k] -= bottom * phi[k - 1]
        inflow = bottom - top
        lateral[k] = inflow * (far[k] if inflow > 0.0 else phi[k])
        rate[k] += lateral[k]
    return rate, lateral


class WeakTemperatureGradientCirculation(Module):
    """Couples a rising and a sinking column (EQ-WTG-OVERTURNING-CIRCULATION)."""

    def __init__(
        self,
        *,
        module_id: str = "circulation",
        ascent_region: str = "ascent",
        descent_region: str = "descent",
        ascent_fraction: float = OBSERVED_ASCENT_FRACTION,
        relaxation_distance: float = OBSERVED_RELAXATION_DISTANCE_M,
        tropopause: float = TROPOPAUSE_PA,
        status: ModuleStatus = ModuleStatus.ACTIVE,
    ) -> None:
        regions = (ascent_region, descent_region)
        per_region = (
            _TEMPERATURE,
            _VAPOUR,
            "atmosphere_layer_mass",
            "atmosphere_layer_pressure",
            "atmosphere_interface_pressure",
            "atmosphere_layer_heat_capacity",
            "layer_radiative_heating",
        )
        contract = ModuleContract(
            module_id=module_id,
            module_name="weak-temperature-gradient overturning circulation",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.REDUCED,
            physical_basis=(
                "Gravity waves remove horizontal temperature differences in the time "
                "they take to cross between rising and sinking air, each vertical mode "
                "at its own speed; the vertical motion that does it sinks where the air "
                "is colder than the mean and rises where it is warmer. Mass continuity "
                "closes the overturning, and energy and water are carried with it."
            ),
            governing_equation_ids=(EQUATION_IDS["wtg_circulation"],),
            assumptions=(
                "one rising and one sinking region stand for the globe",
                "each column relaxes toward the area mean on the gravity-wave crossing time",
                "the modes assume a uniform stratification, the rising column's bulk value",
                "transport is flux-form and upwind in both directions",
            ),
            validity_domain=(
                "Global-mean experiments in which the WTG balance is a fair description "
                "of the overturning; midlatitude eddies are not represented."
            ),
            required_inputs=tuple(f"{r}:{f}" for r in regions for f in per_region),
            produced_outputs=(
                "circulation_subsidence_mass_flux",
                "circulation_subsidence_omega_500",
                "circulation_relaxation_time",
                # Published into each region, so the columns' water and
                # convection modules compose with this step's transport.
                *(
                    f"{r}:{f}"
                    for r in regions
                    for f in ("vapour_advection_rate", "advective_heating_rate")
                ),
            ),
            conserved_quantities=("energy", "water"),
            internal_sources_and_sinks=(),
            provided_interfaces=(CIRCULATION_INTERFACE,),
            uncertainty_outputs=(
                "two_region_structural",
                "relaxation_distance_observational",
                "tropopause_parameter",
                "ascent_fraction_observational",
            ),
            correctness_tests=(
                "test_the_circulation_conserves_energy_and_water_exactly",
                "test_the_gravest_mode_relaxes_on_tau_1",
                "test_each_mode_relaxes_n_times_more_slowly",
                "test_no_temperature_difference_means_no_circulation",
                "test_the_two_columns_carry_equal_and_opposite_mass",
            ),
            known_limitations=(
                "Two regions only; no midlatitude dynamics.",
                "The modes are sines in pressure, exact only for uniform stratification.",
                "The rising column's ascent is uniform across it.",
            ),
            fidelity_justification=(
                "The single column is the null alternative: it cannot represent "
                "subsidence at all, which is the defect this module exists to remove."
            ),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(n_layers)",
                memory_complexity="O(n_layers)",
                scaling_variables=("n_layers",),
                work_units_per_step=0.2,
                benchmark_configuration="two 30-layer columns",
                peak_memory_bytes=65536,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
                lower_fidelity_alternatives=("single_column",),
            ),
        )
        super().__init__(contract, status=status)
        self.ascent_region = ascent_region
        self.descent_region = descent_region
        self.ascent_fraction = float(ascent_fraction)
        self.relaxation_distance = float(relaxation_distance)
        self.tropopause = float(tropopause)
        self.gas_constant = value_of("gas_constant_dry_air")
        self._relaxation_time = 0.0
        self.latent_heat = value_of("water_latent_heat_vaporization_273")
        self.gravity = value_of("earth_standard_gravity")
        self._flux = np.zeros(0)
        self._omega_500 = 0.0
        self._transfers: list[Transfer] = []
        self._rates: dict[str, np.ndarray] = {}

    def validate_configuration(self) -> None:
        super().validate_configuration()
        if not 0.0 < self.ascent_fraction < 1.0:
            raise ConfigurationError("the rising column's area fraction must lie in (0, 1)")
        if not self.relaxation_distance > 0.0:
            raise ConfigurationError("the relaxation distance must be positive")

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        self.mark_initialized()

    def _column(self, region: str) -> dict[str, np.ndarray]:
        def get(name: str, units: str) -> np.ndarray:
            return np.asarray(
                self.input_field(f"{region}:{name}").require_units(units).values, dtype=float
            )

        column = {
            "T": get(_TEMPERATURE, "K"),
            "w": get(_VAPOUR, "kg m^-2"),
            "m": get("atmosphere_layer_mass", "kg m^-2"),
            "p": get("atmosphere_layer_pressure", "Pa"),
            "pi": get("atmosphere_interface_pressure", "Pa"),
            "C": get("atmosphere_layer_heat_capacity", "J m^-2 K^-1"),
            "Q": get("layer_radiative_heating", "W m^-2"),
        }
        z = hydrostatic_heights(column["T"], column["p"], float(column["pi"][0]))
        column["s"] = column["C"] / column["m"] * column["T"] + self.gravity * z
        column["q"] = column["w"] / column["m"]
        return column

    def _troposphere(self, column: dict[str, np.ndarray]) -> tuple[float, float]:
        """The rising column's bulk stratification and gravest-mode relaxation time.

        Bulk, between its lowest layer and its highest layer below the
        tropopause, because the modes are exact only for a uniform buoyancy
        frequency. ``c_1 = N H / pi`` with ``H`` the tropopause height, and
        ``tau_1 = L / c_1``.
        """
        p, T, s = column["p"], column["T"], column["s"]
        inside = np.flatnonzero(p > self.tropopause)
        if inside.size < 2:
            raise ConfigurationError(
                "the troposphere needs at least two layers below the tropopause"
            )
        bottom, top = int(inside[0]), int(inside[-1])
        stratification = float((s[top] - s[bottom]) / (p[top] - p[bottom]))
        z = hydrostatic_heights(T, p, float(column["pi"][0]))
        mass = column["m"][inside]
        mean_temperature = float(np.sum(T[inside] * mass) / np.sum(mass))
        specific_heat = float(np.sum(column["C"][inside]) / np.sum(mass))
        ds_dz = float((s[top] - s[bottom]) / (z[top] - z[bottom]))
        if not ds_dz > 0.0:
            raise ConfigurationError(
                "the rising column's troposphere is not stable in bulk, so no gravity "
                "wave relaxes it"
            )
        buoyancy_frequency = float(
            np.sqrt(self.gravity * ds_dz / (specific_heat * mean_temperature))
        )
        height = float(
            z[top] + self.gas_constant * T[top] / self.gravity * np.log(p[top] / self.tropopause)
        )
        wave_speed = buoyancy_frequency * height / np.pi
        return stratification, self.relaxation_distance / wave_speed

    def compute_tendencies(self) -> list[TendencyContribution]:
        ctx = self.context
        # Scheduled after both columns' radiation and before their water and
        # convection: this step's cooling in, this step's transport out for
        # them to compose with.
        up, down = self._column(self.ascent_region), self._column(self.descent_region)
        if up["T"].size != down["T"].size:
            raise ConfigurationError("the two columns must have the same layers")
        sigma = self.ascent_fraction
        stratification, self._relaxation_time = self._troposphere(up)
        # The sinking column's departure from the area mean; the rising
        # column's is -(1 - sigma) / sigma times it, which is why its flux is
        # the mirror scaled by the area ratio.
        flux_down = spectral_wtg_subsidence(
            temperature_anomaly=sigma * (down["T"] - up["T"]),
            specific_heat=down["C"] / down["m"],
            bulk_stratification=stratification,
            relaxation_time=self._relaxation_time,
            interface_pressure=down["pi"],
            tropopause=self.tropopause,
            gravity=self.gravity,
        )
        flux_up = -(1.0 - sigma) / sigma * flux_down
        self._flux = flux_down
        i500 = int(np.argmin(np.abs(down["pi"] - 5.0e4)))
        # Pressure velocity, Pa s^-1, positive downward (sinking).
        self._omega_500 = float(self.gravity * flux_down[i500])

        energy_down, lateral_s_down = upwind_tendency(flux_down, down["s"], up["s"])
        energy_up, _ = upwind_tendency(flux_up, up["s"], down["s"])
        water_down, lateral_q_down = upwind_tendency(flux_down, down["q"], up["q"])
        water_up, _ = upwind_tendency(flux_up, up["q"], down["q"])

        self._rates = {
            f"{self.ascent_region}:vapour_advection_rate": water_up,
            f"{self.descent_region}:vapour_advection_rate": water_down,
            f"{self.ascent_region}:advective_heating_rate": energy_up,
            f"{self.descent_region}:advective_heating_rate": energy_down,
        }
        # What crosses between the columns, per unit of global area over the step.
        energy_across = (1.0 - sigma) * float(np.sum(lateral_s_down)) * ctx.dt
        water_across = (1.0 - sigma) * float(np.sum(lateral_q_down)) * ctx.dt
        self._transfers = [
            *self._crossing("energy", ConservedQuantity.ENERGY, "J m^-2", energy_across, "s"),
            *self._crossing("water", ConservedQuantity.WATER, "kg m^-2", water_across, "q"),
            *self._crossing(
                "latent", ConservedQuantity.ENERGY, "J m^-2", water_across * self.latent_heat, "L"
            ),
        ]

        def contribution(
            region: str, state: str, value: np.ndarray, units: str, label: str
        ) -> TendencyContribution:
            return TendencyContribution(
                state_id=f"{region}:{state}",
                provider_module_id=self.module_id,
                label=label,
                value=value,
                units=units,
                equation_id=EQUATION_IDS["wtg_circulation"],
            )

        return [
            contribution(
                self.ascent_region,
                _TEMPERATURE,
                energy_up / up["C"],
                "K s^-1",
                "circulation_ascent",
            ),
            contribution(
                self.descent_region,
                _TEMPERATURE,
                energy_down / down["C"],
                "K s^-1",
                "circulation_subsidence",
            ),
            contribution(
                self.ascent_region, _VAPOUR, water_up, "kg m^-2 s^-1", "circulation_ascent"
            ),
            contribution(
                self.descent_region, _VAPOUR, water_down, "kg m^-2 s^-1", "circulation_subsidence"
            ),
        ]

    def _crossing(
        self,
        name: str,
        quantity: ConservedQuantity,
        units: str,
        amount: float,
        tag: str,
    ) -> list[Transfer]:
        """The net lateral exchange of one quantity, sender to receiver.

        Recorded by the circulation for both sides because it is the circulation
        that moves it: the columns' own modules never see the other column.
        """
        ctx = self.context
        owner = "atmosphere" if tag == "s" else "water"
        into_descent = amount >= 0.0
        sender_region = self.ascent_region if into_descent else self.descent_region
        receiver_region = self.descent_region if into_descent else self.ascent_region
        fraction = {
            self.ascent_region: self.ascent_fraction,
            self.descent_region: 1.0 - self.ascent_fraction,
        }
        return [
            Transfer(
                transfer_id=f"{self.module_id}_{name}@{ctx.time:g}",
                sender_module=f"{sender_region}/{owner}",
                receiver_module=f"{receiver_region}/{owner}",
                quantity=quantity,
                amount=abs(amount),
                units=units,
                start_time=ctx.time,
                end_time=ctx.time + ctx.dt,
                spatial_mapping="area_weighted_regions",
                equation_id=EQUATION_IDS["wtg_circulation"],
                channel=f"circulation_{name}",
                sender_area_fraction=fraction[sender_region],
                receiver_area_fraction=fraction[receiver_region],
            )
        ]

    def report_transfers(self) -> tuple[list[Transfer], list[Transfer]]:
        return list(self._transfers), list(self._transfers)

    def report_conservation(self) -> list[ModuleBudget]:
        return []

    def emit_outputs(self) -> dict[str, Field]:
        ctx = self.context
        provenance = Provenance(producer=self.module_id, derivation=EQUATION_IDS["wtg_circulation"])
        common = {
            "quantity_kind": QuantityKind.DIAGNOSTIC,
            "time_semantics": TimeSemantics.INTERVAL_AVERAGE,
            "time_interval": ctx.dt,
            "timestamp": ctx.time,
            "provenance": provenance,
        }
        rates = {
            name: Field(
                field_id=name,
                physical_name=name.split(":", 1)[1].replace("_", " "),
                values=np.array(value),
                units="kg m^-2 s^-1" if "vapour" in name else "W m^-2",
                spatial_support=SpatialSupport.LAYER_AVERAGE,
                dimensions=("layer",),
                **common,  # type: ignore[arg-type]
            )
            for name, value in self._rates.items()
        }
        return {
            **rates,
            "circulation_subsidence_mass_flux": Field(
                field_id="circulation_subsidence_mass_flux",
                physical_name="downward mass flux at the sinking column's interfaces",
                values=np.array(self._flux),
                units="kg m^-2 s^-1",
                spatial_support=SpatialSupport.LAYER_AVERAGE,
                dimensions=("interface",),
                **common,  # type: ignore[arg-type]
            ),
            "circulation_relaxation_time": Field(
                field_id="circulation_relaxation_time",
                physical_name="gravest mode's gravity-wave relaxation time",
                values=np.array(self._relaxation_time),
                units="s",
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                **common,  # type: ignore[arg-type]
            ),
            "circulation_subsidence_omega_500": Field(
                field_id="circulation_subsidence_omega_500",
                physical_name="sinking column's pressure velocity near 500 hPa",
                values=np.array(self._omega_500),
                units="Pa s^-1",
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                **common,  # type: ignore[arg-type]
            ),
        }

    def report_uncertainty(self) -> UncertaintyReport:
        return UncertaintyReport(
            module_id=self.module_id,
            quantified={"omega_500_hPa_day": self._omega_500 * 86400.0 / 100.0},
            method="none propagated",
            unresolved=(
                "two regions stand for a continuum of dynamical regimes",
                "the relaxation distance is one number for a spread of distances",
                "the tropopause the modes vanish at is declared, not derived",
            ),
        )

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "ascent_region": self.ascent_region,
            "descent_region": self.descent_region,
            "ascent_fraction": self.ascent_fraction,
            "ascent_fraction_source": (
                "NCEP/NCAR Reanalysis 1, omega at 500 hPa < 0, 1991-2020 monthly LTM"
            ),
            "relaxation_distance_m": self.relaxation_distance,
            "relaxation_distance_source": (
                "NCEP/NCAR Reanalysis 1 1991-2020 monthly LTM: mean distance from "
                "omega500 > 0 to the nearest omega500 < 0, within 20-30 degrees"
            ),
            "tropopause_Pa": self.tropopause,
        }


class AreaWeightedMeansModule(Module):
    """Publishes area-weighted global means of regional scalar fields (SPA-003).

    The integral over the globe of a per-area quantity is the area-weighted sum
    of the regions' values, so this is the regridding from regions to a global
    scalar that conserves the integral exactly. It lets every global diagnostic
    -- the error report above all -- read a two-column run unchanged.
    """

    def __init__(
        self,
        *,
        regions: Sequence[tuple[str, float]],
        fields: Sequence[tuple[str, str]],
        module_id: str = "global_means",
        status: ModuleStatus = ModuleStatus.DIAGNOSTIC,
    ) -> None:
        total = sum(a for _, a in regions)
        if abs(total - 1.0) > 1e-12:
            raise ConfigurationError(f"region areas {regions} must sum to one")
        contract = ModuleContract(
            module_id=module_id,
            module_name="area-weighted global means of regional fields",
            module_version=V1,
            interface_version=V1,
            fidelity_level=FidelityLevel.STANDARD,
            physical_basis="The global integral of a per-area quantity is the "
            "area-weighted sum over regions that partition the globe.",
            governing_equation_ids=(EQUATION_IDS["area_weighted_mean"],),
            assumptions=("the regions partition the globe",),
            validity_domain="Per-area and intensive scalar fields of regions that "
            "partition the globe.",
            required_inputs=tuple(f"{r}:{f}" for r, _ in regions for f, _ in fields),
            produced_outputs=tuple(f for f, _ in fields),
            conserved_quantities=(),
            correctness_tests=("test_the_global_mean_preserves_the_integral",),
            cost=CostDeclaration(
                module_id=module_id,
                time_complexity="O(n_fields)",
                memory_complexity="O(n_fields)",
                scaling_variables=("n_fields",),
                work_units_per_step=0.01,
                benchmark_configuration=f"{len(fields)} fields over {len(regions)} regions",
                peak_memory_bytes=4096,
                supported_resource_classes=(
                    ResourceClass.LAPTOP,
                    ResourceClass.ENGINEERING_WORKSTATION,
                ),
                lower_fidelity_alternatives=("report_each_region_separately",),
            ),
            fidelity_justification=(
                "Reporting regions separately would leave the global error "
                "undefined; the mean is the one comparison with a global observation."
            ),
        )
        super().__init__(contract, status=status)
        self.regions = tuple(regions)
        self.fields = tuple(fields)

    def initialize(self, initial_state: dict[str, Any], boundary_state: dict[str, Any]) -> None:
        self.mark_initialized()

    def compute_tendencies(self) -> list[TendencyContribution]:
        return []

    def emit_outputs(self) -> dict[str, Field]:
        ctx = self.context
        out: dict[str, Field] = {}
        for name, units in self.fields:
            parts = [(a, self._inputs.get(f"{r}:{name}")) for r, a in self.regions]
            present = [(a, f) for a, f in parts if f is not None]
            if len(present) != len(parts):
                continue
            value = sum(a * float(np.sum(f.require_units(units).values)) for a, f in present)
            out[name] = Field(
                field_id=name,
                physical_name=f"area-weighted global mean of {name}",
                values=np.array(value),
                units=units,
                quantity_kind=QuantityKind.DIAGNOSTIC,
                time_semantics=TimeSemantics.INTERVAL_AVERAGE,
                time_interval=ctx.dt,
                timestamp=ctx.time,
                spatial_support=SpatialSupport.GLOBAL_SCALAR,
                provenance=Provenance(
                    producer=self.module_id, derivation=EQUATION_IDS["area_weighted_mean"]
                ),
            )
        return out

    def report_uncertainty(self) -> UncertaintyReport:
        return UncertaintyReport(module_id=self.module_id, method="exact; no uncertainty added")

    def parameters(self) -> dict[str, Any]:
        return {
            **dict(self.configuration),
            "regions": [list(r) for r in self.regions],
            "fields": [list(f) for f in self.fields],
        }
