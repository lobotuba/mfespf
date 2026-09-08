"""Fundamental and reference constants with declared provenance.

Charter 4.13/4.14: a parameter must be classified, and a fundamental quantity
enters the framework at a declared precision and provenance.  Constants are
therefore records, not bare floats -- a run manifest can then state which value
of the Stefan-Boltzmann constant or the total solar irradiance it used and where
that number came from.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

from mfespf.core.errors import RegistryError

__all__ = ["CONSTANTS", "Constant", "ParameterClass", "get", "value_of"]


class ParameterClass(StrEnum):
    """Charter 4.13.  EQR-004 depends on this distinction being carried, not
    inferred: a calibrated correction must never be reported as a measured
    constant."""

    FUNDAMENTAL = "fundamental"
    MEASURED = "measured"
    DERIVED = "derived"
    PARAMETERIZED = "parameterized"
    CALIBRATED = "calibrated"


@dataclass(frozen=True, slots=True)
class Constant:
    name: str
    value: float
    units: str
    parameter_class: ParameterClass
    source: str
    standard_uncertainty: float | None = None
    notes: str = ""

    @property
    def relative_uncertainty(self) -> float | None:
        if self.standard_uncertainty is None or self.value == 0.0:
            return None
        return abs(self.standard_uncertainty / self.value)


def _c(
    name: str,
    value: float,
    units: str,
    parameter_class: ParameterClass,
    source: str,
    standard_uncertainty: float | None = None,
    notes: str = "",
) -> Constant:
    return Constant(name, value, units, parameter_class, source, standard_uncertainty, notes)


# CODATA 2018 values are exact by SI definition where marked; the astronomical
# and solar values are measured and carry uncertainty.
CONSTANTS: dict[str, Constant] = {
    c.name: c
    for c in (
        _c(
            "speed_of_light",
            2.99792458e8,
            "m s^-1",
            ParameterClass.FUNDAMENTAL,
            "SI definition (exact)",
            0.0,
        ),
        _c(
            "planck",
            6.62607015e-34,
            "J s",
            ParameterClass.FUNDAMENTAL,
            "SI definition (exact)",
            0.0,
        ),
        _c(
            "boltzmann",
            1.380649e-23,
            "J K^-1",
            ParameterClass.FUNDAMENTAL,
            "SI definition (exact)",
            0.0,
        ),
        _c(
            "stefan_boltzmann",
            5.670374419e-8,
            "W m^-2 K^-4",
            ParameterClass.DERIVED,
            "CODATA 2018, derived from exact h, k_B, c",
            0.0,
        ),
        _c(
            "gravitational",
            6.67430e-11,
            "m^3 kg^-1 s^-2",
            ParameterClass.MEASURED,
            "CODATA 2018",
            1.5e-15,
        ),
        _c(
            "astronomical_unit",
            1.495978707e11,
            "m",
            ParameterClass.FUNDAMENTAL,
            "IAU 2012 definition (exact)",
            0.0,
        ),
        _c(
            "earth_radius_mean",
            6.371e6,
            "m",
            ParameterClass.MEASURED,
            "IUGG mean radius R1",
            1.0e3,
            notes="Volumetric mean; a spherical-Earth approximation, declared as such.",
        ),
        _c(
            "earth_orbital_eccentricity_j2000",
            0.0167086,
            "1",
            ParameterClass.MEASURED,
            "VSOP87 / IAU J2000 mean element",
            1.0e-6,
            notes="Slowly varying; treat as time-dependent above millennial timescales.",
        ),
        _c(
            "earth_obliquity_j2000",
            0.409092804,
            "rad",
            ParameterClass.MEASURED,
            "IAU 2006 obliquity at J2000",
            1.0e-9,
        ),
        _c(
            "total_solar_irradiance_reference",
            1361.0,
            "W m^-2",
            ParameterClass.MEASURED,
            "TSIS-1/SORCE-era composite at 1 au, solar-cycle mean",
            0.5,
            notes="Reference value only. A run testing solar variability must use a "
            "time-varying dataset through the solar-source interface, not this.",
        ),
        _c(
            "water_latent_heat_vaporization_273",
            2.501e6,
            "J kg^-1",
            ParameterClass.MEASURED,
            "Standard thermodynamic tables at 273.15 K",
            1.0e3,
        ),
        _c(
            "water_latent_heat_fusion",
            3.337e5,
            "J kg^-1",
            ParameterClass.MEASURED,
            "Standard thermodynamic tables at 273.15 K",
            1.0e2,
        ),
        _c(
            "specific_heat_water_vapour",
            1885.0,
            "J kg^-1 K^-1",
            ParameterClass.MEASURED,
            "Standard thermodynamic tables at constant pressure, near 273 K",
            10.0,
            notes="Needed by Kirchhoff's relation dL/dT = c_pv - c_condensed, which "
            "is what makes the latent heat vary with temperature.",
        ),
        _c(
            "specific_heat_liquid_water",
            4218.0,
            "J kg^-1 K^-1",
            ParameterClass.MEASURED,
            "Standard thermodynamic tables at 273.15 K",
            10.0,
        ),
        _c(
            "specific_heat_ice",
            2106.0,
            "J kg^-1 K^-1",
            ParameterClass.MEASURED,
            "Standard thermodynamic tables at 273.15 K",
            10.0,
        ),
        _c(
            "specific_heat_seawater",
            3985.0,
            "J kg^-1 K^-1",
            ParameterClass.MEASURED,
            "TEOS-10 nominal value at typical surface salinity",
            5.0,
        ),
        _c(
            "density_seawater",
            1025.0,
            "kg m^-3",
            ParameterClass.MEASURED,
            "Nominal surface seawater density",
            2.0,
        ),
        _c(
            "earth_standard_gravity",
            9.80665,
            "m s^-2",
            ParameterClass.FUNDAMENTAL,
            "CGPM 1901 conventional value (exact)",
            0.0,
            notes="Conventional, not local. The gravity module (GRV-001) supersedes it.",
        ),
        _c(
            "specific_heat_dry_air",
            1004.6,
            "J kg^-1 K^-1",
            ParameterClass.MEASURED,
            "Standard thermodynamic tables near 273 K, isobaric",
            0.5,
            notes="Weakly temperature dependent; held constant at this fidelity.",
        ),
        _c(
            "gas_constant_dry_air",
            287.052874,
            "J kg^-1 K^-1",
            ParameterClass.DERIVED,
            "Universal gas constant over the standard dry-air molar mass",
            1.0e-5,
        ),
        _c(
            "standard_surface_pressure",
            101325.0,
            "Pa",
            ParameterClass.FUNDAMENTAL,
            "ISO 2533 standard atmosphere (exact by definition)",
            0.0,
            notes="A reference value, not Earth's global-mean surface pressure.",
        ),
        _c(
            "gas_constant_water_vapour",
            461.523,
            "J kg^-1 K^-1",
            ParameterClass.DERIVED,
            "Universal gas constant over the molar mass of water",
            1.0e-3,
        ),
        _c(
            "triple_point_temperature",
            273.16,
            "K",
            ParameterClass.FUNDAMENTAL,
            "ITS-90 definition of the water triple point (exact)",
            0.0,
        ),
        _c(
            "saturation_vapour_pressure_triple_point",
            611.657,
            "Pa",
            ParameterClass.MEASURED,
            "Measured saturation vapour pressure at the water triple point",
            0.01,
        ),
        _c(
            "radiative_diffusivity_factor",
            1.66,
            "1",
            ParameterClass.PARAMETERIZED,
            "Elsasser diffusivity approximation to the hemispheric transmission integral",
            0.05,
            notes="Replaces the angular integral of transmission with one effective slant "
            "path. Accurate to a few percent near unit optical depth and degrading "
            "outside that range; it is a parameterization, not a measured quantity.",
        ),
        _c(
            "ocean_deep_exchange_coefficient",
            0.7,
            "W m^-2 K^-1",
            ParameterClass.CALIBRATED,
            "Two-layer energy-balance fits to coupled-model step experiments",
            0.3,
            notes="Fitted to model behaviour, not measured on seawater. EQR-004 forbids "
            "presenting it as a measured constant, and the large declared uncertainty is "
            "the honest statement of how well the fit constrains it.",
        ),
    )
}


def get(name: str) -> Constant:
    try:
        return CONSTANTS[name]
    except KeyError:
        raise RegistryError(
            f"constant {name!r} is not registered; add it to mfespf.core.constants "
            f"with a declared parameter class and source before use"
        ) from None


def value_of(name: str) -> float:
    """Numeric value in the constant's declared units."""
    return get(name).value
