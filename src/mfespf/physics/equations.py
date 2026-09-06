"""Equation records for the Stage 0 and Stage 1 kernel.

Charter section 19 stages 0-1, requirements EQR-001 through EQR-005.

Every equation the reference kernel evaluates is registered here before a module
may cite it.  The categories matter as much as the formulae: the Planck law and
Kepler's equation are fundamental, the grey-body emission and slab energy
balance are reduced physical equations with stated approximations, and the
Newton solve and forward Euler step are numerical operators whose error is
accounted separately from physical uncertainty (UNC-003).
"""

from __future__ import annotations

from mfespf.core.versioning import Version
from mfespf.registries.equations import (
    EquationCategory,
    EquationRecord,
    EquationRegistry,
    ValidityDomain,
)

__all__ = ["EQUATION_IDS", "build_equation_registry"]

V1 = Version(1, 0, 0)

EQUATION_IDS = {
    "planck": "EQ-PLANCK-SPECTRAL-RADIANCE",
    "greybody": "EQ-GREYBODY-EMISSION",
    "inverse_square": "EQ-INVERSE-SQUARE-IRRADIANCE",
    "kepler": "EQ-KEPLER-ORBIT-RADIUS",
    "kepler_solver": "EQ-NUM-KEPLER-NEWTON",
    "insolation": "EQ-SPHERICAL-MEAN-INSOLATION",
    "spectral_absorption": "EQ-SPECTRAL-SURFACE-ABSORPTION",
    "slab": "EQ-SLAB-ENERGY-BALANCE",
    "forward_euler": "EQ-NUM-FORWARD-EULER",
    "equilibrium_temperature": "EQ-RADIATIVE-EQUILIBRIUM-TEMPERATURE",
    "boundary_accumulation": "EQ-BOUNDARY-ACCUMULATION",
    "hydrostatic": "EQ-HYDROSTATIC-LAYER-MASS",
    "layer_energy_balance": "EQ-LAYER-ENERGY-BALANCE",
    "schwarzschild": "EQ-SCHWARZSCHILD-TRANSFER",
    "beer_lambert": "EQ-BEER-LAMBERT-EXTINCTION",
    "diffusivity": "EQ-DIFFUSIVITY-APPROXIMATION",
    "grey_greenhouse": "EQ-GREY-SINGLE-LAYER-GREENHOUSE",
    "potential_temperature": "EQ-POTENTIAL-TEMPERATURE",
    "convective_adjustment": "EQ-DRY-CONVECTIVE-ADJUSTMENT",
    "clausius_clapeyron": "EQ-CLAUSIUS-CLAPEYRON-SATURATION",
    "specific_humidity": "EQ-SATURATION-SPECIFIC-HUMIDITY",
    "surface_evaporation": "EQ-BULK-SURFACE-EVAPORATION",
    "condensation": "EQ-SATURATION-ADJUSTMENT-CONDENSATION",
    "latent_energy": "EQ-LATENT-ENERGY-OF-VAPOUR",
    "sensible_heat": "EQ-BULK-SENSIBLE-HEAT-FLUX",
    "cloud_optics": "EQ-CLOUD-OPTICAL-PROPERTIES",
    "cloud_radiative_effect": "EQ-CLOUD-RADIATIVE-EFFECT",
    "two_stream": "EQ-TWO-STREAM-LAYER-COEFFICIENTS",
    "deep_ocean_exchange": "EQ-TWO-RESERVOIR-OCEAN-EXCHANGE",
    "co2_band_absorption": "EQ-CO2-BAND-ABSORPTION",
    "volcanic_aerosol_optics": "EQ-VOLCANIC-AEROSOL-OPTICS",
    "vapour_band_absorption": "EQ-VAPOUR-BAND-ABSORPTION",
    "bulk_plume_humidity": "EQ-BULK-PLUME-HUMIDITY",
    "entraining_lapse_rate": "EQ-ENTRAINING-LAPSE-RATE",
}


def build_equation_registry() -> EquationRegistry:
    """Register the Stage 0/1 equation set and return the registry."""
    reg = EquationRegistry()

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["planck"],
            version=V1,
            name="Planck spectral radiance",
            category=EquationCategory.FUNDAMENTAL_LAW,
            mathematical_form=(
                "B_lambda(T) = (2 h c^2 / lambda^5) / (exp(h c / (lambda k_B T)) - 1)"
            ),
            physical_basis=(
                "Quantised thermal emission from a body in local thermodynamic equilibrium."
            ),
            variables=(
                ("B_lambda", "spectral radiance", "W m^-3"),
                ("lambda", "wavelength", "m"),
                ("T", "thermodynamic temperature", "K"),
            ),
            parameters=(
                ("planck", "fundamental", "SI definition"),
                ("boltzmann", "fundamental", "SI definition"),
                ("speed_of_light", "fundamental", "SI definition"),
            ),
            assumptions=("local thermodynamic equilibrium", "isotropic emission into 2pi sr"),
            validity_domain=ValidityDomain(
                ranges=(("T", 1.0, 1.0e5, "K"), ("lambda", 1.0e-9, 1.0e-2, "m")),
                conditions=("non-degenerate, non-relativistic radiating matter",),
            ),
            limiting_cases=(
                "Rayleigh-Jeans as h c / (lambda k T) -> 0",
                "Wien as h c / (lambda k T) -> infinity",
                "integral over lambda and solid angle recovers sigma T^4",
            ),
            conserved_quantities=("energy",),
            correctness_tests=(
                "test_planck_integrates_to_stefan_boltzmann",
                "test_planck_wien_displacement",
            ),
            known_failure_modes=("numerical overflow in the exponential at small lambda T",),
            source="Planck 1901; CODATA 2018 constants",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["greybody"],
            version=V1,
            name="Grey-body hemispheric emission",
            category=EquationCategory.REDUCED_PHYSICAL_EQUATION,
            mathematical_form="M = epsilon * sigma * T^4",
            physical_basis=(
                "Planck emission integrated over wavelength and hemisphere, with a "
                "wavelength-independent emissivity."
            ),
            variables=(
                ("M", "emitted flux density", "W m^-2"),
                ("epsilon", "hemispheric emissivity", "1"),
                ("T", "surface temperature", "K"),
            ),
            parameters=(("stefan_boltzmann", "derived", "CODATA 2018"),),
            assumptions=("emissivity independent of wavelength and angle",),
            approximations=(
                "Replaces the spectral emission integral with a single grey coefficient; "
                "this is a modelling choice and is refuted whenever a spectral observable "
                "disagrees beyond tolerance.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("T", 1.0, 6000.0, "K"), ("epsilon", 0.0, 1.0, "1")),
            ),
            limiting_cases=("epsilon = 0 gives no emission", "epsilon = 1 gives a black body"),
            conserved_quantities=("energy",),
            permitted_discretizations=("pointwise evaluation",),
            correctness_tests=(
                "test_greybody_matches_planck_integral",
                "test_greybody_zero_emissivity_is_null",
            ),
            known_failure_modes=(
                "systematically wrong outgoing spectrum even when the broadband flux matches",
            ),
            supersedes=(),
            source="Standard radiative transfer; grey approximation stated explicitly",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["inverse_square"],
            version=V1,
            name="Inverse-square irradiance scaling",
            category=EquationCategory.FUNDAMENTAL_LAW,
            mathematical_form="S(d) = S_ref * (d_ref / d)^2",
            physical_basis="Conservation of radiant flux through concentric spheres in vacuum.",
            variables=(
                ("S", "irradiance normal to the beam", "W m^-2"),
                ("d", "distance from the source", "m"),
            ),
            assumptions=("point source", "no intervening absorption or scattering"),
            validity_domain=ValidityDomain(
                ranges=(("d", 1.0e9, 1.0e13, "m"),),
                conditions=("distance large compared with the solar radius",),
            ),
            limiting_cases=("S(d_ref) = S_ref exactly",),
            conserved_quantities=("energy",),
            correctness_tests=("test_inverse_square_scaling", "test_inverse_square_identity"),
            source="Geometric optics",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["kepler"],
            version=V1,
            name="Kepler orbit radius from mean anomaly",
            category=EquationCategory.FUNDAMENTAL_LAW,
            mathematical_form="M = E - e sin E ;  r = a (1 - e cos E)",
            physical_basis="Two-body Newtonian gravitation with the Sun-Earth reduced mass.",
            variables=(
                ("M", "mean anomaly", "rad"),
                ("E", "eccentric anomaly", "rad"),
                ("r", "instantaneous Sun-Earth distance", "m"),
            ),
            parameters=(
                ("earth_orbital_eccentricity_j2000", "measured", "IAU J2000 mean element"),
                ("astronomical_unit", "fundamental", "IAU 2012 definition"),
            ),
            assumptions=(
                "two-body Keplerian orbit",
                "orbital elements fixed over the simulated interval",
            ),
            approximations=(
                "Neglects lunar, planetary and relativistic perturbations; adequate for "
                "annual-mean insolation, not for precise ephemeris work.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("e", 0.0, 0.9, "1"),),
                conditions=("simulated interval short compared with element drift timescales",),
            ),
            limiting_cases=("e = 0 gives r = a for all M",),
            symmetries=("r(M) = r(-M)",),
            conserved_quantities=("energy", "angular_momentum"),
            correctness_tests=(
                "test_kepler_circular_limit",
                "test_kepler_perihelion_aphelion",
                "test_kepler_radius_symmetry",
            ),
            source="Kepler's equation; standard celestial mechanics",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["kepler_solver"],
            version=V1,
            name="Newton-Raphson solution of Kepler's equation",
            category=EquationCategory.NUMERICAL_OPERATOR,
            mathematical_form="E_{n+1} = E_n - (E_n - e sin E_n - M) / (1 - e cos E_n)",
            physical_basis="Iterative root finding; introduces no physics.",
            variables=(("E", "eccentric anomaly", "rad"),),
            assumptions=("e < 1 so the iteration is contractive from E_0 = M",),
            validity_domain=ValidityDomain(ranges=(("e", 0.0, 0.9, "1"),)),
            conserved_quantities=(),
            correctness_tests=("test_kepler_solver_residual_below_tolerance",),
            known_failure_modes=("slow convergence as e approaches 1",),
            source="Standard numerical celestial mechanics",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["insolation"],
            version=V1,
            name="Spherical-mean insolation geometry",
            category=EquationCategory.REDUCED_PHYSICAL_EQUATION,
            mathematical_form="S_toa_mean = S(d) * (pi R^2) / (4 pi R^2) = S(d) / 4",
            physical_basis=(
                "A sphere of radius R intercepts the beam over its projected disc and "
                "distributes it over its whole surface."
            ),
            variables=(
                ("S_toa_mean", "global-mean top-of-atmosphere irradiance", "W m^-2"),
                ("S", "irradiance normal to the beam at the planet", "W m^-2"),
            ),
            assumptions=("spherical planet", "parallel beam", "instantaneous global averaging"),
            approximations=(
                "Discards the diurnal and latitudinal distribution entirely. Any claim about "
                "regional or diurnal behaviour is outside this equation's domain.",
            ),
            validity_domain=ValidityDomain(
                conditions=("zero-dimensional global-mean configurations only",),
            ),
            limiting_cases=("intercepted power equals distributed power exactly",),
            conserved_quantities=("energy",),
            correctness_tests=(
                "test_intercepted_equals_distributed_power",
                "test_global_mean_factor_is_one_quarter",
            ),
            known_failure_modes=(
                "cannot represent polar or seasonal contrast; a claim requiring either must "
                "activate a latitude-resolved implementation",
            ),
            source="Elementary geometry",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["spectral_absorption"],
            version=V1,
            name="Spectral surface absorption",
            category=EquationCategory.DIAGNOSTIC_IDENTITY,
            mathematical_form="A = integral over lambda of (1 - alpha(lambda)) S(lambda) d lambda",
            physical_basis=(
                "Definition of spectral albedo as the reflected fraction per wavelength."
            ),
            variables=(
                ("A", "absorbed flux density", "W m^-2"),
                ("alpha", "spectral albedo", "1"),
                ("S", "incident spectral irradiance", "W m^-2 m^-1"),
            ),
            assumptions=("no transmission through the surface", "no wavelength redistribution"),
            validity_domain=ValidityDomain(ranges=(("alpha", 0.0, 1.0, "1"),)),
            limiting_cases=(
                "grey albedo reduces to A = (1 - alpha) * S_total for any spectral grid",
            ),
            conserved_quantities=("energy",),
            correctness_tests=(
                "test_grey_albedo_matches_broadband",
                "test_absorbed_plus_reflected_equals_incident",
            ),
            source="Definition",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["slab"],
            version=V1,
            name="Slab surface energy balance",
            category=EquationCategory.REDUCED_PHYSICAL_EQUATION,
            mathematical_form="C dT/dt = A - epsilon sigma T^4",
            physical_basis=(
                "First law applied to a well-mixed layer of areal heat capacity C exchanging "
                "energy only by radiation."
            ),
            variables=(
                ("C", "areal heat capacity", "J m^-2 K^-1"),
                ("T", "layer temperature", "K"),
                ("A", "absorbed flux density", "W m^-2"),
            ),
            assumptions=(
                "the layer is isothermal",
                "no conduction to a deeper reservoir",
                "no turbulent or latent exchange",
            ),
            approximations=(
                "Collapses the vertical thermal structure to one temperature; the response "
                "timescale it predicts is therefore a single exponential by construction, "
                "which is exactly the assumption the timescale experiments must test.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("T", 1.0, 1000.0, "K"), ("C", 1.0e3, 1.0e10, "J m^-2 K^-1")),
                conditions=("airless or transparent-atmosphere configurations",),
            ),
            limiting_cases=(
                "steady state gives T_eq = (A / (epsilon sigma))^(1/4)",
                "small perturbations relax with tau = C / (4 epsilon sigma T_eq^3)",
            ),
            conserved_quantities=("energy",),
            permitted_discretizations=("forward Euler", "implicit Euler", "RK4"),
            correctness_tests=(
                "test_slab_reaches_analytic_equilibrium",
                "test_slab_relaxation_timescale",
                "test_slab_energy_closure",
            ),
            known_failure_modes=(
                "forward Euler is unstable for dt greater than about 2 C / (4 epsilon sigma T^3)",
            ),
            source="First law of thermodynamics applied to a lumped reservoir",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["forward_euler"],
            version=V1,
            name="Forward Euler time integration",
            category=EquationCategory.NUMERICAL_OPERATOR,
            mathematical_form="y_{n+1} = y_n + dt * f(t_n, y_n)",
            physical_basis="Explicit first-order integration; introduces no physics.",
            variables=(("y", "integrated state", "varies"), ("dt", "time step", "s")),
            assumptions=("dt within the stability limit of the local Jacobian",),
            approximations=("local truncation error O(dt^2), global O(dt)",),
            validity_domain=ValidityDomain(
                conditions=("dt * |df/dy| < 2 for the linearised problem",),
            ),
            conserved_quantities=(),
            correctness_tests=(
                "test_forward_euler_first_order_convergence",
                "test_forward_euler_exact_energy_bookkeeping",
            ),
            known_failure_modes=("instability and unphysical states beyond the stability limit",),
            source="Standard numerical analysis",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["equilibrium_temperature"],
            version=V1,
            name="Radiative equilibrium temperature of an airless sphere",
            category=EquationCategory.DIAGNOSTIC_IDENTITY,
            mathematical_form="T_eq = ( (1 - alpha) S / (4 epsilon sigma) )^(1/4)",
            physical_basis="Steady state of the slab balance with spherical-mean insolation.",
            variables=(
                ("T_eq", "equilibrium temperature", "K"),
                ("S", "irradiance normal to the beam", "W m^-2"),
                ("alpha", "broadband albedo", "1"),
            ),
            assumptions=("grey albedo and emissivity", "no atmosphere", "steady state"),
            validity_domain=ValidityDomain(ranges=(("alpha", 0.0, 1.0, "1"),)),
            limiting_cases=("alpha -> 1 gives T_eq -> 0",),
            conserved_quantities=("energy",),
            correctness_tests=("test_analytic_equilibrium_matches_integration",),
            source="Closed form of EQ-SLAB-ENERGY-BALANCE at steady state",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["boundary_accumulation"],
            version=V1,
            name="System-boundary energy accumulation",
            category=EquationCategory.DIAGNOSTIC_IDENTITY,
            mathematical_form="E_out(t + dt) = E_out(t) + sum_channels F_channel * dt",
            physical_basis=(
                "Definition of a cumulative reservoir at the modeled system boundary. It "
                "introduces no process; it exists so the top-level budget can be closed."
            ),
            variables=(
                ("E_out", "cumulative energy leaving the system", "J m^-2"),
                ("F_channel", "outgoing flux density per channel", "W m^-2"),
            ),
            assumptions=("no return flux from beyond the boundary at this fidelity",),
            validity_domain=ValidityDomain(
                conditions=("accounting only; carries no predictive content",),
            ),
            conserved_quantities=("energy",),
            correctness_tests=("test_space_accumulates_exactly_what_surface_exports",),
            source="Definition",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["hydrostatic"],
            version=V1,
            name="Hydrostatic layer mass",
            category=EquationCategory.FUNDAMENTAL_LAW,
            mathematical_form="dp = -rho g dz  =>  m_layer = dp / g",
            physical_basis=(
                "Balance of the vertical pressure gradient against gravity in a fluid "
                "with no vertical acceleration."
            ),
            variables=(
                ("p", "pressure", "Pa"),
                ("m_layer", "layer mass per unit area", "kg m^-2"),
                ("g", "gravitational acceleration", "m s^-2"),
            ),
            parameters=(("earth_standard_gravity", "fundamental", "CGPM 1901"),),
            assumptions=(
                "no vertical acceleration",
                "gravity constant over the depth of the column",
            ),
            approximations=(
                "Constant gravity introduces an error of order z / R_earth: about 1.5 "
                "percent at 100 km and negligible through the troposphere.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("p", 0.0, 1.2e5, "Pa"),),
                conditions=("hydrostatic, non-accelerating column",),
            ),
            limiting_cases=("total column mass equals (p_surface - p_top) / g exactly",),
            conserved_quantities=("mass",),
            correctness_tests=("test_column_mass_is_hydrostatic",),
            known_failure_modes=(
                "invalid where vertical accelerations matter, such as deep convection, "
                "which this fidelity does not represent in any case",
            ),
            source="Standard atmospheric thermodynamics",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["layer_energy_balance"],
            version=V1,
            name="Atmospheric layer energy balance",
            category=EquationCategory.REDUCED_PHYSICAL_EQUATION,
            mathematical_form="c_p m_i dT_i/dt = H_i",
            physical_basis=(
                "First law applied to an isothermal layer of fixed mass, with the net "
                "radiative flux divergence as the only heating term."
            ),
            variables=(
                ("T_i", "layer temperature", "K"),
                ("m_i", "layer mass per unit area", "kg m^-2"),
                ("H_i", "net radiative heating of the layer", "W m^-2"),
            ),
            parameters=(("specific_heat_dry_air", "measured", "thermodynamic tables"),),
            assumptions=(
                "the layer is isothermal",
                "dry air with constant specific heat at constant pressure",
                "no advection, convection, latent heating or kinetic energy",
            ),
            approximations=(
                "Omitting convection makes the solution pure radiative equilibrium, which "
                "is convectively unstable in the real lower troposphere. That is a known "
                "wrong answer for the lapse rate and is declared rather than hidden.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("T_i", 100.0, 400.0, "K"),),
                conditions=("radiative-only column experiments",),
            ),
            limiting_cases=(
                "zero heating leaves the layer temperature unchanged",
                "the sum over layers equals the column net radiative gain",
            ),
            conserved_quantities=("energy",),
            permitted_discretizations=("forward Euler", "implicit Euler"),
            correctness_tests=(
                "test_column_energy_budget_closes",
                "test_zero_insolation_column_only_cools",
            ),
            known_failure_modes=(
                "produces a super-adiabatic lapse rate that convection would remove",
            ),
            source="First law of thermodynamics applied to a lumped layer",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["beer_lambert"],
            version=V1,
            name="Beer-Lambert extinction",
            category=EquationCategory.FUNDAMENTAL_LAW,
            mathematical_form="tau = kappa * m ;  transmittance = exp(-tau)",
            physical_basis=(
                "Exponential attenuation along a path through an absorbing medium, from "
                "the independence of successive absorption events."
            ),
            variables=(
                ("tau", "optical depth", "1"),
                ("kappa", "mass absorption coefficient", "m^2 kg^-1"),
                ("m", "absorber mass per unit area", "kg m^-2"),
            ),
            assumptions=("no scattering", "no stimulated emission along the path"),
            validity_domain=ValidityDomain(
                ranges=(("tau", 0.0, 100.0, "1"),),
                conditions=("absorbing, non-scattering medium",),
            ),
            limiting_cases=(
                "tau = 0 gives unit transmittance",
                "tau -> infinity gives zero transmittance",
            ),
            conserved_quantities=("energy",),
            correctness_tests=(
                "test_transparent_column_recovers_the_airless_result",
                "test_incident_equals_absorbed_plus_reflected_plus_transmitted",
            ),
            known_failure_modes=(
                "a band-mean transmittance is not the mean of the transmittances, so "
                "coarse bands bias the flux; the spectral convergence test measures it",
            ),
            source="Standard radiative transfer",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["schwarzschild"],
            version=V1,
            name="Schwarzschild two-stream transfer with emission",
            category=EquationCategory.REDUCED_PHYSICAL_EQUATION,
            mathematical_form=("F(i+1) = F(i) * (1 - e_i) + e_i * B_i ;  H_i = dF_up + dF_down"),
            physical_basis=(
                "Radiative transfer through an absorbing and emitting medium in local "
                "thermodynamic equilibrium, where emissivity equals absorptivity."
            ),
            variables=(
                ("F", "hemispheric flux at an interface", "W m^-2"),
                ("e_i", "layer emissivity in the band", "1"),
                ("B_i", "band-integrated Planck exitance of the layer", "W m^-2"),
                ("H_i", "net radiative heating of the layer", "W m^-2"),
            ),
            assumptions=(
                "plane-parallel, horizontally homogeneous layers",
                "no scattering",
                "local thermodynamic equilibrium, so Kirchhoff's law holds",
                "each layer emits at a single temperature",
            ),
            approximations=(
                "The hemispheric angular integral is replaced by one diffusivity factor.",
                "Layer emission uses the layer-mean temperature rather than integrating "
                "the Planck function across the layer's internal gradient.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("e_i", 0.0, 1.0, "1"), ("T", 100.0, 400.0, "K")),
                conditions=("non-scattering column in local thermodynamic equilibrium",),
            ),
            limiting_cases=(
                "zero opacity transmits the incident field unchanged and emits nothing",
                "one fully absorbing layer reproduces the analytic grey greenhouse",
                "flux divergence summed over layers equals the column net radiative gain",
            ),
            symmetries=("Kirchhoff: absorptivity equals emissivity in every band",),
            conserved_quantities=("energy",),
            permitted_discretizations=("layer-interface sweep, downward then upward",),
            correctness_tests=(
                "test_transparent_column_recovers_the_airless_result",
                "test_single_grey_layer_matches_the_analytic_greenhouse",
                "test_radiation_flux_divergence_closes_identically",
                "test_spectral_convergence_of_the_column_solution",
            ),
            known_failure_modes=(
                "cannot represent scattering at all, so a scattering atmosphere is out of "
                "scope rather than approximated",
                "coarse layers bias emission where the internal gradient is large",
            ),
            source="Schwarzschild 1906; standard two-stream formulation",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["diffusivity"],
            version=V1,
            name="Diffusivity approximation to the hemispheric transmission integral",
            category=EquationCategory.CLOSURE_OR_PARAMETERIZATION,
            mathematical_form=(
                "integral over the hemisphere of exp(-tau/mu) ~ exp(-D tau), D ~ 1.66"
            ),
            physical_basis=(
                "The exact hemispheric integral of exponential transmission is the "
                "exponential integral E3; one effective slant path reproduces it closely "
                "over the optical depths that matter."
            ),
            variables=(("D", "diffusivity factor", "1"), ("tau", "vertical optical depth", "1")),
            parameters=(("radiative_diffusivity_factor", "parameterized", "Elsasser 1942"),),
            assumptions=("isotropic radiance within each hemisphere",),
            approximations=(
                "Accurate to a few percent near unit optical depth, degrading toward both "
                "the optically thin and optically thick limits.",
            ),
            validity_domain=ValidityDomain(ranges=(("tau", 0.01, 10.0, "1"),)),
            limiting_cases=("D approaches 2 as tau -> 0 and 1.66 near unit optical depth",),
            conserved_quantities=("energy",),
            correctness_tests=("test_diffusivity_is_a_declared_parameterization",),
            known_failure_modes=(
                "biases heating rates in optically very thin or very thick layers",
            ),
            source="Elsasser 1942; standard in reduced radiation schemes",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["grey_greenhouse"],
            version=V1,
            name="Analytic single-layer grey greenhouse",
            category=EquationCategory.DIAGNOSTIC_IDENTITY,
            mathematical_form="T_s^4 = 2 T_e^4 / (2 - eps_a) ;  T_a^4 = T_s^4 / 2",
            physical_basis=(
                "Steady state of one isothermal grey layer transparent to sunlight that "
                "absorbs a fraction eps_a of the surface emission and radiates equally up "
                "and down."
            ),
            variables=(
                ("T_s", "surface temperature", "K"),
                ("T_a", "layer temperature", "K"),
                ("T_e", "planetary emission temperature", "K"),
                ("eps_a", "layer longwave emissivity", "1"),
            ),
            assumptions=(
                "a single layer",
                "the layer is transparent to the solar spectrum",
                "grey longwave emissivity",
                "steady state",
            ),
            validity_domain=ValidityDomain(ranges=(("eps_a", 0.0, 1.0, "1"),)),
            limiting_cases=(
                "eps_a = 0 recovers the airless result T_s = T_e",
                "eps_a = 1 gives T_s = 2^(1/4) T_e",
            ),
            conserved_quantities=("energy",),
            correctness_tests=("test_single_grey_layer_matches_the_analytic_greenhouse",),
            source="Closed form of EQ-SCHWARZSCHILD-TRANSFER for one grey layer",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["potential_temperature"],
            version=V1,
            name="Dry potential temperature",
            category=EquationCategory.DIAGNOSTIC_IDENTITY,
            mathematical_form="theta = T (p_ref / p)^(R_d / c_p)",
            physical_basis=(
                "Temperature a dry parcel would have if brought adiabatically to the "
                "reference pressure; conserved under dry adiabatic displacement."
            ),
            variables=(
                ("theta", "potential temperature", "K"),
                ("T", "temperature", "K"),
                ("p", "pressure", "Pa"),
            ),
            parameters=(
                ("gas_constant_dry_air", "derived", "molar mass of dry air"),
                ("specific_heat_dry_air", "measured", "thermodynamic tables"),
            ),
            assumptions=("dry air", "reversible adiabatic displacement"),
            validity_domain=ValidityDomain(
                ranges=(("T", 100.0, 400.0, "K"), ("p", 1.0, 1.2e5, "Pa")),
            ),
            limiting_cases=("theta equals T at the reference pressure",),
            symmetries=("invariant under dry adiabatic displacement",),
            conserved_quantities=(),
            correctness_tests=("test_exner_and_potential_temperature_are_inverse",),
            source="Standard atmospheric thermodynamics",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["convective_adjustment"],
            version=V1,
            name="Dry convective adjustment",
            category=EquationCategory.CLOSURE_OR_PARAMETERIZATION,
            mathematical_form=(
                "where theta decreases with height, set theta = "
                "sum(C_i T_i) / sum(C_i pi_i) over the unstable group; in the "
                "surface-based form the surface is element zero of the profile, with "
                "C_s and pi_s = pi(p_surface)"
            ),
            physical_basis=(
                "A column whose potential temperature decreases with height is unstable "
                "to vertical displacement; convection is assumed to have mixed such "
                "layers to neutrality within the coupling interval. The instability that "
                "drives it is usually between the ground and the air above rather than "
                "between two layers of air, so the surface is admitted to the profile as "
                "its bottom element and the same criterion decides whether it convects."
            ),
            variables=(
                ("theta", "potential temperature", "K"),
                ("C_i", "areal heat capacity of a layer", "J m^-2 K^-1"),
                ("pi_i", "Exner function at the layer", "1"),
            ),
            assumptions=(
                "mixing is instantaneous on the coupling timestep",
                "mixing is complete inside an unstable group and absent outside it",
                "the neutral target is the dry adiabat",
                "in the surface-based form the surface is a well-mixed slab that exchanges "
                "with the air by the same mixing, so no exchange coefficient sets the rate",
            ),
            approximations=(
                "Replaces resolved convective motion with an instantaneous adjustment. "
                "It has no timescale, so the answer does not converge as the timestep "
                "refines -- it is already independent of it, which is a limitation and "
                "not a strength.",
                "Targets the dry adiabat, over-stabilising a saturated troposphere whose "
                "real neutral profile is the smaller moist lapse rate.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("T", 100.0, 400.0, "K"),),
                conditions=("hydrostatic columns without resolved vertical motion",),
            ),
            limiting_cases=(
                "a stable profile is left exactly unchanged",
                "a surface colder than the air above it takes no part in the mixing",
                "a wholly unstable column mixes to a single potential temperature",
                "column enthalpy sum(C_i T_i) is conserved exactly",
            ),
            conserved_quantities=("energy",),
            permitted_discretizations=("iterative group merging from the lowest instability",),
            correctness_tests=(
                "test_convective_adjustment_conserves_column_enthalpy_exactly",
                "test_convective_adjustment_leaves_a_stable_profile",
                "test_an_already_stable_profile_is_left_untouched",
                "test_convection_adds_no_energy_to_the_column",
                "test_the_surface_mixing_conserves_the_enthalpy_of_the_surface_and_the_air",
                "test_an_inversion_between_layers_cannot_hide_a_buoyant_surface",
                "test_a_surface_colder_than_the_air_is_left_alone_and_bad_geometry_refused",
            ),
            known_failure_modes=(
                "predicts a dry lapse rate where the real atmosphere follows a moist one",
                "the layer-only form compares each layer with its neighbour, so an "
                "inversion at the bottom of the column makes it inactive however buoyant "
                "the surface underneath is -- measured as a 50 K potential-temperature "
                "step that the surface-based form removes entirely",
                "the surface-based form mixes completely rather than entraining, so it is "
                "an upper bound on how much the surface moistens the air above it",
                "no convective timescale, so it cannot represent the diurnal cycle of "
                "convection or any transient adjustment",
            ),
            source="Manabe and Strickler 1964; Manabe and Wetherald 1967",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["bulk_plume_humidity"],
            version=V1,
            name="Bulk-plume free-tropospheric relative humidity",
            category=EquationCategory.CLOSURE_OR_PARAMETERIZATION,
            mathematical_form=(
                "RH = (delta + alpha gamma - alpha epsilon) / (delta + gamma - alpha epsilon); "
                "with constant mass flux epsilon = delta and RH = (a + alpha) / (1 + a), "
                "a = epsilon (1 - alpha) / gamma"
            ),
            physical_basis=(
                "The environment's water budget in radiative-convective equilibrium. "
                "Subsidence between clouds carries dry air downward and dries the "
                "environment; detrainment of saturated plume air and the evaporation of "
                "condensate moisten it. The balance of the two fixes the relative "
                "humidity, which is therefore a property of how convection turns over "
                "and not a number that has to be imposed."
            ),
            variables=(
                ("RH", "environmental relative humidity", "1"),
                ("epsilon", "fractional entrainment rate", "m^-1"),
                ("delta", "fractional detrainment rate", "m^-1"),
                ("gamma", "water-vapour lapse rate, -d ln q* / dz", "m^-1"),
                ("alpha", "fraction of condensate that evaporates", "1"),
                ("a", "dimensionless bulk-plume parameter", "1"),
            ),
            parameters=(
                (
                    "bulk_plume_parameter",
                    "parameterized",
                    "0.2 for Earth's tropics, from Romps (2016); not fitted here, and "
                    "never to be fitted against a temperature record (INF-003)",
                ),
                (
                    "evaporating_condensate_fraction",
                    "parameterized",
                    "0.67, from measured precipitation efficiencies of 19 to 38 percent "
                    "(Langhans, Yeo and Romps 2015; Braham 1952); asserted over 0.62 to "
                    "0.81. Replaced an assumed 0.3 over 0.1 to 0.5, a range that "
                    "excluded the measured value rather than bracketing it",
                ),
            ),
            assumptions=(
                "radiative-convective equilibrium, so the environment's water budget is "
                "in steady state",
                "one bulk plume stands for the whole cloud ensemble",
                "the plume is saturated above cloud base",
                "constant mass flux with height in the closed form, which forces "
                "entrainment to equal detrainment",
                "the evaporating fraction is uniform in height",
            ),
            approximations=(
                "Derived for a convecting column and applied here to a global-mean "
                "one. The bias runs the other way from what that phrasing suggests: "
                "the closure reproduces the humidity of a radiative-convective column, "
                "near eighty percent, while the global-mean free troposphere is far "
                "drier because it is dominated by subtropical subsidence that no bulk "
                "plume represents. So this *overstates* the global-mean humidity, by an "
                "amount the framework cannot estimate from within.",
                "The constant-mass-flux form gives one relative humidity for the whole "
                "column, so the observed C shape -- moist near the surface and near the "
                "tropopause, driest in the middle -- is absent. Height-varying rates "
                "recover it and are supported by the general form.",
                "Condensate that evaporates is assumed to do so at the height it "
                "formed, so falling rain moistening the layers beneath it is absent.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("a", 0.0, 2.0, "1"), ("alpha", 0.0, 0.9, "1")),
                conditions=(
                    "convecting columns in radiative-convective equilibrium, where "
                    "cumulus mass flux rather than large-scale advection sets the "
                    "free-tropospheric humidity",
                ),
            ),
            limiting_cases=(
                "undiluted convection (a = 0) leaves the environment at the humidity "
                "evaporating condensate alone supplies, RH = alpha",
                "infinitely dilute convection saturates the environment, because "
                "everything the plume carries is detrained before it can precipitate",
                "the general form reduces exactly to the closed form when entrainment "
                "equals detrainment and tracks the water-vapour lapse rate",
            ),
            conserved_quantities=("water",),
            correctness_tests=(
                "test_the_closed_form_is_the_general_form_at_constant_mass_flux",
                "test_undiluted_convection_leaves_only_the_evaporated_condensate",
                "test_relative_humidity_rises_with_entrainment_and_saturates",
                "test_the_predicted_humidity_lands_in_the_observed_range",
            ),
            known_failure_modes=(
                "applied to a global-mean column it predicts the humidity of the "
                "convecting regions and misses the subtropical dry zones, so such a "
                "column will be too moist",
                "the closed form cannot produce the observed C-shaped profile at all",
                "the entrainment rate it implies sits at 0.198 per kilometre, the "
                "bottom edge of the 0.2 to 0.8 that cloud-resolving simulations "
                "diagnose for bulk plumes rather than comfortably inside it",
            ),
            source=(
                "Romps, D. M. (2014), An analytical model for tropical relative "
                "humidity, J. Climate 27(19), 7432-7449, doi:10.1175/JCLI-D-14-00255.1; "
                "closed form from Romps, D. M. (2016), Clausius-Clapeyron scaling of "
                "CAPE from analytical solutions to RCE, J. Atmos. Sci. 73(9), 3719-3737, "
                "doi:10.1175/JAS-D-15-0327.1"
            ),
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["entraining_lapse_rate"],
            version=V1,
            name="Zero-buoyancy entraining plume lapse rate",
            category=EquationCategory.CLOSURE_OR_PARAMETERIZATION,
            mathematical_form=(
                "Gamma = (g/c_p) [1 + a + q* L/(R T)] / [1 + a + q* L^2/(c_p R_v T^2)]"
            ),
            physical_basis=(
                "A plume that entrains subsaturated environmental air loses moist static "
                "energy as it rises, so it is neutrally buoyant against a profile more "
                "stable than a moist adiabat. The stratification convection leaves "
                "behind therefore depends on how dry the air it entrains is, and the gap "
                "between this profile and the undilute moist adiabat is what appears as "
                "convective available potential energy."
            ),
            variables=(
                ("Gamma", "lapse rate", "K m^-1"),
                ("T", "temperature", "K"),
                ("q_s", "saturation specific humidity", "kg kg^-1"),
                ("a", "dimensionless bulk-plume parameter", "1"),
            ),
            parameters=(
                (
                    "bulk_plume_parameter",
                    "parameterized",
                    "the same a as the humidity closure, which is what makes the pair "
                    "over-determined: fix it against one observable and the other is a "
                    "prediction that can fail",
                ),
            ),
            assumptions=(
                "the plume is exactly neutrally buoyant with respect to its environment",
                "entrainment is fractional and the entrained air is environmental",
                "virtual temperature effects are neglected in the buoyancy",
            ),
            approximations=(
                "Zero buoyancy is a limit, not a statement that convection has no "
                "buoyancy; real plumes carry some, so this slightly overstates the "
                "stability entrainment produces.",
                "A constant latent heat, inherited from the saturation relation.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("T", 200.0, 350.0, "K"), ("a", 0.0, 2.0, "1")),
                conditions=("saturated convecting columns above cloud base",),
            ),
            limiting_cases=(
                "a = 0 reduces exactly to the moist adiabatic lapse rate",
                "a dry column gives the dry adiabat at any entrainment rate, so "
                "entrainment cannot matter where there is no water",
                "the result lies between the moist and dry adiabats for any positive a",
            ),
            conserved_quantities=("energy",),
            correctness_tests=(
                "test_the_entraining_lapse_rate_reduces_to_the_moist_adiabat",
                "test_entrainment_stabilises_between_the_moist_and_dry_adiabats",
                "test_a_dry_column_is_dry_adiabatic_whatever_the_entrainment",
            ),
            known_failure_modes=(
                "it predicts a stratification but not the height convection reaches, so "
                "it cannot on its own say where the troposphere ends",
                "below cloud base there is no saturated plume and the relation does not apply",
            ),
            source=(
                "Singh, M. S. and O'Gorman, P. A. (2013), Influence of entrainment on "
                "the thermal stratification in simulations of radiative-convective "
                "equilibrium, Geophys. Res. Lett. 40, 4398-4403, doi:10.1002/grl.50796; "
                "closed form from Romps (2016), doi:10.1175/JAS-D-15-0327.1"
            ),
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["clausius_clapeyron"],
            version=V1,
            name="Clausius-Clapeyron saturation vapour pressure",
            category=EquationCategory.CONSTITUTIVE_RELATION,
            mathematical_form=("e_s(T) = e_triple exp( (L_v / R_v) (1/T_triple - 1/T) )"),
            physical_basis=(
                "Equilibrium between liquid water and its vapour along the phase "
                "boundary, integrated at constant latent heat."
            ),
            variables=(
                ("e_s", "saturation vapour pressure over liquid", "Pa"),
                ("T", "temperature", "K"),
            ),
            parameters=(
                ("water_latent_heat_vaporization_273", "measured", "thermodynamic tables"),
                ("gas_constant_water_vapour", "derived", "molar mass of water"),
                (
                    "saturation_vapour_pressure_triple_point",
                    "measured",
                    "measured triple-point value",
                ),
            ),
            assumptions=(
                "the latent heat of vaporisation is constant over the range used",
                "the vapour behaves as an ideal gas",
                "the liquid volume is negligible beside the vapour volume",
                "saturation is over liquid water, not ice",
            ),
            approximations=(
                "Holding L constant costs about a percent between 250 and 300 K and more "
                "outside it. Over ice below the triple point the true saturation pressure "
                "is lower than this returns, which matters for any claim about cirrus or "
                "polar humidity.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("T", 230.0, 340.0, "K"),),
                conditions=("saturation over liquid water",),
            ),
            limiting_cases=("returns the measured value exactly at the triple point",),
            conserved_quantities=(),
            correctness_tests=(
                "test_saturation_vapour_pressure_at_the_triple_point",
                "test_saturation_vapour_pressure_matches_reference_values",
            ),
            known_failure_modes=(
                "systematically wrong over ice",
                "extrapolates poorly far from the triple point",
            ),
            source="Clausius-Clapeyron relation, integrated form",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["specific_humidity"],
            version=V1,
            name="Saturation specific humidity",
            category=EquationCategory.DIAGNOSTIC_IDENTITY,
            mathematical_form="q_s = eps e_s / (p - (1 - eps) e_s),  eps = R_d / R_v",
            physical_basis="Definition of specific humidity applied at saturation.",
            variables=(
                ("q_s", "saturation specific humidity", "kg kg^-1"),
                ("e_s", "saturation vapour pressure", "Pa"),
                ("p", "total pressure", "Pa"),
            ),
            assumptions=("ideal-gas mixture of dry air and vapour",),
            approximations=(
                "The thin form q_s = eps e_s / p is not used: the two differ by about "
                "2 percent at 300 K and sea level.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("p", 1.0, 1.2e5, "Pa"),),
                conditions=("e_s well below the total pressure",),
            ),
            limiting_cases=("q_s -> 0 as the temperature falls",),
            conserved_quantities=(),
            correctness_tests=(
                "test_saturation_specific_humidity_is_monotone_and_bounded",
                "test_specific_humidity_uses_the_full_form_not_the_thin_approximation",
            ),
            known_failure_modes=("undefined once e_s approaches the total pressure",),
            source="Definition",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["surface_evaporation"],
            version=V1,
            name="Bulk surface evaporation",
            category=EquationCategory.CLOSURE_OR_PARAMETERIZATION,
            mathematical_form="E = k_E * max(0, q_s(T_s) - q_1)",
            physical_basis=(
                "Turbulent transfer of vapour down the humidity gradient between a "
                "saturated surface and the air above it."
            ),
            variables=(
                ("E", "evaporative mass flux", "kg m^-2 s^-1"),
                ("k_E", "bulk exchange coefficient", "kg m^-2 s^-1"),
                ("q_1", "specific humidity of the lowest layer", "kg kg^-1"),
            ),
            parameters=(("k_E", "parameterized", "free parameter of the reduced model"),),
            assumptions=(
                "the surface is saturated and has unlimited water available",
                "the exchange coefficient absorbs wind speed and stability entirely",
                "no condensation onto the surface, so the flux is one-way",
            ),
            approximations=(
                "A real bulk formula scales with wind speed and a stability-dependent "
                "transfer coefficient. Collapsing both into one constant means this "
                "cannot represent any wind or stability dependence of evaporation, and "
                "no claim about either is in scope.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("T_s", 230.0, 340.0, "K"),),
                conditions=("saturated surface with unlimited water",),
            ),
            limiting_cases=(
                "zero when the air above is already saturated at the surface temperature",
            ),
            conserved_quantities=("water", "energy"),
            correctness_tests=(
                "test_evaporation_conserves_water_and_carries_its_latent_energy",
                "test_evaporation_stops_when_the_air_is_saturated",
            ),
            known_failure_modes=(
                "no wind or stability dependence",
                "cannot dry out: the surface store is not allowed to limit the flux at "
                "this fidelity, so the run must declare an unlimited reservoir",
            ),
            source="Reduced bulk aerodynamic formula",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["condensation"],
            version=V1,
            name="Saturation adjustment condensation",
            category=EquationCategory.CLOSURE_OR_PARAMETERIZATION,
            mathematical_form="C_i = max(0, w_i - q_s(T_i, p_i) m_i)",
            physical_basis=(
                "Vapour in excess of saturation condenses; the latent heat released warms "
                "the layer and the condensate leaves as precipitation."
            ),
            variables=(
                ("C_i", "condensed mass in a layer", "kg m^-2"),
                ("w_i", "vapour mass in a layer", "kg m^-2"),
                ("m_i", "dry air mass of the layer", "kg m^-2"),
            ),
            assumptions=(
                "condensation is instantaneous once saturation is exceeded",
                "all condensate precipitates immediately, so no cloud water is stored",
                "no supersaturation and no nucleation barrier",
            ),
            approximations=(
                "Removing all condensate at once means there is no cloud, and therefore "
                "no cloud radiative effect at all. A claim about clouds is outside this "
                "module's scope until the Stage 4 cloud state exists.",
                "The saturation check uses the temperature before latent release rather "
                "than iterating to a consistent post-release state.",
            ),
            validity_domain=ValidityDomain(
                conditions=(
                    "columns where cloud storage and its radiative effect can be neglected",
                ),
            ),
            limiting_cases=("no condensation in subsaturated air",),
            conserved_quantities=("water", "energy"),
            correctness_tests=(
                "test_condensation_conserves_water_and_energy_together",
                "test_moist_column_budgets_close_for_water_and_energy",
            ),
            known_failure_modes=(
                "no cloud, so no cloud radiative effect",
                "instantaneous adjustment cannot represent precipitation timing",
            ),
            source="Standard saturation adjustment",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["latent_energy"],
            version=V1,
            name="Latent energy content of atmospheric vapour",
            category=EquationCategory.DIAGNOSTIC_IDENTITY,
            mathematical_form="E_latent = L_v * sum_i w_i",
            physical_basis=(
                "Energy that was supplied to evaporate the vapour and is returned when it "
                "condenses. Treating it as a stored quantity is what lets the water and "
                "energy budgets close simultaneously between the two events."
            ),
            variables=(
                ("E_latent", "latent energy per unit area", "J m^-2"),
                ("w_i", "vapour mass in a layer", "kg m^-2"),
            ),
            parameters=(
                ("water_latent_heat_vaporization_273", "measured", "thermodynamic tables"),
            ),
            assumptions=("constant latent heat of vaporisation",),
            approximations=(
                "The sensible heat carried by the water itself is neglected, both in the "
                "vapour and in the precipitation that returns to the surface.",
            ),
            validity_domain=ValidityDomain(ranges=(("T", 230.0, 340.0, "K"),)),
            limiting_cases=("zero in a dry column",),
            conserved_quantities=("energy",),
            correctness_tests=("test_moist_column_budgets_close_for_water_and_energy",),
            source="Definition",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["sensible_heat"],
            version=V1,
            name="Bulk sensible heat flux",
            category=EquationCategory.CLOSURE_OR_PARAMETERIZATION,
            mathematical_form="H = k_H (T_s - T_1)",
            physical_basis=(
                "Turbulent transport of sensible heat down the temperature gradient "
                "between the surface and the air immediately above it."
            ),
            variables=(
                ("H", "sensible heat flux", "W m^-2"),
                ("k_H", "bulk exchange coefficient", "W m^-2 K^-1"),
                ("T_1", "temperature of the lowest layer", "K"),
            ),
            parameters=(("k_H", "parameterized", "free parameter of the reduced model"),),
            assumptions=(
                "the exchange coefficient absorbs wind speed and stability entirely",
                "the flux is linear in the temperature difference",
            ),
            approximations=(
                "A real bulk formula scales with wind speed and a stability-dependent "
                "coefficient, and is strongly asymmetric between stable and unstable "
                "conditions. One constant can represent none of that, so no claim about "
                "the wind or stability dependence of the surface energy balance is in "
                "scope.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("T_s", 150.0, 400.0, "K"),),
                conditions=("columns without resolved boundary-layer turbulence",),
            ),
            limiting_cases=("zero flux when the surface and the air above are isothermal",),
            conserved_quantities=("energy",),
            correctness_tests=(
                "test_sensible_heat_flux_vanishes_when_isothermal",
                "test_sensible_heat_is_transferred_not_created",
            ),
            known_failure_modes=(
                "no wind or stability dependence",
                "linear in the temperature difference, where the real flux is not",
            ),
            source="Reduced bulk aerodynamic formula",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["two_stream"],
            version=V1,
            name="Two-stream layer reflectance, transmittance and emissivity",
            category=EquationCategory.REDUCED_PHYSICAL_EQUATION,
            mathematical_form=(
                "r = g2 sinh(k tau) / (k cosh(k tau) + g1 sinh(k tau)) ; "
                "t = k / (k cosh(k tau) + g1 sinh(k tau)) ; e = 1 - r - t"
            ),
            physical_basis=(
                "Radiative transfer through a homogeneous slab that absorbs, scatters "
                "and emits, with the radiance reduced to one upward and one downward "
                "hemispheric stream."
            ),
            variables=(
                ("r", "layer reflectance", "1"),
                ("t", "layer transmittance", "1"),
                ("e", "layer emissivity", "1"),
                ("tau", "optical depth", "1"),
                ("omega", "single-scattering albedo", "1"),
                ("g", "asymmetry parameter", "1"),
            ),
            parameters=(("radiative_diffusivity_factor", "parameterized", "Elsasser 1942"),),
            assumptions=(
                "isotropic radiance within each hemisphere",
                "homogeneous layer properties",
                "local thermodynamic equilibrium, so emissivity is one minus reflectance "
                "minus transmittance",
            ),
            approximations=(
                "Two streams cannot represent an anisotropic radiance field. Delta "
                "scaling removes the forward-scattering peak first, without which a "
                "cloud of asymmetry 0.85 would reflect far too much.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("tau", 0.0, 1000.0, "1"), ("omega", 0.0, 1.0, "1")),
            ),
            limiting_cases=(
                "omega = 0 gives r = 0, t = exp(-D tau), e = 1 - t: the absorption-only "
                "solution it replaces",
                "omega = 1 gives r + t = 1 and e = 0: a purely scattering layer neither "
                "absorbs nor emits",
                "tau = 0 gives a transparent layer",
            ),
            symmetries=("Kirchhoff: emissivity equals absorptivity",),
            conserved_quantities=("energy",),
            permitted_discretizations=("adding method over layers",),
            correctness_tests=(
                "test_two_stream_reduces_to_pure_absorption",
                "test_conservative_scattering_neither_absorbs_nor_emits",
                "test_scattering_solver_reproduces_the_absorption_only_results",
            ),
            known_failure_modes=(
                "two streams misrepresent strongly anisotropic fields even after delta scaling",
                "no polarisation and no angular detail",
            ),
            source="Standard two-stream transfer; Joseph, Wiscombe and Weinman 1976",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["cloud_optics"],
            version=V1,
            name="Cloud optical properties from water path",
            category=EquationCategory.CLOSURE_OR_PARAMETERIZATION,
            mathematical_form=(
                "shortwave: tau_sca = 3 W / (2 rho r_eff) ; longwave: tau_abs = k_lw W"
            ),
            physical_basis=(
                "Geometric-optics extinction by droplets large compared with the solar "
                "wavelength, and absorption by droplets comparable with the infrared "
                "wavelength."
            ),
            variables=(
                ("W", "condensed water path", "kg m^-2"),
                ("r_eff", "effective radius", "m"),
                ("tau_sca", "shortwave scattering optical depth", "1"),
                ("tau_abs", "longwave absorption optical depth", "1"),
            ),
            parameters=(
                ("k_lw", "parameterized", "reduced longwave mass absorption coefficient"),
                ("r_eff", "parameterized", "one prescribed effective radius per phase"),
                ("g", "parameterized", "droplet asymmetry parameter"),
            ),
            assumptions=(
                "one effective radius per phase, with no size distribution",
                "shortwave scattering is conservative; droplet absorption is neglected",
                "longwave scattering is neglected",
                "ice differs from liquid only through a larger effective radius",
            ),
            approximations=(
                "The geometric-optics limit is good to a few percent for droplets in the "
                "solar spectrum and poor for ice crystals, whose habit dominates their "
                "scattering and is absent here.",
                "A single effective radius replaces the distribution. The shortwave "
                "optical depth is inversely proportional to it, so that choice sets the "
                "cloud albedo as strongly as the water path does.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("r_eff", 1.0e-6, 1.0e-4, "m"), ("W", 0.0, 10.0, "kg m^-2")),
                conditions=("liquid or simple ice clouds, no mixed phase",),
            ),
            limiting_cases=(
                "zero water path gives a transparent cloud and no radiative effect",
                "optical depth is linear in water path at fixed effective radius",
            ),
            conserved_quantities=("energy",),
            correctness_tests=(
                "test_cloud_optical_depth_scales_with_water_path",
                "test_no_cloud_control_matches_the_clear_column",
                "test_a_thin_cloud_is_longwave_dominated_and_a_thick_one_shortwave",
            ),
            known_failure_modes=(
                "ice optical properties are a crude analogy to liquid",
                "no sub-grid variability, so a broken cloud field would be "
                "systematically too reflective",
            ),
            source="Standard reduced cloud optics",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["cloud_radiative_effect"],
            version=V1,
            name="Cloud radiative effect",
            category=EquationCategory.DIAGNOSTIC_IDENTITY,
            mathematical_form=(
                "CRE_sw = F_sw,net(all sky) - F_sw,net(clear) ; "
                "CRE_lw = OLR(clear) - OLR(all sky) ; CRE = CRE_sw + CRE_lw"
            ),
            physical_basis=(
                "Difference between the radiation budget of the column as it is and of "
                "the same column with the cloud removed, evaluated at the top."
            ),
            variables=(
                ("CRE_sw", "shortwave cloud radiative effect", "W m^-2"),
                ("CRE_lw", "longwave cloud radiative effect", "W m^-2"),
            ),
            assumptions=(
                "the clear-sky column differs from the all-sky one only by the cloud",
                "maximum overlap, with all-sky formed as a linear mix of a clear and an "
                "overcast column weighted by the cloud fraction",
            ),
            approximations=(
                "Independent-column mixing ignores horizontal photon transport between "
                "the clear and cloudy parts of a partially covered column.",
            ),
            validity_domain=ValidityDomain(ranges=(("f_c", 0.0, 1.0, "1"),)),
            limiting_cases=(
                "zero cloud fraction gives exactly zero in both components",
                "the shortwave component is negative and the longwave positive for a "
                "cloud that reflects sunlight and absorbs terrestrial radiation",
            ),
            conserved_quantities=("energy",),
            correctness_tests=(
                "test_no_cloud_control_matches_the_clear_column",
                "test_cloud_shortwave_effect_cools_and_longwave_effect_warms",
                "test_cloud_radiative_effect_is_reported_before_any_net",
            ),
            known_failure_modes=(
                "reporting only the net hides which component is wrong, since the two "
                "are of opposite sign and comparable magnitude",
            ),
            source="Standard definition; Ramanathan et al. 1989",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["deep_ocean_exchange"],
            version=V1,
            name="Two-reservoir ocean heat exchange",
            category=EquationCategory.CLOSURE_OR_PARAMETERIZATION,
            mathematical_form="H = gamma (T_s - T_d) ; C_d dT_d/dt = H",
            physical_basis=(
                "Bulk representation of every process that carries heat from the ocean "
                "mixed layer into the water below it -- subduction, diffusion and "
                "overturning -- as a single rate proportional to the temperature "
                "difference between two well-mixed reservoirs."
            ),
            variables=(
                ("H", "downward heat flux from the mixed layer", "W m^-2"),
                ("T_s", "mixed-layer temperature", "K"),
                ("T_d", "deep-reservoir temperature", "K"),
                ("C_d", "areal heat capacity of the deep reservoir", "J m^-2 K^-1"),
            ),
            parameters=(
                ("ocean_deep_exchange_coefficient", "calibrated", "fitted, not measured"),
                ("density_seawater", "measured", "nominal surface seawater density"),
                ("specific_heat_seawater", "measured", "TEOS-10 nominal value"),
            ),
            assumptions=(
                "both reservoirs are isothermal and well mixed",
                "the exchange coefficient does not depend on the state or on time",
                "no salinity, circulation or meridional transport",
            ),
            approximations=(
                "A linear law replaces a vertical structure with two numbers, so the "
                "vertical profile of heat uptake is absent rather than approximated.",
                "The exchange coefficient is calibrated against the behaviour of models "
                "and observations, not derived; EQR-004 requires that be carried with it.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("T_s", 200.0, 340.0, "K"), ("T_d", 200.0, 340.0, "K")),
                conditions=(
                    "global-mean heat uptake on timescales from the mixed-layer "
                    "adjustment time to a few deep-reservoir times",
                ),
            ),
            limiting_cases=(
                "equal reservoir temperatures give exactly zero exchange",
                "an infinite deep capacity fixes T_d and recovers a constant-temperature "
                "lower boundary",
                "a zero exchange coefficient recovers the single-reservoir slab exactly",
            ),
            conserved_quantities=("energy",),
            permitted_discretizations=("forward_euler",),
            correctness_tests=(
                "test_equal_temperatures_exchange_no_heat",
                "test_two_reservoir_timescales_match_the_analytic_modes",
                "test_disabling_the_exchange_recovers_the_single_reservoir_result",
            ),
            known_failure_modes=(
                "constant exchange efficiency, so the model cannot represent a change in "
                "ocean heat uptake efficiency under sustained warming",
                "forward Euler is unstable if the step approaches the fast mode's "
                "timescale, which is set by the mixed layer and not by this equation",
            ),
            source="Gregory 2000; Held et al. 2010 two-layer energy-balance formulation",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["co2_band_absorption"],
            version=V1,
            name="Carbon dioxide band absorption coefficient",
            category=EquationCategory.CONSTITUTIVE_RELATION,
            mathematical_form=("kappa(v) = kappa_bg + k0 (M_CO2/M_air) x exp(-|v - v0| / D)"),
            physical_basis=(
                "A vibration-rotation band whose line strengths fall away from the band "
                "centre, represented as a single exponential envelope in wavenumber. The "
                "coefficient is per kilogram of air, so it multiplies the layer mass the "
                "radiation module already carries."
            ),
            variables=(
                ("kappa", "mass absorption coefficient of air", "m^2 kg^-1"),
                ("v", "wavenumber", "m^-1"),
                ("x", "carbon dioxide mole fraction", "1"),
            ),
            parameters=(
                ("co2_band_centre_absorption", "parameterized", "set for a saturated core"),
                ("co2_wing_decay", "calibrated", "fitted to a published forcing per doubling"),
                ("infrared_background_absorption", "parameterized", "stands for water vapour"),
            ),
            assumptions=(
                "carbon dioxide uniformly mixed through the column",
                "only the 15 micron band; the 4.3 micron band is omitted",
                "a single exponential wing envelope in place of a line-by-line spectrum",
            ),
            approximations=(
                "The wing envelope is an envelope, not a spectrum: it reproduces how the "
                "opaque width grows with concentration and not the line structure inside "
                "it, so nothing about spectrally resolved outgoing radiance is in scope.",
                "The background coefficient is held fixed, so there is no water-vapour "
                "feedback and a change between two runs is attributable to CO2 alone.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("x", 5.0e-5, 2.0e-3, "1"),),
                conditions=(
                    "concentrations within roughly a factor of four of present day, where "
                    "the wings stay inside the resolved spectral range",
                ),
            ),
            limiting_cases=(
                "zero concentration leaves only the background coefficient",
                "a saturated core makes the forcing depend on the logarithm of "
                "concentration rather than on concentration",
            ),
            conserved_quantities=("energy",),
            correctness_tests=(
                "test_forcing_per_doubling_is_nearly_constant",
                "test_zero_concentration_leaves_only_the_background",
                "test_the_band_is_confined_to_the_infrared",
            ),
            known_failure_modes=(
                "a grid that does not resolve the wings gives linear growth or complete "
                "saturation instead of a logarithm",
                "extrapolating far outside the calibrated range moves the wings outside "
                "the resolved band, which weakens the logarithm for a numerical reason",
            ),
            source=(
                "Standard band-model treatment; the wing decay is calibrated against a "
                "published 3.7 W m^-2 forcing per doubling and achieves 3.698 at the "
                "background opacity and spectral resolution the configurations use, "
                "both of which are constants shared with the calibration"
            ),
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["volcanic_aerosol_optics"],
            version=V1,
            name="Stratospheric volcanic aerosol optical properties",
            category=EquationCategory.CONSTITUTIVE_RELATION,
            mathematical_form=(
                "tau_scat(lambda) = tau_550 (lambda / 550 nm)^-alpha ; tau_abs = 0 ; g = g_s"
            ),
            physical_basis=(
                "Sulfate droplets in the stratosphere scatter sunlight with a "
                "single-scattering albedo very near one. The record supplies the optical "
                "depth at a reference wavelength directly, so only the spectral shape and "
                "the vertical placement are constructed."
            ),
            variables=(
                ("tau_scat", "scattering optical depth", "1"),
                ("tau_550", "optical depth at 550 nm, from the record", "1"),
                ("lambda", "wavelength", "m"),
                ("g", "scattering asymmetry", "1"),
            ),
            parameters=(
                ("sulfate_angstrom_exponent", "parameterized", "typical of stratospheric sulfate"),
                ("sulfate_asymmetry", "parameterized", "typical of stratospheric sulfate"),
            ),
            assumptions=(
                "a pure scatterer in the shortwave; single-scattering albedo taken as one",
                "no longwave effect, so no stratospheric warming after an eruption",
                "globally uniform and spread evenly through the declared layers",
            ),
            approximations=(
                "A single Angstrom exponent replaces a size distribution, so the spectral "
                "shape is an envelope rather than a spectrum.",
                "Vertical placement is crude in a column whose top layer is not a "
                "stratosphere; it matters for the omitted longwave effect and much less "
                "for the shortwave reflection this represents.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("tau_550", 0.0, 0.5, "1"),),
                conditions=("global-mean shortwave reflection by stratospheric aerosol",),
            ),
            limiting_cases=(
                "zero optical depth scatters nothing and changes no flux",
                "the scattering is confined to the shortwave and to the declared layers",
            ),
            conserved_quantities=("energy",),
            correctness_tests=(
                "test_zero_optical_depth_scatters_nothing",
                "test_the_aerosol_is_confined_to_the_shortwave_and_the_stratosphere",
                "test_a_volcanic_layer_reflects_and_cools",
            ),
            known_failure_modes=(
                "representing the aerosol as partly absorbing would put energy into the "
                "layer that in reality leaves the planet",
                "an asymmetry taken from cloud droplets rather than sulfate would halve "
                "the reflected fraction for the same optical depth",
            ),
            source="Sato et al. 1993 record; standard Mie parameters for stratospheric sulfate",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["vapour_band_absorption"],
            version=V1,
            name="Water vapour band and continuum absorption",
            category=EquationCategory.CONSTITUTIVE_RELATION,
            mathematical_form=(
                "kappa(v) = kappa_cont + k0 [ exp(-max(v - v_rot, 0) / D) "
                "+ exp(-|v - v_vib| / D) ];  tau_i(v) = kappa(v) w_i"
            ),
            physical_basis=(
                "A pure rotational band that makes the far infrared opaque, a bending "
                "band at 6.3 microns, and a continuum across the window between them. "
                "The coefficient is per kilogram of vapour and multiplies the vapour "
                "path in each layer, because water vapour is not uniformly mixed and "
                "cannot enter as a well-mixed coefficient times the layer mass."
            ),
            variables=(
                ("kappa", "mass absorption coefficient of water vapour", "m^2 kg^-1"),
                ("v", "wavenumber", "m^-1"),
                ("w", "vapour path in a layer", "kg m^-2"),
                ("tau", "absorption optical depth of a layer", "1"),
            ),
            parameters=(
                ("vapour_band_absorption", "parameterized", "set for saturated cores"),
                (
                    "vapour_wing_decay",
                    "calibrated",
                    "fitted to the CERES all-sky outgoing longwave flux, about "
                    "240 W m^-2, which is what a column carrying the all-sky "
                    "albedo of 0.3 must emit in equilibrium",
                ),
                (
                    "vapour_window_continuum",
                    "calibrated",
                    "fitted to a published clear-sky window transmittance near 0.65 at 10 microns",
                ),
            ),
            assumptions=(
                "local thermodynamic equilibrium, so Kirchhoff's law gives the emission",
                "the vapour distribution is supplied by the water module, not assumed here",
                "one wing decay serves both bands",
            ),
            approximations=(
                "The continuum is flat in wavelength and linear in vapour amount. The "
                "real one strengthens towards longer wavelengths and is quadratic, "
                "because it is dominated by self-broadening -- so this understates how "
                "fast the window closes as the atmosphere moistens.",
                "Near-infrared vapour absorption is omitted entirely, so the shortwave "
                "heating of moist air is absent and only the longwave effect is present.",
                "Calibrated against the all-sky outgoing longwave while the column "
                "carries no cloud, so this opacity stands in for the cloud longwave "
                "effect as well as the vapour's own. A configuration that also ran "
                "the cloud module would double-count it.",
                "The envelope reproduces where the wings cross unit optical depth, not "
                "the line structure inside them; spectrally resolved radiance is out of "
                "scope.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("w", 0.0, 1.0e2, "kg m^-2"),),
                conditions=(
                    "column vapour amounts within roughly a factor of three of the "
                    "present-day global mean, where the band wings stay inside the "
                    "resolved spectral range",
                ),
            ),
            limiting_cases=(
                "a dry column leaves no absorption at all and the window is fully open",
                "saturated band cores make the far infrared and the 6.3 micron band "
                "opaque at any ordinary vapour amount, so the response comes from the "
                "window and the wings",
                "holding the vapour fixed removes the feedback and leaves the Planck "
                "response, which is what the null control does",
            ),
            conserved_quantities=("energy",),
            correctness_tests=(
                "test_a_dry_column_absorbs_nothing_and_the_window_is_open",
                "test_the_window_closes_as_the_column_moistens",
                "test_the_bands_are_saturated_and_confined_to_the_infrared",
                "test_the_outgoing_longwave_matches_the_flux_the_column_absorbs",
                "test_letting_the_vapour_follow_the_temperature_amplifies_the_response",
            ),
            known_failure_modes=(
                "a grid that does not resolve the window separately from the bands loses "
                "the continuum's effect and with it the feedback",
                "extrapolating to a much moister column moves the wings outside the "
                "resolved range and weakens the feedback for a numerical reason",
            ),
            source=(
                "Standard band-model treatment; the wing decay and the window continuum "
                "are calibrated against published clear-sky radiation observations "
                "(CERES outgoing longwave, window transmittance) and never against a "
                "temperature record"
            ),
        )
    )

    return reg
