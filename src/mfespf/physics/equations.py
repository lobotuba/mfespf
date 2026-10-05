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
    "relaxed_convective_adjustment": "EQ-RELAXED-CONVECTIVE-ADJUSTMENT",
    "solar_vapour_absorption": "EQ-SOLAR-VAPOUR-K-DISTRIBUTION",
    "observed_cloud_field": "EQ-ISCCP-CLOUD-FIELD",
    "wtg_circulation": "EQ-WTG-OVERTURNING-CIRCULATION",
    "area_weighted_mean": "EQ-AREA-WEIGHTED-GLOBAL-MEAN",
    "vapour_eddy_diffusion": "EQ-VAPOUR-EDDY-DIFFUSION",
    "ozone_solar_absorption": "EQ-OZONE-SOLAR-ABSORPTION",
    "measured_solar_spectrum": "EQ-MEASURED-SOLAR-SPECTRUM",
    "brewer_dobson_upwelling": "EQ-BREWER-DOBSON-UPWELLING",
    "lateral_heat_export": "EQ-LATERAL-HEAT-EXPORT",
    "band_insolation": "EQ-BAND-MEAN-INSOLATION",
    "co2_15_micron_band": "EQ-CO2-15-MICRON-K-DISTRIBUTION",
    "ozone_9_6_micron_band": "EQ-O3-9.6-MICRON-K-DISTRIBUTION",
    "fractional_top_convective_adjustment": "EQ-FRACTIONAL-TOP-CONVECTIVE-ADJUSTMENT",
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
                "gravity is whatever the configured GravityField says it is; the "
                "default ConstantGravity holds it fixed over the depth of the column",
            ),
            approximations=(
                "Constant gravity introduces an error of order 2 z / R_earth. This "
                "was recorded as 'about 1.5 percent at 100 km and negligible through "
                "the troposphere', which is true of a troposphere and wrong about "
                "this column: its top at 1000 Pa is 19 to 27 km up depending on the "
                "resolution. Measured against NewtonianAltitudeGravity, the constant "
                "puts 0.21 percent too little mass in the column and 0.6 to 0.85 "
                "percent too little in the topmost layer -- the layer that sets the "
                "outgoing longwave -- and the equilibrated column lands 0.037 K "
                "cooler with 0.085 W m^-2 more outgoing longwave, at 8 and at 16 "
                "layers alike.",
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
                ranges=(("tau", 0.0, float("inf"), "1"),),
                conditions=(
                    "absorbing, non-scattering medium",
                    # The ceiling was 100 until a column that reaches 1110 in the
                    # carbon dioxide band centre was measured against it. The
                    # record was already contradicting itself: it declares the
                    # limiting case "tau -> infinity gives zero transmittance"
                    # and then refused to be evaluated past 100. Exponential
                    # attenuation is exact at every non-negative optical depth,
                    # and beyond about 709 the transmittance underflows to
                    # exactly zero in double precision, which is the correct
                    # value and not a failure. Negative tau remains excluded,
                    # because it is a defect rather than a thick layer.
                    "no upper bound on tau: the relation is exact wherever it is "
                    "evaluated, and a real band centre is genuinely this opaque",
                ),
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
                "Two different errors, and they behave oppositely. The absolute error in "
                "the transmittance is worst near unit optical depth -- 0.032 at tau = "
                "1.28 -- and vanishes at both limits. The relative error in the "
                "absorptance is worst in the thin limit, where it tends to the fixed "
                "ratio D/2 = 0.83, a 17 percent underestimate however thin the layer is.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("tau", 0.0, float("inf"), "1"),),
                conditions=(
                    # The ceiling of 10 went the same way as the floor of 0.01,
                    # and for the same reason. Measured against the exact
                    # hemispheric integral, the absolute error in the
                    # transmittance is bounded by 0.0316 at *every* optical
                    # depth, peaks at tau = 1.28 in the middle of the old range,
                    # and falls to 7e-6 by tau = 10 and below double precision
                    # by tau = 30. The thick limit is where this approximation
                    # is at its best, not its worst: a saturated layer emits as
                    # a black body whatever slant path is assumed for it.
                    "no upper bound on tau: the absolute transmittance error is "
                    "bounded by 0.0316 everywhere and vanishes in the thick limit",
                    # The floor was 0.01 until it was measured. The approximation is
                    # an order of magnitude *more* accurate below that floor than
                    # inside the declared range -- 0.003 maximum transmittance error
                    # on [0, 0.01] against 0.032 on [0.01, 10] -- so the floor
                    # excluded the region where the relation is at its best. It came
                    # from reading the relative error as though it were the absolute
                    # one. The 17 percent thin-limit bias is real, but it is a bias
                    # on a vanishing quantity: it belongs to the uncertainty this
                    # relation carries, which the radiation module already declares
                    # as diffusivity_approximation, and not to the question of
                    # whether the relation may be evaluated at all.
                    "the accuracy is a property of the path, not of a layer: the "
                    "approximate transmittance is exactly multiplicative, so splitting a "
                    "path into layers leaves the approximate answer unchanged while "
                    "moving every per-layer optical depth, and sixteen layers of tau = "
                    "0.08 each -- comfortably inside this range -- carry the path optical "
                    "depth 1.28 at which the approximation is at its worst. A per-layer "
                    "range is necessary, not sufficient.",
                ),
            ),
            limiting_cases=("D approaches 2 as tau -> 0 and 1.66 near unit optical depth",),
            conserved_quantities=("energy",),
            correctness_tests=(
                "test_diffusivity_is_a_declared_parameterization",
                "test_the_diffusivity_approximation_is_most_accurate_in_the_thin_limit",
            ),
            known_failure_modes=(
                "underestimates the absorptance of an optically thin layer by 17 percent, "
                "and so its heating rate, however thin the layer is -- the flux error "
                "vanishes with tau but the relative heating-rate error does not",
                "biases heating rates in optically very thick layers",
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
            record_id=EQUATION_IDS["wtg_circulation"],
            version=V1,
            name="Relaxed spectral weak-temperature-gradient overturning between two columns",
            category=EquationCategory.REDUCED_PHYSICAL_EQUATION,
            mathematical_form=(
                "sinking column: f_k = c_p sigma (T_d - T_a)_k / (g tau_1 dS/dp), "
                "A_n = 2 sum_k f_k sin(n pi zeta_k) dzeta_k, "
                "F(zeta) = sum_(n=1..M) (A_n / n) sin(n pi zeta), "
                "zeta = (p_s - p) / (p_s - p_t), F = 0 above p_t; "
                "tau_1 = L / c_1, c_1 = N H / pi, N^2 = g (dS/dz) / (c_p T_bar), with "
                "dS/dp and dS/dz the rising column's bulk tropospheric values; "
                "rising column: -(1 - sigma) / sigma F; "
                "d(phi m)_k/dt from upwind vertical fluxes and lateral exchange "
                "(F_k - F_(k+1)) phi_upstream for phi in (s, q)"
            ),
            physical_basis=(
                "Horizontal temperature differences in the low-latitude free "
                "troposphere are removed by gravity waves in the time they take to "
                "cross between rising and sinking air (Raymond and Zeng 2005); each "
                "vertical mode travels at c_1 / n and so relaxes n times more slowly "
                "(Herman and Raymond 2014). Mass continuity closes the overturning; "
                "the air carries its dry static energy and water."
            ),
            variables=(
                ("F", "downward mass flux per unit area of the sinking column", "kg m^-2 s^-1"),
                ("T", "layer temperature of the sinking (d) and rising (a) column", "K"),
                ("S", "dry static energy, c_p T + g z", "J kg^-1"),
                ("sigma", "area fraction of the rising region", "1"),
                ("tau_1", "gravest mode's relaxation time", "s"),
            ),
            parameters=(
                (
                    "sigma",
                    "measured",
                    "0.476, the area-weighted fraction of the globe with upward "
                    "monthly-mean omega at 500 hPa, NCEP/NCAR Reanalysis 1, 1991-2020 "
                    "long-term means (NOAA PSL)",
                ),
                (
                    "L",
                    "measured",
                    "750 km, the area-weighted mean distance from descending "
                    "(omega500 > 0) to the nearest ascending grid point within 20-30 "
                    "degrees of the equator, same data; monthly 589 to 853 km",
                ),
                (
                    "p_t",
                    "parameterized",
                    "100 hPa, the tropopause of the GASS WTG intercomparison and SCAM6",
                ),
                ("M", "numerical", "modes kept: the layers below p_t, at most 20"),
            ),
            assumptions=(
                "two regions, one rising and one sinking",
                "each relaxes toward the area mean on the gravity-wave crossing time",
                "flux-form upwind transport",
            ),
            approximations=(
                "Midlatitude eddy transport is absent.",
                "The modes are sines in pressure, exact only for a uniform "
                "stratification, which is why the bulk value is used.",
            ),
            validity_domain=ValidityDomain(
                conditions=(
                    "a troposphere stable in bulk below the declared tropopause; the "
                    "module refuses to run otherwise",
                ),
            ),
            limiting_cases=(
                "equal column temperatures give no circulation and two independent columns",
                "a temperature difference in the gravest mode relaxes on tau_1",
                "energy and water leaving one column arrive in the other exactly",
                "equal and opposite mass fluxes weighted by area at every interface",
            ),
            conserved_quantities=("energy", "water", "mass"),
            permitted_discretizations=("flux-form upwind", "truncated sine series"),
            correctness_tests=(
                "test_the_circulation_conserves_energy_and_water_exactly",
                "test_the_gravest_mode_relaxes_on_tau_1",
                "test_each_mode_relaxes_n_times_more_slowly",
                "test_no_temperature_difference_means_no_circulation",
                "test_the_two_columns_carry_equal_and_opposite_mass",
            ),
            known_failure_modes=(
                "a single sinking region cannot represent the spread of subsidence rates",
                "the strict balance this replaced diverged near the cold top of convection "
                "(docs/protocols/two-column-circulation.md)",
            ),
            source=(
                "Herman and Raymond 2014, doi:10.1002/2014MS000359; Raymond and Zeng "
                "2005; Cohen et al. 2024, doi:10.1029/2023MS003866; Sobel and "
                "Bretherton 2000; Nilsson and Emanuel 1999"
            ),
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["area_weighted_mean"],
            version=V1,
            name="Area-weighted global mean of partitioning regions",
            category=EquationCategory.DIAGNOSTIC_IDENTITY,
            mathematical_form="X_global = sum_r a_r X_r, with sum_r a_r = 1",
            physical_basis=(
                "The global integral of a per-area quantity over regions that "
                "partition the globe is the area-weighted sum of their values."
            ),
            variables=(("a_r", "area fraction of region r", "1"),),
            assumptions=("the regions partition the globe",),
            approximations=(),
            validity_domain=ValidityDomain(conditions=("regions that partition the globe",)),
            limiting_cases=("one region of area one returns its own value",),
            conserved_quantities=("energy", "water"),
            permitted_discretizations=("exact",),
            correctness_tests=("test_the_global_mean_preserves_the_integral",),
            source="definition",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["observed_cloud_field"],
            version=V1,
            name="Observed global cloud field as disjoint sub-columns",
            category=EquationCategory.EMPIRICAL_RELATION,
            mathematical_form=(
                "F = (1 - sum_i a_i) F_clear + sum_i a_i F_i, where sub-column i holds "
                "one cloud in the layer containing its top pressure, with solar "
                "scattering optical depth tau_c, asymmetry 0.85, and thermal absorption "
                "optical depth tau_i (f_ice / 2.13 + (1 - f_ice) / 2.56)"
            ),
            physical_basis=(
                "Clouds reflect sunlight from their tops and emit at the temperature "
                "there. Satellite amounts are seen from above and so partition the "
                "sky; the independent-column approximation solves each part and "
                "weights it, which is exact for fluxes because the transfer is linear "
                "in its sources."
            ),
            variables=(
                ("a_i", "amount of cloud type i seen from above", "1"),
                (
                    "tau_c",
                    "radiatively weighted mean visible optical thickness of a cloud type",
                    "1",
                ),
                ("p_top", "cloud top pressure", "Pa"),
                ("f_ice", "ice share of the type", "1"),
            ),
            parameters=(
                (
                    "amounts, top pressures, optical thicknesses, ice shares",
                    "measured",
                    "Rossow and Schiffer (1999), BAMS 80, Table 4 (D2 global annual "
                    "means) and Table 5 (phase split); infrared ratios 2.56 (liquid, "
                    "Mie) and 2.13 (ice, Minnis et al. 1993a) from their section 3",
                ),
            ),
            assumptions=(
                "cloud types do not overlap and nothing lies beneath a visible cloud",
                "each type occupies one model layer",
                "no scattering of thermal radiation",
            ),
            approximations=(
                "Logarithmic-mean optical thicknesses stand for each type's spread; "
                "ISCCP chose them to reproduce area-averaged radiation.",
                "One asymmetry parameter for liquid and ice.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("tau_c", 0.0, 400.0, "1"),),
                conditions=("global-mean column; the field is imposed, not predicted",),
            ),
            limiting_cases=(
                "zero amounts recover the clear column exactly",
                "the amount-weighted fluxes conserve energy because each sub-column does",
            ),
            conserved_quantities=("energy",),
            permitted_discretizations=("independent sub-columns",),
            correctness_tests=(
                "test_the_cloud_types_partition_the_sky",
                "test_each_cloud_sits_in_the_layer_containing_its_top",
                "test_the_infrared_optical_thickness_follows_the_phase",
                "test_an_empty_field_is_the_clear_column",
            ),
            known_failure_modes=(
                "multi-layer cloud, which ISCCP cannot see and CloudSat does, is absent",
            ),
            source="Rossow and Schiffer 1999, BAMS 80, 2261-2287",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["solar_vapour_absorption"],
            version=V1,
            name="Absorption of sunlight by water vapour (k-distribution)",
            category=EquationCategory.CONSTITUTIVE_RELATION,
            mathematical_form=(
                "T(w) = sum_i h_i exp(-k_i w) / sum_i h_i per band, with the scaled path "
                "w = u (p / 300 hPa)^0.8 [1 + 0.00135 (T - 240 K)]; carried as one "
                "sub-band per term, sized to the term's share of the band's sunlight"
            ),
            physical_basis=(
                "Water vapour absorbs sunlight in near-infrared vibration-rotation "
                "bands at 0.94, 1.14, 1.38, 1.87 and 2.7 microns, and weakly in the "
                "visible. A k-distribution sorts each band by absorption strength; the "
                "one-parameter scaling folds the pressure broadening of line wings and "
                "the weak temperature dependence of line strength into the path."
            ),
            variables=(
                ("u", "water vapour path of a layer", "kg m^-2"),
                ("p", "layer pressure", "Pa"),
                ("T", "layer temperature", "K"),
                ("k_i", "absorption coefficient of term i at 300 hPa and 240 K", "m^2 kg^-1"),
                ("h_i", "fraction of the solar flux in term i", "1"),
            ),
            parameters=(
                (
                    "k_i, h_i, scaling exponent 0.8, temperature coefficient 0.00135",
                    "published",
                    "Chou and Suarez (1999), NASA/TM-1999-104606 vol. 15, Tables 2 and 3 "
                    "and equations 3.2-3.5; fitted by its authors to line-by-line "
                    "calculations and not adjusted here",
                ),
            ),
            assumptions=(
                "the k-distribution's ordering holds at every level (the correlation "
                "assumption, carried by the one-parameter scaling)",
                "within a source band the surface albedo does not vary with wavelength",
                "vapour absorbs and does not scatter sunlight",
            ),
            approximations=(
                "Sunlight is solved as diffuse radiation with the solver's diffusivity "
                "factor, not as a beam at a mean solar zenith angle.",
                "The 2.27-10 micron band's k-distribution is applied only up to the "
                "grid's first edge beyond 2.27 microns; the thermal bands above keep "
                "their own vapour representation.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("T", 150.0, 330.0, "K"),),
                conditions=("gas absorption only; no ozone, oxygen or carbon dioxide",),
            ),
            limiting_cases=(
                "no vapour transmits every solar band unchanged",
                "a uniform path transmits sum h exp(-k w) / sum h, the source's eq. 3.4",
                "at 300 hPa and 240 K the scaling factor is exactly one",
            ),
            conserved_quantities=("energy",),
            permitted_discretizations=("one sub-band per k-term",),
            correctness_tests=(
                "test_the_sub_bands_carry_the_source_weights",
                "test_a_uniform_path_transmits_the_source_k_distribution",
                "test_the_scaling_is_one_at_the_reference_state",
                "test_a_dry_column_on_the_rebuilt_grid_changes_nothing",
            ),
            known_failure_modes=(
                "with reflection placed at the surface, reflected sunlight crosses the "
                "absorbing vapour twice, lowering the planetary albedo",
            ),
            source="Chou and Suarez 1999, NASA/TM-1999-104606 vol. 15",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["ozone_solar_absorption"],
            version=V1,
            name="Absorption of sunlight by ozone (effective band coefficients)",
            category=EquationCategory.CONSTITUTIVE_RELATION,
            mathematical_form=(
                "tau_band = k_band x per layer, x the layer's ozone in cm-atm at STP; "
                "k_band from Chou and Suarez (1999) Table 3 in eight bands from 0.175 to "
                "0.7 microns, band 2 as two disjoint ranges"
            ),
            physical_basis=(
                "Ozone absorbs ultraviolet sunlight in the Hartley and Huggins bands and "
                "visible sunlight in the Chappuis band. The ultraviolet is taken up "
                "where the ozone is, in the stratosphere, so the air there is heated "
                "from above and warms upward."
            ),
            variables=(
                ("x", "ozone amount of a layer", "cm-atm"),
                ("k_band", "effective ozone absorption coefficient of a band", "cm-atm^-1"),
            ),
            parameters=(
                (
                    "k_band",
                    "published",
                    "Chou and Suarez (1999) Table 3, derived by its authors from WMO "
                    "(1986) cross sections over 127 sub-bands; band 8 includes the "
                    "near-infrared absorption by their eq. 3.14; not adjusted here",
                ),
                (
                    "ozone profile",
                    "measured",
                    "SWOOSH v2.73 combined, anomaly-filled, 1991-2020 global area-weighted "
                    "mean, 316 to 1 hPa; held at the 316 hPa value below; the column "
                    "between the model top and 1 hPa placed in the top layer",
                ),
            ),
            assumptions=(
                "one coefficient per band, with no temperature or pressure dependence, "
                "as in the source",
                "ozone absorbs and does not scatter",
                "the ozone amount is prescribed and does not respond to the column",
            ),
            approximations=(
                "Ozone's 9.6-micron longwave band is not represented, so its absorption "
                "of upwelling longwave in the lower stratosphere is absent.",
                "Heating that really happens above the model top is deposited in the top layer.",
                "Sunlight is solved as diffuse radiation with the solver's diffusivity "
                "factor, not as a beam at a mean solar zenith angle.",
            ),
            validity_domain=ValidityDomain(
                conditions=("sunlight between 0.175 and 0.7 microns",),
            ),
            limiting_cases=(
                "no ozone transmits every band unchanged",
                "the column total is the profile's integral, about 290 dobson units",
            ),
            conserved_quantities=("energy",),
            permitted_discretizations=("one grid band per source band",),
            correctness_tests=(
                "test_the_ozone_bands_carry_the_source_coefficients",
                "test_the_ozone_column_is_the_profile_integral",
                "test_ozone_absorbs_the_ultraviolet_and_little_else",
            ),
            known_failure_modes=(
                "without the longwave band the lower stratosphere is biased cold",
            ),
            source=(
                "Chou and Suarez 1999, NASA/TM-1999-104606 vol. 15; Davis et al. 2016, "
                "doi:10.5194/essd-8-461-2016"
            ),
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["measured_solar_spectrum"],
            version=V1,
            name="Measured solar spectral irradiance (TSIS-1 HSRS v2)",
            category=EquationCategory.EMPIRICAL_RELATION,
            mathematical_form=(
                "S_band = S_total x integral over the band of E(lambda) / integral of "
                "E(lambda), E the TSIS-1 HSRS v2 full-spectrum irradiance"
            ),
            physical_basis=(
                "The sun is not a black body: its ultraviolet is depleted by "
                "absorption in its own atmosphere, so a 5772 K black body puts half "
                "again too much sunlight below 0.4 microns."
            ),
            variables=(("E", "measured spectral irradiance at 1 au", "W m^-2 nm^-1"),),
            parameters=(
                (
                    "E",
                    "measured",
                    "Coddington et al. 2023, doi:10.1029/2022EA002637, 115.5 nm to 200 "
                    "microns; committed as its cumulative at 0.5 percent logarithmic "
                    "steps, within 0.25 percent of the full file on every band tested",
                ),
            ),
            assumptions=("the spectral shape is constant; only the total is configured",),
            approximations=(
                "The reduced table interpolates the cumulative linearly between nodes "
                "0.5 percent apart.",
            ),
            validity_domain=ValidityDomain(
                conditions=("wavelengths from 115.5 nm to 200 microns",),
            ),
            limiting_cases=(
                "the bands of a grid spanning the support sum to the total",
                "Chou and Suarez's band fractions are reproduced to within a few percent",
            ),
            conserved_quantities=("energy",),
            permitted_discretizations=("band integrals of the cumulative",),
            correctness_tests=(
                "test_the_measured_spectrum_integrates_to_its_total",
                "test_the_measured_spectrum_reproduces_chou_and_suarez_band_fractions",
            ),
            known_failure_modes=(
                "a fixed shape carries no solar-cycle variation in the ultraviolet",
            ),
            source="Coddington et al. 2023, Earth and Space Science 10",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["band_insolation"],
            version=V1,
            name="Annual-mean insolation over a latitude band",
            category=EquationCategory.REDUCED_PHYSICAL_EQUATION,
            mathematical_form=(
                "Q / S0 = < (h0 sin(phi) sin(delta) + cos(phi) cos(delta) sin(h0)) / pi >, "
                "averaged over solar longitude (sin(delta) = sin(epsilon) sin(lambda), "
                "circular orbit) and over the band by area; cos(h0) = -tan(phi) tan(delta)"
            ),
            physical_basis=(
                "A latitude band receives the daily-mean sunlight its latitudes see as the "
                "declination moves through the year; near the equator that is more than "
                "the spherical mean."
            ),
            variables=(
                ("phi", "latitude", "rad"),
                ("delta", "solar declination", "rad"),
                ("h0", "sunset hour angle", "rad"),
            ),
            parameters=(
                (
                    "epsilon",
                    "measured",
                    "mean obliquity at J2000, 84381.406 arcseconds (IAU 2006 precession)",
                ),
            ),
            assumptions=("circular orbit", "no atmosphere between the boundary and the band"),
            approximations=(
                "The band is treated as one column, so its latitudes' different sunlight is "
                "averaged before it is absorbed.",
            ),
            validity_domain=ValidityDomain(conditions=("annual-mean experiments",)),
            limiting_cases=("the whole globe (90 degrees) gives exactly one quarter",),
            conserved_quantities=("energy",),
            correctness_tests=("test_the_band_insolation_is_geometry",),
            source="Standard solar geometry; IAU 2006 obliquity",
        )
    )

    reg.register(
        EquationRecord(
            record_id="EQ-O3-9.6-MICRON-K-DISTRIBUTION",
            version=V1,
            name="Ozone's 9.6-micron band as a pressure-tabulated k-distribution",
            category=EquationCategory.CONSTITUTIVE_RELATION,
            mathematical_form=(
                "tau_n = k_n(p) u [1 + alpha (T - 250) + beta (T - 250)^2], six terms in "
                "980-1100 cm^-1 with weights g_n shared by every pressure; log k_n "
                "interpolated in log p between 0.1 and 1000 hPa, held above 0.1 hPa and "
                "extrapolated log-linearly below 1000 hPa; one sub-band per term, sized to "
                "its weight of the Planck flux at 250 K"
            ),
            physical_basis=(
                "Ozone's asymmetric-stretch band absorbs the infrared in the atmospheric "
                "window. Its lines narrow with pressure until Doppler broadening takes "
                "over, so its absorption per molecule depends on the pressure it is at."
            ),
            variables=(
                ("u", "layer ozone amount", "cm-atm"),
                ("p", "layer pressure", "Pa"),
                ("T", "layer temperature", "K"),
            ),
            parameters=(
                (
                    "g_n, k_n(p), alpha, beta",
                    "derived",
                    "least-squares exponential sum fitted to Chou, Suarez, Liang and Yan "
                    "2001 band-5 line-by-line transmission tables (section 4.2, Table 8), "
                    "as carried by WRF module_ra_goddard.F; RMS 0.0015 at 250 K, 0.0018 at "
                    "held-out pressures (ozone-infrared-band.md attempt 2)",
                ),
            ),
            assumptions=(
                "the k-distribution's ordering holds at every level",
                "ozone and the smooth vapour absorption within a sub-band multiply "
                "their transmittances",
            ),
            approximations=(
                "Chou's tables are band means; splitting the band by weight at 250 K is "
                "the k-distribution assumption.",
                "Fitted on homogeneous paths at the effective pressure and temperature "
                "of Chou's scaling.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("T", 170.0, 330.0, "K"),),
                conditions=("the temperature fit's range, as for Chou's tables",),
            ),
            limiting_cases=(
                "no ozone transmits every band unchanged",
                "above 0.1 hPa the absorption per molecule stops depending on pressure "
                "(the Doppler limit)",
            ),
            conserved_quantities=("energy",),
            permitted_discretizations=("one sub-band per k-term",),
            correctness_tests=(
                "test_the_band_reproduces_chou_s_tabulated_transmittance",
                "test_the_terms_fill_the_band_by_their_planck_weights",
                "test_no_ozone_absorbs_nothing_and_thin_air_holds_the_doppler_limit",
            ),
            known_failure_modes=(
                "a single power-law pressure scaling of the same band fails Chou's tables "
                "(attempt 1); the pressure table is what makes it fit",
            ),
            source="Chou, Suarez, Liang and Yan 2001, NASA/TM-2001-104606 vol. 19",
        )
    )

    reg.register(
        EquationRecord(
            record_id="EQ-CO2-15-MICRON-K-DISTRIBUTION",
            version=V1,
            name="The 15-micron CO2 band as a scaled k-distribution",
            category=EquationCategory.CONSTITUTIVE_RELATION,
            mathematical_form=(
                "tau_n = (k_n / mu) / 1.66 * u (p / p_r)^m [1 + alpha (T - 250) + beta "
                "(T - 250)^2], k_n = k_1 8^(n-1), six terms each for the centre (620-720 "
                "cm^-1; p_r 30 hPa, m 0.85) and the wings (540-620, 720-800 cm^-1; p_r 300 "
                "hPa, m 0.5); one sub-band per term, sized to its Planck-weighted fraction "
                "at 250 K"
            ),
            physical_basis=(
                "CO2's bending-mode band at 15 microns is made of lines whose widths "
                "shrink with pressure, so the band's absorption per molecule depends on "
                "the pressure it is at; the scaled amount carries that and the band's "
                "weak temperature dependence."
            ),
            variables=(
                ("u", "layer CO2 amount", "cm-atm"),
                ("p", "layer pressure", "Pa"),
                ("T", "layer temperature", "K"),
            ),
            parameters=(
                (
                    "k_1, fractions, p_r, m, alpha, beta",
                    "published",
                    "Chou, Suarez, Liang and Yan 2001, NASA/TM-2001-104606 vol. 19, "
                    "Tables 3 and 10, fitted by the authors to line-by-line calculations "
                    "with HITRAN 1996; not adjusted here",
                ),
            ),
            assumptions=(
                "the k-distribution's ordering holds at every level, carried by the scaling",
                "CO2 and the smooth vapour absorption within a sub-band multiply their "
                "transmittances",
            ),
            approximations=(
                "Chou et al. note the scaled k-distribution underestimates the cooling "
                "above 10 hPa, where Doppler broadening matters; their accurate method "
                "there is pre-computed transmission tables.",
                "Ozone's 9.6-micron band, methane, N2O, CFCs and the minor CO2 bands are "
                "not represented (about -4.9 W m^-2 of MLS outgoing longwave).",
            ),
            validity_domain=ValidityDomain(
                ranges=(("T", 170.0, 330.0, "K"),),
                conditions=("the temperature fit's range, as for Chou's tables",),
            ),
            limiting_cases=(
                "no CO2 transmits every band unchanged",
                "the centre's fractions sum to the Planck share of 620-720 cm^-1 at 250 K",
            ),
            conserved_quantities=("energy",),
            permitted_discretizations=("one sub-band per k-term",),
            correctness_tests=(
                "test_the_band_terms_carry_the_source_coefficients",
                "test_the_centre_fractions_are_the_planck_share",
                "test_thin_air_absorbs_less_per_molecule",
            ),
            known_failure_modes=(
                "above 10 hPa the scaled k-distribution underestimates the cooling",
            ),
            source="Chou, Suarez, Liang and Yan 2001, NASA/TM-2001-104606 vol. 19",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["brewer_dobson_upwelling"],
            version=V1,
            name="Brewer-Dobson upwelling through the tropical tropopause",
            category=EquationCategory.REDUCED_PHYSICAL_EQUATION,
            mathematical_form=(
                "d(phi m)_k/dt = F_k (phi_(k-1) - phi_k), phi in (c_p T + g z, q); F zero "
                "below 100 hPa, F0 from 100 to 70 hPa, linear in pressure to zero at the "
                "top; air enters and leaves laterally with the local properties"
            ),
            physical_basis=(
                "Wave drag in the stratosphere drives a mean meridional circulation whose "
                "tropical branch rises through the tropopause and whose air leaves "
                "poleward above it; rising air cools adiabatically."
            ),
            variables=(
                ("F", "upward mass flux per unit area of the column", "kg m^-2 s^-1"),
                ("phi", "dry static energy or specific humidity", "J kg^-1 or kg kg^-1"),
            ),
            parameters=(
                (
                    "F0",
                    "measured",
                    "5.9e9 kg s^-1 at 70 hPa between the turnaround latitudes, ERA-Interim "
                    "1989-2009 (Seviour et al. 2012), over the band's area",
                ),
                ("vertical shape", "parameterized", "declared in tropical-column.md"),
            ),
            assumptions=(
                "the upwelling is prescribed and does not respond to the column",
                "flux-form upwind transport",
            ),
            approximations=(
                "The vertical shape of the mass flux above 70 hPa is declared, not measured.",
            ),
            validity_domain=ValidityDomain(
                conditions=("a column standing for the band between the turnaround latitudes",),
            ),
            limiting_cases=(
                "no flux changes nothing",
                "in an isothermal layer rising air cools at g / c_p per metre",
                "mass is conserved layer by layer",
            ),
            conserved_quantities=("energy", "water", "mass"),
            permitted_discretizations=("flux-form upwind",),
            correctness_tests=(
                "test_the_upwelling_conserves_mass_layer_by_layer",
                "test_rising_air_cools_at_the_dry_adiabatic_rate_in_an_isothermal_layer",
                "test_no_exchange_leaves_the_column_unchanged",
            ),
            known_failure_modes=("a prescribed circulation cannot respond to warming",),
            source="Seviour, Butchart and Hardiman 2012, QJRMS; Butchart 2014, Rev. Geophys.",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["lateral_heat_export"],
            version=V1,
            name="A band's poleward heat export as a lateral boundary flux",
            category=EquationCategory.EMPIRICAL_RELATION,
            mathematical_form=(
                "H_k = -E m_k / sum(m below 100 hPa) for layers below 100 hPa, zero above; "
                "E the band's observed net TOA surplus"
            ),
            physical_basis=(
                "In equilibrium a band exports the radiation it absorbs in excess of what "
                "it emits; the atmosphere and ocean carry it poleward."
            ),
            variables=(("E", "heat export per unit area of the band", "W m^-2"),),
            parameters=(
                (
                    "E",
                    "measured",
                    "NCEP/NCAR R1 1991-2020 band-mean net TOA minus its global mean: 42.51, "
                    "38.464 and 32.95 W m^-2 for 20, 25 and 30 degrees",
                ),
            ),
            assumptions=(
                "the export is prescribed and does not respond to the column",
                "it is taken from the atmosphere below 100 hPa by layer mass",
            ),
            approximations=(
                "Part of the real export is oceanic and leaves from the surface.",
                "NCEP's global net TOA bias (-12.87 W m^-2) is removed as if uniform.",
            ),
            validity_domain=ValidityDomain(conditions=("annual-mean experiments",)),
            limiting_cases=("the export removed is exactly E, from below 100 hPa only",),
            conserved_quantities=("energy",),
            correctness_tests=(
                "test_no_exchange_leaves_the_column_unchanged",
                "test_the_exchange_budget_closes_through_the_lateral_boundary",
            ),
            source="Kalnay et al. 1996 (NCEP/NCAR R1); Trenberth and Caron 2001 for scale",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["fractional_top_convective_adjustment"],
            version=V1,
            name="Convective adjustment with a top inside a layer",
            category=EquationCategory.CLOSURE_OR_PARAMETERIZATION,
            mathematical_form=(
                "above a mixed group, s(m) = s_k + (s_{k+1} - s_k)(m - 1/2) over the "
                "layer's mass fraction m; the group mixes the fraction f below the "
                "level where s(f) equals its mixed value (N + f w mbar(f)) / "
                "(D + f w), and the remainder keeps its reconstructed mean"
            ),
            physical_basis=(
                "A convecting parcel rises to its level of neutral buoyancy, and that "
                "level lies wherever the environment's stability says, not at a model "
                "layer edge. Reconstructing the stability inside the layer above the "
                "group lets the mixing stop there, so the depth of convection varies "
                "continuously with the forcing instead of in whole layers."
            ),
            variables=(
                ("T", "layer temperature", "K"),
                (
                    "s",
                    "stability variable: potential temperature or the lapse-rate equivalent",
                    "K",
                ),
                ("f", "mass fraction of the layer above the group that is mixed", "1"),
            ),
            assumptions=(
                "the stability variable is linear in mass inside a layer, with the "
                "gradient it has to the layer above, and integrates to the layer mean",
                "the unmixed part of a partially penetrated layer keeps its "
                "reconstructed profile's mean",
                "vapour is uniform in specific humidity within a layer, so a mass "
                "fraction of the layer carries the same fraction of its vapour",
            ),
            approximations=(
                "The reconstruction is exact for a linear profile; at a layer interface "
                "consecutive reconstructions disagree by half the second difference of "
                "the stability, which vanishes as the grid refines.",
                "Onset is not continuous: the first mixing of a newly unstable pair "
                "already reaches part-way into the layer above. Convection that is "
                "active continuously, as it is in the equilibrium column, does not "
                "cross onset.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("T", 100.0, 400.0, "K"),),
                conditions=(
                    "hydrostatic columns without resolved vertical motion",
                    "a stable gradient above the group; where it is not, the layer "
                    "joins or not as a whole, as in the whole-layer scheme",
                ),
            ),
            limiting_cases=(
                "f = 0 everywhere is the whole-layer adjustment",
                "f = 1 is the whole layer joining, and the top then continues into "
                "the next layer continuously",
                "column enthalpy sum(C_i T_i) is conserved exactly",
                "the adjusted profile is stable",
            ),
            conserved_quantities=("energy", "water"),
            permitted_discretizations=("bisection for f; exact linear reconstruction",),
            correctness_tests=(
                "test_the_fractional_top_conserves_enthalpy_and_leaves_a_stable_profile",
                "test_the_fractional_top_moves_continuously_with_the_forcing",
                "test_without_partial_tops_the_adjustment_is_the_whole_layer_one",
                "test_the_fractional_top_carries_water_in_the_same_fraction",
            ),
            known_failure_modes=("onset of convection is still a jump of about half a layer",),
            source="this framework; the neutral-buoyancy top is standard (Emanuel 1994)",
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["relaxed_convective_adjustment"],
            version=V1,
            name="Convective adjustment relaxed over a timescale",
            category=EquationCategory.CLOSURE_OR_PARAMETERIZATION,
            mathematical_form=(
                "T_new = T* + (1 - exp(-dt / tau)) (T_adj - T*), where T* is the "
                "provisional profile and T_adj the instantaneous adjustment of "
                "EQ-DRY-CONVECTIVE-ADJUSTMENT applied to it; the same factor scales "
                "the surface's share and the vapour mixed over the same groups"
            ),
            physical_basis=(
                "Convection removes an instability over a finite time rather than "
                "within one coupling interval. Relaxing toward the neutral profile at "
                "a rate 1 / tau is the Betts-Miller form; integrated exactly over a "
                "step it applies the fraction 1 - exp(-dt / tau) of the adjustment, "
                "so the result does not depend on the step for a steady target."
            ),
            variables=(
                ("T", "layer temperature", "K"),
                ("tau", "convective relaxation timescale", "s"),
                ("dt", "coupling timestep", "s"),
            ),
            parameters=(
                (
                    "convective_relaxation_timescale",
                    "parameterized",
                    "declared over 1 h to 16 h, nominal 2 h (Frierson 2007, who "
                    "found an aquaplanet remarkably insensitive across that range; "
                    "Betts 1986). Chosen under the frozen rule in "
                    "docs/protocols/convective-relaxation-timescale.md, never "
                    "against a temperature record (INF-003)",
                ),
            ),
            assumptions=(
                "the instantaneous adjustment's target and groups are the right "
                "destination; only the rate of reaching it changes",
                "one timescale serves every layer and the surface",
            ),
            approximations=(
                "A single relaxation rate for all depths of convection, where real "
                "shallow and deep convection adjust on different timescales.",
                "The unstable groups are still whole layers, so the relaxation smooths "
                "the mixing in time and does not by itself remove the grid from where "
                "the top of the mixing lies.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("T", 100.0, 400.0, "K"),),
                conditions=("hydrostatic columns without resolved vertical motion",),
            ),
            limiting_cases=(
                "tau -> 0 recovers the instantaneous adjustment exactly",
                "tau -> infinity leaves the provisional profile unchanged",
                "column enthalpy sum(C_i T_i) is conserved exactly for every tau",
            ),
            conserved_quantities=("energy", "water"),
            permitted_discretizations=("exact exponential integration over the step",),
            correctness_tests=(
                "test_a_zero_timescale_is_the_instantaneous_adjustment_exactly",
                "test_the_relaxed_adjustment_applies_the_exact_exponential_fraction",
                "test_the_relaxed_adjustment_conserves_energy_and_water_every_step",
            ),
            known_failure_modes=(
                "if the timescale is long against the radiative destabilisation, the "
                "column carries a standing instability and the lapse rate departs "
                "from the target it was meant to reach",
            ),
            source="Betts 1986; Betts and Miller 1986; Frierson 2007",
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
                ranges=(("T", 150.0, 350.0, "K"), ("a", 0.0, 2.0, "1")),
                conditions=(
                    "saturated convecting columns above cloud base",
                    "below about 200 K the saturation humidity is negligible and the "
                    "relation degenerates to the dry adiabat, exactly and smoothly; "
                    "the floor is where that limit is verified rather than where the "
                    "convecting premise still holds",
                ),
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
                "test_the_entraining_lapse_rate_degenerates_to_dry_as_the_column_cools",
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
            mathematical_form=(
                "e_s(T) = e_triple (T/T_triple)^(dC/R_v) "
                "exp( (L_0 - dC T_triple)/R_v (1/T_triple - 1/T) ), "
                "with L(T) = L_0 + dC (T - T_triple), dC = c_pv - c_condensed, "
                "and the condensed phase ice below the triple point, liquid above"
            ),
            physical_basis=(
                "Equilibrium between water and its vapour along the phase boundary, "
                "integrated exactly with the latent heat varying as Kirchhoff's "
                "relation requires, over whichever condensed phase is stable."
            ),
            variables=(
                ("e_s", "saturation vapour pressure over the stable phase", "Pa"),
                ("T", "temperature", "K"),
            ),
            parameters=(
                ("water_latent_heat_vaporization_273", "measured", "thermodynamic tables"),
                ("water_latent_heat_fusion", "measured", "thermodynamic tables"),
                ("specific_heat_water_vapour", "measured", "thermodynamic tables"),
                ("specific_heat_liquid_water", "measured", "thermodynamic tables"),
                ("specific_heat_ice", "measured", "thermodynamic tables"),
                ("gas_constant_water_vapour", "derived", "molar mass of water"),
                (
                    "saturation_vapour_pressure_triple_point",
                    "measured",
                    "measured triple-point value",
                ),
            ),
            assumptions=(
                "the latent heat varies linearly with temperature, as Kirchhoff's "
                "relation gives for constant heat capacities",
                "the vapour behaves as an ideal gas",
                "the condensed volume is negligible beside the vapour volume",
                "below the triple point the stable condensed phase is ice, so the "
                "relevant latent heat is that of sublimation",
            ),
            approximations=(
                "The heat capacities are held constant, so the latent heat is linear "
                "in temperature rather than curved. Against Murphy and Koop (2005) "
                "the error is below 0.6 percent from 240 to 330 K, below 2 percent to "
                "200 K, 6.3 percent at 172 K and 13.6 percent by 150 K. The large "
                "relative errors are where the absolute amounts are radiatively "
                "negligible -- 6 percent of 1e-3 Pa -- but they are the reason the "
                "declared floor is where it is.",
                "At the hot end, against IAPWS SR1-86(1992) eq. 1, the error is -0.52 "
                "percent at 330 K, -0.71 at 340, -0.95 at 350 and -1.24 at 360. The "
                "ceiling was 340 K, just past Murphy and Koop's 332 K reach, and "
                "unmeasured above it; it is 350 K now, the last temperature at which "
                "the declared 1 percent holds. The measurement was prompted by a top "
                "layer reaching 342 K during spin-up (docs/protocols/model-top.md).",
                "Supercooled liquid is not represented: below the triple point this "
                "returns the ice value, so a mixed-phase cloud is out of scope.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("T", 170.0, 350.0, "K"),),
                conditions=(
                    "saturation over the stable phase: ice below the triple point, liquid above",
                    "accuracy better than 7 percent throughout and better than 1 "
                    "percent above 220 K, against Murphy and Koop (2005)",
                ),
            ),
            limiting_cases=(
                "returns the measured value exactly at the triple point, from both sides",
                "reduces to the constant-latent-heat form as the heat capacities converge",
            ),
            conserved_quantities=(),
            correctness_tests=(
                "test_saturation_vapour_pressure_at_the_triple_point",
                "test_saturation_vapour_pressure_matches_reference_values",
                "test_saturation_follows_clausius_clapeyron_with_a_varying_latent_heat",
                "test_saturation_over_ice_is_lower_than_over_supercooled_liquid",
                "test_saturation_is_continuous_across_the_triple_point",
            ),
            known_failure_modes=(
                "degrades below 180 K as the constant-heat-capacity assumption tires",
                "no supercooled liquid branch, so mixed-phase conditions are out of scope",
            ),
            source=(
                "Clausius-Clapeyron integrated with Kirchhoff's relation; validated "
                "against Murphy, D. M. and Koop, T. (2005), Review of the vapour "
                "pressures of ice and supercooled water for atmospheric applications, "
                "Q. J. R. Meteorol. Soc. 131, 1539-1565, and above 330 K against IAPWS "
                "SR1-86(1992), Revised Supplementary Release on Saturation Properties "
                "of Ordinary Water Substance"
            ),
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
            record_id=EQUATION_IDS["vapour_eddy_diffusion"],
            version=V1,
            name="Vertical eddy diffusion of water vapour",
            category=EquationCategory.CLOSURE_OR_PARAMETERIZATION,
            mathematical_form=(
                "F_i = -rho_i K (q_k - q_(k-1)) / dz_i upward across the interface between "
                "layers k-1 and k, rho_i = p_i / (R T_i), dz_i = (R T_i / g_i) "
                "ln(p_(k-1) / p_k), g = dp / m per layer; dw_k/dt = F_k - F_(k+1), "
                "F = 0 at the surface and the top"
            ),
            physical_basis=(
                "Mixing too small or slow for convection to represent carries trace "
                "gases across the tropopause and through the stratosphere; its net "
                "effect on a horizontal mean is down the mean gradient, and its rate "
                "is inferred from observed tracer profiles."
            ),
            variables=(
                ("F", "upward vapour flux at an interface", "kg m^-2 s^-1"),
                ("K", "eddy diffusivity", "m^2 s^-1"),
                ("q", "specific humidity", "kg kg^-1"),
            ),
            parameters=(
                (
                    "K",
                    "measured",
                    "0.39 m^2 s^-1, the minimum of Massie and Hunten's (1981) "
                    "tracer-derived profile (3.9e3 cm^2 s^-1 at 17 km, Hu et al. 2012 "
                    "Table 3); their profile rises to 11 m^2 s^-1 near 86 km",
                ),
            ),
            assumptions=(
                "unresolved mixing acts down the mean gradient",
                "one diffusivity at every level",
                "convective mixing is represented separately and is not in K",
            ),
            approximations=(
                "A horizontal-mean column has no mean ascent, so the stratospheric "
                "circulation's transport is folded into K along with small-scale mixing.",
                "The sensible enthalpy carried by the vapour is neglected.",
            ),
            validity_domain=ValidityDomain(
                conditions=(
                    "explicit steps with dt (D_below + D_above) / m <= 1 in every layer, "
                    "D = rho K / dz; the diffusion refuses a step that breaks it",
                ),
            ),
            limiting_cases=(
                "a uniform specific humidity gives no flux",
                "a linear profile on a uniform grid gives the analytic flux -rho K dq/dz",
                "the column total is unchanged",
            ),
            conserved_quantities=("water",),
            permitted_discretizations=("explicit flux form",),
            correctness_tests=(
                "test_eddy_diffusion_conserves_the_column",
                "test_a_uniform_humidity_does_not_diffuse",
                "test_the_flux_is_minus_rho_k_dq_dz",
                "test_an_unstable_step_is_refused",
            ),
            known_failure_modes=(
                "a global-mean column's tropopause is warmer than the tropical cold "
                "point that dries the stratosphere, so a 1-D cold trap can over-predict "
                "stratospheric water by orders of magnitude (Hu et al. 2012)",
            ),
            source=(
                "Massie and Hunten 1981, JGR 86, 9859; Hu, Seager and Bains 2012, "
                "doi:10.1088/0004-637X/761/2/166"
            ),
        )
    )

    reg.register(
        EquationRecord(
            record_id=EQUATION_IDS["latent_energy"],
            version=V1,
            name="Latent energy content of atmospheric vapour",
            category=EquationCategory.DIAGNOSTIC_IDENTITY,
            mathematical_form=(
                "E_latent = L_v * sum_i w_i, with the release per layer "
                "L_v + L_f where T_i < T_triple and L_v above it"
            ),
            physical_basis=(
                "Energy that was supplied to evaporate the vapour and is returned when it "
                "condenses. Treating it as a stored quantity is what lets the water and "
                "energy budgets close simultaneously between the two events. Below the "
                "triple point the vapour deposits as ice rather than condensing as "
                "liquid, so it returns the latent heat of sublimation; the fusion part of "
                "that is paid by the surface when the ice lands and melts, which is what "
                "makes the cycle close without a temperature-dependent store."
            ),
            variables=(
                ("E_latent", "latent energy per unit area", "J m^-2"),
                ("w_i", "vapour mass in a layer", "kg m^-2"),
                ("T_i", "layer temperature, which selects the condensed phase", "K"),
            ),
            parameters=(
                ("water_latent_heat_vaporization_273", "measured", "thermodynamic tables"),
                ("water_latent_heat_fusion", "measured", "thermodynamic tables"),
                ("triple_point_temperature", "measured", "thermodynamic tables"),
            ),
            assumptions=(
                "the latent heat is held at its triple-point value within each phase",
                "the condensed phase is ice below the triple point and liquid above, so "
                "the energy released follows the phase the saturation formula already "
                "assumed",
            ),
            approximations=(
                "The sensible heat carried by the water itself is neglected, both in the "
                "vapour and in the precipitation that returns to the surface. Holding the "
                "latent heat constant within a phase is what that neglect buys: a store "
                "that moved with temperature would change when a layer warmed with no "
                "water going anywhere, and no transfer would balance it.",
                "Against the Kirchhoff-varying latent heat the error is below 1 percent "
                "from 170 to 284 K, below 2 percent to 294 K, 2.6 percent at 300 K and "
                "6.7 percent at 340 K. The declared range used to be [230, 340] with no "
                "stated condition, and it was inverted: the error below the triple point "
                "is now the smallest anywhere, 0.3 percent at 230 K, while the hot end it "
                "admitted is the worst.",
                "The hot end moved from 340 to 350 K with the saturation relation's, "
                "whose ceiling it shares: there the error is 7.7 percent against the "
                "Kirchhoff-varying latent heat (8.0 against IAPWS SR1-86), so the "
                "declared accuracy below is restated to that rather than left claiming "
                "6.7 over a wider range.",
            ),
            validity_domain=ValidityDomain(
                ranges=(("T", 170.0, 350.0, "K"),),
                conditions=(
                    "the latent heat released follows the condensed phase, which is ice "
                    "below the triple point",
                    "accuracy better than 7.7 percent throughout and better than 2 "
                    "percent below 294 K, against the Kirchhoff-varying latent heat",
                    "the floor is the saturation relation's, because below it no vapour "
                    "is placed and there is nothing for this relation to act on",
                ),
            ),
            limiting_cases=(
                "zero in a dry column",
                "reduces to the constant-L_v form in a column that never reaches the triple point",
            ),
            conserved_quantities=("energy",),
            correctness_tests=(
                "test_moist_column_budgets_close_for_water_and_energy",
                "test_depositing_as_ice_releases_the_heat_of_sublimation",
            ),
            known_failure_modes=(
                "The fusion energy is settled at the surface on the same step the ice "
                "forms, so a snowpack that persists is not represented: the melting is "
                "immediate by construction, as the precipitation already is.",
            ),
            source="Definition, with the phase taken from the stable-phase saturation",
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
                ranges=(("tau", 0.0, float("inf"), "1"), ("omega", 0.0, 1.0, "1")),
                conditions=(
                    # The ceiling of 1000 was an unmeasured guess and a real
                    # column passes it -- 1110 in the carbon dioxide band
                    # centre. Measured instead: reflectance, transmittance and
                    # emissivity sum to one to the last bit at every optical
                    # depth from zero to 1e8, at every single-scattering albedo
                    # up to and including conservative scattering, and with and
                    # without forward asymmetry. The transmittance underflows
                    # smoothly to zero and the layer saturates at its correct
                    # limits. The solver carries an explicit large-tau branch;
                    # it was written for this and the ceiling denied it.
                    "no upper bound on tau: energy closes exactly and the layer "
                    "saturates correctly at every optical depth measured, to 1e8",
                ),
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
