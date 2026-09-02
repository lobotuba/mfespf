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

    return reg
