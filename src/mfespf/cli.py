"""Command-line entry point.

``mfespf run <configuration>`` executes a reference configuration and prints its
manifest, conservation report and cost estimate.  The output is the run record,
not a summary of it: whatever the terminal shows is what the manifest holds.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from mfespf import __version__
from mfespf.core.errors import MFESPFError
from mfespf.evaluation.reconstruction import reconstruct_budget
from mfespf.evaluation.timescales import (
    ExponentialFit,
    TimescaleMeasurement,
    TimescaleReport,
    efolding_time,
    fit_exponential_modes,
)
from mfespf.experiments.forcing import (
    Forcing,
    ImpulseForcing,
    NoForcing,
    PeriodicForcing,
    PulseForcing,
    StepForcing,
)
from mfespf.run import column_configurations, configurations, historical

_CONFIGURATIONS = {
    # Stage 1: the airless radiating body.
    "airless": configurations.airless_radiating_body,
    "airless-banded": configurations.banded_variant,
    "null-control": configurations.null_control,
    # Stage 2: the radiating atmospheric column.
    "transparent": column_configurations.transparent_atmosphere,
    "absorbing": column_configurations.absorbing_emitting_atmosphere,
    "window": column_configurations.window_atmosphere,
    "grey-layer": column_configurations.single_grey_layer,
    # Stage 3: thermodynamics, convection and water.
    "radiative-convective": column_configurations.radiative_convective,
    "moist": column_configurations.moist_radiative_convective,
    # Stage 4: clouds.
    "no-cloud": column_configurations.no_cloud_control,
    "cloud": column_configurations.prescribed_cloud,
    # Stage 5: thermal storage and response timescales.
    "surface-storage": configurations.surface_thermal_storage,
    "ocean": configurations.reduced_ocean_heat_storage,
    "column-ocean": column_configurations.column_over_reduced_ocean,
    # Stage 6: historical replay of a declared record.
    "replay": historical.historical_replay,
    "replay-observed": historical.observed_replay,
}

# Configurations whose surface temperature is worth analysing for response
# timescales. Running the analysis on an unforced configuration would fit modes
# to a spin-up transient and label them response times, which is exactly the
# confusion TIM-006 is about.
_TIMESCALE_CONFIGURATIONS = frozenset({"surface-storage", "ocean", "column-ocean"})


def _build_forcing(args: argparse.Namespace) -> Forcing | None:
    """Assemble the declared forcing experiment (TIM-007)."""
    if args.forcing is None:
        return None
    amplitude = args.forcing_amplitude
    onset = args.forcing_onset
    if args.forcing == "none":
        return NoForcing()
    if args.forcing == "step":
        return StepForcing(amplitude=amplitude, onset=onset)
    if args.forcing == "pulse":
        return PulseForcing(amplitude=amplitude, onset=onset, duration=args.forcing_duration)
    if args.forcing == "impulse":
        # Energy rather than flux: an impulse one step wide has a flux that
        # depends on the step, and only the energy is reproducible.
        return ImpulseForcing(energy=amplitude * args.forcing_duration, onset=onset)
    return PeriodicForcing(amplitude=amplitude, period=args.forcing_period, onset=onset)


def _report_timescales(result: object, regime: str) -> None:
    """Measure and print the response timescales, per mechanism (TIM-008)."""
    steps = result.steps  # type: ignore[attr-defined]
    times = np.array([s.time for s in steps])
    values = np.array([s.outputs["surface_temperature"].scalar() for s in steps])
    if times.size < 32 or float(np.std(values)) == 0.0:
        return

    report = TimescaleReport()
    subsystem = "mixed_layer_and_deep_reservoir" if _has_ocean(steps) else "surface_slab"
    fits: dict[int, ExponentialFit] = {}
    for n_modes in (1, 2):
        try:
            fits[n_modes] = fit_exponential_modes(times, values, n_modes=n_modes)
        except MFESPFError:
            continue
    if not fits:
        return

    # Which model to believe is decided by the residuals, not by which
    # configuration was run: a second mode is reported only when admitting it
    # improves the fit by more than an order of magnitude.
    best = 1
    if 2 in fits and 1 in fits and fits[2].relative_residual * 10.0 < fits[1].relative_residual:
        best = 2
    fit = fits[best]
    mechanisms = (
        ("radiative_relaxation_of_the_mixed_layer", "deep_reservoir_heat_uptake")
        if best == 2
        else ("radiative_relaxation_of_the_mixed_layer",)
    )
    for amplitude, tau, mechanism, constrained in zip(
        fit.amplitudes, fit.timescales, mechanisms, fit.constrained, strict=True
    ):
        report.add(
            TimescaleMeasurement(
                observable="surface_temperature",
                subsystem=subsystem,
                mechanism=mechanism,
                forcing_regime=regime,
                timescale_s=tau,
                method=f"{best}-mode exponential fit to the {regime} response",
                amplitude=amplitude,
                relative_residual=fit.relative_residual,
                constrained=constrained,
            )
        )

    print()
    print(report.describe())
    for n_modes, candidate in sorted(fits.items()):
        print(f"  {n_modes}-mode relative residual: {candidate.relative_residual:.4g}")
    if any(not m.constrained for m in report.measurements):
        print(
            f"  the record spans {fit.record_span:.4g} s; a timescale longer than that "
            f"is an extrapolation under the fitted model, not a measurement. Run longer "
            f"before quoting it."
        )
    try:
        print(f"  e-folding time (no model assumed): {efolding_time(times, values):.6g} s")
    except MFESPFError as exc:
        print(f"  e-folding time unavailable: {exc}")


def _has_ocean(steps: list) -> bool:
    return bool(steps) and "ocean_deep_temperature" in steps[-1].outputs


def _run(args: argparse.Namespace) -> int:
    builder = _CONFIGURATIONS[args.configuration]
    kwargs: dict[str, object] = {}
    if args.steps_per_solar_exchange is not None:
        kwargs["solar_step_multiple"] = args.steps_per_solar_exchange
    if args.eccentricity is not None:
        kwargs["eccentricity"] = args.eccentricity
    if args.heat_capacity is not None:
        kwargs["areal_heat_capacity"] = args.heat_capacity
    if args.master_dt is not None:
        kwargs["master_dt"] = args.master_dt
    forcing = _build_forcing(args)
    if forcing is not None:
        kwargs["forcing"] = forcing
    try:
        kernel = builder(**kwargs)  # type: ignore[arg-type]
    except MFESPFError as exc:
        # Assembling a configuration can fail before anything runs: a replay
        # whose record has not been fetched is the common case, and the message
        # says what to do about it.
        print(f"cannot assemble {args.configuration!r}: {exc}", file=sys.stderr)
        return 2

    adapter = getattr(kernel, "adapter", None)
    steps = args.steps
    if adapter is not None:
        # The record's coverage, not a step count, is what bounds a replay: a run
        # asked for a boundary condition past the end of the record fails rather
        # than extrapolating, and saying so here is friendlier than the failure.
        available = historical.replay_steps(
            adapter, "total_solar_irradiance", kernel.schedule.master_dt
        )
        if steps > available:
            print(
                f"the record supports {available} master steps at "
                f"{kernel.schedule.master_dt:g} s; running {available} instead of {steps}"
            )
            steps = available

    try:
        result = kernel.run(steps, assert_conservation=not args.allow_residual)
    except MFESPFError as exc:
        print(f"run failed: {exc}", file=sys.stderr)
        return 2

    print(result.describe())
    if result.final_conservation is not None:
        print()
        print(result.final_conservation.describe())
    if result.steps:
        outputs = result.steps[-1].outputs
        surface = outputs.get("surface_temperature")
        if surface is not None:
            print()
            print(f"final surface temperature: {surface.scalar():.4f} K")
        layers = outputs.get("atmosphere_layer_temperature")
        if layers is not None:
            profile = ", ".join(f"{v:.2f}" for v in layers.values)
            print(f"layer temperatures (surface first): {profile} K")
        cre_sw = outputs.get("cloud_shortwave_radiative_effect")
        cre_lw = outputs.get("cloud_longwave_radiative_effect")
        if cre_sw is not None and cre_lw is not None:
            # SPC-010: the components first, the net only after them.
            print(
                f"cloud radiative effect  shortwave: {cre_sw.scalar():+.4f} W m^-2   "
                f"longwave: {cre_lw.scalar():+.4f} W m^-2   "
                f"net: {cre_sw.scalar() + cre_lw.scalar():+.4f} W m^-2"
            )
        sensible = outputs.get("surface_sensible_heat_flux")
        latent = outputs.get("surface_latent_heat_flux")
        if sensible is not None:
            line = f"sensible heat flux: {sensible.scalar():.4f} W m^-2"
            if latent is not None:
                line += f"   latent heat flux: {latent.scalar():.4f} W m^-2"
            print(line)
        water = outputs.get("atmosphere_precipitable_water")
        precipitation = outputs.get("precipitation_flux")
        if water is not None and precipitation is not None:
            # Reported in mm/day, the unit a hydrological claim would use.
            print(
                f"precipitable water: {water.scalar():.4f} kg m^-2   "
                f"precipitation: {precipitation.scalar() * 86400.0:.4f} mm/day"
            )
        deep = outputs.get("ocean_deep_temperature")
        uptake = outputs.get("ocean_heat_uptake")
        if deep is not None and uptake is not None:
            print(
                f"deep reservoir temperature: {float(deep.values):.4f} K   "
                f"ocean heat uptake: {float(uptake.values):+.4f} W m^-2"
            )
        anomaly = outputs.get("solar_irradiance_anomaly")
        if anomaly is not None:
            print(f"imposed irradiance anomaly at the final step: {anomaly.scalar():+.4f} W m^-2")
        olr = outputs.get("outgoing_longwave_radiation")
        reflected = outputs.get("reflected_shortwave_radiation")
        if olr is not None and reflected is not None:
            print(
                f"outgoing longwave: {olr.scalar():.4f} W m^-2   "
                f"reflected shortwave: {reflected.scalar():.4f} W m^-2"
            )

    if args.configuration in _TIMESCALE_CONFIGURATIONS and forcing is not None:
        _report_timescales(result, forcing.regime)

    if adapter is not None:
        print()
        print(adapter.describe())

    if result.conservation_reports and (adapter is not None or args.reconstruct):
        # The per-interval checks have already passed; this asks the different
        # question of whether the budget closes when summed over the whole run,
        # where a bias too small to see in one interval would have accumulated.
        print()
        print(reconstruct_budget(result.conservation_reports).describe())

    print()
    print(f"  {result.operating_mode.describe()}")

    if args.manifest:
        path = Path(args.manifest)
        result.manifest.write(path)
        print(f"\nmanifest written to {path}")

    if not result.may_support_a_claim:
        print(
            "\nThis run is not admissible as evidence for or against a claim.",
            file=sys.stderr,
        )
        return 1
    return 0


_EXPERIMENTS = {
    "solar-only": "mfespf.experiments.solar_only:run_solar_only_experiment",
    "solar-and-co2": "mfespf.experiments.solar_and_co2:run_solar_and_co2_experiment",
}


def _experiment(args: argparse.Namespace) -> int:
    """Run a declared experiment and print its record, verdict and limits."""
    import importlib

    module_name, function_name = _EXPERIMENTS[args.experiment].split(":")
    try:
        runner = getattr(importlib.import_module(module_name), function_name)
        result = runner()
    except MFESPFError as exc:
        print(f"cannot run {args.experiment!r}: {exc}", file=sys.stderr)
        return 2

    print(result.describe())
    print()
    print(result.outcome.describe())
    print()
    print(f"whole-run conservation residual: {result.worst_relative_residual:.3e}")

    # A refuted claim is a result, not an error: it is what the framework is for,
    # and exiting non-zero would make a successful refutation look like a broken
    # run to anything reading the exit code.
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="mfespf", description=__doc__)
    parser.add_argument("--version", action="version", version=f"mfespf {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="execute a reference configuration")
    run.add_argument("configuration", choices=sorted(_CONFIGURATIONS))
    run.add_argument("--steps", type=int, default=2000, help="number of master steps")
    run.add_argument("--manifest", help="path to write the run manifest to")
    run.add_argument(
        "--steps-per-solar-exchange",
        type=int,
        help="master steps between solar/orbit exchanges (multirate demonstration)",
    )
    run.add_argument("--eccentricity", type=float, help="orbital eccentricity")
    run.add_argument(
        "--heat-capacity",
        type=float,
        dest="heat_capacity",
        help="areal heat capacity of the surface slab in J m^-2 K^-1; lowering it "
        "shortens the response time and so the run needed to reach equilibrium",
    )
    run.add_argument(
        "--master-dt",
        type=float,
        dest="master_dt",
        help="master time step in seconds",
    )
    run.add_argument(
        "--forcing",
        choices=("none", "step", "pulse", "impulse", "periodic"),
        help="declared forcing experiment applied to total solar irradiance (TIM-007)",
    )
    run.add_argument(
        "--forcing-amplitude",
        type=float,
        default=20.0,
        help="forcing amplitude in W m^-2 of irradiance at 1 au; for an impulse this "
        "is multiplied by --forcing-duration to give the delivered energy",
    )
    run.add_argument(
        "--forcing-onset", type=float, default=0.0, help="forcing onset time in seconds"
    )
    run.add_argument(
        "--forcing-duration",
        type=float,
        default=8.64e6,
        help="pulse duration in seconds, or the impulse's energy divisor",
    )
    run.add_argument(
        "--forcing-period",
        type=float,
        default=3.456e7,
        help="period of a periodic forcing in seconds",
    )
    run.add_argument(
        "--reconstruct",
        action="store_true",
        help="reconstruct the energy budget over the whole run, not per interval",
    )
    run.add_argument(
        "--allow-residual",
        action="store_true",
        help="continue past a conservation failure and report it rather than stopping",
    )
    run.set_defaults(func=_run)

    experiment = sub.add_parser("experiment", help="run a declared experiment")
    experiment.add_argument("experiment", choices=sorted(_EXPERIMENTS))
    experiment.set_defaults(func=_experiment)

    listing = sub.add_parser("configurations", help="list reference configurations")
    listing.set_defaults(func=lambda _a: (print("\n".join(sorted(_CONFIGURATIONS))), 0)[1])

    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
