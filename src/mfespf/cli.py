"""Command-line entry point.

``mfespf run <configuration>`` executes a reference configuration and prints its
manifest, conservation report and cost estimate.  The output is the run record,
not a summary of it: whatever the terminal shows is what the manifest holds.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from mfespf import __version__
from mfespf.core.errors import MFESPFError
from mfespf.run import column_configurations, configurations

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
}


def _run(args: argparse.Namespace) -> int:
    builder = _CONFIGURATIONS[args.configuration]
    kwargs: dict[str, object] = {}
    if args.steps_per_solar_exchange is not None:
        kwargs["solar_step_multiple"] = args.steps_per_solar_exchange
    if args.eccentricity is not None:
        kwargs["eccentricity"] = args.eccentricity
    if args.heat_capacity is not None:
        kwargs["areal_heat_capacity"] = args.heat_capacity
    kernel = builder(**kwargs)  # type: ignore[arg-type]

    try:
        result = kernel.run(args.steps, assert_conservation=not args.allow_residual)
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
        olr = outputs.get("outgoing_longwave_radiation")
        reflected = outputs.get("reflected_shortwave_radiation")
        if olr is not None and reflected is not None:
            print(
                f"outgoing longwave: {olr.scalar():.4f} W m^-2   "
                f"reflected shortwave: {reflected.scalar():.4f} W m^-2"
            )

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
        "--allow-residual",
        action="store_true",
        help="continue past a conservation failure and report it rather than stopping",
    )
    run.set_defaults(func=_run)

    listing = sub.add_parser("configurations", help="list reference configurations")
    listing.set_defaults(func=lambda _a: (print("\n".join(sorted(_CONFIGURATIONS))), 0)[1])

    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
