"""Reconstructing a run's energy budget over its whole length.

Charter section 19 stage 6 and 8.5; requirements CON-009, CON-010, KER-013,
ACC-004.

The coupler already checks conservation at every synchronization point, and
those checks pass to round-off. This does something different and asks a
question they cannot: over the *entire* replay period, does the change in the
system's stored energy equal the net energy that crossed its boundary?

The distinction matters because a per-interval residual of 1e-16 accumulated
over fifty thousand intervals is still 1e-16 relative to the interval, and
whether it is negligible against the *run* depends on how the errors add. A
drift that biases every interval the same way grows linearly; one that is
random grows as the square root. Neither is visible in a per-interval report,
and the second stage of the ladder onward is where a slow drift would first be
mistaken for a trend.

Nothing here recomputes physics. Every number comes from records the modules
already published -- the coupling residuals for what crossed the boundary, the
module budgets for what is stored -- so a disagreement is a disagreement between
the run and itself, which is the only kind worth reporting.
"""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field as dc_field

from mfespf.core.errors import ConfigurationError
from mfespf.coupling.conservation import ConservationReport
from mfespf.coupling.transfers import ConservedQuantity

__all__ = ["EnergyBudgetReconstruction", "reconstruct_budget"]


@dataclass(frozen=True)
class EnergyBudgetReconstruction:
    """A conserved quantity's whole-run budget, assembled from the run's own records."""

    quantity: ConservedQuantity
    units: str
    interval: tuple[float, float]
    n_intervals: int
    system_modules: tuple[str, ...]
    boundary_modules: tuple[str, ...]
    entered: float
    left: float
    stored_change: float
    physical_sources: float = 0.0
    numerical_sources: float = 0.0
    source_labels: dict[str, float] = dc_field(default_factory=dict)

    @property
    def duration(self) -> float:
        return self.interval[1] - self.interval[0]

    @property
    def net_boundary_flow(self) -> float:
        return self.entered - self.left

    @property
    def residual(self) -> float:
        """What the system stored, less what crossed its boundary and was made."""
        return self.stored_change - (
            self.net_boundary_flow + self.physical_sources + self.numerical_sources
        )

    @property
    def scale(self) -> float:
        """The largest term the residual is judged against, never the net change.

        A run whose inflow and outflow nearly cancel has a small stored change
        and large flows; scaling by the change would make a residual that is a
        percent of the flows look like a rounding error.
        """
        return max(
            abs(self.entered),
            abs(self.left),
            abs(self.stored_change),
            abs(self.physical_sources),
            abs(self.numerical_sources),
            1e-300,
        )

    @property
    def relative_residual(self) -> float:
        return abs(self.residual) / self.scale

    @property
    def mean_imbalance(self) -> float:
        """The residual expressed as an average flux over the run, in W m^-2.

        The form a reader can compare against something physical. A whole-run
        residual of a few gigajoules per square metre means nothing until it is
        divided by the duration; a drift of 0.01 W m^-2 is immediately
        recognisable as comparable with the imbalances the model is used to
        study.
        """
        if self.duration <= 0.0 or self.quantity is not ConservedQuantity.ENERGY:
            return 0.0
        return self.residual / self.duration

    def passed(self, tolerance: float = 1e-10) -> bool:
        return self.relative_residual <= tolerance

    def describe(self) -> str:
        lines = [
            f"{self.quantity.value} reconstruction over [{self.interval[0]:g}, "
            f"{self.interval[1]:g}] s ({self.n_intervals} coupling intervals)",
            f"  system:     {', '.join(self.system_modules)}",
            f"  boundary:   {', '.join(self.boundary_modules) or 'none'}",
            f"  entered     {self.entered:+.8g} {self.units}",
            f"  left        {self.left:+.8g} {self.units}",
            f"  net across  {self.net_boundary_flow:+.8g} {self.units}",
            f"  stored      {self.stored_change:+.8g} {self.units}",
        ]
        if self.physical_sources or self.numerical_sources:
            lines.append(
                f"  sources     physical {self.physical_sources:+.8g}, "
                f"numerical {self.numerical_sources:+.8g}"
            )
        for label, amount in sorted(self.source_labels.items()):
            lines.append(f"    {label}: {amount:+.8g}")
        lines.append(f"  residual    {self.residual:+.8g} (relative {self.relative_residual:.3e})")
        if self.quantity is ConservedQuantity.ENERGY and self.duration > 0.0:
            lines.append(f"  as a mean imbalance: {self.mean_imbalance:+.4g} W m^-2")
        return "\n".join(lines)


def reconstruct_budget(
    reports: list[ConservationReport],
    *,
    quantity: ConservedQuantity = ConservedQuantity.ENERGY,
    boundary_modules: tuple[str, ...] = ("space",),
) -> EnergyBudgetReconstruction:
    """Sum a run's per-interval records into one budget for the whole run.

    ``boundary_modules`` names the stores that are outside the system. That is a
    modelling decision and not something to infer: the space boundary
    accumulates energy that has left, so counting it inside the system would
    make every run conserve trivially and check nothing. Requiring it to be
    named is the point.
    """
    if not reports:
        raise ConfigurationError(
            "no conservation reports to reconstruct from; run a whole number of coupling "
            "periods so at least one interval is reconciled (CON-009)"
        )
    boundary = set(boundary_modules)

    modules = {b.module_id for r in reports for b in r.budgets_for(quantity)}
    system = tuple(sorted(modules - boundary))
    if not system:
        raise ConfigurationError(
            f"every module holding {quantity.value} is named as a boundary; there is no system "
            f"left to reconstruct"
        )
    units = next(
        (b.units for r in reports for b in r.budgets_for(quantity) if b.module_id in system),
        quantity.si_units,
    )

    entered = left = stored_change = 0.0
    physical = numerical = 0.0
    labels: dict[str, float] = {}
    for report in reports:
        for budget in report.budgets_for(quantity):
            if budget.module_id in boundary:
                continue
            stored_change += budget.stored_change
            physical += budget.net_physical
            numerical += budget.net_numerical
            for label, amount in budget.source_labels.items():
                labels[label] = labels.get(label, 0.0) + amount
        for residual in report.coupling_residuals:
            if residual.quantity is not quantity:
                continue
            inside_sender = residual.sender in system
            inside_receiver = residual.receiver in system
            if inside_receiver and not inside_sender:
                # The receiver's own record, not the sender's: a mismatch between
                # the two is the coupler's business and is already reported there.
                entered += residual.imported
            elif inside_sender and not inside_receiver:
                left += residual.exported

    return EnergyBudgetReconstruction(
        quantity=quantity,
        units=units,
        interval=(reports[0].interval[0], reports[-1].interval[1]),
        n_intervals=len(reports),
        system_modules=system,
        boundary_modules=tuple(sorted(boundary)),
        entered=entered,
        left=left,
        stored_change=stored_change,
        physical_sources=physical,
        numerical_sources=numerical,
        source_labels=labels,
    )
