"""Spectral grid descriptor, negotiation and conservative regridding.

Charter section 10 and requirements SPC-001 through SPC-010.

The charter's central spectral rule is that no module may assume a fixed
spectrum (SPC-002).  Two consequences are implemented here.  First, a spectral
quantity always travels with the grid it is defined on, so a receiver can regrid
rather than assume.  Second, regridding conserves the integrated quantity to a
declared tolerance and logs its residual (SPC-005, SPC-006) -- an energy budget
that closes only because a regrid quietly lost 2 W m-2 is not a closed budget.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

import numpy as np

from mfespf.core.errors import ConfigurationError, DimensionalityError
from mfespf.core.units import parse_unit

__all__ = [
    "IntegrationConvention",
    "RegridResult",
    "SpectralCoordinate",
    "SpectralGrid",
    "negotiate_spectral_grid",
    "regrid_spectral_density",
]


class SpectralCoordinate(StrEnum):
    """SPC-003."""

    WAVELENGTH = "wavelength"
    FREQUENCY = "frequency"
    WAVENUMBER = "wavenumber"
    PHOTON_ENERGY = "photon_energy"
    BAND_ID = "band_id"


class IntegrationConvention(StrEnum):
    """How a stored value relates to the integral over its bin (SPC-004)."""

    BIN_INTEGRATED = "bin_integrated"  # value is the integral over the bin
    SPECTRAL_DENSITY = "spectral_density"  # value is per unit spectral coordinate


@dataclass(frozen=True)
class SpectralGrid:
    """A declared spectral discretisation (charter 10.2 / SPC-004).

    Bins are defined by their edges, so bin widths, adjacency and total coverage
    are all derivable rather than separately asserted.  ``resolution_identifier``
    is what a run manifest records and what SPC-008 compares against a request.
    """

    grid_id: str
    coordinate: SpectralCoordinate
    coordinate_units: str
    edges: np.ndarray
    integration_convention: IntegrationConvention = IntegrationConvention.SPECTRAL_DENSITY
    density_units: str | None = None
    interpolation_policy: str = "conservative_rebin"
    extrapolation_policy: str = "zero_outside_declared_support"
    resolution_identifier: str = ""

    def __post_init__(self) -> None:
        edges = np.asarray(self.edges, dtype=float)
        object.__setattr__(self, "edges", edges)
        if edges.ndim != 1 or edges.size < 2:
            raise ConfigurationError(f"spectral grid {self.grid_id!r} needs at least two bin edges")
        if not np.all(np.diff(edges) > 0):
            raise ConfigurationError(
                f"spectral grid {self.grid_id!r} edges must be strictly increasing; "
                f"declare ordering explicitly rather than relying on the consumer to sort"
            )
        if self.coordinate is not SpectralCoordinate.BAND_ID:
            parse_unit(self.coordinate_units)

    @property
    def n_bins(self) -> int:
        return int(self.edges.size - 1)

    @property
    def widths(self) -> np.ndarray:
        return np.diff(self.edges)

    @property
    def centers(self) -> np.ndarray:
        return 0.5 * (self.edges[:-1] + self.edges[1:])

    @property
    def support(self) -> tuple[float, float]:
        return (float(self.edges[0]), float(self.edges[-1]))

    def integrate(self, values: np.ndarray) -> float:
        """Integrate a spectral quantity over the whole grid.

        Interpretation follows ``integration_convention``: bin-integrated values
        are summed, densities are multiplied by bin width first.  Getting this
        wrong by a factor of the bin width is the classic spectral bug, so it is
        resolved once here rather than at every call site.
        """
        v = np.asarray(values, dtype=float)
        if v.shape[-1] != self.n_bins:
            raise ConfigurationError(
                f"values of shape {v.shape} do not match grid {self.grid_id!r} "
                f"with {self.n_bins} bins"
            )
        if self.integration_convention is IntegrationConvention.BIN_INTEGRATED:
            return float(np.sum(v, axis=-1))
        return float(np.sum(v * self.widths, axis=-1))

    def describe(self) -> str:
        lo, hi = self.support
        return (
            f"{self.grid_id} [{self.coordinate.value} {self.coordinate_units}] "
            f"{self.n_bins} bins over {lo:g}-{hi:g}, {self.integration_convention.value}"
        )


@dataclass(frozen=True, slots=True)
class RegridResult:
    """The regridded values together with the residual SPC-006 requires logging."""

    values: np.ndarray
    source_integral: float
    target_integral: float
    unmatched_source_integral: float
    method: str

    @property
    def absolute_residual(self) -> float:
        return abs(self.target_integral - self.source_integral)

    @property
    def relative_residual(self) -> float:
        if self.source_integral == 0.0:
            return 0.0
        return self.absolute_residual / abs(self.source_integral)

    def describe(self) -> str:
        return (
            f"regrid[{self.method}]: source={self.source_integral:.6g} "
            f"target={self.target_integral:.6g} residual={self.relative_residual:.3e}"
            + (
                f" (unmatched source support carried {self.unmatched_source_integral:.6g})"
                if self.unmatched_source_integral
                else ""
            )
        )


def regrid_spectral_density(
    values: np.ndarray,
    source: SpectralGrid,
    target: SpectralGrid,
    *,
    tolerance: float = 1e-12,
) -> RegridResult:
    """Rebin a spectral quantity conservatively (SPC-005).

    The method is exact overlap integration: each target bin receives the portion
    of each source bin's integral that falls inside it.  Within the overlap of
    the two supports the integral is conserved to round-off.  Source content
    outside the target's support is reported in ``unmatched_source_integral``
    rather than silently discarded -- that number is what tells a caller its
    target grid does not cover the physics it was handed.
    """
    if source.coordinate is not target.coordinate:
        raise DimensionalityError(
            f"cannot regrid from {source.coordinate.value} to {target.coordinate.value}; "
            f"convert the spectral coordinate at a declared boundary first"
        )
    if source.coordinate_units != target.coordinate_units:
        raise DimensionalityError(
            f"spectral grids use different coordinate units "
            f"({source.coordinate_units!r} vs {target.coordinate_units!r}); "
            f"convert at a declared boundary (DAT-004)"
        )

    v = np.asarray(values, dtype=float)
    if v.shape[-1] != source.n_bins:
        raise ConfigurationError(
            f"values of shape {v.shape} do not match source grid with {source.n_bins} bins"
        )

    # Work in bin-integrated space, where rebinning is a linear redistribution.
    if source.integration_convention is IntegrationConvention.SPECTRAL_DENSITY:
        source_integrals = v * source.widths
    else:
        source_integrals = v.copy()

    src_lo, src_hi = source.edges[:-1], source.edges[1:]
    tgt_lo, tgt_hi = target.edges[:-1], target.edges[1:]

    overlap = np.clip(
        np.minimum(src_hi[None, :], tgt_hi[:, None]) - np.maximum(src_lo[None, :], tgt_lo[:, None]),
        0.0,
        None,
    )
    fraction = overlap / source.widths[None, :]
    target_integrals = fraction @ source_integrals

    total_source = float(np.sum(source_integrals))
    matched = float(np.sum(fraction.sum(axis=0) * source_integrals))
    unmatched = total_source - matched

    if target.integration_convention is IntegrationConvention.SPECTRAL_DENSITY:
        out = target_integrals / target.widths
    else:
        out = target_integrals

    result = RegridResult(
        values=out,
        source_integral=total_source,
        target_integral=float(np.sum(target_integrals)),
        unmatched_source_integral=unmatched,
        method="exact_overlap_rebin",
    )

    # Conservation is checked over the shared support: content the target grid
    # does not cover is a coverage statement, not a conservation failure.
    covered_residual = abs(result.target_integral - matched)
    scale = abs(matched) or 1.0
    if covered_residual / scale > max(tolerance, 1e-12):
        raise DimensionalityError(
            f"spectral regrid lost {covered_residual:.6g} of {matched:.6g} within the shared "
            f"support (relative {covered_residual / scale:.3e} > tolerance {tolerance:.3e})"
        )
    return result


def negotiate_spectral_grid(
    *,
    provider_grids: list[SpectralGrid],
    requested_range: tuple[float, float] | None = None,
    requested_min_bins: int = 0,
    coordinate: SpectralCoordinate | None = None,
) -> tuple[SpectralGrid, tuple[str, ...]]:
    """Choose a provider grid for a consumer's request (SPC-008).

    Returns the selected grid and any *unmet* parts of the request.  An unmet
    request is returned rather than raised because the charter's rule is that the
    shortfall be recorded with its associated uncertainty, not that the run stop:
    a broadband module answering a request for line-by-line resolution is a
    legitimate low-fidelity configuration, provided the report says so.
    """
    if not provider_grids:
        raise ConfigurationError("no spectral grids offered by the provider")

    candidates = provider_grids
    if coordinate is not None:
        matching = [g for g in candidates if g.coordinate is coordinate]
        if matching:
            candidates = matching

    def covers(grid: SpectralGrid) -> bool:
        if requested_range is None:
            return True
        lo, hi = grid.support
        return lo <= requested_range[0] and hi >= requested_range[1]

    covering = [g for g in candidates if covers(g)]
    pool = covering or candidates
    # Prefer the coarsest grid that still meets the bin request: SPC/CMP both
    # push toward the cheapest configuration that answers the question.
    adequate = [g for g in pool if g.n_bins >= requested_min_bins]
    selected = (
        min(adequate, key=lambda g: g.n_bins) if adequate else max(pool, key=lambda g: g.n_bins)
    )

    unmet: list[str] = []
    if coordinate is not None and selected.coordinate is not coordinate:
        unmet.append(
            f"requested coordinate {coordinate.value}, received {selected.coordinate.value}"
        )
    if requested_range is not None and not covers(selected):
        lo, hi = selected.support
        unmet.append(
            f"requested support [{requested_range[0]:g}, {requested_range[1]:g}], "
            f"received [{lo:g}, {hi:g}]"
        )
    if selected.n_bins < requested_min_bins:
        unmet.append(f"requested >={requested_min_bins} bins, received {selected.n_bins}")
    return selected, tuple(unmet)
