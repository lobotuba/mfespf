"""KER-009 and KER-010: interchangeable spectral resolutions and implementations.

The charter's variable-spectral requirement is that no radiation subsystem may
assume a fixed spectrum.  The operational test of that is whether the same
physics can be run at two spectral resolutions and give the same answer when --
and only when -- the physics is spectrally flat.  A framework that agrees under a
grey albedo but cannot represent a wavelength-dependent one has not met the
requirement; it has merely hidden the spectrum somewhere else.
"""

from __future__ import annotations

import numpy as np
import pytest

from mfespf.contracts.module import ModuleStatus
from mfespf.core.constants import value_of
from mfespf.physics.solar import ConstantSpectrumSolarSource, banded_grid, broadband_grid
from mfespf.run.configurations import airless_radiating_body, banded_variant

pytestmark = pytest.mark.stage(2)

TSI = value_of("total_solar_irradiance_reference")
C_SMALL = 1.0e6


def _equilibrium(kernel, steps=12000):
    result = kernel.run(steps)
    assert result.may_support_a_claim
    return result.steps[-1].outputs["surface_temperature"].scalar()


@pytest.mark.req("SPC-002", "SOL-001")
def test_spectral_integral_matches_total_irradiance():
    """Whatever the grid, the spectrum integrates to the declared total."""
    for grid in (broadband_grid(), banded_grid(4), banded_grid(32)):
        solar = ConstantSpectrumSolarSource(grid=grid, status=ModuleStatus.PRESCRIBED)
        solar.initialize({}, {})
        integrated = float(np.sum(solar._spectral_irradiance * grid.widths))
        assert integrated == pytest.approx(TSI, rel=1e-12)


@pytest.mark.req("KER-009", "SPC-007")
def test_grey_albedo_gives_the_same_answer_at_two_spectral_resolutions():
    """With a spectrally flat albedo, band count must not change the result."""
    broadband = _equilibrium(
        airless_radiating_body(
            albedo=0.3, areal_heat_capacity=C_SMALL, master_dt=600.0, solar_step_multiple=1
        )
    )
    banded = _equilibrium(
        banded_variant(
            16, albedo=0.3, areal_heat_capacity=C_SMALL, master_dt=600.0, solar_step_multiple=1
        )
    )
    assert banded == pytest.approx(broadband, rel=1e-9)


@pytest.mark.req("KER-009", "SPC-009")
def test_a_wavelength_dependent_albedo_changes_the_result():
    """Equal integrated power with a different spectrum need not give the same response.

    SPC-009 forbids assuming otherwise. Here the same total irradiance meets a
    surface that reflects strongly at short wavelengths, and the equilibrium
    temperature differs from the grey case -- which is the physical content the
    broadband configuration is structurally unable to represent.
    """
    n_bands = 16
    grid = banded_grid(n_bands)
    # Reflective below 1 micron, absorbing above: a crude snow-like surface.
    albedo = np.where(grid.centers < 1.0e-6, 0.8, 0.1)

    spectral = _equilibrium(
        banded_variant(
            n_bands,
            albedo=albedo,
            areal_heat_capacity=C_SMALL,
            master_dt=600.0,
            solar_step_multiple=1,
        )
    )
    grey = _equilibrium(
        banded_variant(
            n_bands,
            albedo=float(np.mean(albedo)),
            areal_heat_capacity=C_SMALL,
            master_dt=600.0,
            solar_step_multiple=1,
        )
    )
    # The spectral case is colder: most solar energy arrives where this surface
    # reflects, so a band-mean albedo materially misstates the absorption.
    assert spectral < grey - 1.0


@pytest.mark.req("KER-010", "SYS-001")
def test_two_solar_implementations_satisfy_one_interface():
    """The replaceability demonstration: same contract, different fidelity."""
    broadband = ConstantSpectrumSolarSource(grid=broadband_grid())
    banded = ConstantSpectrumSolarSource(grid=banded_grid(16))

    for module in (broadband, banded):
        module.validate_configuration()
        assert module.contract.provided_interfaces[0].interface_id == "solar_source"
        assert module.contract.provided_interfaces[0].version == banded.contract.interface_version

    assert broadband.contract.fidelity_level != banded.contract.fidelity_level
    assert (
        broadband.contract.supported_spectral_representations
        != banded.contract.supported_spectral_representations
    )


@pytest.mark.req("SPC-001", "SPC-004")
def test_spectral_fields_travel_with_their_grid():
    """A receiver must never have to infer the grid from context."""
    kernel = banded_variant(8, areal_heat_capacity=C_SMALL, master_dt=600.0, solar_step_multiple=1)
    result = kernel.run(2)
    field = result.steps[-1].outputs["toa_spectral_irradiance"]
    grid = field.metadata["spectral_grid"]
    assert grid.n_bins == 8
    assert field.values.shape == (8,)
    assert field.spectral_coordinate == grid.grid_id
