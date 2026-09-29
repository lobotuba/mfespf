"""The error report's process quantities count every channel they claim to."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from mfespf.evaluation.bias import surface_temperature_bias

pytestmark = pytest.mark.stage(3)


def _result(**series: float):
    """A run of ten identical steps publishing the given scalar fields."""
    base = {
        "surface_temperature": 290.0,
        "toa_incident_irradiance": 340.0,
        "reflected_shortwave_radiation": 100.0,
        "outgoing_longwave_radiation": 240.0,
        **series,
    }
    step = SimpleNamespace(
        outputs={k: SimpleNamespace(values=np.array(v)) for k, v in base.items()}
    )
    return SimpleNamespace(steps=[step] * 10)


@pytest.mark.req("VAL-005")
def test_the_bowen_ratio_counts_the_heat_convection_carries_off_the_surface():
    """Surface-based convection heats the air by mixing as well as by the bulk flux.

    Both are sensible heat. Counting only the bulk flux read the moist-convective
    column's Bowen ratio as 0.01 when it was 0.47, because the mixing had already
    carried the heat up and left the bulk flux little to do.
    """
    bulk_only = surface_temperature_bias(
        _result(surface_sensible_heat_flux=2.0, surface_latent_heat_flux=100.0), label="bulk"
    )
    both = surface_temperature_bias(
        _result(
            surface_sensible_heat_flux=2.0,
            surface_convective_heat_flux=48.0,
            surface_latent_heat_flux=100.0,
        ),
        label="both",
    )
    assert bulk_only.bowen_ratio == pytest.approx(0.02)
    assert both.bowen_ratio == pytest.approx(0.50)
