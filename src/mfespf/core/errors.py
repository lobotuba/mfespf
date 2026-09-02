"""Framework error taxonomy.

The charter distinguishes failure classes (section 15) and requires that a run
fail safely rather than continue past an incompatible interface (SYS-006).  The
exception hierarchy mirrors that taxonomy so a caught error already carries its
scientific classification.
"""

from __future__ import annotations

__all__ = [
    "ConfigurationError",
    "ConservationError",
    "CouplingError",
    "DimensionalityError",
    "InterfaceVersionError",
    "InvalidForScientificInferenceError",
    "MFESPFError",
    "OwnershipError",
    "ProvenanceError",
    "RegistryError",
    "ResourceBudgetError",
    "UnitError",
]


class MFESPFError(Exception):
    """Base class for every framework-raised error."""


class ConfigurationError(MFESPFError):
    """A configuration failed validation before initialisation (MOD-005)."""


class UnitError(MFESPFError):
    """A unit symbol is unknown or unusable on an interface."""


class DimensionalityError(UnitError):
    """Two quantities that must share a dimension do not (DAT-005)."""


class InterfaceVersionError(MFESPFError):
    """Interface version negotiation failed (VER-006, A6.5)."""


class OwnershipError(MFESPFError):
    """A state-ownership rule was violated (STA-001, STA-004, MOD-010)."""


class ConservationError(MFESPFError):
    """A conserved-quantity budget failed to reconcile (CON-004, CON-010)."""


class CouplingError(MFESPFError):
    """A coupling schedule or transfer pairing is invalid (CON-002, TIM-003)."""


class ProvenanceError(MFESPFError):
    """A dataset role or information-separation rule was violated (INF-003)."""


class RegistryError(MFESPFError):
    """A registry lookup or registration is invalid."""


class ResourceBudgetError(MFESPFError):
    """Estimated cost exceeds the declared resource budget (CMP-005)."""


class InvalidForScientificInferenceError(MFESPFError):
    """Raised when a run is asked to support a claim after failing a gate.

    A run may fail a mandatory gate and still be retained for debugging
    (AUD-007); this error marks the boundary it must not cross.
    """
