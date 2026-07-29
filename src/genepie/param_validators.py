"""Parameter validation utilities for GENESIS Python interface."""

from typing import List, Optional, Union

from .constants import (
    FITTING_METHODS as FITTING_METHODS,
    PBC_CORRECTS as PBC_CORRECTS,
    TRJ_FORMATS as TRJ_FORMATS,
    TRJ_TYPES as TRJ_TYPES,
)
from .exceptions import GenesisValidationError


# Enum definitions matching GENESIS control file options
FORCEFIELDS = ["CHARMM", "AMBER", "GROAMBER", "GROMACS"]
ELECTROSTATICS = ["PME", "CUTOFF"]
ENSEMBLES = ["NVE", "NVT", "NPT", "NPAT", "NPgT"]
TPCONTROLS = ["NO", "LANGEVIN", "BUSSI", "BERENDSEN"]
INTEGRATORS = ["VVER", "LEAP", "VVER_CG"]
MINIMIZERS = ["SD", "LBFGS"]
BOUNDARY_TYPES = ["PBC", "NOBC"]
IMPLICIT_SOLVENTS = ["NONE", "GBSA", "EEF1", "IMM1"]
# FITTING_METHODS, TRJ_FORMATS, TRJ_TYPES and PBC_CORRECTS are re-exported from
# constants, which is the single source of truth shared with the Fortran layer.


def validate_enum(
    value: Optional[str],
    allowed: List[str],
    name: str,
    allow_none: bool = True
) -> None:
    """Validate that a string value is in the allowed list.

    Args:
        value: Value to validate
        allowed: List of allowed values (case-insensitive)
        name: Parameter name for error messages
        allow_none: If True, None is a valid value

    Raises:
        GenesisValidationError: If value is not in allowed list
    """
    if value is None:
        if not allow_none:
            raise GenesisValidationError(f"{name} is required")
        return

    upper_value = value.upper()
    upper_allowed = [a.upper() for a in allowed]

    if upper_value not in upper_allowed:
        raise GenesisValidationError(
            f"{name}: Invalid value '{value}'. "
            f"Allowed values: {', '.join(allowed)}"
        )


def validate_positive(
    value: Optional[Union[int, float]],
    name: str,
    allow_none: bool = True
) -> None:
    """Validate that a value is positive (> 0).

    Args:
        value: Value to validate
        name: Parameter name for error messages
        allow_none: If True, None is a valid value

    Raises:
        GenesisValidationError: If value is not positive
    """
    if value is None:
        if not allow_none:
            raise GenesisValidationError(f"{name} is required")
        return

    if value <= 0:
        raise GenesisValidationError(
            f"{name}: Must be positive (> 0), got {value}"
        )


def validate_non_negative(
    value: Optional[Union[int, float]],
    name: str,
    allow_none: bool = True
) -> None:
    """Validate that a value is non-negative (>= 0).

    Args:
        value: Value to validate
        name: Parameter name for error messages
        allow_none: If True, None is a valid value

    Raises:
        GenesisValidationError: If value is negative
    """
    if value is None:
        if not allow_none:
            raise GenesisValidationError(f"{name} is required")
        return

    if value < 0:
        raise GenesisValidationError(
            f"{name}: Must be non-negative (>= 0), got {value}"
        )


def validate_range(
    value: Optional[Union[int, float]],
    min_val: Union[int, float],
    max_val: Union[int, float],
    name: str,
    allow_none: bool = True
) -> None:
    """Validate that a value is within a range.

    Args:
        value: Value to validate
        min_val: Minimum allowed value (inclusive)
        max_val: Maximum allowed value (inclusive)
        name: Parameter name for error messages
        allow_none: If True, None is a valid value

    Raises:
        GenesisValidationError: If value is out of range
    """
    if value is None:
        if not allow_none:
            raise GenesisValidationError(f"{name} is required")
        return

    if value < min_val or value > max_val:
        raise GenesisValidationError(
            f"{name}: Must be in range [{min_val}, {max_val}], got {value}"
        )


def validate_distance_ordering(
    switchdist: Optional[float],
    cutoffdist: Optional[float],
    pairlistdist: Optional[float]
) -> None:
    """Validate switchdist <= cutoffdist <= pairlistdist.

    Args:
        switchdist: Switch distance for force smoothing
        cutoffdist: Cutoff distance for interactions
        pairlistdist: Pair list distance

    Raises:
        GenesisValidationError: If distance ordering is invalid

    Note:
        switchdist == cutoffdist is allowed (disables force switching).
    """
    if switchdist is not None and cutoffdist is not None:
        if switchdist > cutoffdist:
            raise GenesisValidationError(
                f"switchdist ({switchdist}) must be <= "
                f"cutoffdist ({cutoffdist})"
            )

    if cutoffdist is not None and pairlistdist is not None:
        if pairlistdist < cutoffdist:
            raise GenesisValidationError(
                f"pairlistdist ({pairlistdist}) must be >= "
                f"cutoffdist ({cutoffdist})"
            )


def validate_pme_params(
    electrostatic: Optional[str],
    pme_alpha: Optional[float] = None,
    pme_ngrid_x: Optional[int] = None,
    pme_ngrid_y: Optional[int] = None,
    pme_ngrid_z: Optional[int] = None,
    pme_nspline: Optional[int] = None
) -> None:
    """Validate PME parameters when electrostatic='PME'.

    Args:
        electrostatic: Electrostatic method
        pme_alpha: Ewald coefficient
        pme_ngrid_x/y/z: PME grid dimensions
        pme_nspline: B-spline order

    Raises:
        GenesisValidationError: If PME parameters are invalid
    """
    if electrostatic is None or electrostatic.upper() != "PME":
        return

    # PME alpha should be positive if specified
    if pme_alpha is not None:
        validate_positive(pme_alpha, "pme_alpha")

    # Grid dimensions should be positive if specified
    for name, val in [("pme_ngrid_x", pme_ngrid_x),
                      ("pme_ngrid_y", pme_ngrid_y),
                      ("pme_ngrid_z", pme_ngrid_z)]:
        if val is not None:
            validate_positive(val, name)

    # Spline order should be >= 4
    if pme_nspline is not None and pme_nspline < 4:
        raise GenesisValidationError(
            f"pme_nspline: Must be >= 4, got {pme_nspline}"
        )


def validate_shake_params(
    rigid_bond: Optional[bool],
    shake_iteration: Optional[int] = None,
    shake_tolerance: Optional[float] = None
) -> None:
    """Validate SHAKE parameters when rigid_bond=True.

    Args:
        rigid_bond: Whether SHAKE is enabled
        shake_iteration: Maximum SHAKE iterations
        shake_tolerance: SHAKE convergence tolerance

    Raises:
        GenesisValidationError: If SHAKE parameters are invalid
    """
    if not rigid_bond:
        return

    if shake_iteration is not None:
        validate_positive(shake_iteration, "shake_iteration")

    if shake_tolerance is not None:
        validate_positive(shake_tolerance, "shake_tolerance")


def validate_ensemble_params(
    ensemble: Optional[str],
    temperature: Optional[float] = None,
    pressure: Optional[float] = None,
    tpcontrol: Optional[str] = None
) -> None:
    """Validate ensemble-specific parameters.

    Args:
        ensemble: Ensemble type (NVE, NVT, NPT, etc.)
        temperature: Target temperature in K
        pressure: Target pressure
        tpcontrol: Temperature/pressure control method

    Raises:
        GenesisValidationError: If ensemble parameters are invalid
    """
    if ensemble is None:
        return

    validate_enum(ensemble, ENSEMBLES, "ensemble")
    ensemble_upper = ensemble.upper()

    # NVT/NPT require temperature
    if ensemble_upper in ["NVT", "NPT", "NPAT", "NPgT"]:
        if temperature is not None:
            validate_non_negative(temperature, "temperature")

    # NPT requires pressure
    if ensemble_upper in ["NPT", "NPAT", "NPgT"]:
        if pressure is not None:
            validate_positive(pressure, "pressure")

    # Validate tpcontrol if specified
    if tpcontrol is not None:
        validate_enum(tpcontrol, TPCONTROLS, "tpcontrol")


def validate_pbc_params(
    boundary_type: Optional[str],
    box_size_x: Optional[float] = None,
    box_size_y: Optional[float] = None,
    box_size_z: Optional[float] = None
) -> None:
    """Validate PBC parameters.

    Args:
        boundary_type: Boundary type (PBC, NOBC)
        box_size_x/y/z: Box dimensions

    Raises:
        GenesisValidationError: If PBC parameters are invalid
    """
    if boundary_type is None or boundary_type.upper() != "PBC":
        return

    # Box sizes should be positive if specified
    for name, val in [("box_size_x", box_size_x),
                      ("box_size_y", box_size_y),
                      ("box_size_z", box_size_z)]:
        if val is not None:
            validate_positive(val, name)


def validate_timestep(
    timestep: Optional[float],
    rigid_bond: Optional[bool] = None
) -> None:
    """Validate timestep is reasonable for MD simulation.

    Args:
        timestep: Time step in ps
        rigid_bond: Whether SHAKE/RATTLE is used

    Raises:
        GenesisValidationError: If timestep is invalid
    """
    if timestep is None:
        return

    validate_positive(timestep, "timestep")

    # Reject excessively large timesteps
    if timestep > 0.005:
        typical = "0.002 ps with SHAKE/RATTLE" if rigid_bond else "0.001 ps without SHAKE/RATTLE"
        raise GenesisValidationError(
            f"timestep ({timestep} ps) is very large. "
            f"Typical value is {typical}."
        )
