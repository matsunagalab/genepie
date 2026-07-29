"""WHAM and MBAR free energy estimates from umbrella sampling data."""

import ctypes
import io
import os
import re
import tempfile
from typing import (
    Iterable,
    Optional,
)

from ..libgenesis import LibGenesis
from .. import ctrl_files
from .. import c2py_util
from ..exceptions import (
    ErrorCode,
    GenesisFortranNotSupportedError,
    GenesisValidationError,
)
from .._fortran import (
    ctrl_to_bytes,
    fortran_status,
)


def _check_free_energy_input(tool_name: str,
                            cvfile: Optional[str],
                            dcdfile: Optional[str]) -> None:
    """Reject the input modes the free energy wrappers cannot service.

    Unlike the CLI setup, wham_c_mod.fpp/mbar_c_mod.fpp never call
    define_molecules, so a DCD trajectory would be analysed against an empty
    molecule and silently produce garbage. Only precomputed CV files work.
    """
    if dcdfile:
        raise GenesisFortranNotSupportedError(
            f"{tool_name}: dcdfile input is not supported by the Python "
            "interface. Precompute the reaction coordinates (for example with "
            "trj_analysis) and pass them via cvfile.",
            code=ErrorCode.ERROR_NOT_SUPPORTED,
        )
    if not cvfile:
        raise GenesisValidationError(
            f"{tool_name}: cvfile is required. Pass a filename pattern whose "
            "placeholder expands to the replica index."
        )


# Matches the first ``{...}`` placeholder, the same span get_replicate_name1()
# replaces with the replica index in the Fortran analyze code.
_REPLICA_PLACEHOLDER = re.compile(r"\{[^}]*\}")


def _infer_n_replica(*sequences: Optional[Iterable]) -> Optional[int]:
    """Best-effort replica count from the per-replica restraint arrays.

    ``constant``/``reference`` are given as one sequence per dimension, each
    holding one value per replica, so the length of the first inner sequence is
    the replica count for a 1-D reaction coordinate. Returns ``None`` when it
    cannot be determined without consuming a lazy iterable.
    """
    for seq in sequences:
        if isinstance(seq, (list, tuple)) and seq:
            first = seq[0]
            if isinstance(first, (list, tuple)):
                return len(first)
    return None


def _validate_cvfiles_exist(tool_name: str, cvfile: str,
                            n_replica: Optional[int]) -> None:
    """Turn a missing CV file into a catchable exception before Fortran runs.

    ``open_file(..., IOFileInput)`` calls ``error_msg`` -> ``exit(1)`` when the
    file does not exist, which kills the host Python process instead of raising.
    We expand the replica placeholder (mirroring ``get_replicate_name1``) and
    check existence up front. When the replica count is unknown we still check
    the first name, which is the one the Fortran side opens first, so an obvious
    wrong path is always caught.
    """
    if _REPLICA_PLACEHOLDER.search(cvfile):
        count = n_replica if (n_replica and n_replica > 0) else 1
        names = [_REPLICA_PLACEHOLDER.sub(str(i), cvfile, count=1)
                 for i in range(1, count + 1)]
    else:
        names = [cvfile]

    missing = [name for name in names if not os.path.isfile(name)]
    if missing:
        listed = ", ".join(missing)
        raise GenesisValidationError(
            f"{tool_name}: cvfile not found: {listed}. Check the path and that "
            "the placeholder expands to the replica index."
        )


def wham_analysis(
        psffile: Optional[str] = None,
        prmtopfile: Optional[str] = None,
        ambcrdfile: Optional[str] = None,
        grotopfile: Optional[str] = None,
        grocrdfile: Optional[str] = None,
        pdbfile: Optional[str] = None,
        dcdfile: Optional[str] = None,
        cvfile: Optional[str] = None,
        check_only: Optional[bool] = None,
        allow_backup: Optional[bool] = None,
        dimension: Optional[int] = None,
        nblocks: Optional[int] = None,
        temperature: Optional[float] = None,
        tolerance: Optional[float] = None,
        rest_function: Optional[Iterable[int]] = None,
        grids: Optional[Iterable[tuple[float, float, int]]] = None,
        selection_group: Optional[Iterable[str]] = None,
        selection_mole_name: Optional[Iterable[str]] = None,
        function: Optional[Iterable[str]] = None,
        select_index: Iterable[Iterable[int]] = None,
        constant: Iterable[Iterable[float]] = None,
        reference: Iterable[Iterable[float]] = None,
        is_periodic: Iterable[bool] = None,
        box_size: Iterable[float] = None,
        ):
    """
    Computes a free energy profile with WHAM from umbrella sampling data.

    Only precomputed reaction coordinates (``cvfile``) are supported; see
    :func:`_check_free_energy_input` for why ``dcdfile`` is rejected.

    Args:
        cvfile: Filename pattern for the CV files, where the placeholder
            expands to the replica index (e.g. ``"run{}.dis"``).
        dimension: Number of reaction coordinates (1 or 2).
        nblocks: Number of blocks for block averaging.
        temperature: Simulation temperature in K.
        tolerance: Convergence threshold of the WHAM iteration.
        rest_function: Restraint function indices used as reaction coordinates.
        grids: ``(min, max, num_grids)`` per dimension.
        constant: Restraint force constants per replica.
        reference: Restraint centers per replica.
        is_periodic: Whether each reaction coordinate is periodic.
        box_size: Period of each periodic reaction coordinate.

    Returns:
        PMF as an array of shape ``(n_bins, n_columns)``, where the first
        column holds the bin centers.
    """
    # === INPUT VALIDATION ===
    from ..file_validators import validate_file_exists, validate_file_pattern
    from ..param_validators import validate_positive, validate_range

    _check_free_energy_input("wham_analysis", cvfile, dcdfile)

    # Validate input files
    validate_file_exists(psffile, "psffile", required=False)
    validate_file_exists(prmtopfile, "prmtopfile", required=False)
    validate_file_exists(grotopfile, "grotopfile", required=False)
    validate_file_pattern(cvfile, "cvfile", required=False)
    _validate_cvfiles_exist(
        "wham_analysis", cvfile, _infer_n_replica(reference, constant))

    # Validate parameters
    if dimension is not None:
        validate_range(dimension, 1, 2, "dimension")
    if temperature is not None:
        validate_positive(temperature, "temperature")
    if nblocks is not None:
        validate_positive(nblocks, "nblocks")
    # === END VALIDATION ===

    result_pmf_c = ctypes.c_void_p(None)
    n_bins = ctypes.c_int(0)
    n_bin_x = ctypes.c_int(0)
    try:
        ctrl = io.BytesIO()
        ctrl_files.write_ctrl_input(
                ctrl,
                psffile=psffile,
                prmtopfile=prmtopfile,
                ambcrdfile=ambcrdfile,
                grotopfile=grotopfile,
                grocrdfile=grocrdfile,
                pdbfile=pdbfile,
                dcdfile=dcdfile,
                cvfile=cvfile,
                )
        ctrl_files.write_ctrl_output(
                ctrl,
                pmffile="dummy.pmf")
        ctrl.write(b'[WHAM]\n')
        ctrl_files.write_kwargs(
                ctrl,
                check_only=check_only,
                allow_backup=allow_backup,
                dimension=dimension,
                nblocks=nblocks,
                temperature=temperature,
                tolerance=tolerance,
                rest_function=ctrl_files.NumberingData(rest_function),
                grids=ctrl_files.NumberingData(grids),
                )
        ctrl_files.write_ctrl_selection(
                ctrl, selection_group, selection_mole_name)
        ctrl.write(b'[RESTRAINTS]\n')
        ctrl_files.write_kwargs(
                ctrl,
                function=ctrl_files.NumberingData(function),
                select_index=ctrl_files.NumberingData(select_index),
                constant=ctrl_files.NumberingData(constant),
                reference=ctrl_files.NumberingData(reference),
                is_periodic=ctrl_files.NumberingData(is_periodic),
                box_size=ctrl_files.NumberingData(box_size),
                )

        ctrl_bytes, ctrl_len = ctrl_to_bytes(ctrl)
        with fortran_status() as (status, msgbuf, msglen):
            LibGenesis().lib.wa_analysis_c(
                    ctrl_bytes,
                    ctrl_len,
                    ctypes.byref(result_pmf_c),
                    ctypes.byref(n_bins),
                    ctypes.byref(n_bin_x),
                    ctypes.byref(status),
                    msgbuf,
                    ctypes.c_int(msglen),
                    )

        result_pmf = c2py_util.conv_double_ndarray(
                result_pmf_c, [n_bins.value, n_bin_x.value])

        return result_pmf
    finally:
        if result_pmf_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_pmf_c),
                    ctypes.byref(n_bins), ctypes.byref(n_bin_x))


def mbar_analysis(
        psffile: Optional[str] = None,
        prmtopfile: Optional[str] = None,
        ambcrdfile: Optional[str] = None,
        grotopfile: Optional[str] = None,
        grocrdfile: Optional[str] = None,
        pdbfile: Optional[str] = None,
        dcdfile: Optional[str] = None,
        cvfile: Optional[str] = None,
        check_only: Optional[bool] = None,
        allow_backup: Optional[bool] = None,
        nreplica: Optional[int] = None,
        input_type: Optional[str] = None,
        dimension: Optional[int] = None,
        temperature: Optional[float] = None,
        target_temperature: Optional[float] = None,
        nblocks: Optional[int] = None,
        tolerance: Optional[float] = None,
        rest_function: Optional[Iterable[int]] = None,
        grids: Optional[Iterable[tuple[float, float, int]]] = None,
        selection_group: Optional[Iterable[str]] = None,
        selection_mole_name: Optional[Iterable[str]] = None,
        constant: Iterable[Iterable[float]] = None,
        reference: Iterable[Iterable[float]] = None,
        is_periodic: Iterable[bool] = None,
        box_size: Iterable[float] = None,
        ):
    """
    Computes relative free energies with MBAR from multi-state sampling data.

    Only precomputed reaction coordinates (``cvfile``) are supported; see
    :func:`_check_free_energy_input` for why ``dcdfile`` is rejected.

    Args:
        cvfile: Filename pattern for the CV files, where the placeholder
            expands to the replica index (e.g. ``"run{}.dat"``).
        nreplica: Number of replicas (states).
        input_type: Sampling type, e.g. ``"US"`` for umbrella sampling.
        dimension: Number of reaction coordinates (1 or 2).
        temperature: Simulation temperature in K.
        target_temperature: Temperature at which the result is reported.
        nblocks: Number of blocks for block averaging.
        tolerance: Convergence threshold of the MBAR iteration.
        rest_function: Restraint function indices used as reaction coordinates.
        grids: ``(min, max, num_grids)`` per dimension.
        constant: Restraint force constants per replica.
        reference: Restraint centers per replica.
        is_periodic: Whether each reaction coordinate is periodic.
        box_size: Period of each periodic reaction coordinate.

    Returns:
        Relative free energies as an array of shape ``(n_replica, n_blocks)``.
    """
    # === INPUT VALIDATION ===
    from ..file_validators import validate_file_exists, validate_file_pattern
    from ..param_validators import validate_positive, validate_range

    _check_free_energy_input("mbar_analysis", cvfile, dcdfile)

    # Validate input files
    validate_file_exists(psffile, "psffile", required=False)
    validate_file_exists(prmtopfile, "prmtopfile", required=False)
    validate_file_exists(grotopfile, "grotopfile", required=False)
    validate_file_pattern(cvfile, "cvfile", required=False)
    _validate_cvfiles_exist("mbar_analysis", cvfile, nreplica)

    # Validate parameters
    if dimension is not None:
        validate_range(dimension, 1, 2, "dimension")
    if temperature is not None:
        validate_positive(temperature, "temperature")
    if nreplica is not None:
        validate_positive(nreplica, "nreplica")
    if nblocks is not None:
        validate_positive(nblocks, "nblocks")
    # === END VALIDATION ===

    result_fene_c = ctypes.c_void_p(None)
    n_replica = ctypes.c_int(0)
    n_blocks = ctypes.c_int(0)
    # The Fortran writer fills the returned array only when fenefile is set, so
    # it cannot be left empty. Point it at a scratch directory so that repeated
    # calls neither pollute the cwd nor collide with an existing fene.dat.
    scratch = tempfile.TemporaryDirectory()
    try:
        ctrl = io.BytesIO()
        ctrl_files.write_ctrl_input(
                ctrl,
                psffile=psffile,
                prmtopfile=prmtopfile,
                ambcrdfile=ambcrdfile,
                grotopfile=grotopfile,
                grocrdfile=grocrdfile,
                pdbfile=pdbfile,
                dcdfile=dcdfile,
                cvfile=cvfile,
                )
        ctrl_files.write_ctrl_output(
                ctrl,
                fenefile=os.path.join(scratch.name, "fene.dat"))
        ctrl.write(b'[MBAR]\n')
        ctrl_files.write_kwargs(
                ctrl,
                check_only=check_only,
                allow_backup=allow_backup,
                nreplica=nreplica,
                input_type=input_type,
                dimension=dimension,
                temperature=temperature,
                target_temperature=target_temperature,
                nblocks=nblocks,
                tolerance=tolerance,
                rest_function=ctrl_files.NumberingData(rest_function),
                grids=ctrl_files.NumberingData(grids),
                )
        ctrl_files.write_ctrl_selection(
                ctrl, selection_group, selection_mole_name)
        ctrl.write(b'[RESTRAINTS]\n')
        ctrl_files.write_kwargs(
                ctrl,
                constant=ctrl_files.NumberingData(constant),
                reference=ctrl_files.NumberingData(reference),
                is_periodic=ctrl_files.NumberingData(is_periodic),
                box_size=ctrl_files.NumberingData(box_size),
                )

        ctrl_bytes, ctrl_len = ctrl_to_bytes(ctrl)
        with fortran_status() as (status, msgbuf, msglen):
            LibGenesis().lib.mbar_analysis_c(
                    ctrl_bytes,
                    ctrl_len,
                    ctypes.byref(result_fene_c),
                    ctypes.byref(n_replica),
                    ctypes.byref(n_blocks),
                    ctypes.byref(status),
                    msgbuf,
                    ctypes.c_int(msglen),
                    )

        result_fene = c2py_util.conv_double_ndarray(
                result_fene_c, [n_replica.value, n_blocks.value])
        return result_fene
    finally:
        scratch.cleanup()
        if result_fene_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_fene_c),
                    ctypes.byref(n_replica), ctypes.byref(n_blocks))
