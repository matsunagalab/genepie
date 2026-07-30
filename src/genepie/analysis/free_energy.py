"""WHAM and MBAR free energy estimates from umbrella sampling data."""

import ctypes
import io
import os
import re
import tempfile
from collections import namedtuple
from typing import (
    Iterable,
    Optional,
    Union,
)

import numpy as np

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
from ._common import run_analysis_isolated


# Result bundle returned by ``mbar_analysis(..., return_weights=True)``.
#   fene:      relative free energies, shape (n_replica, n_blocks)
#   weights:   per-sample MBAR weights at the target temperature,
#              shape (n_replica, n_step); each row sums with the others to 1.
#   n_replica: number of states (replicas)
#   n_step:    number of samples per state
MbarResult = namedtuple(
    "MbarResult", ["fene", "weights", "n_replica", "n_step"])


# Result of a 1-D :func:`pmf_analysis`.
#   cv:           bin-center coordinates, shape (n_bins,)
#   pmf:          free energy from the standard histogram estimator (kcal/mol)
#   pmf_gaussian: free energy from the Gaussian-kernel estimator (kcal/mol)
# Both PMF curves have their minimum shifted to zero.
Pmf1DResult = namedtuple("Pmf1DResult", ["cv", "pmf", "pmf_gaussian"])


# Result of a 2-D :func:`pmf_analysis`.
#   cv1: bin-center coordinates along the first reaction coordinate, (n_x,)
#   cv2: bin-center coordinates along the second reaction coordinate, (n_y,)
#   pmf: Gaussian-kernel free energy matrix, shape (n_x, n_y); ``pmf[i, j]`` is
#        the free energy at ``(cv1[i], cv2[j])``, minimum shifted to zero.
Pmf2DResult = namedtuple("Pmf2DResult", ["cv1", "cv2", "pmf"])


# MBAR input_type keywords accepted by the Fortran reader (mbar_option_str.fpp).
MBAR_INPUT_TYPES = (
    "CV", "US", "ENESINGLE", "REMD", "ENEPAIR", "FEP", "ENEALL", "REST", "MBGO",
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
        temperature: Optional[Union[float, Iterable[float]]] = None,
        target_temperature: Optional[float] = None,
        nblocks: Optional[int] = None,
        tolerance: Optional[float] = None,
        self_iteration: Optional[int] = None,
        newton_iteration: Optional[int] = None,
        rest_function: Optional[Iterable[int]] = None,
        grids: Optional[Iterable[tuple[float, float, int]]] = None,
        selection_group: Optional[Iterable[str]] = None,
        selection_mole_name: Optional[Iterable[str]] = None,
        constant: Iterable[Iterable[float]] = None,
        reference: Iterable[Iterable[float]] = None,
        is_periodic: Iterable[bool] = None,
        box_size: Iterable[float] = None,
        return_weights: bool = False,
        ):
    """
    Computes relative free energies with MBAR from multi-state sampling data.

    Only precomputed reaction coordinates (``cvfile``) are supported; see
    :func:`_check_free_energy_input` for why ``dcdfile`` is rejected.

    Args:
        cvfile: Filename pattern for the CV files, where the placeholder
            expands to the replica index (e.g. ``"run{}.dat"``). For
            ``input_type="EneSingle"``/``"REMD"`` this is the per-state
            potential-energy time series (one ``time value`` pair per line).
        nreplica: Number of replicas (states).
        input_type: Sampling type. ``"US"``/``"CV"`` for umbrella sampling,
            ``"EneSingle"`` (equivalently ``"REMD"``) for temperature REMD.
        dimension: Number of reaction coordinates (1 or 2).
        temperature: Simulation temperature in K. May be a single value (used
            for every state) or one value per state, e.g. the T-REMD ladder.
        target_temperature: Temperature at which the result is reported.
        nblocks: Number of blocks for block averaging. Weight output requires
            ``nblocks == 1`` (the default).
        tolerance: Convergence threshold of the MBAR iteration.
        self_iteration: Number of self-consistent iterations.
        newton_iteration: Number of Newton-Raphson iterations.
        rest_function: Restraint function indices used as reaction coordinates.
        grids: ``(min, max, num_grids)`` per dimension.
        constant: Restraint force constants per replica.
        reference: Restraint centers per replica.
        is_periodic: Whether each reaction coordinate is periodic.
        box_size: Period of each periodic reaction coordinate.
        return_weights: When True, also compute the per-sample MBAR weights at
            ``target_temperature`` and return an :class:`MbarResult` instead of
            just the free-energy array. Requires ``nreplica`` and
            ``nblocks == 1``. Weights are returned directly from Fortran memory;
            no weight files are created.

    Returns:
        By default, relative free energies as an array of shape
        ``(n_replica, n_blocks)``. When ``return_weights=True``, an
        :class:`MbarResult` namedtuple ``(fene, weights, n_replica, n_step)``
        where ``weights`` has shape ``(n_replica, n_step)`` and the whole set of
        weights sums to 1 (they are the unbiased ensemble weights at the target
        temperature, ready for resampling).
    """
    # === INPUT VALIDATION ===
    from ..file_validators import validate_file_exists, validate_file_pattern
    from ..param_validators import validate_enum, validate_positive, validate_range

    _check_free_energy_input("mbar_analysis", cvfile, dcdfile)

    # Validate input files
    validate_file_exists(psffile, "psffile", required=False)
    validate_file_exists(prmtopfile, "prmtopfile", required=False)
    validate_file_exists(grotopfile, "grotopfile", required=False)
    validate_file_pattern(cvfile, "cvfile", required=False)
    _validate_cvfiles_exist("mbar_analysis", cvfile, nreplica)

    # Validate parameters
    validate_enum(input_type, MBAR_INPUT_TYPES, "input_type")
    if dimension is not None:
        validate_range(dimension, 1, 2, "dimension")
    if temperature is not None:
        if isinstance(temperature, (list, tuple)):
            for t in temperature:
                validate_positive(t, "temperature")
        else:
            validate_positive(temperature, "temperature")
    if nreplica is not None:
        validate_positive(nreplica, "nreplica")
    if nblocks is not None:
        validate_positive(nblocks, "nblocks")

    if return_weights:
        if nreplica is None:
            raise GenesisValidationError(
                "mbar_analysis: return_weights=True requires nreplica so the "
                "returned weight array dimensions can be validated."
            )
        if nblocks not in (None, 1):
            raise GenesisValidationError(
                "mbar_analysis: return_weights=True requires nblocks == 1 "
                f"(got nblocks={nblocks}). Per-sample weights are only defined "
                "for a single block."
            )
        if input_type is not None and input_type.upper() in ("ENEPAIR", "FEP"):
            raise GenesisFortranNotSupportedError(
                "mbar_analysis: per-sample target-ensemble weights are not "
                f"available for input_type={input_type}.",
                code=ErrorCode.ERROR_NOT_SUPPORTED,
            )
    # === END VALIDATION ===

    result_fene_c = ctypes.c_void_p(None)
    n_replica = ctypes.c_int(0)
    n_blocks = ctypes.c_int(0)
    result_weights_c = ctypes.c_void_p(None)
    n_weight_replica = ctypes.c_int(0)
    n_weight_step = ctypes.c_int(0)
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
                self_iteration=self_iteration,
                newton_iteration=newton_iteration,
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
                    ctypes.c_int(int(return_weights)),
                    ctypes.byref(result_fene_c),
                    ctypes.byref(n_replica),
                    ctypes.byref(n_blocks),
                    ctypes.byref(result_weights_c),
                    ctypes.byref(n_weight_replica),
                    ctypes.byref(n_weight_step),
                    ctypes.byref(status),
                    msgbuf,
                    ctypes.c_int(msglen),
                    )

        result_fene = c2py_util.conv_double_ndarray(
                result_fene_c, [n_replica.value, n_blocks.value])

        if not return_weights:
            return result_fene

        weights = c2py_util.conv_double_ndarray(
                result_weights_c,
                [n_weight_replica.value, n_weight_step.value])
        return MbarResult(
            fene=result_fene,
            weights=weights,
            n_replica=n_weight_replica.value,
            n_step=n_weight_step.value,
        )
    finally:
        scratch.cleanup()
        if result_fene_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_fene_c),
                    ctypes.byref(n_replica), ctypes.byref(n_blocks))
        if result_weights_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_weights_c),
                    ctypes.byref(n_weight_replica),
                    ctypes.byref(n_weight_step))


def _bin_centers(grid):
    """Reproduce the bin centers pmf_analysis places at each grid cell.

    Mirrors ``setup_option`` in ``pm_option.fpp``: ``num_grids`` grid edges give
    ``num_grids - 1`` bins whose centers are offset by half a bin width.
    """
    grid_min, grid_max, num_grids = grid
    num_grids = int(num_grids)
    delta = (grid_max - grid_min) / (num_grids - 1)
    return grid_min + 0.5 * delta + np.arange(num_grids - 1) * delta


def _as_per_dim(value, name, dimension):
    """Normalize a scalar or per-dimension sequence into a length-``dimension`` list."""
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        seq = list(value)
    else:
        seq = [value]
    if len(seq) == 1 and dimension > 1:
        seq = seq * dimension
    if len(seq) != dimension:
        raise GenesisValidationError(
            f"pmf_analysis: {name} must have one entry per dimension "
            f"(expected {dimension}, got {len(seq)})."
        )
    return seq


def _write_cv_weight_files(directory, cv, weight):
    """Serialize in-memory CV (and weight) samples to the pmf_analysis format.

    Returns ``(cvfile_pattern, weightfile_pattern)`` using a single ``{}``
    replica so the wrapper can drive the same file-based Fortran routine used
    for the CLI. The first column is a sample index (the Fortran reader treats
    it as the time stamp and ignores it).
    """
    cv = np.asarray(cv, dtype=np.float64)
    if cv.ndim == 1:
        cv = cv[:, None]
    n_sample, ndim = cv.shape

    cvfile = os.path.join(directory, "{}.cv")
    idx = np.arange(1, n_sample + 1)
    np.savetxt(cvfile.format(1), np.column_stack([idx, cv]))

    weightfile = None
    if weight is not None:
        weight = np.asarray(weight, dtype=np.float64).ravel()
        if weight.shape[0] != n_sample:
            raise GenesisValidationError(
                "pmf_analysis: weight length "
                f"({weight.shape[0]}) does not match the number of CV samples "
                f"({n_sample})."
            )
        weightfile = os.path.join(directory, "{}.weight")
        np.savetxt(weightfile.format(1), np.column_stack([idx, weight]))

    return cvfile, weightfile, ndim


def pmf_analysis(
        cv: Optional[Iterable] = None,
        weight: Optional[Iterable] = None,
        cvfile: Optional[str] = None,
        weightfile: Optional[str] = None,
        distfile: Optional[str] = None,
        grids: Optional[Iterable[tuple[float, float, int]]] = None,
        band_width: Optional[Iterable[float]] = None,
        dimension: Optional[int] = None,
        nreplica: Optional[int] = None,
        temperature: float = 300.0,
        cutoff: Optional[float] = None,
        is_periodic: Optional[Iterable[bool]] = None,
        box_size: Optional[Iterable[float]] = None,
        check_only: Optional[bool] = None,
        allow_backup: Optional[bool] = None,
        ):
    """
    Estimate a potential of mean force (PMF) from reaction-coordinate samples.

    This wraps the GENESIS ``pmf_analysis`` tool, which builds the PMF directly
    from collective-variable samples and optional per-sample weights (for
    example MBAR weights) using a histogram and a Gaussian-kernel estimator.

    There are two ways to provide the data:

    * **In-memory arrays** (convenient in notebooks): pass ``cv`` (and, for a
      reweighted PMF, ``weight``). ``cv`` is ``(n_sample,)`` for a 1-D PMF or
      ``(n_sample, 2)`` for a 2-D PMF. The arrays are written to temporary
      files internally.
    * **File patterns** (same as the CLI): pass ``cvfile`` (and optionally
      ``weightfile`` / ``distfile``) as filename patterns whose ``{}``
      placeholder expands to the replica index, e.g. ``"run{}.pathcv"``.

    Args:
        cv: In-memory CV samples, shape ``(n_sample,)`` or ``(n_sample, ndim)``.
        weight: In-memory per-sample weights, shape ``(n_sample,)``. When
            omitted every sample is weighted equally.
        cvfile: CV filename pattern (alternative to ``cv``).
        weightfile: Weight filename pattern (alternative to ``weight``).
        distfile: Optional path-CV distance filename pattern; samples whose
            distance is >= ``cutoff`` are discarded.
        grids: ``(min, max, num_grids)`` per dimension. ``num_grids`` grid
            edges produce ``num_grids - 1`` bins.
        band_width: Gaussian-kernel sigma per dimension.
        dimension: Number of reaction coordinates (1 or 2). Inferred from
            ``cv``/``grids`` when omitted.
        nreplica: Number of replicas when using ``cvfile`` patterns (default 1).
        temperature: Temperature in K used for ``-kT log P``.
        cutoff: Distance cutoff paired with ``distfile`` (0 disables filtering).
        is_periodic: Whether each reaction coordinate is periodic.
        box_size: Period of each periodic reaction coordinate.

    Returns:
        :class:`Pmf1DResult` for a 1-D PMF or :class:`Pmf2DResult` for a 2-D
        PMF (see their docstrings for the exact fields).
    """
    if cv is not None and cvfile is not None:
        raise GenesisValidationError(
            "pmf_analysis: pass either cv (in-memory) or cvfile (pattern), "
            "not both."
        )
    if cv is None and cvfile is None:
        raise GenesisValidationError(
            "pmf_analysis: provide the reaction coordinate via cv (array) or "
            "cvfile (filename pattern)."
        )
    if grids is None:
        raise GenesisValidationError("pmf_analysis: grids is required.")
    if band_width is None:
        raise GenesisValidationError("pmf_analysis: band_width is required.")

    # Normalize grids to a list of (min, max, num_grids) tuples so both the
    # inferred dimension and the returned bin centers are unambiguous.
    grids = list(grids)
    if grids and not isinstance(grids[0], (list, tuple)):
        grids = [tuple(grids)]
    grids = [tuple(g) for g in grids]

    scratch = tempfile.TemporaryDirectory()
    result_pmf_c = ctypes.c_void_p(None)
    n_out1 = ctypes.c_int(0)
    n_out2 = ctypes.c_int(0)
    try:
        if cv is not None:
            cvfile, weightfile, ndim = _write_cv_weight_files(
                scratch.name, cv, weight)
            nreplica = 1
            if dimension is None:
                dimension = ndim
        if dimension is None:
            dimension = len(grids)
        if len(grids) != dimension:
            raise GenesisValidationError(
                f"pmf_analysis: grids must have one (min, max, num_grids) entry "
                f"per dimension (expected {dimension}, got {len(grids)})."
            )

        band_width = _as_per_dim(band_width, "band_width", dimension)
        is_periodic = _as_per_dim(is_periodic, "is_periodic", dimension)
        box_size = _as_per_dim(box_size, "box_size", dimension)

        _validate_cvfiles_exist("pmf_analysis", cvfile, nreplica)

        ctrl = io.BytesIO()
        ctrl_files.write_ctrl_input(
                ctrl,
                cvfile=cvfile,
                weightfile=weightfile,
                distfile=distfile,
                )
        ctrl_files.write_ctrl_output(
                ctrl,
                pmffile=os.path.join(scratch.name, "pmf.dat"))
        ctrl.write(b'[OPTION]\n')
        ctrl_files.write_kwargs(
                ctrl,
                check_only=check_only,
                allow_backup=allow_backup,
                nreplica=nreplica,
                dimension=dimension,
                temperature=temperature,
                cutoff=cutoff,
                grids=ctrl_files.NumberingData(grids),
                band_width=ctrl_files.NumberingData(band_width),
                is_periodic=ctrl_files.NumberingData(is_periodic),
                box_size=ctrl_files.NumberingData(box_size),
                )

        ctrl_bytes, ctrl_len = ctrl_to_bytes(ctrl)
        with fortran_status() as (status, msgbuf, msglen):
            LibGenesis().lib.pmf_analysis_c(
                    ctrl_bytes,
                    ctrl_len,
                    ctypes.byref(result_pmf_c),
                    ctypes.byref(n_out1),
                    ctypes.byref(n_out2),
                    ctypes.byref(status),
                    msgbuf,
                    ctypes.c_int(msglen),
                    )

        result = c2py_util.conv_double_ndarray(
                result_pmf_c, [n_out1.value, n_out2.value])

        if dimension == 1:
            # Columns: bin center, standard PMF, Gaussian-kernel PMF.
            return Pmf1DResult(
                cv=result[:, 0],
                pmf=result[:, 1],
                pmf_gaussian=result[:, 2],
            )

        return Pmf2DResult(
            cv1=_bin_centers(grids[0]),
            cv2=_bin_centers(grids[1]),
            pmf=result,
        )
    finally:
        scratch.cleanup()
        if result_pmf_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_pmf_c),
                    ctypes.byref(n_out1), ctypes.byref(n_out2))


def pmf_analysis_isolated(timeout: Optional[float] = None, **kwargs):
    """Run :func:`pmf_analysis` in an isolated subprocess.

    Like the WHAM/MBAR isolated variants, this runs the estimate in a throwaway
    subprocess so accumulated Fortran state cannot destabilize the host kernel
    across many sequential calls.

    Args:
        timeout: Maximum time in seconds to wait for completion (None = no limit).
        **kwargs: All arguments passed to :func:`pmf_analysis`.

    Returns:
        The :class:`Pmf1DResult` or :class:`Pmf2DResult`, identical to what
        :func:`pmf_analysis` returns.
    """
    return run_analysis_isolated(
        func_name="pmf_analysis",
        timeout=timeout,
        task_description="pmf_analysis",
        **kwargs,
    )


def wham_analysis_isolated(timeout: Optional[float] = None, **kwargs):
    """Run :func:`wham_analysis` in an isolated subprocess.

    The WHAM/MBAR solvers accumulate global Fortran state, so running four or
    more of them in the same interpreter crashes the kernel. This variant runs
    the solve in a throwaway subprocess, giving every call clean state so any
    number of estimates can be computed back to back in one kernel session.

    Args:
        timeout: Maximum time in seconds to wait for completion (None = no limit).
        **kwargs: All arguments passed to :func:`wham_analysis`.

    Returns:
        The PMF array, identical to what :func:`wham_analysis` returns.

    Raises:
        GenesisValidationError: If input validation fails.
        GenesisFortranNotSupportedError: If an unsupported input mode is used.
        GenesisFortranError: If the Fortran solver returns an error.
        TimeoutError: If the solve exceeds ``timeout``.
        RuntimeError: If the subprocess fails unexpectedly.
    """
    return run_analysis_isolated(
        func_name="wham_analysis",
        timeout=timeout,
        task_description="wham_analysis",
        **kwargs,
    )


def mbar_analysis_isolated(timeout: Optional[float] = None, **kwargs):
    """Run :func:`mbar_analysis` in an isolated subprocess.

    The WHAM/MBAR solvers accumulate global Fortran state, so running four or
    more of them in the same interpreter crashes the kernel. This variant runs
    the solve in a throwaway subprocess, giving every call clean state so any
    number of estimates can be computed back to back in one kernel session.

    Args:
        timeout: Maximum time in seconds to wait for completion (None = no limit).
        **kwargs: All arguments passed to :func:`mbar_analysis`.

    Returns:
        The relative free energy array, identical to what
        :func:`mbar_analysis` returns.

    Raises:
        GenesisValidationError: If input validation fails.
        GenesisFortranNotSupportedError: If an unsupported input mode is used.
        GenesisFortranError: If the Fortran solver returns an error.
        TimeoutError: If the solve exceeds ``timeout``.
        RuntimeError: If the subprocess fails unexpectedly.
    """
    return run_analysis_isolated(
        func_name="mbar_analysis",
        timeout=timeout,
        task_description="mbar_analysis",
        **kwargs,
    )
