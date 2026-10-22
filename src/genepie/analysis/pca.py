"""Principal component analysis of Cartesian coordinates.

The GENESIS PCA pipeline is four command-line tools, each of which is wrapped
here by a function that shares its Fortran core with the CLI program:

1. :func:`~genepie.analysis.avecrd.avecrd_analysis` -- iterative average
   structure (fit every frame onto the current average, average, repeat).
2. :func:`flccrd_analysis` -- fit every frame onto the converged average and
   accumulate the variance-covariance matrix of the (mass-weighted)
   fluctuations.
3. :func:`eigmat_analysis` -- diagonalise the matrix (LAPACK ``DSYEV``).
4. :func:`prjcrd_analysis` -- fit every frame onto the average again and
   project the fluctuations onto the principal modes.

:func:`pca_analysis` runs the whole pipeline. The fitting removes overall
translation and rotation before the covariance is accumulated, exactly as the
"superpose onto the iteratively refined mean structure, then PCA" recipe does
with mdtraj + scikit-learn; see the notes in :func:`pca_analysis` for the two
small conventions that differ.
"""

import ctypes
import io
from collections import namedtuple
from typing import Optional

import numpy as np

from ..libgenesis import LibGenesis
from ..s_molecule import SMolecule
from ..s_trajectories import STrajectories
from ..exceptions import GenesisValidationError
from .. import ctrl_files
from .._fortran import (
    ctrl_to_bytes,
    fortran_status,
    molecule_c,
)
from ._common import _SourceArgs
from .avecrd import avecrd_analysis
from .converter import selection


_VCV_MATRIX_MODES = ("GLOBAL", "LOCAL")


FlccrdAnalysisResult = namedtuple(
        'FlccrdAnalysisResult',
        ['vcv_matrix', 'rmsf'])
FlccrdAnalysisResult.__doc__ = """Result of :func:`flccrd_analysis`.

Attributes:
    vcv_matrix: Variance-covariance matrix of the fitted (mass-weighted)
        fluctuations, shape ``(3*nsel, 3*nsel)``; element ``3*i + d`` is
        coordinate ``d`` of the ``i``-th analysed atom. Normalised by the
        number of frames, as ``flccrd_analysis`` does.
    rmsf: Root-mean-square fluctuation of each analysed atom in angstrom,
        shape ``(nsel,)`` (``rmsfile`` column 4 of the CLI tool).
"""

EigmatAnalysisResult = namedtuple(
        'EigmatAnalysisResult',
        ['eigenvalues', 'eigenvectors', 'contribution'])
EigmatAnalysisResult.__doc__ = """Result of :func:`eigmat_analysis`.

Attributes:
    eigenvalues: Eigenvalues in descending order, shape ``(n,)``.
    eigenvectors: Eigenvectors, shape ``(n, n)``; ``eigenvectors[k]`` is the
        ``k``-th principal mode (the ``k``-th column of the CLI ``vecfile``).
    contribution: Fraction of the total variance carried by each mode,
        ``eigenvalues / eigenvalues.sum()`` (the explained variance ratio).
"""

PrjcrdAnalysisResult = namedtuple(
        'PrjcrdAnalysisResult',
        ['projection'])
PrjcrdAnalysisResult.__doc__ = """Result of :func:`prjcrd_analysis`.

Attributes:
    projection: Projection of every analysed frame onto the first ``num_pca``
        modes, shape ``(nframe, num_pca)`` in angstrom (``prjfile`` columns
        2.. of the CLI tool).
"""

PcaAnalysisResult = namedtuple(
        'PcaAnalysisResult',
        ['projection', 'eigenvalues', 'eigenvectors', 'contribution',
         'rmsf', 'ave_coord', 'vcv_matrix'])
PcaAnalysisResult.__doc__ = """Result of :func:`pca_analysis`.

Attributes:
    projection: ``(nframe, num_pca)`` projections onto the principal modes.
    eigenvalues: ``(3*nsel,)`` eigenvalues of the covariance matrix,
        descending.
    eigenvectors: ``(3*nsel, 3*nsel)``; ``eigenvectors[k]`` is mode ``k``.
    contribution: Explained variance ratio of every mode.
    rmsf: ``(nsel,)`` RMSF of the analysed atoms about the average, angstrom.
    ave_coord: ``(natom, 3)`` converged average structure of every atom of
        ``molecule`` (the reference the frames were fitted onto).
    vcv_matrix: ``(3*nsel, 3*nsel)`` variance-covariance matrix.
"""


def _check_vcv_mode(vcv_matrix: str) -> str:
    mode = str(vcv_matrix).upper()
    if mode not in _VCV_MATRIX_MODES:
        raise GenesisValidationError(
            f"Invalid vcv_matrix: {vcv_matrix}. "
            f"Valid values: {list(_VCV_MATRIX_MODES)}"
        )
    return mode


def _ave_coord_fortran(molecule: SMolecule, ave_coord) -> np.ndarray:
    """Return ``ave_coord`` as the C-contiguous (natom, 3) float64 array that
    aliases the Fortran (3, natom) buffer."""
    arr = np.ascontiguousarray(ave_coord, dtype=np.float64)
    if arr.shape != (molecule.num_atoms, 3):
        raise GenesisValidationError(
            f"ave_coord must have shape ({molecule.num_atoms}, 3), "
            f"got {arr.shape}"
        )
    if not np.all(np.isfinite(arr)):
        raise GenesisValidationError("ave_coord contains non-finite values")
    return arr


def _pca_groups(analysis_selection: str, fitting_selection: Optional[str]):
    """Selection groups and the 1-based group indices of a PCA control text."""
    if not analysis_selection:
        raise GenesisValidationError("analysis_selection must not be empty")
    groups = [analysis_selection]
    if fitting_selection is None or fitting_selection == analysis_selection:
        fitting_group = 1
    else:
        groups.append(fitting_selection)
        fitting_group = 2
    return groups, 1, fitting_group


def _write_pca_ctrl(
        analysis_selection: str,
        fitting_selection: Optional[str],
        fitting_method: str,
        mass_weighted: bool,
        vcv_matrix: str,
        num_pca: Optional[int] = None,
        ) -> io.BytesIO:
    groups, analysis_group, fitting_group = _pca_groups(
        analysis_selection, fitting_selection)
    ctrl = io.BytesIO()
    ctrl_files.write_ctrl_selection(ctrl, groups)
    ctrl_files.write_ctrl_fitting(
        ctrl, fitting_method=fitting_method, fitting_atom=fitting_group,
        mass_weight=mass_weighted)
    ctrl.write(b"[OPTION]\n")
    ctrl_files.write_kwargs(
        ctrl,
        check_only=False,
        analysis_atom=analysis_group,
        vcv_matrix=vcv_matrix,
        num_pca=num_pca,
    )
    return ctrl


def _analysed_atoms(molecule: SMolecule, analysis_selection: str,
                    fitting_selection: Optional[str], vcv_matrix: str):
    """1-based indices of the atoms that span the PC space."""
    if vcv_matrix == "GLOBAL" or fitting_selection is None:
        return selection(molecule, analysis_selection)
    return selection(molecule, fitting_selection)


def flccrd_analysis(
        molecule: SMolecule,
        trajs: STrajectories,
        ave_coord,
        analysis_selection: str,
        fitting_selection: Optional[str] = None,
        fitting_method: str = "TR+ROT",
        mass_weighted: bool = False,
        vcv_matrix: str = "GLOBAL",
        ana_period: int = 1,
        ) -> FlccrdAnalysisResult:
    """
    Variance-covariance matrix of the fluctuations about an average structure
    (``flccrd_analysis``).

    Every frame is fitted onto ``ave_coord`` (translation and rotation removed
    with ``fitting_method`` on the ``fitting_selection`` atoms), then the
    mass-weighted deviations of the analysed atoms from ``ave_coord`` are
    accumulated into the covariance matrix.

    Args:
        molecule: Molecular structure (must match ``trajs``)
        trajs: Trajectory, in memory or lazy
        ave_coord: Average structure of every atom, shape ``(natom, 3)``;
            typically ``avecrd_analysis(...).coord``
        analysis_selection: Atoms spanning the PC space (``GLOBAL`` mode),
            e.g. ``"an:CA"``
        fitting_selection: Atoms used for fitting; default: the analysis
            atoms
        fitting_method: ``NO``, ``TR+ROT`` (default), ``TR``, ``TR+ZROT``,
            ``XYTR`` or ``XYTR+ZROT``
        mass_weighted: Mass-weighted fitting and fluctuations. Masses are
            assigned from the atom names, as the CLI tool does.
        vcv_matrix: ``GLOBAL`` (analysis atoms span the PC space, default) or
            ``LOCAL`` (fitting atoms span it)
        ana_period: Analysis period (default: 1)

    Returns:
        FlccrdAnalysisResult(vcv_matrix, rmsf)
    """
    vcv_mode = _check_vcv_mode(vcv_matrix)
    src = _SourceArgs(trajs, ana_period)
    ave = _ave_coord_fortran(molecule, ave_coord)
    nsel = len(_analysed_atoms(molecule, analysis_selection,
                               fitting_selection, vcv_mode))
    nvcv = 3 * nsel

    # Symmetric, so the C-order (nvcv, nvcv) buffer aliases the Fortran one.
    vcv = np.zeros((nvcv, nvcv), dtype=np.float64)
    rmsf = np.zeros(nsel, dtype=np.float64)
    nstru_out = ctypes.c_int(0)

    ctrl = _write_pca_ctrl(analysis_selection, fitting_selection,
                           fitting_method, mass_weighted, vcv_mode)
    ctrl_bytes, ctrl_len = ctrl_to_bytes(ctrl)

    with molecule_c(molecule) as mol_c:
        with fortran_status() as (status, msg, msglen):
            LibGenesis().lib.flccrd_analysis_c(
                ctypes.byref(mol_c),
                *src.args,
                ctrl_bytes,
                ctrl_len,
                ave.ctypes.data_as(ctypes.c_void_p),
                vcv.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(nvcv),
                rmsf.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(nsel),
                ctypes.byref(nstru_out),
                ctypes.byref(status),
                msg,
                ctypes.c_int(msglen),
            )

    return FlccrdAnalysisResult(vcv, rmsf)


def eigmat_analysis(vcv_matrix) -> EigmatAnalysisResult:
    """
    Eigen decomposition of a variance-covariance matrix (``eigmat_analysis``).

    Args:
        vcv_matrix: Symmetric ``(n, n)`` matrix, e.g.
            ``flccrd_analysis(...).vcv_matrix``

    Returns:
        EigmatAnalysisResult(eigenvalues, eigenvectors, contribution) with the
        eigenvalues in descending order and ``eigenvectors[k]`` the ``k``-th
        mode.

    Note:
        The CLI ``cntfile`` reports contributions relative to the sum of the
        first ``n - 7`` eigenvalues; ``contribution`` here is relative to the
        sum of all ``n``. After rigid-body fitting the trailing eigenvalues
        are ~0, so the two agree to rounding.
    """
    vcv = np.ascontiguousarray(vcv_matrix, dtype=np.float64)
    if vcv.ndim != 2 or vcv.shape[0] != vcv.shape[1] or vcv.shape[0] == 0:
        raise GenesisValidationError(
            f"vcv_matrix must be a non-empty square matrix, got {vcv.shape}"
        )
    if not np.all(np.isfinite(vcv)):
        raise GenesisValidationError("vcv_matrix contains non-finite values")
    if not np.allclose(vcv, vcv.T, rtol=1e-8, atol=1e-10):
        raise GenesisValidationError("vcv_matrix must be symmetric")
    n = vcv.shape[0]

    eigval = np.zeros(n, dtype=np.float64)
    # Fortran (j, k) = component j of mode k aliases Python eigvec[k, j].
    eigvec = np.zeros((n, n), dtype=np.float64)

    with fortran_status() as (status, msg, msglen):
        LibGenesis().lib.eigmat_analysis_c(
            vcv.ctypes.data_as(ctypes.c_void_p),
            ctypes.c_int(n),
            eigval.ctypes.data_as(ctypes.c_void_p),
            eigvec.ctypes.data_as(ctypes.c_void_p),
            ctypes.byref(status),
            msg,
            ctypes.c_int(msglen),
        )

    total = eigval.sum()
    contribution = eigval / total if total != 0.0 else np.zeros_like(eigval)
    return EigmatAnalysisResult(eigval, eigvec, contribution)


def prjcrd_analysis(
        molecule: SMolecule,
        trajs: STrajectories,
        ave_coord,
        eigenvectors,
        analysis_selection: str,
        num_pca: Optional[int] = None,
        fitting_selection: Optional[str] = None,
        fitting_method: str = "TR+ROT",
        mass_weighted: bool = False,
        vcv_matrix: str = "GLOBAL",
        ana_period: int = 1,
        ) -> PrjcrdAnalysisResult:
    """
    Projection of a trajectory onto principal modes (``prjcrd_analysis``).

    Every frame is fitted onto ``ave_coord`` exactly as in
    :func:`flccrd_analysis`, and the mass-weighted deviation of the analysed
    atoms is projected onto the first ``num_pca`` rows of ``eigenvectors``.
    ``trajs`` need not be the trajectory the modes were derived from, so this
    also projects new data onto an existing PC space.

    Args:
        molecule: Molecular structure (must match ``trajs``)
        trajs: Trajectory, in memory or lazy
        ave_coord: Average structure of every atom, shape ``(natom, 3)``
        eigenvectors: ``(nmode, 3*nsel)`` array whose row ``k`` is mode ``k``
            (``eigmat_analysis(...).eigenvectors``)
        analysis_selection: Atoms spanning the PC space (``GLOBAL`` mode)
        num_pca: Number of modes to project onto (default: all rows of
            ``eigenvectors``)
        fitting_selection, fitting_method, mass_weighted, vcv_matrix,
        ana_period: As in :func:`flccrd_analysis`; use the same values the
            covariance matrix was computed with.

    Returns:
        PrjcrdAnalysisResult(projection) with shape ``(nframe, num_pca)``
    """
    vcv_mode = _check_vcv_mode(vcv_matrix)
    src = _SourceArgs(trajs, ana_period)
    ave = _ave_coord_fortran(molecule, ave_coord)
    nsel = len(_analysed_atoms(molecule, analysis_selection,
                               fitting_selection, vcv_mode))
    nvcv = 3 * nsel

    # Python (nmode, nvcv) C-order aliases Fortran (nvcv, nmode).
    vec = np.ascontiguousarray(eigenvectors, dtype=np.float64)
    if vec.ndim != 2 or vec.shape[1] != nvcv or vec.shape[0] == 0:
        raise GenesisValidationError(
            f"eigenvectors must have shape (nmode, {nvcv}), got {vec.shape}"
        )
    nmode = vec.shape[0]
    if num_pca is None:
        num_pca = nmode
    if not isinstance(num_pca, (int, np.integer)) or num_pca < 1 \
            or num_pca > nmode:
        raise GenesisValidationError(
            f"num_pca must be an integer in 1..{nmode}, got {num_pca}"
        )
    num_pca = int(num_pca)

    # Python (nframe, num_pca) C-order aliases Fortran (num_pca, nframe).
    proj = np.zeros((src.n_frames, num_pca), dtype=np.float64)
    nstru_out = ctypes.c_int(0)

    ctrl = _write_pca_ctrl(analysis_selection, fitting_selection,
                           fitting_method, mass_weighted, vcv_mode,
                           num_pca=num_pca)
    ctrl_bytes, ctrl_len = ctrl_to_bytes(ctrl)

    with molecule_c(molecule) as mol_c:
        with fortran_status() as (status, msg, msglen):
            LibGenesis().lib.prjcrd_analysis_c(
                ctypes.byref(mol_c),
                *src.args,
                ctrl_bytes,
                ctrl_len,
                ave.ctypes.data_as(ctypes.c_void_p),
                vec.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(nvcv),
                ctypes.c_int(nmode),
                proj.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_int(num_pca),
                ctypes.c_int(src.n_frames),
                ctypes.byref(nstru_out),
                ctypes.byref(status),
                msg,
                ctypes.c_int(msglen),
            )

    return PrjcrdAnalysisResult(proj[:nstru_out.value])


def pca_analysis(
        molecule: SMolecule,
        trajs: STrajectories,
        analysis_selection: str = "an:CA",
        fitting_selection: Optional[str] = None,
        fitting_method: str = "TR+ROT",
        mass_weighted: bool = False,
        num_pca: Optional[int] = None,
        num_iterations: int = 10,
        vcv_matrix: str = "GLOBAL",
        ana_period: int = 1,
        ) -> PcaAnalysisResult:
    """
    Principal component analysis of Cartesian coordinates.

    Runs the GENESIS pipeline ``avecrd_analysis`` -> ``flccrd_analysis`` ->
    ``eigmat_analysis`` -> ``prjcrd_analysis``:

    1. Iteratively fit every frame onto the current average structure and
       re-average (``num_iterations`` times, starting from the coordinates of
       ``molecule``), as the GENESIS ``avecrd_analysis`` tool does. This
       removes overall translation and rotation. The change of the average in
       each cycle is reported in the Fortran log (``check the convergence``).
    2. Fit every frame onto the converged average and accumulate the
       variance-covariance matrix of the analysed atoms.
    3. Diagonalise it (LAPACK ``DSYEV``).
    4. Fit every frame onto the average once more and project it onto the
       first ``num_pca`` modes.

    Lazy trajectories are supported; the DCD file is read once per pass.

    Args:
        molecule: Molecular structure matching ``trajs`` (for lazy
            trajectories, the subset molecule returned by ``crd_convert``);
            its coordinates seed the iterative average
        trajs: Trajectory, in memory or lazy
        analysis_selection: Atoms spanning the PC space (default: ``"an:CA"``)
        fitting_selection: Atoms used for fitting (default: the analysis
            atoms)
        fitting_method: ``TR+ROT`` (default) or another GENESIS fitting method
        mass_weighted: Mass-weighted fitting, averaging and fluctuations
        num_pca: Number of modes to project onto (default: all ``3*nsel``)
        num_iterations: Fit-and-average iterations of step 1 (default: 10;
            the average of a well-sampled trajectory typically converges to
            machine precision within a few cycles)
        vcv_matrix: ``GLOBAL`` (default) or ``LOCAL``, see
            :func:`flccrd_analysis`
        ana_period: Analysis period (default: 1)

    Returns:
        PcaAnalysisResult(projection, eigenvalues, eigenvectors,
        contribution, rmsf, ave_coord, vcv_matrix)

    Note:
        Compared with ``PCA().fit_transform`` on mdtraj-superposed
        coordinates, (a) the covariance is normalised by ``nframe`` rather
        than ``nframe - 1`` (eigenvalues scale by ``nframe/(nframe-1)``;
        eigenvectors and projections do not change) and (b) the fluctuations
        are taken about the converged average of step 1 rather than about the
        mean of the finally superposed frames; the two coincide once the
        iteration has converged.

    Examples:
        >>> trajs, mol = crd_convert(mol, ["traj.dcd"], selection="an:CA")
        >>> pca = pca_analysis(mol, trajs[0], analysis_selection="an:CA",
        ...                    num_pca=2)
        >>> pc1, pc2 = pca.projection[:, 0], pca.projection[:, 1]
        >>> print(pca.contribution[:2])
    """
    vcv_mode = _check_vcv_mode(vcv_matrix)
    if not isinstance(num_iterations, (int, np.integer)) or num_iterations < 1:
        raise GenesisValidationError(
            f"num_iterations must be a positive integer, got {num_iterations}"
        )
    groups, analysis_group, fitting_group = _pca_groups(
        analysis_selection, fitting_selection)

    # 1. iterative average structure (all atoms of molecule)
    ave = avecrd_analysis(
        molecule, trajs,
        ana_period=ana_period,
        selection_group=groups,
        fitting_method=fitting_method,
        fitting_atom=fitting_group,
        mass_weight=mass_weighted,
        num_iterations=int(num_iterations),
        analysis_atom=analysis_group,
    )

    # 2. covariance of the fitted fluctuations
    flc = flccrd_analysis(
        molecule, trajs, ave.coord, analysis_selection,
        fitting_selection=fitting_selection,
        fitting_method=fitting_method,
        mass_weighted=mass_weighted,
        vcv_matrix=vcv_mode,
        ana_period=ana_period,
    )

    # 3. principal modes
    eig = eigmat_analysis(flc.vcv_matrix)

    # 4. projections
    if num_pca is None:
        num_pca = eig.eigenvalues.shape[0]
    prj = prjcrd_analysis(
        molecule, trajs, ave.coord, eig.eigenvectors, analysis_selection,
        num_pca=num_pca,
        fitting_selection=fitting_selection,
        fitting_method=fitting_method,
        mass_weighted=mass_weighted,
        vcv_matrix=vcv_mode,
        ana_period=ana_period,
    )

    return PcaAnalysisResult(
        prj.projection, eig.eigenvalues, eig.eigenvectors, eig.contribution,
        flc.rmsf, ave.coord, flc.vcv_matrix)
