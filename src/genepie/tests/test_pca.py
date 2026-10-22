"""Tests of the PCA pipeline: flccrd, eigmat, prjcrd and pca_analysis.

The CLI reference data in ``tests/regression_test/test_analysis`` were produced
by the command-line tools, so the first group of tests pins the Python wrappers
to the CLI numerics. The remaining tests check lazy == memory, the end-to-end
``pca_analysis`` against the mdtraj + scikit-learn recipe, and the validation
errors.
"""
import copy

import numpy as np
import pytest

from .conftest import (
    ANALYSIS_REF_ROOT, BPTI_DCD, BPTI_PDB, BPTI_PSF,
    has_module, requires_mdtraj, requires_regression_data,
)
from ..s_molecule import SMolecule
from ..exceptions import GenesisValidationError, GenesisFortranError
from .. import genesis_exe

FLCCRD_REF = ANALYSIS_REF_ROOT / "test_flccrd_analysis" / "BPTI"
PRJCRD_REF = ANALYSIS_REF_ROOT / "test_prjcrd_analysis" / "BPTI"
EIGMAT_REF = ANALYSIS_REF_ROOT / "test_eigmat_analysis" / "APP"

requires_sklearn = pytest.mark.skipif(
    not has_module("sklearn"), reason="scikit-learn is not installed")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _load_bpti(lazy=False, selection="an:CA"):
    """BPTI C-alpha trajectory (58 atoms) and the matching subset molecule.

    The bundled system is solvated (37k atoms); the CA subset keeps every
    covariance matrix small while exercising the same code paths.
    """
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF)
    trajs, sub = genesis_exe.crd_convert(
        mol, trj_files=[str(BPTI_DCD)], trj_format="DCD",
        trj_type="COOR+BOX", selection=selection, lazy=lazy)
    return sub, trajs[0]


def _pdb_coords(path):
    """Coordinates of the ATOM records of a PDB file, shape (n, 3)."""
    xyz = []
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                xyz.append([float(line[30:38]), float(line[38:46]),
                            float(line[46:54])])
    return np.array(xyz)


def _cli_ave_coord(mol, ca_pdb):
    """Average structure of every atom of ``mol`` from the CLI pdb_avefile
    (which holds the CA atoms); non-CA atoms, if any, are left at zero."""
    ca_idx = genesis_exe.selection(mol, "an:CA") - 1   # 0-based
    ave = np.zeros((mol.num_atoms, 3))
    ave[ca_idx] = _pdb_coords(ca_pdb)
    return ave


# ---------------------------------------------------------------------------
# CLI reference data
# ---------------------------------------------------------------------------
@requires_regression_data
def test_flccrd_matches_cli_reference():
    """RMSF column of the CLI rmsfile; the VCV matrix must be symmetric and
    reproduce the RMSF from its diagonal."""
    mol, traj = _load_bpti()
    ave = _cli_ave_coord(mol, FLCCRD_REF / "input_ave.pdb")
    res = genesis_exe.flccrd_analysis(mol, traj, ave, "an:CA")

    ref = np.loadtxt(FLCCRD_REF / "ref", usecols=(0, 1, 3, 4))
    nsel = ref.shape[0]
    assert res.vcv_matrix.shape == (3 * nsel, 3 * nsel)
    assert res.rmsf.shape == (nsel,)
    np.testing.assert_allclose(res.rmsf, ref[:, 2], rtol=0, atol=1e-3)
    np.testing.assert_allclose(res.vcv_matrix, res.vcv_matrix.T,
                               rtol=0, atol=1e-12)
    diag_rmsf = np.sqrt(np.diag(res.vcv_matrix).reshape(nsel, 3).sum(axis=1))
    np.testing.assert_allclose(diag_rmsf, res.rmsf, rtol=1e-12, atol=0)
    bfactor = (8.0 / 3.0) * np.pi ** 2 * res.rmsf ** 2
    np.testing.assert_allclose(bfactor, ref[:, 3], rtol=0, atol=2e-3)


@requires_regression_data
def test_eigmat_matches_cli_reference():
    """Eigenvalues of the CLI valfile; eigenvectors orthonormal and
    consistent with the matrix."""
    with open(EIGMAT_REF / "input.pca") as f:
        natom = int(f.readline())
        values = np.loadtxt(f)
    n = 3 * natom
    vcv = values.reshape(n, n)

    res = genesis_exe.eigmat_analysis(vcv)
    ref = np.loadtxt(EIGMAT_REF / "ref")
    np.testing.assert_allclose(res.eigenvalues, ref, rtol=1e-8, atol=1e-8)
    assert np.all(np.diff(res.eigenvalues) <= 1e-9)        # descending
    np.testing.assert_allclose(res.eigenvectors @ res.eigenvectors.T,
                               np.eye(n), rtol=0, atol=1e-10)
    # vcv v_k = lambda_k v_k
    recon = res.eigenvectors @ vcv
    np.testing.assert_allclose(recon, res.eigenvalues[:, None] * res.eigenvectors,
                               rtol=0, atol=1e-8)
    np.testing.assert_allclose(res.contribution.sum(), 1.0, rtol=0, atol=1e-12)


@requires_regression_data
def test_prjcrd_matches_cli_reference():
    """Projections of the CLI prjfile from the CLI vecfile/valfile inputs."""
    mol, traj = _load_bpti()
    ave = _cli_ave_coord(mol, PRJCRD_REF / "input_ave.pdb")
    nvcv = 3 * len(genesis_exe.selection(mol, "an:CA"))
    # vecfile order: for mode k, for component j -> eigvec(j, k)
    vec = np.loadtxt(PRJCRD_REF / "input.vec").reshape(nvcv, nvcv)

    res = genesis_exe.prjcrd_analysis(mol, traj, ave, vec, "an:CA", num_pca=3)
    ref = np.loadtxt(PRJCRD_REF / "ref")
    assert res.projection.shape == (traj.nframe, 3)
    np.testing.assert_allclose(res.projection, ref[:, 1:], rtol=0, atol=1e-4)


# ---------------------------------------------------------------------------
# End-to-end pipeline
# ---------------------------------------------------------------------------
def test_pca_analysis_pipeline_consistency():
    """pca_analysis must equal the four steps run by hand, and the
    projections must have the covariance the eigenvalues predict."""
    mol, traj = _load_bpti(selection="an:CA")
    pca = genesis_exe.pca_analysis(mol, traj, analysis_selection="all",
                                   num_pca=4, num_iterations=10)

    ave = genesis_exe.avecrd_analysis(
        mol, traj, selection_group=["all"], fitting_method="TR+ROT",
        fitting_atom=1, num_iterations=10, analysis_atom=1)
    np.testing.assert_allclose(pca.ave_coord, ave.coord, rtol=0, atol=1e-12)
    flc = genesis_exe.flccrd_analysis(mol, traj, ave.coord, "all")
    eig = genesis_exe.eigmat_analysis(flc.vcv_matrix)
    prj = genesis_exe.prjcrd_analysis(mol, traj, ave.coord, eig.eigenvectors,
                                      "all", num_pca=4)
    np.testing.assert_allclose(pca.vcv_matrix, flc.vcv_matrix, rtol=0, atol=1e-12)
    np.testing.assert_allclose(pca.eigenvalues, eig.eigenvalues, rtol=0, atol=1e-12)
    np.testing.assert_allclose(pca.projection, prj.projection, rtol=0, atol=1e-12)

    # the average structure is the origin of the PC space ...
    assert pca.projection.shape == (traj.nframe, 4)
    np.testing.assert_allclose(pca.projection.mean(axis=0), 0.0, rtol=0, atol=1e-8)
    # ... and the variance along mode k is its eigenvalue (1/N normalisation)
    var = (pca.projection ** 2).mean(axis=0)
    np.testing.assert_allclose(var, pca.eigenvalues[:4], rtol=1e-8, atol=1e-10)
    # six rigid-body modes are removed by the fitting
    assert np.all(np.abs(pca.eigenvalues[-6:]) < 1e-8)


def test_pca_fitting_selection_differs_from_analysis():
    """Fit on the first 30 CA atoms, analyse all CA atoms: GLOBAL spans the
    analysis atoms, LOCAL spans the fitting atoms."""
    mol, traj = _load_bpti(selection="an:CA")
    fit_sel = "resno:1-30"
    n_all = mol.num_atoms
    n_fit = len(genesis_exe.selection(mol, fit_sel))
    assert 0 < n_fit < n_all

    g = genesis_exe.pca_analysis(mol, traj, analysis_selection="all",
                                 fitting_selection=fit_sel, num_pca=2,
                                 num_iterations=3)
    assert g.eigenvalues.shape == (3 * n_all,)
    assert g.rmsf.shape == (n_all,)
    assert g.projection.shape == (traj.nframe, 2)

    l = genesis_exe.pca_analysis(mol, traj, analysis_selection="all",
                                 fitting_selection=fit_sel, num_pca=2,
                                 num_iterations=3, vcv_matrix="LOCAL")
    assert l.eigenvalues.shape == (3 * n_fit,)
    assert l.rmsf.shape == (n_fit,)


@pytest.mark.parametrize("ana_period", [1, 2])
def test_pca_lazy_vs_memory(ana_period):
    """Lazy DCD input must give the same PCA as the in-memory trajectory."""
    mol_m, traj_m = _load_bpti(lazy=False)
    mol_l, traj_l = _load_bpti(lazy=True)
    assert traj_l.is_lazy

    kw = dict(analysis_selection="an:CA", num_pca=3,
              ana_period=ana_period)
    mem = genesis_exe.pca_analysis(mol_m, traj_m, **kw)
    lazy = genesis_exe.pca_analysis(mol_l, traj_l, **kw)

    assert mem.projection.shape == (traj_m.nframe // ana_period, 3)
    np.testing.assert_allclose(lazy.ave_coord, mem.ave_coord, rtol=0, atol=1e-10)
    np.testing.assert_allclose(lazy.vcv_matrix, mem.vcv_matrix, rtol=0, atol=1e-10)
    np.testing.assert_allclose(lazy.eigenvalues, mem.eigenvalues, rtol=0, atol=1e-10)
    np.testing.assert_allclose(lazy.projection, mem.projection, rtol=0, atol=1e-8)


@requires_mdtraj
@requires_sklearn
def test_pca_matches_mdtraj_sklearn_recipe():
    """Same result as 'superpose onto the iteratively refined mean, then
    sklearn PCA' (tutorial_md/analyze/05_md_dimensionalreduction.ipynb)."""
    import mdtraj as md
    from sklearn.decomposition import PCA

    mol, traj = _load_bpti(selection="an:CA")
    pca = genesis_exe.pca_analysis(mol, traj, analysis_selection="all",
                                   num_pca=3)

    w = md.load(str(BPTI_DCD), top=str(BPTI_PDB))
    w.restrict_atoms([a.index for a in w.topology.atoms if a.name == "CA"])
    ref = copy.copy(w[0])
    for _ in range(10):
        w.superpose(ref)
        ref.xyz = np.mean(w.xyz, axis=0, keepdims=True)
    w.superpose(ref, 0)
    X = w.xyz.reshape(w.n_frames, -1) * 10.0          # nm -> angstrom
    sk = PCA().fit(X)
    skp = sk.transform(X)
    nf = w.n_frames

    # sklearn normalises by N-1, GENESIS by N
    np.testing.assert_allclose(pca.eigenvalues[:3],
                               sk.explained_variance_[:3] * (nf - 1) / nf,
                               rtol=1e-4, atol=0)
    np.testing.assert_allclose(pca.contribution[:3],
                               sk.explained_variance_ratio_[:3],
                               rtol=1e-4, atol=0)
    for k in range(3):                                 # sign is arbitrary
        a, b = pca.projection[:, k], skp[:, k]
        sign = 1.0 if np.dot(a, b) >= 0 else -1.0
        np.testing.assert_allclose(a, sign * b, rtol=0, atol=1e-4)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
def test_pca_validation_errors():
    mol, traj = _load_bpti(selection="an:CA")
    ave = np.zeros((mol.num_atoms, 3))

    with pytest.raises(GenesisValidationError):
        genesis_exe.eigmat_analysis(np.zeros((3, 4)))
    with pytest.raises(GenesisValidationError):
        genesis_exe.eigmat_analysis(np.array([[0.0, 1.0], [2.0, 0.0]]))
    with pytest.raises(GenesisValidationError):
        genesis_exe.flccrd_analysis(mol, traj, np.zeros((5, 3)), "all")
    with pytest.raises(GenesisValidationError):
        genesis_exe.flccrd_analysis(mol, traj, ave, "all", vcv_matrix="BOTH")
    with pytest.raises(GenesisValidationError):
        genesis_exe.flccrd_analysis(mol, traj, ave, "all", ana_period=0)
    n = 3 * mol.num_atoms
    with pytest.raises(GenesisValidationError):
        genesis_exe.prjcrd_analysis(mol, traj, ave, np.eye(n), "all",
                                    num_pca=n + 1)
    with pytest.raises(GenesisValidationError):
        genesis_exe.prjcrd_analysis(mol, traj, ave, np.eye(n - 1), "all")
    with pytest.raises(GenesisValidationError):
        genesis_exe.pca_analysis(mol, traj, num_iterations=0)
    with pytest.raises(GenesisValidationError):
        genesis_exe.pca_analysis(mol, traj, analysis_selection="")


def test_pca_fortran_errors_are_typed():
    """Failures inside the Fortran core surface as typed exceptions, not
    as a process exit."""
    mol, traj = _load_bpti(selection="an:CA")
    ave = np.zeros((mol.num_atoms, 3))
    # the molecule has only CA atoms: an empty analysis selection
    with pytest.raises((GenesisFortranError, GenesisValidationError)):
        genesis_exe.flccrd_analysis(mol, traj, ave, "an:CB")
    # trajectory/molecule atom-count mismatch
    full = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF)
    with pytest.raises((GenesisFortranError, GenesisValidationError)):
        genesis_exe.flccrd_analysis(full, traj, np.zeros((full.num_atoms, 3)),
                                    "an:CA")
