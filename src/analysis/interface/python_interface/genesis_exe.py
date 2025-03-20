import ctypes
from collections import namedtuple
import os
import tempfile
from typing import Iterable, Optional, TextIO
import numpy as np
import numpy.typing as npt
from libgenesis import LibGenesis
from s_molecule import SMolecule
from s_trajectories import STrajectories, STrajectoriesArray
import ctrl_files
import c2py_util
import py2c_util


def crd_convert(molecule: SMolecule,
                ctrl_filename: str | bytes | os.PathLike) -> STrajectoriesArray:
    """
    Executes crd_convert.

    Args:
        molecule:
        ctrl_filename:

    Returns:
        Array of STrajectories
    """
    buf = ctypes.c_void_p(None)
    num_trajs_c = ctypes.c_int(0)
    mol_c = molecule.to_SMoleculeC()
    LibGenesis().lib.crd_convert_c(
            ctypes.byref(mol_c),
            py2c_util.pathlike_to_byte(ctrl_filename),
            ctypes.byref(buf),
            ctypes.byref(num_trajs_c))
    LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))
    return STrajectoriesArray(buf, num_trajs_c.value)


TrjAnalysisResult = namedtuple(
        'TrjAnalysisResult',
        ['distance',
         'angle',
         'torsion',
         'com_distance',
         'com_angle',
         'com_torsion'])


def trj_analysis(molecule: SMolecule, trajs :STrajectories,
                 ana_period: Optional[int] = 1,
                 selection_group: Optional[Iterable[str]] = None,
                 selection_mole_name: Optional[Iterable[str]] = None,
                 check_only: Optional[bool] = None,
                 distance: Optional[Iterable[str]] = [],
                 dist_weight: Optional[Iterable[str]] = [],
                 angle: Optional[Iterable[str]] = [],
                 torsion: Optional[Iterable[str]] = [],
                 com_distance: Optional[Iterable[str]] = [],
                 com_angle: Optional[Iterable[str]] = [],
                 com_torsion: Optional[Iterable[str]] = [],
                 ) -> TrjAnalysisResult:
    """
    Executes trj_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
        selection_group:
        selection_mole_name:
    Returns:
        (distance, angle, torsion, cdis, cang, ctor)
    """
    num_distance = ctypes.c_int(0)
    result_distance_c = ctypes.c_void_p(None)
    num_angle = ctypes.c_int(0)
    result_angle_c = ctypes.c_void_p(None)
    num_torsion = ctypes.c_int(0)
    result_torsion_c = ctypes.c_void_p(None)
    num_cdis = ctypes.c_int(0)
    result_cdis_c = ctypes.c_void_p(None)
    num_cang = ctypes.c_int(0)
    result_cang_c = ctypes.c_void_p(None)
    num_ctor = ctypes.c_int(0)
    result_ctor_c = ctypes.c_void_p(None)
    ana_period_c = ctypes.c_int(ana_period)
    mol_c = ctypes.c_void_p(None)
    try:
        mol_c = molecule.to_SMoleculeC()
        with tempfile.NamedTemporaryFile(dir=os.getcwd(), delete=True) as ctrl:
            ctrl_files.write_ctrl_output(
                    ctrl,
                    disfile="dummy.dis",
                    angfile="dummy.ang",
                    torfile="dummy.tor",
                    comdisfile="dummy.comdis",
                    comangfile="dummy.comangr",
                    comtorfile="dummy.comtor",
                    qntfile="dummy.qnt",
                    )
            ctrl_files.write_ctrl_selection(
                    ctrl, selection_group, selection_mole_name)
            ctrl.write(b"[OPTION]\n")
            ctrl_files.write_kwargs(
                    ctrl,
                    check_only = check_only,
                    distance = distance,
                    dist_weight = dist_weight,
                    angle = angle,
                    torsion = torsion,
                    com_distance = com_distance,
                    com_angle = com_angle,
                    com_torsion = com_torsion,
                    )

            ctrl.seek(0)
            LibGenesis().lib.trj_analysis_c(
                    ctypes.byref(mol_c),
                    ctypes.byref(trajs.get_c_obj()),
                    ctypes.byref(ana_period_c),
                    py2c_util.pathlike_to_byte(ctrl.name),
                    ctypes.byref(result_distance_c),
                    ctypes.byref(num_distance),
                    ctypes.byref(result_angle_c),
                    ctypes.byref(num_angle),
                    ctypes.byref(result_torsion_c),
                    ctypes.byref(num_torsion),
                    ctypes.byref(result_cdis_c),
                    ctypes.byref(num_cdis),
                    ctypes.byref(result_cang_c),
                    ctypes.byref(num_cang),
                    ctypes.byref(result_ctor_c),
                    ctypes.byref(num_ctor),
                    )
        n_frame_c = ctypes.c_int(int(trajs.nframe / ana_period))
        result_distance = (c2py_util.conv_double_ndarray(
            result_distance_c, [n_frame_c.value, num_distance.value])
                           if result_distance_c else None)
        result_angle = (c2py_util.conv_double_ndarray(
                result_angle_c, [n_frame_c.value, num_angle.value])
                        if result_angle_c else None)
        result_torsion = (c2py_util.conv_double_ndarray(
                result_torsion_c, [n_frame_c.value, num_torsion.value])
                        if result_torsion_c else None)
        result_cdis = (c2py_util.conv_double_ndarray(
            result_cdis_c, [n_frame_c.value, num_cdis.value])
                           if result_distance_c else None)
        result_cang = (c2py_util.conv_double_ndarray(
                result_cang_c, [n_frame_c.value, num_cang.value])
                        if result_angle_c else None)
        result_ctor = (c2py_util.conv_double_ndarray(
                result_ctor_c, [n_frame_c.value, num_ctor.value])
                        if result_ctor_c else None)
        return TrjAnalysisResult(
                result_distance, result_angle, result_torsion,
                result_cdis, result_cang, result_ctor)
    finally:
        if result_ctor_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_ctor_c),
                    ctypes.byref(n_frame_c), ctypes.byref(num_ctor))
        if result_cang_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_cang_c),
                    ctypes.byref(n_frame_c), ctypes.byref(num_cang))
        if result_cdis_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_cdis_c),
                    ctypes.byref(n_frame_c), ctypes.byref(num_cdis))
        if result_torsion_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_torsion_c),
                    ctypes.byref(n_frame_c), ctypes.byref(num_torsion))
        if result_angle_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_angle_c),
                    ctypes.byref(n_frame_c), ctypes.byref(num_angle))
        if result_distance_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_distance_c),
                    ctypes.byref(n_frame_c), ctypes.byref(num_distance))
        if mol_c:
            LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))


RgAnalysisResult = namedtuple(
        'RgAnalysisResult',
        ['rg'])


def rg_analysis(molecule: SMolecule, trajs :STrajectories,
        ana_period: Optional[int] = 1,
        selection_group: Optional[Iterable[str]] = None,
        selection_mole_name: Optional[Iterable[str]] = None,
        fitting_method: Optional[str] = None,
        fitting_atom: Optional[int] = None,
        check_only: Optional[bool] = None,
        analysis_atom: Optional[int] = None,
        mass_weighted: Optional[bool] = None,
        ) -> RgAnalysisResult:
    """
    Executes rg_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
    Returns:
        rg
    """
    mol_c = molecule.to_SMoleculeC()
    ana_period_c = ctypes.c_int(ana_period)
    result_rg_c = ctypes.c_void_p()
    try:
        mol_c = molecule.to_SMoleculeC()
        with tempfile.NamedTemporaryFile(dir=os.getcwd(), delete=True) as ctrl:
            ctrl_files.write_ctrl_output(
                    ctrl,
                    rgfile = "dummy.rg")
            ctrl_files.write_ctrl_selection(
                    ctrl, selection_group, selection_mole_name)
            ctrl_files.write_ctrl_fitting(
                    ctrl, fitting_method, fitting_atom)
            ctrl.write(b"[OPTION]\n")
            ctrl_files.write_value(ctrl, "check_only", check_only)
            ctrl_files.write_value(ctrl, "analysis_atom", analysis_atom)
            ctrl_files.write_value(ctrl, "mass_weighted", mass_weighted)

            ctrl.seek(0)
            LibGenesis().lib.rg_analysis_c(
                    ctypes.byref(mol_c),
                    ctypes.byref(trajs.get_c_obj()),
                    ctypes.byref(ana_period_c),
                    py2c_util.pathlike_to_byte(ctrl.name),
                    ctypes.byref(result_rg_c),
                    )
        n_frame_c = ctypes.c_int(int(trajs.nframe / ana_period))
        result_rg = (c2py_util.conv_double_ndarray(
            result_rg_c, n_frame_c.value)
                    if result_rg_c else None)
        return RgAnalysisResult(
                result_rg)
    finally:
        if result_rg_c:
            LibGenesis().lib.deallocate_double(
                    ctypes.byref(result_rg_c),
                    ctypes.byref(n_frame_c))
        if mol_c:
            LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))


RmsdAnalysisResult = namedtuple(
        'RmsdAnalysisResult',
        ['rmsd'])


def rmsd_analysis(
        molecule: SMolecule, trajs :STrajectories,
        ana_period: Optional[int] = 1,
        selection_group: Optional[Iterable[str]] = None,
        selection_mole_name: Optional[Iterable[str]] = None,
        fitting_method: Optional[str] = None,
        fitting_atom: Optional[int] = None,
        check_only: Optional[bool] = None,
        analysis_atom: Optional[int] = None,
        ) -> RmsdAnalysisResult:
    """
    Executes rmsd_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
    Returns:
        rmsd
    """
    mol_c = None
    ana_period_c = ctypes.c_int(ana_period)
    result_rmsd_c = ctypes.c_void_p()
    try:
        mol_c = molecule.to_SMoleculeC()
        with tempfile.NamedTemporaryFile(dir=os.getcwd(), delete=True) as ctrl:
            ctrl_files.write_ctrl_output(
                    ctrl,
                    rmsfile = "dummy.rms")
            ctrl_files.write_ctrl_selection(
                    ctrl, selection_group, selection_mole_name)
            ctrl_files.write_ctrl_fitting(
                    ctrl, fitting_method, fitting_atom)
            ctrl.write(b"[OPTION]\n")
            ctrl_files.write_value(ctrl, "check_only", check_only)
            ctrl_files.write_value(ctrl, "analysis_atom", analysis_atom)

            ctrl.seek(0)
            LibGenesis().lib.ra_analysis_c(
                    ctypes.byref(mol_c),
                    ctypes.byref(trajs.get_c_obj()),
                    ctypes.byref(ana_period_c),
                    py2c_util.pathlike_to_byte(ctrl.name),
                    ctypes.byref(result_rmsd_c),
                    )
        n_frame_c = ctypes.c_int(int(trajs.nframe / ana_period))
        result_rmsd = (c2py_util.conv_double_ndarray(
            result_rmsd_c, n_frame_c.value)
                    if result_rmsd_c else None)
        return RmsdAnalysisResult(
                result_rmsd)
    finally:
        if result_rmsd_c:
            LibGenesis().lib.deallocate_double(
                    ctypes.byref(result_rmsd_c),
                    ctypes.byref(n_frame_c))
        if mol_c:
            LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))


DrmsAnalysisResult = namedtuple(
        'DrmsAnalysisResult',
        ['drms'])


def drms_analysis(
        molecule: SMolecule, trajs :STrajectories,
        ana_period: Optional[int] = 1,
        selection_group: Optional[Iterable[str]] = None,
        selection_mole_name: Optional[Iterable[str]] = None,
        fitting_method: Optional[str] = None,
        fitting_atom: Optional[int] = None,
        check_only: Optional[bool] = None,
        contact_groups: Optional[int] = None,
        ignore_hydrogen: Optional[bool] = None,
        two_states: Optional[bool] = None,
        avoid_bonding: Optional[bool] = None,
        exclude_residues: Optional[int] = None,
        minimum_distance: Optional[float] = None,
        maximum_distance: Optional[float] = None,
        pbc_correct: Optional[bool] = None,
        verbose: Optional[bool] = None,
        ) -> DrmsAnalysisResult:
    """
    Executes drms_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
    Returns:
        drms
    """
    mol_c = None
    ana_period_c = ctypes.c_int(ana_period)
    result_drms_c = ctypes.c_void_p()
    try:
        mol_c = molecule.to_SMoleculeC()
        with tempfile.NamedTemporaryFile(dir=os.getcwd(), delete=True) as ctrl:
            ctrl_files.write_ctrl_output(
                    ctrl,
                    rmsfile = "dummy.rms")
            ctrl_files.write_ctrl_selection(
                    ctrl, selection_group, selection_mole_name)
            ctrl_files.write_ctrl_fitting(
                    ctrl, fitting_method, fitting_atom)
            ctrl.write(b"[OPTION]\n")
            ctrl_files.write_value(ctrl, "check_only", check_only)
            ctrl_files.write_value(ctrl, "contact_groups", contact_groups)
            ctrl_files.write_value(ctrl, "ignore_hydrogen", ignore_hydrogen)
            ctrl_files.write_value(ctrl, "two_states", two_states)
            ctrl_files.write_value(ctrl, "avoid_bonding", avoid_bonding)
            ctrl_files.write_value(ctrl, "exclude_residues", exclude_residues)
            ctrl_files.write_value(ctrl, "minimum_distance", minimum_distance)
            ctrl_files.write_value(ctrl, "maximum_distance", maximum_distance)
            ctrl_files.write_value(ctrl, "pbc_correct", pbc_correct)
            ctrl_files.write_value(ctrl, "verbose", verbose)

            ctrl.seek(0)
            LibGenesis().lib.dr_analysis_c(
                    ctypes.byref(mol_c),
                    ctypes.byref(trajs.get_c_obj()),
                    ctypes.byref(ana_period_c),
                    py2c_util.pathlike_to_byte(ctrl.name),
                    ctypes.byref(result_drms_c),
                    )
        n_frame_c = ctypes.c_int(int(trajs.nframe / ana_period))
        result_drms = (c2py_util.conv_double_ndarray(
            result_drms_c, n_frame_c.value)
                    if result_drms_c else None)
        return DrmsAnalysisResult(
                result_drms)
    finally:
        if result_drms_c:
            LibGenesis().lib.deallocate_double(
                    ctypes.byref(result_drms_c),
                    ctypes.byref(n_frame_c))
        if mol_c:
            LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))


MsdAnalysisResult = namedtuple(
        'MsdAnalysisResult',
        ['msd'])


def msd_analysis(
        molecule: SMolecule, trajs :STrajectories,
        ana_period: Optional[int] = 1,
        selection_group: Optional[Iterable[str]] = None,
        selection_mole_name: Optional[Iterable[str]] = None,
        selection: Optional[Iterable[str]] = None,
        mode: Optional[Iterable[int]] = None,
        check_only: Optional[bool] = None,
        oversample: Optional[bool] = None,
        delta: Optional[int] = None,
        ) -> MsdAnalysisResult:
    """
    Executes msd_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
    Returns:
        msd
    """
    mol_c = None
    ana_period_c = ctypes.c_int(ana_period)
    num_analysis_mols_c = ctypes.c_int(0)
    num_delta_c = ctypes.c_int(0)
    result_msd_c = ctypes.c_void_p()
    try:
        mol_c = molecule.to_SMoleculeC()
        with tempfile.NamedTemporaryFile(dir=os.getcwd(), delete=True) as ctrl:
            ctrl_files.write_ctrl_output(
                    ctrl,
                    msdfile = "dummy.msd")
            ctrl_files.write_ctrl_selection(
                    ctrl, selection_group, selection_mole_name)
            ctrl_files.write_ctrl_molecule_selection(
                    ctrl, selection, mode)
            ctrl.write(b"[OPTION]\n")
            ctrl_files.write_kwargs(
                    ctrl,
                    check_only = check_only,
                    oversample = oversample,
                    delta = delta,
                    )

            ctrl.seek(0)
            LibGenesis().lib.ma_analysis_c(
                    ctypes.byref(mol_c),
                    ctypes.byref(trajs.get_c_obj()),
                    ctypes.byref(ana_period_c),
                    py2c_util.pathlike_to_byte(ctrl.name),
                    ctypes.byref(result_msd_c),
                    ctypes.byref(num_analysis_mols_c),
                    ctypes.byref(num_delta_c),
                    )
        result_msd = c2py_util.conv_double_ndarray(
            result_msd_c, [num_delta_c.value, num_analysis_mols_c.value])
        return MsdAnalysisResult(
                result_msd)
    finally:
        if result_msd_c:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(result_msd_c),
                    ctypes.byref(num_delta_c), ctypes.byref(num_analysis_mols_c))
        if mol_c:
            LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))


def hb_analysis(molecule: SMolecule, trajs :STrajectories,
                ana_period: int,
                ctrl_path: str | bytes | os.PathLike
                ) -> str:
    """
    Executes hb_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
        ctrl_path:

    Returnsu:
        TODO
    """
    mol_c = molecule.to_SMoleculeC()
    ana_period_c = ctypes.c_int(ana_period)
    result = ctypes.c_void_p()
    LibGenesis().lib.hb_analysis_c(
            ctypes.byref(mol_c),
            ctypes.byref(trajs.get_c_obj()),
            ctypes.byref(ana_period_c),
            py2c_util.pathlike_to_byte(ctrl_path),
            ctypes.byref(result),
            )
    s = c2py_util.conv_string(result)
    LibGenesis().lib.deallocate_c_string(ctypes.byref(result))
    return s


def diffusion_analysis(msd_data: npt.NDArray[np.float64],
                       ctrl_path: str | bytes | os.PathLike,
                       ) -> npt.NDArray[np.float64]:
    """
    Executes diffusion_analysis.

    Args:
        ctrl_path:

    Returns:
        diffusion
    """
    c_msd = None
    c_out = ctypes.c_void_p(0)
    if msd_data.ndim != 2:
        raise Exception
    d0 = ctypes.c_int(msd_data.shape[0])
    d1 = ctypes.c_int(msd_data.shape[1])
    try:
        c_msd = ctypes.c_void_p(
                LibGenesis().lib.allocate_c_double_array2(
                ctypes.byref(d0),
                ctypes.byref(d1)))
        py2c_util.write_double_ndarray(msd_data, c_msd)
        LibGenesis().lib.diffusion_analysis_c(
                ctypes.byref(c_msd),
                ctypes.byref(d0),
                ctypes.byref(d1),
                py2c_util.pathlike_to_byte(ctrl_path),
                ctypes.byref(c_out),
                )
        return c2py_util.conv_double_ndarray(
                c_out, (msd_data.shape[0], msd_data.shape[1] * 2 - 1))
    finally:
        if c_msd:
            LibGenesis().lib.deallocate_double2(
                    ctypes.byref(c_msd), ctypes.byref(d0), ctypes.byref(d1))
        if c_out:
            ci = ctypes.c_int(msd_data.shape[0] * (msd_data.shape[1] * 2 - 1))
            LibGenesis().lib.deallocate_double(
                    ctypes.byref(c_out), ctypes.byref(ci))


def avecrd_analysis(molecule: SMolecule, trajs :STrajectories,
                ana_period: int,
                ctrl_path: str | bytes | os.PathLike
                ):
    """
    Executes aa_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
        ctrl_path:

    Returns:
        TODO
    """
    mol_c = None
    pdb_ave_c = ctypes.c_void_p()
    try:
        mol_c = molecule.to_SMoleculeC()
        ana_period_c = ctypes.c_int(ana_period)
        LibGenesis().lib.aa_analysis_c(
                ctypes.byref(mol_c),
                ctypes.byref(trajs.get_c_obj()),
                ctypes.byref(ana_period_c),
                py2c_util.pathlike_to_byte(ctrl_path),
                ctypes.byref(pdb_ave_c),
                )
        if pdb_ave_c:
            s = c2py_util.conv_string(pdb_ave_c)
            LibGenesis().lib.deallocate_c_string(ctypes.byref(pdb_ave_c))
            print(s)
        else:
            s = None
    finally:
        if mol_c:
            LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))
        if pdb_ave_c:
            LibGenesis().lib.deallocate_c_string(ctypes.byref(pdb_ave_c))
    return


def wham_analysis(ctrl_path: str | bytes | os.PathLike
                ):
    """
    Executes wham_analysis.

    Args:
        ctrl_path:

    Returns:
        pmf
    """
    result_pmf_c = ctypes.c_void_p(None)
    n_bins = ctypes.c_int(0)
    n_bin_x = ctypes.c_int(0)
    LibGenesis().lib.wa_analysis_c(
            py2c_util.pathlike_to_byte(ctrl_path),
            ctypes.byref(result_pmf_c),
            ctypes.byref(n_bins),
            ctypes.byref(n_bin_x),
            )
    result_pmf = c2py_util.conv_double_ndarray(
            result_pmf_c, [n_bins.value, n_bin_x.value])
    LibGenesis().lib.deallocate_double2(
            ctypes.byref(result_pmf_c),
            ctypes.byref(n_bins), ctypes.byref(n_bin_x))
    return result_pmf


def mbar_analysis(ctrl_path: str | bytes | os.PathLike
                ):
    """
    Executes mbar_analysis.

    Args:
        ctrl_path:

    Returns:
        fene
    """
    result_fene_c = ctypes.c_void_p(None)
    n_replica = ctypes.c_int(0)
    n_blocks = ctypes.c_int(0)
    LibGenesis().lib.mbar_analysis_c(
            py2c_util.pathlike_to_byte(ctrl_path),
            ctypes.byref(result_fene_c),
            ctypes.byref(n_replica),
            ctypes.byref(n_blocks),
            )
    result_fene = c2py_util.conv_double_ndarray(
             result_fene_c, [n_replica.value, n_blocks.value])
    LibGenesis().lib.deallocate_double2(
            ctypes.byref(result_fene_c),
            ctypes.byref(n_replica), ctypes.byref(n_blocks))
    return result_fene


KmeansClusteringResult = namedtuple(
        'KmeansClusteringResult',
        ['pdb_str',
         'cluster_idxs'])


def kmeans_clustering(
        molecule: SMolecule, trajs :STrajectories,
        ana_period: Optional[int] = 1,
        selection_group: Optional[Iterable[str]] = None,
        selection_mole_name: Optional[Iterable[str]] = None,
        fitting_method: Optional[str] = None,
        fitting_atom: Optional[int] = None,
        zrot_ngrid: Optional[int] = None,
        zrot_grid_size: Optional[float] = None,
        mass_weight: Optional[bool] = None,
        check_only: Optional[bool] = None,
        allow_backup: Optional[bool] = None,
        analysis_atom: Optional[int] = None,
        num_clusters: Optional[int] = None,
        max_iteration: Optional[int] = None,
        stop_threshold: Optional[float] = None,
        num_iterations: Optional[int] = None,
        trjout_atom: Optional[int] = None,
        trjout_format: Optional[str] = None,
        trjout_type: Optional[str] = None,
        iseed: Optional[int] = None,
        ) -> KmeansClusteringResult:
    """
    Executes kmeans_clustering.

    Args:
        molecule:
        trajs:
        ana_period:
    Returns:
        (pdb string, cluster indices)
    """
    mol_c = None
    pdb_c = ctypes.c_void_p()
    cluster_idxs_c = ctypes.c_void_p()
    ana_period_c = ctypes.c_int(ana_period)
    cluster_size = ctypes.c_int(0)
    try:
        mol_c = molecule.to_SMoleculeC()
        with tempfile.NamedTemporaryFile(dir=os.getcwd(), delete=True) as ctrl:
            ctrl_files.write_ctrl_output(
                    ctrl,
                    indexfile = "dummy.idx",
                    pdbfile = "dummy_{}.pdb",
                    trjfile = "dummy{}.trj")
            ctrl_files.write_ctrl_selection(
                    ctrl, selection_group, selection_mole_name)
            ctrl_files.write_ctrl_fitting(
                    ctrl, fitting_method, fitting_atom, zrot_ngrid, zrot_grid_size, mass_weight)
            ctrl.write(b"[OPTION]\n")
            ctrl_files.write_kwargs(
                    ctrl,
                    check_only = check_only,
                    allow_backup = allow_backup,
                    analysis_atom = analysis_atom,
                    num_clusters = num_clusters,
                    max_iteration = max_iteration,
                    stop_threshold = stop_threshold,
                    num_iterations = num_iterations,
                    trjout_atom = trjout_atom,
                    trjout_format = trjout_format,
                    trjout_type = trjout_type,
                    iseed = iseed,
                    )

            ctrl.seek(0)
            LibGenesis().lib.kc_analysis_c(
                    ctypes.byref(mol_c),
                    ctypes.byref(trajs.get_c_obj()),
                    ctypes.byref(ana_period_c),
                    py2c_util.pathlike_to_byte(ctrl.name),
                    ctypes.byref(pdb_c),
                    ctypes.byref(cluster_idxs_c),
                    ctypes.byref(cluster_size),
                    )
            if pdb_c:
                pdb = c2py_util.conv_string(pdb_c)
                LibGenesis().lib.deallocate_c_string(ctypes.byref(pdb_c))
            else:
                pdb = None
            cluster_idxs = (c2py_util.conv_int_ndarray(
                    cluster_idxs_c, cluster_size.value)
                            if cluster_idxs_c else None)
            return KmeansClusteringResult(pdb, cluster_idxs)
    finally:
        if cluster_idxs_c:
            LibGenesis().lib.deallocate_int(
                    ctypes.byref(cluster_idxs_c), ctypes.byref(cluster_size))
        if pdb_c:
            LibGenesis().lib.deallocate_c_string(ctypes.byref(pdb_c))
        if mol_c:
            LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))
    return None
