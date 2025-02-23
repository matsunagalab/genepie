import ctypes
import os
import numpy as np
import numpy.typing as npt
from libgenesis import LibGenesis
from s_molecule import SMolecule, py2c_s_molecule
from s_trajectories import STrajectories, STrajectoriesArray
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
    mol_c = py2c_s_molecule(molecule)
    LibGenesis().lib.crd_convert_c(
            ctypes.byref(mol_c),
            py2c_util.pathlike_to_byte(ctrl_filename),
            ctypes.byref(buf),
            ctypes.byref(num_trajs_c))
    LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))
    return STrajectoriesArray(buf, num_trajs_c.value)


def trj_analysis(molecule: SMolecule, trajs :STrajectories,
                 ana_period: int,
                 ctrl_path: str | bytes | os.PathLike
                 ) -> tuple[npt.NDArray[np.float64],
                            npt.NDArray[np.float64],
                            npt.NDArray[np.float64],
                            npt.NDArray[np.float64],
                            npt.NDArray[np.float64],
                            npt.NDArray[np.float64]]:
    """
    Executes trj_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
        ctrl_path:

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
        mol_c = py2c_s_molecule(molecule)
        LibGenesis().lib.trj_analysis_c(
                ctypes.byref(mol_c),
                ctypes.byref(trajs.get_c_obj()),
                ctypes.byref(ana_period_c),
                py2c_util.pathlike_to_byte(ctrl_path),
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
        return (result_distance, result_angle, result_torsion,
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


def rg_analysis(molecule: SMolecule, trajs :STrajectories,
                ana_period: int,
                ctrl_path: str | bytes | os.PathLike
                ) -> npt.NDArray[np.float64]:
    """
    Executes rg_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
        ctrl_path:

    Returns:
        rg
    """
    mol_c = py2c_s_molecule(molecule)

    ana_period_c = ctypes.c_int(ana_period)
    result_rg_c = ctypes.c_void_p(None)

    LibGenesis().lib.rg_analysis_c(
            ctypes.byref(mol_c),
            ctypes.byref(trajs.get_c_obj()),
            ctypes.byref(ana_period_c),
            py2c_util.pathlike_to_byte(ctrl_path),
            ctypes.byref(result_rg_c),
            )

    n_frame_c = ctypes.c_int(int(trajs.nframe / ana_period))
    result_rg = c2py_util.conv_double_ndarray(
            result_rg_c, n_frame_c.value)
    LibGenesis().lib.deallocate_double(
            ctypes.byref(result_rg_c),
            ctypes.byref(n_frame_c))
    return result_rg


def rmsd_analysis(molecule: SMolecule, trajs :STrajectories,
                ana_period: int,
                ctrl_path: str | bytes | os.PathLike
                ) -> npt.NDArray[np.float64]:
    """
    Executes rmsd_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
        ctrl_path:

    Returns:
        rmsd
    """
    mol_c = py2c_s_molecule(molecule)

    ana_period_c = ctypes.c_int(ana_period)
    result_ra_c = ctypes.c_void_p(None)

    LibGenesis().lib.ra_analysis_c(
            ctypes.byref(mol_c),
            ctypes.byref(trajs.get_c_obj()),
            ctypes.byref(ana_period_c),
            py2c_util.pathlike_to_byte(ctrl_path),
            ctypes.byref(result_ra_c),
            )

    n_frame_c = ctypes.c_int(int(trajs.nframe / ana_period))
    result_ra = c2py_util.conv_double_ndarray(
            result_ra_c, n_frame_c.value)
    LibGenesis().lib.deallocate_double(
            ctypes.byref(result_ra_c),
            ctypes.byref(n_frame_c))
    return result_ra


def drms_analysis(molecule: SMolecule, trajs :STrajectories,
                ana_period: int,
                ctrl_path: str | bytes | os.PathLike
                ) -> npt.NDArray[np.float64]:
    """
    Executes drms_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
        ctrl_path:

    Returns:
        drms
    """
    mol_c = py2c_s_molecule(molecule)

    ana_period_c = ctypes.c_int(ana_period)
    result_dr_c = ctypes.c_void_p(None)

    LibGenesis().lib.dr_analysis_c(
            ctypes.byref(mol_c),
            ctypes.byref(trajs.get_c_obj()),
            ctypes.byref(ana_period_c),
            py2c_util.pathlike_to_byte(ctrl_path),
            ctypes.byref(result_dr_c),
            )

    n_frame_c = ctypes.c_int(int(trajs.nframe / ana_period))
    result_dr = c2py_util.conv_double_ndarray(
            result_dr_c, n_frame_c.value)
    LibGenesis().lib.deallocate_double(
            ctypes.byref(result_dr_c),
            ctypes.byref(n_frame_c))
    return result_dr


def msd_analysis(molecule: SMolecule, trajs :STrajectories,
                ana_period: int,
                ctrl_path: str | bytes | os.PathLike
                ) -> tuple(npt.NDArray[np.float64]):
    """
    Executes msd_analysis.

    Args:
        molecule:
        trajs:
        ana_period:
        ctrl_path:

    Returns:
        msd
    """
    mol_c = py2c_s_molecule(molecule)

    ana_period_c = ctypes.c_int(ana_period)
    num_analysis_mols_c = ctypes.c_int(0)
    num_delta_c = ctypes.c_int(0)
    result_msd_c = ctypes.c_void_p(None)

    LibGenesis().lib.ma_analysis_c(
            ctypes.byref(mol_c),
            ctypes.byref(trajs.get_c_obj()),
            ctypes.byref(ana_period_c),
            py2c_util.pathlike_to_byte(ctrl_path),
            ctypes.byref(result_msd_c),
            ctypes.byref(num_analysis_mols_c),
            ctypes.byref(num_delta_c),
            )

    result_msd = c2py_util.conv_double_ndarray(
            result_msd_c, [num_delta_c.value, num_analysis_mols_c.value])
    LibGenesis().lib.deallocate_double2(
            ctypes.byref(result_msd_c),
            ctypes.byref(num_delta_c), ctypes.byref(num_analysis_mols_c))
    return result_msd


def hb_analysis(molecule: SMolecule, trajs :STrajectories,
                ana_period: int,
                ctrl_path: str | bytes | os.PathLike
                ):
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
    mol_c = py2c_s_molecule(molecule)
    ana_period_c = ctypes.c_int(ana_period)
    LibGenesis().lib.hb_analysis_c(
            ctypes.byref(mol_c),
            ctypes.byref(trajs.get_c_obj()),
            ctypes.byref(ana_period_c),
            py2c_util.pathlike_to_byte(ctrl_path),
            )
    return


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
    mol_c = py2c_s_molecule(molecule)
    ana_period_c = ctypes.c_int(ana_period)
    LibGenesis().lib.aa_analysis_c(
            ctypes.byref(mol_c),
            ctypes.byref(trajs.get_c_obj()),
            ctypes.byref(ana_period_c),
            py2c_util.pathlike_to_byte(ctrl_path),
            )
    return


def wham_analysis(n_bins: int, n_bin_x: int,
                ctrl_path: str | bytes | os.PathLike
                ):
    """
    Executes wham_analysis.

    Args:
        n_bins:
        n_bin_x
        ctrl_path:

    Returns:
        pmf
    """
    result_pmf_c = ctypes.c_void_p(None)
    n_bins_c = ctypes.c_int(n_bins)
    n_bin_x_c = ctypes.c_int(n_bin_x)
    LibGenesis().lib.wa_analysis_c(
            py2c_util.pathlike_to_byte(ctrl_path),
            ctypes.byref(result_pmf_c),
            )
    result_pmf = c2py_util.conv_double_ndarray(
            result_pmf_c, [n_bins_c.value, n_bin_x_c.value])
    LibGenesis().lib.deallocate_double2(
            ctypes.byref(result_pmf_c),
            ctypes.byref(n_bins_c), ctypes.byref(n_bin_x_c))
    return result_pmf


def mbar_analysis(n_replica: int, n_blocks: int,
                ctrl_path: str | bytes | os.PathLike
                ):
    """
    Executes mbar_analysis.

    Args:
        n_replica:
        n_blocks:
        ctrl_path:

    Returns:
        fene
    """
    result_fene_c = ctypes.c_void_p(None)
    n_replica_c = ctypes.c_int(n_replica)
    n_blocks_c = ctypes.c_int(n_blocks)
    LibGenesis().lib.mbar_analysis_c(
            py2c_util.pathlike_to_byte(ctrl_path),
            ctypes.byref(result_fene_c),
            )
    result_fene = c2py_util.conv_double_ndarray(
             result_fene_c, [n_replica_c.value, n_blocks_c.value])
    LibGenesis().lib.deallocate_double2(
            ctypes.byref(result_fene_c),
            ctypes.byref(n_replica_c), ctypes.byref(n_blocks_c))
    return result_fene


def kmeans_clustering(molecule: SMolecule, trajs :STrajectories,
                ana_period: int,
                ctrl_path: str | bytes | os.PathLike
                ):
    """
    Executes kmeans_clustering.

    Args:
        molecule:
        trajs:
        ana_period:
        ctrl_path:

    Returns:
        TODO
    """
    mol_c = py2c_s_molecule(molecule)
    ana_period_c = ctypes.c_int(ana_period)
    LibGenesis().lib.kc_analysis_c(
            ctypes.byref(mol_c),
            ctypes.byref(trajs.get_c_obj()),
            ctypes.byref(ana_period_c),
            py2c_util.pathlike_to_byte(ctrl_path),
            )
    return
