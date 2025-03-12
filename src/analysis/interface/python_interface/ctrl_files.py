from dataclasses import dataclass
from typing import TextIO, Optional, List


def write_ctrl_input(
    dst: TextIO,
    amb_crd_file: Optional[str] = None,
    amb_ref_file: Optional[str] = None,
    atp_file: Optional[str] = None,
    coor_file: Optional[str] = None,
    cv_file: Optional[str] = None,
    target_file: Optional[str] = None,
    dcd_file: Optional[str] = None,
    dcd_vel_file: Optional[str] = None,
    exc_file: Optional[str] = None,
    gpr_file: Optional[str] = None,
    gro_ref_file: Optional[str] = None,
    gro_crd_file: Optional[str] = None,
    gro_top_file: Optional[str] = None,
    mt_file: Optional[str] = None,
    index_file: Optional[str] = None,
    log_file: Optional[str] = None,
    ene_file: Optional[str] = None,
    msd_file: Optional[str] = None,
    path_file: Optional[str] = None,
    path_cv_file: Optional[str] = None,
    pca_file: Optional[str] = None,
    pdb_file: Optional[str] = None,
    pdb_tgt_file: Optional[str] = None,
    pdb_ave_file: Optional[str] = None,
    pdb_aft_file: Optional[str] = None,
    pdb_sph_file: Optional[str] = None,
    pdb_wbx_file: Optional[str] = None,
    prm_top_file: Optional[str] = None,
    psf_file: Optional[str] = None,
    rad_file: Optional[str] = None,
    ref_ene_file: Optional[str] = None,
    ref_file: Optional[str] = None,
    fit_file: Optional[str] = None,
    rem_file: Optional[str] = None,
    rst_file: Optional[str] = None,
    rtp_file: Optional[str] = None,
    top_file: Optional[str] = None,
    val_file: Optional[str] = None,
    vec_file: Optional[str] = None,
    vel_file: Optional[str] = None,
    weight_file: Optional[str] = None,
    xsc_file: Optional[str] = None,
    dist_file: Optional[str] = None
) -> None:
    """
    Write genesis control file input information.

    Args:
        file: The file object to which the output will be written.
        All other arguments : File paths to be written. Optional.
    """
    mapping = {
        "Ambcrdfile": amb_crd_file,
        "Ambreffile": amb_ref_file,
        "Atpfile": atp_file,
        "Coorfile": coor_file,
        "Cvfile": cv_file,
        "Targetfile": target_file,
        "Dcdfile": dcd_file,
        "Dcdvelfile": dcd_vel_file,
        "excfile": exc_file,
        "gprfile": gpr_file,
        "Groreffile": gro_ref_file,
        "Grocrdfile": gro_crd_file,
        "Grotopfile": gro_top_file,
        "Mtfile": mt_file,
        "Indexfile": index_file,
        "Logfile": log_file,
        "Enefile": ene_file,
        "msdfile": msd_file,
        "pathfile": path_file,
        "pathcvfile": path_cv_file,
        "pcafile": pca_file,
        "Pdbfile": pdb_file,
        "Pdb_tgtfile": pdb_tgt_file,
        "Pdb_avefile": pdb_ave_file,
        "Pdb_aftfile": pdb_aft_file,
        "Pdb_sphfile": pdb_sph_file,
        "Pdb_wbxfile": pdb_wbx_file,
        "Prmtopfile": prm_top_file,
        "Psffile": psf_file,
        "Radfile": rad_file,
        "Refenefile": ref_ene_file,
        "Reffile": ref_file,
        "Fitfile": fit_file,
        "Remfile": rem_file,
        "Rstfile": rst_file,
        "Rtpfile": rtp_file,
        "Topfile": top_file,
        "Valfile": val_file,
        "Vecfile": vec_file,
        "Velfile": vel_file,
        "Weightfile": weight_file,
        "Xscfile": xsc_file,
        "Distfile": dist_file,
    }

    dst.write("b[INPUT]\n")
    for ctrl_name, value in mapping.items():
        if value is not None:
            dst.write(f"{ctrl_name} = {value}\n".encode('utf-8'))


def write_ctrl_output(
    dst: TextIO,
    amb_crd_file: Optional[str] = None,
    ang_file: Optional[str] = None,
    cnt_file: Optional[str] = None,
    comang_file: Optional[str] = None,
    comdis_file: Optional[str] = None,
    comtor_file: Optional[str] = None,
    coor_file: Optional[str] = None,
    crd_file: Optional[str] = None,
    crs_file: Optional[str] = None,
    dis_file: Optional[str] = None,
    ene_file: Optional[str] = None,
    exc_file: Optional[str] = None,
    fenefile: Optional[str] = None,
    gpr_file: Optional[str] = None,
    hb_list_file: Optional[str] = None,
    map_file: Optional[str] = None,
    msd_file: Optional[str] = None,
    gro_top_file: Optional[str] = None,
    gro_crd_file: Optional[str] = None,
    gro_crd_tgt_file: Optional[str] = None,
    index_file: Optional[str] = None,
    log_file: Optional[str] = None,
    out_file: Optional[str] = None,
    par_file: Optional[str] = None,
    path_cv_file: Optional[str] = None,
    pca_file: Optional[str] = None,
    pdb_file: Optional[str] = None,
    pdb_ave_file: Optional[str] = None,
    pdb_aft_file: Optional[str] = None,
    pmf_file: Optional[str] = None,
    pml_file: Optional[str] = None,
    prj_file: Optional[str] = None,
    prob_file: Optional[str] = None,
    qmmm_crd_file: Optional[str] = None,
    qmmm_psf_file: Optional[str] = None,
    qmmm_pdb_file: Optional[str] = None,
    qnt_file: Optional[str] = None,
    rdf_file: Optional[str] = None,
    rms_file: Optional[str] = None,
    rg_file: Optional[str] = None,
    rst_file: Optional[str] = None,
    txt_file: Optional[str] = None,
    top_file: Optional[str] = None,
    tor_file: Optional[str] = None,
    trj_file: Optional[str] = None,
    trr_file: Optional[str] = None,
    val_file: Optional[str] = None,
    vcv_file: Optional[str] = None,
    vec_file: Optional[str] = None,
    vel_file: Optional[str] = None,
    voronoi_file: Optional[str] = None,
    vmd_file: Optional[str] = None,
    weight_file: Optional[str] = None,
    xsc_file: Optional[str] = None,
    tbl_file: Optional[str] = None,
    morph_file: Optional[str] = None
) -> None:
    """
    Write genesis control file output information.

    Args:
        file: The file object to which the output will be written.
        All other arguments : File paths to be written. Optional.
    """
    mapping = {
        "Ambcrdfile": amb_crd_file,
        "Angfile": ang_file,
        "Cntfile": cnt_file,
        "Comangfile": comang_file,
        "Comdisfile": comdis_file,
        "Comtorfile": comtor_file,
        "Coorfile": coor_file,
        "Crdfile": crd_file,
        "Crsfile": crs_file,
        "Disfile": dis_file,
        "Enefile": ene_file,
        "Excfile": exc_file,
        "Fenefile": fenefile,
        "Gprfile": gpr_file,
        "HB_listfile": hb_list_file,
        "Mapfile": map_file,
        "Msdfile": msd_file,
        "Grotopfile": gro_top_file,
        "Grocrdfile": gro_crd_file,
        "Grocrd_tgtfile": gro_crd_tgt_file,
        "Indexfile": index_file,
        "Logfile": log_file,
        "Outfile": out_file,
        "Parfile": par_file,
        "Pathcvfile": path_cv_file,
        "Pcafile": pca_file,
        "Pdbfile": pdb_file,
        "Pdb_avefile": pdb_ave_file,
        "Pdb_aftfile": pdb_aft_file,
        "Pmffile": pmf_file,
        "Pmlfile": pml_file,
        "Prjfile": prj_file,
        "Probfile": prob_file,
        "qmmm_crdfile": qmmm_crd_file,
        "qmmm_psffile": qmmm_psf_file,
        "qmmm_pdbfile": qmmm_pdb_file,
        "Qntfile": qnt_file,
        "Rdffile": rdf_file,
        "Rmsfile": rms_file,
        "Rgfile": rg_file,
        "Rstfile": rst_file,
        "Txtfile": txt_file,
        "Topfile": top_file,
        "Torfile": tor_file,
        "Trjfile": trj_file,
        "Trrfile": trr_file,
        "Valfile": val_file,
        "Vcvfile": vcv_file,
        "Vecfile": vec_file,
        "Velfile": vel_file,
        "voronoifile": voronoi_file,
        "Vmdfile": vmd_file,
        "weightfile": weight_file,
        "Xscfile": xsc_file,
        "Tblfile": tbl_file,
        "Morphfile": morph_file,
    }

    dst.write(b"[OUTPUT]\n")
    for fortran_key, value in mapping.items():
        if value is not None:
            dst.write(f"{fortran_key} = {value}\n".encode('utf-8'))


@dataclass
class TrajectoryParameters:
    trjfile: Optional[str] = None
    md_step: Optional[int] = None
    mdout_period: Optional[int] = None
    ana_period: Optional[int] = None
    repeat: Optional[int] = None


def write_trajectory_info(
    dst: TextIO,
    trajectories: List[TrajectoryParameters],
    trj_format: Optional[str] = None,
    trj_type: Optional[str] = None,
    trj_natom: Optional[int] = None
) -> None:
    """
    Write trajectory-related information to a file.

    Args:
        dst: 出力先のファイルオブジェクト
        trajectories: TrajectoryParametersのリスト
        trj_format: トラジェクトリフォーマット
        trj_type: トラジェクトリタイプ
        trj_natom: 原子数
    """
    dst.write(b"[TRAJECTORY]\n")
    for idx, traj in enumerate(trajectories, 1):
        if traj.trjfile is not None:
            dst.write(f"trjfile{idx} = {traj.trjfile}\n".encode('utf-8'))
        if traj.md_step is not None:
            dst.write(f"md_step{idx} = {traj.md_step}\n".encode('utf-8'))
        if traj.mdout_period is not None:
            dst.write(f"mdout_period{idx} = {traj.mdout_period}\n".encode('utf-8'))
        if traj.ana_period is not None:
            dst.write(f"ana_period{idx} = {traj.ana_period}\n".encode('utf-8'))
        if traj.repeat is not None:
            dst.write(f"repeat{idx} = {traj.repeat}\n".encode('utf-8'))
    if trj_format is not None:
        dst.write(f"trj_format = {trj_format}\n".encode('utf-8'))
    if trj_type is not None:
        dst.write(f"trj_type = {trj_type}\n".encode('utf-8'))
    if trj_natom is not None:
        dst.write(f"trj_natom = {trj_natom}\n".encode('utf-8'))


def write_ctrl_selection(dst: TextIO, group: Optional[List[str]] = None,
                         mole_name: Optional[List[str]] = None) -> None:
    """
    Write selection information to a file.
    """
    if group is None:
        group = []
    if mole_name is None:
        mole_name = []
    dst.write(b"[SELECTION]\n")
    for i, g in enumerate(group, start=1):
        dst.write(f"group{i} = {g}\n".encode('utf-8'))
    for i, m in enumerate(mole_name, start=1):
        dst.write(f"mole_name{i} = {m}\n".encode('utf-8'))
