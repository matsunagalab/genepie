from dataclasses import dataclass
from typing import Iterable, Optional, TextIO, Union


def write_ctrl_input(
    dst: TextIO,
    ambcrdfile: Optional[str] = None,
    ambreffile: Optional[str] = None,
    atpfile: Optional[str] = None,
    coorfile: Optional[str] = None,
    cvfile: Optional[str] = None,
    targetfile: Optional[str] = None,
    dcdfile: Optional[str] = None,
    dcdvelfile: Optional[str] = None,
    excfile: Optional[str] = None,
    gprfile: Optional[str] = None,
    groreffile: Optional[str] = None,
    grocrdfile: Optional[str] = None,
    grotopfile: Optional[str] = None,
    mtfile: Optional[str] = None,
    indexfile: Optional[str] = None,
    logfile: Optional[str] = None,
    enefile: Optional[str] = None,
    msdfile: Optional[str] = None,
    pathfile: Optional[str] = None,
    pathcvfile: Optional[str] = None,
    pcafile: Optional[str] = None,
    pdbfile: Optional[str] = None,
    pdb_tgtfile: Optional[str] = None,
    pdb_avefile: Optional[str] = None,
    pdb_aftfile: Optional[str] = None,
    pdb_sphfile: Optional[str] = None,
    pdb_wbxfile: Optional[str] = None,
    prmtopfile: Optional[str] = None,
    psffile: Optional[str] = None,
    radfile: Optional[str] = None,
    refenefile: Optional[str] = None,
    reffile: Optional[str] = None,
    fitfile: Optional[str] = None,
    remfile: Optional[str] = None,
    rstfile: Optional[str] = None,
    rtpfile: Optional[str] = None,
    topfile: Optional[str] = None,
    valfile: Optional[str] = None,
    vecfile: Optional[str] = None,
    velfile: Optional[str] = None,
    weightfile: Optional[str] = None,
    xscfile: Optional[str] = None,
    distfile: Optional[str] = None
) -> None:
    """
    Write genesis control file input information.

    Args:
        file: The file object to which the output will be written.
        All other arguments : File paths to be written. Optional.
    """
    mapping = {
        "Ambcrdfile": ambcrdfile,
        "Ambreffile": ambreffile,
        "Atpfile": atpfile,
        "Coorfile": coorfile,
        "Cvfile": cvfile,
        "Targetfile": targetfile,
        "Dcdfile": dcdfile,
        "Dcdvelfile": dcdvelfile,
        "excfile": excfile,
        "gprfile": gprfile,
        "Groreffile": groreffile,
        "Grocrdfile": grocrdfile,
        "Grotopfile": grotopfile,
        "Mtfile": mtfile,
        "Indexfile": indexfile,
        "Logfile": logfile,
        "Enefile": enefile,
        "msdfile": msdfile,
        "pathfile": pathfile,
        "pathcvfile": pathcvfile,
        "pcafile": pcafile,
        "Pdbfile": pdbfile,
        "Pdb_tgtfile": pdb_tgtfile,
        "Pdb_avefile": pdb_avefile,
        "Pdb_aftfile": pdb_aftfile,
        "Pdb_sphfile": pdb_sphfile,
        "Pdb_wbxfile": pdb_wbxfile,
        "Prmtopfile": prmtopfile,
        "Psffile": psffile,
        "Radfile": radfile,
        "Refenefile": refenefile,
        "Reffile": reffile,
        "Fitfile": fitfile,
        "Remfile": remfile,
        "Rstfile": rstfile,
        "Rtpfile": rtpfile,
        "Topfile": topfile,
        "Valfile": valfile,
        "Vecfile": vecfile,
        "Velfile": velfile,
        "Weightfile": weightfile,
        "Xscfile": xscfile,
        "Distfile": distfile,
    }

    dst.write("b[INPUT]\n")
    for ctrl_name, value in mapping.items():
        if value is not None:
            dst.write(f"{ctrl_name} = {value}\n".encode('utf-8'))


def write_ctrl_output(
    dst: TextIO,
    ambcrdfile: Optional[str] = None,
    angfile: Optional[str] = None,
    cntfile: Optional[str] = None,
    comangfile: Optional[str] = None,
    comdisfile: Optional[str] = None,
    comtorfile: Optional[str] = None,
    coorfile: Optional[str] = None,
    crdfile: Optional[str] = None,
    crsfile: Optional[str] = None,
    disfile: Optional[str] = None,
    enefile: Optional[str] = None,
    excfile: Optional[str] = None,
    fenefile: Optional[str] = None,
    gprfile: Optional[str] = None,
    hb_listfile: Optional[str] = None,
    mapfile: Optional[str] = None,
    msdfile: Optional[str] = None,
    grotopfile: Optional[str] = None,
    grocrdfile: Optional[str] = None,
    grocrd_tgtfile: Optional[str] = None,
    indexfile: Optional[str] = None,
    logfile: Optional[str] = None,
    outfile: Optional[str] = None,
    parfile: Optional[str] = None,
    pathcvfile: Optional[str] = None,
    pcafile: Optional[str] = None,
    pdbfile: Optional[str] = None,
    pdb_avefile: Optional[str] = None,
    pdb_aftfile: Optional[str] = None,
    pmffile: Optional[str] = None,
    pmlfile: Optional[str] = None,
    prjfile: Optional[str] = None,
    probfile: Optional[str] = None,
    qmmm_crdfile: Optional[str] = None,
    qmmm_psffile: Optional[str] = None,
    qmmm_pdbfile: Optional[str] = None,
    qntfile: Optional[str] = None,
    rdffile: Optional[str] = None,
    rmsfile: Optional[str] = None,
    rgfile: Optional[str] = None,
    rstfile: Optional[str] = None,
    txtfile: Optional[str] = None,
    topfile: Optional[str] = None,
    torfile: Optional[str] = None,
    trjfile: Optional[str] = None,
    trrfile: Optional[str] = None,
    valfile: Optional[str] = None,
    vcvfile: Optional[str] = None,
    vecfile: Optional[str] = None,
    velfile: Optional[str] = None,
    voronoifile: Optional[str] = None,
    vmdfile: Optional[str] = None,
    weightfile: Optional[str] = None,
    xscfile: Optional[str] = None,
    tblfile: Optional[str] = None,
    morphfile: Optional[str] = None
) -> None:
    """
    Write genesis control file output information.

    Args:
        file: The file object to which the output will be written.
        All other arguments : File paths to be written. Optional.
    """
    mapping = {
        "Ambcrdfile": ambcrdfile,
        "Angfile": angfile,
        "Cntfile": cntfile,
        "Comangfile": comangfile,
        "Comdisfile": comdisfile,
        "Comtorfile": comtorfile,
        "Coorfile": coorfile,
        "Crdfile": crdfile,
        "Crsfile": crsfile,
        "Disfile": disfile,
        "Enefile": enefile,
        "Excfile": excfile,
        "Fenefile": fenefile,
        "Gprfile": gprfile,
        "HB_listfile": hb_listfile,
        "Mapfile": mapfile,
        "Msdfile": msdfile,
        "Grotopfile": grotopfile,
        "Grocrdfile": grocrdfile,
        "Grocrd_tgtfile": grocrd_tgtfile,
        "Indexfile": indexfile,
        "Logfile": logfile,
        "Outfile": outfile,
        "Parfile": parfile,
        "Pathcvfile": pathcvfile,
        "Pcafile": pcafile,
        "Pdbfile": pdbfile,
        "Pdb_avefile": pdb_avefile,
        "Pdb_aftfile": pdb_aftfile,
        "Pmffile": pmffile,
        "Pmlfile": pmlfile,
        "Prjfile": prjfile,
        "Probfile": probfile,
        "qmmm_crdfile": qmmm_crdfile,
        "qmmm_psffile": qmmm_psffile,
        "qmmm_pdbfile": qmmm_pdbfile,
        "Qntfile": qntfile,
        "Rdffile": rdffile,
        "Rmsfile": rmsfile,
        "Rgfile": rgfile,
        "Rstfile": rstfile,
        "Txtfile": txtfile,
        "Topfile": topfile,
        "Torfile": torfile,
        "Trjfile": trjfile,
        "Trrfile": trrfile,
        "Valfile": valfile,
        "Vcvfile": vcvfile,
        "Vecfile": vecfile,
        "Velfile": velfile,
        "voronoifile": voronoifile,
        "Vmdfile": vmdfile,
        "weightfile": weightfile,
        "Xscfile": xscfile,
        "Tblfile": tblfile,
        "Morphfile": morphfile,
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
    trajectories: Iterable [TrajectoryParameters],
    trj_format: Optional[str] = None,
    trj_type: Optional[str] = None,
    trj_natom: Optional[int] = None
) -> None:
    """
    Write trajectory-related information to a file.

    Args:
        dst: output destination
        trajectories: Set of TrajectoryParameters
        trj_format: Trajectory format
        trj_type: Trajectory type
        trj_natom: Number of atoms in trajectories
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


def write_ctrl_selection(dst: TextIO, group: Optional[Iterable[str]] = None,
                         mole_name: Optional[Iterable[str]] = None) -> None:
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


def write_ctrl_fitting(
    dst: TextIO,
    fitting_method: Optional[str] = None,
    fitting_atom: Optional[int] = None,
    zrot_ngrid: Optional[int] = None,
    zrot_grid_size: Optional[float] = None,
    mass_weight: Optional[bool] = None,
) -> None:
    """
    Write fitting information to a file.

    Args:
        dst: output destination
        fitting_method:
        fitting_atom:
        zrot_ngrid:
        zrot_grid_size:
        mass_weight:
    """
    dst.write(b"[FITTING]\n")
    if fitting_method is not None:
        dst.write(f"fitting_method = {fitting_method}\n".encode('utf-8'))
    if fitting_atom is not None:
        dst.write(f"fitting_atom = {fitting_atom}\n".encode('utf-8'))
    if zrot_ngrid is not None:
        dst.write(f"zrot_ngrid = {zrot_ngrid}\n".encode('utf-8'))
    if zrot_grid_size is not None:
        dst.write(f"zrot_grid_size = {zrot_grid_size:.6D}\n".encode('utf-8'))
    write_yes_no(dst, "mass_weight", mass_weight)


def write_string(dst: TextIO, name: str, v: Optional[str]) -> None:
    if v is not None:
        dst.write(f"{name} = {v}\n".encode('utf-8'))


def write_yes_no(dst: TextIO, name: str, v: Optional[bool]) -> None:
    if v is not None:
        dst.write("{} = {}\n".format(name, "YES" if v else "NO")
                   .encode('utf-8'))


def write_int(dst: TextIO, name: str, v: Optional[int]) -> None:
    if v is not None:
        dst.write(f"{name} = {v}\n".encode('utf-8'))


def write_float(dst: TextIO, name: str, v: Optional[float]) -> None:
    if v is not None:
        dst.write(f"{name} = {v:.12E}\n".encode('utf-8'))


def write_value(dst: TextIO, name: str, v: Optional[Union[bool, int, float, str]]) -> None:
    if v is None:
        pass
    elif isinstance(v, bool):
        write_yes_no(dst, name, v)
    elif isinstance(v, int):
        write_int(dst, name, v)
    elif isinstance(v, float):
        write_float(dst, name, v)
    elif isinstance(v, str):
        write_string(dst, name, v)
    else:
        raise TypeError(f"Unsupported type for value: {type(v)}")


def write_iterable(
        dst: TextIO, name: str,
        vals: Optional[Iterable[Optional[Union[bool, int, float, str]]]]
        ) -> None:
    for idx, v in enumerate(vals, 1):
        write_value(dst, f"{name}{idx}", v)
