from typing import Any, Iterable, Optional, TextIO


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
    parfile: Optional[str] = None,
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
    strfile: Optional[str] = None,
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
        "Parfile": parfile,
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
        "Strfile": strfile,
        "Topfile": topfile,
        "Valfile": valfile,
        "Vecfile": vecfile,
        "Velfile": velfile,
        "Weightfile": weightfile,
        "Xscfile": xscfile,
        "Distfile": distfile,
    }

    dst.write(b"[INPUT]\n")
    for ctrl_name, value in mapping.items():
        if value is not None:
            dst.write(f"{ctrl_name} = {value}\n".encode('utf-8'))


def write_ctrl_output(
    dst: TextIO,
    ambcrdfile: Optional[str] = None,
    angfile: Optional[str] = None,
    cntfile: Optional[str] = None,
    dcdfile: Optional[str] = None,
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
        "Dcdfile": dcdfile,
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


def write_ctrl_selection(dst: TextIO, group: Optional[Iterable[str]] = None,
                         mole_name: Optional[Iterable[str]] = None) -> None:
    """
    Write selection information to a file.
    """
    dst.write(b"[SELECTION]\n")
    write_kwargs(dst,
                 group = NumberingData(group),
                 mole_name = NumberingData(mole_name),
                 )


def write_ctrl_molecule_selection(dst: TextIO, selection: Optional[Iterable[str]] = None,
                                  mode: Optional[Iterable[int]] = None) -> None:
    """
    Write selection information to a file.
    """
    dst.write(b"[MOLECULE_SELECTION]\n")
    write_kwargs(dst,
                 selection = NumberingData(selection),
                 mode = NumberingData(mode),
                 )


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
    write_kwargs(dst,
                 fitting_method = fitting_method,
                 fitting_atom = fitting_atom,
                 zrot_ngrid = zrot_ngrid,
                 zrot_grid_size = zrot_grid_size,
                 mass_weight = mass_weight,
                 )


def bool_to_yes_no(v: bool) -> str:
    return "YES" if v else "NO"


def float_to_str(v: float) -> str:
    return f"{v:.6E}"


def iterable_to_str(vals: Iterable[Any]) -> str:
    return ' '.join(map(value_to_str_auto, vals))



def value_to_str_auto(v: Any) -> Optional[str]:
    if isinstance(v, bool):
        return bool_to_yes_no(v)
    elif isinstance(v, float):
        return float_to_str(v)
    elif isinstance(v, Iterable) and not isinstance(v, (str, bytes)):
        return iterable_to_str(v)
    else:
        return str(v)


def write_value(dst: TextIO, name: str,
                v: Optional[Any]):
    if v is not None:
        dst.write("{} = {}\n".format(name, value_to_str_auto(v))
                   .encode('utf-8'))


def write_numbering_values(
        dst: TextIO, name: str,
        vals: Optional[Iterable[Any]]
        ) -> None:
    """
    name1 = vals[0]
    name2 = vals[1]
    name3 = vals[2]
    """
    if vals is not None:
        for idx, v in enumerate(vals, 1):
            write_value(dst, f"{name}{idx}", v)


class NumberingData:
    def __init__(self, src: Iterable[Any]):
        self.src = src


def write_values(
        dst: TextIO,
        values: Iterable[tuple[str, Optional[Any]]]) -> None:
    for name, v in values:
        if isinstance(v, NumberingData):
            write_numbering_values(dst, name, v.src)
        else:
            write_value(dst, name, v)


def write_kwargs(dst: TextIO, **kwargs) -> None:
    write_values(dst, kwargs.items())


# ============================================================================
# ATDYN Control Sections
# ============================================================================

def write_ctrl_energy(
    dst: TextIO,
    forcefield: Optional[str] = None,
    electrostatic: Optional[str] = None,
    switchdist: Optional[float] = None,
    cutoffdist: Optional[float] = None,
    pairlistdist: Optional[float] = None,
    dielec_const: Optional[float] = None,
    pme_alpha: Optional[float] = None,
    pme_ngrid_x: Optional[int] = None,
    pme_ngrid_y: Optional[int] = None,
    pme_ngrid_z: Optional[int] = None,
    pme_nspline: Optional[int] = None,
    vdw_force_switch: Optional[bool] = None,
    vdw_shift: Optional[bool] = None,
    implicit_solvent: Optional[str] = None,
    gbsa_eps_solvent: Optional[float] = None,
    gbsa_eps_solute: Optional[float] = None,
    gbsa_salt_cons: Optional[float] = None,
    gbsa_surf_tens: Optional[float] = None,
    table_density: Optional[float] = None,
    output_style: Optional[str] = None,
    dispersion_corr: Optional[str] = None,
    contact_check: Optional[bool] = None,
    vacuum: Optional[bool] = None,
) -> None:
    """Write [ENERGY] section for atdyn."""
    dst.write(b"[ENERGY]\n")
    write_kwargs(dst,
                 forcefield=forcefield,
                 electrostatic=electrostatic,
                 switchdist=switchdist,
                 cutoffdist=cutoffdist,
                 pairlistdist=pairlistdist,
                 dielec_const=dielec_const,
                 pme_alpha=pme_alpha,
                 pme_ngrid_x=pme_ngrid_x,
                 pme_ngrid_y=pme_ngrid_y,
                 pme_ngrid_z=pme_ngrid_z,
                 pme_nspline=pme_nspline,
                 vdw_force_switch=vdw_force_switch,
                 vdw_shift=vdw_shift,
                 implicit_solvent=implicit_solvent,
                 gbsa_eps_solvent=gbsa_eps_solvent,
                 gbsa_eps_solute=gbsa_eps_solute,
                 gbsa_salt_cons=gbsa_salt_cons,
                 gbsa_surf_tens=gbsa_surf_tens,
                 table_density=table_density,
                 output_style=output_style,
                 dispersion_corr=dispersion_corr,
                 contact_check=contact_check,
                 vacuum=vacuum,
                 )


def write_ctrl_dynamics(
    dst: TextIO,
    integrator: Optional[str] = None,
    nsteps: Optional[int] = None,
    timestep: Optional[float] = None,
    eneout_period: Optional[int] = None,
    crdout_period: Optional[int] = None,
    velout_period: Optional[int] = None,
    rstout_period: Optional[int] = None,
    stoptr_period: Optional[int] = None,
    nbupdate_period: Optional[int] = None,
    iseed: Optional[int] = None,
    initial_time: Optional[float] = None,
    annealing: Optional[bool] = None,
    anneal_period: Optional[int] = None,
    dtemperature: Optional[float] = None,
    verbose: Optional[bool] = None,
    target_md: Optional[bool] = None,
    steered_md: Optional[bool] = None,
) -> None:
    """Write [DYNAMICS] section for atdyn."""
    dst.write(b"[DYNAMICS]\n")
    write_kwargs(dst,
                 integrator=integrator,
                 nsteps=nsteps,
                 timestep=timestep,
                 eneout_period=eneout_period,
                 crdout_period=crdout_period,
                 velout_period=velout_period,
                 rstout_period=rstout_period,
                 stoptr_period=stoptr_period,
                 nbupdate_period=nbupdate_period,
                 iseed=iseed,
                 initial_time=initial_time,
                 annealing=annealing,
                 anneal_period=anneal_period,
                 dtemperature=dtemperature,
                 verbose=verbose,
                 target_md=target_md,
                 steered_md=steered_md,
                 )


def write_ctrl_minimize(
    dst: TextIO,
    method: Optional[str] = None,
    nsteps: Optional[int] = None,
    eneout_period: Optional[int] = None,
    crdout_period: Optional[int] = None,
    rstout_period: Optional[int] = None,
    nbupdate_period: Optional[int] = None,
    force_scale_init: Optional[float] = None,
    force_scale_max: Optional[float] = None,
    verbose: Optional[bool] = None,
    tol_rmsg: Optional[float] = None,
    tol_maxg: Optional[float] = None,
) -> None:
    """Write [MINIMIZE] section for atdyn."""
    dst.write(b"[MINIMIZE]\n")
    write_kwargs(dst,
                 method=method,
                 nsteps=nsteps,
                 eneout_period=eneout_period,
                 crdout_period=crdout_period,
                 rstout_period=rstout_period,
                 nbupdate_period=nbupdate_period,
                 force_scale_init=force_scale_init,
                 force_scale_max=force_scale_max,
                 verbose=verbose,
                 tol_rmsg=tol_rmsg,
                 tol_maxg=tol_maxg,
                 )


def write_ctrl_boundary(
    dst: TextIO,
    type: Optional[str] = None,
    box_size_x: Optional[float] = None,
    box_size_y: Optional[float] = None,
    box_size_z: Optional[float] = None,
    domain_x: Optional[int] = None,
    domain_y: Optional[int] = None,
    domain_z: Optional[int] = None,
) -> None:
    """Write [BOUNDARY] section for atdyn."""
    dst.write(b"[BOUNDARY]\n")
    write_kwargs(dst,
                 type=type,
                 box_size_x=box_size_x,
                 box_size_y=box_size_y,
                 box_size_z=box_size_z,
                 domain_x=domain_x,
                 domain_y=domain_y,
                 domain_z=domain_z,
                 )


def write_ctrl_ensemble(
    dst: TextIO,
    ensemble: Optional[str] = None,
    tpcontrol: Optional[str] = None,
    temperature: Optional[float] = None,
    pressure: Optional[float] = None,
    gamma_t: Optional[float] = None,
    gamma_p: Optional[float] = None,
    isotropy: Optional[str] = None,
) -> None:
    """Write [ENSEMBLE] section for atdyn."""
    dst.write(b"[ENSEMBLE]\n")
    write_kwargs(dst,
                 ensemble=ensemble,
                 tpcontrol=tpcontrol,
                 temperature=temperature,
                 pressure=pressure,
                 gamma_t=gamma_t,
                 gamma_p=gamma_p,
                 isotropy=isotropy,
                 )


def write_ctrl_constraints(
    dst: TextIO,
    rigid_bond: Optional[bool] = None,
    shake_iteration: Optional[int] = None,
    shake_tolerance: Optional[float] = None,
    water_model: Optional[str] = None,
    hydrogen_type: Optional[str] = None,
) -> None:
    """Write [CONSTRAINTS] section for atdyn."""
    dst.write(b"[CONSTRAINTS]\n")
    write_kwargs(dst,
                 rigid_bond=rigid_bond,
                 shake_iteration=shake_iteration,
                 shake_tolerance=shake_tolerance,
                 water_model=water_model,
                 hydrogen_type=hydrogen_type,
                 )
