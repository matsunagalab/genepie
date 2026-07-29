"""K-means clustering of trajectory frames."""

import ctypes
import io
import os
import tempfile
import numpy as np
import numpy.typing as npt
from typing import (
    Iterable,
    NamedTuple,
    Optional,
)

from ..libgenesis import LibGenesis
from ..s_molecule import SMolecule
from ..s_trajectories import STrajectories
from .. import ctrl_files
from .. import c2py_util
from .._fortran import (
    ctrl_to_bytes,
    fortran_status,
    molecule_c,
)
from ._common import extract_model_blocks


class KmeansClusteringResult(NamedTuple):
    mols_from_pdb: SMolecule
    cluster_idxs: npt.NDArray[np.int64]


def kmeans_clustering(
        molecule: SMolecule, trajs: STrajectories,
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
    pdb_c = ctypes.c_void_p()
    cluster_idxs_c = ctypes.c_void_p()
    ana_period_c = ctypes.c_int(ana_period)
    cluster_size = ctypes.c_int(0)
    try:
        with molecule_c(molecule) as mol_c:
            ctrl = io.BytesIO()
            ctrl_files.write_ctrl_output(
                    ctrl,
                    indexfile="dummy.idx",
                    pdbfile="dummy_{}.pdb",
                    trjfile="dummy{}.trj")
            ctrl_files.write_ctrl_selection(
                    ctrl, selection_group, selection_mole_name)
            ctrl_files.write_ctrl_fitting(
                    ctrl, fitting_method, fitting_atom,
                    zrot_ngrid, zrot_grid_size, mass_weight)
            ctrl.write(b"[OPTION]\n")
            ctrl_files.write_kwargs(
                    ctrl,
                    check_only=check_only,
                    allow_backup=allow_backup,
                    analysis_atom=analysis_atom,
                    num_clusters=num_clusters,
                    max_iteration=max_iteration,
                    stop_threshold=stop_threshold,
                    num_iterations=num_iterations,
                    trjout_atom=trjout_atom,
                    trjout_format=trjout_format,
                    trjout_type=trjout_type,
                    iseed=iseed,
                    )

            ctrl_bytes, ctrl_len = ctrl_to_bytes(ctrl)
            with fortran_status() as (status, msgbuf, msglen):
                LibGenesis().lib.kc_analysis_c(
                        ctypes.byref(mol_c),
                        ctypes.byref(trajs.get_c_obj()),
                        ctypes.byref(ana_period_c),
                        ctrl_bytes,
                        ctrl_len,
                        ctypes.byref(pdb_c),
                        ctypes.byref(cluster_idxs_c),
                        ctypes.byref(cluster_size),
                        ctypes.byref(status),
                        msgbuf,
                        ctypes.c_int(msglen),
                        )

        if pdb_c:
            pdb_str = c2py_util.conv_string(pdb_c)
            pdb_mols = []
            for pdb_block in extract_model_blocks(pdb_str):
                with tempfile.NamedTemporaryFile(
                        dir=os.getcwd(), delete=True) as pdb_file:
                    pdb_file.write(pdb_block.encode())
                    pdb_file.seek(0)
                    pdb_mols.append(SMolecule.from_file(pdb=pdb_file.name))
        else:
            pdb_mols = None
        cluster_idxs = (c2py_util.conv_int_ndarray(
                cluster_idxs_c, cluster_size.value)
                        if cluster_idxs_c else None)
        return KmeansClusteringResult(pdb_mols, cluster_idxs)
    finally:
        if cluster_idxs_c:
            LibGenesis().lib.deallocate_int(
                    ctypes.byref(cluster_idxs_c), ctypes.byref(cluster_size))
        if pdb_c:
            LibGenesis().lib.deallocate_c_string(ctypes.byref(pdb_c))
