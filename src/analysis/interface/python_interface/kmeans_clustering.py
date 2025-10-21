# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent))
    __package__ = "python_interface"
# --------------------------------------------
import os
import pathlib
from .ctrl_files import TrajectoryParameters
from .s_molecule import SMolecule
from . import genesis_exe


def test_kmeans_clustering():
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")

    mol = SMolecule.from_file(pdb=pdb_path, psf=psf_path)
    trajs, subset_mol =  genesis_exe.crd_convert(
            mol,
            traj_params = [
                TrajectoryParameters(
                    trjfile = "BPTI_run.dcd",
                    md_step = 10,
                    mdout_period = 1,
                    ana_period = 1,
                    repeat = 1,
                    ),
                ],
            trj_format = "DCD",
            trj_type = "COOR+BOX",
            trj_natom = 0,
            selection_group = ["all", ],
            fitting_method = "NO",
            fitting_atom = 1,
            check_only = False,
            pbc_correct = "NO",
    )

    _ = subset_mol

    try: 
        for t in trajs:
            ret = genesis_exe.kmeans_clustering(
                    mol, t,
                    selection_group = ["an:CA", ],
                    fitting_method = "TR+ROT",
                    fitting_atom = 1,
                    check_only = False,
                    allow_backup    = False,
                    analysis_atom   = 1,
                    num_clusters    = 2,
                    max_iteration   = 100,
                    stop_threshold  = 98.0,
                    num_iterations  = 5,
                    trjout_atom     = 1,
                    trjout_format   = "DCD",
                    trjout_type     = "COOR",
                    iseed           = 3141592,
                    )
            for mol in ret.mols_from_pdb:
                print("num_atoms = ", mol.num_atoms)
            print(ret.cluster_idxs)
    finally:
        if hasattr(trajs, "close"):
            trajs.close()
    

def main():
    if os.path.exists("dummy.trj"):
        os.remove("dummy.trj")
    test_kmeans_clustering()


if __name__ == "__main__":
    main()
