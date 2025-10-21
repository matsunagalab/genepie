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


def test_drms_analysis():
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")
    ref_path = pathlib.Path("BPTI_ionize.pdb")

    mol = SMolecule.from_file(pdb=pdb_path, psf=psf_path, ref=ref_path)
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
            drms, = genesis_exe.drms_analysis(
                    mol, t,
                    selection_group = ["an: CA", ],
                    check_only = False,
                    contact_groups = 1,
                    ignore_hydrogen  = False,
                    two_states       = False,
                    avoid_bonding    = True,
                    exclude_residues = 4,
                    minimum_distance = 1.0,
                    maximum_distance = 6.0,
                    pbc_correct      = False,
                    verbose          = True,
                    )
            print(drms)
    finally:
        if hasattr(trajs, "close"):
            trajs.close()


def main():
    if os.path.exists("dummy.trj"):
        os.remove("dummy.trj")
    test_drms_analysis()


if __name__ == "__main__":
    main()
