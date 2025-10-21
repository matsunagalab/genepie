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


def test_crd():
    # 関数を呼び出す
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
            selection_group = ["an:CA", ],
            fitting_method = "TR+ROT",
            fitting_atom = 1,
            check_only = False,
            centering = True,
            centering_atom = 1,
            center_coord = (0.0, 0.0, 0.0),
            pbc_correct = "molecule",
            rename_res = ["HSE HIS", "HSE HIS"],
        )

    _ = subset_mol

    try:
        for t in trajs:
           pass

    finally:
        if hasattr(trajs, "close"):
            trajs.close()


def main():
    if os.path.exists("dummy.trj"):
        os.remove("dummy.trj")
    test_crd()


if __name__ == "__main__":
    main()
