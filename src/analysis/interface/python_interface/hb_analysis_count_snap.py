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


def test_hb_analysis_Count_atom():
    pdb_path = pathlib.Path("RALP_DPPC_run.pdb")
    psf_path = pathlib.Path("RALP_DPPC.psf")

    mol = SMolecule.from_file(pdb=pdb_path, psf=psf_path)
    trajs, subset_mol =  genesis_exe.crd_convert(
            mol,
            traj_params = [
                TrajectoryParameters(
                    trjfile = "RALP_DPPC_run.dcd",
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
            centering      = True,
            centering_atom = 1,
            center_coord   = (0.0,0.0,0.0),
            pbc_correct = "NO",
            rename_res = ["HSE HIS","HSD HIS",],
    ) 

    _ = subset_mol

    try: 

        for t in trajs:
            d = genesis_exe.hb_analysis(
                    mol, t,
                    selection_group = ["sid:PROA",
                                       "resname:DPPC & (an:O11 | an:O12 | an:O13 | an:O14)", ],
                    check_only = False,
                    output_type = "Count_Snap",
                    solvent_list  = "DPPC",
                    analysis_atom = 1,
                    target_atom   = 2,
                    boundary_type = "PBC",
                    hb_distance   = 3.4,
                    dha_angle     = 120.0,
                    hda_angle     = 30.0,
                    )
            print(d, flush=True)
    finally:
        if hasattr(trajs, "close"):
            trajs.close()


def main():
    if os.path.exists("dummy.trj"):
        os.remove("dummy.trj")
    test_hb_analysis_Count_atom()


if __name__ == "__main__":
    main()
