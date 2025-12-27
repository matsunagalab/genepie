# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------
import os
import pathlib
from ..ctrl_files import TrajectoryParameters
from ..s_molecule import SMolecule
from .. import genesis_exe


def test_msd_analysis():
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")

    mol = SMolecule.from_file(pdb=pdb_path, psf=psf_path)
    trajs, subset_mol =  genesis_exe.crd_convert(
            mol,
            traj_params = [
                TrajectoryParameters(
                    trjfile = "BPTI_run.dcd",
                    md_step = 10000,
                    mdout_period = 1000,
                    ana_period = 1,
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
            d = genesis_exe.msd_analysis(
                    mol, t,
                    selection_group = ["rnam:TIP3", ],
                    selection = [1, ],
                    mode = ["ALL", ],
                    oversample = True,
                    delta = 9,
                    )
            # Validate MSD results
            assert d.msd is not None, "MSD result should not be None"
            assert d.msd.shape[0] > 0, "MSD result should have at least one time point"
            assert d.msd.shape[1] > 0, "MSD result should have at least one molecule"
            # MSD should be non-negative
            assert (d.msd >= 0).all(), "MSD values should be non-negative"
            print(f"MSD shape: {d.msd.shape}, min={d.msd.min():.3f}, max={d.msd.max():.3f}", flush=True)
    finally:
        if hasattr(trajs, "close"):
            trajs.close()


def main():
    if os.path.exists("dummy.trj"):
        os.remove("dummy.trj")
    test_msd_analysis()


if __name__ == "__main__":
    main()
