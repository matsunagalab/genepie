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


def test_rmsd_analysis():
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
            d = genesis_exe.rmsd_analysis(
                    mol, t,
                    selection_group = ["sid:BPTI and an:CA", ],
                    fitting_method = "TR+ROT",
                    fitting_atom = 1,
                    check_only = False,
                    analysis_atom  = 1,
                    )
            # Validate RMSD results
            assert d.rmsd is not None, "RMSD result should not be None"
            assert len(d.rmsd) > 0, "RMSD result should have at least one frame"
            assert all(r >= 0 for r in d.rmsd), "RMSD values should be non-negative"
            # RMSD should be in reasonable range (0-50 Angstroms for proteins)
            assert all(r < 50.0 for r in d.rmsd), "RMSD values should be reasonable (< 50 Å)"
            print(f"RMSD values (n={len(d.rmsd)}): min={min(d.rmsd):.3f}, max={max(d.rmsd):.3f}", flush=True)

    finally:
        if hasattr(trajs, "close"):
            trajs.close()


def main():
    if os.path.exists("dummy.trj"):
        os.remove("dummy.trj")
    test_rmsd_analysis()


if __name__ == "__main__":
    main()
