# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------
import os
from .conftest import BPTI_PDB, BPTI_PSF, BPTI_DCD
from ..ctrl_files import TrajectoryParameters
from ..s_molecule import SMolecule
from .. import genesis_exe


def test_rg_analysis():
    """Test RG analysis with zerocopy interface."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF)
    trajs, subset_mol = genesis_exe.crd_convert(
            mol,
            traj_params=[
                TrajectoryParameters(
                    trjfile=str(BPTI_DCD),
                    md_step=10,
                    mdout_period=1,
                    ana_period=1,
                    repeat=1,
                ),
            ],
            trj_format="DCD",
            trj_type="COOR+BOX",
            trj_natom=0,
            selection_group=["all"],
            fitting_method="NO",
            fitting_atom=1,
            check_only=False,
            pbc_correct="NO",
    )
    _ = subset_mol

    try:
        for t in trajs:
            result = genesis_exe.rg_analysis(
                mol, t,
                analysis_selection="all",
                ana_period=1,
                mass_weighted=True,
            )

            # Validate Rg results
            assert result.rg is not None, "Rg result should not be None"
            assert len(result.rg) > 0, "Rg result should have at least one frame"
            assert all(r > 0 for r in result.rg), "Rg values should be positive"
            # Rg should be in reasonable range for proteins (5-100 Angstroms)
            assert all(r < 100.0 for r in result.rg), "Rg values should be reasonable (< 100 A)"

            print(f"Rg values (n={len(result.rg)}): min={min(result.rg):.3f}, max={max(result.rg):.3f}")
    finally:
        if hasattr(trajs, "close"):
            trajs.close()


def main():
    if os.path.exists("dummy.trj"):
        os.remove("dummy.trj")
    try:
        test_rg_analysis()
        print("\nAll RG tests passed!")
    except Exception as e:
        print(f"\nRG test FAILED: {e}")
        raise


if __name__ == "__main__":
    main()
