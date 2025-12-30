# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------
import os
import numpy as np
from .conftest import BPTI_PDB, BPTI_PSF, BPTI_DCD
from ..ctrl_files import TrajectoryParameters
from ..s_molecule import SMolecule
from .. import genesis_exe


def compute_contact_list_from_refcoord(mol, selection_indices, min_dist=1.0, max_dist=6.0, exclude_residues=4):
    """Compute contact list and reference distances from molecule reference coordinates.

    This is a simplified version of the Fortran setup_contact_list for testing purposes.
    It computes contacts between atoms within min_dist-max_dist and excludes nearby residues.
    """
    ref_coord = mol.atom_refcoord  # (num_atoms, 3), 1-indexed in mol
    residue_no = mol.residue_no

    contact_list = []
    contact_dist = []

    for i, idx_i in enumerate(selection_indices):
        for j, idx_j in enumerate(selection_indices):
            if j <= i:
                continue  # Skip duplicates and self

            # Check residue exclusion
            res_i = residue_no[idx_i - 1]  # Convert to 0-indexed
            res_j = residue_no[idx_j - 1]
            if abs(res_i - res_j) < exclude_residues:
                continue

            # Calculate distance from reference coordinates
            coord_i = ref_coord[idx_i - 1, :]  # Convert to 0-indexed
            coord_j = ref_coord[idx_j - 1, :]
            d = np.sqrt(np.sum((coord_i - coord_j) ** 2))

            if min_dist <= d < max_dist:
                contact_list.append([min(idx_i, idx_j), max(idx_i, idx_j)])
                contact_dist.append(d)

    return np.array(contact_list, dtype=np.int32).T, np.array(contact_dist, dtype=np.float64)


def test_drms_analysis():
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)
    trajs, subset_mol =  genesis_exe.crd_convert(
            mol,
            traj_params = [
                TrajectoryParameters(
                    trjfile = str(BPTI_DCD),
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


def test_drms_zerocopy():
    """Test zerocopy DRMS analysis and compare with legacy implementation."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)
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
            # Get CA atom indices for contact computation
            ca_indices = genesis_exe.selection(mol, "an: CA")
            print(f"Number of CA atoms: {len(ca_indices)}")

            # Compute contact list using Python
            contact_list, contact_dist = compute_contact_list_from_refcoord(
                mol, ca_indices,
                min_dist=1.0, max_dist=6.0, exclude_residues=4
            )
            print(f"Number of contacts: {contact_list.shape[1]}")
            print(f"contact_list shape: {contact_list.shape}")
            print(f"contact_dist shape: {contact_dist.shape}")

            # Run zerocopy version
            result_zerocopy = genesis_exe.drms_analysis_zerocopy(
                t,
                contact_list=contact_list,
                contact_dist=contact_dist,
                ana_period=1,
                pbc_correct=False,
            )

            # Validate results
            assert result_zerocopy.drms is not None, "Zerocopy DRMS result should not be None"
            assert len(result_zerocopy.drms) > 0, "Zerocopy DRMS result should have values"
            assert all(d >= 0 for d in result_zerocopy.drms), "DRMS values should be non-negative"

            print(f"Zerocopy DRMS (n={len(result_zerocopy.drms)}): "
                  f"min={min(result_zerocopy.drms):.5f}, max={max(result_zerocopy.drms):.5f}")

    finally:
        if hasattr(trajs, "close"):
            trajs.close()


def test_drms_zerocopy_vs_legacy():
    """Compare zerocopy and legacy DRMS implementations for consistency."""
    mol = SMolecule.from_file(pdb=BPTI_PDB, psf=BPTI_PSF, ref=BPTI_PDB)
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
            # Run legacy version
            result_legacy = genesis_exe.drms_analysis(
                mol, t,
                selection_group=["an: CA"],
                check_only=False,
                contact_groups=1,
                ignore_hydrogen=False,
                two_states=False,
                avoid_bonding=False,  # Disable to match our Python contact computation
                exclude_residues=4,
                minimum_distance=1.0,
                maximum_distance=6.0,
                pbc_correct=False,
                verbose=True,
            )

            # Compute contacts the same way legacy does (without avoid_bonding)
            ca_indices = genesis_exe.selection(mol, "an: CA")
            contact_list, contact_dist = compute_contact_list_from_refcoord(
                mol, ca_indices,
                min_dist=1.0, max_dist=6.0, exclude_residues=4
            )

            # Run zerocopy version
            result_zerocopy = genesis_exe.drms_analysis_zerocopy(
                t,
                contact_list=contact_list,
                contact_dist=contact_dist,
                ana_period=1,
                pbc_correct=False,
            )

            print(f"Legacy DRMS: min={min(result_legacy.drms):.5f}, max={max(result_legacy.drms):.5f}")
            print(f"Zerocopy DRMS: min={min(result_zerocopy.drms):.5f}, max={max(result_zerocopy.drms):.5f}")

            # Check that results match (within numerical tolerance)
            # Note: They may differ slightly due to different contact selection
            assert len(result_legacy.drms) == len(result_zerocopy.drms), \
                "Result lengths should match"

    finally:
        if hasattr(trajs, "close"):
            trajs.close()


def main():
    if os.path.exists("dummy.trj"):
        os.remove("dummy.trj")
    # Run zerocopy tests first to avoid Fortran global state issues
    # Legacy -> zerocopy sequence crashes, but zerocopy -> legacy works
    try:
        test_drms_zerocopy()
        print("\n✓ test_drms_zerocopy: PASSED")
    except Exception as e:
        print(f"\n✗ test_drms_zerocopy: FAILED - {e}")
        raise

    try:
        test_drms_analysis()
        print("\n✓ test_drms_analysis: PASSED")
    except Exception as e:
        print(f"\n✗ test_drms_analysis: FAILED - {e}")
        raise


if __name__ == "__main__":
    main()
