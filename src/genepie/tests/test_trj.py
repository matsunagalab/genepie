# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------
import unittest
import numpy as np
from .. import genesis_exe
from ..custom_test_case import CustomTestCase


class TestTrjAnalysis(CustomTestCase):
    # Note: Test methods are ordered alphabetically by unittest.
    # Zerocopy tests must run BEFORE legacy tests to avoid Fortran state issues.
    # So we use 'test_a_*' prefix for zerocopy and 'test_z_*' for legacy.

    def test_z_trj_analysis(self):
        trajs, mol = self.create_traj_by_genesis(
                self.TRJ_PATH, pdb=self.PDB_PATH, psf=self.PSF_PATH)
        for t in trajs:
            r = genesis_exe.trj_analysis(
                    mol, t,
                    distance = ["BPTI:1:ARG:CA  BPTI:2:PRO:CA",
                                "BPTI:2:PRO:CA  BPTI:3:ASP:CA"],
                    angle = ["BPTI:1:ARG:CA  BPTI:2:PRO:CA  BPTI:3:ASP:CA",],
                    torsion = ["BPTI:1:ARG:CA  BPTI:2:PRO:CA  "
                               + "BPTI:3:ASP:CA   BPTI:4:PHE:CA", ])
            ta_test_root = self.TEST_ROOT / "test_analysis/test_trj_analysis"
            ref_dist = np.loadtxt(ta_test_root / "Distance/ref")
            self.assertAlmostEqual(ref_dist[:,1:], r.distance, places=3)
            ref_ang = np.loadtxt(ta_test_root / "Angle/ref")
            self.assertAlmostEqual(ref_ang[:,1:], r.angle, places=3)
            ref_tor = np.loadtxt(ta_test_root / "Dihedral/ref")
            self.assertAlmostEqual(ref_tor[:,1:], r.torsion, places=3)

    def test_a_trj_analysis_zerocopy_com(self):
        """Test trj_analysis_zerocopy_com with COM calculations."""
        trajs, mol = self.create_traj_by_genesis(
                self.TRJ_PATH, pdb=self.PDB_PATH, psf=self.PSF_PATH)
        for t in trajs:
            # Test COM distance using selection to get atom indices
            # First, get indices for residue 1 and residue 2
            res1_idx = genesis_exe.selection(mol, "rno:1")  # 1-indexed
            res2_idx = genesis_exe.selection(mol, "rno:2")
            res3_idx = genesis_exe.selection(mol, "rno:3")

            # Test COM distance between residue 1 and residue 2
            cdis_groups = [
                (list(res1_idx), list(res2_idx)),
            ]

            # Test COM angle between residue 1, 2, and 3
            cang_groups = [
                (list(res1_idx), list(res2_idx), list(res3_idx)),
            ]

            result = genesis_exe.trj_analysis_zerocopy_com(
                mol, t,
                cdis_groups=cdis_groups,
                cang_groups=cang_groups,
            )

            # Verify we got results
            self.assertIsNotNone(result.cdis)
            self.assertIsNotNone(result.cang)

            # Check shapes
            n_frames = t.nframe
            self.assertEqual(result.cdis.shape, (n_frames, 1))
            self.assertEqual(result.cang.shape, (n_frames, 1))

            # Verify distance is positive
            self.assertTrue(np.all(result.cdis > 0))

            # Verify angle is in valid range (0-180 degrees)
            self.assertTrue(np.all(result.cang >= 0))
            self.assertTrue(np.all(result.cang <= 180))

            # Verify values are reasonable (not NaN or Inf)
            self.assertTrue(np.all(np.isfinite(result.cdis)))
            self.assertTrue(np.all(np.isfinite(result.cang)))


if __name__ == "__main__":
    unittest.main()
