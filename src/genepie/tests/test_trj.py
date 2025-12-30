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

    def test_a_trj_analysis_zerocopy_full(self):
        """Test trj_analysis_zerocopy_full (pre-allocated result arrays)."""
        trajs, mol = self.create_traj_by_genesis(
                self.TRJ_PATH, pdb=self.PDB_PATH, psf=self.PSF_PATH)
        for t in trajs:
            # Test with same atom pairs as zerocopy
            ca1_idx = genesis_exe.selection(mol, "rno:1 and an:CA")
            ca2_idx = genesis_exe.selection(mol, "rno:2 and an:CA")
            ca3_idx = genesis_exe.selection(mol, "rno:3 and an:CA")
            ca4_idx = genesis_exe.selection(mol, "rno:4 and an:CA")

            dist_pairs = np.array([
                [ca1_idx[0], ca2_idx[0]],
                [ca2_idx[0], ca3_idx[0]],
            ], dtype=np.int32)
            angle_triplets = np.array([
                [ca1_idx[0], ca2_idx[0], ca3_idx[0]],
            ], dtype=np.int32)
            torsion_quads = np.array([
                [ca1_idx[0], ca2_idx[0], ca3_idx[0], ca4_idx[0]],
            ], dtype=np.int32)

            # Run zerocopy version (for comparison)
            result_zerocopy = genesis_exe.trj_analysis_zerocopy(
                t,
                distance_pairs=dist_pairs,
                angle_triplets=angle_triplets,
                torsion_quadruplets=torsion_quads,
            )

            # Run zerocopy_full version
            result_full = genesis_exe.trj_analysis_zerocopy_full(
                t,
                distance_pairs=dist_pairs,
                angle_triplets=angle_triplets,
                torsion_quadruplets=torsion_quads,
            )

            # Verify shapes match
            self.assertEqual(result_full.distance.shape, result_zerocopy.distance.shape)
            self.assertEqual(result_full.angle.shape, result_zerocopy.angle.shape)
            self.assertEqual(result_full.torsion.shape, result_zerocopy.torsion.shape)

            # Verify values match exactly
            np.testing.assert_allclose(result_full.distance, result_zerocopy.distance,
                                       rtol=1e-10, atol=1e-10)
            np.testing.assert_allclose(result_full.angle, result_zerocopy.angle,
                                       rtol=1e-10, atol=1e-10)
            np.testing.assert_allclose(result_full.torsion, result_zerocopy.torsion,
                                       rtol=1e-10, atol=1e-10)

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

    def test_a_trj_analysis_zerocopy_full_com(self):
        """Test trj_analysis_zerocopy_full_com (pre-allocated result arrays with COM)."""
        trajs, mol = self.create_traj_by_genesis(
                self.TRJ_PATH, pdb=self.PDB_PATH, psf=self.PSF_PATH)
        for t in trajs:
            # Get atom indices for residues
            res1_idx = genesis_exe.selection(mol, "rno:1")
            res2_idx = genesis_exe.selection(mol, "rno:2")
            res3_idx = genesis_exe.selection(mol, "rno:3")

            # COM distance between residue 1 and residue 2
            cdis_groups = [
                (list(res1_idx), list(res2_idx)),
            ]

            # COM angle between residue 1, 2, and 3
            cang_groups = [
                (list(res1_idx), list(res2_idx), list(res3_idx)),
            ]

            # Run zerocopy version for comparison
            result_zerocopy = genesis_exe.trj_analysis_zerocopy_com(
                mol, t,
                cdis_groups=cdis_groups,
                cang_groups=cang_groups,
            )

            # Run zerocopy_full version
            result_full = genesis_exe.trj_analysis_zerocopy_full_com(
                mol, t,
                cdis_groups=cdis_groups,
                cang_groups=cang_groups,
            )

            # Verify shapes match
            self.assertEqual(result_full.cdis.shape, result_zerocopy.cdis.shape)
            self.assertEqual(result_full.cang.shape, result_zerocopy.cang.shape)

            # Verify values match exactly
            np.testing.assert_allclose(result_full.cdis, result_zerocopy.cdis,
                                       rtol=1e-10, atol=1e-10)
            np.testing.assert_allclose(result_full.cang, result_zerocopy.cang,
                                       rtol=1e-10, atol=1e-10)


if __name__ == "__main__":
    unittest.main()
