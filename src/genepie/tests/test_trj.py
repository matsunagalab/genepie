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
from ..s_molecule import SMolecule
from ..custom_test_case import CustomTestCase


class TestTrjAnalysis(CustomTestCase):

    def test_trj_analysis_atoms(self):
        """Test trj_analysis with atom-based measurements."""
        trajs, mol = self.create_traj_by_genesis(
                self.TRJ_PATH, pdb=self.PDB_PATH, psf=self.PSF_PATH)
        for t in trajs:
            # Get CA atom indices for first 4 residues
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

            result = genesis_exe.trj_analysis(
                t,
                distance_pairs=dist_pairs,
                angle_triplets=angle_triplets,
                torsion_quadruplets=torsion_quads,
            )

            # Validate results
            self.assertIsNotNone(result.distance)
            self.assertIsNotNone(result.angle)
            self.assertIsNotNone(result.torsion)

            # Check shapes
            n_frames = t.nframe
            self.assertEqual(result.distance.shape, (n_frames, 2))
            self.assertEqual(result.angle.shape, (n_frames, 1))
            self.assertEqual(result.torsion.shape, (n_frames, 1))

            # Check values are reasonable
            self.assertTrue(np.all(result.distance > 0))  # Distances positive
            self.assertTrue(np.all(result.angle >= 0))    # Angles 0-180
            self.assertTrue(np.all(result.angle <= 180))

            # Compare with reference data
            ta_test_root = self.TEST_ROOT / "test_analysis/test_trj_analysis"
            ref_dist = np.loadtxt(ta_test_root / "Distance/ref")
            self.assertAlmostEqual(ref_dist[:, 1:], result.distance, places=3)
            ref_ang = np.loadtxt(ta_test_root / "Angle/ref")
            self.assertAlmostEqual(ref_ang[:, 1:], result.angle, places=3)
            ref_tor = np.loadtxt(ta_test_root / "Dihedral/ref")
            self.assertAlmostEqual(ref_tor[:, 1:], result.torsion, places=3)

    def test_trj_lazy_vs_memory(self):
        """Lazy trj_analysis must match the memory-based result."""
        mol = SMolecule.from_file(
            pdb=self.PDB_PATH, psf=self.PSF_PATH, ref=self.PDB_PATH)

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

        common = dict(
            trj_files=[str(self.TRJ_PATH)],
            trj_format="DCD",
            trj_type="COOR+BOX",
            selection="all",
        )
        for ana_period in (1, 2):
            mem_trajs, _ = genesis_exe.crd_convert(mol, lazy=False, **common)
            lazy_trajs, _ = genesis_exe.crd_convert(mol, lazy=True, **common)
            self.assertTrue(lazy_trajs[0].is_lazy)

            mem = genesis_exe.trj_analysis(
                mem_trajs[0],
                distance_pairs=dist_pairs,
                angle_triplets=angle_triplets,
                torsion_quadruplets=torsion_quads,
                ana_period=ana_period,
            )
            lazy = genesis_exe.trj_analysis(
                lazy_trajs[0],
                distance_pairs=dist_pairs,
                angle_triplets=angle_triplets,
                torsion_quadruplets=torsion_quads,
                ana_period=ana_period,
            )

            self.assertEqual(mem.distance.shape, lazy.distance.shape)
            np.testing.assert_allclose(
                lazy.distance, mem.distance, rtol=1e-4, atol=1e-5)
            np.testing.assert_allclose(
                lazy.angle, mem.angle, rtol=1e-4, atol=1e-5)
            np.testing.assert_allclose(
                lazy.torsion, mem.torsion, rtol=1e-4, atol=1e-5)

    def test_trj_com_lazy_vs_memory(self):
        """Lazy COM trj_analysis must match the memory-based result."""
        mol = SMolecule.from_file(
            pdb=self.PDB_PATH, psf=self.PSF_PATH, ref=self.PDB_PATH)

        res1_idx = genesis_exe.selection(mol, "rno:1")
        res2_idx = genesis_exe.selection(mol, "rno:2")
        res3_idx = genesis_exe.selection(mol, "rno:3")
        res4_idx = genesis_exe.selection(mol, "rno:4")

        ca1_idx = genesis_exe.selection(mol, "rno:1 and an:CA")
        ca2_idx = genesis_exe.selection(mol, "rno:2 and an:CA")

        # Mix atom-based and COM-based measurements to exercise every stream.
        dist_pairs = np.array([[ca1_idx[0], ca2_idx[0]]], dtype=np.int32)
        cdis_groups = [(list(res1_idx), list(res2_idx))]
        cang_groups = [(list(res1_idx), list(res2_idx), list(res3_idx))]
        ctor_groups = [
            (list(res1_idx), list(res2_idx), list(res3_idx), list(res4_idx)),
        ]

        common = dict(
            trj_files=[str(self.TRJ_PATH)],
            trj_format="DCD",
            trj_type="COOR+BOX",
            selection="all",
        )
        for ana_period in (1, 2):
            mem_trajs, _ = genesis_exe.crd_convert(mol, lazy=False, **common)
            lazy_trajs, _ = genesis_exe.crd_convert(mol, lazy=True, **common)
            self.assertTrue(lazy_trajs[0].is_lazy)

            kwargs = dict(
                distance_pairs=dist_pairs,
                cdis_groups=cdis_groups,
                cang_groups=cang_groups,
                ctor_groups=ctor_groups,
                molecule=mol,
                ana_period=ana_period,
            )
            mem = genesis_exe.trj_analysis(mem_trajs[0], **kwargs)
            lazy = genesis_exe.trj_analysis(lazy_trajs[0], **kwargs)

            self.assertEqual(mem.distance.shape, lazy.distance.shape)
            np.testing.assert_allclose(
                lazy.distance, mem.distance, rtol=1e-4, atol=1e-5)
            np.testing.assert_allclose(
                lazy.com_distance, mem.com_distance, rtol=1e-4, atol=1e-5)
            np.testing.assert_allclose(
                lazy.com_angle, mem.com_angle, rtol=1e-4, atol=1e-5)
            np.testing.assert_allclose(
                lazy.com_torsion, mem.com_torsion, rtol=1e-4, atol=1e-5)

    def test_trj_analysis_com(self):
        """Test trj_analysis with COM-based measurements."""
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

            result = genesis_exe.trj_analysis(
                t,
                cdis_groups=cdis_groups,
                cang_groups=cang_groups,
                molecule=mol,
            )

            # Verify we got results
            self.assertIsNotNone(result.com_distance)
            self.assertIsNotNone(result.com_angle)

            # Check shapes
            n_frames = t.nframe
            self.assertEqual(result.com_distance.shape, (n_frames, 1))
            self.assertEqual(result.com_angle.shape, (n_frames, 1))

            # Verify distance is positive
            self.assertTrue(np.all(result.com_distance > 0))

            # Verify angle is in valid range (0-180 degrees)
            self.assertTrue(np.all(result.com_angle >= 0))
            self.assertTrue(np.all(result.com_angle <= 180))

            # Verify values are reasonable (not NaN or Inf)
            self.assertTrue(np.all(np.isfinite(result.com_distance)))
            self.assertTrue(np.all(np.isfinite(result.com_angle)))

    def test_trj_analysis_mixed(self):
        """Test trj_analysis with both atom-based and COM-based measurements."""
        trajs, mol = self.create_traj_by_genesis(
                self.TRJ_PATH, pdb=self.PDB_PATH, psf=self.PSF_PATH)
        for t in trajs:
            # Atom-based measurements
            ca1_idx = genesis_exe.selection(mol, "rno:1 and an:CA")
            ca2_idx = genesis_exe.selection(mol, "rno:2 and an:CA")
            dist_pairs = np.array([
                [ca1_idx[0], ca2_idx[0]],
            ], dtype=np.int32)

            # COM-based measurements
            res1_idx = genesis_exe.selection(mol, "rno:1")
            res2_idx = genesis_exe.selection(mol, "rno:2")
            cdis_groups = [
                (list(res1_idx), list(res2_idx)),
            ]

            result = genesis_exe.trj_analysis(
                t,
                distance_pairs=dist_pairs,
                cdis_groups=cdis_groups,
                molecule=mol,
            )

            # Verify we got both atom-based and COM-based results
            self.assertIsNotNone(result.distance)
            self.assertIsNotNone(result.com_distance)

            # Check shapes
            n_frames = t.nframe
            self.assertEqual(result.distance.shape, (n_frames, 1))
            self.assertEqual(result.com_distance.shape, (n_frames, 1))

            # Verify values are reasonable
            self.assertTrue(np.all(result.distance > 0))
            self.assertTrue(np.all(result.com_distance > 0))
            self.assertTrue(np.all(np.isfinite(result.distance)))
            self.assertTrue(np.all(np.isfinite(result.com_distance)))


if __name__ == "__main__":
    unittest.main()
