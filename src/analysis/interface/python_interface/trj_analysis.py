import unittest
import numpy as np
from . import genesis_exe
from .custom_test_case import CustomTestCase


class TestTrjAnalysis(CustomTestCase):
    def test_trj_analysis(self):
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


if __name__ == "__main__":
    unittest.main()
