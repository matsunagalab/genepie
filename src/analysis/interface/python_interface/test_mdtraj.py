import unittest
import mdtraj
from .s_trajectories import STrajectories
from .custom_test_case import CustomTestCase


class TestMDTraj(CustomTestCase):

    def test_from_mdtraj_trajectory(self):
        mdt = mdtraj.load(self.TRJ_PATH, top=self.PDB_PATH)

        trj, mol = STrajectories.from_mdtraj_trajectory(mdt)
        gtrajs, gmol = self.create_traj_by_genesis(
                self.TRJ_PATH, pdb=self.PDB_PATH)
        with trj, gtrajs:
            self.assertAlmostEqual(gtrajs[0], trj)
            # self.assertAlmostEqual(gmol, mol)

    def test_to_mdtraj_trajectory(self):
        strajs, smol = self.create_traj_by_genesis(
                self.TRJ_PATH, pdb=self.PDB_PATH, psf=self.PSF_PATH)
        with strajs:
            for t in strajs:
                mdt = t.to_mdtraj_trajectory(smol)
                gt, gm = STrajectories.from_mdtraj_trajectory(mdt)
                with gt:
                    self.assertAlmostEqual(t, gt)
                    # self.assertAlmostEqual(smol, gm)


if __name__ == "__main__":
    unittest.main()
