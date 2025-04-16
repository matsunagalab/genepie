import pathlib
import unittest
from ctrl_files import TrajectoryParameters
import genesis_exe
from s_molecule import SMolecule
from s_trajectories import STrajectories
import mdtraj

class TestMDTraj(unittest.TestCase):

    def test_from_mdtraj_trajectory(self):
        pdb_path = pathlib.Path("BPTI_ionize.pdb")
        trj_path = pathlib.Path("BPTI_run.dcd")
        mdt = mdtraj.load(trj_path, top=pdb_path)

        trj, mol = STrajectories.from_mdtraj_trajectory(mdt)
        d = genesis_exe.trj_analysis(
                mol, trj,
                distance = ["BPTI:1:ARG:CA  BPTI:2:PRO:CA",
                            "BPTI:2:PRO:CA  BPTI:3:ASP:CA"],
                angle =    ["BPTI:1:ARG:CA  BPTI:2:PRO:CA  BPTI:3:ASP:CA",],
                torsion =  ["BPTI:1:ARG:CA  BPTI:2:PRO:CA  BPTI:3:ASP:CA   BPTI:4:PHE:CA", ])
        print(d.distance)
        print(d.angle)
        print(d.torsion)

    def test_to_mdtraj_trajectory(self):
        pdb_path = pathlib.Path("BPTI_ionize.pdb")
        with SMolecule.from_file(pdb=pdb_path) as mol:
            with genesis_exe.crd_convert(
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
                    ) as trajs:
                for t in trajs:
                    _ = t.to_mdtraj_trajectory(mol)


if __name__ == "__main__":
    unittest.main()
