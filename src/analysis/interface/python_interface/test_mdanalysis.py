import pathlib
import unittest
from s_trajectories import STrajectories
import MDAnalysis as mda
from ctrl_files import TrajectoryParameters
from s_molecule import SMolecule
import genesis_exe


class TestMdAnalysis(unittest.TestCase):

    def test_from_mdanalysis_universe(self):
        pdb_path = pathlib.Path("BPTI_ionize.pdb")
        trj_file = pathlib.Path("BPTI_run.dcd")
        uni = mda.Universe(pdb_path, trj_file)
        trj, mol = STrajectories.from_mdanalysis_universe(uni)
        with trj:
            d = genesis_exe.trj_analysis(
                    mol, trj,
                    distance = ["BPTI:1:ARG:CA  BPTI:2:PRO:CA",
                                "BPTI:2:PRO:CA  BPTI:3:ASP:CA"],
                    angle =    ["BPTI:1:ARG:CA  BPTI:2:PRO:CA  BPTI:3:ASP:CA",],
                    torsion =  ["BPTI:1:ARG:CA  BPTI:2:PRO:CA  BPTI:3:ASP:CA   BPTI:4:PHE:CA", ])
            print(d.distance)
            print(d.angle)
            print(d.torsion)
        # with SMolecule.from_file(pdb=pdb_path) as gmol:
            # self.assertEqual(mol.num_deg_freedom, gmol.num_deg_freedom)
            # self.assertEqual(mol.num_atoms, gmol.num_atoms)
            # self.assertEqual(mol.num_bonds, gmol.num_bonds)
            # self.assertEqual(mol.num_enm_bonds, gmol.num_enm_bonds)
            # self.assertEqual(mol.num_angles, gmol.num_angles)
            # self.assertEqual(mol.num_dihedrals, gmol.num_dihedrals)
            # self.assertEqual(mol.num_impropers, gmol.num_impropers)

    @unittest.skip
    def test_to_mdanalysis_universe(self):
        pdb_path = pathlib.Path("BPTI_ionize.pdb")
        psf_path = pathlib.Path("BPTI_ionize.psf")
        with SMolecule.from_file(pdb=pdb_path, psf=psf_path) as mol:
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
                    selection_group = ["an:CA", ],
                    fitting_method = "TR+ROT",
                    fitting_atom = 1,
                    check_only = False,
                    centering = True,
                    centering_atom = 1,
                    center_coord = (0.0, 0.0, 0.0),
                    pbc_correct = "molecule",
                    rename_res = ["HSE HIS", "HSE HIS"],
                    ) as trajs:
                for t in trajs:
                    _ = t.to_mdanalysis_universe(mol)
                    # uni.atoms.write("test.pdb")

if __name__ == "__main__":
    unittest.main()
