import pathlib
import unittest
from s_molecule import SMolecule


class TestMDTraj(unittest.TestCase):

    def test_to_mdtraj_topology(self):
        pdb_path = pathlib.Path("BPTI_ionize.pdb")
        with SMolecule.from_file(pdb=pdb_path) as mol:
            topo = mol.to_mdtraj_topology()
            SMolecule.from_mdtraj_topology(topo)


if __name__ == "__main__":
    unittest.main()
