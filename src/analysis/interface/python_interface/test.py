import ctypes
import pathlib
from .libgenesis import LibGenesis
from .s_molecule import SMolecule


def test():
    # 関数を呼び出す
    pdb_filename = pathlib.Path("molecule.pdb")
    mol = SMolecule.from_file(pdb=pdb_filename)
    # 結果を処理する
    print("num_atoms = ", mol.num_atoms)
    for i in range(max(0, mol.num_atoms - 5), mol.num_atoms):
        print(mol.atom_coord[i])
        print(mol.atom_no[i], mol.segment_name[i], mol.atom_name[i])

    print("num_atoms = ", mol.num_atoms)
    mol_c = mol.to_SMoleculeC()
    LibGenesis().lib.test_conv_c2f(ctypes.byref(mol_c))


def main():
    test()


if __name__ == "__main__":
    main()
