import ctypes
from libgenesis import LibGenesis
from s_molecule import SMoleculeC, c2py_s_molecule, py2c_s_molecule


def test():
    # 関数を呼び出す
    pdb_filename = b"molecule.pdb"
    mol_c = SMoleculeC()

    LibGenesis().lib.define_molecule_from_pdb(
            pdb_filename,
            ctypes.byref(mol_c))
    mol = c2py_s_molecule(mol_c)

    # 結果を処理する
    print("num_atoms = ", mol.num_atoms)
    for i in range(max(0, mol.num_atoms - 5), mol.num_atoms):
        print(mol.atom_coord[i])
        print(mol.atom_no[i], mol.segment_name[i], mol.atom_name[i])

    # メモリを解放する
    LibGenesis().lib.deallocate_s_molecule_c(ctypes.byref(mol_c))

    print("num_atoms = ", mol.num_atoms)
    mol_c = py2c_s_molecule(mol)
    LibGenesis().lib.test_conv_c2f(ctypes.byref(mol_c))

if __name__ == "__main__":
    test()
