import os
import pathlib
from s_molecule import SMolecule
import genesis_exe


def test_crd():
    # 関数を呼び出す
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")
    ctrl_path = pathlib.Path("test_crd_inp")
    with SMolecule.from_file(pdb=pdb_path, psf=psf_path) as mol:
        with genesis_exe.crd_convert(mol, ctrl_path) as trajs:
            for t in trajs:
                pass


def main():
    if os.path.exists("out"):
        os.remove("out")
    test_crd()


if __name__ == "__main__":
    main()
