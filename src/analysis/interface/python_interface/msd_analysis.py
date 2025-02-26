import os
import pathlib
from s_molecule import SMolecule
import genesis_exe


def test_msd_analysis():
    # 関数を呼び出す
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")
    crd_ctrl_path = pathlib.Path("test_no_crd_msd_analysis_inp")
    msd_analysis_ctrl_path = pathlib.Path("test_msd_analysis_inp")

    with SMolecule.from_pdb_psf_file(pdb_path, psf_path) as mol:
        with genesis_exe.crd_convert(mol, crd_ctrl_path) as trajs:
            for t in trajs:
                msd = genesis_exe.msd_analysis(mol, t, 1, msd_analysis_ctrl_path)
                print(msd)


def main():
    if os.path.exists("out"):
        os.remove("out")
    if os.path.exists("out1"):
        os.remove("out1")
    test_msd_analysis()


if __name__ == "__main__":
    main()
