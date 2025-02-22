import os
import pathlib
from s_molecule import SMolecule
import genesis_exe


def test_trj_analysis_distance():
    # 関数を呼び出す
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")
    crd_ctrl_path = pathlib.Path("test_no_crd_inp")
    trj_analysis_distance_ctrl_path = pathlib.Path("test_trj_analysis_Distance_inp")

    with SMolecule.from_pdb_psf_file(pdb_path, psf_path) as mol:
        with genesis_exe.crd_convert(mol, crd_ctrl_path) as trajs:
            for t in trajs:
                d = genesis_exe.trj_analysis(
                        mol, t, 1, trj_analysis_distance_ctrl_path)
                print(d, flush=True)


def main():
    if os.path.exists("out"):
        os.remove("out")
    test_trj_analysis_distance()


if __name__ == "__main__":
    main()
