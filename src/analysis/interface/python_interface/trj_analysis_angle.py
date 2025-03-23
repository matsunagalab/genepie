import os
import pathlib
from s_molecule import SMolecule
import genesis_exe


def test_trj_analysis_angle():
    # 関数を呼び出す
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")
    crd_ctrl_path = pathlib.Path("test_no_crd_inp")

    with SMolecule.from_file(pdb=pdb_path, psf=psf_path) as mol:
        with genesis_exe.crd_convert(mol, crd_ctrl_path) as trajs:
            for t in trajs:
                d = genesis_exe.trj_analysis(
                        mol, t,
                        angle = ["BPTI:1:ARG:CA  BPTI:2:PRO:CA  BPTI:3:ASP:CA",])
                print(d.angle, flush=True)


def main():
    if os.path.exists("out"):
        os.remove("out")
    test_trj_analysis_angle()


if __name__ == "__main__":
    main()
