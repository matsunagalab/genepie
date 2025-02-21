import os
import ctypes
import pathlib
from libgenesis import LibGenesis
from s_molecule import SMolecule, py2c_s_molecule
import genesis_exe


def test_hb_analysis_Count_atom():
    # 関数を呼び出す
    pdb_path = pathlib.Path("RALP_DPPC_run.pdb")
    psf_path = pathlib.Path("RALP_DPPC.psf")
    crd_ctrl_path = pathlib.Path("test_no_crd_hb_analysis_inp")
    hb_analysis_ctrl_path = pathlib.Path("test_hb_analysis_Count_atom_inp")

    with SMolecule.from_pdb_psf_file(pdb_path, psf_path) as mol:
        with genesis_exe.crd_convert(mol, crd_ctrl_path) as trajs:
            for t in trajs:
                genesis_exe.hb_analysis(
                        mol, t, 1, hb_analysis_ctrl_path)


def test_hb_analysis_Count_snap():
    # 関数を呼び出す
    pdb_path = pathlib.Path("RALP_DPPC_run.pdb")
    psf_path = pathlib.Path("RALP_DPPC.psf")
    crd_ctrl_path = pathlib.Path("test_no_crd_hb_analysis_inp")
    hb_analysis_ctrl_path = pathlib.Path("test_hb_analysis_Count_snap_inp")

    with SMolecule.from_pdb_psf_file(pdb_path, psf_path) as mol:
        with genesis_exe.crd_convert(mol, crd_ctrl_path) as trajs:
            for t in trajs:
                genesis_exe.hb_analysis(
                        mol, t, 1, hb_analysis_ctrl_path)


def main():
    if os.path.exists("out"):
        os.remove("out")
    if os.path.exists("out1"):
        os.remove("out1")
    test_hb_analysis_Count_atom()
    if os.path.exists("out"):
        os.remove("out")
    if os.path.exists("out2"):
        os.remove("out2")
    # test_hb_analysis_Count_snap()


if __name__ == "__main__":
    main()
