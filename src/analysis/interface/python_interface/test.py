import os
import ctypes
import pathlib
from libgenesis import LibGenesis
from s_molecule import SMolecule, py2c_s_molecule
import genesis_exe
import msd_reader


def test():
    # 関数を呼び出す
    pdb_filename = pathlib.Path("molecule.pdb")
    with SMolecule.from_pdb_file(pdb_filename) as mol:
        # 結果を処理する
        print("num_atoms = ", mol.num_atoms)
        for i in range(max(0, mol.num_atoms - 5), mol.num_atoms):
            print(mol.atom_coord[i])
            print(mol.atom_no[i], mol.segment_name[i], mol.atom_name[i])

        print("num_atoms = ", mol.num_atoms)
        mol_c = py2c_s_molecule(mol)
        LibGenesis().lib.test_conv_c2f(ctypes.byref(mol_c))


def test_crd():
    # 関数を呼び出す
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")
    ctrl_path = pathlib.Path("test_crd_inp")
    with SMolecule.from_pdb_psf_file(pdb_path, psf_path) as mol:
        with genesis_exe.crd_convert(mol, ctrl_path) as trajs:
            for t in trajs:
                pass


def test_trj_analysis():
    # 関数を呼び出す
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")
    crd_ctrl_path = pathlib.Path("test_no_crd_inp")
    trj_analysis_ctrl_path = pathlib.Path("test_trj_analysis_inp")

    with SMolecule.from_pdb_psf_file(pdb_path, psf_path) as mol:
        with genesis_exe.crd_convert(mol, crd_ctrl_path) as trajs:
            for t in trajs:
                d = genesis_exe.trj_analysis(
                        mol, t, 1, trj_analysis_ctrl_path)
                print(d)


def test_rg_analysis():
    # 関数を呼び出す
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")
    crd_ctrl_path = pathlib.Path("test_no_crd_inp")
    rg_analysis_ctrl_path = pathlib.Path("test_rg_analysis_inp")

    with SMolecule.from_pdb_psf_file(pdb_path, psf_path) as mol:
        with genesis_exe.crd_convert(mol, crd_ctrl_path) as trajs:
            for t in trajs:
                rg = genesis_exe.rg_analysis(mol, t, 1, rg_analysis_ctrl_path)
                print(rg)


def test_rmsd_analysis():
    # 関数を呼び出す
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")
    crd_ctrl_path = pathlib.Path("test_no_crd_inp")
    rmsd_analysis_ctrl_path = pathlib.Path("test_rmsd_analysis_inp")
    with SMolecule.from_pdb_psf_file(pdb_path, psf_path) as mol:
        with genesis_exe.crd_convert(mol, crd_ctrl_path) as trajs:
            for t in trajs:
                rmsd = genesis_exe.rmsd_analysis(mol, t, 1, rmsd_analysis_ctrl_path)
                print(rmsd)


def test_drms_analysis():
    # 関数を呼び出す
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")
    ref_path = pathlib.Path("BPTI_ionize.pdb")
    crd_ctrl_path = pathlib.Path("test_no_crd_inp")
    drms_analysis_ctrl_path = pathlib.Path("test_drms_analysis_inp")

    with SMolecule.from_pdb_psf_ref_file(pdb_path, psf_path, ref_path) as mol:
        with genesis_exe.crd_convert(mol, crd_ctrl_path) as trajs:
            for t in trajs:
                drms = genesis_exe.drms_analysis(mol, t, 1, drms_analysis_ctrl_path)
                print(drms)


def test_msd_analysis():
    # 関数を呼び出す
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")
    crd_ctrl_path = pathlib.Path("test_no_crd_inp")
    msd_analysis_ctrl_path = pathlib.Path("test_msd_analysis_inp")

    with SMolecule.from_pdb_psf_file(pdb_path, psf_path) as mol:
        with genesis_exe.crd_convert(mol, crd_ctrl_path) as trajs:
            for t in trajs:
                msd = genesis_exe.msd_analysis(mol, t, 1, msd_analysis_ctrl_path)
                print(msd)


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


def test_diffusion_analysis():
    ctrl_path = pathlib.Path("test_da_analysis_inp")
    msd = msd_reader.read_msd_from_file("msd.data")
    ret = genesis_exe.diffusion_analysis(msd, ctrl_path)
    print(ret)


def test_avecrd_analysis():
    # 関数を呼び出す
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")
    crd_ctrl_path = pathlib.Path("test_no_crd_inp")
    aa_analysis_ctrl_path = pathlib.Path("test_avecrd_analysis_inp")

    with SMolecule.from_pdb_psf_file(pdb_path, psf_path) as mol:
        with genesis_exe.crd_convert(mol, crd_ctrl_path) as trajs:
            for t in trajs:
                genesis_exe.avecrd_analysis(
                        mol, t, 1, aa_analysis_ctrl_path)


def main():
    if os.path.exists("out"):
      os.remove("out")
    # test()
    # test_crd()
    if os.path.exists("out"):
        os.remove("out")
    # test_trj_analysis()
    if os.path.exists("out"):
        os.remove("out")
    # test_rg_analysis()
    if os.path.exists("out"):
        os.remove("out")
    # test_rmsd_analysis()
    if os.path.exists("out"):
        os.remove("out")
    # test_drms_analysis()
    if os.path.exists("out"):
        os.remove("out")
    # test_msd_analysis()
    if os.path.exists("out"):
        os.remove("out")
    # test_hb_analysis_Count_atom()
    if os.path.exists("out"):
        os.remove("out")
    # test_hb_analysis_Count_snap()
    if os.path.exists("out"):
        os.remove("out")
    # test_diffusion_analysis()
    if os.path.exists("out"):
        os.remove("out")
    test_avecrd_analysis()


if __name__ == "__main__":
    main()
