import os
import pathlib
from s_molecule import SMolecule
import genesis_exe


def test_drms_analysis():
    # 関数を呼び出す
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")
    ref_path = pathlib.Path("BPTI_ionize.pdb")
    crd_ctrl_path = pathlib.Path("test_no_crd_inp")

    with SMolecule.from_pdb_psf_ref_file(pdb_path, psf_path, ref_path) as mol:
        with genesis_exe.crd_convert(mol, crd_ctrl_path) as trajs:
            for t in trajs:
                drms, = genesis_exe.drms_analysis(
                        mol, t,
                        selection_group = ["an: CA", ],
                        check_only = False,
                        contact_groups = 1,
                        ignore_hydrogen  = False,
                        two_states       = False,
                        avoid_bonding    = True,
                        exclude_residues = 4,
                        minimum_distance = 1.0,
                        maximum_distance = 6.0,
                        pbc_correct      = False,
                        verbose          = True,
                        )
                print(drms)


def main():
    if os.path.exists("out"):
        os.remove("out")
    test_drms_analysis()


if __name__ == "__main__":
    main()
