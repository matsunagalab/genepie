import os
import pathlib
from s_molecule import SMolecule
import genesis_exe


def test_hb_analysis_Count_atom():
    # 関数を呼び出す
    pdb_path = pathlib.Path("RALP_DPPC_run.pdb")
    psf_path = pathlib.Path("RALP_DPPC.psf")
    crd_ctrl_path = pathlib.Path("test_no_crd_hb_analysis_inp")

    with SMolecule.from_file(pdb=pdb_path, psf=psf_path) as mol:
        with genesis_exe.crd_convert(mol, crd_ctrl_path) as trajs:
            for t in trajs:
                d = genesis_exe.hb_analysis(
                        mol, t,
                        selection_group = ["sid:PROA",
                                           "resname:DPPC & (an:O11 | an:O12 | an:O13 | an:O14)", ],
                        check_only = False,
                        output_type = "Count_atom",
                        solvent_list  = "DPPC",
                        analysis_atom = 1,
                        target_atom   = 2,
                        boundary_type = "PBC",
                        hb_distance   = 3.4,
                        dha_angle     = 120.0,
                        hda_angle     = 30.0,
                        )
                print(d, flush=True)


def main():
    if os.path.exists("out"):
        os.remove("out")
    test_hb_analysis_Count_atom()


if __name__ == "__main__":
    main()
