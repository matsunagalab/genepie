import os
import pathlib
from s_molecule import SMolecule
import genesis_exe


def test_kmeans_clustering():
    # 関数を呼び出す
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")
    crd_ctrl_path = pathlib.Path("test_no_crd_inp")

    with SMolecule.from_pdb_psf_file(pdb_path, psf_path) as mol:
        with genesis_exe.crd_convert(mol, crd_ctrl_path) as trajs:
            for t in trajs:
                pdb, cluster_idx = genesis_exe.kmeans_clustering(
                        mol, t,
                        selection_group = ["an:CA", ],
                        fitting_method = "TR+ROT",
                        fitting_atom = 1,
                        check_only = False,
                        allow_backup    = False,
                        analysis_atom   = 1,
                        num_clusters    = 2,
                        max_iteration   = 100,
                        stop_threshold  = 98.0,
                        num_iterations  = 5,
                        trjout_atom     = 1,
                        trjout_format   = "DCD",
                        trjout_type     = "COOR",
                        iseed           = 3141592,
                        )
                print(pdb)
                print(cluster_idx)


def main():
    if os.path.exists("out"):
        os.remove("out")
    if os.path.exists("out1"):
        os.remove("out1")
    if os.path.exists("output_1.pdb"):
        os.remove("output_1.pdb")
    if os.path.exists("output_2.pdb"):
        os.remove("output_2.pdb")
    if os.path.exists("output1.trj"):
        os.remove("output1.trj")
    if os.path.exists("output2.trj"):
        os.remove("output2.trj")
    test_kmeans_clustering()


if __name__ == "__main__":
    main()
