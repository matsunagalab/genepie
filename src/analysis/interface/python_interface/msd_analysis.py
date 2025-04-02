import os
import pathlib
from ctrl_files import TrajectoryParameters
from s_molecule import SMolecule
import genesis_exe


def test_msd_analysis():
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")

    with SMolecule.from_file(pdb=pdb_path, psf=psf_path) as mol:
        with genesis_exe.crd_convert(
                mol,
                traj_params = [
                    TrajectoryParameters(
                        trjfile = "BPTI_run.dcd",
                        md_step = 10000,
                        mdout_period = 1000,
                        ana_period = 1,
                        ),
                    ],
                trj_format = "DCD",
                trj_type = "COOR+BOX",
                trj_natom = 0,
                selection_group = ["all", ],
                fitting_method = "NO",
                fitting_atom = 1,
                check_only = False,
                pbc_correct = "NO",
                ) as trajs:
            for t in trajs:
                d = genesis_exe.msd_analysis(
                        mol, t,
                        selection_group = ["rnam:TIP3", ],
                        selection = [1, ],
                        mode = ["ALL", ],
                        oversample = True,
                        delta = 9,
                        )
                print(d.msd, flush=True)


def main():
    if os.path.exists("dummy.trj"):
        os.remove("dummy.trj")
    test_msd_analysis()


if __name__ == "__main__":
    main()
