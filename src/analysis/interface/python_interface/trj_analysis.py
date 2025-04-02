import os
import pathlib
from ctrl_files import TrajectoryParameters
from s_molecule import SMolecule
import genesis_exe


def test_trj_analysis_distance():
    pdb_path = pathlib.Path("BPTI_ionize.pdb")
    psf_path = pathlib.Path("BPTI_ionize.psf")

    with SMolecule.from_file(pdb=pdb_path, psf=psf_path) as mol:
        with genesis_exe.crd_convert(
                mol,
                traj_params = [
                    TrajectoryParameters(
                        trjfile = "BPTI_run.dcd",
                        md_step = 10,
                        mdout_period = 1,
                        ana_period = 1,
                        repeat = 1,
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
                d = genesis_exe.trj_analysis(
                        mol, t,
                        distance = ["BPTI:1:ARG:CA  BPTI:2:PRO:CA",
                                    "BPTI:2:PRO:CA  BPTI:3:ASP:CA"],
                        angle =    ["BPTI:1:ARG:CA  BPTI:2:PRO:CA  BPTI:3:ASP:CA",],
                        torsion =  ["BPTI:1:ARG:CA  BPTI:2:PRO:CA  BPTI:3:ASP:CA   BPTI:4:PHE:CA", ])
                print(d.distance, flush=True)
                print(d.angle, flush=True)
                print(d.torsion, flush=True)


def main():
    if os.path.exists("dummy.trj"):
        os.remove("dummy.trj")
    test_trj_analysis_distance()


if __name__ == "__main__":
    main()
