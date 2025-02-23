import os
import pathlib
import genesis_exe


def test_mbar_analysis():
    ctrl_path = pathlib.Path("test_mbar_analysis_umbrella_1d_inp")
    fene = genesis_exe.mbar_analysis(61, 1, ctrl_path)
    print(fene)


def main():
    if os.path.exists("out"):
        os.remove("out")
    if os.path.exists("fene.dat"):
        os.remove("fene.dat")
    test_mbar_analysis()


if __name__ == "__main__":
    main()
