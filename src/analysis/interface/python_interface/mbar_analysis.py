import os
import ctypes
import pathlib
from libgenesis import LibGenesis
from s_molecule import SMolecule, py2c_s_molecule
import genesis_exe


def test_mbar_analysis_umbrella_1d():
    ctrl_path = pathlib.Path("test_mbar_analysis_umbrella_1d_inp")
    genesis_exe.mbar_analysis(ctrl_path)


def test_mbar_analysis_umbrella_block():
    ctrl_path = pathlib.Path("test_mbar_analysis_umbrella_block_inp")
    genesis_exe.mbar_analysis(ctrl_path)


def main():
    if os.path.exists("out"):
        os.remove("out")
    if os.path.exists("fene.dat"):
        os.remove("fene.dat")
    test_mbar_analysis_umbrella_1d()
    if os.path.exists("out"):
        os.remove("out")
    if os.path.exists("fene.dat"):
        os.remove("fene.dat")
    test_mbar_analysis_umbrella_block()


if __name__ == "__main__":
    main()
