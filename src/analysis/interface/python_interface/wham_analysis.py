import os
import pathlib
import genesis_exe


def test_wham_analysis():
    ctrl_path = pathlib.Path("test_wham_analysis_inp")
    pmf = genesis_exe.wham_analysis(300, 2, ctrl_path)
    print(pmf)


def main():
    if os.path.exists("out"):
        os.remove("out")
    if os.path.exists("pmf"):
        os.remove("pmf")
    test_wham_analysis()


if __name__ == "__main__":
    main()
