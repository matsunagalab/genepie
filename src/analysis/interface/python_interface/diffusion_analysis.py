import os
import pathlib
import genesis_exe
import msd_reader


def test_diffusion_analysis():
    msd = msd_reader.read_msd_from_file("msd.data")
    ret = genesis_exe.diffusion_analysis(
            msd,
            time_step = 2,
            start = "20 %"
            )
    print(ret)


def main():
    if os.path.exists("out"):
        os.remove("out")
    test_diffusion_analysis()


if __name__ == "__main__":
    main()
