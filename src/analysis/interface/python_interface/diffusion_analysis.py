from . import genesis_exe
import numpy as np


def test_diffusion_analysis():
    msd = np.loadtxt("msd.data", dtype=np.float64)
    print(msd)
    ret = genesis_exe.diffusion_analysis(
            msd,
            time_step = 2,
            start = "20 %"
            )
    print(ret)


def main():
    test_diffusion_analysis()


if __name__ == "__main__":
    main()
