# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent))
    __package__ = "python_interface"
# --------------------------------------------
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
