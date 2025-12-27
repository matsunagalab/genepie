# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------
from .. import genesis_exe
import numpy as np


def test_diffusion_analysis():
    msd = np.loadtxt("msd.data", dtype=np.float64)
    assert msd.ndim == 2, "MSD data should be 2D array"
    print(f"MSD data shape: {msd.shape}")

    ret = genesis_exe.diffusion_analysis(
            msd,
            time_step = 2,
            start = "20 %"
            )

    # Validate diffusion coefficient results
    assert ret is not None, "Diffusion result should not be None"
    assert len(ret) > 0, "Diffusion result should have at least one value"
    # Diffusion coefficients should be non-negative
    assert all(d >= 0 for d in ret), "Diffusion coefficients should be non-negative"
    print(f"Diffusion coefficients (n={len(ret)}): {ret}")


def main():
    test_diffusion_analysis()


if __name__ == "__main__":
    main()
