# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------
import unittest
import numpy as np
from .. import genesis_exe
from ..custom_test_case import CustomTestCase
from ..exceptions import (
    GenesisFortranNotSupportedError,
    GenesisValidationError,
)


class TestWhamAnalysis(CustomTestCase):
    """Regression tests for wham_analysis against the CLI reference output."""

    # Same tolerance the CLI regression harness uses for this test case, see
    # tests/regression_test/test_analysis/test_wham_analysis/triala_cv/config.ini
    TOLERANCE = 1.0e-8

    # 14 umbrella windows along an end-to-end distance of trialanine.
    CONSTANT = (1.2,) * 14
    REFERENCE = (1.80, 2.72, 3.64, 4.56, 5.48, 6.40, 7.32,
                 8.24, 9.16, 10.08, 11.00, 11.92, 12.84, 13.76)

    def _cvfile(self):
        return str(self.TEST_ROOT / "test_analysis" / "trajectories"
                   / "triala_cv" / "{}.dis")

    def _run(self, **overrides):
        kwargs = dict(
            cvfile=self._cvfile(),
            dimension=1,
            nblocks=1,
            temperature=300.0,
            tolerance=10E-08,
            rest_function=(1,),
            grids=((0.0, 15.0, 301),),
            constant=(self.CONSTANT,),
            reference=(self.REFERENCE,),
            is_periodic=(False,),
        )
        kwargs.update(overrides)
        return genesis_exe.wham_analysis(**kwargs)

    def test_wham_analysis_matches_cli_reference(self):
        """PMF must reproduce the value the CLI wham_analysis writes."""
        pmf = self._run()

        ref = np.loadtxt(
            self.TEST_ROOT / "test_analysis" / "test_wham_analysis"
            / "triala_cv" / "ref")

        self.assertEqual(ref.shape, pmf.shape)
        self.assertTrue(np.all(np.isfinite(pmf)))
        # Column 0 holds the bin centers, column 1 the free energy.
        self.assertAlmostEqual(ref, pmf, delta=self.TOLERANCE)

    def test_wham_analysis_rejects_dcd_input(self):
        """DCD input is unsupported because molecules are never defined."""
        with self.assertRaises(GenesisFortranNotSupportedError):
            self._run(dcdfile="whatever.dcd")

    def test_wham_analysis_requires_cvfile(self):
        with self.assertRaises(GenesisValidationError):
            self._run(cvfile=None)


if __name__ == "__main__":
    unittest.main()
