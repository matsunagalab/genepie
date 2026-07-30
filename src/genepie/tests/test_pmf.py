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
from ..exceptions import GenesisValidationError


class TestPmfAnalysis(CustomTestCase):
    """Regression tests for pmf_analysis against the CLI reference output."""

    # Same tolerance the CLI regression harness uses for this test case, see
    # tests/regression_test/test_analysis/test_pmf_analysis/pathcv/config.ini
    TOLERANCE = 0.1

    def _traj(self, name):
        return str(self.TEST_ROOT / "test_analysis" / "trajectories"
                   / "umbrella_path" / name)

    def _run(self, **overrides):
        kwargs = dict(
            cvfile=self._traj("{}.pathcv"),
            weightfile=self._traj("{}.weight"),
            distfile=self._traj("{}.pathdist"),
            nreplica=16,
            dimension=1,
            temperature=300.0,
            cutoff=0.04,
            grids=((1.0, 16.0, 50),),
            band_width=(0.3,),
            is_periodic=(False,),
        )
        kwargs.update(overrides)
        return genesis_exe.pmf_analysis(**kwargs)

    def test_pmf_analysis_matches_cli_reference(self):
        """PMF must reproduce the values the CLI pmf_analysis writes."""
        res = self._run()

        ref = np.loadtxt(
            self.TEST_ROOT / "test_analysis" / "test_pmf_analysis"
            / "pathcv" / "pmf.dat.ref")

        # ref columns: bin center, standard PMF, Gaussian-kernel PMF
        self.assertEqual(ref.shape[0], res.cv.shape[0])
        self.assertTrue(np.all(np.isfinite(res.cv)))
        self.assertTrue(np.all(np.isfinite(res.pmf)))
        self.assertTrue(np.all(np.isfinite(res.pmf_gaussian)))
        self.assertAlmostEqual(ref[:, 0], res.cv, delta=self.TOLERANCE)
        self.assertAlmostEqual(ref[:, 1], res.pmf, delta=self.TOLERANCE)
        self.assertAlmostEqual(ref[:, 2], res.pmf_gaussian, delta=self.TOLERANCE)

    def test_pmf_analysis_memory_matches_file(self):
        """In-memory arrays must give the same PMF as the equivalent files.

        Uses a single replica (no dist cutoff) so the memory path and the
        file path see identical samples.
        """
        cv = np.loadtxt(self._traj("1.pathcv"))[:, 1]
        weight = np.loadtxt(self._traj("1.weight"))[:, 1]

        file_res = genesis_exe.pmf_analysis(
            cvfile=self._traj("{}.pathcv"),
            weightfile=self._traj("{}.weight"),
            nreplica=1,
            dimension=1,
            temperature=300.0,
            grids=((1.0, 16.0, 50),),
            band_width=(0.3,),
            is_periodic=(False,),
        )
        mem_res = genesis_exe.pmf_analysis(
            cv=cv,
            weight=weight,
            temperature=300.0,
            grids=((1.0, 16.0, 50),),
            band_width=(0.3,),
            is_periodic=(False,),
        )
        self.assertAlmostEqual(file_res.cv, mem_res.cv, delta=1e-8)
        self.assertAlmostEqual(file_res.pmf, mem_res.pmf, delta=1e-6)
        self.assertAlmostEqual(
            file_res.pmf_gaussian, mem_res.pmf_gaussian, delta=1e-6)

    def test_pmf_analysis_2d_runs(self):
        """A 2-D reweighted PMF returns a finite matrix with matching axes."""
        rng = np.random.default_rng(0)
        phi = rng.uniform(-170, 170, size=2000)
        psi = rng.uniform(-170, 170, size=2000)
        weight = rng.uniform(0.5, 1.5, size=2000)

        res = genesis_exe.pmf_analysis(
            cv=np.column_stack([phi, psi]),
            weight=weight,
            temperature=300.0,
            grids=((-180.0, 180.0, 37), (-180.0, 180.0, 37)),
            band_width=(15.0, 15.0),
            is_periodic=(True, True),
            box_size=(360.0, 360.0),
        )
        self.assertEqual(res.pmf.shape, (36, 36))
        self.assertEqual(res.cv1.shape[0], 36)
        self.assertEqual(res.cv2.shape[0], 36)
        self.assertTrue(np.all(np.isfinite(res.pmf)))
        self.assertAlmostEqual(float(res.pmf.min()), 0.0, delta=1e-9)

    def test_pmf_analysis_requires_cv(self):
        with self.assertRaises(GenesisValidationError):
            genesis_exe.pmf_analysis(
                grids=((1.0, 16.0, 50),), band_width=(0.3,))

    def test_pmf_analysis_rejects_both_cv_and_cvfile(self):
        with self.assertRaises(GenesisValidationError):
            genesis_exe.pmf_analysis(
                cv=np.zeros(10), cvfile="x{}.cv",
                grids=((1.0, 16.0, 50),), band_width=(0.3,))


if __name__ == "__main__":
    unittest.main()
