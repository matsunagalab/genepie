# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------
import os
import unittest
import numpy as np
from .. import genesis_exe
from ..custom_test_case import CustomTestCase
from ..exceptions import (
    GenesisFortranNotSupportedError,
    GenesisValidationError,
)

# 61 umbrella windows spaced 3 degrees apart along a periodic dihedral.
NREPLICA = 61
CONSTANT = (0.06092,) * NREPLICA
REFERENCE = tuple(3.0 * i for i in range(NREPLICA))


def load_fene_reference(path):
    """Load a fene reference file as a 2-D (n_replica, n_blocks) array."""
    ref = np.loadtxt(path)
    if ref.ndim == 1:
        ref = ref[:, np.newaxis]
    return ref


class MbarUmbrellaMixin:
    """Shared regression checks for the umbrella_1d MBAR dataset.

    Kept separate from CustomTestCase so that subclassing it from another test
    module does not make unittest collect the same tests twice.
    """

    # The MBAR iteration converges to 1e-8 and the CLI regression harness
    # accepts 0.01 across platforms, so 1e-6 catches regressions while leaving
    # room for BLAS-dependent differences in the solver.
    TOLERANCE = 1.0e-6
    NBLOCKS = None
    REFERENCE_DIR = "umbrella_1d"

    def reference_path(self):
        return (self.TEST_ROOT / "test_analysis" / "test_mbar_analysis"
                / self.REFERENCE_DIR / "fene.dat.ref")

    def run_mbar(self, **overrides):
        kwargs = dict(
            cvfile=str(self.TEST_ROOT / "test_analysis" / "trajectories"
                       / "umbrella_1d" / "{}.dat"),
            nreplica=NREPLICA,
            input_type="US",
            dimension=1,
            temperature=300.0,
            target_temperature=300.0,
            tolerance=1E-08,
            rest_function=(1,),
            grids=((-1.0, 181.0, 81),),
            constant=(CONSTANT,),
            reference=(REFERENCE,),
            is_periodic=(True,),
            box_size=(360.0,),
        )
        if self.NBLOCKS is not None:
            kwargs["nblocks"] = self.NBLOCKS
        kwargs.update(overrides)
        return genesis_exe.mbar_analysis(**kwargs)

    def test_mbar_analysis_matches_cli_reference(self):
        """Free energies must reproduce what the CLI mbar_analysis writes."""
        fene = self.run_mbar()

        ref = load_fene_reference(self.reference_path())

        self.assertEqual(ref.shape, fene.shape)
        self.assertTrue(np.all(np.isfinite(fene)))
        self.assertAlmostEqual(ref, fene, delta=self.TOLERANCE)

    def test_mbar_analysis_leaves_cwd_clean(self):
        """The Fortran writer must not drop fene.dat into the cwd."""
        self.assertFalse(os.path.exists("fene.dat"))
        self.run_mbar()
        self.assertFalse(os.path.exists("fene.dat"))

    def test_mbar_analysis_rejects_dcd_input(self):
        """DCD input is unsupported because molecules are never defined."""
        with self.assertRaises(GenesisFortranNotSupportedError):
            self.run_mbar(dcdfile="whatever.dcd")

    def test_mbar_analysis_requires_cvfile(self):
        with self.assertRaises(GenesisValidationError):
            self.run_mbar(cvfile=None)

    def test_mbar_analysis_rejects_missing_cvfile(self):
        """A nonexistent cvfile must raise, not abort the whole process.

        The Fortran open_file() calls exit(1) on a missing input file, so the
        wrapper has to check existence up front.
        """
        with self.assertRaises(GenesisValidationError):
            self.run_mbar(cvfile="/no/such/dir/{}.dat")


class TestMbarAnalysis1D(MbarUmbrellaMixin, CustomTestCase):
    pass


if __name__ == "__main__":
    unittest.main()
