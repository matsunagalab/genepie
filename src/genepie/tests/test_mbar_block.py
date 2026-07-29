# --------------------------------------------
if __name__ == "__main__" and __package__ is None:
    import sys, pathlib
    pkg_dir = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(pkg_dir.parent.parent))
    __package__ = "genepie.tests"
# --------------------------------------------
import unittest
from ..custom_test_case import CustomTestCase
from .test_mbar_1d import MbarUmbrellaMixin


class TestMbarAnalysisBlock(MbarUmbrellaMixin, CustomTestCase):
    """Same dataset as the 1D case, split into blocks for error estimation."""

    NBLOCKS = 5
    REFERENCE_DIR = "umbrella_block"


if __name__ == "__main__":
    unittest.main()
