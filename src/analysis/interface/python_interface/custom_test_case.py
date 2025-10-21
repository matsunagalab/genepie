import os
import pathlib
import unittest
from typing import Optional, Union
import numpy as np
from .ctrl_files import TrajectoryParameters
from . import genesis_exe
from .s_molecule import SMolecule
from .s_trajectories import STrajectories, STrajectoriesArray


class CustomTestCase(unittest.TestCase):
    """"""
    TEST_ROOT = pathlib.Path("../../../../tests/regression_test")
    PDB_PATH = pathlib.Path("BPTI_ionize.pdb")
    PSF_PATH = pathlib.Path("BPTI_ionize.psf")
    TRJ_PATH = pathlib.Path("BPTI_run.dcd")

    def assertAlmostEqualNumpyNDArray(
            self, expected, actual,
            rtol: Optional[float] = None,
            atol: Optional[float] = None,
            msg: Optional[str] = None,
            ):
        if rtol is None:
            rtol = 1.e-7
        if atol is None:
            atol = 0
        if not np.allclose(expected, actual, rtol=rtol, atol=atol):
            standard_msg = (f"{expected}.coords != {actual}.coords "
                            + f"within rtol={rtol}, atol={atol}")
            self.fail(self._formatMessage(msg, standard_msg))

    def assertObjectsAlmostEqual(self, expected, actual,
                                 places: Optional[int] = None,
                                 msg: Optional[str] = None,
                                 delta: Optional[float] = None,
                                 except_attr: set[str] = set()):
        edict = expected.__dict__
        adict = actual.__dict__

        if msg is not None:
            msg = msg + " : "
        else:
            msg = ""
        self.assertEqual(set(edict.keys()), set(adict.keys()),
                         msg + "Object attribute names do not match.")
        for key in edict:
            if key in except_attr:
                continue
            ev = edict[key]
            av = adict[key]
            if (
                    (isinstance(ev, np.ndarray)
                        and (ev.dtype == np.float64)
                        and isinstance(av, np.ndarray)
                        and (av.dtype == np.float64))
                    or (isinstance(ev, float) and isinstance(av, float))
             ):
                self.assertAlmostEqual(ev, av, places,
                                       msg + f"{key} is not equal.", delta)
            else:
                self.assertEqual(ev, av, msg + f"{key} is not equal.")

    def assertAlmostEqualSTrajectories(
            self, e_trj: STrajectories, a_trj: STrajectories,
            places: Optional[int] = None,
            msg: Optional[str] = None,
            delta: Optional[float] = None):
        if (places is None) and (delta is None):
            places = 4
        if (
                not isinstance(e_trj, STrajectories)
                or not isinstance(a_trj, STrajectories)
        ):
            self.fail("Both arguments must be instances of STrajectories.")
        if msg is not None:
            msg = msg + " : Strajectories"
        else:
            msg = "Strajectories"
        self.assertObjectsAlmostEqual(e_trj, a_trj, places, msg, delta,
                                      except_attr={"c_obj", "_mem_owner"})

    def assertAlmostEqualSMolecule(
            self, e_mol: STrajectories, a_mol: STrajectories,
            places: Optional[int] = None,
            msg: Optional[str] = None,
            delta: Optional[float] = None):
        if (places is None) and (delta is None):
            places = 4
        if (
                not isinstance(e_mol, SMolecule)
                or not isinstance(a_mol, SMolecule)
        ):
            self.fail("Both arguments must be instances of SMolecule.")
        if msg is not None:
            msg = msg + " : SMolecule"
        else:
            msg = "SMolecule"
        self.assertObjectsAlmostEqual(
                e_mol, a_mol, places, msg, delta,
                except_attr={
                    "atom_coord", "num_atoms_fep", "num_bonds_fep",
                    "num_angles_fep", "num_dihedrals_fep", "num_impropers_fep",
                    "num_cmaps_fep", "bond_list_fep", "angl_list_fep",
                    "dihe_list_fep", "impr_list_fep", "cmap_list_fep",
                    "id_singleA", "id_singleB", "fepgrp", "fepgrp_bond",
                    "fepgrp_angl", "fepgrp_dihe", "fepgrp_cmap"})

    def assertAlmostEqual(self, expected, actual,
                          places: Optional[int] = None,
                          msg: Optional[str] = None,
                          delta: Optional[float] = None):
        if isinstance(expected, np.ndarray) and isinstance(actual, np.ndarray):
            if delta is None:
                if places is not None:
                    atol = 10**(-places)
                else:
                    atol = None
            else:
                atol = delta
            if (expected.dtype == np.float64) and (actual.dtype == np.float64):
                self.assertAlmostEqualNumpyNDArray(
                        expected, actual, rtol=0.0, atol=atol, msg=msg)
                return
        elif (isinstance(expected, STrajectories)
              and isinstance(actual, STrajectories)):
            self.assertAlmostEqualSTrajectories(
                    expected, actual, places, msg, delta)
            return
        elif (isinstance(expected, SMolecule)
              and isinstance(actual, SMolecule)):
            self.assertAlmostEqualSMolecule(
                    expected, actual, places, msg, delta)
            return
        if (places is None) and (delta is None):
            places = 7
        super().assertAlmostEqual(
                expected, actual, places, msg, delta)

    def assertEqualNdarray(self, expected, actual, msg: str = None):
        self.assertTrue(np.array_equal(expected, actual), msg=msg)

    def assertEqual(self, expected, actual, msg: str = None):
        if isinstance(expected, np.ndarray) and isinstance(actual, np.ndarray):
            self.assertEqualNdarray(expected, actual, msg)
        else:
            super().assertEqual(expected, actual, msg)

    @staticmethod
    def create_traj_by_genesis(
            dcd: Union[str, bytes, os.PathLike],
            pdb: Union[str, bytes, os.PathLike] = '',
            top: Union[str, bytes, os.PathLike] = '',
            gpr: Union[str, bytes, os.PathLike] = '',
            psf: Union[str, bytes, os.PathLike] = '',
            ref: Union[str, bytes, os.PathLike] = '',
            fit: Union[str, bytes, os.PathLike] = '',
            prmtop: Union[str, bytes, os.PathLike] = '',
            ambcrd: Union[str, bytes, os.PathLike] = '',
            ambref: Union[str, bytes, os.PathLike] = '',
            grotop: Union[str, bytes, os.PathLike] = '',
            grocrd: Union[str, bytes, os.PathLike] = '',
            groref: Union[str, bytes, os.PathLike] = '',
            ) -> tuple[STrajectoriesArray, SMolecule]:
        mol = SMolecule.from_file(
                pdb=pdb, top=top, gpr=gpr, psf=psf, ref=ref, fit=fit,
                prmtop=prmtop, ambcrd=ambcrd, ambref=ambref,
                grotop=grotop, grocrd=grocrd, groref=groref)
        trajs, subset_mol = genesis_exe.crd_convert(
                mol,
                traj_params=[
                    TrajectoryParameters(
                        trjfile=dcd,
                        md_step=10,
                        mdout_period=1,
                        ana_period=1,
                        repeat=1,
                        ),
                    ],
                trj_format="DCD",
                trj_type="COOR+BOX",
                trj_natom=0,
                selection_group=["all", ],
                fitting_method="NO",
                fitting_atom=1,
                check_only=False,
                pbc_correct="NO",
                )
        _ = subset_mol
        return (trajs, mol)
