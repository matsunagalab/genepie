"""Enumeration constants shared with the GENESIS Fortran layer.

The integer values must stay in sync with the Fortran parameters they mirror,
so this module is the single place that defines them. The ``*_MAP`` dictionaries
translate the control-file spellings accepted by the CLI into those integers,
and the name lists used for validation are derived from the maps so that the two
cannot drift apart.
"""


# Mirrors src/analysis/libana/trajectory_str.fpp
class TrjFormat:
    PDB = 1
    AMBER = 2
    DCD = 3
    GROMACS = 4
    CHARMM_RST = 5
    NAMD_RST = 6


# Mirrors src/analysis/libana/trajectory_str.fpp
class TrjType:
    COOR = 1
    COOR_BOX = 2


# Mirrors src/lib/fitting_str.fpp (FittingMethodTypes)
class FittingMethod:
    NO = 1
    TR_ROT = 2
    TR = 3
    TR_ZROT = 4
    XYTR = 5
    XYTR_ZROT = 6


# Mirrors src/analysis/libana/pbc_correct.fpp
class PBCCMode:
    NO = 1
    MOLECULE = 2


TRJ_FORMAT_MAP = {
    "PDB": TrjFormat.PDB,
    "AMBER": TrjFormat.AMBER,
    "DCD": TrjFormat.DCD,
    "GROMACS": TrjFormat.GROMACS,
    "CHARMM_RST": TrjFormat.CHARMM_RST,
    "NAMD_RST": TrjFormat.NAMD_RST,
}

TRJ_TYPE_MAP = {
    "COOR": TrjType.COOR,
    "COOR+BOX": TrjType.COOR_BOX,
}

FITTING_METHOD_MAP = {
    "NO": FittingMethod.NO,
    "TR+ROT": FittingMethod.TR_ROT,
    "TR": FittingMethod.TR,
    "TR+ZROT": FittingMethod.TR_ZROT,
    "XYTR": FittingMethod.XYTR,
    "XYTR+ZROT": FittingMethod.XYTR_ZROT,
}

PBCC_MODE_MAP = {
    "NO": PBCCMode.NO,
    "MOLECULE": PBCCMode.MOLECULE,
}

TRJ_FORMATS = list(TRJ_FORMAT_MAP)
TRJ_TYPES = list(TRJ_TYPE_MAP)
FITTING_METHODS = list(FITTING_METHOD_MAP)
PBC_CORRECTS = list(PBCC_MODE_MAP)
