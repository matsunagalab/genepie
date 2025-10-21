"""
GENESIS Python Interface

This package provides Python bindings for GENESIS analysis tools.
"""

from .s_molecule import SMolecule
from .s_trajectories import STrajectories, STrajectoriesArray
from . import genesis_exe
from .genesis_exe import LibGenesis
from . import ctrl_files 

__version__ = "1.0.0"
__all__ = [
    "SMolecule",
    "STrajectories", 
    "STrajectoriesArray",
    "genesis_exe", 
    "LibGenesis", 
    "ctrl_files",
]
